
#include "sched_tools.h"
using std::cout ;
using std::cerr ;
using std::endl ;
#include "distribute.h"

namespace Loci {
#define MAX_CHUNKS 100
#define RTS_HLP    9930
#define RTR_HLP    9920
#define RTS_RES    9910
#define RTR_RES    9900
unsigned char *buf;  /*Send and receive data*/
int position=0;        
int size=0;            /*Size of buffer*/


  //global variables
  rule_implP rp1;
  rule_implP local_comp1; 
  fact_db *local_facts1;
  variableSet inputs1;
  variableSet outputs1;
  entitySet exec_set1;
  fact_db *facts1;
  MPI_Comm procGrp;
 
  /*To allocate inputs and outputs over the iterate space*/
void Allocate_func(){
    
    for(variableSet::const_iterator vi=inputs1.begin();vi!=inputs1.end();++vi) {
      storeRepP sp = local_facts1->get_variable(*vi) ;
      sp->allocate(interval(0,exec_set1.size()-1)) ;     
    }
    for(variableSet::const_iterator vi=outputs1.begin();vi!=outputs1.end();++vi) {
      storeRepP sp = local_facts1->get_variable(*vi) ;
      sp->allocate(interval(0,exec_set1.size()-1)) ;
    }     
}
  /*To deallocate inputs and outputs*/
void Deallocate_func(){
  
    
    //deallocate the temporaries
    for(variableSet::const_iterator vi=inputs1.begin();vi!=inputs1.end();++vi) {
      storeRepP sp = local_facts1->get_variable(*vi) ;
      sp->allocate(EMPTY) ;
    } 
    for(variableSet::const_iterator vi=outputs1.begin();vi!=outputs1.end();++vi) {
      storeRepP sp = local_facts1->get_variable(*vi) ;
      sp->allocate(EMPTY) ;
    }
}

  // Transfer inputs to dest 
void SendInput (int tStart, int tSize, int dest, int tag,MPI_Comm procGrp) {  
  
    // Pack inputs from facts
    position = 0 ;
    size = 0 ;
    for(variableSet::const_iterator vi=inputs1.begin();vi!=inputs1.end();++vi) {
      storeRepP s_ptr = facts1->get_variable(*vi) ;
      size += s_ptr->pack_size(interval(tStart,tStart+tSize-1)) ;
    }
   
    buf = new unsigned char[size] ;
    entitySet myent=interval(tStart,tStart+tSize-1);
    for(variableSet::const_iterator vi=inputs1.begin();vi!=inputs1.end();++vi) {
      storeRepP s_ptr = facts1->get_variable(*vi) ;
      s_ptr->pack(buf,position,size,myent) ;
    } 
    MPI_Status tStatus;
    //Send Inputs
    int send_array[2];
    for(int k=0;k<=1;k++){
      send_array[k]=0;
    } 
    send_array[0]=size;
    send_array[1]=tSize;    
    MPI_Sendrecv(send_array, 2, MPI_INT, dest,RTS_HLP,NULL,0,MPI_INT,dest,RTR_HLP, procGrp,&tStatus);    
    MPI_Rsend(buf,size,MPI_UNSIGNED_CHAR,dest,tag,procGrp);   
    delete [] buf;
}

  // Receive inputs from sender 
void ReceiveInput (int rcvStart,int *rcvSize,int src,int tag,MPI_Comm procGrp) { 
 
    MPI_Status tStatus;
    int recvArray[2];
    for(int l=0;l<=1;l++){
      recvArray[l]=0;
    }    
    int size=0;
    MPI_Recv(recvArray, 2, MPI_INT, src, RTS_HLP, procGrp,&tStatus);  
    MPI_Send(NULL,0,MPI_INT,src,RTR_HLP,procGrp); 
    size=recvArray[0];
    *rcvSize=recvArray[1];
    buf = new unsigned char[size];
   
    MPI_Recv (buf,size,MPI_UNSIGNED_CHAR,src, tag, procGrp,&tStatus);   
    
    // unpack inputs into local facts
    position = 0 ;
    for(variableSet::const_iterator vi=inputs1.begin();vi!=inputs1.end();++vi) {     
      storeRepP s_ptr = local_facts1->get_variable(*vi) ;
      s_ptr->unpack(buf,position,size,sequence(interval(rcvStart,rcvStart+*rcvSize-1))) ;
    }
    
    // Delete allocated temporaries
    delete[] buf ;
 
}
//Send outputs to dest
void SendOutput (int tStart, int tSize, int dest,int tag,MPI_Comm procGrp) { 

    // Pack outputs
    position = 0 ;
    size = 0 ;
    for(variableSet::const_iterator vi=outputs1.begin();vi!=outputs1.end();++vi){
      storeRepP s_ptr = local_facts1->get_variable(*vi) ;
      size += s_ptr->pack_size(interval(0,tSize-1)) ;
    }
  
    buf = new unsigned char[size] ;
    entitySet myent2=interval(0,tSize-1);
    for(variableSet::const_iterator vi=outputs1.begin();vi!=outputs1.end();++vi) {
      storeRepP s_ptr = local_facts1->get_variable(*vi) ;
      s_ptr->pack(buf,position,size,myent2) ;
    }
    MPI_Status tStatus;
    //Send outputs 
    MPI_Sendrecv(&size, 1, MPI_INT, dest,RTS_RES,NULL,0,MPI_INT,dest,RTR_RES,procGrp,&tStatus);
    MPI_Rsend(buf,size,MPI_UNSIGNED_CHAR, dest, tag, procGrp);

    delete [] buf;

  
}
  //Receive outputs from src
void ReceiveOutput (int rcvStart, int rcvSize, int src, int tag,MPI_Comm procGrp) { 
 
    MPI_Status tStatus;
    size=0;
    MPI_Recv(&size, 1, MPI_INT, src, RTS_RES, procGrp, &tStatus);
    MPI_Send(NULL,0,MPI_INT,src,RTR_RES,procGrp);
    buf = new unsigned char[size];
   
    MPI_Recv (buf, size, MPI_UNSIGNED_CHAR, src, tag, procGrp, &tStatus);
    // unpack outputs into facts
    position = 0 ;     
   
    for(variableSet::const_iterator vi=outputs1.begin();vi!=outputs1.end();++vi) {
      storeRepP s_ptr = facts1->get_variable(*vi) ;
      s_ptr->unpack(buf,position,size,sequence(interval(rcvStart,rcvStart+rcvSize-1))) ; 
    }
 
    // Delete allocated temporaries
    delete[] buf ;
}
  
void workCompute(int start, int size,int Signal){
  rule_implP rp2;
  rp2=rp1; 
  //Compute from local facts
  if(Signal==1){   
       rp2=local_comp1;  
  }  
  rp2->compute(sequence(interval(start,start+size-1)));
}

void ExecuteLoop (void (*workCompute) (int,int,int),void (*SendInput) (int,int,int,int,MPI_Comm),void (*ReceiveInput) (int,int *,int,int,MPI_Comm),void (*SendOutput) (int,int,int,int,MPI_Comm),void (*ReceiveOutput) (int,int,int,int,MPI_Comm),void (*Allocate_func) (),void (*Deallocate_func) (),int method,int *yMap,double *stats,int *chunkMap,MPI_Comm procGrp) ;


dynamic_schedule_rule::dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds)  {

  rp = fi.get_rule_implP() ; //get rule from rule database 
  rule_tag = fi ;            //store rule tag in rule_tag
  local_compute1 = rp->new_rule_impl() ; //another instance of rule 
  entitySet in = rule_tag.sources() ; //inputs from rule rhs
  outputs = rule_tag.targets() ;      //outputs as in rhs  
  exec_set = eset ;   
  
  // Setup local facts input variables (types only no allocation)
  for(variableSet::const_iterator vi=in.begin();vi!=in.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    if((store_ptr != 0) && store_ptr->RepType() == Loci::STORE) {	  
      inputs += *vi ;
      // we make it an intensional fact, since it is local fact_db
      local_facts.create_intensional_fact(*vi,store_ptr->new_store(EMPTY)) ;
    } else {
      local_facts.create_intensional_fact(*vi,facts.get_variable(*vi)) ;
    }
  }
  // Setup local facts output variables 
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    local_facts.create_intensional_fact(*vi,store_ptr->new_store(EMPTY)) ;
  }
  
  // Initialize both functions for remote and local execution.
  local_compute1->initialize(local_facts) ;
  rp->initialize(facts) ;  
  
  // MPI_Comm procGrp=MPI_COMM_WORLD;
  MPI_Comm_dup(MPI_COMM_WORLD, &procGrp);
}


void dynamic_schedule_rule::execute(fact_db &facts) { 
   
    //Assignment to global variables
  rp1=rp;
  local_comp1=local_compute1;
  exec_set1=exec_set;
  facts1=&facts;
  local_facts1=&local_facts;
  inputs1=inputs;
  outputs1=outputs; 
  extern int method; 
  
  int *yMap=new int[2*Loci::MPI_processes]; 
  for(int i = 0; i < Loci::MPI_processes; i++) {
    yMap[2*i] = 0;
    yMap[2*i+1] = 0;
  }
   //Broadcast start,size to all procs
  int sendstart=0;
  sendstart = exec_set.Min();
  int *rstart = new int[Loci::MPI_processes];
  for(int q = 0; q < Loci::MPI_processes; q++) {
    rstart[q] = 0;
  }
   MPI_Allgather(&sendstart, 1, MPI_INT, rstart, 1, MPI_INT, MPI_COMM_WORLD);
  int sendsize=0;
  sendsize = exec_set.Max()-sendstart+1;
  int *rsize = new int[Loci::MPI_processes];
  for(int n = 0; n < Loci::MPI_processes; n++) {
    rsize[n] = 0;
  }
   MPI_Allgather(&sendsize, 1, MPI_INT, rsize, 1, MPI_INT, MPI_COMM_WORLD);
 
  for(int i = 0; i < Loci::MPI_processes; i++) {
    yMap[2*i] = rstart[i];
    yMap[2*i+1] = rsize[i];
  }
  delete [] rsize;
  delete [] rstart;
  
  double *stats=new double[4*Loci::MPI_processes];
  for(int i = 0; i < Loci::MPI_processes; i++) {
   stats[4*i+0] = 0;
   stats[4*i+1] =0;
   stats[4*i+2] = 0;
   stats[4*i+3] =0;
  }
  int *chunkMap=new int[3*MAX_CHUNKS];
  for(int i = 0; i < MAX_CHUNKS; i++) {
   chunkMap[3*i+0]=0;
   chunkMap[3*i+1]=0;
   chunkMap[3*i+2]=0;
  }     
  double t1=0.0;
  double t2=0.0;
  double t3=0.0;
    
  Loci::ExecuteLoop(&workCompute,&SendInput,&ReceiveInput,&SendOutput,&ReceiveOutput,&Allocate_func,&Deallocate_func,method,yMap,stats,chunkMap,procGrp);
  
  delete [] yMap;        
  delete [] stats;
  delete [] chunkMap;

}
   

void dynamic_schedule_rule::Print(std::ostream &s) const {

 s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl ; 
  
}


}

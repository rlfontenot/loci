#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <algorithm>
#include <math.h>
#include <mpi.h>

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;
using std::cin ; 

#define DIM 3

/*
 * retained for historical reasons: the number of bits in an attribute value:
 * effectively the order of a curve
 */
#define		NUMBITS			32

/*
 * the number of bits in a word used to store an hcode (or in an element of
 * an array that's used)
 */
#define		WORDBITS		32

typedef struct 
{
  double x;
  double y;
  double z;
  double node_num;
}Node;

typedef unsigned long U_long;

typedef struct 
{
  U_long int_x;
  U_long int_y;
  U_long int_z;
  U_long node_num;
}Node_int;

typedef struct 
{
  U_long key1; 
  U_long key2;
  U_long key3;
  U_long key_num;
}Key;

typedef struct
{
  int nodes;
  int cell_left;
  int cell_right;
}Cell;


typedef struct {
  U_long hcode[DIM];
}Hcode;

typedef Hcode Point;

// Mask for 3-Dimensions
const U_long g_mask[] = {4,2,1};


void dec2bin(unsigned long t , int b[]);
unsigned long bin2dec(int c[]);
void Converter(Node_int *int_cod , int x, Key *new2old, int counter);
void Rand_key(int n, Key *new2old );

// functions to compute the Hilbert Order 
Hcode H_encode(Point pt);


inline bool compare(unsigned long *a1,unsigned long *a2)
{
  return a1[0]<a2[0] || (a1[0]==a2[0] && a1[1]<a2[1]) || ( a1[0]==a2[0] &&a1[1]==a2[1] && a1[2]<a2[2]);
}


inline bool compare_key(Key a1,Key a2)
{
  return a1.key1 < a2.key1 || (a1.key1 == a2.key1 && a1.key2 < a2.key2) || (a1.key1 == a2.key1 && a1.key2 == a2.key2 && a1.key3 < a2.key3 );  }


inline bool compare_splitter(Key *a1,unsigned long *a2)
{
  return a1->key1 < a2[0] || (a1->key1 == a2[0] && a1->key2 < a2[1]) || (a1->key1 == a2[0] && a1->key2 == a2[1] && a1->key3 < a2[2] ); 
}


inline bool Gcompare(Key *a1,unsigned long *a2)
{
  return a1->key1 > a2[0] || (a1->key1 == a2[0] && a1->key2 > a2[1]) || ( a1->key1==a2[0] && a1->key2 == a2[1] && a1->key3 > a2[2]) ||( a1->key1 == a2[0] && a1->key2 == a2[1] && a1->key3 == a2[2]) ;
}


inline bool Ecompare_splitter(Key *a1,unsigned long *a2)
{
  return a1->key1 < a2[0] || (a1->key1 == a2[0] && a1->key2 < a2[1]) || (a1->key1 == a2[0] && a1->key2 == a2[1] && a1->key3 < a2[2]) || ( a1->key1 == a2[0] && a1->key2 == a2[1] && a1->key3 == a2[2] ); 
}


inline bool compare_node_num(Node a1,Node a2)
{
  return a1.node_num < a2.node_num ;
}


inline bool compare_coordkey(double *a1,double *a2)
{
  return a1[3]<a2[3] ;
}

inline bool min_compare(int *a1,int *a2)
{
  return a1[1]<a2[1];
}


int Server(int ac, char *av[],int , int );
void Client(int , int ,int ,char *av[]);

/***************************************************/
int main(int ac, char *av[]) {
 
  int myId ;
  int numProcessors ; 

  MPI_Init(&ac,&av) ;
 

  /* Get the number of processors and my processor identification */
  MPI_Comm_size( MPI_COMM_WORLD, &numProcessors) ;
  MPI_Comm_rank( MPI_COMM_WORLD, &myId) ;

  if(myId == 0) {
    // Processor 0 runs the server code
    Server(ac,av, myId, numProcessors) ;    
  } else {
    // all other processors run the client code.
    Client(myId, numProcessors,ac,av);
  }

  // All MPI programs must call this before exit  
  MPI_Finalize() ;
  return 0 ; 
}  



 
int Server(int ac, char *av[], int Id, int numProcessors)
{
   
  if(ac <= 1) {
    cerr << "usage:  xdropt <xdrfile>" << endl;
    exit(-1);
  }
  
  int myrank, npes; 
  MPI_Status stat;
  myrank = Id;
  npes = numProcessors;
  
  int ln;
  ln = strlen(av[1]);
  char *new_filename = new char[ln+5];
  
  strcpy(new_filename,"new_");
  strcat(new_filename,av[1]);
  
  
  FILE *FP = fopen(av[1], "r") ;
  if(FP == NULL) {
    cerr << "can't open " << av[1] <<  endl ;
    return(-1);
  }
   
   
  int Info[8];
   
  XDR xdr_handle ;
  xdrstdio_create(&xdr_handle, FP, XDR_DECODE) ;
  
  int nfaces,ncells , num_nodes;
   
  //Reading the grid info 
  for(int i = 0; i < 8; i++)
    xdr_int(&xdr_handle, &Info[i]) ;
   
  num_nodes = Info[3] ;
  nfaces = Info[4];
  ncells = Info[5];


  MPI_Bcast(&Info[3], 3, MPI_INT, 0, MPI_COMM_WORLD);
   
  cout<<"intializing......"<<endl;
  cout<<"nodes "<<num_nodes<<endl;
  cout<<"nfaces "<<nfaces<<endl;
  cout<<"ncells "<<ncells<<endl;
  

  // Computes number of nodes to send to each PE
   
  int nodes_send,nodes_sendlast,num,ctr,temp,count,flag=0;  
  num = num_nodes/npes;
   
  if( num_nodes % npes == 0)
    {
      nodes_send = num;
      nodes_sendlast = num;
    }
  else
    {
      nodes_send = num;
      nodes_sendlast = num_nodes - num*(npes-1);
    } 
   
  // Computes number of faces to send to each PE
  int nfaces_send, nfaces_sendfirst;
  num = nfaces/npes;
   
  if( nfaces % npes == 0)
    {
      nfaces_send = num;
      nfaces_sendfirst = num;
    }
  else
    {
      nfaces_send = num;
      nfaces_sendfirst = nfaces - num*(npes-1);
    } 
   
  
  double  max = 0, min = 0;
  int counter = 0;  
  double *temp_ptmp = new double[nodes_sendlast*3];
  Node *ptmp = new Node[nodes_send];   

  //Reading in nodes for the source PE
  // Reading in the xyz coordinates & finding max & min
  for(int i = 0; i < nodes_send ; i++)
    {
       
      xdr_double(&xdr_handle,&ptmp[i].x) ;
      xdr_double(&xdr_handle,&ptmp[i].y) ;
      xdr_double(&xdr_handle,&ptmp[i].z) ;      
       
      if( ptmp[i].x > max )
        max=ptmp[i].x;
       
      if( ptmp[i].y > max )
        max=ptmp[i].y;   
       
      if( ptmp[i].z > max )
        max=ptmp[i].z;
       
      if(i == 0)
        min = ptmp[i].x;
                       
      if( ptmp[i].x < min )
        min = ptmp[i].x;
       
      if( ptmp[i].y < min )
        min = ptmp[i].y;   
       
      if( ptmp[i].z < min )
        min = ptmp[i].z;    
     
      ptmp[i].node_num = counter + i;
    }
 
  int tag1 = 1;
  //  int tp = nodes_sendlast*3;
  ctr = 0;
  count = 0;
   
  /*********************************************************/
  // SENDING NODES TO THE PE AND FINDING MAX & MIN
   
  if(npes > 1)
    while(ctr<npes-1)
      {
        if(ctr==npes-2)
          {
            temp=nodes_sendlast;
            flag = 1;
          }
        else if(ctr==0)
          {
            temp=nodes_send;
            flag=0;
          }
        else
          {
            temp=nodes_send;
            flag=1;
          }
	 
        for(int i = 0; i < temp; i++)
          {
            for(int j = 0; j < 3; j++)
              {
                xdr_double(&xdr_handle,&temp_ptmp[i*3 + j]) ;
                if( temp_ptmp[i*3 + j] > max )
                  max = temp_ptmp[i*3 + j];  
		 
                if( temp_ptmp[i*3 + j] < min)
                  min = temp_ptmp[i*3 + j];           		 
              }
          }      
	 
        count = temp * 3;
	 
        MPI_Send(&tag1, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD );       
        MPI_Send(&temp_ptmp[0], count, MPI_DOUBLE, ctr+1, 0, MPI_COMM_WORLD );   
        ctr++ ;  
        count = 0;
      
      }
    
  delete temp_ptmp;


  // BROADCASTING MAX & MIN

  MPI_Bcast(&min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
 
  // CONVERTING THE X,Y,Z COORDINATES TO INTEGERS
  
  double diff = max - min;
  double result;
  typedef unsigned long length[4]; 
  Node_int *int_cod = new Node_int[nodes_send]; 
   
  counter = 0;
  unsigned long max_int = ~0L ; 
  
  
  for(int i = 0; i < nodes_send; i++)
    {
      result = (ptmp[i].x - min) / diff ;
      int_cod[i].int_x = (unsigned long)(result * max_int);
       
      result = (ptmp[i].y - min) / diff ;
      int_cod[i].int_y = (unsigned long)(result * max_int);
	
      result = (ptmp[i].z - min) / diff ;
      int_cod[i].int_z = (unsigned long)(result * max_int);
	
      int_cod[i].node_num = i + counter;
    }




  // Declares memory for the morton , hilbert or random key 
  Key *new2old = new Key[nodes_send];
    
    
  // Selects the morton or the random function 
    
  if(ac > 2 && !strcmp("random",av[2]))
    {
      cout<<"Random ordering selected........."<<endl;
      for(int i=0;i<nodes_send;i++)
        Rand_key(i+counter,&new2old[i]);
    }
  else if(ac > 2 && !strcmp("morton",av[2])) 
    {
      // Performs the bit interleaving and o/p in new array
      cout<<"morton ordering selected........."<<endl;
      // new2old contains the key and counter gives an index to the key
      for(int i = 0; i < nodes_send; i++)
        Converter(&int_cod[i], i, &new2old[i], counter + i);
    }
  else 
    {
      // Performs the bit interleaving and o/p in new array
      // new2old contains the key and counter gives an index to the key
      cout<<"hilbert ordering selected........."<<endl;
   
     
      Point pt;
      Hcode ht;
      for(int i = 0; i < nodes_send; i++)
        {
          pt.hcode[0] = int_cod[i].int_x;      
          pt.hcode[1] = int_cod[i].int_y; 
          pt.hcode[2] = int_cod[i].int_z;     
      
          ht = H_encode(pt); 
      
          new2old[i].key3 = ht.hcode[0]  ;      
          new2old[i].key2 = ht.hcode[1]  ;      
          new2old[i].key1 = ht.hcode[2]  ;      
          new2old[i].key_num = counter + i;
        }
    }

  delete int_cod;

  /*******************************************************************/
  //          START OF SAMPLE SORT FOR THE NODES
  // Locally sorting the new2old array
   
  std::sort(new2old, new2old + nodes_send, compare_key);
  cout<<"Re-ordering the nodes......."<<endl; 

  // FINDING THE SPLITTERS

  count = npes * (npes-1);
  length *splitter = new length[npes-1];  
  length *bcast_splitter = new length[npes];   
  unsigned long *allpick = new unsigned long[count*4];

  for(int i=1;i<npes;i++)
    {
      
      splitter[i-1][0]= new2old[i*(nodes_send/npes)].key1;
      splitter[i-1][1]= new2old[i*(nodes_send/npes)].key2;
      splitter[i-1][2]= new2old[i*(nodes_send/npes)].key3;
      splitter[i-1][3]= new2old[i*(nodes_send/npes)].key_num;
       
    }
  
  count = (npes-1) * npes;
  unsigned long **gather;
  gather = new unsigned long*[count];
 
  for(int i=0;i<count;i++)
    gather[i] = new unsigned long[4];
    
 
  // Collecting the splitters
  count = (npes-1)*4;
  MPI_Gather( &splitter[0][0], count, MPI_UNSIGNED_LONG, &allpick[0], count, MPI_UNSIGNED_LONG,0, MPI_COMM_WORLD );
 

  count =  npes*(npes-1);
  int t=0; 

  for(int i=0;i<count;i++)
    {
      for(int j=0;j<4;j++)
        {
          gather[i][j]=allpick[t++];
        }      
    }
  
  // Sorting the splitters
  std::sort( gather, gather + count , compare);

 
  for(int i=1;i<npes;i++)
    {
      for(int j=0;j<4;j++)
        {
          bcast_splitter[i-1][j]= gather[i*npes-1][j];
        }
    }
  
  unsigned long x=0;
  
  for(int j = 0;j<4;j++)
    {
      bcast_splitter[npes-1][j] = ~(x<<31);
    }
  
 
  // SENDING THE SPLITTERS TO ALL the PE 
  ctr=0;
  count = npes*4; 
  
  while(ctr < npes-1)
    {
      
      MPI_Send(&bcast_splitter[0], count, MPI_UNSIGNED_LONG, npes-ctr-1, 0, MPI_COMM_WORLD );  
      ctr++;
    }

 
  //Compute the number of elements that belong to each bucket  
  int *scounts = new int[npes];
  
  for(int i=0;i<npes;i++)
    scounts[i] = 0;

  
  for(int i=0;i<nodes_send;i++)
    {
      for(int j=0; j< npes; j++)
        {   
          if(j==0)
            {
              if(compare_splitter(&new2old[i],bcast_splitter[j])) 
                {
                  scounts[j]++; 
                  break; 
                }
            } 
          else 
            {   
              if( j!=npes-1 && compare_splitter(&new2old[i],bcast_splitter[j]) && Gcompare(&new2old[i], bcast_splitter[j-1]))
                {
                  scounts[j]++;
                  break;    
                }
              else if(Ecompare_splitter(&new2old[i],bcast_splitter[j]) && Gcompare(&new2old[i], bcast_splitter[j-1]))
                {
                  scounts[j]++;
                  break;  
                }
            }
        }
    }
     
  
  int  y = 0; 
  // Multiplying by 4 because each key has x,y,z co-ordinate and index
  
  for(int i=0;i<npes;i++)
    scounts[i]= scounts[i] * 4;
  
  // Determine the starting location of each buckets element 
 
  int *sdisp = new int[npes];  
  sdisp[0] = 0;
  
  for(int i=1; i<npes; i++)
    {
      sdisp[i]= sdisp[i-1]+scounts[i-1];
    }
  
  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive
   
  int *rcounts = new int[npes];
  MPI_Alltoall(&scounts[0], 1, MPI_INT, &rcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
  
  
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
  
  int *rdisp = new int[npes];
  rdisp[0] = 0;
  
  for(int i=1; i<npes; i++)
    rdisp[i]= rdisp[i-1]+rcounts[i-1];
    
  int total; 
  // In case of one processor total = nodes_send    
   
  if(npes > 1)
    total = (rdisp[npes-1]+rcounts[npes-1])/4;
  else
    total = nodes_send;
 
 
  unsigned long *elemn = new unsigned long[nodes_send*4];
  unsigned long *sort_elemn = new unsigned long[total*4];
    
  for(int i=0;i<total*4;i++)
    sort_elemn[i]=0;

  t = 0;
  for(int i=0;i<nodes_send;i++)
    {       
      elemn[t++]= new2old[i].key1;  
      elemn[t++]= new2old[i].key2;  
      elemn[t++]= new2old[i].key3;  
      elemn[t++]= new2old[i].key_num;  
    }
  
 
  // Each PE Sends & recieve element using MPI_Alltoallv
  
  MPI_Alltoallv(&elemn[0], &scounts[0], &sdisp[0], MPI_UNSIGNED_LONG, &sort_elemn[0], &rcounts[0], &rdisp[0], MPI_UNSIGNED_LONG, MPI_COMM_WORLD );

  // THIS ARRAY CONTAINS THE KEY AFTER THE SAMPLE SORT 
  
  Key *PE_key = new Key[total];
   
  for(int i=0;i<total;i++)
    {   
      PE_key[i].key1 = sort_elemn[i*4 + 0];
      PE_key[i].key2 = sort_elemn[i*4 + 1];
      PE_key[i].key3 = sort_elemn[i*4 + 2];
      PE_key[i].key_num = sort_elemn[i*4 + 3];
    }
  
  // LOCAL SORT OF THE KEY 
  std::sort( PE_key, PE_key  + total , compare_key);

  delete elemn;
  delete sort_elemn;
  delete bcast_splitter;
  delete splitter;
  delete gather;

  /********************************************************************/
  // Gathering the new node number from all the other processors and making an old2new array
  
  // array new_nodes contain the index after the sample sort   
  unsigned long *new_nodes = new  unsigned long[total];   
  
  for(int i=0;i<total;i++)
    {
      new_nodes[i] = PE_key[i].key_num;
    }

  int *recvbuff = new int[npes];
  MPI_Allgather(&total,1, MPI_INT, &recvbuff[0], 1, MPI_INT, MPI_COMM_WORLD );
 
  int *rdispla = new int[npes];
  rdispla[0] = 0;
  
  for(int i=1;i<npes;i++)
    rdispla[i]= rdispla[i-1]+recvbuff[i-1];
  
  int sum = rdispla[npes - 1] + recvbuff[npes - 1];

  if(npes == 1)
    sum = num_nodes;

  
  // array new_node_num serves as the new numbering of the nodes 

  unsigned long *new_node_num = new  unsigned long[sum];  
  
  MPI_Allgatherv(&new_nodes[0],total, MPI_UNSIGNED_LONG, &new_node_num[0], &recvbuff[0], &rdispla[0], MPI_UNSIGNED_LONG, MPI_COMM_WORLD );

  delete new_nodes;
 
  // This generates a old2new numbering of the nodes 
  int *old2new = new int[sum];
 
  // In case of one processor new2old is the new node numbering 
  if(npes > 1)
    {
      for(int i=0;i<num_nodes;i++)
        old2new[new_node_num[i]] = i;
    }
  else
    {
      // In Case of one processor
      for(int i=0;i<num_nodes;i++)
        old2new[new2old[i].key_num] = i; 
    }

  // Reordering the xyz coordinates index acc to old2ew numbering 
  
  for(int i=0;i<nodes_send;i++)
    ptmp[i].node_num = (double)old2new[(int)ptmp[i].node_num];
     
  std::sort(ptmp, ptmp + nodes_send, compare_node_num);   
 
  /*****************************************************************/
  // A final SAMPLE SORT is required based on the new node numbering 
  // to align the nodes in selected Ordering 
   
  y=0;
  int *split = new int[npes];
  int *scount = new int[npes];
  
  for(int i=0;i<npes;i++)
    scount[i] = 0 ;
  
   
  for(int i=0;i<npes-1;i++)    
    split[i] = num_nodes/npes * (i+1) ;
    
  split[npes-1] = num_nodes ;



  for(int i=0;i<nodes_send;i++)
    {
      for(int j=0; j< npes; j++)
        {   
          if(j==0)
            {
              if(ptmp[i].node_num < split[j])
                {
                  scount[j]++;  
                }
            } 
          else 
            {
              if( (ptmp[i].node_num < split[j]) && ((ptmp[i].node_num > split[j-1]) ||(ptmp[i].node_num == split[j-1]))) 
                {
                  scount[j]++; 
                }
            }
        }
    }


  for(int i=0;i<npes;i++)
    scount[i]= scount[i]*4;
    
  // Determine the starting location of each buckets element 

  int *sdispl = new int[npes];
  sdispl[0] = 0;
  
  for(int i=1; i<npes; i++)
    sdispl[i]= sdispl[i-1]+scount[i-1];
     
  
  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive
   
  int *rcount = new int[npes];  
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);
  

  //Based on rcounts determine where in the local array the data
  //from each PE will be stored  
 
  int *rdispl = new int[npes];
  rdispl[0] = 0;

  for(int i=1; i<npes; i++)
    rdispl[i]= rdispl[i-1]+rcount[i-1];
     

  if(npes>1)
    total = (rdispl[npes-1]+rcount[npes-1])/4;
  else
    total = nodes_send;

  t=0;
  Node *new_ptmp = new Node[total];
 
  double *seq = new  double[nodes_send*4];
  double *temp_seq = new  double[total*4];

  for(int i=0;i<nodes_send;i++)
    {       
      seq[i*4 + 0] = ptmp[i].x;  
      seq[i*4 + 1] = ptmp[i].y;  
      seq[i*4 + 2] = ptmp[i].z;  
      seq[i*4 + 3] = ptmp[i].node_num;  
    }
   
  // Each PE Sends & recieve element using MPI_Alltoallv
   
  MPI_Alltoallv(&seq[0], &scount[0], &sdispl[0], MPI_DOUBLE, &temp_seq[0], &rcount[0], &rdispl[0], MPI_DOUBLE, MPI_COMM_WORLD );
  cout<<"where r we "<<endl;
   
  for(int i=0;i<total;i++)
    {       
      new_ptmp[i].x = temp_seq[i*4 + 0];  
      new_ptmp[i].y = temp_seq[i*4 + 1];  
      new_ptmp[i].z = temp_seq[i*4 + 2];  
      new_ptmp[i].node_num = temp_seq[i*4 + 3];  
    }
   
  // This variable is used later to send the new_ptmp elements to the PE 0 where the coordinates are written back to the xdr file 
 
  int num_coord = total;
   
  // Local sort of the nodes which aligns nodes in Morton ordering        
  std::sort(new_ptmp  ,new_ptmp + total ,compare_node_num );
  /*
    FILE *TFP = fopen("testing_hilbert.dat", "w") ;
    if(TFP == NULL) {
    cerr << "can't open " << "testing_grid.xls" <<  endl ;
    return(-1);
    }
     
    for(int i=0; i<total ;i++)
    fprintf( TFP, " %f %f %f \n" , new_ptmp[i].x, new_ptmp[i].y, new_ptmp[i].z);
   
  */

  delete temp_seq;
  delete seq;
  delete ptmp;
   
  /************************************************************************/  
  // Reading in the cellnum info for PE 0
  // If there is only one PE 
   
  cout<<"Reading in  the cellnum......."<<endl; 

  int off=0;
  if(npes > 1)
    off = 0;
  else
    off = 1;

  // Contains one more element in case there is only one processor
  int *cellnum_off = new int[ nfaces_send*3 + off];
  t=0;
   
  for(int i = 0; i < nfaces_send ; i++)
    {
      for(int j = 0; j < 3; j++)
        {
          xdr_int(&xdr_handle,&cellnum_off[t++]) ;      
        }
    }
   
  if(off == 1)
    xdr_int(&xdr_handle,&cellnum_off[t]) ;      

  //SENDING THE CELLS TO OTHER PROCESSORS
  
  int *cellnum = new int[nfaces_sendfirst*3+1];
   
  ctr=0;
  t=0;
  temp=0;
       
  while(ctr < npes-1)
    {
      if(ctr == npes-2)       
        temp = nfaces_sendfirst;
      else
        temp = nfaces_send;
       
      for(int i = 0; i < temp ; i++)
        {
          for(int j = 0; j < 3; j++)
            {
              xdr_int(&xdr_handle,&cellnum[t++]) ;      
            }
        }
     
      if(ctr == npes-2)  
        {
          xdr_int(&xdr_handle,&cellnum[t]);        
          count=temp*3+1;
        }
      else
        count = temp*3;

      MPI_Send(&cellnum[0], count, MPI_INT, ctr+1, 0, MPI_COMM_WORLD );  
      
      ctr++;  
      count = 0;
      t = 0; 
    }
   
  delete cellnum;
   
  int offset;
 
  // Receving the first offset so that next PE can calculate its last offset
  if( npes > 1 )
    MPI_Recv(&offset, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD ,&stat ); 


  Cell *size_cellnum = new Cell[nfaces_send];
  
  // Tells how many nodes in a face and the cellnum
  for(int i = 0; i <nfaces_send  ; i++)
    {
      size_cellnum[i].cell_left = cellnum_off[i*3 + 1];	
      size_cellnum[i].cell_right = cellnum_off[i*3 + 2];	            

      if(i != nfaces_send-1)
        size_cellnum[i].nodes = cellnum_off[(i+1)*3] - cellnum_off[i*3] ;
      else
        {
          if(npes > 1)
            size_cellnum[i].nodes = offset - cellnum_off[i*3] ; 
          else
            size_cellnum[i].nodes = cellnum_off[(i+1)*3] - cellnum_off[i*3] ; 
        }
    }
   


  // Total number of nodes in all the faces of PE 0
  int num_facenodes;
   
  if(npes > 1)
    num_facenodes = offset - cellnum_off[0];
  else
    num_facenodes = cellnum_off[nfaces_send * 3] - cellnum_off[0];

  int *recvcount = new int[npes-1]; 
  ctr=0;
   
  // REceiving the count of faces each PE has to recieve
  
  if(npes > 1)
    while(ctr<npes-1) 
      {      
        MPI_Recv(&recvcount[ctr],1 , MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);      
        ctr++;
      }
    
  /********************************************************************/
  // READING IN THE FACES FOR THE PE 0 
  cout<<"reading in faces........"<<endl;
   
  int *temp_node = new int[num_facenodes];
   
  for(int i = 0; i < num_facenodes ; i++)
    {
      xdr_int(&xdr_handle,&temp_node[i]) ;      
    }
   
  int maximum = recvcount[0];
 
  for(int i=1;i<npes-1;i++)
    if( recvcount[i]>maximum)
      maximum =recvcount[i];

  
  // SENDING THE FACES TO THE PE
   
  int *temp_nodes = new int[maximum];
  ctr=0;
  t=0;
  while(ctr<npes-1)
    {
      for(int i = 0; i < recvcount[ctr] ; i++)
        {
          xdr_int(&xdr_handle,&temp_nodes[t++]) ;      	   
        }
       
      count= recvcount[ctr];
       
      MPI_Send(&temp_nodes[0], count, MPI_INT, ctr+1, 0, MPI_COMM_WORLD );  
      ctr++;  
      count=0;
      t=0; 
    }
  
  delete temp_nodes;

 
  // changing the node number to match the new node numbering
  for(int i=0; i<num_facenodes; i++)
    {
      temp_node[i]= old2new[temp_node[i]];
    }
   
   
  // Contains the faces acc to the new node numbering 
  int **new_faces;
  new_faces = new int*[nfaces_send];
   
  for(int i=0;i< nfaces_send ;i++)
    new_faces[i] = new int[size_cellnum[i].nodes];
   
  t=0;
   
  // Preparing the new faces 
  for(int i=0;i<nfaces_send ;i++)
    {
      for(int j=0;j<size_cellnum[i].nodes;j++)
        {
          new_faces[i][j] = temp_node[t++];
        }
    }
   
   
  int minimum; 
  int **min_faces;
  min_faces = new int*[nfaces_send];
   
  for(int i = 0; i <nfaces_send ; i++) 
    min_faces[i]= new int[2];
   
     
  /****************************************************************/
  //    Reordering the faces in the morton order
  
  cout<<"Re-ordering the faces......."<<endl;   


  if(ac > 3 && !strcmp("boundary", av[3]))
    {      
      //finding the minimum of the face and storing it into the min_faces
      for(int i=0;i<nfaces_send  ;i++)
        {
          minimum = new_faces[i][0];
          for(int j=1;j<size_cellnum[i].nodes;j++)
            {
              if( new_faces[i][j] < minimum )
                minimum =new_faces[i][j]; 	
            }
	   
          int tagval = size_cellnum[i].cell_left;
          if(tagval >  size_cellnum[i].cell_right)
            tagval =  size_cellnum[i].cell_right ;
          if(tagval < 0)
            minimum = tagval * num_nodes + minimum ;
	   
          min_faces[i][0] = i;
          min_faces[i][1] = minimum;
	   
        }
    }
  else
    {
      //finding the minimum of the face and storing it into the min_faces
      cout<<"boundary "<<endl;      
      for(int i=0;i<nfaces_send  ;i++)
        {
          minimum = new_faces[i][0];
          for(int j=1;j<size_cellnum[i].nodes;j++)
            {
              if( new_faces[i][j] < minimum )
                minimum =new_faces[i][j]; 	
            }
	  
          min_faces[i][0] = i;
          min_faces[i][1] = minimum;
	   
        }
    }    

   
  // Sorting the faces acc to the minimum node number
  std::sort(min_faces, min_faces + nfaces_send , min_compare);
   
  // Allocating the memory to the changed faces 
  int **changed_faces;
  changed_faces = new int*[nfaces_send];
   
  for(int i=0;i<nfaces_send;i++)
    changed_faces[i] = new int[size_cellnum[min_faces[i][0]].nodes];
   
   
  // Re - Aligns faces
  for(int i=0;i<nfaces_send;i++)
    {
      for(int j = 0;j < size_cellnum[min_faces[i][0]].nodes ; j++)
        {
          changed_faces[i][j] = new_faces[min_faces[i][0]][j];
        }
    }
  
  // SAMPLE SORT ON THE FACES 
  count = npes*(npes-1);
   
  // Selecting the splitters for the Sample sort 
   
  int *splitter1 = new int[npes-1];  
  int *bcast_splitter1 = new int[npes];   
  int *allpick1 = new int[count];
   
   
  for(int i=1;i<npes;i++)
    splitter1[i-1]= min_faces[ i * (nfaces_send / npes)][1];

  count = npes-1;
   
  MPI_Gather( &splitter1[0], count, MPI_INT, &allpick1[0], count, MPI_INT,0, MPI_COMM_WORLD );
   
   
  count =  npes*(npes-1);
      
  std::sort(allpick1 ,allpick1  + count );
 
   
  for(int i=1;i<npes;i++)
    {
      bcast_splitter1[i-1] = allpick1[i*npes-1];     
    }
   
  bcast_splitter1[npes-1] = num_nodes;

  // Broadcasting splitters
 
  MPI_Bcast(&bcast_splitter1[0], npes, MPI_INT, 0, MPI_COMM_WORLD);
   
   
  //Compute the number of elements that belong to each bucket  
   
  int *scount1 = new int[npes];
   
  for(int i=0;i<npes;i++)
    scount1[i]=0;
     

  for(int i=0;i<nfaces_send;i++)
    {
      for(int j=0; j< npes; j++)
        {   
          if(j==0)
            {
              if(min_faces[i][1] < bcast_splitter1[j]) 
                {
                  scount1[j]++;  
                }
            } 
          else 
            {
              if( (min_faces[i][1]<bcast_splitter1[j]) &&( min_faces[i][1]> bcast_splitter1[j-1] || min_faces[i][1] == bcast_splitter1[j-1]))
                {
                  scount1[j]++; 
                }
            }
        }
    }
   
  
  // Sending the sorted min_faces to other PE's so that they can perform local sort on basis of that
  // Determine the starting location of each buckets element 
   
  int *sdisp_min = new int[npes];
  sdisp_min[0] = 0;
   
  for(int i=1; i<npes; i++)    
    sdisp_min[i]= sdisp_min[i-1]+scount1[i-1];
 
   
  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive
   
  int *rcount_min = new int[npes];
   
  MPI_Alltoall(&scount1[0], 1, MPI_INT, &rcount_min[0], 1, MPI_INT, MPI_COMM_WORLD);

  
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored

  int *rdisp_min = new int[npes];   
  rdisp_min[0] = 0;
   
  for(int i=1; i<npes; i++) 
    rdisp_min[i]= rdisp_min[i-1]+rcount_min[i-1];
     
   
  total = 0; 
  total = (rdisp_min[npes-1] + rcount_min[npes-1]);
   
  int *pick_minfaces = new int[total];
  int *temp_min = new int[nfaces_send];
   
  for(int i=0; i < nfaces_send; i++)
    temp_min[i] = min_faces[i][1];
     
   
  // Each PE Sends & recieve element using MPI_Alltoallv
  MPI_Alltoallv(&temp_min[0], &scount1[0], &sdisp_min[0], MPI_INT, &pick_minfaces[0], &rcount_min[0], &rdisp_min[0], MPI_INT, MPI_COMM_WORLD );
   
  // gving an index to the collected minfaces so that the faces can be re-ordered acc top that
   
  int **gather_minfaces; 
  gather_minfaces = new int*[total];
  for(int i=0;i<total;i++)
    gather_minfaces[i] = new int[2];  
   
  for(int i=0;i<total;i++)
    {
      gather_minfaces[i][0]=i;  
      gather_minfaces[i][1]=pick_minfaces[i];  
    } 
   
   
  //          LOCAL SORT
  // Sorting the faces acc to the minimum node number  
  std::sort(gather_minfaces,gather_minfaces + total , min_compare);
   
  delete temp_min;
  delete pick_minfaces;
   
 
  // SENDING THE FACES TO THE PROCESSORS ACC TO THE MIN_FACES ARRAY                             
  t=0;
  total=0;
  int *total1 = new int[npes];
  for(int i=0;i<npes;i++)
    {
      total1[i]=0;
    }
     
  // total1[] gives the total number of nodes in the faces that are to be sent  
     
  for(int i=0;i<npes;i++)
    {
      for(int j=0;j<scount1[i];j++)
        { 
          total1[i] = total1[i] + size_cellnum[min_faces[t++][0]].nodes;
        }
    }
     
          
  // Determine the starting location of each buckets element 
     
  int *sdisp1 = new int[npes];
  sdisp1[0] = 0;
     
  for(int i=1; i<npes; i++)
    sdisp1[i]= sdisp1[i-1] + total1[i-1];
       
     
  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive
     
  int *rcount1 = new int[npes];     
  MPI_Alltoall(&total1[0], 1, MPI_INT, &rcount1[0], 1, MPI_INT, MPI_COMM_WORLD);
     
  
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
     
     
  int *rdisp1 = new int[npes];
  rdisp1[0] = 0;
     
  for(int i=1; i<npes; i++)      
    rdisp1[i]= rdisp1[i-1]+rcount1[i-1];
        
     
  total = (rdisp1[npes-1]+rcount1[npes-1]);
     
  int *elemn1 = new int[num_facenodes];
  int *sort_elemn1 = new int[total];
     
  // This variable is used later to send the faces to the first PE to write the faces back to the xdr file
  int num_finalfaces = total;
     
     
  t=0;
     
  for(int i=0;i<nfaces_send;i++)
    {
      for(int j= 0;j < size_cellnum[min_faces[i][0]].nodes ;j++)
        {
          elemn1[t++]= changed_faces[i][j] ;   
        }
    }
      
     
  // Each PE Sends & recieve element using MPI_Alltoallv
     
  MPI_Alltoallv(&elemn1[0], &total1[0], &sdisp1[0], MPI_INT, &sort_elemn1[0], &rcount1[0], &rdisp1[0], MPI_INT, MPI_COMM_WORLD );
     
     
  delete elemn1;
     
     
  //  SENDING THE SIZE AND CELLNUM OF THE FACES  
     
  for(int i=0;i<npes;i++)
    {
      scount1[i]= scount1[i]*3;	 
    }
     
  // For sending the face size to the PE's since faces can have different sizes, rcount2 tells the number of faces received from each PE
     
  int *rcount2 = new int[npes];
     
  MPI_Alltoall(&scount1[0], 1, MPI_INT, &rcount2[0], 1, MPI_INT, MPI_COMM_WORLD);
     
     
  // size array has the size of the faces , which is needed by othere PE's 
     
  int *size = new int[nfaces_send*3];
     
  t=0;
  for(int i=0;i<nfaces_send;i++)
    {
      size[t++]= size_cellnum[min_faces[i][0]].nodes;
      size[t++]= size_cellnum[min_faces[i][0]].cell_left;
      size[t++]= size_cellnum[min_faces[i][0]].cell_right;
    }
     
     
  sdisp1[0] = 0;
  for(int i=1; i<npes; i++)
    {
      sdisp1[i]= 0;
    }
     
     
  for(int i=1; i<npes; i++)
    {
      sdisp1[i]= sdisp1[i-1]+scount1[i-1];
    }
   
  rdisp1[0] = 0;
  for(int i=1; i<npes; i++)
    {
      rdisp1[i]= 0;
    } 
     
  for(int i=1; i<npes; i++)
    {
      rdisp1[i]= rdisp1[i-1]+rcount2[i-1];
    } 
     
     
  total = (rdisp1[npes-1]+rcount2[npes-1]);
     
  int *recv_size = new int[total];
         
  // Each PE Sends & recieve element using MPI_Alltoallv
     
  MPI_Alltoallv(&size[0], &scount1[0], &sdisp1[0], MPI_INT, &recv_size[0], &rcount2[0], &rdisp1[0], MPI_INT, MPI_COMM_WORLD );
     
  // Re-aligning the cellnum and size array acc to the new local sort of the gather_minface
  // faces[][] contains the faces received from other PE's but we need to do the local sort on min_faces again
     
  int **faces;
  faces = new int*[total/3];
     
  for(int i=0;i<total/3;i++)
    {
      faces[i] = new int[recv_size[i*3]];
    }
     
     
  t=0;
   
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[i*3] ;j++)
        {
          faces[i][j] = sort_elemn1[t++];
        }
    }
     
     
  // Allocating memory to faces acc to the local sort of the minfaces
     
  int **final_faces;
  final_faces = new int*[total/3];
     
  for(int i=0;i<total/3;i++)
    {
      final_faces[i] = new int[recv_size[(gather_minfaces[i][0])*3]];
    }
     
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[(gather_minfaces[i][0])*3] ;j++)
        {
          final_faces[i][j] = faces[(gather_minfaces[i][0])][j];
        }
    }
     
   
  t=0;
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[(gather_minfaces[i][0])*3] ;j++)
        {
          sort_elemn1[t++]=final_faces[i][j];
        }
    }
     
  delete final_faces;
  delete faces;

  // Sending the first offset to the next processor for the offset that is along with cellnum
     
  sum=0;
        
  for(int i=0;i<total/3;i++)
    sum = sum + recv_size[i*3];
       
  // sum_all is used to calculate the offset
 
  int *sum_all = new int[npes];
  MPI_Allgather(&sum, 1, MPI_INT, &sum_all[0], 1, MPI_INT, MPI_COMM_WORLD);
    
  // Giving new offset to the cellnum acc to the Morton ordering
  cout<<"reordering the cell number "<<endl;
   
  if(npes == 1)
    total = nfaces_send * 3 + 1;       
     
  int off_set = 0;
  int *new_array = new int[total ];
     
  t=0;
     
  // Contains the cellnumber in the morton ordering with the new offset
  if(npes > 1)
    {
      for(int i = 0;i < total/3; i++)
        {
          new_array[t++] =  off_set ;     
          new_array[t++] =  recv_size[(gather_minfaces[i][0])*3+1] ;
          new_array[t++] =  recv_size[(gather_minfaces[i][0])*3+2] ;
          off_set = off_set + recv_size[(gather_minfaces[i][0])*3] ;
        }
    }
  else
    {	
      // For single processor  
      for(int i = 0;i < nfaces_send; i++)
        { 
          new_array[t++] =  off_set ;     
          new_array[t++] =  size_cellnum[min_faces[i][0]].cell_left;
          new_array[t++] =  size_cellnum[min_faces[i][0]].cell_right;
          off_set = off_set + size_cellnum[min_faces[i][0]].nodes;
        }
    }


  if(npes == 1)
    new_array[t] = off_set;
  
  /****************************************************************/
  //       Reordering & renumbering the cell number 
  /***************************************************************/
   
  ctr=1;
  flag=1;
  t=0;
   
  int *new_cellnum = new int[ncells];

  for(int i=0;i<ncells;i++)
    {
      new_cellnum[i]=-1;
    }
   

  // Cells reordered acc to their first visit by faces
  
  for(int i=0; i < total/3; i++)
    {    
      if(new_array[i*3+1]>0)
        {
          if(new_cellnum[new_array[i*3+1]-1]==-1)
            new_cellnum[new_array[i*3+1]-1]=ctr++;
        }  
      
      if(new_array[i*3+2]>0)
        {
          if(new_cellnum[new_array[i*3+2]-1]==-1)
            new_cellnum[new_array[i*3+2]-1]=ctr++;
        }  
    }

  // sending the counter and new_cellnum array which is updated by each PE
 
  if(npes > 1)
    {
      MPI_Send(&ctr, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD ); 
       
      MPI_Send(&new_cellnum[0], ncells, MPI_INT, myrank+1, 0, MPI_COMM_WORLD );  
      MPI_Barrier( MPI_COMM_WORLD);   
       
      MPI_Bcast( &new_cellnum[0], ncells, MPI_INT, npes-1, MPI_COMM_WORLD);  
    }
  
  

  for(int i=0; i < total/3; i++)
    {    
      if(new_array[i*3+1]>0)
        {
          new_array[i*3+1] = new_cellnum[new_array[i*3+1]-1];
        }
      
      if(new_array[i*3+2]>0)
        {
          new_array[i*3+2] = new_cellnum[new_array[i*3+2]-1];
        }
      
    }
  
   
  
  cout<<"writing out to the new xdr file...."<<endl;
  
  // Writing out into new xdr file 
  
  FILE *NFP = fopen(new_filename, "w") ;
  if(NFP == NULL) {
    cerr << "can't open " << "new.xdr" <<  endl ;
    return(-1);
  }
  
  
  XDR xdr_handler2 ;
  xdrstdio_create(&xdr_handler2, NFP, XDR_ENCODE) ;
  
 
  for(int i = 0; i < 8; i++)
    xdr_int(&xdr_handler2, &Info[i]) ;
    
  cout<<"writing the coord........"<<endl; 
  for(int i=0;i<num_coord;i++)
    {
      /*    for(int j=0; j < 3 ;j++)
            {
            xdr_double(&xdr_handler2,&new_ptmp[i][j]) ;
            // cout<<new_ptmp[i][j]<<endl;
            }*/

      xdr_double(&xdr_handler2, &new_ptmp[i].x) ;
      xdr_double(&xdr_handler2, &new_ptmp[i].y) ;
      xdr_double(&xdr_handler2, &new_ptmp[i].z) ;

    } 


  // Receiving the count of the number of co-ordinates that other PE's are sending to write back to the xdr file

  ctr=0;
  int  *recv_numcoord = new int[npes-1];   
 
  if(npes > 1)
    {
      while(ctr<npes-1)
        {
     
          MPI_Recv(&recv_numcoord[ctr],1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  

          ctr++;        
        }


      maximum = recv_numcoord[0];

  
      for(int i=1;i<npes-1;i++)
        if( recv_numcoord[i] > maximum )
          maximum =  recv_numcoord[i];
  
      // Writing the co-ordinates back to the xdr file 
  
      double *temp_newptmp = new double[maximum];  
      ctr=0;
      tag1=0;

  
      while(ctr<npes-1)
        {
          
          MPI_Recv(&tag1 , 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  
     
          if(tag1 == 1) 
            MPI_Recv(&temp_newptmp[0] ,recv_numcoord[ctr], MPI_DOUBLE, ctr+1, 0, MPI_COMM_WORLD, &stat);  
      
          for(int i=0; i <recv_numcoord[ctr] ;i++)
            {
              xdr_double(&xdr_handler2,&temp_newptmp[i]);
              // cout<<temp_newptmp[i]<<endl; 
            }

          ctr++;        
        }
  


      delete temp_newptmp;
    }  

  // Writing the cellnum's of the PE 0  
  cout<<"writing the cellnum........"<<endl; 

  for(int i=0;i < total ;i++)
    {
      xdr_int(&xdr_handler2, &new_array[i]) ;
      //	 cout<<new_array[i]<<endl; 
    }
  
  // Receiving the count of the number of cellnum that other PE's are sending to write back to the xdr file
  
 
  ctr=0;
  int  *recv_cellnum = new int[npes-1];   


  while(ctr<npes-1)
    {
     
      MPI_Recv(&recv_cellnum[ctr],1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  

      ctr++;        
    }

  if(npes > 1)
    {

      maximum = recv_cellnum[0];

  
      for(int i=1;i<npes-1;i++)
        if( recv_cellnum[i] > maximum )
          maximum =  recv_cellnum[i];
    }
  // Writing the cellnum back to the xdr file
  
  int *temp_cellnum = new int[maximum];  
  ctr=0;
  tag1=0;
    
 
  while(ctr<npes-1)
    {

      
      MPI_Recv(&tag1 , 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  
     
      if(tag1 == 1)            
        MPI_Recv(&temp_cellnum[0] ,recv_cellnum[ctr], MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  

    
      for(int i=0; i <recv_cellnum[ctr] ;i++)
        {
          xdr_int(&xdr_handler2,&temp_cellnum[i]);
          //	cout<<temp_cellnum[i]<<endl; 
        }

      ctr++;        
      tag1 = 0;
    }
  
  delete temp_cellnum;
  delete recv_cellnum;
  delete recv_numcoord;  
       
  cout<<"writing the faces........"<<endl;

  if(npes == 1)
    {
      for(int i=0;i < nfaces_send; i++)
        {
          for(int j=0;j<size_cellnum[min_faces[i][0]].nodes;j++)
            {
              xdr_int(&xdr_handler2, &changed_faces[i][j]);
            }
        }
    }
  else
    {
      for(int i = 0; i < num_finalfaces; i++)
        {
          xdr_int(&xdr_handler2, &sort_elemn1[i]) ;
          // 	cout<<sort_elemn1[i]<<endl; 
        }
    }

  // Receiving the count of the number of faces that other PE's are sending to write back to the xdr file

  ctr=0;
  int  *recv_numfaces = new int[npes-1];   

  
  while(ctr<npes-1)
    {
     
      MPI_Recv(&recv_numfaces[ctr],1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  
      ctr++;        
    }


  if(npes >1)
    {
      maximum = recv_numfaces[0];
      
      for(int i=1;i<npes-1;i++)
        if( recv_numfaces[i] > maximum )
          maximum =  recv_numfaces[i];
    }
  // Writing the faces back to the xdr file
  
  int *temp_faces = new int[maximum];  
  ctr=0;
  tag1=0;


  while(ctr<npes-1)
    {
       
      MPI_Recv(&tag1 , 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  
     
      if(tag1 == 1)            
        MPI_Recv(&temp_faces[0] ,recv_numfaces[ctr], MPI_INT, ctr+1, 0, MPI_COMM_WORLD, &stat);  
    
      for(int i=0; i <recv_numfaces[ctr] ;i++)
        {
          xdr_int(&xdr_handler2,&temp_faces[i]);
          //  cout<<temp_faces[i]<<endl; 
        }
      ctr++;
      tag1=0;        
    }
  
   
  cout<<"closing xdr handle........thanks"<<endl;

  //   fclose(TFP);    
  xdr_destroy(&xdr_handle) ;
  xdr_destroy(&xdr_handler2) ;
  fclose(NFP) ;
  fclose(FP) ;
  return 0;
  
}


/************************************************************/
/********************START OF CLIENT*************************/
/************************************************************/

void Client(int id, int numProcessors,int ac,char *av[]) {
  
  int myrank, npes; 
  MPI_Status stat;
  int Info[3];
  
  npes = numProcessors;
  myrank = id;
  
  int num, nodes_send ,nodes_sendlast,num_nodes;
  
  MPI_Bcast( Info, 3, MPI_INT, 0, MPI_COMM_WORLD);  
  num_nodes = Info[0];

 
  int ORDER = 1;
  int k=1;
  while(k < num_nodes || k == num_nodes ) 
    {
      k = k << 1;
      ORDER++;
    }  





  // Calculating the number of nodes to send & receive to each PE 

  num = Info[0]/npes;
  
  if( Info[0] % npes == 0)
    {
      nodes_send = num;
      nodes_sendlast = num;
    }
  else 
    {
      nodes_send = num;
      // nodes_sendlast = Info[0] - num*(npes-1);
      nodes_sendlast = num;
      if(myrank==npes-1)   
	nodes_send = Info[0] - num*(npes-1);
    } 

  int nfaces_send, nfaces_sendfirst;
  int nfaces = Info[1]; 

 
  num = nfaces/npes;

  if( nfaces % npes == 0)
    {
      nfaces_send = num;
      nfaces_sendfirst = num;
    }
  else
    {
      nfaces_send = num;
      if(myrank==npes-1)
	nfaces_send  = nfaces - num*(npes-1);

      //      nfaces_sendfirst = nfaces - num*(npes-1);
    } 

 
  double  max=0,min=0;
  typedef double triple[3];
  double *temp_ptmp = new double[nodes_send*3];
  Node *ptmp = new Node[nodes_send];

  int tag1=0;
  int tp = 0;
  tp = nodes_send*3;
   
  // Receving the nodes from PE 0
  MPI_Recv( &tag1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);  

  if(tag1 == 1)
    MPI_Recv(&temp_ptmp[0],tp, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);     

  MPI_Bcast(&min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

  // Converting the x,y,z coordinates to integers
  
  double diff = max - min;
  
  typedef unsigned long length[4];
  // length *int_cod = new length[nodes_send]; 
 
  Node_int *int_cod = new Node_int[nodes_send];  
  Key *new2old = new Key[nodes_send];
  

  int counter;  
  double result;

  if(myrank==npes-1)  
    counter = myrank * nodes_sendlast;
  else
    counter = myrank * nodes_send;
  

  for(int i=0;i<nodes_send ;i++)
    {
    
      ptmp[i].x =  temp_ptmp[i*3 + 0];   
      ptmp[i].y =  temp_ptmp[i*3 + 1];   
      ptmp[i].z =  temp_ptmp[i*3 + 2];   
      ptmp[i].node_num =  i + counter;   
    }

  delete temp_ptmp;
 
  unsigned long max_int = ~0L ; 

  for(int i = 0; i < nodes_send; i++)
    {
      result = (ptmp[i].x - min) / diff ;
      int_cod[i].int_x = (unsigned long)(result * max_int);
      
      result = (ptmp[i].y - min) / diff ;
      int_cod[i].int_y = (unsigned long)(result * max_int);
      
      result = (ptmp[i].z - min) / diff ;
      int_cod[i].int_z = (unsigned long)(result * max_int);
	
      int_cod[i].node_num = i + counter;
    }


  if(ac>2 && !strcmp("random",av[2]))
    {
     
      for(int i=0;i<nodes_send;i++)
	Rand_key(counter+i,&new2old[i]);
    }
  else if(ac > 2 && !strcmp("morton",av[2]))  
    {
      // Performs the bit interleaving and o/p in new array
     
   
      for(int i=0;i<nodes_send;i++)
	Converter(&int_cod[i],i,&new2old[i], counter + i);
    }
  else 
    {
      // Performs the bit interleaving and o/p in new array
     
      Point pt;
      Hcode ht;
      for(int i = 0; i < nodes_send; i++)
	{
	  pt.hcode[0] = int_cod[i].int_x;      
	  pt.hcode[1] = int_cod[i].int_y; 
	  pt.hcode[2] = int_cod[i].int_z;     
      
	  ht = H_encode(pt); 
      
	  new2old[i].key3 = ht.hcode[0]  ;      
	  new2old[i].key2 = ht.hcode[1]  ;      
	  new2old[i].key1 = ht.hcode[2]  ;      
	  new2old[i].key_num = counter + i;
	}
    }

  delete int_cod;

  /**********************************************************/
  // PERFORMING THE  SAMPLE SORT ON THE NODES

  
  std::sort(new2old, new2old + nodes_send, compare_key);
 
  // Performing the sample sort on the node numbering according to morton     ordering  
  // first step splitter selection

  int  count = (npes-1)*npes;
 
  length *splitter = new length[npes-1];  
  length *bcast_splitter = new length[npes];    
  
  unsigned long *allpick = new unsigned long[count*4];
 
  for(int i=1;i<npes;i++)
    {       
      splitter[i-1][0]= new2old[i*(nodes_send/npes)].key1;
      splitter[i-1][1]= new2old[i*(nodes_send/npes)].key2;
      splitter[i-1][2]= new2old[i*(nodes_send/npes)].key3;
      splitter[i-1][3]= new2old[i*(nodes_send/npes)].key_num;    
    }
 

  count = (npes-1)*4;
  MPI_Gather( &splitter[0][0], count, MPI_UNSIGNED_LONG, &allpick[0], count, MPI_UNSIGNED_LONG,0, MPI_COMM_WORLD );
 

  count= npes*4;
  MPI_Recv(&bcast_splitter[0],count, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &stat);

 

  //Compute the number of elements that belong to each bucket  
  int *scounts = new int[npes];

  for(int i=0;i<npes;i++)
    scounts[i]=0;
    
  
  for(int i=0;i < nodes_send; i++)
    {
      for(int j=0; j< npes; j++)
	{   
	  if(j==0)
	    {
	      if(compare_splitter(&new2old[i],bcast_splitter[j])) 
		{
		  scounts[j]++;
		  break;  
		}
	    } 
	  else 
	    {
	      if(compare_splitter(&new2old[i],bcast_splitter[j]) && Gcompare(&new2old[i], bcast_splitter[j-1]))
		{
		  scounts[j]++;
		  break; 
		}
	    }
	}
    }

  
  int  y=0;
  // Multiplying by 4 because there are x,y,z coordinates and the key number
  
  for(int i=0;i<npes;i++)
    scounts[i]= scounts[i]*4;
    
  // Determine the starting location of each buckets element 

  int *sdisp = new int[npes];
  sdisp[0] = 0;

  for(int i=1; i<npes; i++)
    sdisp[i]= sdisp[i-1]+scounts[i-1];
    

  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive

  int *rcounts = new int[npes];
  MPI_Alltoall(&scounts[0], 1, MPI_INT, &rcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
 
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
  
  int *rdisp = new int[npes];
  rdisp[0] = 0;

  for(int i=1; i<npes; i++)
    rdisp[i]= rdisp[i-1]+rcounts[i-1];
    
  int total = (rdisp[npes-1]+rcounts[npes-1])/4;
  int t=0;
  
  unsigned long *elemn = new unsigned long[nodes_send*4];
  unsigned long *sort_elemn = new unsigned long[total*4];
 
  for(int i=0;i<total*4;i++)
    sort_elemn[i]=0;
    
 
  for(int i=0;i<nodes_send;i++)
    {
     
      elemn[t++]= new2old[i].key1;  
      elemn[t++]= new2old[i].key2;  
      elemn[t++]= new2old[i].key3;  
      elemn[t++]= new2old[i].key_num;  
    }

 
  // Each PE Sends & recieve element using MPI_Alltoallv
 
  MPI_Alltoallv(&elemn[0], &scounts[0], &sdisp[0], MPI_UNSIGNED_LONG, &sort_elemn[0], &rcounts[0], &rdisp[0], MPI_UNSIGNED_LONG, MPI_COMM_WORLD );

 
  Key *PE_key = new Key[total];
  
  for(int i=0;i<total;i++)
    {   
      PE_key[i].key1 = sort_elemn[i*4 + 0];
      PE_key[i].key2 = sort_elemn[i*4 + 1];
      PE_key[i].key3 = sort_elemn[i*4 + 2];
      PE_key[i].key_num = sort_elemn[i*4 + 3];
    }
  
 
  // LOCAL SORT OF THE KEY 
  std::sort( PE_key, PE_key  + total , compare_key);
  
  
  delete scounts;
  delete elemn;
  delete sort_elemn;
  delete bcast_splitter;
  delete splitter; 
  
  
  unsigned long *new_nodes = new  unsigned long[total];  
  
  for(int i=0;i<total;i++)
    new_nodes[i] = PE_key[i].key_num;
    
  int *recvbuff = new int[npes];
  MPI_Allgather(&total,1, MPI_INT, &recvbuff[0], 1, MPI_INT, MPI_COMM_WORLD );

  int r=0;
  for(int i=0;i<npes;i++)
    r = r + recvbuff[i];
 
  int *rdispla = new int[npes];
  rdispla[0] = 0;
   
  for(int i=1;i<npes;i++)
    rdispla[i]= rdispla[i-1] + recvbuff[i-1];
  
  int sum = rdispla[npes-1] + recvbuff[npes-1];  
  unsigned long *new_node_num = new  unsigned long[sum];
 
  // GATHERS THE NEW NODE NUMBERING IN THE new_node_num array
  MPI_Allgatherv(&new_nodes[0],total, MPI_UNSIGNED_LONG, &new_node_num[0], &recvbuff[0], &rdispla[0], MPI_UNSIGNED_LONG, MPI_COMM_WORLD );
  
  int *old2new = new int[sum];

  for(int i=0; i < sum ;i++)
    old2new[(int)new_node_num[i]] = i; 
  
  // Reordering the xyz coordinates index acc to old2ew numbering 
  
  for(int i=0;i<nodes_send;i++)
    ptmp[i].node_num = (double)old2new[(int) ptmp[i].node_num];

  std::sort(ptmp  , ptmp + nodes_send, compare_node_num);
 
  /**********************************************************/
  // PERFORMING THE  SAMPLE SORT ON THE NODES BECAUSE OF THE NEW NODE NUMBERING 
 
   
  int *split = new int[npes];
  int *scount = new int[npes];
  
  for(int i=0;i<npes;i++)    
    scount[i] = 0 ;
    

  
  for(int i=0;i<npes-1;i++)
    split[i] = Info[0]/npes * (i+1) ;
    
  split[npes-1] = Info[0] ;



  for(int i=0;i<nodes_send;i++)
    {
      for(int j=0; j< npes; j++)
	{   
	  if(j==0)
	    {
	      if(ptmp[i].node_num < split[j])
		{
		  scount[j]++;  
		}
	    } 
	  else 
	    {
	      if( (ptmp[i].node_num < split[j]) && ((ptmp[i].node_num > split[j-1]) ||(ptmp[i].node_num == split[j-1]))) 
		{
		  scount[j]++; 
		}
	    }
	}
    }


  for(int i=0;i<npes;i++)
    scount[i]= scount[i]*4;
  
  
  // Determine the starting location of each buckets element 
  
  int *sdispl = new int[npes];  
  sdispl[0] = 0;
  
  for(int i=1; i<npes; i++)   
    sdispl[i]= sdispl[i-1]+scount[i-1];

  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive

  int *rcount = new int[npes];
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);
 
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
  
  int *rdispl = new int[npes];
  rdispl[0] = 0;

  for(int i=1; i<npes; i++)
    {
      rdispl[i]= rdispl[i-1]+rcount[i-1];
    }

  total = (rdispl[npes-1]+rcount[npes-1])/4;
  t=0;

 
  Node *new_ptmp = new Node[total];

  double *seq = new double[nodes_send*4];
  double *temp_seq = new double[total*4];
    
  for(int i=0;i<nodes_send;i++)
    {     
      seq[i*4 + 0] = ptmp[i].x;  
      seq[i*4 + 1] = ptmp[i].y;  
      seq[i*4 + 2] = ptmp[i].z;  
      seq[i*4 + 3] = ptmp[i].node_num;  
    }
   
  // Each PE Sends & recieve element using MPI_Alltoallv
   
  MPI_Alltoallv(&seq[0], &scount[0], &sdispl[0], MPI_DOUBLE, &temp_seq[0], &rcount[0], &rdispl[0], MPI_DOUBLE, MPI_COMM_WORLD );
   
  for(int i = 0; i < total; i++)
    {
     
      new_ptmp[i].x = temp_seq[i*4 + 0];  
      new_ptmp[i].y = temp_seq[i*4 + 1];  
      new_ptmp[i].z = temp_seq[i*4 + 2];  
      new_ptmp[i].node_num = temp_seq[i*4 + 3];  
      
    }
      
  
  std::sort(new_ptmp  , new_ptmp + total, compare_node_num);
   
  // This variable is used later to send the new_ptmp elements to the PE 0 where the coordinates are written back to the xdr file 
 
  int num_coord = total;
   
  delete temp_seq;
  delete seq;
  delete ptmp;
      
     
  // recvbuff[] gives the information abt the number of faces each PE
  // is going to receive and also abt the number of cellnum pair and offset they are going to receive
  // cellnum_off contains the offset and the cellnum of a face
  int size_array;  
 
  if(myrank==npes-1)
    size_array = nfaces_send*3+1;
  else
    size_array = nfaces_send*3;     

  int *cellnum_off = new int[size_array];
  int offset;
  
  // Receiving cellnum from the PE 0   

  MPI_Recv(&cellnum_off[0],size_array , MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

  // Sending and Receiving the first value(offset) of the cellnum to compute the size of the face by subtracting it from the offset
 
  if(myrank%2==1)
    {     
      MPI_Send(&cellnum_off[0],1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD ); 
      if(myrank!=npes-1)
	MPI_Recv(&offset , 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD, &stat);
    }
  else
    {
      if(myrank!=npes-1)
	MPI_Recv(&offset , 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD, &stat);       
      if(myrank!=0)
	MPI_Send(&cellnum_off[0],1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD );
    }
   
 
  Cell *size_cellnum = new Cell[nfaces_send];

  // Tells how many nodes in a face and the cellnum
  for(int i = 0; i < nfaces_send ; i++)
    {
      size_cellnum[i].cell_left = cellnum_off[i*3 + 1];
      size_cellnum[i].cell_right = cellnum_off[i*3 + 2];	
	       
      if( i != nfaces_send - 1 )
	size_cellnum[i].nodes = cellnum_off[(i+1)*3] - cellnum_off[i*3] ;
      else
	{
	  if(myrank != npes-1 )
	    size_cellnum[i].nodes = offset - cellnum_off[i*3] ;
	  else
	    size_cellnum[i].nodes = cellnum_off[(i+1)*3] - cellnum_off[i*3] ;
	}
    } 
 
  // tells abt the total number of nodes in all the faces on a pe
  int num_facenodes; 

  if(myrank!=npes-1)
    num_facenodes = offset - cellnum_off[0];
  else
    num_facenodes = cellnum_off[nfaces_send * 3] - cellnum_off[0];

  // Sending the faces that is to be received 
  MPI_Send(&num_facenodes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
   
  int *temp_nodes = new int[num_facenodes];

  MPI_Recv(&temp_nodes[0], num_facenodes , MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
 

  // changing the node number to match the new node numbering
  for(int i=0; i<num_facenodes; i++)
    temp_nodes[i]= old2new[temp_nodes[i]];
     
  // size_cellnum[i][0] contains the number of nodes in a face
  // Contains the faces acc to the new node numbering 
   
  int **new_faces;
  new_faces = new int*[nfaces_send];
   
  for(int i = 0; i < nfaces_send ;i++)
    new_faces[i] = new int[size_cellnum[i].nodes];
   
  t=0;

  // Preparing the new faces 
  for(int i=0;i<nfaces_send ;i++)
    {
      for(int j=0;j<size_cellnum[i].nodes;j++)
	{
	  new_faces[i][j]=temp_nodes[t++];
	}
    }

  int minimum; 
 
  int **min_faces;
  min_faces = new int*[nfaces_send];

  for(int i = 0; i <nfaces_send ; i++)
    min_faces[i]= new int[2];
  
     
  /****************************************************************/
  //    Reordering the faces in the morton order

  
  if(ac > 3 && !strcmp("boundary", av[3]))
    {  
      //finding the minimum of the face and storing it into the min_faces
      for(int i=0;i<nfaces_send ;i++)
	{
	  minimum = new_faces[i][0];
	  for(int j = 1; j < size_cellnum[i].nodes ; j++)
	    {
	      if( new_faces[i][j] < minimum )
		minimum = new_faces[i][j]; 	
	    }
	  
	  int tagval = size_cellnum[i].cell_left;
	  if(tagval >  size_cellnum[i].cell_right)
	    tagval =  size_cellnum[i].cell_right ;
	  if(tagval < 0)
	    minimum = tagval * Info[0] + minimum ;
	  
	  min_faces[i][0]=i;
	  min_faces[i][1]=minimum;	  
	}
    }
  else
    {
      //finding the minimum of the face and storing it into the min_faces
      for(int i=0;i<nfaces_send ;i++)
	{
	  minimum = new_faces[i][0];
	  for(int j = 1; j < size_cellnum[i].nodes ; j++)
	    {
	      if( new_faces[i][j] < minimum )
		minimum = new_faces[i][j]; 	
	    }
	  
	  min_faces[i][0]=i;
	  min_faces[i][1]=minimum;	  
	}
    }


 
  // Sorting the faces acc to the minimum node number
  
  std::sort(min_faces, min_faces + nfaces_send , min_compare);
  
  // Allocating the memory to the changed faces 
  int **changed_faces;
  changed_faces = new int*[nfaces_send];
  
  for(int i=0;i<nfaces_send;i++)
    changed_faces[i] = new int[size_cellnum[min_faces[i][0]].nodes];
    
  
  // Aligns faces in morton ordering
  for(int i=0;i < nfaces_send; i++)
    {
      for(int j=0; j<size_cellnum[min_faces[i][0]].nodes;j++)
        {
          changed_faces[i][j] = new_faces[min_faces[i][0]][j];
        }
    }
 
  count = npes*(npes-1);
  
  
  int *splitter1 = new int[npes-1];  
  int *bcast_splitter1 = new int[npes];   
  
  int *allpick1 = new int[count];
  
  
  for(int i=1;i<npes;i++)
    splitter1[i-1]= min_faces[i*(nfaces_send/npes)][1];
    
    
  count = npes-1;
  
  MPI_Gather( &splitter1[0], count, MPI_INT, &allpick1[0], count, MPI_INT,0, MPI_COMM_WORLD );
  MPI_Bcast(&bcast_splitter1[0], npes, MPI_INT, 0, MPI_COMM_WORLD);


  //Compute the number of elements that belong to each bucket  
  int *scount1 = new int[npes];
  
  for(int i=0;i<npes;i++)
    scount1[i]=0;
       
  y=0;
  
  for(int i=0;i<nfaces_send;i++)
    {
      for(int j=0; j< npes; j++)
	{   
	  if(j==0)
	    {
	      if(min_faces[i][1]<bcast_splitter1[j]) 
		{
		  scount1[j]++;  
		}
	    } 
	  else 
	    {
	      if( (min_faces[i][1]<bcast_splitter1[j]) &&( min_faces[i][1]> bcast_splitter1[j-1] || min_faces[i][1] == bcast_splitter1[j-1]))
		{
		  scount1[j]++; 
		}
	    }
	}
    }
   

   
  // Sending the sorted min_faces to other PE's so that they can perform local sort on basis of that
  // Determine the starting location of each buckets element 
   
  int *sdisp_min =  new int[npes];   
  sdisp_min[0] = 0;
   
  for(int i=1; i<npes; i++)
    sdisp_min[i]= sdisp_min[i-1]+scount1[i-1];
     
   
  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive
   
  int *rcount_min = new int[npes];
   
  MPI_Alltoall(&scount1[0], 1, MPI_INT, &rcount_min[0], 1, MPI_INT, MPI_COMM_WORLD);
   
 
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
   
  int *rdisp_min = new int[npes];
  rdisp_min[0] = 0;
   
  for(int i=1; i<npes; i++) 
    rdisp_min[i]= rdisp_min[i-1]+rcount_min[i-1];
      
   
  total = 0; 
  total = (rdisp_min[npes-1]+rcount_min[npes-1]);
   
  int *pick_minfaces = new int[total];
  int *temp_min = new int[nfaces_send];
   
  for(int i=0; i < nfaces_send; i++)
    {
      temp_min[i] = min_faces[i][1];
    }
   
  
  // Each PE Sends & recieve element using MPI_Alltoallv   
  MPI_Alltoallv(&temp_min[0], &scount1[0], &sdisp_min[0], MPI_INT, &pick_minfaces[0], &rcount_min[0], &rdisp_min[0], MPI_INT, MPI_COMM_WORLD );
   
  // gving an index to the collected minfaces
   
  int **gather_minfaces; 
  gather_minfaces = new int*[total];
  
  for(int i=0;i<total;i++)
    gather_minfaces[i] = new int[2];  
   
   
  for(int i=0;i<total;i++)
    {
      gather_minfaces[i][0]=i;  
      gather_minfaces[i][1]=pick_minfaces[i];  
    }
   
 
  delete temp_min;
  delete pick_minfaces; 
  

  //          LOCAL SORT
  // Sorting the faces acc to the minimum node number  
  std::sort(gather_minfaces,gather_minfaces + total , min_compare);
   
   
  //APPLYING THE BUCKET SORT ON THE FACES                              
                           
  t=0;
  total=0;
  int *total1 = new int[npes];
     
  for(int i=0;i<npes;i++)
    {
      total1[i]=0;
    }
     
  for(int i=0;i<npes;i++)
    {
      for(int j=0;j<scount1[i];j++)
	{ 
	  total1[i] = total1[i] + size_cellnum[min_faces[t++][0]].nodes;
	}
    }
     
     
     
  // Determine the starting location of each buckets element 

  int *sdisp1 = new int[npes];
  sdisp1[0] = 0;

  for(int i=1; i<npes; i++)
    {
      sdisp1[i]= sdisp1[i-1]+total1[i-1];
    }

  // Perform an all-all to inform the corresponding processes
  // the number of elements they are going to receive

  int *rcount1 = new int[npes];

  MPI_Alltoall(&total1[0], 1, MPI_INT, &rcount1[0], 1, MPI_INT, MPI_COMM_WORLD);

 
  //Based on rcounts determine where in the local array the data
  //from each PE will be stored
 
  int *rdisp1 = new int[npes];
  rdisp1[0] = 0;

  for(int i=1; i<npes; i++)
    rdisp1[i]= rdisp1[i-1]+rcount1[i-1];
     

  total = (rdisp1[npes-1]+rcount1[npes-1]);
  
  int *elemn1 = new int[num_facenodes];
  int *sort_elemn1 = new int[total];
 
  // This variable is used later to send the faces to the first PE to write the faces back to the xdr file;

  int num_finalfaces = total;
  t=0;

  for(int i=0;i<nfaces_send;i++)
    {
      for(int j = 0;j < size_cellnum[min_faces[i][0]].nodes ;j++)
	elemn1[t++]= changed_faces[i][j] ;   
    }
    
  
  // Each PE Sends & recieve element using MPI_Alltoallv
     
  MPI_Alltoallv(&elemn1[0], &total1[0], &sdisp1[0], MPI_INT, &sort_elemn1[0], &rcount1[0], &rdisp1[0], MPI_INT, MPI_COMM_WORLD);

  delete elemn1;

  for(int i=0;i<npes;i++)     
    scount1[i]= scount1[i]*3;     
     
 
  // For sending the face size to the PE's since faces can have different sizes, rcount2 tells the number of faces received from each PE
 
  int *rcount2 = new int[npes];
     
  MPI_Alltoall(&scount1[0], 1, MPI_INT, &rcount2[0], 1, MPI_INT, MPI_COMM_WORLD);


  // size array has the size of the faces only, which is needed by othere PE's 
  // Preparing cellnum and size array to send to other PE 

  int *size = new int[nfaces_send*3];
  t=0;
  
  for(int i=0;i<nfaces_send;i++)
    {
      size[t++]= size_cellnum[min_faces[i][0]].nodes;
      size[t++]= size_cellnum[min_faces[i][0]].cell_left;
      size[t++]= size_cellnum[min_faces[i][0]].cell_right;
    }
  
  
  sdisp1[0] = 0;
  for(int i=1; i<npes; i++)
    sdisp1[i]= 0;
    
   
  for(int i=1; i<npes; i++)
    sdisp1[i]= sdisp1[i-1] + scount1[i-1];
     
   
  rdisp1[0] = 0;
  for(int i=1; i<npes; i++)
    rdisp1[i]= 0;
       
   
  for(int i=1; i<npes; i++)
    rdisp1[i]= rdisp1[i-1]+rcount2[i-1];
         
    
  total = (rdisp1[npes-1]+rcount2[npes-1]);     
  int *recv_size = new int[total];
    
  // Each PE Sends & recieve element using MPI_Alltoallv   
  MPI_Alltoallv(&size[0], &scount1[0], &sdisp1[0], MPI_INT, &recv_size[0], &rcount2[0], &rdisp1[0], MPI_INT, MPI_COMM_WORLD );

  // Allocating memory TO THE NEW FACES 
  int **faces;
  faces = new int*[total/3];
   
  for(int i=0;i<total/3;i++)
    faces[i] = new int[recv_size[i*3]];
     
  t=0;
   
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[i*3] ;j++)
	{
	  faces[i][j] = sort_elemn1[t++];
	}
    }


  // Allocating memory to faces acc to the local sort of the minfaces
 
  int **final_faces;
  final_faces = new int*[total/3];
 
  for(int i=0;i<total/3;i++)
    final_faces[i] = new int[recv_size[(gather_minfaces[i][0])*3]];
     
 
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[(gather_minfaces[i][0])*3] ;j++)	
	final_faces[i][j] = faces[(gather_minfaces[i][0])][j];	 
    }
    
  t=0;
  for(int i=0;i<total/3;i++)
    {
      for(int j=0; j<recv_size[(gather_minfaces[i][0])*3] ;j++)
	sort_elemn1[t++]= final_faces[i][j];
    }
    
  delete final_faces;
  delete faces;
   
  // Sending the first offset to the next processor
  sum=0;
   
  for(int i=0;i<total/3;i++)
    sum = sum + recv_size[i*3];
   
   
  // sum_all[] is used to give the offset to the pe's for cellnum
   
  int *sum_all = new int[npes];
  MPI_Allgather(&sum, 1, MPI_INT, &sum_all[0], 1, MPI_INT, MPI_COMM_WORLD);
   
  offset=0;
   
  for(int i=0;i<myrank;i++)     
    offset = offset + sum_all[i];
     
   
  // Giving new offset to the cellnum acc to the Morton ordering
  // off is used because last PE has one more element i.e last offset than others 
   
  int off = 0;
   
  if(myrank == npes-1)
    off=1;
  else
    off=0;
   
  int *new_array = new int[total + off];
    
  //This variable is used later while sending the cellnum to the PE 0 for writing it back to the xdr file

  int num_cellnum = total + off;

  t=0;
  // Contains the cellnumber in the morton ordering with the new offset
   
  for(int i = 0;i < total/3; i++)
    {
      new_array[t++] =  offset ;     
      new_array[t++] =  recv_size[(gather_minfaces[i][0])*3+1] ;
      new_array[t++] =  recv_size[(gather_minfaces[i][0])*3+2] ;
      offset = offset + recv_size[(gather_minfaces[i][0])*3] ;
    }
   
  if(myrank == npes-1)
    new_array[t]=offset;
   
  int ncells = Info[2];
  int ctr=0;
  t=0;
   
 
  // ctr gets value from myrank-1 PE which is used as counter in renumbering the cells 
  // according to first visit by cellnum
  
  int *new_cellnum = new int[ncells];
   
  MPI_Recv(&ctr , 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, &stat);   
 
  MPI_Recv(&new_cellnum[0] ,ncells, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, &stat);

  // Cells reordered acc to their first visit by cellnum
   
  for(int i=0; i < total/3; i++)
    {    
      if(new_array[i*3+1]>0)
        {
          if(new_cellnum[new_array[i*3+1]-1]==-1)
            new_cellnum[new_array[i*3+1]-1]=ctr++;
        }  

      if(new_array[i*3+2]>0)
        {
          if(new_cellnum[new_array[i*3+2]-1]==-1)
            new_cellnum[new_array[i*3+2]-1]=ctr++;
        }  
    }

  // Sending the ctr and new_cellnum array to next PE to give a new number to the cell
    
  if(myrank != npes-1)
    {
      MPI_Send(&ctr, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD ); 
      MPI_Send(&new_cellnum[0], ncells, MPI_INT, myrank+1, 0, MPI_COMM_WORLD );  
    }

  // new_cellnum contains the new cellnumber for each cell

  MPI_Barrier( MPI_COMM_WORLD);
  MPI_Bcast( &new_cellnum[0], ncells, MPI_INT, npes-1, MPI_COMM_WORLD);  

  // Re-ordering the cells acc to their first visit by cellnum

  for(int i=0; i < total/3; i++)
    {    
      if(new_array[i*3+1]>0)
	new_array[i*3+1] = new_cellnum[new_array[i*3+1]-1];
       
      if(new_array[i*3+2]>0)	
	new_array[i*3+2] = new_cellnum[new_array[i*3+2]-1];    
    }
   
   
  delete new_cellnum;
  
  // Preparing the nodes to send to the source processor

  num_coord = num_coord * 3;
  double *send_coord = new double[num_coord];
  
  for(int i=0;i< num_coord/3;i++)
    {       
      send_coord[i*3+0] = new_ptmp[i].x;  
      send_coord[i*3+1] = new_ptmp[i].y;  
      send_coord[i*3+2] = new_ptmp[i].z;  
    }

  tag1=1;
  
  /***************** sending the nodes information *****************/
  
  MPI_Send(&num_coord, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );  
  MPI_Send(&tag1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );  
  MPI_Send(&send_coord[0], num_coord, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );  
  
  /************** sending the cellnum information *******************/
  
  MPI_Send(&num_cellnum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );  
  MPI_Send(&tag1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );  
  MPI_Send(&new_array[0], num_cellnum, MPI_INT, 0, 0, MPI_COMM_WORLD ); 

 
  /************** sending the faces information *******************/
  
  MPI_Send(&num_finalfaces, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ); 
  MPI_Send(&tag1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );  
  MPI_Send(&sort_elemn1[0], num_finalfaces, MPI_INT, 0, 0, MPI_COMM_WORLD );  
 

}



// Function for decimal to binary conversion

void dec2bin(unsigned long t , int b[])
{

  
  unsigned long a = t;  
  
  for(int i=31;i>=0;i--)
    { 
      if(a>0 )
	{
	  b[i] = a % 2;
	  a = a>>1;
	}
      else
	b[i] = 0;
    }
   
}
  
// Function for binary to decimal conversion

unsigned long bin2dec(int c[])
{
  unsigned long num=0;
  unsigned long  n = 2;
  int t=31;
  
  for(int i=0;i<32;i++)
    {
      if(i==0)
	num = num+c[t--];
      else  
	num = num+(unsigned long)((n<<(i-1))*c[t--]);
    } 
  
  return num;
}


// Bit interleaving routine

void Converter(Node_int *int_cod,int x, Key *new2old,int counter)
{
  
  int b1[32],b2[32],b3[32];
  int c1[32],c2[32],c3[32];    

  dec2bin(int_cod->int_x, b1);
  dec2bin(int_cod->int_y, b2);
  dec2bin(int_cod->int_z, b3);

  // Arranging the bits so that the bits are interleaved

  for(int i=0;i<10;i++)
    {
      c1[i*3]=b1[i];
      c1[i*3+1]=b2[i];
      c1[i*3+2]=b3[i];
    }
  
  c1[30]=b1[10];
  c1[31]=b2[10];
  
  c2[0]=b3[10];

  int t=1;
  
  for(int i=11;i<21;i++)
    {
      c2[t++] = b1[i];
      c2[t++] = b2[i];
      c2[t++] = b3[i];
    }

  t=0;  

  for(int i=21;i<32;i++)
    {
      if(i==21)
	{
	  c2[31] = b1[i];
	  c3[t++] = b2[i];
	  c3[t++] = b3[i];
	}
      else
	{
	  c3[t++] = b1[i];
	  c3[t++] = b2[i];
	  c3[t++] = b3[i];
	}
    }

   
  new2old->key1= bin2dec(c1);
  new2old->key2= bin2dec(c2);
  new2old->key3= bin2dec(c3);
  new2old->key_num= counter;
   

}




void Rand_key(int n, Key *new2old )
{
 
  srand48((unsigned int)time((time_t *)NULL));
   
  new2old->key1 = (unsigned long)lrand48();
  new2old->key2 = (unsigned long)lrand48();
  new2old->key3 = (unsigned long)lrand48();
  new2old->key_num = n;
   
}   


/*============================================================================*/
/*                            H_encode					      */
/*============================================================================*/
/*
 * given the coordinates of a point, it finds the sequence number of the point
 * on the Hilbert Curve
 */
Hcode H_encode(Point p)
{
  U_long	mask = (U_long)1 << WORDBITS - 1, element, temp1, temp2,
    A, W = 0, S, tS, T, tT, J, P = 0, xJ;

  Hcode	h;
  h.hcode[0] = 0;
  h.hcode[1] = 0;
  h.hcode[2] = 0;
	
  int	i = NUMBITS * DIM - DIM, j;
	
  for (j = A = 0; j < DIM; j++)
    if (p.hcode[j] & mask)
      A |= g_mask[j];
	
  S = tS = A;
	
  P |= S & g_mask[0];
  for (j = 1; j < DIM; j++)
    if( S & g_mask[j] ^ (P >> 1) & g_mask[j])
      P |= g_mask[j];
	
  /* add in DIM bits to hcode */
  element = i / WORDBITS;
  if (i % WORDBITS > WORDBITS - DIM)
    {
      h.hcode[element] |= P << i % WORDBITS;
      h.hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
    }
  else
    h.hcode[element] |= P << i - element * WORDBITS;

  J = DIM;
  for (j = 1; j < DIM; j++)
    if ((P >> j & 1) == (P & 1))
      continue;
    else
      break;
  if (j != DIM)
    J -= j;
  xJ = J - 1;
	
  if (P < 3)
    T = 0;
  else
    if (P % 2)
      T = (P - 1) ^ (P - 1) / 2;
    else
      T = (P - 2) ^ (P - 2) / 2;
  tT = T;

  for (i -= DIM, mask >>= 1; i >=0; i -= DIM, mask >>= 1)
    {
      for (j = A = 0; j < DIM; j++)
        if (p.hcode[j] & mask)
          A |= g_mask[j];

      W ^= tT;
      tS = A ^ W;
      if (xJ % DIM != 0)
        {
          temp1 = tS << xJ % DIM;
          temp2 = tS >> DIM - xJ % DIM;
          S = temp1 | temp2;
          S &= ((U_long)1 << DIM) - 1;
        }
      else
        S = tS;

      P = S & g_mask[0];
      for (j = 1; j < DIM; j++)
        if( S & g_mask[j] ^ (P >> 1) & g_mask[j])
          P |= g_mask[j];

      /* add in DIM bits to hcode */
      element = i / WORDBITS;
      if (i % WORDBITS > WORDBITS - DIM)
        {
          h.hcode[element] |= P << i % WORDBITS;
          h.hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
        }
      else
        h.hcode[element] |= P << i - element * WORDBITS;

      if (i > 0)
        {
          if (P < 3)
            T = 0;
          else
            if (P % 2)
              T = (P - 1) ^ (P - 1) / 2;
            else
              T = (P - 2) ^ (P - 2) / 2;

          if (xJ % DIM != 0)
            {
              temp1 = T >> xJ % DIM;
              temp2 = T << DIM - xJ % DIM;
              tT = temp1 | temp2;
              tT &= ((U_long)1 << DIM) - 1;
            }
          else
            tT = T;

          J = DIM;
          for (j = 1; j < DIM; j++)
            if ((P >> j & 1) == (P & 1))
              continue;
            else
              break;
          if (j != DIM)
            J -= j;

          xJ += J - 1;
          /*	J %= DIM;*/
        }
    }
  return h;
}


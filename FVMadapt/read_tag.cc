//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
//******************************************************************************************
// this file read in the refinement plan 
//moified 03/12/07, read in cell tags before face tags
//******************************************************************************************
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include <Tools/tools.h>
#include "sciTypes.h"
#include "defines.h"
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using namespace Loci;
int currentMem(void);

//functions from distribute_io
namespace Loci{
  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm );
  vector<int> simplePartitionVec(int mn, int mx, int p) {
    vector<int> nums(p+1) ;
    int n = mx-mn+1 ;
    int dn = n/p ; // divisor
    int rn = n%p ; // remainder
    int start = mn ;
    nums[0] = start ;
    for(int i=0;i<p;++i) {
      start += dn+((i<rn)?1:0) ;
      nums[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return nums ;
  }
  vector<sequence> transposeSeq(const vector<sequence> sv) {
     vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = sv[i].num_intervals()*2 ;
    vector<int> recv_sz(MPI_processes) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }
    //    outRep->allocate(new_alloc) ;
    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<sv[i].num_intervals();++j) {
        send_store[send_displacement[i]+j*2] = sv[i][j].first ;
        send_store[send_displacement[i]+j*2+1] = sv[i][j].second ;
      }
    
    
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  

    vector<sequence> sv_t(MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        sv_t[i] += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return sv_t ;
  }
  
// result allocated over local numbering 
// input is in file numbering
// Routine redistributes input to be arranged in the local numbering
// and placed in result.
// input is assumed to be partitioned by simplePartion
  void distribute_reorder_store(Loci::storeRepP &result, entitySet resultSet,
                                Loci::storeRepP input, fact_db &facts ){
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    
    dMap g2f ;
    g2f = df->g2f.Rep() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    
    Map newnum ;
    newnum.allocate(resultSet) ;
    
    FORALL(resultSet,i) {
      newnum[i] = g2f[l2g[i]] ;
    } ENDFORALL ;

    const int p = MPI_processes ;
    // Get input (file) distribution
    int imx = GLOBAL_MAX(input->domain().Max()) ;
    int imn = GLOBAL_MIN(input->domain().Min()) ;
    vector<int> fptn = simplePartitionVec(imn,imx,p) ;

    // Get distribution plan
    vector<vector<pair<int,int> > > dist_plan(p) ;
    int np = (imx-imn+1)/p ; // number of elements per processor
    FORALL(resultSet,i) {
      int fn = newnum[i] ;
      if(fn < imn || fn > imx) {
        cerr << "Problem with distribute_reorder_store, index out of bounds"
             << endl ;
        cerr << "fn = " << fn << "imx = " << imx << "imn = " << imn << endl ;
        Loci::Abort() ;
      }
      // processor that contains this value
      int r = min((fn-imn)/np,p-1) ; // Guess which processor
      for(;;) { // search from guess
        if(fn >= fptn[r] && fn < fptn[r+1])
          break ;
        r+= (fn < fptn[r])?-1:1 ;
        if(r < 0 || r >= p) {
          debugout << "bad r in processor search " << r << endl ;
        }
        FATAL(r >= p) ;
        FATAL(r < 0) ;
      }
      dist_plan[r].push_back(pair<int,int>(fn,i)) ;
    } ENDFORALL ;
    // Compute recv requests from distribution plan
    vector<sequence> recv_seq(p),send_req(p) ;
    for(int i=0;i<p;++i) {
      std::sort(dist_plan[i].begin(),dist_plan[i].end()) ;
      sequence s1,s2 ;
      int psz = dist_plan[i].size() ;
      for(int j=0;j<psz;++j) {
        s1 +=dist_plan[i][j].first ;
        s2 +=dist_plan[i][j].second ;
      }
      send_req[i] = s1 ;
      recv_seq[i] = s2 ;
    }

    // Transpose the send requests to get the sending sequences
    // from this processor
    vector<sequence> send_seq = transposeSeq(send_req) ;
    vector<entitySet> send_sets(p) ;
    for(int i=0;i<p;++i)
      send_sets[i] = entitySet(send_seq[i]) ;
    
    int *send_sizes = new int[p] ;
    int *recv_sizes = new int[p] ;


    for(int i=0;i<p;++i)
      send_sizes[i] = input->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int *send_dspl = new int[p] ;
    int *recv_dspl = new int[p] ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    unsigned char *send_store = new unsigned char[send_sz] ;
    unsigned char *recv_store = new unsigned char[recv_sz] ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      input->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(send_store, &send_sizes[0], send_dspl, MPI_PACKED,
		  recv_store, &recv_sizes[0], recv_dspl, MPI_PACKED,
		  MPI_COMM_WORLD) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      result->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                     recv_seq[i]) ;
    }
    delete[] recv_store ;
    delete[] send_store ;
    delete[] recv_dspl ;
    delete[] send_dspl ;
    delete[] recv_sizes ;
    delete[] send_sizes ;
  
  }
}

class get_postag : public pointwise_rule{
  const_param<std::string> tagfile_par;
  store<char> posTag;
  
public:
  get_postag(){
    name_store("tagfile_par", tagfile_par);
    name_store("posTag", posTag);
    input("tagfile_par");
    output("posTag");
    constraint("pos");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    ifstream inFile;
    int nprocs = Loci::MPI_processes;
    
    //process 0 open tagFile
    if(Loci::MPI_rank == 0){
      inFile.open((*tagfile_par).c_str());
      if(!inFile){
        cerr <<"can not open " << *tagfile_par << " for input" << endl;
        Loci::Abort();
      }
      //  cout << "Reading tagfile " << *tagfile_par << " for posTag" << endl;
    }

    entitySet dom = entitySet(seq);
  
    //serial version
    if(nprocs == 1){
      char tag; 
      FORALL(dom, cc){
        if(inFile >> tag) posTag[cc] = char(tag - '0');
        else posTag[cc] = 0;
      
      }ENDFORALL;
      //close tagfile
      inFile.close();
      //      cout << "Finish reading posTag" << endl;
      return;
    }
    //parallel version
    
    fact_db::distribute_infoP d = (Loci::exec_current_fact_db)->get_distribute_info() ;
        
    //map local node dom to global node dom,
    // and check if the size of  two dom is the same 
    Map l2g ;
    l2g = d->l2g.Rep() ;
    MapRepP l2gP = MapRepP(l2g.Rep()) ;
    entitySet dom_global = l2gP->image(dom) ;
    FATAL(dom.size() != dom_global.size()) ;
    
    //get map from global to file
    dMap g2f ;
    g2f = d->g2f.Rep() ;
    
    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[l2g[i]] ;
    } ENDFORALL ;
    
    //local min and max of file numbering of node dom
    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;
    
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;
    
    //find the min and max value over all processes
    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;
  
    // simple patition the whole file numbering node domain 
    vector<entitySet> out_ptn = Loci::simplePartition(imn,imx,MPI_COMM_WORLD) ;
    //end of code from collect_reorder_store()  
    
    //temp_posTag is allocated on the simple partition of file numbering of nodes
    store<char> temp_posTag;
    temp_posTag.allocate(out_ptn[Loci::MPI_rank]);
    
    //buf size, i.e., maximum  num of local nodes of simple partition
    int buf_size = 0;
    for(int i = 0; i < nprocs; i++){
      buf_size = max(buf_size, out_ptn[i].size());
    }
    
    char* buf = new char[buf_size] ;
        
    if(Loci::MPI_rank == 0){
      char tag;
      //process 0 read in its local tag 
      FORALL(out_ptn[0], cc){
        if(inFile >> tag) temp_posTag[cc] = char(tag - '0');
        else temp_posTag[cc] = 0; 
      }ENDFORALL;
      
      //process 0 read the tags of other processes into a buf and then send the buf
      for(int i = 1; i < nprocs; i++){
        int send_size = out_ptn[i].size();
        for(int j = 0; j < send_size; j++)  
          {
            if(inFile >> tag) buf[j] = char(tag - '0');
            else buf[j] = 0; 
          }
        MPI_Send(buf, send_size, MPI_CHAR, i, 20, MPI_COMM_WORLD);
      }          
    }else{ //other processes recv the buf and unpack it
   
    
      MPI_Status status;
      MPI_Recv(buf, buf_size, MPI_CHAR, 0, 20, MPI_COMM_WORLD, &status);
      int ptr = 0;
      FORALL(out_ptn[Loci::MPI_rank], cc){
        temp_posTag[cc] = buf[ptr++];
      }ENDFORALL;
    }
    //finish read in temp_posTag;

        
    //redistribute temp_posTag to local node dom
    
    Loci::storeRepP posTagRepP = posTag.getRep();
    Loci::distribute_reorder_store(posTagRepP, dom, temp_posTag.getRep(),
                                   *(Loci::exec_current_fact_db));
            
    delete [] buf;
    //process 0 close file
    if(Loci::MPI_rank == 0){
      inFile.close();
      //  cout << "Finish reading  posTag " << endl;
    }
    
    //   cerr << currentMem() << " after reading posTag " << Loci::MPI_rank << endl; 
  }
};
register_rule<get_postag> register_get_postag;

class get_nodetag : public pointwise_rule{
  const_param<std::string> tagfile_par;
  const_param<int> num_original_nodes;
  const_store<int> num_inner_nodes;
  store<std::vector<char> > nodeTag;
public:
  get_nodetag(){
    name_store("tagfile_par", tagfile_par);
    name_store("num_original_nodes", num_original_nodes);
    name_store("nodeTag", nodeTag);
    name_store("num_inner_nodes", num_inner_nodes);
    input("tagfile_par");
    input("num_original_nodes");
    input("num_inner_nodes");
    output("nodeTag");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
 
    
    ifstream inFile;
    int nprocs = Loci::MPI_processes;
    //process 0 open tagFile
    if(Loci::MPI_rank == 0){
      inFile.open((*tagfile_par).c_str());
      if(!inFile){
        cerr <<"can not open " << *tagfile_par << " for input" << endl;
        Loci::Abort();
      }
   
    }
    //find the length of the file
    
   
    
    //serial version
    if(nprocs == 1){
      char tag;
      for(int i = 0; i < *num_original_nodes; i++) inFile>>tag;
      
      Loci::constraint edges, geom_cells, faces;
      Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
      *edges = e2n->domain();
      faces = (Loci::exec_current_fact_db)->get_variable("faces");
      geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");
        
      FORALL(*edges, cc){
        if(num_inner_nodes[cc] != 0){
          std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
          for(int i = 0; i < num_inner_nodes[cc]; i++){
            if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
            else nodeTag[cc][i] = 0;
          }
        }else{
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();
        }
        
      }ENDFORALL; 
      
      
      
      FORALL(*geom_cells, cc){
        if(num_inner_nodes[cc] != 0){
          std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
          for(int i = 0; i < num_inner_nodes[cc]; i++){
            if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
            else nodeTag[cc][i] = 0;
          }
        }else{
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();
        }
      }ENDFORALL;
        
      FORALL(*faces, cc){
        if(num_inner_nodes[cc] != 0){
          std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
          for(int i = 0; i < num_inner_nodes[cc]; i++){
            if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
            else nodeTag[cc][i] = 0;
          }
        }else{
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();
        }
      }ENDFORALL;
      
      
      //close tagfile
      inFile.close();
      
      return;
    }
    //parallel version
    
    
    fact_db::distribute_infoP d = (Loci::exec_current_fact_db)->get_distribute_info() ;
    Loci::constraint my_entities ; 
    my_entities = d->my_entities ;
    
    
    //start working on nodeTag
    
    Loci::constraint edges, geom_cells, faces;
    Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
    *edges = e2n->domain(); 
    faces = (Loci::exec_current_fact_db)->get_variable("faces");
    geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");
    
    entitySet local_edges, local_faces, local_geom_cells;
    //don't know if it's necessray
    local_edges = (*my_entities) & (*edges) ;
    local_faces = (*my_entities) & (*faces);
    local_geom_cells = (*my_entities)&(*geom_cells);
    
    

    
    
    //compute num of local inner nodes on each process
    int num_local_inner_nodes = 0;
    
    FORALL(local_edges, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    //write out cell_nodes first, then write face_nodes
    FORALL(local_geom_cells, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    
    FORALL(local_faces, ei){
      num_local_inner_nodes += num_inner_nodes[ei];
    }ENDFORALL;
    
    
    //compute buf_size on each process
    unsigned int  buf_size = 0;
    MPI_Allreduce(&num_local_inner_nodes, &buf_size, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    
    if(buf_size == 0){

      FORALL(local_edges, cc){
        nodeTag[cc].clear();
      }ENDFORALL;
      FORALL(local_geom_cells, cc){
        nodeTag[cc].clear();
      }ENDFORALL;
      
      FORALL(local_faces, cc){
        nodeTag[cc].clear();
      }ENDFORALL;
      
      return;
    }

      
      //process 0 find out size of buffer for each process
      int* size_buf = new int[nprocs];
      MPI_Gather(&num_local_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      char *nbuf = new char[buf_size]; //nodeTag buf
      
      if(Loci::MPI_rank == 0){
        char tag;
        
        for(int i = 0; i < *num_original_nodes; i++) inFile>>tag;
        //process 0 read in its local nodeTag
        FORALL(local_edges, cc){
          if(num_inner_nodes[cc] != 0){
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
              else nodeTag[cc][i] = 0;
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
          
        }ENDFORALL; 
        
    
        
        FORALL(local_geom_cells, cc){
          if(num_inner_nodes[cc] != 0){
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
              else nodeTag[cc][i] = 0;
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
        }ENDFORALL;
        
        FORALL(local_faces, cc){
          if(num_inner_nodes[cc] != 0){
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
              else nodeTag[cc][i] = 0;
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
        }ENDFORALL;
        
        //process 0 read in nodeTags into nbuf and send nbuf to other processes
        for(int i = 1; i < nprocs; i++){
          
          for(int j = 0; j < size_buf[i]; j++)  
            {
              if(inFile >> tag) nbuf[j] = char(tag - '0');
              else nbuf[j] = 0; 
            }
          if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_CHAR, i, 11, MPI_COMM_WORLD);
        }
      }else{ // all other process recv nbuf from process 0 and unpack nodeTag
        int recv_size = num_local_inner_nodes;
        if(recv_size != 0){
          MPI_Status status;
          MPI_Recv(&nbuf[0], buf_size, MPI_CHAR, 0, 11, MPI_COMM_WORLD, &status);
        }  
        int ptr = 0;
        FORALL(local_edges, cc){
          if(num_inner_nodes[cc] != 0){
            
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              nodeTag[cc][i] = nbuf[ptr++];
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
        }ENDFORALL; 
        
        FORALL(local_geom_cells, cc){
          if(num_inner_nodes[cc]!= 0){
            
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              nodeTag[cc][i] = nbuf[ptr++];
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
        }ENDFORALL;
        
        FORALL(local_faces, cc){
          if(num_inner_nodes[cc] != 0){
            std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
            for(int i = 0; i < num_inner_nodes[cc]; i++){
              nodeTag[cc][i] = nbuf[ptr++];
            }
          }else{
             std::vector<char>(1).swap(nodeTag[cc]);
            nodeTag[cc].clear();
          }
        }ENDFORALL;
        
      }
        
  
         
   
      //clean up
       
      delete [] nbuf;
      delete [] size_buf;
    

     //process 0 close file
      if(Loci::MPI_rank == 0){
        inFile.close();
        // cout << "Finish reading  nodeTag " << endl;
      }
      //   cerr << currentMem() << " After reading in nodeTag " << Loci::MPI_rank << endl; 
  }
  
};

register_rule<get_nodetag> register_get_nodetag;

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
#include <cstdlib>
#include <string>
#include <Loci.h>
#include <vector>
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


// result allocated over local numbering 
// input is in file numbering
// Routine redistributes input to be arranged in the local numbering
// and placed in result.
// input is assumed to be partitioned by simplePartion
//void distribute_reorder_store(Loci::storeRepP &result, entitySet resultSet,
//                              Loci::storeRepP input, fact_db &facts );
  

  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm);

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
    Loci::File2LocalOrder(posTagRepP, dom, temp_posTag.getRep(), 0,
                          Loci::exec_current_fact_db->get_distribute_info(),
                          MPI_COMM_WORLD) ;
    //    Loci::distribute_reorder_store(posTagRepP, dom, temp_posTag.getRep(),
    //                                   *(Loci::exec_current_fact_db));
            
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

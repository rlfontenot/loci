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
#include <hdf5.h>

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


  std::vector<int> all_collect_sizes(int size);
// Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);
  
  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm);
  
}


// class get_postag : public pointwise_rule{
//   const_param<std::string> tagfile_par;
//   store<char> posTag;
  
// public:
//   get_postag(){
//     name_store("tagfile_par", tagfile_par);
//     name_store("posTag", posTag);
//     input("tagfile_par");
//     output("posTag");
//     constraint("pos");
//     disable_threading();
//   }
//   virtual void compute(const sequence &seq){
//     ifstream inFile;
//     int nprocs = Loci::MPI_processes;
    
//     //process 0 open tagFile
//     if(Loci::MPI_rank == 0){
//       inFile.open((*tagfile_par).c_str());
//       if(!inFile){
//         cerr <<"can not open " << *tagfile_par << " for input" << endl;
//         Loci::Abort();
//       }
//     }
    
//     entitySet dom = entitySet(seq);
    
//     //serial version
//     if(nprocs == 1){
//       char tag; 
//       FORALL(dom, cc){
//         if(inFile >> tag) posTag[cc] = char(tag - '0');
//         else posTag[cc] = 0;
      
//       }ENDFORALL;
//       //close tagfile
//       inFile.close();
//       return;
//     }
//     //parallel version
    
//     fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
//     //first create store in the order of file numbering      
//     int offset = 0;
//     store<char> temp_posTag;
//     temp_posTag = Local2FileOrder(posTag.Rep(),dom,offset,dist,MPI_COMM_WORLD) ;;
//     int num_nodes = temp_posTag.domain().size();
    
//     int* size_buf = new int[nprocs];
    
//     //compute buf_size on each process
//     unsigned int  buf_size = 0;
//     MPI_Allreduce(&num_nodes, &buf_size, 1, MPI_INT,
//                   MPI_MAX, MPI_COMM_WORLD);
//     //process 0 find out size of buffer for each process
     
//     MPI_Gather(&num_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
//     char *buf = new char[buf_size]; //posTag buf
    
//     if(Loci::MPI_rank == 0){
//       char tag;
//       //process 0 read in its local tag 
//       FORALL(temp_posTag.domain(), cc){
//         if(inFile >> tag) temp_posTag[cc] = char(tag - '0');
//         else temp_posTag[cc] = 0; 
//       }ENDFORALL;
      
//       //process 0 read the tags of other processes into a buf and then send the buf
//       for(int i = 1; i < nprocs; i++){
//         int send_size = size_buf[i];
//         for(int j = 0; j < send_size; j++)  
//           {
//             if(inFile >> tag) buf[j] = char(tag - '0');
//             else buf[j] = 0; 
//           }
//         MPI_Send(buf, send_size, MPI_CHAR, i, 20, MPI_COMM_WORLD);
//       }          
//     }else{ //other processes recv the buf and unpack it
   
    
//       MPI_Status status;
//       MPI_Recv(buf, buf_size, MPI_CHAR, 0, 20, MPI_COMM_WORLD, &status);
//       int ptr = 0;
//       FORALL(temp_posTag.domain(), cc){
//         temp_posTag[cc] = buf[ptr++];
//       }ENDFORALL;
//     }
    
//     //finish read in temp_posTag;

        
//     //redistribute temp_posTag to local node dom
//     storeRepP localVar=posTag.Rep();
//     File2LocalOrder(localVar, dom,
//                     temp_posTag.Rep(), offset,
//                     dist,
//                     MPI_COMM_WORLD);
  
               
//     delete [] buf;
//     delete [] size_buf;
//     //process 0 close file
//     if(Loci::MPI_rank == 0){
//       inFile.close();
//       //  cout << "Finish reading  posTag " << endl;
//     }
    
    
//   }
// };
//   register_rule<get_postag> register_get_postag;
  
// class get_nodetag : public pointwise_rule{
//   const_param<std::string> tagfile_par;
//   const_param<int> num_original_nodes;
//   const_store<int> num_inner_nodes;
//   store<std::vector<char> > nodeTag;
// public:
//   get_nodetag(){
//     name_store("tagfile_par", tagfile_par);
//     name_store("num_original_nodes", num_original_nodes);
//     name_store("nodeTag", nodeTag);
//     name_store("num_inner_nodes", num_inner_nodes);
//     input("tagfile_par");
//     input("num_original_nodes");
//     input("num_inner_nodes");
//     output("nodeTag");
//     disable_threading();
//   }
//   virtual void compute(const sequence &seq){
 
    
//     ifstream inFile;
//     int nprocs = Loci::MPI_processes;
//     //process 0 open tagFile
//     if(Loci::MPI_rank == 0){
//       inFile.open((*tagfile_par).c_str());
//       if(!inFile){
//         cerr <<"can not open " << *tagfile_par << " for input" << endl;
//         Loci::Abort();
//       }
   
//     }
//     //find the length of the file
    
   
    
//     //serial version
//     if(nprocs == 1){
//       char tag;
//       for(int i = 0; i < *num_original_nodes; i++) inFile>>tag;
      
//       Loci::constraint edges, geom_cells, faces;
//       Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
//       *edges = e2n->domain();
//       faces = (Loci::exec_current_fact_db)->get_variable("faces");
//       geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");
        
//       FORALL(*edges, cc){
//         if(num_inner_nodes[cc] != 0){
//           std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
//           for(int i = 0; i < num_inner_nodes[cc]; i++){
//             if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
//             else nodeTag[cc][i] = 0;
//           }
//         }else{
//           std::vector<char>(1).swap(nodeTag[cc]);//to avoid default allocated size for vector
//           nodeTag[cc].clear();
//         }
        
//       }ENDFORALL; 
      
      
      
//       FORALL(*geom_cells, cc){
//         if(num_inner_nodes[cc] != 0){
//           std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
//           for(int i = 0; i < num_inner_nodes[cc]; i++){
//             if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
//             else nodeTag[cc][i] = 0;
//           }
//         }else{
//           std::vector<char>(1).swap(nodeTag[cc]);
//           nodeTag[cc].clear();
//         }
//       }ENDFORALL;
        
//       FORALL(*faces, cc){
//         if(num_inner_nodes[cc] != 0){
//           std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
//           for(int i = 0; i < num_inner_nodes[cc]; i++){
//             if(inFile >> tag) nodeTag[cc][i] = char(tag - '0');
//             else nodeTag[cc][i] = 0;
//           }
//         }else{
//           std::vector<char>(1).swap(nodeTag[cc]);
//           nodeTag[cc].clear();
//         }
//       }ENDFORALL;
      
      
//       //close tagfile
//       inFile.close();
      
//       return;
//     }
//     //parallel version
    
    
//     fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
//     Loci::constraint my_entities ; 
//     my_entities = dist->my_entities ;
    
    
//     //start working on nodeTag
    
//     Loci::constraint edges, geom_cells, faces;
//     Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
//     *edges = e2n->domain(); 
//     faces = (Loci::exec_current_fact_db)->get_variable("faces");
//     geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");
    
//     entitySet local_edges, local_faces, local_geom_cells;
//     //don't know if it's necessray
//     local_edges = (*my_entities) & (*edges) ;
//     local_faces = (*my_entities) & (*faces);
//     local_geom_cells = (*my_entities)&(*geom_cells);


//     int offset = 0;
//     unsigned int buf_size = 0;
//     int* size_buf = new int[MPI_processes];
//     //put in a block to release memory when out of block
//     {
//       // file numbering, the domain of this container shifted by offset
//       // is the actual file numbering.
//       store<int> edge_num_inner_nodes;
//       edge_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
//       //for each processor, 
//       int num_edge_inner_nodes = 0;
//       FORALL(edge_num_inner_nodes.domain(), ei){
//         num_edge_inner_nodes += edge_num_inner_nodes[ei];
//       }ENDFORALL;
   
//       buf_size = 0;
//       MPI_Allreduce(&num_edge_inner_nodes, &buf_size, 1, MPI_INT,
//                     MPI_MAX, MPI_COMM_WORLD);
      
//       if(buf_size == 0){
//         FORALL(local_edges, cc){
//           std::vector<char>(1).swap(nodeTag[cc]);
//           nodeTag[cc].clear();  
//         }ENDFORALL;
        
//       }else{
//         //process 0 find out size of buffer for each process
        
        
//         MPI_Gather(&num_edge_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
//         char *nbuf = new char[buf_size]; //nodeTag buf
        
//         //nodeTag in the order of file numbering
//         store<std::vector<char> > edge_nodeTag;
//         edge_nodeTag.allocate(edge_num_inner_nodes.domain());
        
//         if(Loci::MPI_rank == 0){
//           char tag;
//           for(int i = 0; i < *num_original_nodes; i++) inFile>>tag;
//           //process 0 read in its local nodeTag
//           FORALL(edge_num_inner_nodes.domain(), cc){
//             if(edge_num_inner_nodes[cc] != 0){
//               std::vector<char>(int(edge_num_inner_nodes[cc])).swap(edge_nodeTag[cc]);
//               for(int i = 0; i < edge_num_inner_nodes[cc]; i++){
//                 if(inFile >> tag) edge_nodeTag[cc][i] = char(tag - '0');
//                 else edge_nodeTag[cc][i] = 0;
//               }
//             }else{
//               std::vector<char>(1).swap(edge_nodeTag[cc]);
//               edge_nodeTag[cc].clear();
//             }
            
//           }ENDFORALL; 
          
//           //process 0 read in nodeTags into nbuf and send nbuf to other processes
//           for(int i = 1; i < nprocs; i++){
            
//             for(int j = 0; j < size_buf[i]; j++)  
//               {
//                 if(inFile >> tag) nbuf[j] = char(tag - '0');
//                 else nbuf[j] = 0; 
//               }
//             if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_CHAR, i, 11, MPI_COMM_WORLD);
//           }
//         }else{ // all other process recv nbuf from process 0 and unpack nodeTag
//           int recv_size = num_edge_inner_nodes;
//           if(recv_size != 0){
//             MPI_Status status;
//             MPI_Recv(&nbuf[0], buf_size, MPI_CHAR, 0, 11, MPI_COMM_WORLD, &status);
//           }  
//           int ptr = 0;
//           FORALL(edge_num_inner_nodes.domain(), cc){
//             if(edge_num_inner_nodes[cc] != 0){
              
//               std::vector<char>(int(edge_num_inner_nodes[cc])).swap(edge_nodeTag[cc]);
//               for(int i = 0; i < edge_num_inner_nodes[cc]; i++){
//                 edge_nodeTag[cc][i] = nbuf[ptr++];
//               }
//             }else{
//               std::vector<char>(1).swap(edge_nodeTag[cc]);
//               edge_nodeTag[cc].clear();
//             }
//           }ENDFORALL; 
//         }
        
//       storeRepP localVar=nodeTag.Rep();
//       File2LocalOrder(localVar, local_edges,
//                       edge_nodeTag.Rep(), offset,
//                       dist,
//                       MPI_COMM_WORLD);
      
//       delete [] nbuf;
//       }
//     }//finish edge nodeTag
//     //start with cell nodeTag
//      offset = 0;
//      {
//        // file numbering, the domain of this container shifted by offset
//        // is the actual file numbering.
//        store<int> cell_num_inner_nodes;
//        cell_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
//        //for each processor, 
//        int num_cell_inner_nodes = 0;
//        FORALL(cell_num_inner_nodes.domain(), ei){
//          num_cell_inner_nodes += cell_num_inner_nodes[ei];
//        }ENDFORALL;
       
//        buf_size = 0;
//        MPI_Allreduce(&num_cell_inner_nodes, &buf_size, 1, MPI_INT,
//                      MPI_MAX, MPI_COMM_WORLD);
       
//        if(buf_size == 0){
//          FORALL(local_geom_cells, cc){
//            std::vector<char>(1).swap(nodeTag[cc]);
//         nodeTag[cc].clear();  
//          }ENDFORALL;
         
//        }else{
//          //process 0 find out size of buffer for each process

         
//          MPI_Gather(&num_cell_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
         
//          char *nbuf = new char[buf_size]; //nodeTag buf
        
//          //nodeTag in the order of file numbering
//          store<std::vector<char> > cell_nodeTag;
//          cell_nodeTag.allocate(cell_num_inner_nodes.domain());
         
//       if(Loci::MPI_rank == 0){
//         char tag;
        
//        //process 0 read in its local nodeTag
//         FORALL(cell_num_inner_nodes.domain(), cc){
//           if(cell_num_inner_nodes[cc] != 0){
//             std::vector<char>(int(cell_num_inner_nodes[cc])).swap(cell_nodeTag[cc]);
//             for(int i = 0; i < cell_num_inner_nodes[cc]; i++){
//               if(inFile >> tag) cell_nodeTag[cc][i] = char(tag - '0');
//               else cell_nodeTag[cc][i] = 0;
//             }
//           }else{
//             std::vector<char>(1).swap(cell_nodeTag[cc]);
//             cell_nodeTag[cc].clear();
//           }
          
//         }ENDFORALL; 
        
//        //process 0 read in nodeTags into nbuf and send nbuf to other processes
//         for(int i = 1; i < nprocs; i++){
          
//           for(int j = 0; j < size_buf[i]; j++)  
//            {
//              if(inFile >> tag) nbuf[j] = char(tag - '0');
//              else nbuf[j] = 0; 
//            }
//          if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_CHAR, i, 12, MPI_COMM_WORLD);
//         }
//       }else{ // all other process recv nbuf from process 0 and unpack nodeTag
//         int recv_size = num_cell_inner_nodes;
//         if(recv_size != 0){
//             MPI_Status status;
//             MPI_Recv(&nbuf[0], buf_size, MPI_CHAR, 0, 12, MPI_COMM_WORLD, &status);
//         }  
//           int ptr = 0;
//           FORALL(cell_num_inner_nodes.domain(), cc){
//             if(cell_num_inner_nodes[cc] != 0){
              
//               std::vector<char>(int(cell_num_inner_nodes[cc])).swap(cell_nodeTag[cc]);
//               for(int i = 0; i < cell_num_inner_nodes[cc]; i++){
//                 cell_nodeTag[cc][i] = nbuf[ptr++];
//               }
//             }else{
//               std::vector<char>(1).swap(cell_nodeTag[cc]);
//               cell_nodeTag[cc].clear();
//             }
//           }ENDFORALL; 
//       }
      
//       storeRepP localVar=nodeTag.Rep();
//       File2LocalOrder(localVar, local_geom_cells,
//                       cell_nodeTag.Rep(), offset,
//                       dist,
//                       MPI_COMM_WORLD);
     
//       delete [] nbuf;
//        }
//      }//finsh cell nodeTag
     
//        //start for faces
//      offset= 0;
//      {
//        // Create container vardist that is ordered across processors in the
//        // file numbering, the domain of this container shifted by offset
//        // is the actual file numbering.
//        store<int> face_num_inner_nodes;
//        face_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_faces,offset,dist,MPI_COMM_WORLD) ;
       
       
//        int num_face_inner_nodes = 0;
//        FORALL(face_num_inner_nodes.domain(), ei){
//          num_face_inner_nodes += face_num_inner_nodes[ei];
//        }ENDFORALL;
       
//        buf_size = 0;
//        MPI_Allreduce(&num_face_inner_nodes, &buf_size, 1, MPI_INT,
//                      MPI_MAX, MPI_COMM_WORLD);
       
//        if(buf_size == 0){
//          FORALL(local_faces, cc){
//            std::vector<char>(1).swap(nodeTag[cc]);
//            nodeTag[cc].clear();  
//          }ENDFORALL;
         
//        }else{
         
         
//          //process 0 find out size of buffer for each process
      
//          MPI_Gather(&num_face_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
         
//          char* nbuf = new char[buf_size]; //nodeTag buf
         
//          store<std::vector<char> > face_nodeTag;
//          face_nodeTag.allocate(face_num_inner_nodes.domain());
         
//          if(Loci::MPI_rank == 0){
//            char tag;
           
//            //process 0 read in its local nodeTag
//            FORALL(face_num_inner_nodes.domain(), cc){
//              if(face_num_inner_nodes[cc] != 0){
//             std::vector<char>(int(face_num_inner_nodes[cc])).swap(face_nodeTag[cc]);
//             for(int i = 0; i < face_num_inner_nodes[cc]; i++){
//               if(inFile >> tag) face_nodeTag[cc][i] = char(tag - '0');
//               else face_nodeTag[cc][i] = 0;
//             }
//              }else{
//                std::vector<char>(1).swap(face_nodeTag[cc]);
//                face_nodeTag[cc].clear();
//              }
             
//            }ENDFORALL; 
           
//            //process 0 read in nodeTags into nbuf and send nbuf to other processes
//            for(int i = 1; i < nprocs; i++){
             
//              for(int j = 0; j < size_buf[i]; j++)  
//                {
//               if(inFile >> tag) nbuf[j] = char(tag - '0');
//               else nbuf[j] = 0; 
//             }
//              if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_CHAR, i, 13, MPI_COMM_WORLD);
//            }
//          }else{ // all other process recv nbuf from process 0 and unpack nodeTag
//            int recv_size = num_face_inner_nodes;
//            if(recv_size != 0){
//           MPI_Status status;
//           MPI_Recv(&nbuf[0], buf_size, MPI_CHAR, 0, 13, MPI_COMM_WORLD, &status);
//            }  
//            int ptr = 0;
//            FORALL(face_num_inner_nodes.domain(), cc){
//           if(face_num_inner_nodes[cc] != 0){
            
//             std::vector<char>(int(face_num_inner_nodes[cc])).swap(face_nodeTag[cc]);
//             for(int i = 0; i < face_num_inner_nodes[cc]; i++){
//               face_nodeTag[cc][i] = nbuf[ptr++];
//             }
//           }else{
//             std::vector<char>(1).swap(face_nodeTag[cc]);
//             face_nodeTag[cc].clear();
//           }
//            }ENDFORALL; 
//          }
         
//          storeRepP localVar=nodeTag.Rep();
//          File2LocalOrder(localVar, local_faces,
//                          face_nodeTag.Rep(), offset,
//                          dist,
//                          MPI_COMM_WORLD);
         
         
         
         
//          //clean up
         
//          delete [] nbuf;
//        }
//      }
  
//      delete [] size_buf;
     
//      //process 0 close file
//      if(Loci::MPI_rank == 0){
//        inFile.close();
       
//      }
     
//   }
  
// };

//   register_rule<get_nodetag> register_get_nodetag;
  

// class get_postag : public pointwise_rule{
//   const_param<std::string> tagfile_par;
//   store<char> posTag;
  
// public:
//   get_postag(){
//     name_store("tagfile_par", tagfile_par);
//     name_store("posTag", posTag);
//     input("tagfile_par");
//     output("posTag");
//     constraint("pos");
//     disable_threading();
//   }
//   virtual void compute(const sequence &seq){
//     hid_t file_id;
//     entitySet dom = entitySet(seq);
//     file_id = H5Fopen((*tagfile_par).c_str(), H5F_ACC_RDONLY,
//                       H5P_DEFAULT);
       
//     Loci::readContainer(file_id,"ref",posTag.Rep(),dom) ;
//     H5Fclose(file_id);


//     // string filename = "output/" + sname + "_sca." + iter + "_" + casename ;
    
  
//   }
//  };
// register_rule<get_postag> register_get_postag;




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
   
    int nprocs = Loci::MPI_processes;
    entitySet dom = entitySet(seq);
       
    //serial version
    if(nprocs == 1){
      hid_t file_id;
      file_id = H5Fopen((*tagfile_par).c_str(), H5F_ACC_RDONLY,
                        H5P_DEFAULT);
      hid_t group_id = 0;
      group_id = H5Gopen(file_id, "ref") ;
      //process 0 read in its local tag
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int rank = 1 ;
      hsize_t dimension = dom.size() ;
      count = dimension ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      float *buf = new float[dom.size()]; 
      H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
               H5P_DEFAULT, buf);
      H5Sclose(memspace) ;
      
      int ptr = 0;
      FORALL(dom, cc){
        posTag[cc] = buf[ptr++]==1?1:0;
      }ENDFORALL;  
      
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Gclose(group_id) ;
      H5Fclose(file_id);
      delete [] buf;
      return;
    }
   
    
    
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    //first create store in the order of file numbering      
    int offset = 0;
    store<char> temp_posTag;
    temp_posTag = Local2FileOrder(posTag.Rep(),dom,offset,dist,MPI_COMM_WORLD) ;;
    int num_nodes = temp_posTag.domain().size();
    
    int* size_buf = new int[nprocs];
    
    //compute buf_size on each process
    unsigned int  buf_size = 0;
    MPI_Allreduce(&num_nodes, &buf_size, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    //process 0 find out size of buffer for each process
     
    MPI_Gather(&num_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
    float *buf = new float[buf_size]; //posTag buf
    
    if(Loci::MPI_rank == 0){
      hid_t file_id;
      file_id = H5Fopen((*tagfile_par).c_str(), H5F_ACC_RDONLY,
                        H5P_DEFAULT);
      hid_t group_id = 0;
      group_id = H5Gopen(file_id, "ref") ;
       
      //process 0 read in its local tag
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int rank = 1 ;
      hsize_t dimension = size_buf[0] ;
      count = dimension ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;

      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
               H5P_DEFAULT, buf);
      H5Sclose(memspace) ;
      start += count;
      int ptr = 0;
      FORALL(temp_posTag.domain(), cc){
        temp_posTag[cc] = (buf[ptr++]==1?1:0);
      }ENDFORALL;
      
      //process 0 read the tags of other processes into a buf and then send the buf
      for(int i = 1; i < nprocs; i++){
        int send_size = size_buf[i];
        hsize_t dimension = size_buf[i] ;
        count = dimension ;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
        H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                 H5P_DEFAULT, buf);
        H5Sclose(memspace) ;
        start += count;
        MPI_Send(buf, send_size, MPI_FLOAT, i, 20, MPI_COMM_WORLD);
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Gclose(group_id) ;
      H5Fclose(file_id);
       
    }else{ //other processes recv the buf and unpack it
   
    
      MPI_Status status;
      MPI_Recv(buf, buf_size, MPI_FLOAT, 0, 20, MPI_COMM_WORLD, &status);
      int ptr = 0;
      FORALL(temp_posTag.domain(), cc){
        temp_posTag[cc] = (buf[ptr++]==1?1:0);
      }ENDFORALL;
    }
    
    //finish read in temp_posTag;

        
    //redistribute temp_posTag to local node dom
    storeRepP localVar=posTag.Rep();
    File2LocalOrder(localVar, dom,
                    temp_posTag.Rep(), offset,
                    dist,
                    MPI_COMM_WORLD);
  
               
    delete [] buf;
    delete [] size_buf;
    //process 0 close file
    //  cout << "Finish reading  posTag " << endl;
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
     int nprocs = Loci::MPI_processes;
    //process 0 open tagFile
     
    
    //serial version
    if(nprocs == 1){
      hid_t file_id;
      file_id = H5Fopen((*tagfile_par).c_str(), H5F_ACC_RDONLY,
                        H5P_DEFAULT);
      hid_t group_id = 0;
      group_id = H5Gopen(file_id, "ref") ;
      //process 0 read in its local tag
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int rank = 1 ;
      
      start += *num_original_nodes; //skip the posTag part
      hsize_t dimension = 0;
      entitySet dom = entitySet(seq);
      FORALL(dom, cc){
        dimension += num_inner_nodes[cc];
      }ENDFORALL; 
      count = dimension ;
     
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      float *buf = new float[count]; 
      H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
               H5P_DEFAULT, buf);
      H5Sclose(memspace) ;
      
            
      Loci::constraint edges, geom_cells, faces;
      Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
      *edges = e2n->domain();
      faces = (Loci::exec_current_fact_db)->get_variable("faces");
      geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");   
      
      int ptr = 0;
      FORALL(*edges, cc){
        if(num_inner_nodes[cc] != 0){
          std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
          for(int i = 0; i < num_inner_nodes[cc]; i++){
            nodeTag[cc][i] = (buf[ptr++]==1?1:0);
          }
        }else{
          std::vector<char>(1).swap(nodeTag[cc]);//to avoid default allocated size for vector
          nodeTag[cc].clear();
        }
        
      }ENDFORALL; 
      
      
      
      FORALL(*geom_cells, cc){
        if(num_inner_nodes[cc] != 0){
          std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
          for(int i = 0; i < num_inner_nodes[cc]; i++){
            nodeTag[cc][i] = (buf[ptr++]==1?1:0);
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
            nodeTag[cc][i] = (buf[ptr++]==1?1:0);
          }
        }else{
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();
        }
      }ENDFORALL;
      
      
      //close tagfile
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Gclose(group_id) ;
      H5Fclose(file_id);
      delete [] buf; 
      return;
    }
    //parallel version
    
    
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    Loci::constraint my_entities ; 
    my_entities = dist->my_entities ;
    
    
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

#ifdef H5_INTERFACE_1_6_4
    hsize_t start = 0 ;
#else
    hssize_t start = 0 ;
#endif
    hid_t file_id;
    hid_t group_id = 0;
    hid_t dataset;
    hid_t dataspace;
    int rank = 1 ;
    hsize_t stride = 1 ;
    
    int offset = 0;
    unsigned int buf_size = 0;
    int* size_buf = new int[MPI_processes];
    //put in a block to release memory when out of block
    {
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      store<int> edge_num_inner_nodes;
      edge_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
      //for each processor, 
      int num_edge_inner_nodes = 0;
      FORALL(edge_num_inner_nodes.domain(), ei){
        num_edge_inner_nodes += edge_num_inner_nodes[ei];
      }ENDFORALL;
   
      buf_size = 0;
      MPI_Allreduce(&num_edge_inner_nodes, &buf_size, 1, MPI_INT,
                    MPI_MAX, MPI_COMM_WORLD);
      
      if(buf_size == 0){
        FORALL(local_edges, cc){
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();  
        }ENDFORALL;
        
      }else{
        //process 0 find out size of buffer for each process
        
        
        MPI_Gather(&num_edge_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        float *nbuf = new float[buf_size]; //nodeTag buf
        
        //nodeTag in the order of file numbering
        store<std::vector<char> > edge_nodeTag;
        edge_nodeTag.allocate(edge_num_inner_nodes.domain());
        
        if(Loci::MPI_rank == 0){
         
          file_id = H5Fopen((*tagfile_par).c_str(), H5F_ACC_RDONLY,
                            H5P_DEFAULT);
          
          group_id = H5Gopen(file_id, "ref") ;
          //process 0 read in its local tag
          dataset =  H5Dopen(group_id, "data") ;
          dataspace = H5Dget_space(dataset) ;
          start += *num_original_nodes; //skip the posTag part
         
          hsize_t count = 0 ;
          hsize_t dimension = size_buf[0];
          count = dimension ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
          H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                   H5P_DEFAULT, nbuf);
          H5Sclose(memspace) ;
          start += count;
          
          //process 0 read in its local nodeTag
          int ptr = 0;
          FORALL(edge_num_inner_nodes.domain(), cc){
            if(edge_num_inner_nodes[cc] != 0){
              std::vector<char>(int(edge_num_inner_nodes[cc])).swap(edge_nodeTag[cc]);
              for(int i = 0; i < edge_num_inner_nodes[cc]; i++){
                edge_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              }
            }else{
              std::vector<char>(1).swap(edge_nodeTag[cc]);
              edge_nodeTag[cc].clear();
            }
          }ENDFORALL; 
          
          //process 0 read in nodeTags into nbuf and send nbuf to other processes
          for(int i = 1; i < nprocs; i++){
            dimension = size_buf[i];
            count = dimension ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                     H5P_DEFAULT, nbuf);
            H5Sclose(memspace) ;
             start += count;
            if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_FLOAT, i, 11, MPI_COMM_WORLD);
          }
    
         
        }else{ // all other process recv nbuf from process 0 and unpack nodeTag
          int recv_size = num_edge_inner_nodes;
          if(recv_size != 0){
            MPI_Status status;
            MPI_Recv(&nbuf[0], buf_size, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, &status);
          }  
          int ptr = 0;
          FORALL(edge_num_inner_nodes.domain(), cc){
            if(edge_num_inner_nodes[cc] != 0){
              
              std::vector<char>(int(edge_num_inner_nodes[cc])).swap(edge_nodeTag[cc]);
              for(int i = 0; i < edge_num_inner_nodes[cc]; i++){
                edge_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              }
            }else{
              std::vector<char>(1).swap(edge_nodeTag[cc]);
              edge_nodeTag[cc].clear();
            }
          }ENDFORALL; 
        }
        
        storeRepP localVar=nodeTag.Rep();
        File2LocalOrder(localVar, local_edges,
                        edge_nodeTag.Rep(), offset,
                        dist,
                        MPI_COMM_WORLD);
        
        delete [] nbuf;
      }
    }//finish edge nodeTag

    //start with cell nodeTag
    offset = 0;
    {
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      store<int> cell_num_inner_nodes;
      cell_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
      //for each processor, 
      int num_cell_inner_nodes = 0;
      FORALL(cell_num_inner_nodes.domain(), ei){
        num_cell_inner_nodes += cell_num_inner_nodes[ei];
      }ENDFORALL;
       
      buf_size = 0;
      MPI_Allreduce(&num_cell_inner_nodes, &buf_size, 1, MPI_INT,
                    MPI_MAX, MPI_COMM_WORLD);
       
      if(buf_size == 0){
        FORALL(local_geom_cells, cc){
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();  
        }ENDFORALL;
         
      }else{
        //process 0 find out size of buffer for each process

         
        MPI_Gather(&num_cell_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
         
        float *nbuf = new float[buf_size]; //nodeTag buf
        
        //nodeTag in the order of file numbering
        store<std::vector<char> > cell_nodeTag;
        cell_nodeTag.allocate(cell_num_inner_nodes.domain());
         
        if(Loci::MPI_rank == 0){
         
          hsize_t count = 0 ;
          hsize_t dimension = size_buf[0];
          count = dimension ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
          H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                   H5P_DEFAULT, nbuf);
          H5Sclose(memspace) ;
           start += count; 
          //process 0 read in its local nodeTag
          int ptr = 0;
          FORALL(cell_num_inner_nodes.domain(), cc){
            if(cell_num_inner_nodes[cc] != 0){
              std::vector<char>(int(cell_num_inner_nodes[cc])).swap(cell_nodeTag[cc]);
              for(int i = 0; i < cell_num_inner_nodes[cc]; i++){
                cell_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              }
            }else{
              std::vector<char>(1).swap(cell_nodeTag[cc]);
              cell_nodeTag[cc].clear();
            }
          
          }ENDFORALL; 
        
          //process 0 read in nodeTags into nbuf and send nbuf to other processes
          for(int i = 1; i < nprocs; i++){
            hsize_t count = 0 ;
            hsize_t dimension = size_buf[i];
            count = dimension ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                     H5P_DEFAULT, nbuf);
            H5Sclose(memspace) ;
            start += count;
            if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_FLOAT, i, 12, MPI_COMM_WORLD);
          }
          
        }else{ // all other process recv nbuf from process 0 and unpack nodeTag
          int recv_size = num_cell_inner_nodes;
          if(recv_size != 0){
            MPI_Status status;
            MPI_Recv(&nbuf[0], buf_size, MPI_FLOAT, 0, 12, MPI_COMM_WORLD, &status);
          }  
          int ptr = 0;
          FORALL(cell_num_inner_nodes.domain(), cc){
            if(cell_num_inner_nodes[cc] != 0){
            
              std::vector<char>(int(cell_num_inner_nodes[cc])).swap(cell_nodeTag[cc]);
              for(int i = 0; i < cell_num_inner_nodes[cc]; i++){
                cell_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              }
            }else{
              std::vector<char>(1).swap(cell_nodeTag[cc]);
              cell_nodeTag[cc].clear();
            }
          }ENDFORALL; 
        }
      
        storeRepP localVar=nodeTag.Rep();
        File2LocalOrder(localVar, local_geom_cells,
                        cell_nodeTag.Rep(), offset,
                        dist,
                        MPI_COMM_WORLD);
     
        delete [] nbuf;
      }
    }//finsh cell nodeTag
     
    //start for faces
    offset= 0;
    {
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      store<int> face_num_inner_nodes;
      face_num_inner_nodes = Local2FileOrder(num_inner_nodes.Rep(),local_faces,offset,dist,MPI_COMM_WORLD) ;
       
       
      int num_face_inner_nodes = 0;
      FORALL(face_num_inner_nodes.domain(), ei){
        num_face_inner_nodes += face_num_inner_nodes[ei];
      }ENDFORALL;
       
      buf_size = 0;
      MPI_Allreduce(&num_face_inner_nodes, &buf_size, 1, MPI_INT,
                    MPI_MAX, MPI_COMM_WORLD);
       
      if(buf_size == 0){
        FORALL(local_faces, cc){
          std::vector<char>(1).swap(nodeTag[cc]);
          nodeTag[cc].clear();  
        }ENDFORALL;
         
      }else{
         
         
        //process 0 find out size of buffer for each process
      
        MPI_Gather(&num_face_inner_nodes, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
         
        float* nbuf = new float[buf_size]; //nodeTag buf
         
        store<std::vector<char> > face_nodeTag;
        face_nodeTag.allocate(face_num_inner_nodes.domain());
         
        if(Loci::MPI_rank == 0){
          hsize_t count = 0 ;
          hsize_t dimension = size_buf[0];
          count = dimension ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
          H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                   H5P_DEFAULT, nbuf);
          H5Sclose(memspace) ;
            start += count;
          int ptr = 0;
          //process 0 read in its local nodeTag
          FORALL(face_num_inner_nodes.domain(), cc){
            if(face_num_inner_nodes[cc] != 0){
              std::vector<char>(int(face_num_inner_nodes[cc])).swap(face_nodeTag[cc]);
              for(int i = 0; i < face_num_inner_nodes[cc]; i++){
                face_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              
              }
            }else{
              std::vector<char>(1).swap(face_nodeTag[cc]);
              face_nodeTag[cc].clear();
            }
             
          }ENDFORALL; 
           
          //process 0 read in nodeTags into nbuf and send nbuf to other processes
          for(int i = 1; i < nprocs; i++){
            hsize_t count = 0 ;
            hsize_t dimension = size_buf[i];
            count = dimension ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                     H5P_DEFAULT, nbuf);
            H5Sclose(memspace) ;
             start += count;            
            if(size_buf[i] != 0) MPI_Send(&nbuf[0], size_buf[i], MPI_FLOAT, i, 13, MPI_COMM_WORLD);
          }
        }else{ // all other process recv nbuf from process 0 and unpack nodeTag
          int recv_size = num_face_inner_nodes;
          if(recv_size != 0){
            MPI_Status status;
            MPI_Recv(&nbuf[0], buf_size, MPI_FLOAT, 0, 13, MPI_COMM_WORLD, &status);
          }  
          int ptr = 0;
          FORALL(face_num_inner_nodes.domain(), cc){
            if(face_num_inner_nodes[cc] != 0){
               
              std::vector<char>(int(face_num_inner_nodes[cc])).swap(face_nodeTag[cc]);
              for(int i = 0; i < face_num_inner_nodes[cc]; i++){
                face_nodeTag[cc][i] = (nbuf[ptr++]==1?1:0);
              }
            }else{
              std::vector<char>(1).swap(face_nodeTag[cc]);
              face_nodeTag[cc].clear();
            }
          }ENDFORALL; 
        }
         
        storeRepP localVar=nodeTag.Rep();
        File2LocalOrder(localVar, local_faces,
                        face_nodeTag.Rep(), offset,
                        dist,
                        MPI_COMM_WORLD);
         
         
         
         
        //clean up
         
        delete [] nbuf;
      }
       
    }
  
    delete [] size_buf;
     
    //process 0 close file
    if(Loci::MPI_rank == 0){
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Gclose(group_id) ;
      H5Fclose(file_id);
    }
     
  }
  
};

register_rule<get_nodetag> register_get_nodetag;
  

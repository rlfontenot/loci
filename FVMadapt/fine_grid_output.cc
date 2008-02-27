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
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include <rpc/xdr.h>
#include <rpc/rpc.h>
#include "sciTypes.h"
#include "defines.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;

//typedef Loci::vector3d<double> vect3d;
namespace Loci{
storeRepP collect_reorder_store(storeRepP &sp, fact_db &facts);
std::vector<int> all_collect_sizes(int size);
}



class node_output_file : public pointwise_rule {
 
  store<bool> node_output ;
  const_param<string> outfile_par ;
  const_store<Loci::FineNodes> inner_nodes;
 
  
  
public:
  node_output_file(){
    name_store("node_output", node_output);
    name_store("outfile_par", outfile_par);
    name_store("inner_nodes", inner_nodes);
 
 
    
    input("outfile_par");
    input("inner_nodes");
 
    
    output("node_output");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    
    FILE *out;
    XDR xdr_handle;
    
   
   
    if(num_procs == 1){

      out = fopen((*outfile_par).c_str(),"a");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;
     
      Loci::constraint faces, geom_cells;
      faces = Loci::exec_current_fact_db->get_variable("faces");
      geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
      entitySet local_edges = Loci::exec_current_fact_db->get_variable("edge2node")->domain();
      //out inner nodes
       FORALL(local_edges, cc){
         for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
           vect3d p = inner_nodes[cc][i];
           xdr_double(&xdr_handle, &p.x) ; 
           xdr_double(&xdr_handle, &p.y) ;
           xdr_double(&xdr_handle, &p.z) ;
         }
       }ENDFORALL;
       FORALL(*geom_cells, cc){
         for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
           vect3d p = inner_nodes[cc][i];
           xdr_double(&xdr_handle, &p.x) ; 
           xdr_double(&xdr_handle, &p.y) ;
           xdr_double(&xdr_handle, &p.z) ;
         }
       }ENDFORALL; 
       FORALL(*faces, cc){
         for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
           vect3d p = inner_nodes[cc][i];
           xdr_double(&xdr_handle, &p.x) ; 
           xdr_double(&xdr_handle, &p.y) ;
           xdr_double(&xdr_handle, &p.z) ;
         }
       }ENDFORALL;
      
       fclose (out);
       xdr_destroy(&xdr_handle) ;
       return;
    }
    fact_db::distribute_infoP d =  Loci::exec_current_fact_db->get_distribute_info() ;
    Loci::constraint my_entities, faces, geom_cells, interior_faces, boundary_faces; 
    my_entities = d->my_entities ;
    faces = Loci::exec_current_fact_db->get_variable("faces");
    geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
    interior_faces =  Loci::exec_current_fact_db->get_variable("interior_faces");
    boundary_faces =  Loci::exec_current_fact_db->get_variable("boundary_faces");
    
    entitySet local_edges =*my_entities &
      (Loci::exec_current_fact_db->get_variable("edge2node"))->domain(); 
    entitySet local_faces = *my_entities & *faces;
    entitySet local_cells =  *my_entities & *geom_cells;
    entitySet local_interior_faces = *my_entities & *interior_faces;
    entitySet local_boundary_faces = *my_entities & *boundary_faces;

 //process 0 open the file     
    if(my_id == 0){
      out = fopen((*outfile_par).c_str(),"a");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;
    }
    
   
    //output inner nodes
    //each proc compute its num_local_inner_nodes
    int num_local_inner_nodes = 0;
    FORALL(local_edges, cc){
      num_local_inner_nodes += inner_nodes[cc].size();
    }ENDFORALL;
    FORALL(local_cells, cc){
      num_local_inner_nodes += inner_nodes[cc].size();
    }ENDFORALL; 
    FORALL(local_faces, cc){
      num_local_inner_nodes += inner_nodes[cc].size();
    }ENDFORALL;
    
    std::vector<int> inner_nodes_sizes;
    inner_nodes_sizes = Loci::all_collect_sizes(num_local_inner_nodes);
   
    int inner_nodes_buf_size =  (*std::max_element(inner_nodes_sizes.begin(), inner_nodes_sizes.end()))*3;
    double* inner_nodes_buf = new double[inner_nodes_buf_size];
     cerr << currentMem() << " node_output " <<  Loci::MPI_rank << endl;
    //first process 0 write out its inner_nodes
    if(my_id == 0){
    
      FORALL(local_edges, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          xdr_double(&xdr_handle, &p.x) ; 
          xdr_double(&xdr_handle, &p.y) ;
          xdr_double(&xdr_handle, &p.z) ;
          
        }
      }ENDFORALL;
      FORALL(local_cells, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          xdr_double(&xdr_handle, &p.x) ; 
          xdr_double(&xdr_handle, &p.y) ;
          xdr_double(&xdr_handle, &p.z) ;
          
        }
      }ENDFORALL; 
      FORALL(local_faces, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          xdr_double(&xdr_handle, &p.x) ; 
          xdr_double(&xdr_handle, &p.y) ;
          xdr_double(&xdr_handle, &p.z) ;
          
        }
      }ENDFORALL;
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 12, MPI_COMM_WORLD) ;
        MPI_Recv(inner_nodes_buf,inner_nodes_buf_size,
                 MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <inner_nodes_sizes[i]*3; j++)
          xdr_double(&xdr_handle, &inner_nodes_buf[j]);
      }
      fclose (out);
      xdr_destroy(&xdr_handle) ;
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      FORALL(local_edges,cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          inner_nodes_buf[index++] = p.x;
          inner_nodes_buf[index++] = p.y;
          inner_nodes_buf[index++] = p.z;
         
        }
      }ENDFORALL;
      FORALL(local_cells,cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          inner_nodes_buf[index++] = p.x;
          inner_nodes_buf[index++] = p.y;
          inner_nodes_buf[index++] = p.z;
          
        }
      }ENDFORALL;
      FORALL(local_faces,cc){
       for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          vect3d p = inner_nodes[cc][i];
          inner_nodes_buf[index++] = p.x;
          inner_nodes_buf[index++] = p.y;
          inner_nodes_buf[index++] = p.z;
          
        }  
      }ENDFORALL;
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 12, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(inner_nodes_buf, inner_nodes_sizes[my_id]*3, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD) ;
      }  
    }//finish output inner nodes

    delete [] inner_nodes_buf;

   
  }
} ;
register_rule<node_output_file> register_node_output_file;
  

class face_output_file : public pointwise_rule {
 
  store<bool> face_output ;
  const_param<string> outfile_par ;
  const_store<Loci::FineFaces> fine_faces;

  
public:
  face_output_file(){
    name_store("face_output", face_output);
    name_store("outfile_par", outfile_par);
    name_store("fine_faces", fine_faces);

    
    input("outfile_par");
    input("fine_faces");

    output("face_output");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    
    FILE *out;
    XDR xdr_handle;
    
   
   
    if(num_procs == 1){

      out = fopen((*outfile_par).c_str(),"a");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;

      Loci::constraint faces, geom_cells;
      faces = Loci::exec_current_fact_db->get_variable("faces");
      geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
     
      //output offset c1, c2 
      int offset = 0;
      int c1, c2;
      FORALL(*geom_cells, cc){
      
        for(unsigned int f = 0; f < fine_faces[cc].size(); f++){
      
          c1 = fine_faces[cc][f][0];
          c2 = fine_faces[cc][f][1];
          xdr_int(&xdr_handle, &offset) ; 
           xdr_int(&xdr_handle, &c1) ;
           xdr_int(&xdr_handle, &c2) ;
           offset += fine_faces[cc][f].size()-2;
        }
        
      }ENDFORALL; 
      FORALL(*faces, cc){
        for(unsigned int f = 0; f < fine_faces[cc].size(); f++){
          c1 = fine_faces[cc][f][0];
          c2 = fine_faces[cc][f][1];
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &c1) ;
          xdr_int(&xdr_handle, &c2) ;
          offset += fine_faces[cc][f].size()-2;
        }
      }ENDFORALL;
      xdr_int(&xdr_handle, &offset) ;
      
      //output face2node
      int nindex = 0;
      
      FORALL(*geom_cells, cc){
         for(unsigned f = 0; f < fine_faces[cc].size(); f++){
           for(unsigned int i = 2; i < fine_faces[cc][f].size(); i++){
             nindex = fine_faces[cc][f][i]-1;
             xdr_int(&xdr_handle, &nindex) ;
           }
         }
       }ENDFORALL; 
       FORALL(*faces, cc){
         for(unsigned f = 0; f < fine_faces[cc].size(); f++){
           for(unsigned int i = 2; i < fine_faces[cc][f].size(); i++){
             nindex = fine_faces[cc][f][i]-1;
             xdr_int(&xdr_handle, &nindex) ;
           }
         }
       }ENDFORALL;
       
       fclose (out);
       xdr_destroy(&xdr_handle) ;
       return;
    }
      
  
    //next, output pos

    fact_db::distribute_infoP d =  Loci::exec_current_fact_db->get_distribute_info() ;
    Loci::constraint my_entities, faces, geom_cells, interior_faces, boundary_faces; 
    my_entities = d->my_entities ;
    faces = Loci::exec_current_fact_db->get_variable("faces");
    geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
    interior_faces =  Loci::exec_current_fact_db->get_variable("interior_faces");
    boundary_faces =  Loci::exec_current_fact_db->get_variable("boundary_faces");
    
   
    entitySet local_faces = *my_entities & *faces;
    entitySet local_cells =  *my_entities & *geom_cells;
    entitySet local_interior_faces = *my_entities & *interior_faces;
    entitySet local_boundary_faces = *my_entities & *boundary_faces;

    //process 0 open the file     
    if(my_id == 0){
      out = fopen((*outfile_par).c_str(),"a");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;
    }
    
    //start output offset, c1, c2
    //first output interior faces
    int offset = 0;
    int c1, c2;
    int num_local_interior_faces = 0;
    
    FORALL(local_cells, cc){
      num_local_interior_faces += fine_faces[cc].size();
    }ENDFORALL; 
    FORALL(local_interior_faces, cc){
      num_local_interior_faces += fine_faces[cc].size();
    }ENDFORALL;
    
    
    std::vector<int> interior_neib_sizes;
    interior_neib_sizes = Loci::all_collect_sizes(num_local_interior_faces);
    int interior_buf_size =  (*std::max_element(interior_neib_sizes.begin(), interior_neib_sizes.end()))*3;
    int* interior_buf = new int[interior_buf_size];

    cerr << currentMem() << " face_output " <<  Loci::MPI_rank << endl;
    //first process 0 write out its neib
    if(my_id == 0){
      
      FORALL(local_cells, cc){
        for(unsigned int f = 0; f < fine_faces[cc].size(); f++){
          c1 = fine_faces[cc][f][0];
          c2 = fine_faces[cc][f][1];
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &c1) ;
          xdr_int(&xdr_handle, &c2) ;
         
          offset += fine_faces[cc][f].size()-2;
        }
        
      }ENDFORALL; 
      FORALL(local_interior_faces, cc){
        for(unsigned int f = 0; f < fine_faces[cc].size(); f++){
          c1 = fine_faces[cc][f][0];
          c2 = fine_faces[cc][f][1];
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &c1) ;
          xdr_int(&xdr_handle, &c2) ;
          offset += fine_faces[cc][f].size()-2;
        }
      }ENDFORALL;
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 14, MPI_COMM_WORLD) ;
        MPI_Recv(interior_buf,interior_buf_size,
                 MPI_INT, i, 15, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <interior_neib_sizes[i]; j++){
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &interior_buf[3*j+1]);
          xdr_int(&xdr_handle, &interior_buf[3*j+2]);
          offset += interior_buf[3*j];
        }
      }

    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      
      FORALL(local_cells,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          interior_buf[index++] = fine_faces[cc][f].size()-2;
          interior_buf[index++] = fine_faces[cc][f][0];
          interior_buf[index++] = fine_faces[cc][f][1];
        }
      }ENDFORALL;
      FORALL(local_interior_faces,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          interior_buf[index++] = fine_faces[cc][f].size()-2;
          interior_buf[index++] = fine_faces[cc][f][0];
          interior_buf[index++] = fine_faces[cc][f][1];
        } 
      }ENDFORALL;
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 14, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(interior_buf, interior_neib_sizes[my_id]*3, MPI_INT, 0, 15, MPI_COMM_WORLD) ;
      }  
    }//finish output interior offset, c1, c2

    delete [] interior_buf;
    //then output boundary faces

    int num_local_boundary_faces = 0;
   
    
    FORALL(local_boundary_faces, cc){
      num_local_boundary_faces += fine_faces[cc].size();
    }ENDFORALL;
    
    
    std::vector<int> boundary_neib_sizes;
    
    boundary_neib_sizes = Loci::all_collect_sizes(num_local_boundary_faces);
    int boundary_buf_size =  (*std::max_element(boundary_neib_sizes.begin(), boundary_neib_sizes.end()))*3;
    int* boundary_buf = new int[boundary_buf_size];


    //first process 0 write out its neib
    if(my_id == 0){
      
      FORALL(local_boundary_faces, cc){
        for(unsigned int f = 0; f < fine_faces[cc].size(); f++){
          c1 = fine_faces[cc][f][0];
          c2 = fine_faces[cc][f][1];
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &c1) ;
          xdr_int(&xdr_handle, &c2) ;
          offset += fine_faces[cc][f].size()-2;
        }
      }ENDFORALL;
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 16, MPI_COMM_WORLD) ;
        MPI_Recv(boundary_buf,boundary_buf_size,
                 MPI_INT, i, 17, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <boundary_neib_sizes[i]; j++){
          xdr_int(&xdr_handle, &offset) ; 
          xdr_int(&xdr_handle, &boundary_buf[3*j+1]);
          xdr_int(&xdr_handle, &boundary_buf[3*j+2]);
          offset += boundary_buf[3*j];
        }
      }
      xdr_int(&xdr_handle, &offset) ; 
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      FORALL(local_boundary_faces,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          boundary_buf[index++] = fine_faces[cc][f].size()-2;
          boundary_buf[index++] = fine_faces[cc][f][0];
          boundary_buf[index++] = fine_faces[cc][f][1];
        } 
      }ENDFORALL;
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 16, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(boundary_buf, boundary_neib_sizes[my_id]*3, MPI_INT, 0, 17, MPI_COMM_WORLD) ;
      }  
    }//finish output boundary offset, c1, c2


    delete [] boundary_buf;

    
    //start interior face2node,
    int f2n_vec_size = 0;
    
    FORALL(local_cells, cc){
      for(unsigned int f = 0; f < fine_faces[cc].size(); f++) 
        f2n_vec_size += (fine_faces[cc][f].size()-2);
    }ENDFORALL; 
    FORALL(local_interior_faces, cc){
      for(unsigned int f = 0; f < fine_faces[cc].size(); f++)
        f2n_vec_size += (fine_faces[cc][f].size()-2);
    }ENDFORALL;
    

    std::vector<int> f2n_sizes;
    f2n_sizes = Loci::all_collect_sizes(f2n_vec_size);
   
    int f2n_buf_size =  *std::max_element(f2n_sizes.begin(), f2n_sizes.end());
    int* f2n_buf = new int[f2n_buf_size];
   
    //first process 0 write out its f2n
    if(my_id == 0){
      int nindex = 0;
      
      FORALL(local_cells, cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int i = 2; i < fine_faces[cc][f].size(); i++){
            nindex = fine_faces[cc][f][i]-1;
            xdr_int(&xdr_handle, &nindex) ; 
          }
        }
      }ENDFORALL; 
      FORALL(local_interior_faces, cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int i = 2; i < fine_faces[cc][f].size(); i++){
            nindex = fine_faces[cc][f][i]-1;
            xdr_int(&xdr_handle, &nindex) ; 
          }
        }
      }ENDFORALL;
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 18, MPI_COMM_WORLD) ;
        MPI_Recv(f2n_buf,f2n_buf_size,
                 MPI_INT, i, 19, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <f2n_sizes[i]; j++){
          xdr_int(&xdr_handle, &f2n_buf[j]) ; 
        }
      }
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      
      FORALL(local_cells,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int j = 2; j < fine_faces[cc][f].size(); j++){ 
            f2n_buf[index++] = fine_faces[cc][f][j]-1;
          }
        }
      }ENDFORALL;
      FORALL(local_interior_faces,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int j = 2; j < fine_faces[cc][f].size(); j++){ 
            f2n_buf[index++] = fine_faces[cc][f][j]-1;
          }
        }
      }ENDFORALL;
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 18, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(f2n_buf, f2n_sizes[my_id], MPI_INT, 0, 19, MPI_COMM_WORLD) ;
      }  
    }//finish output interior face2node

    delete [] f2n_buf;
    //start boundary face2node,

    int f2n_boundary_size = 0;
    FORALL(local_boundary_faces, cc){
      for(unsigned int f = 0; f < fine_faces[cc].size(); f++)
        f2n_boundary_size += fine_faces[cc][f].size()-2;
    }ENDFORALL;
    
    
    std::vector<int> f2n_boundary_sizes;
    f2n_boundary_sizes = Loci::all_collect_sizes(f2n_boundary_size);
   

    int f2n_boundary_buf_size =  *std::max_element(f2n_boundary_sizes.begin(), f2n_boundary_sizes.end());
    int* f2n_boundary_buf = new int[f2n_boundary_buf_size];
   
    //first process 0 write out its f2n
    if(my_id == 0){
      int nindex = 0;
      FORALL(local_boundary_faces, cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int i = 2; i < fine_faces[cc][f].size(); i++){
            nindex = fine_faces[cc][f][i]-1;
            xdr_int(&xdr_handle, &nindex) ; 
          }
        }
      }ENDFORALL;
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 20, MPI_COMM_WORLD) ;
        MPI_Recv(f2n_boundary_buf,f2n_boundary_buf_size,
                 MPI_INT, i, 21, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <f2n_boundary_sizes[i]; j++){
          xdr_int(&xdr_handle, &f2n_boundary_buf[j]) ; 
        }
      }
      fclose (out);
      xdr_destroy(&xdr_handle) ;
      
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      FORALL(local_boundary_faces,cc){
        for(unsigned f = 0; f < fine_faces[cc].size(); f++){
          for(unsigned int j = 2; j < fine_faces[cc][f].size(); j++){ 
            f2n_boundary_buf[index++] = fine_faces[cc][f][j]-1;
          }
        }
      }ENDFORALL;
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(f2n_boundary_buf, f2n_boundary_sizes[my_id], MPI_INT, 0, 21, MPI_COMM_WORLD) ;
      }  
    }//finish output boundary face2node
    
    delete [] f2n_boundary_buf;
  }
} ;
register_rule<face_output_file> register_face_output_file;




class pos_output_file : public pointwise_rule {
 
  store<bool> pos_output ;
  const_param<string> outfile_par ;
  const_param<int> npnts_par;
  const_param<int> ncells_par;
  const_param<int> nfaces_par;
  const_store<vect3d> pos;
  
 
  
public:
  pos_output_file(){
    name_store("pos_output", pos_output);
    name_store("outfile_par", outfile_par);
    name_store("npnts", npnts_par);
    name_store("ncells", ncells_par);
    name_store("nfaces", nfaces_par);
    name_store("pos", pos);
 
    
    input("outfile_par");
    input("npnts");
    input("ncells");
    input("nfaces");
    input("pos");
    
    output("pos_output");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
   
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    
    FILE *out;
    XDR xdr_handle;
    
   
   
    if(num_procs == 1){

      out = fopen((*outfile_par).c_str(),"w");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;
      // int dummy = 0;
      int ndm = 3;
      int nzones = 1;
      int npatch = 6;
      int maxppf = 0;
      int maxfpc = 0;
      int npnts, ncells, nfaces;
      
      
     
      
      npnts = *npnts_par + pos.domain().size();
      ncells = *ncells_par;
      nfaces = *nfaces_par;
      cout << " npnts, nfaces, ncells " << npnts << " " << nfaces <<" " << ncells << endl;
      //output header
      xdr_int(&xdr_handle, &ndm) ;
      xdr_int(&xdr_handle, &nzones) ;
      xdr_int(&xdr_handle, &npatch) ;
      xdr_int(&xdr_handle, &npnts) ; //num of Nodes
      xdr_int(&xdr_handle, &nfaces) ;
      xdr_int(&xdr_handle, &ncells) ;
      xdr_int(&xdr_handle, &maxppf) ;
      xdr_int(&xdr_handle, &maxfpc) ;

      //output original nodes
      FORALL(pos.domain(), cc){
        vect3d p = pos[cc];
        xdr_double(&xdr_handle, &p.x) ; 
        xdr_double(&xdr_handle, &p.y) ;
        xdr_double(&xdr_handle, &p.z) ;
      }ENDFORALL;
      
      fclose (out);
      xdr_destroy(&xdr_handle) ;
      return;
    }
      
  
    //next, output pos

    fact_db::distribute_infoP d =  Loci::exec_current_fact_db->get_distribute_info() ;
    Loci::constraint my_entities, faces, geom_cells, interior_faces, boundary_faces; 
    my_entities = d->my_entities ;
    
    //reorder store first, from local to io entities
    storeRepP posRep = pos.Rep();
    store<vect3d> pos_io;
    pos_io =  collect_reorder_store(posRep, *(Loci::exec_current_fact_db));
    entitySet nodes = pos_io.domain();
    //compute the size of buf
    
    int local_size = nodes.size();
    std::vector<int> sort_max ;
    sort_max = Loci::all_collect_sizes(local_size) ;

    int num_original_nodes = 0;
    for(unsigned int i = 0; i < sort_max.size(); i++)num_original_nodes += sort_max[i];
    
    int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
    total_size = total_size*3; //x, y, z value
    
 //process 0 open the file and write out xdr header    
    if(my_id == 0){
      out = fopen((*outfile_par).c_str(),"w");
      if(out == NULL) {
        cerr << "can't open " << *outfile_par << " for writing" << endl ;
        exit(-1) ;
      }
      xdrstdio_create(&xdr_handle,out,XDR_ENCODE) ;
      //      int dummy = 0;
      int ndm = 3;
      int nzones = 1;
      int npatch = 6;
      int maxppf = 0;
      int maxfpc = 0;
      int npnts, ncells, nfaces;
     
      npnts = *npnts_par + num_original_nodes;
      ncells = *ncells_par;
      nfaces = *nfaces_par;
      cout << " npnts, nfaces, ncells " << npnts << " " << nfaces <<" " << ncells << endl;
      xdr_int(&xdr_handle, &ndm) ;
      xdr_int(&xdr_handle, &nzones) ;
      xdr_int(&xdr_handle, &npatch) ;
      xdr_int(&xdr_handle, &npnts) ; //num of Nodes
      xdr_int(&xdr_handle, &nfaces) ;
      xdr_int(&xdr_handle, &ncells) ;
      xdr_int(&xdr_handle, &maxppf) ;
      xdr_int(&xdr_handle, &maxfpc) ;
    }
    
   
    //allocate buf
    double* buf = new double[total_size];
    cerr << currentMem() << " pos_output " <<  Loci::MPI_rank << endl;
    
    //first process 0 write out its pos_io
    if(my_id == 0){
      FORALL(nodes, cc){
        xdr_double(&xdr_handle, &pos_io[cc].x) ; 
        xdr_double(&xdr_handle, &pos_io[cc].y) ;
        xdr_double(&xdr_handle, &pos_io[cc].z) ;
       
      }ENDFORALL;
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
       	  MPI_Send(&flag, 1, MPI_INT, i, 10, MPI_COMM_WORLD) ;
          MPI_Recv(buf,total_size, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &status) ;
          //write the buffer out
          for(int j = 0; j <sort_max[i]*3; j++) {
            xdr_double(&xdr_handle, &buf[j]);
          }
          
      }
      fclose(out);
      xdr_destroy(&xdr_handle) ; 
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      FORALL(nodes,cc){
        buf[index++] = pos_io[cc].x;
        buf[index++] = pos_io[cc].y;
        buf[index++] = pos_io[cc].z;
      }ENDFORALL;
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status) ;
      if(flag) { 
        MPI_Send(buf, sort_max[my_id]*3, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD) ;
      }  
    }//finish output pos
    
    delete [] buf;
  }
} ;
register_rule<pos_output_file> register_pos_output_file;
  

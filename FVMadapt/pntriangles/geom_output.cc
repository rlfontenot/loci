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

#include "coeff.h"

using std::vector;
using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using Loci::storeRepP;
using Loci::constraint;

//typedef Loci::vector3d<double> vect3d;
namespace Loci{
  //storeRepP collect_reorder_store(storeRepP &sp, fact_db &facts);
std::vector<int> all_collect_sizes(int size);
}



class bnode_output_file : public pointwise_rule {
  
  store<bool> bnode_output ;
  const_param<string> outfile_par ;
  const_store<vect3d> pos;
  
  
  
public:
  bnode_output_file(){
    name_store("bnode_output", bnode_output);
    name_store("outfile_par", outfile_par);
    name_store("pos", pos);
    input("outfile_par");
    input("pos");
    output("bnode_output");
    constraint("boundary_nodes");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
        
    if(num_procs == 1){
      ofstream ofile((*outfile_par).c_str(), ios::out);
      ofile.precision(10);
      ofile << seq.size()<< endl;
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        ofile << pos[*ei].x <<' ' << pos[*ei].y <<' '<< pos[*ei].z  << endl;
      }
      
      ofile.close();
      return;
    }
    
    
    
    std::vector<int> boundary_nodes_sizes;
    boundary_nodes_sizes = Loci::all_collect_sizes(seq.size());
    
    int boundary_nodes_buf_size =  (*std::max_element(boundary_nodes_sizes.begin(), boundary_nodes_sizes.end()))*3;
    int num_boundary_nodes = 0;
    for(unsigned int i = 0; i < boundary_nodes_sizes.size(); i++){
      num_boundary_nodes += boundary_nodes_sizes[i];
    }
    
    double* boundary_nodes_buf = new double[boundary_nodes_buf_size];
    //process 0 open the file     
    if(my_id == 0){
      ofstream ofile((*outfile_par).c_str(), ios::out);
      ofile.precision(10);
      //write out num_boundary_nodes
      ofile << num_boundary_nodes << endl;
      
      //write out my own boundary nodes
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        ofile << pos[*ei].x <<' ' << pos[*ei].y <<' '<< pos[*ei].z   << endl;
      }
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 12, MPI_COMM_WORLD) ;
        MPI_Recv(boundary_nodes_buf,boundary_nodes_buf_size,
                 MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <boundary_nodes_sizes[i]; j++){
          ofile << boundary_nodes_buf[3*j] << ' '
                << boundary_nodes_buf[3*j+1] << ' ' <<boundary_nodes_buf[3*j+2]<< endl;
        }
      }
      
      ofile.close();
      
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        
        boundary_nodes_buf[index++] = pos[*ei].x;
        boundary_nodes_buf[index++] = pos[*ei].y;
        boundary_nodes_buf[index++] = pos[*ei].z;
      }
      
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 12, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(boundary_nodes_buf, boundary_nodes_sizes[my_id]*3, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD) ;
      }  
    }//finish output boundary nodes
    
    delete [] boundary_nodes_buf;
    
    
  }
} ;
register_rule<bnode_output_file> register_bnode_output_file;


class bface_output_file : public pointwise_rule {
  
  store<bool> bface_output ;
  const_param<string> outfile_par ;
  const_store<int> boundary_node_index;
  const_store<std::vector<Loci::GeomCoeff> > geom;
  const_multiMap face2node;
  const_store<vect3d> pos;
  
public:
  bface_output_file(){
    name_store("bface_output", bface_output);
    name_store("outfile_par", outfile_par);
    name_store("geom", geom);
    name_store("boundary_node_index", boundary_node_index);
    name_store("face2node", face2node);
    name_store("pos", pos);
    input("outfile_par");
    input("geom");
    input("face2node->(boundary_node_index,pos)");
    output("bface_output");
    constraint("geometry_BC");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    
    //compute how many triangle boundary faces each process has
    int num_local_boundary_faces = 0;
    for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
      num_local_boundary_faces += face2node[*ei].size()-2;
    }

    
    std::vector<int> boundary_faces_sizes;
    boundary_faces_sizes = Loci::all_collect_sizes(num_local_boundary_faces);

    //total number of triangle boundary faces on all processes
    int num_boundary_faces = 0;
    for(unsigned int i = 0; i < boundary_faces_sizes.size(); i++){
      num_boundary_faces += boundary_faces_sizes[i];
    }
    
    
    int boundary_faces_buf_size =  (*std::max_element(boundary_faces_sizes.begin(), boundary_faces_sizes.end()))*3;
    int* boundary_faces_buf = new int[boundary_faces_buf_size];
    
    ofstream ofile; 
    
    
    if(num_procs == 1){
      
      ofile.open((*outfile_par).c_str(), ios::app);
      ofile << num_boundary_faces << endl;
      
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for( int triId = 1; triId <= (face2node[*ei].size()-2); triId++){
          ofile << boundary_node_index[face2node[*ei][0]]<<' '
                << boundary_node_index[face2node[*ei][triId]]<<' '
                <<boundary_node_index[face2node[*ei][triId+1]] << endl;
        }
      }
      ofile.precision(10);
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for(unsigned int f = 0; f < geom[*ei].size(); f++){
          Loci::GeomCoeff g = geom[*ei][f];
        
          ofile << g.b300<<' '
                << g.b030<<' '
                << g.b003<<' '
                << g.b210<<' '
                << g.b120<<' '
                << g.b021<<' '
                << g.b012<<' '
                << g.b102<<' '
                << g.b201<<' '
                << g.b111<<endl;
        }
      }


     //  for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
         
//        for( int triId = 1; triId <= (face2node[*ei].size()-2); triId++){
//          cout<<"face2node: " <<boundary_node_index[face2node[*ei][0]]<<' '
//              << boundary_node_index[face2node[*ei][triId]]<<' '
//              <<boundary_node_index[face2node[*ei][triId+1]] << endl;
//          cout<< pos[face2node[*ei][0]]<<endl;
//          cout  << pos[face2node[*ei][triId]]<<endl;
//          cout   <<pos[face2node[*ei][triId+1]] << endl;
       
//          Loci::GeomCoeff g = geom[*ei][triId-1];
//          cout << "geom: " << endl;
//          cout  << g.b300<<endl;
//          cout   << g.b030<<endl;
//          cout      << g.b003<<endl;
         
//        }
//        cout <<endl<<endl;
//       }
      ofile.close();
      return;
    }
    
    
    // output face2node
    
    
    //process 0 open the file     
    if(my_id == 0){
      ofile.open((*outfile_par).c_str(),ios::app);
      ofile << num_boundary_faces << endl;
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for( int triId = 1; triId <= (face2node[*ei].size()-2); triId++){  
          ofile << boundary_node_index[face2node[*ei][0]]<<' '
                <<boundary_node_index[face2node[*ei][triId]]<<' '
                <<boundary_node_index[face2node[*ei][triId+1]] << endl;
        }
      }
      
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 14, MPI_COMM_WORLD) ;
        MPI_Recv(boundary_faces_buf,boundary_faces_buf_size,
                 MPI_INT, i, 15, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <boundary_faces_sizes[i]; j++){
          ofile << boundary_faces_buf[3*j] <<' '
                << boundary_faces_buf[3*j+1] <<' '
                << boundary_faces_buf[3*j+2] <<endl;
        }
        
        
      }
      
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for( int triId = 1; triId <= (face2node[*ei].size()-2); triId++){ 
          boundary_faces_buf[index++] =  boundary_node_index[face2node[*ei][0]];
          boundary_faces_buf[index++] = boundary_node_index[face2node[*ei][triId]] ;
          boundary_faces_buf[index++] =  boundary_node_index[face2node[*ei][triId+1]];
        }
      }
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 14, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(boundary_faces_buf, boundary_faces_sizes[my_id]*3, MPI_INT, 0, 15, MPI_COMM_WORLD) ;
      }  
    }//finish output face2node
    
    delete [] boundary_faces_buf;
    
    //then output boundary faces coefficient
    int boundary_face_buf_size =  (*std::max_element(boundary_faces_sizes.begin(), boundary_faces_sizes.end()))*30;
    double* boundary_face_buf = new double[boundary_face_buf_size];
    
    if(my_id == 0){
       ofile.precision(10);
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for(unsigned int f = 0; f < geom[*ei].size(); f++){
          Loci::GeomCoeff g = geom[*ei][f];
        
          ofile << g.b300<<' '
                << g.b030<<' '
                << g.b003<<' '
                << g.b210<<' '
                << g.b120<<' '
                << g.b021<<' '
                << g.b012<<' '
                << g.b102<<' '
                << g.b201<<' '
                << g.b111<<' ' <<endl;
        }
      }
      
      
      //then process 0 recv data from other processes
      for(int i = 1; i < num_procs; i++){
        MPI_Status status;
        int flag = 1;
        MPI_Send(&flag, 1, MPI_INT, i, 16, MPI_COMM_WORLD) ;
        MPI_Recv(boundary_face_buf,boundary_face_buf_size,
                 MPI_DOUBLE, i, 17, MPI_COMM_WORLD, &status) ;
        //write the buffer out
        for(int j = 0; j <boundary_faces_sizes[i]; j++){
          for(int k = 0; k < 10; k++){
            vect3d p;
            p.x =  boundary_face_buf[30*j+k*3];
            p.y =  boundary_face_buf[30*j+k*3+1];
            p.z =  boundary_face_buf[30*j+k*3+2];
            ofile << p << ' ';
          }
          ofile << endl;
        }
        
        ofile.close();
      }
    }else{
      //process [1,numprocs-1] send data to process 0
      MPI_Status status ;
      //pack the buf
      int index=0;
      
      for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
        for(unsigned int f = 0; f < geom[*ei].size(); f++){
          Loci::GeomCoeff g = geom[*ei][f];
          boundary_face_buf[index++] = g.b300.x;
          boundary_face_buf[index++] = g.b300.y;
          boundary_face_buf[index++] = g.b300.z;
          
          boundary_face_buf[index++] = g.b030.x;
          boundary_face_buf[index++] = g.b030.y;
          boundary_face_buf[index++] = g.b030.z;
          
          boundary_face_buf[index++] = g.b003.x;
          boundary_face_buf[index++] = g.b003.y;
          boundary_face_buf[index++] = g.b003.z;
          
          boundary_face_buf[index++] = g.b210.x;
          boundary_face_buf[index++] = g.b210.y;
          boundary_face_buf[index++] = g.b210.z;

          boundary_face_buf[index++] = g.b120.x;
          boundary_face_buf[index++] = g.b120.y;
          boundary_face_buf[index++] = g.b120.z;
          
          boundary_face_buf[index++] = g.b021.x;
          boundary_face_buf[index++] = g.b021.y;
          boundary_face_buf[index++] = g.b021.z;
          
          boundary_face_buf[index++] = g.b012.x;
          boundary_face_buf[index++] = g.b012.y;
          boundary_face_buf[index++] = g.b012.z;
          
          boundary_face_buf[index++] = g.b102.x;
          boundary_face_buf[index++] = g.b102.y;
          boundary_face_buf[index++] = g.b102.z;
          
          boundary_face_buf[index++] = g.b201.x;
          boundary_face_buf[index++] = g.b201.y;
          boundary_face_buf[index++] = g.b201.z;
          
          boundary_face_buf[index++] = g.b111.x;
          boundary_face_buf[index++] = g.b111.y;
          boundary_face_buf[index++] = g.b111.z;
          
        }
      }
      int flag = 0 ;
      MPI_Recv(&flag, 1, MPI_INT, 0, 16, MPI_COMM_WORLD, &status) ;
      if(flag) {
        MPI_Send(boundary_face_buf, boundary_faces_sizes[my_id]*30, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD) ;
      }  
    }//finish output face geometry
    
    delete [] boundary_face_buf;
    
    
    
    
  }
} ;
register_rule<bface_output_file> register_bface_output_file;

  

//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#ifndef DISTRIBUTE_IO_H
#define DISTRIBUTE_IO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <parSampleSort.h>

#include <store_rep.h>
#include <DMap.h>
#include <fact_db.h>
#include <mpi.h>
#include <store.h>
#include <MapVec.h>
#include <multiStore.h>

namespace Loci {

  void redistribute_write_container(hid_t file_id, std::string vname,
                                    Loci::storeRepP var, fact_db &facts) ;
  void read_container_redistribute(hid_t file_id, std::string vname,
                                   Loci::storeRepP var, entitySet read_set,
                                   fact_db &facts) ;

  inline hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id) {
    hid_t file_id = 0 ;
    if(Loci::MPI_rank==0) {
      file_id = H5Fcreate(name,flags,create_id,access_id) ;
      if(file_id < 0) {
	cerr << "H5Fcreate unable to create file '" << name << "'" << endl ;
	Loci::Abort() ;
      }
    }
    return file_id ;
  }

  inline hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fcreate(name,flags,create_id,access_id) ;
    else
      return 0 ;
  }

  inline hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id) {
    if(Loci::MPI_rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
  }

  inline hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id,
                            MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
  }

  inline herr_t hdf5CloseFile(hid_t file_id) {
    if(Loci::MPI_rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }

  inline herr_t hdf5CloseFile(hid_t file_id, MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }
    
  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var, fact_db &facts) {

    redistribute_write_container(file_id,vname,var,facts) ;
  }
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet, fact_db &facts) {
    read_container_redistribute(file_id,vname,var,readSet, facts) ;
  }

  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::writeContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      redistribute_write_container(file_id,vname,var,
                                   *Loci::exec_current_fact_db) ;
  }
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::readContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      read_container_redistribute(file_id,vname,var,readSet,
                                  *Loci::exec_current_fact_db) ;
  }

  void writeContainerRAW(hid_t file_id, std::string vname,
                         storeRepP var, MPI_Comm comm) ;

  void readContainerRAW(hid_t file_id, std::string vname,
                        storeRepP var, MPI_Comm comm ) ;

  template<class T> void writeUnorderedVector(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v,
                                              MPI_Comm comm) {
    int my_rank = 0 ;
    MPI_Comm_rank(comm,&my_rank) ;
    int procs = 1 ;
    MPI_Comm_size(comm,&procs) ;

    //serial version
    if(procs==1){
      
      hsize_t array_size_combined = v.size() ;
      if(array_size_combined == 0)
        return ;
     
      int rank = 1 ;
      hsize_t dimension = array_size_combined ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
     
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = v.size() ;
      hid_t datatype = dp->get_hdf5_type() ;
     
#ifdef H5_USE_16_API
      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT) ;
#else
      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
     
      if(count != 0) {
     
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
     
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
     
        H5Dwrite(dataset,datatype,memspace,dataspace,
                 H5P_DEFAULT, &v[0]) ;
     
        H5Sclose(memspace) ;
     
      }
     
      
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      return;
    }

    size_t local_size = v.size() ;
    std::vector<size_t> recv_sizes(procs) ;
    MPI_Gather(&local_size,sizeof(size_t),MPI_BYTE,
               &recv_sizes[0],sizeof(size_t),MPI_BYTE,0,comm) ;

    if(my_rank == 0) {
      hsize_t array_size = 0 ;
      for(int i=0;i<procs;++i)
        array_size += recv_sizes[i] ;
      if(array_size == 0)
        return ;
      int rank = 1 ;
      hsize_t dimension = array_size ;

      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = recv_sizes[0] ;
      hid_t datatype = dp->get_hdf5_type() ;
#ifdef H5_USE_16_API
      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT) ;
#else
      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
      if(count != 0) {
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,datatype,memspace,dataspace,
                 H5P_DEFAULT, &v[0]) ;
        H5Sclose(memspace) ;
      }
      for(int i=1;i<procs;++i) {
        start += recv_sizes[i-1] ;
        if(recv_sizes[i] == 0)
          continue ;
        int flag = 0 ;
        MPI_Send(&flag,1,MPI_INT,i,0,comm) ;
        std::vector<T> rv(recv_sizes[i]) ;
        MPI_Status mstat ;
        MPI_Recv(&rv[0],sizeof(T)*recv_sizes[i],MPI_BYTE,i,0,comm,
                 &mstat) ;
        count = recv_sizes[i] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,datatype,memspace,dataspace,
                 H5P_DEFAULT, &rv[0]) ;
        H5Sclose(memspace) ;
      }

      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
    } else {
      if(local_size == 0)
        return ;

      int flag = 0;
      MPI_Status mstat ;
      MPI_Recv(&flag,1,MPI_INT,0,0,comm,&mstat) ;
      MPI_Send(&v[0],sizeof(T)*local_size,MPI_BYTE,0,0,comm) ;
    }
  }
 
   
  template<class T> void writeUnorderedVector(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v) {
    writeUnorderedVector(group_id,element_name,v,MPI_COMM_WORLD) ;
  }
  
 
  
  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) ;
  
  hid_t createUnorderedFile(const char * filename, entitySet set,
                            fact_db &facts) ;

  inline hid_t createUnorderedFile(const char * filename, entitySet set) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::createUnorderedFile()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
      return -1 ;
    } else
      return createUnorderedFile(filename, set, *Loci::exec_current_fact_db) ;
  }

  void closeUnorderedFile(hid_t file_id) ;

  template<class T> void writeUnorderedStore(hid_t file_id,
                                             const_store<T> &s, entitySet set,
                                             const char *name) {
    std::vector<T> v(set.size()) ;
    size_t c = 0 ;
    FORALL(set,ii) {
      v[c++] = s[ii] ;
    } ENDFORALL ;
    writeUnorderedVector(file_id,name,v) ;
  }

  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells,
                                 fact_db &facts) ;

  inline
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::parallelWriteGridTopology()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      parallelWriteGridTopology(filename, upperRep, lowerRep, boundary_mapRep,
                                face2nodeRep, refRep, bnamesRep, posRep,
                                localCells, *Loci::exec_current_fact_db) ;
  }
  
  //open /output/$bc_name/$file_name
  hid_t open_boundary_file(std::string bc_name,
                           std::string file_name
                           );
    
  //get boundary faces that belong to a boundary surface current_bc
  entitySet get_boundary_faces(std::string current_bc,//boundary name
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset //all boundary faces 
                               );

  //get boundary nodes that belong to a boundary surface current_bc
  entitySet get_boundary_nodes(std::string current_bc,//boundary name
                               storeRepP face2nodeRep,
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset, //all boundary faces 
                               fact_db &facts );

  void writeBoundaryTopo(hid_t file_id, //file_id of this boudnary surface
                         storeRepP face2nodeRep,
                         entitySet bfaces, //boundary faces belong to this surface 
                         fact_db &facts ); 
  


 
  //return value at cut position 
  template<class T>  inline T interpolate_val( double t, //weight of first node
                                               T a, //value at first node
                                               T b //value at second node
                                               ){
    
    return t*a + (1.0 - t)*b ;
  }
 
  
 
  
  struct CutPlane {
    storeRepP edgesWeight; //it's a store<double> allocated on edgesCut, containing the weight for interpoplation for each edge in edgesCut
    std::vector<std::vector<int > > faceLoops;  //loops formed, the values stored are edge ids, which is either local edge entity or negated index to inner_edges
    
    CutPlane( storeRepP ew, 
              std::vector<std::vector<int > >& fl){
      edgesWeight = ew;
      faceLoops = fl;
    }
    
    CutPlane(const CutPlane& cp){
      edgesWeight = cp.edgesWeight ;
      faceLoops =  cp.faceLoops;
    }
    
    CutPlane(){}
    
  };
  
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP posRep,
                       entitySet localCells,//all geom_cells
                       fact_db &facts);

  inline
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP posRep,
                       entitySet localCells//all geom_cells
                       ) {

    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::getCutPlane()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
      return CutPlane() ;
    } else return getCutPlane( upperRep,
                               lowerRep,
                               boundary_mapRep,
                               face2nodeRep,
                               face2edgeRep,
                               edge2nodeRep,
                               posRep,
                               localCells,//all geom_cells
                               *Loci::exec_current_fact_db);
  }

  // Convert container from local numbering to output file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // fact_db pointer  (facts)
  // MPI Communicator
  storeRepP Local2FileOrder_output(storeRepP sp, entitySet dom,
                                   fact_db& facts, MPI_Comm comm);
 
  template<class T>   void writeCutPlaneNodalVal(hid_t file_id,
                                                 std::string element_name,
                                                 storeRepP face2nodeRep,
                                                 storeRepP edge2nodeRep,
                                                 const_store<T> & pos,
                                                 const Loci::CutPlane &cp,
                                                 fact_db &facts){
    
   
    
    const_multiMap face2node(face2nodeRep) ;
    const_MapVec<2> edge2node(edge2nodeRep);
    const_store<double> edgesWeight(cp.edgesWeight); //the weight for interpoplation for each edgesCut, allocated on edgesCut
    entitySet edgesCut = edgesWeight.domain();
        
    
    //check the domain
    if((edgesCut-edge2node.domain())!=EMPTY){
      debugout<< "ERROR: the domain of edge2node is smaller than cp.edgesCut"<<endl;
    }
    
    //compute the cutting positions of edges 
    store<T> edge_pos;
    edge_pos.allocate(edgesCut);
    FORALL(edgesCut, e){
      double w =edgesWeight[e];
      T a = pos[edge2node[e][0]];
      T b = pos[edge2node[e][1]];
      T p = interpolate_val(w, a, b);
      edge_pos[e] = p;
    }ENDFORALL;
    
   
    //transform the store into output order
    store<T> gedge_pos;
    storeRepP geposRep =  Local2FileOrder_output(edge_pos.Rep(),  edgesCut, 
                                                 facts, MPI_COMM_WORLD);
       
    if(geposRep == NULL){
      gedge_pos .allocate(EMPTY);
    }else{
      gedge_pos = geposRep;
    }
   
   
    //get positions std::vector
    entitySet local_edges_cut = gedge_pos.domain();
    int num_edge_nodes = local_edges_cut.size();
    std::vector<T>  vpos(num_edge_nodes);
    int cnt = 0 ;
    entitySet::const_iterator ei ;
    for(ei=local_edges_cut.begin();ei!=local_edges_cut.end();++ei)
      vpos[cnt++] = gedge_pos[*ei];
      
    //write out the vector
    writeUnorderedVector(file_id, element_name.c_str(), vpos) ;
  }
  
  template<class T>   void writeCutPlaneNodalVal(hid_t file_id,
                                                 std::string element_name,
                                                 storeRepP face2nodeRep,
                                                 storeRepP edge2nodeRep,
                                                 const_store<T> & pos,
                                                 const Loci::CutPlane &cp
                                                 ){
    
    writeCutPlaneNodalVal( file_id,
                           element_name,
                           face2nodeRep,
                           edge2nodeRep,
                           pos,
                           cp,
                           *Loci::exec_current_fact_db );
  }
  
                          
  void writeCutPlaneTopo(hid_t bc_id,
                         const CutPlane& cp,
                         fact_db &facts) ;
   
}
#endif

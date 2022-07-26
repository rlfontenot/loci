//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#define io_performance
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi_containerIO.h>
#include <execute.h>
#include <dist_internal.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include <distribute_io.h>
#include "dist_tools.h"
#include <fact_db.h>
#include <constraint.h>
#include <string>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


#ifndef MPI_STUBB

namespace Loci {
  extern entitySet BcastEntitySet(entitySet set, int root, MPI_Comm comm);
  extern std::vector<int> simplePartitionVec(int mn, int mx, int p);
  extern std::vector<Loci::entitySet> simplePartition(int mn, int mx, MPI_Comm comm);
  
  typedef struct {
    unsigned long vec_size;
    unsigned long dom_size;
    unsigned long data_size;
  } store_header;

  static void handle_error(int errcode, string str)
  {
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    cerr <<  str << " : " << msg << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
#define MPI_CHECK(fn) { int errcode; errcode = (fn); if (errcode != MPI_SUCCESS) handle_error(errcode, #fn ); }
  
  // similar to hdf5_ReadDomain in hdf5_readwrite.cc

  //only process 0 will call this function to read in header and domain
  //Notice:: this function need to be called right after file is opened, fh not updated
  //MPI_file_read() is used, it uses individual file pointer, blocking, non-collective
  void pmpi_ReadDomain( MPI_File fh,  store_header& header,  entitySet &eset)
  {
    // read in the whole header and domain
   
    MPI_CHECK(MPI_File_read(fh, &header, sizeof(header), MPI_BYTE, MPI_STATUS_IGNORE));
    
    unsigned long    dimension = header.dom_size;
    eset = EMPTY;
    if(dimension ==0)return;
    
    int *data = new int[dimension];

    
    MPI_CHECK(MPI_File_read(fh, data, dimension, MPI_INT, MPI_STATUS_IGNORE));
    eset = EMPTY;
    for(size_t i=0;i< dimension;i+=2){
      eset |= interval(data[i],data[i+1]);
    }
    delete [] data;
  }


  //only process 0 will call this function to write header and domain,
  //Notice:: this function need to be called right after file is opened, fh not updated
  //MPI_file_write() is used, it uses individual file pointer, blocking, non-collective
  void pmpi_WriteDomain(MPI_File fh, const entitySet &en, const store_header& header )
  {

   
    MPI_CHECK(MPI_File_write(fh,
                             &header, sizeof(header), MPI_BYTE,
                             MPI_STATUS_IGNORE) );
     
    int num_intervals = en.num_intervals(); 
    if( num_intervals < 1) return;
    unsigned long  dimension = num_intervals*2; //size of 1D Array
    
    interval *it = new interval[num_intervals];
    
    int *data = new int[num_intervals*2];      
    
    for(int i=0;i<num_intervals;i++){
      it[i]       = en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }
    MPI_CHECK(MPI_File_write(fh,
                             data, dimension, MPI_INT,
                             MPI_STATUS_IGNORE) );
    delete [] data;
    delete [] it;
  }



  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);
  

  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm) ;

  int getMinFileNumberFromLocal(entitySet read_set,
                                fact_db::distribute_infoP dist ) ;
  
  void pmpi_write_store(MPI_File fh, Loci::storeRepP qrep, entitySet dom, int offset, MPI_Comm comm, int xfer_type) {
#ifdef io_performance 
    MPI_Barrier(MPI_COMM_WORLD);
    stopWatch s;
    s.start();
#endif

    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;
 
    // Shift domain by offset
    entitySet dom_file = dom >> offset ;

    // Compute overall domain across processors
    std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
    entitySet q_dom;
    for(int i = 0; i < np; i++)
      q_dom += dom_vector[i];

    // If nothing to write, don't proceed
    if(q_dom == EMPTY)
      return ;

    //write frame info
    Loci::frame_info fi = qrep->get_frame_info() ;
    
    //get arr_sizes and total_arr_size
    int pack_size = qrep->pack_size(dom);
    std::vector<int> arr_sizes = Loci::all_collect_sizes(pack_size,comm) ;
    size_t tot_arr_size = 0 ;
    for(int i = 0; i < np; ++i)
      tot_arr_size += size_t(max(0,arr_sizes[i])) ;

    store_header header;
    
    header.vec_size = fi.size;
    header.data_size = tot_arr_size;  //count in bytes
    header.dom_size = 2*q_dom.num_intervals(); //how many int

    
    if(prank == 0)
      pmpi_WriteDomain(fh, q_dom, header);

    //need set moffset
    MPI_Offset moffset;
    moffset =(MPI_Offset) sizeof(header); //skip header
    moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain
    for(int i = 0; i < prank; i++)moffset += (MPI_Offset)arr_sizes[i]; 
 
    std::vector<unsigned char> buffer(pack_size);
    int loc_pack = 0; 
    //fill the buffer
    qrep->pack((void*)&buffer[0], loc_pack, pack_size, dom);

    if(xfer_type == DXFER_COLLECTIVE_IO){
      MPI_CHECK( MPI_File_write_at_all(fh, moffset, &buffer[0], pack_size, MPI_BYTE, 
                                       MPI_STATUS_IGNORE));
    }else{
      MPI_CHECK(MPI_File_write_at(fh, moffset, &buffer[0], pack_size, MPI_BYTE, 
                                  MPI_STATUS_IGNORE));
    }
#ifdef io_performance  
    MPI_Barrier(MPI_COMM_WORLD);
    double wall_time = s.stop();
    if(prank == 0) std::cerr << "                                                    mpi parallel time to write_store : "  << wall_time << std::endl; 
#endif
  }


  


  void pmpi_redistribute_write_container( std::string& filename,
                                          Loci::storeRepP var, Loci::fact_db &facts, int xfer_type) {

    if(var->RepType() == Loci::PARAMETER){
      if(MPI_rank==0) cerr << " mpi container io doesn't provide parameter io" << endl;
      return ;
    }
    
    MPI_File fh;
    MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                            MPI_MODE_CREATE|MPI_MODE_WRONLY, Loci::PHDF5_MPI_Info, &fh));
    int offset = 0 ;         
    Loci::fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist == 0) {
      pmpi_write_store(fh,var,var->domain(),offset,MPI_COMM_WORLD, xfer_type);
      return ;
    }
    
  
    // Redistribute container to map from local to global numbering
    if(var->RepType() != Loci::PARAMETER && Loci::MPI_processes != 1) {
      // parallel store write.. reorder to file numbering then write out
      // reorder from local to file numbering
      
      entitySet dom = var->domain() ;
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      Loci::storeRepP vardist = Local2FileOrder(var,dom,offset,dist,MPI_COMM_WORLD) ;
      // Write out container that has been distributed in the file numbering
      pmpi_write_store(fh,vardist,vardist->domain(),offset,MPI_COMM_WORLD, xfer_type) ;//group_id is replaced by fh
    }else if(var->RepType() != Loci::PARAMETER && Loci::MPI_processes == 1){
      pmpi_write_store(fh,var, var->domain(),offset,MPI_COMM_WORLD, xfer_type);
    }
    
    MPI_File_close(&fh);
  }

 

  void pmpi_read_store(MPI_File fh, Loci::storeRepP qrep, int &offset, MPI_Comm comm, int xfer_type) {
#ifdef io_performance
    MPI_Barrier(MPI_COMM_WORLD);
    Loci::stopWatch s;
    s.start();
#endif
    offset = 0 ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    store_header header;
    // Here we read in a store container.  First lets read in the domain
    entitySet q_dom ;
    if(prank == 0) pmpi_ReadDomain(fh, header, q_dom) ; //this function is called right after file is opened
    q_dom = BcastEntitySet(q_dom,0,comm) ;
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, comm);

    //each process do partition and compute dom
    std::vector<int> interval_sizes ;
    entitySet dom ;
    if(q_dom != EMPTY) {
      std::vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
      for(int i=0;i<np;++i) {
        entitySet qset = ptn[i] &q_dom ;
        interval_sizes.push_back(qset.size()) ;
      }
      dom = ptn[prank] &q_dom ;
    } else
      for(int i=0;i<np;++i)
        interval_sizes.push_back(0) ;

    if(q_dom==EMPTY) {
      qrep->allocate(q_dom) ;
      return ;
    }
    offset = dom.Min() ;
    dom <<= offset ;
  
    //allocate qrep
    qrep->allocate(dom) ;
    if(header.vec_size > 1)
      qrep->set_elem_size(header.vec_size) ;
    //find pack_size
    int pack_size = qrep->pack_size(dom);
    std::vector<int> arr_sizes = Loci::all_collect_sizes(pack_size,comm) ;

    //set moffset
    MPI_Offset moffset;
    moffset =(MPI_Offset) sizeof(header); //skip header
    moffset += (MPI_Offset) sizeof(int)*header.dom_size; //skip domain
    //first compute  start
    for(int p = 0; p < prank; ++p){
      moffset += (MPI_Offset) arr_sizes[p];
    }
    
    //read in buffer
    int loc_unpack = 0;
    vector<unsigned char> buffer(pack_size) ;
    if(xfer_type == DXFER_COLLECTIVE_IO){
      MPI_CHECK(MPI_File_read_at_all(fh, moffset, &buffer[0], pack_size, MPI_BYTE, MPI_STATUS_IGNORE));
    }else{
      MPI_CHECK(MPI_File_read_at(fh, moffset, &buffer[0], pack_size, MPI_BYTE, MPI_STATUS_IGNORE));
    }
    //unpack
    sequence tmp_seq = sequence(dom) ;
    qrep->unpack(&buffer[0], loc_unpack, pack_size, tmp_seq) ;
  
#ifdef io_performance
    MPI_Barrier(MPI_COMM_WORLD);
    double wall_time = s.stop();
    if(prank == 0) std::cerr << "                                                    mpiio parallel time to read_store: "  << wall_time << endl; 
#endif
  }

  void pmpi_read_container_redistribute(std::string& filename,
                                        Loci::storeRepP var, Loci::entitySet read_set,
                                        Loci::fact_db &facts, int xfer_type) {
    
    
    
    
    if(var->RepType() == Loci::PARAMETER) {
      if(MPI_rank==0) cerr << " mpi container io doesn't provide parameter io" << endl;
      return ;
    }

    MPI_File fh;
    MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                            MPI_MODE_RDONLY, Loci::PHDF5_MPI_Info, &fh));
   
    // Read in store in file numbering
    int offset = 0 ;
    Loci::storeRepP new_store = var->new_store(EMPTY) ;
    pmpi_read_store( fh, new_store,offset,MPI_COMM_WORLD, xfer_type) ;
    MPI_File_close(&fh);

    
    // map from file number to local numbering
    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist != 0) {
      // Correct offset if file numbering changes.  Assume read_set is being
      // read in over the same set
      int minID = offset ;
      MPI_Bcast(&minID,1,MPI_INT,0,MPI_COMM_WORLD) ;
      const int minIDf = getMinFileNumberFromLocal(read_set,dist) ;
      const int correct = minIDf - minID ;
      offset += correct  ;

      // Allocate space for reordered container
      Loci::storeRepP result = var->new_store(read_set) ;
      Loci::File2LocalOrder(result,read_set,new_store,offset,dist,MPI_COMM_WORLD) ;
      // Copy results into container
      if(read_set == EMPTY) {
        read_set = result->domain() ;
        var->allocate(read_set) ;
      }
      var->copy(result,read_set) ;
    } else {
      if(read_set == EMPTY) {
        read_set = new_store->domain() ;
        var->allocate(read_set) ;
      } else {
        offset = read_set.Min() - new_store->domain().Min() ;
        if(offset != 0) {
          // shift new store by offset to correct alignment
          new_store->shift(offset) ;
        }
      }
      var->copy(new_store,read_set) ;
    }

  }
}


#endif






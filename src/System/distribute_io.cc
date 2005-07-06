#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;

#include <algorithm>
using std::sort;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>

namespace Loci {
  
  //This is a generalized routine for writing out storeRepP's. Has
  //been tested for stores, storeVecs and multiStores. The initial
  //store in the local numbering is first redistributed such that it
  //ends up in a blocked partitioning format in the global numbering
  //across all the processors. The qrep passed to this routine is in
  //the chunked partitioning format in global numbering.   
  void write_container(hid_t group_id, storeRepP qrep) {
    if(qrep->RepType() == PARAMETER) {
      if(MPI_rank==0) {
        frame_info fi = qrep->write_frame_info(group_id) ;
        int array_size = 0 ;
        if(fi.size) 
          if(fi.is_stat)   
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              array_size += *vi ;
          else
            array_size = fi.size  ;
        else
          if(fi.is_stat) 
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              array_size += *vi ;
          else
            for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
              array_size += *fvi ;

        if(array_size == 0)
          array_size = 1 ;
        hsize_t dimension = array_size ;
        int rank = 1 ;
        hssize_t start = 0 ; 
        hsize_t stride = 1 ;
        hsize_t count = array_size ;

        hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
        DatatypeP dp = qrep->getType() ;
        hid_t datatype = dp->get_hdf5_type() ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        dimension = count ;
        start += dimension ;
        hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace, H5P_DEFAULT) ;
        entitySet dom = qrep->domain() ;
        qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
      return ;
    }
    entitySet dom = qrep->domain() ;
    std::vector<entitySet> dom_vector = all_collect_vectors(dom);
    entitySet q_dom;
    for(int i = 0; i < Loci::MPI_processes; i++) 
      q_dom += dom_vector[i];

    unsigned char* tmp_send_buf ;
    std::vector<int> sort_max ;
    int local_size = qrep->pack_size(dom) ;
    sort_max = all_collect_sizes(local_size) ;
    int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
    tmp_send_buf = new unsigned char[total_size] ;
    if(Loci::MPI_rank == 0)
      Loci::HDF5_WriteDomain(group_id, q_dom);
    frame_info fi = qrep->write_frame_info(group_id) ;
    if(q_dom == EMPTY)
      return ;
    int array_size = 0 ;
    if(fi.size) 
      if(fi.is_stat)   
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	array_size = fi.size * dom.size() ;
    else
      if(fi.is_stat) 
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	  array_size += *fvi ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size) ;
    int tot_arr_size = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_arr_size += arr_sizes[i] ;
    if(Loci::MPI_rank != 0) {
      MPI_Status status ;
      int send_size_buf ;
      send_size_buf = qrep->pack_size(dom) ;
      int tot_size = send_size_buf ;
      int loc_pack = 0 ;
      qrep->pack(tmp_send_buf, loc_pack, total_size, dom) ;
      int flag = 0 ;
      MPI_Recv(&flag,1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status) ;
      if(flag) {
	MPI_Send(&tot_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
	MPI_Send(tmp_send_buf, tot_size, MPI_PACKED, 0, 12, MPI_COMM_WORLD) ;
      }
    } else {
      int rank = 1 ;
      hsize_t dimension = 1 ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = arr_sizes[0] ;
      dimension =  tot_arr_size ;
      if(dimension != 0) {
	hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	dimension = count ;
	start += dimension ;
	hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace, H5P_DEFAULT) ;
	qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
	H5Dclose(dataset) ;
	
	for(int i = 1; i < Loci::MPI_processes; ++i) {
	  MPI_Status status ;
	  int recv_total_size ;
	  entitySet tmpset = dom_vector[i];
	  Loci::storeRepP t_qrep = qrep->new_store(tmpset) ;
	  int loc_unpack = 0 ;
	  int flag = 1 ;
	  MPI_Send(&flag, 1, MPI_INT, i, 10, MPI_COMM_WORLD) ;
	  MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, MPI_COMM_WORLD, &status) ;
	  MPI_Recv(tmp_send_buf, recv_total_size, MPI_PACKED, i, 12, MPI_COMM_WORLD, &status) ;
	  Loci::sequence tmp_seq = Loci::sequence(tmpset) ;
	  t_qrep->unpack(tmp_send_buf, loc_unpack, total_size, tmp_seq) ;
	  dimension = arr_sizes[i] ;
	  count = dimension ; 
          
	  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
	  start += count ;
	  dataset = H5Dopen(group_id, "data") ;
	  t_qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", tmpset) ;
	  H5Dclose(dataset) ;
	}
	H5Sclose(dataspace) ;
      }
    }
    delete [] tmp_send_buf ;
  }

  //This routine reads storeReps from .hdf5 file. If dom specified is 
  //EMPTY then it is found by dividing total entities by processors.
  //dom entities are allocated to qrep on local processor.
  void read_container(hid_t group_id, storeRepP qrep, entitySet &dom) {
    if(qrep->RepType() == PARAMETER) {
      int pack_size = 0 ;
      if(MPI_rank==0) {
        int array_size = 0 ;
        frame_info fi = qrep->read_frame_info(group_id) ;
        if(fi.size)
          if(fi.is_stat) {
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              array_size += *vi ;
          } else {
            if(fi.size > 1)
              qrep->set_elem_size(fi.size) ;
            array_size = fi.size ;
          }
        else { 
          if(fi.is_stat) {
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi) 
              array_size += *vi ;
          }
          else {
            for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi) 
              array_size += *fvi ;
          }
        }

        hid_t dimension = array_size ;
        hid_t dataset =  H5Dopen(group_id, "data") ;
        hid_t dataspace = H5Dget_space(dataset) ;

        qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;

        pack_size = qrep->pack_size(dom) ;

      }

      // Now broadcast the result to other processors
      if(MPI_processes > 1) {
        MPI_Bcast(&pack_size,1,MPI_INT,0,MPI_COMM_WORLD) ;
        unsigned char *pack_buf = new unsigned char[pack_size] ;
        if(MPI_rank == 0) {
          int loc_pack = 0 ;
          int sz = pack_size ;
          qrep->pack(pack_buf,loc_pack,sz,dom) ;
        }
        MPI_Bcast(pack_buf,pack_size,MPI_PACKED,0,MPI_COMM_WORLD) ;
        if(MPI_rank != 0) {
          int loc_pack = 0 ;
          int sz = pack_size ;
          qrep->unpack(pack_buf,loc_pack,sz,dom) ;
        }
      }


      return ;
    }

    entitySet q_dom ;
    if(Loci::MPI_rank == 0)
      Loci::HDF5_ReadDomain(group_id, q_dom) ;
    unsigned char* tmp_buf = 0;
    std::vector<int> interval_sizes ;
    q_dom = all_collect_entitySet(q_dom) ;
    entitySet global_dom = all_collect_entitySet(dom);
    if(global_dom == EMPTY) {
      if(q_dom != EMPTY) {
	int min = q_dom.Min() ;
	int max = q_dom.Max() ;
	int sz = (max - min + 1)/Loci::MPI_processes;
	for(int i = 0; i < MPI_processes; i++) {
	  entitySet tmpEnt;
	  if(i < MPI_processes - 1) {
	    if(sz != 0)
	      tmpEnt = interval(min + i*sz, min + (i+1)*sz-1);
	  }
	  else 
	    tmpEnt = interval(min + i*sz, max);
	  tmpEnt &= q_dom;
	  interval_sizes.push_back(tmpEnt.size());
	  if(i == MPI_rank)
	    dom = tmpEnt;
	}
      }
    }
    else {
      interval_sizes = all_collect_sizes(dom.size()) ;
      entitySet tset = all_collect_entitySet(dom) ;
      if(tset != q_dom) {
	cerr << "The total domain of the container and the sum of domains across the processors does not match" << endl ;
	cerr << "q_dom = " << q_dom << endl ;
	cerr << "tset = " << tset << endl ;
      }
    }
    if(qrep->domain() == EMPTY)
      qrep->allocate(dom) ;
    FATAL(qrep->domain() != dom) ;
    frame_info fi = qrep->read_frame_info(group_id) ;
    int array_size = 0 ;
    int vec_size = 0 ;
    if(fi.size)
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
	vec_size = fi.second_level.size() ;
      }
      else {
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	array_size = fi.size * dom.size() ;
      }
    else { 
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi) 
	  array_size += *vi ;
	vec_size = fi.second_level.size() + dom.size() ;
      }
      else {
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi) 
	  array_size += *fvi ;
	vec_size = dom.size() ;
      }
    }
    std::vector<int> tmp_sizes = all_collect_sizes(vec_size) ;
    int max_tmp_size = *std::max_element(tmp_sizes.begin(), tmp_sizes.end()) ;
    int max_eset_size = *std::max_element(interval_sizes.begin(), interval_sizes.end()) ;
    int* tmp_int  ;
    tmp_int = new int[max_tmp_size] ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size) ;
    int tot_arr_size = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_arr_size += arr_sizes[i] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      int t = 0 ;
      if(fi.size) { 
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	if(fi.is_stat) 
	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	    tmp_int[t++] = *vi ;
      }
      else {
	if(fi.is_stat) {
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;
	  
	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi) 
	    tmp_int[t++] = *vi ;
	}
	else
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;
      }
      if(tmp_sizes[Loci::MPI_rank])
	MPI_Send(tmp_int, tmp_sizes[Loci::MPI_rank], MPI_INT, 0, 10, MPI_COMM_WORLD) ;
      int total_size = 0 ;
      MPI_Recv(&total_size, 1, MPI_INT, 0, 11,
	       MPI_COMM_WORLD, &status) ;  
      tmp_buf = new unsigned char[total_size] ;
      MPI_Recv(tmp_buf, total_size, MPI_PACKED, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      Loci::sequence tmp_seq = Loci::sequence(dom) ;
      int loc_unpack = 0 ;
      if(qrep->domain() == EMPTY)
	qrep->allocate(dom) ;
      FATAL(dom != qrep->domain()) ;
      qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
    } else {
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int curr_indx = 0 ;
      int total_size = 0 ;
      int tmp_total_size = 0 ;
      entitySet max_set;
      if(max_eset_size > 0)
	max_set = interval(0, max_eset_size-1) ;
      storeRepP tmp_sp ;
      if(fi.size) 
	tmp_sp = qrep->new_store(max_set) ; 
      
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	entitySet local_set;
	if(interval_sizes[p] > 0) 
	  local_set = entitySet(interval(curr_indx, interval_sizes[p]+curr_indx-1)) ;
	curr_indx += interval_sizes[p] ;
	hsize_t dimension = arr_sizes[p] ;
	count = dimension ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	entitySet tmp_set;
	if(local_set.size())
	  tmp_set = interval(0, local_set.size()-1) ;
	if(p && tmp_sizes[p]) {
	  MPI_Recv(tmp_int, tmp_sizes[p], MPI_INT, p, 10,
		   MPI_COMM_WORLD, &status) ;
	  std::vector<int> vint, fvint ;
	  int t = 0 ;
	  if(fi.size) {
	    if(fi.is_stat) {
	      for(int i = 0; i < tmp_sizes[p]; ++i)
		vint.push_back(tmp_int[t++]) ;
	      fi.second_level = vint ;
	    } 
	  }
	  else {
	    for(int i = 0; i < local_set.size(); ++i)
	      fvint.push_back(tmp_int[t++]) ;
	    for(int i = 0; i < tmp_sizes[p]-local_set.size(); ++i)
	      vint.push_back(tmp_int[t++]) ;
	    fi.first_level = fvint ;
	    fi.second_level = vint ;
	  }
	}
	storeRepP t_sp ;
	int t = 0 ;
	if(p == 0) 
	  if(!fi.size)
	    for(std::vector<int>::const_iterator vi = fi.first_level.begin(); vi != fi.first_level.end(); ++vi)
	      tmp_int[t++] = *vi ;
	if(fi.size) {
	  tmp_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ; 
	  tmp_total_size = tmp_sp->pack_size(tmp_set) ;
	}
	else {
	  t_sp = qrep->new_store(tmp_set, tmp_int) ;
	  t_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ; 
	  tmp_total_size = t_sp->pack_size(tmp_set) ;
	}
	if(tmp_total_size > total_size) {
	  total_size = tmp_total_size ;
	  if(p)
	    delete [] tmp_buf ;
	  tmp_buf = new unsigned char[total_size] ;
	}
	start += count ;
	int loc = 0 ;
	if(fi.size)
	  tmp_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
	else
	  t_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
	if(p == 0) {
	  int loc_unpack = 0 ;
	  Loci::sequence tmp_seq = Loci::sequence(dom) ;
	  if(fi.size) 
	    if(fi.size > 1) 
	      qrep->set_elem_size(fi.size) ;
	  if(qrep->domain() == EMPTY)
	    qrep->allocate(dom) ;
	  FATAL(dom != qrep->domain()) ;
	  qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
	} else { 
	  MPI_Send(&total_size, 1, MPI_INT, p, 11, MPI_COMM_WORLD) ;
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, p, 12, MPI_COMM_WORLD) ;
	}
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
    delete [] tmp_buf ;
    delete [] tmp_int ; 
  }

  //This routine returns a storeRep which is created by using passed
  //storeRep sp and map remap.  It is used whenever a storeRep needs to 
  //be converted in global numbering to write out to a file.  
  //remap is mapping from io entities(read from grid file)
  //to loci entities(global numbering of Loci)
  storeRepP collect_reorder_store(storeRepP &sp, dMap& remap, fact_db &facts) {
    entitySet dom = sp->domain() ;
    dom =  collect_entitySet(dom) ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    entitySet loci_dom = init_ptn[MPI_rank] & dom ;
    
    entitySet global_loci_dom  =  all_collect_entitySet(dom) ;
    entitySet global_io_dom = remap.preimage(global_loci_dom).first;
    global_io_dom = all_collect_entitySet(global_io_dom);

    if(global_loci_dom.size() != global_io_dom.size()) {
      cerr << "remap is not able to provide information of all entities.";
      cerr << "loci_entities: " << global_loci_dom << " has total " << global_loci_dom.size() << endl; 
      cerr << "mapped io_entities : " << global_io_dom << " has total "  << global_io_dom.size() << endl;
    }

    dMap reverse ; //Mapping from loci entities to io entities 
    {
      dmultiMap d_remap ;
      if(facts.is_distributed_start())
        distributed_inverseMap(d_remap, remap, global_loci_dom, global_io_dom, init_ptn) ;
      else
        inverseMap(d_remap, remap, global_loci_dom, global_loci_dom) ;
      
      FORALL(loci_dom, ri) {
        if(d_remap[ri].size() == 1)
          reverse[ri] = d_remap[ri][0] ;
        else
          if(d_remap[ri].size() > 1)
            cerr << "inverse of remap has multiple entries!" << endl ;
      } ENDFORALL ;
    }
    
    Map l2g ;
    l2g = facts.get_variable("l2g") ;

    entitySet local_dom =  MapRepP(l2g.Rep())->preimage(loci_dom).first ;

    if(loci_dom.size() != local_dom.size()) {
      cerr << "Something fishy going on with l2g, which does not seem to give one to one and onto mapping.";
      cerr << "loci_dom: " << loci_dom << " has total " << loci_dom.size() << endl;
      cerr << "local_dom: " << local_dom << " has total " << local_dom.size() << endl;
    }

    Map local2io ;
    local2io.allocate(local_dom) ;
    FORALL(local_dom, ri) {
      local2io[ri] = l2g[ri] ;
    } ENDFORALL ;
    
    local_dom = MapRepP(local2io.Rep())->preimage(reverse.domain()).first ;
    MapRepP(local2io.Rep())->compose(reverse, local_dom) ;
    entitySet io_owned =  MapRepP(local2io.Rep())->image(local_dom) ;

    std::vector<entitySet> chop_ptn(Loci::MPI_processes);
    if(global_io_dom != EMPTY) {
      int min = global_io_dom.Min();
      int max = global_io_dom.Max();
      int sz = (max - min + 1)/Loci::MPI_processes;
      for(int i = 0; i < MPI_processes; i++) {
	entitySet tmpEnt;
	if(i < MPI_processes - 1) {
	  if(sz != 0)
	    tmpEnt = interval(min + i*sz, min + (i+1)*sz-1);
	}
	else 
	  tmpEnt = interval(min + i*sz, max);
	tmpEnt &= global_io_dom;
	chop_ptn[i] = tmpEnt;
      }
    }

    entitySet io_dom = chop_ptn[MPI_rank]; 
    entitySet owned_io_out_of_dom, filled_io_entities, filled_local_entities ;
    
    FORALL(local_dom, ii) {
      if(io_dom.inSet(local2io[ii])) { 
	filled_io_entities += local2io[ii] ;
	filled_local_entities += ii ;
      }
    } ENDFORALL ;

    owned_io_out_of_dom = io_owned - filled_io_entities ;

    storeRepP local2io_sp = MapRepP(local2io.Rep())->thaw() ;
    dMap d_local2io(local2io_sp) ; 

    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(io_dom) ;
    if(filled_local_entities != EMPTY)
      qcol_rep->scatter(d_local2io, sp, filled_local_entities) ;

    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    int size_send = 0 ;
    std::vector<entitySet> send_clone( MPI_processes) ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    std::vector<storeRepP> scl( MPI_processes), rcl( MPI_processes) ;
    
    for(int i = 0; i <  MPI_processes; ++i) {
      send_clone[i] = owned_io_out_of_dom & chop_ptn[i] ;
      send_dom[i] = MapRepP(local2io.Rep())->preimage(send_clone[i]).first ;
      entitySet tmp = send_dom[i];
      store<int> tmp_store ;
      entitySet tmpEnt;
      if(tmp.size())
	tmpEnt = interval(0, tmp.size() - 1);
      tmp_store.allocate(tmpEnt) ;
      send_count[i] = tmp_store.Rep()->pack_size(tmpEnt) ;
      ei = tmp.begin() ;
      for(int j = 0; j < tmp.size(); ++j) {
	tmp_store[j] = local2io[*ei] ;
	++ei ;
      }
      scl[i] = tmp_store.Rep() ;
      size_send += send_count[i];
    }
    unsigned char *send_buf = new unsigned char[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      scl[i]->pack(send_buf, loc_pack, size_send, scl[i]->domain()) ; 
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    unsigned char *recv_buf = new unsigned char[size_send] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_PACKED,
		  recv_buf, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    std::vector< sequence> recv_dom( MPI_processes) ;
    int loc_unpack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      int tmp_size = recv_count[i] / sizeof(int) ;
      entitySet tmp ;
      if(tmp_size > 0)
	tmp = interval(0, tmp_size-1) ;
      recv_dom[i] =  sequence(tmp) ;
      store<int> tmp_store ;
      tmp_store.allocate(tmp) ;
      tmp_store.Rep()->unpack(recv_buf, loc_unpack, size_send, recv_dom[i]) ; 
      rcl[i] = tmp_store.Rep() ;
    }
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    }
    
    unsigned char *send_store = new unsigned char[size_send] ;
    std::vector<storeRepP> tmp_sp( MPI_processes) ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    int size_recv = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp[i] = sp->new_store(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    }
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_unpack = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      tmp_sp[i]->unpack(recv_store, loc_unpack, size_recv, recv_dom[i]) ; 
    }
    for(int i = 0; i < MPI_processes; ++i) {	
      dMap m;
      m.allocate(rcl[i]->domain()) ;
      store<int> RCL(rcl[i]) ;
      FORALL(m.domain(), mi) {
	m[mi] = RCL[mi] ;
      } ENDFORALL ;

      if(tmp_sp[i]->domain() != EMPTY)
	qcol_rep->scatter(m, tmp_sp[i], tmp_sp[i]->domain()) ;  
    }

    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return qcol_rep ;
  } 

  //This routine fills new_sp storeRep for allocated entities using passed 
  //storeRep sp_init and map remap.  It is used whenever a storeRep needs to 
  //be converted in local numbering after it is read from a file. 
  void distribute_reorder_store(storeRepP &new_sp, storeRepP sp_init, dMap& remap, fact_db &facts) {
    entitySet our_dom = new_sp->domain() ;
    entitySet remap_image  = Loci::collect_entitySet(our_dom) ;
    remap_image = all_collect_entitySet(remap_image) ;
    Map l2g ; 
    l2g = facts.get_variable("l2g") ;
    entitySet read_set = remap.preimage(remap_image).first;
    read_set = all_collect_entitySet(read_set);

    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    entitySet local_own = d->g2l.domain() & remap.image(read_set) ;

#ifdef DEBUG
    entitySet global_own = all_collect_entitySet(local_own);
    if(global_own != remap_image) {
      cerr << "global_own should be same as remap_image." << endl;
      cerr << "g2l does not define " << remap_image - global_own << endl;
    }
#endif

    entitySet t_dom = remap.preimage(local_own).first ;
    dMap tmp_remap ; 
    FORALL(t_dom, ti) { 
      tmp_remap[ti] = remap[ti] ;
    } ENDFORALL ; 
    Loci::MapRepP(tmp_remap.Rep())->compose(d->g2l, t_dom) ;

    Loci::constraint my_entities ; 
    my_entities = d->my_entities ; 
    entitySet local_req = *my_entities & our_dom ;
    entitySet global_req = Loci::MapRepP(l2g.Rep())->image(local_req) ;
    global_req = Loci::MapRepP(remap.Rep())->preimage(global_req).first ;

    entitySet tmp_dom = sp_init->domain() ;
    std::vector<entitySet> sp_init_domain = all_collect_vectors(tmp_dom);
    entitySet global_sp_init_domain;
    for(int i = 0; i < Loci::MPI_processes; i++) {
      global_sp_init_domain += sp_init_domain[i];
    }
    
    FATAL((global_req - t_dom).size() != 0); 
    FATAL((global_req - global_sp_init_domain).size() != 0);
    
    unsigned char *tmp_buf, *recv_buf ;
    int total_size = 0 ;
    total_size = sp_init->pack_size(tmp_dom) ;
    total_size = GLOBAL_MAX(total_size);
    tmp_buf = new unsigned char[total_size] ;
    recv_buf = new unsigned char[total_size] ;

    int loc = 0 ;
    sp_init->pack(tmp_buf, loc, total_size, tmp_dom) ; 
    MPI_Status status, stat_1 ;
    if(Loci::MPI_rank != 0) {
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	MPI_Recv(recv_buf, total_size, MPI_PACKED, Loci::MPI_rank-1, 12,
		 MPI_COMM_WORLD, &status) ;  
	if(Loci::MPI_rank < Loci::MPI_processes-1) { 
	  MPI_Send(recv_buf, total_size, MPI_PACKED, Loci::MPI_rank+1, 12, MPI_COMM_WORLD) ;
	}
	if(p == 0)
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, 0, 12, MPI_COMM_WORLD) ;
	entitySet local_set = sp_init_domain[p] ;
	Loci::sequence tmp_seq = Loci::sequence(local_set) ;
	int loc_unpack = 0 ;
	storeRepP sp = sp_init->new_store(local_set) ;
	sp->unpack(recv_buf, loc_unpack, total_size, tmp_seq) ;
	local_set &= global_req ;
	new_sp->scatter(tmp_remap, sp, local_set) ; 
      }
    } else {
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	if(p > 0) 
	  MPI_Recv(recv_buf, total_size, MPI_PACKED, p, 12,
		   MPI_COMM_WORLD, &stat_1) ;
	if(p == 0)
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, 1, 12, MPI_COMM_WORLD) ;
	else
	  MPI_Send(recv_buf, total_size, MPI_PACKED, 1, 12, MPI_COMM_WORLD) ;
	entitySet local_set = sp_init_domain[p] ;  
	int loc_unpack = 0 ;
	Loci::sequence tmp_seq = Loci::sequence(local_set) ;
	storeRepP sp = sp_init->new_store(local_set) ;
	if(p != 0)
	  sp->unpack(recv_buf, loc_unpack, total_size, tmp_seq) ;
	local_set &= global_req ;
	if(p == 0)
	  new_sp->scatter(tmp_remap, sp_init, local_set) ; 
	else 
	  new_sp->scatter(tmp_remap, sp, local_set) ; 
      }
    }
    delete [] recv_buf ;
    delete [] tmp_buf ;
  }

  void redistribute_write_container(hid_t file_id, std::string vname,
                               Loci::storeRepP var, fact_db &facts) {
    hid_t group_id = 0 ;
    if(Loci::MPI_rank == 0)
      group_id = H5Gcreate(file_id, vname.c_str(), 0) ;
    if(Loci::MPI_processes == 1 || var->RepType() == PARAMETER) {
      Loci::write_container(group_id, var) ;
    } else {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      
      dMap remap = df->remap ;
      Loci::storeRepP vardist = collect_reorder_store(var,remap,facts) ;
      Loci::write_container(group_id,vardist) ;
    }
    if(Loci::MPI_rank == 0)
      H5Gclose(group_id) ;
  }

  void read_container_redistribute(hid_t file_id, std::string vname,
                                   Loci::storeRepP var, entitySet read_set,
                                   fact_db &facts) {
    hid_t group_id = 0;
    if(Loci::MPI_rank == 0)
      group_id = H5Gopen(file_id, vname.c_str()) ;
    if(var->RepType() == PARAMETER) {
      Loci::read_container(group_id, var, read_set) ;
      if(Loci::MPI_rank == 0)
        H5Gclose(group_id) ;
      return ;
    }
    if(Loci::MPI_processes == 1) {
      Loci::read_container(group_id, var, read_set) ;
    } else {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      
      dMap remap = df->remap ;

      entitySet alloc_dom = EMPTY;
      Loci::storeRepP new_store = var->new_store(EMPTY) ;
      Loci::read_container(group_id, new_store , alloc_dom) ;
      Loci::storeRepP result = var->new_store(read_set) ;
      Loci::distribute_reorder_store(result,new_store,remap,facts) ;
      var->copy(result,result->domain()&read_set) ;
    }
    if(Loci::MPI_rank == 0)
      H5Gclose(group_id) ;
    
  }

    
}

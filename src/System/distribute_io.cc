#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;

#include <algorithm>
using std::sort;
using std::pair ;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>

namespace Loci {

  entitySet BcastEntitySet(entitySet set, int root, MPI_Comm comm) {
    
      // Now lets share the domain with all other processors ;
    int sz = set.num_intervals() ;
    MPI_Bcast(&sz,1,MPI_INT,root,comm) ;
    vector<interval> vlist(sz) ;
    if(MPI_rank == 0) {
      for(int i=0;i<sz;++i)
        vlist[i] = set[i] ;
    }
    MPI_Bcast(&vlist[0],sz*2,MPI_INT,root,comm) ;
    set = EMPTY ;
    for(int i = 0;i<sz;++i)
      set += vlist[i] ;
    return set ;
  }

  vector<int> simplePartition(int mn, int mx, int p) {
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
      
  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    vector<int> pl = simplePartition(mn,mx,p) ;
    vector<entitySet> ptn(p) ;
    for(int i=0;i<p;++i)
      ptn[i] = interval(pl[i],pl[i+1]-1) ;
    return ptn ;
  }

      
     

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
#ifdef H5_INTERFACE_1_6_4
        hsize_t start = 0 ;
#else
        hssize_t start = 0 ;
#endif
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
    for(int i = 0; i < MPI_processes; i++) 
      q_dom += dom_vector[i];

    unsigned char* tmp_send_buf ;
    std::vector<int> sort_max ;
    int local_size = qrep->pack_size(dom) ;
    sort_max = all_collect_sizes(local_size) ;
    int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
    tmp_send_buf = new unsigned char[total_size] ;
    if(MPI_rank == 0)
      HDF5_WriteDomain(group_id, q_dom);
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

    if(MPI_rank != 0) {
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
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
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

	for(int i = 1; i < MPI_processes; ++i) {
	  MPI_Status status ;
	  int recv_total_size ;
	  entitySet tmpset = dom_vector[i];

	  storeRepP t_qrep = qrep->new_store(tmpset) ;

	  int loc_unpack = 0 ;
	  int flag = 1 ;
	  MPI_Send(&flag, 1, MPI_INT, i, 10, MPI_COMM_WORLD) ;
	  MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, MPI_COMM_WORLD, &status) ;
	  MPI_Recv(tmp_send_buf, recv_total_size, MPI_PACKED, i, 12, MPI_COMM_WORLD, &status) ;

	  sequence tmp_seq = sequence(tmpset) ;
	  t_qrep->unpack(tmp_send_buf, loc_unpack, total_size, tmp_seq) ;
	  dimension = arr_sizes[i] ;
	  count = dimension ; 

	  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
	  start += count ;

	  dataset = H5Dopen(group_id, "data") ;
	  t_qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", tmpset) ;
          t_qrep->allocate(EMPTY) ;

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
        delete[] pack_buf ;
      }


      return ;
    }

    // Here we read in a store container.  First lets read in the domain
    entitySet q_dom ;
    if(MPI_rank == 0)
      HDF5_ReadDomain(group_id, q_dom) ;

    // Now lets share the domain with all other processors ;
    q_dom = BcastEntitySet(q_dom,0,MPI_COMM_WORLD) ;

    unsigned char* tmp_buf = 0;
    std::vector<int> interval_sizes ;

    int sz = dom.size() ;
    int total_sz = 0 ;
    MPI_Allreduce(&sz,&total_sz,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(total_sz == 0) {
      if(q_dom != EMPTY) {
        vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),
                                                MPI_COMM_WORLD) ;
        for(int i=0;i<MPI_processes;++i) {
          entitySet qset = ptn[i] &q_dom ;
          interval_sizes.push_back(qset.size()) ;
          if(i == MPI_rank) 
            dom = qset ;
        }
      } else
        for(int i=0;i<MPI_processes;++i) 
          interval_sizes.push_back(0) ;
    } else {
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
    if(MPI_rank != 0) {
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
      if(tmp_sizes[MPI_rank])
	MPI_Send(tmp_int, tmp_sizes[MPI_rank], MPI_INT, 0, 10, MPI_COMM_WORLD) ;
      int total_size = 0 ;
      MPI_Recv(&total_size, 1, MPI_INT, 0, 11,
	       MPI_COMM_WORLD, &status) ;  
      tmp_buf = new unsigned char[total_size] ;
      MPI_Recv(tmp_buf, total_size, MPI_PACKED, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      sequence tmp_seq = sequence(dom) ;
      int loc_unpack = 0 ;
      if(qrep->domain() == EMPTY)
	qrep->allocate(dom) ;
      FATAL(dom != qrep->domain()) ;
      qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
    } else {
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
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
      
      for(int p = 0; p < MPI_processes; ++p) {
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
	  sequence tmp_seq = sequence(dom) ;
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

  entitySet findBoundingSet(entitySet in) {
    Entity max_val = in.Max() ;
    Entity min_val = in.Min() ;
    Entity gmin_val = min_val;
    Entity gmax_val = max_val ;
    MPI_Allreduce(&min_val,&gmin_val,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
    MPI_Allreduce(&max_val,&gmax_val,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

    return entitySet(interval(gmin_val,gmax_val)) ;
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
  
  //This routine returns a storeRep which is created by using passed
  //storeRep sp and map remap.  It is used whenever a storeRep needs to 
  //be converted in global numbering to write out to a file.  
  //remap is mapping from io entities(read from grid file)
  //to loci entities(global numbering of Loci)
  storeRepP collect_reorder_store(storeRepP &sp, fact_db &facts) {

    entitySet dom = sp->domain() ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    constraint my_entities ; 
    my_entities = d->my_entities ;
    dom = *my_entities & dom ;


    Map l2g ;
    l2g = d->l2g.Rep() ;
    MapRepP l2gP = MapRepP(l2g.Rep()) ;

    entitySet dom_global = l2gP->image(dom) ;

    FATAL(dom.size() != dom_global.size()) ;

    dMap g2f ;
    g2f = d->g2f.Rep() ;
    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[l2g[i]] ;
      debugout << "newnum["<< i << "] = " << newnum[i] << endl ;
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;

    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;
    
    int p = MPI_processes ;

    
    vector<entitySet> out_ptn = simplePartition(imn,imx,MPI_COMM_WORLD) ;

    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;
    for(int i=0;i<p;++i) {
      send_sets[i] = newnum.preimage(out_ptn[i]).first ;
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;

    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    debugout << "file_dom = " << file_dom << endl ;
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;


    int *send_sizes = new int[p] ;
    int *recv_sizes = new int[p] ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

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
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(send_store, &send_sizes[0], send_dspl, MPI_PACKED,
		  recv_store, &recv_sizes[0], recv_dspl, MPI_PACKED,
		  MPI_COMM_WORLD) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    delete[] recv_store ;
    delete[] send_store ;
    delete[] recv_dspl ;
    delete[] send_dspl ;
    delete[] recv_sizes ;
    delete[] send_sizes ;

    return qcol_rep ;
  } 

  
  // result allocated over local numbering 
  // input is in file numbering
  // Routine redistributes input to be arranged in the local numbering
  // and placed in result.
  // input is assumed to be partitioned by simplePartion
  
  void distribute_reorder_store(storeRepP &result, entitySet resultSet,
                                storeRepP input, fact_db &facts ) {
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
    vector<int> fptn = simplePartition(imn,imx,p) ;

    // Get distribution plan
    vector<vector<pair<int,int> > > dist_plan(p) ;
    FORALL(resultSet,i) {
      int fn = newnum[i] ;
      if(fn < imn || fn > imx) {
        cerr << "Problem with distribute_reorder_store, index out of bounds"
             << endl ;
        cerr << "fn = " << fn << "imx = " << imx << "imn = " << imn << endl ;
        Loci::Abort() ;
      }
      // processor that contains this value
      int r = p*(fn-imn)/(imx-imn+1) ;
      for(;;) { // refine
        if(fn >= fptn[r] && fn < fptn[r+1])
          break ;
        r+= (fn < fptn[r])?-1:1 ;
        FATAL(r >= p) ;
        FATAL(r < 0) ;
      }
      dist_plan[r].push_back(pair<int,int>(fn,i)) ;
    } ENDFORALL ;

    // Compute recv requests from distribution plan
    vector<sequence> recv_seq(p),send_req(p) ;
    for(int i=0;i<p;++i) {
      sort(dist_plan[i].begin(),dist_plan[i].end()) ;
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

  void redistribute_write_container(hid_t file_id, std::string vname,
                               storeRepP var, fact_db &facts) {
    hid_t group_id = 0 ;
    if(MPI_rank == 0)
      group_id = H5Gcreate(file_id, vname.c_str(), 0) ;
    if(MPI_processes == 1 || var->RepType() == PARAMETER) {
      write_container(group_id, var) ;
    } else {
      storeRepP vardist = collect_reorder_store(var,facts) ;
      write_container(group_id,vardist) ;
    }
    if(MPI_rank == 0)
      H5Gclose(group_id) ;
  }

  void read_container_redistribute(hid_t file_id, std::string vname,
                                   storeRepP var, entitySet read_set,
                                   fact_db &facts) {
    hid_t group_id = 0;
    if(MPI_rank == 0)
      group_id = H5Gopen(file_id, vname.c_str()) ;
    if(var->RepType() == PARAMETER) {
      read_container(group_id, var, read_set) ;
      if(MPI_rank == 0)
        H5Gclose(group_id) ;
      return ;
    }
    if(MPI_processes == 1) {
      read_container(group_id, var, read_set) ;
    } else {
      entitySet alloc_dom = EMPTY;
      storeRepP new_store = var->new_store(EMPTY) ;
      read_container(group_id, new_store , alloc_dom) ;
      storeRepP result = var->new_store(read_set) ;
      distribute_reorder_store(result,read_set,new_store,facts) ;
      var->copy(result,read_set) ;
    }
    if(MPI_rank == 0)
      H5Gclose(group_id) ;
    
  }

  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) {
    vector<int> ids(local_set.size()) ;
      
      
    int c = 0 ;
    if(MPI_processes > 1) {
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      FORALL(local_set,ii) {
        ids[c++] = l2g[ii] ;
      } ENDFORALL ;
    } else {
      FORALL(local_set,ii) {
        ids[c++] = ii ;
      } ENDFORALL ;
    }
    writeUnorderedVector(file_id,"entityIds",ids) ;
  }

  hid_t createUnorderedFile(const char * filename, entitySet set, fact_db &facts) {
    hid_t file_id = 0;
    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
      group_id = H5Gcreate(file_id,"dataInfo",0) ;
    }
    writeSetIds(group_id,set,facts) ;
    if(MPI_rank == 0) 
      H5Gclose(group_id) ;
    return file_id ;
  }

  void closeUnorderedFile(hid_t file_id) {
    if(MPI_rank == 0)
      H5Fclose(file_id) ;
  }
    
}

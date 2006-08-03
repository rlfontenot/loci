#include <vector>
using std::vector;
#include <set>
using std::set;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>

namespace Loci {

  // Pass in a set of entitySets one for each processor to broadcast
  // to other processors.  Return with a vector of entitySets as sent to
  // you by each other processor.
  vector<entitySet> Alltoall_entitySet(vector<entitySet> v) {
    WARN(int(v.size()) != MPI_processes) ;

    if(MPI_processes == 1)
      return v ;

    const int p = v.size() ;
    vector<int> ivals(p) ;
    for(int i=0;i<p;++i)
      ivals[i] = v[i].num_intervals() ;
    vector<int> ivalr(p) ;
    MPI_Alltoall(&(ivals[0]),1,MPI_INT,&(ivalr[0]),1,MPI_INT,MPI_COMM_WORLD) ;
    vector<int> sdispls(p),rdispls(p) ;
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+ivals[i-1]*2 ;
      rdispls[i] = rdispls[i-1]+ivalr[i-1]*2 ;
    }
    vector<int> sbuf(sdispls[p-1]+ivals[p-1]*2) ;
    vector<int> rbuf(rdispls[p-1]+ivalr[p-1]*2) ;
    vector<int> scounts(p),rcounts(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = 2*ivals[i] ;
      rcounts[i] = 2*ivalr[i] ;
      for(int j=0;j<ivals[i];++j) {
        sbuf[sdispls[i]+j*2] = v[i][j].first ;
        sbuf[sdispls[i]+j*2+1] = v[i][j].second ;
      }
    }
    MPI_Alltoallv(&(sbuf[0]),&(scounts[0]),&(sdispls[0]),MPI_INT,
                  &(rbuf[0]),&(rcounts[0]),&(rdispls[0]),MPI_INT,
                  MPI_COMM_WORLD) ;

    vector<entitySet> retv(p) ;
    for(int i=0;i<p;++i) {
      for(int j=0;j<ivalr[i];++j) {
        retv[i] += interval(rbuf[rdispls[i]+j*2],rbuf[rdispls[i]+j*2+1]) ;
      }
    }
    return retv ;
  }

  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) {
    if(MPI_processes == 1)
      return m ;

    const int p = MPI_processes ;
    entitySet dom = m.domain() ;
    std::vector<entitySet> send_slices(p) ;
    for(int i=0;i<p;++i) {
      send_slices[i] = dom & init_ptn[i] ;
      dom -= init_ptn[i] ;
    }
    WARN(dom != EMPTY) ;

    vector<entitySet> recv_slices = Alltoall_entitySet(send_slices) ;

    vector<int> scounts(p),rcounts(p), sdispls(p),rdispls(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = send_slices[i].size() ;
      rcounts[i] = recv_slices[i].size() ;
    }
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    vector<int> sbuf(sdispls[p-1]+scounts[p-1]) ;
    vector<int> rbuf(rdispls[p-1]+rcounts[p-1]) ;
    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = send_slices[i].begin();
          ei != send_slices[i].end();++ei,++k) {
        sbuf[sdispls[i]+k] = m[*ei] ;
      }
    }

    MPI_Alltoallv(&(sbuf[0]),&(scounts[0]),&(sdispls[0]),MPI_INT,
                  &(rbuf[0]),&(rcounts[0]),&(rdispls[0]),MPI_INT,
                  MPI_COMM_WORLD) ;

    dMap ret_map ;

    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = recv_slices[i].begin();
          ei != recv_slices[i].end();++ei,++k) {
        ret_map[*ei] = rbuf[rdispls[i]+k] ;
      }
    }

    return ret_map ;

  }

  storeRepP distribute_store(storeRepP &sp, fact_db &facts) {
    if(!facts.isDistributed()) {
      return sp ;
    }
    storeRepP nsp ;
    if(facts.isDistributed()) {
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      d = facts.get_distribute_info() ;
      l2g = facts.get_variable("l2g") ;
      if(d->myid == 0) {
	std::vector<entitySet> ent(MPI_processes-1) ;
	MPI_Status *status, *size_status, *store_status;
	MPI_Request *recv_request, *size_request, *store_request ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	sequence te ;
	entitySet temp, me ;
	entitySet my ;

	for(ti = d->my_entities.begin(); ti != d->my_entities.end(); ++ti) {
	  temp += l2g[*ti] ;
	}
	me = temp & sp->domain() ;
	for(ti = me.begin(); ti != me.end(); ++ti) {
	  te += d->g2l[*ti] ;
	  my += d->g2l[*ti] ;
	}

	recv_size = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,1, MPI_COMM_WORLD, &size_request[k]);
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	recv_buffer = new int*[MPI_processes-1] ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  recv_buffer[i] = new int[recv_size[i]] ;
	}
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;

	for(k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,2, MPI_COMM_WORLD, &recv_request[k] );

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(k = 0; k < MPI_processes-1; ++k) {
	  entitySet re ;
	  for(int i = 0 ; i < recv_size[k]; ++i)
	    re += recv_buffer[k][i] ;

	  ent[k] = re ;
	}
	nsp = sp->new_store(my) ;
	//nsp = sp->new_store(EMPTY) ;
        //nsp->allocate(my) ;

	int sz = 0 ;
	int my_sz = sp->pack_size(me) ;
	int my_pack = 0 ;
	int my_unpack = 0 ;
	int *s_size = new int[MPI_processes-1] ;

	for(int i = 0; i < MPI_processes-1; ++i) {
	  s_size[i] = sp->pack_size(ent[i]) ;
	  sz += sp->pack_size(ent[i]) ;
	}
	unsigned char **send_ptr = new unsigned char*[MPI_processes-1] ;
	unsigned char* my_stuff = new unsigned char[my_sz] ;
	store_request = new MPI_Request[MPI_processes-1] ;
	store_status = new MPI_Status[MPI_processes-1] ;
	sp->pack(my_stuff, my_pack,my_sz,me) ;
	send_ptr[0] = new unsigned char[sz] ;
	for(int i = 1; i < MPI_processes-1; i++)
	  send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
	for(int i = 0; i < MPI_processes-1; i++) {
	  int loc_pack = 0 ;
	  sp->pack(send_ptr[i], loc_pack, s_size[i], ent[i]) ;
	  MPI_Isend(send_ptr[i], s_size[i], MPI_PACKED, i+1, 3,
		    MPI_COMM_WORLD, &store_request[i]) ;
	}
#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, store_request, store_status) ;
	FATAL(err != MPI_SUCCESS) ;
	nsp->unpack(my_stuff, my_unpack, my_sz, te) ;
	delete [] recv_size ;
	delete [] recv_buffer ;
	delete [] status ;
	delete [] recv_request ;
	delete [] size_request ;
	delete [] size_status ;
	delete [] store_request ;
	delete [] store_status ;
	delete [] send_ptr ;

      }
      else {
	int *send_buffer ;
	int send_size ;
	int loc_unpack = 0 ;
	unsigned char *recv_ptr ;
	MPI_Status stat ;
	entitySet re ;
	sequence tempseq ;
	entitySet my = sp->domain() & d->my_entities ;
	send_size = my.size() ;
	send_buffer = new int[send_size] ;
	int sz = sp->pack_size(my) ;

	recv_ptr = new unsigned char[sz] ;
	nsp = sp->new_store(EMPTY) ;
        nsp->allocate(my) ;
	int j = 0 ;

	for(ti = my.begin(); ti != my.end(); ++ti) {
	  re += l2g[*ti] ;
	}
	for(ti = re.begin(); ti != re.end(); ++ti) {
	  send_buffer[j] = *ti  ;
	  ++j ;
	}
	for(ti = re.begin(); ti != re.end(); ++ti)
	  tempseq += d->g2l[*ti] ;

	MPI_Send(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD) ;

	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD) ;
	MPI_Recv(recv_ptr, sz, MPI_PACKED, 0, 3, MPI_COMM_WORLD, &stat) ;
	nsp->unpack(recv_ptr, loc_unpack, sz,tempseq) ;
	delete [] send_buffer ;
	delete [] recv_ptr ;

      }
    }
    return nsp ;

  }
  namespace {
    inline bool fieldSort2(const std::pair<Entity,Entity> &p1,
                          const std::pair<Entity,Entity> &p2) {
      return p1.second < p2.second ;
    }
  }

  void distributed_inverseMap(multiMap &result,
                              vector<pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn) {
    // Sort input according to second field
    sort(input.begin(),input.end(),fieldSort2) ;

    // Now count what we will be sending
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = 0 ;
    int current_p = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;

      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_sz[current_p] += 2 ;
    }

    // transfer recv sizes
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

    current_p = 0 ;
    vector<int> offsets(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      offsets[i] = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;
      int from = input[i].first ;
      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_store[send_displacement[current_p]+offsets[current_p]++] = to ;
      send_store[send_displacement[current_p]+offsets[current_p]++] = from ;
    }
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  

    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    store<int> sizes ;

    sizes.allocate(local_input_image) ;
    FORALL(local_input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      sizes[indx]++ ;
    }
    result.allocate(sizes) ;

    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      int refer = recv_store[i] ;
      //      if(input_image.inSet(indx)) {
        sizes[indx]-- ;
        result[indx][sizes[indx]] = refer ;
        //      }
    }
#ifdef DEBUG
    FORALL(local_input_image,i) {
      WARN(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
    
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;
  }

  void distributed_inverseMap(dmultiMap &result,
                              vector<pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn) {
    // Sort input according to second field
    sort(input.begin(),input.end(),fieldSort2) ;

    // Now count what we will be sending
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = 0 ;
    int current_p = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;

      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_sz[current_p] += 2 ;
    }

    // transfer recv sizes
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

    current_p = 0 ;
    vector<int> offsets(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      offsets[i] = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;
      int from = input[i].first ;
      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_store[send_displacement[current_p]+offsets[current_p]++] = to ;
      send_store[send_displacement[current_p]+offsets[current_p]++] = from ;
    }
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  

    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    store<int> sizes ;
    sizes.allocate(local_input_image) ;
    FORALL(local_input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      //      if(input_image.inSet(indx))
        sizes[indx]++ ;
    }
    
    FORALL(local_input_image,i) {
      vector<int,malloc_alloc<int> > tmp(sizes[i]) ;
      tmp.swap(result[i]) ;
    }ENDFORALL ;

    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      int refer = recv_store[i] ;
      sizes[indx]-- ;
      result[indx][sizes[indx]] = refer ;
    }
#ifdef DEBUG
    FORALL(local_input_image,i) {
      WARN(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
    
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;
  }
 
  void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      for(size_t k=0;k<input_map[i].size();++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  
  
  void distributed_inverseMap(dmultiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int sz = input_map.end(i)-input_map.begin(i) ;
      for(int k=0;k<sz;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }
  // non dynamic multiMap
  void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(multiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      for(size_t k=0;k<input_map[i].size();++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  
  
  void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int sz = input_map.end(i)-input_map.begin(i) ;
      for(int k=0;k<sz;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }
  
  /* This routine collects the store variables into processor
     0. Finally we have a single store allocated over the entities in
     global numbering in processor 0 */
  storeRepP collect_store(storeRepP &sp, fact_db &facts) {
    storeRepP nsp = sp ;
    if(facts.isDistributed())  {
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      if(d->myid == 0) {
	std::vector<sequence> vseq(MPI_processes-1) ;
	MPI_Status *status, *size_status, *store_status ;
	MPI_Request *recv_request, *size_request, *store_request ;
	int **recv_buffer ;
	int *recv_size ;
	entitySet re ;
	sequence te ;
	entitySet temp = sp->domain() & d->my_entities ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  te += l2g[*ti] ;
	re += entitySet(te) ;
	recv_size = new int[MPI_processes-1] ;
	int *recv_size_bytes = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(int k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,11, MPI_COMM_WORLD, &size_request[k]);

#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;

	recv_buffer = new int*[MPI_processes-1] ;
        int recv_size_total = 0 ;
        for(int k=0;k<MPI_processes-1;++k) {
	  recv_size_total += recv_size[k] ;
	}
	recv_buffer[0] = new int[recv_size_total] ;
	for(int i = 1; i < MPI_processes-1; ++i)
	  recv_buffer[i] = recv_buffer[i-1] + recv_size[i-1] ;
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;

	for(int k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,12, MPI_COMM_WORLD, &recv_request[k] );

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(int k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_size_bytes[k],1,MPI_INT, k+1,14,
		    MPI_COMM_WORLD, &size_request[k]) ;

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;

	for(int k = 0; k < MPI_processes-1; ++k) {
	  sequence tempseq ;
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    re += recv_buffer[k][i] ;
	    tempseq += recv_buffer[k][i] ;
	  }
	  vseq[k] = tempseq ;
	}
	int my_sz= sp->pack_size(temp) ;
	int my_unpack = 0 ;
	int loc_unpack = 0 ;
	int loc_pack = 0 ;
	int *r_size = new int[MPI_processes-1] ;
	store_request = new MPI_Request[MPI_processes-1] ;
	store_status = new MPI_Status[MPI_processes-1] ;

	int sz = 0 ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  r_size[i] = recv_size_bytes[i] ;
	  sz += r_size[i] ;
	}

	unsigned char **recv_ptr = new unsigned char*[MPI_processes-1] ;
	unsigned char* my_stuff = new unsigned char[my_sz] ;
	sp->pack(my_stuff, loc_pack,my_sz,temp) ;
	nsp = sp->new_store(re) ;
	recv_ptr[0] = new unsigned char[sz] ;
	for(int i = 1; i < MPI_processes-1; i++)
	  recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;

	for(int i = 0; i < MPI_processes-1; i++)
	  MPI_Irecv(recv_ptr[i],r_size[i] , MPI_PACKED, i+1, 13,
		    MPI_COMM_WORLD, &store_request[i]) ;

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, store_request, store_status) ;
	FATAL(err != MPI_SUCCESS) ;
	nsp->unpack(my_stuff, my_unpack, my_sz, te) ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  loc_unpack = 0 ;
	  nsp->unpack(recv_ptr[i], loc_unpack, r_size[i],
		      vseq[i]) ;
	}

	delete [] recv_size ;
        delete [] recv_size_bytes ;
        delete [] recv_buffer[0] ;
	delete [] recv_buffer ;
	delete [] status ;
	delete [] recv_request ;
	delete [] size_request ;
	delete [] size_status ;
	delete [] store_request ;
	delete [] store_status ;
        delete [] r_size ;
        delete [] my_stuff ;
        delete [] recv_ptr[0] ;
	delete [] recv_ptr ;
      }
      else {

        entitySet dom = sp->domain() ;
	entitySet temp = dom & d->my_entities ;

	int send_size = temp.size() ;
	int *send_buffer = new int[send_size] ;
	int sz = sp->pack_size(temp) ;
	unsigned char *send_ptr = new unsigned char[sz] ;

	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  send_buffer[j++] = l2g[*ti] ;

	int loc_pack = 0;
	sp->pack(send_ptr, loc_pack, sz, temp) ;

	MPI_Send(&send_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 12, MPI_COMM_WORLD) ;
        MPI_Send(&sz,1,MPI_INT,0,14,MPI_COMM_WORLD) ;
	MPI_Send(send_ptr, sz, MPI_PACKED, 0, 13, MPI_COMM_WORLD) ;
	delete [] send_buffer ;
	delete [] send_ptr ;
      }
    }
    return nsp ;

  }

  //This routine is used to collect the entire store on processor
  //0. The store is initially distributed across all the processors in
  //their initial global numbering. This is not a scalable option for
  //large containers. Written for use in the couple_grids routine in
  //fvm_coupler.cc
  storeRepP collect_global_store(storeRepP &sp) {
    storeRepP nsp = sp ;
    if(Loci::MPI_rank == 0) {
      std::vector<sequence> vseq(MPI_processes-1) ;
      MPI_Status *status, *size_status, *store_status ;
      MPI_Request *recv_request, *size_request, *store_request ;
      int **recv_buffer ;
      int *recv_size ;
      entitySet re ;
      entitySet temp = sp->domain() ;
      sequence te = sequence(temp) ;
      re += temp ;
      recv_size = new int[MPI_processes-1] ;
      int *recv_size_bytes = new int[MPI_processes-1] ;
      size_request = new MPI_Request[MPI_processes-1] ;
      size_status = new MPI_Status[MPI_processes-1] ;
      for(int k = 0; k < MPI_processes-1; k++)
	MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,11, MPI_COMM_WORLD, &size_request[k]);

#ifdef DEBUG
      int err =
#endif
	MPI_Waitall(MPI_processes-1, size_request, size_status) ;
      FATAL(err != MPI_SUCCESS) ;

      recv_buffer = new int*[MPI_processes-1] ;
      int recv_size_total = 0 ;
      for(int k=0;k<MPI_processes-1;++k)
	recv_size_total += recv_size[k] ;

      recv_buffer[0] = new int[recv_size_total] ;
      for(int i = 1; i < MPI_processes-1; ++i)
	recv_buffer[i] = recv_buffer[i-1] + recv_size[i-1] ;
      recv_request = new MPI_Request[MPI_processes-1] ;
      status = new MPI_Status[MPI_processes-1] ;

      for(int k = 0; k < MPI_processes-1; k++)
	MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,12, MPI_COMM_WORLD, &recv_request[k] );

#ifdef DEBUG
      err =
#endif
	MPI_Waitall(MPI_processes-1, recv_request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int k = 0; k < MPI_processes-1; k++)
	MPI_Irecv(&recv_size_bytes[k],1,MPI_INT, k+1,14,
	 	  MPI_COMM_WORLD, &size_request[k]) ;

#ifdef DEBUG
      err =
#endif
 	MPI_Waitall(MPI_processes-1, size_request, size_status) ;
      FATAL(err != MPI_SUCCESS) ;

      for(int k = 0; k < MPI_processes-1; ++k) {
	sequence tempseq ;
	for(int i = 0 ; i < recv_size[k]; ++i) {
	  re += recv_buffer[k][i] ;
	  tempseq += recv_buffer[k][i] ;
	}
	vseq[k] = tempseq ;
      }
      int my_sz= sp->pack_size(temp) ;
      int my_unpack = 0 ;
      int loc_unpack = 0 ;
      int loc_pack = 0 ;
      int *r_size = new int[MPI_processes-1] ;
      store_request = new MPI_Request[MPI_processes-1] ;
      store_status = new MPI_Status[MPI_processes-1] ;

      int sz = 0 ;
      for(int i = 0; i < MPI_processes-1; ++i) {
	r_size[i] = recv_size_bytes[i] ;
	sz += r_size[i] ;
      }

      unsigned char **recv_ptr = new unsigned char*[MPI_processes-1] ;
      unsigned char* my_stuff = new unsigned char[my_sz] ;
      sp->pack(my_stuff, loc_pack,my_sz,temp) ;
      nsp = sp->new_store(re) ;
      recv_ptr[0] = new unsigned char[sz] ;
      for(int i = 1; i < MPI_processes-1; i++)
	recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;

      for(int i = 0; i < MPI_processes-1; i++)
	MPI_Irecv(recv_ptr[i],r_size[i] , MPI_PACKED, i+1, 13,
		  MPI_COMM_WORLD, &store_request[i]) ;

#ifdef DEBUG
      err =
#endif
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
      FATAL(err != MPI_SUCCESS) ;
      nsp->unpack(my_stuff, my_unpack, my_sz, te) ;
      for(int i = 0; i < MPI_processes-1; ++i) {
	loc_unpack = 0 ;
	nsp->unpack(recv_ptr[i], loc_unpack, r_size[i],
		    vseq[i]) ;
      }

      delete [] recv_size ;
      delete [] recv_size_bytes ;
      delete [] recv_buffer[0] ;
      delete [] recv_buffer ;
      delete [] status ;
      delete [] recv_request ;
      delete [] size_request ;
      delete [] size_status ;
      delete [] store_request ;
      delete [] store_status ;
      delete [] r_size ;
      delete [] my_stuff ;
      delete [] recv_ptr[0] ;
      delete [] recv_ptr ;
    }
    else {
      entitySet dom = sp->domain() ;
      entitySet temp = dom ;
      int send_size = temp.size() ;
      int *send_buffer = new int[send_size] ;
      int sz = sp->pack_size(temp) ;
      unsigned char *send_ptr = new unsigned char[sz] ;

      int j = 0 ;
      for(entitySet::const_iterator ti = temp.begin(); ti != temp.end(); ++ti)
	send_buffer[j++] = *ti ;

      int loc_pack = 0;
      sp->pack(send_ptr, loc_pack, sz, temp) ;

      MPI_Send(&send_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
      MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 12, MPI_COMM_WORLD) ;
      MPI_Send(&sz,1,MPI_INT,0,14,MPI_COMM_WORLD) ;
      MPI_Send(send_ptr, sz, MPI_PACKED, 0, 13, MPI_COMM_WORLD) ;
      delete [] send_buffer ;
      delete [] send_ptr ;
    }
    return nsp ;
  }

}

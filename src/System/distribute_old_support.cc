#include "distribute_old_support.h"
#include <entitySet.h>

#include "Tools/debugger.h"

#include <vector>
using std::vector ;

#include <string>
using std::string ;

#include <mpi.h>

#include <iostream>
#include <fstream>
using std::cout ;
using std::cerr ; 
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::ofstream ;
using std::ostringstream ;

#include <algorithm>
using std::sort ;

namespace Loci {

  //This routine was first written to read in a partition of
  //entities. No longer needed if we are relying completely on the
  //scalable version.  
  vector<entitySet> read_partition(const char *fname,int num_partitions) {
    vector<entitySet> ptn ;
    ifstream infile ;
    ostringstream oss ;
    int part ;
    oss << fname << "." << num_partitions ;
    string filename = oss.str() ; 
    infile.open(filename.c_str(), ios::in) ;
    if(infile.fail()) {
      cerr << "File " << filename <<  "   not found \n First create the file using -exit option \n " << endl ;
      Finalize() ;
      exit(0) ;
    }
    infile >> part ;
    FATAL(part != num_partitions) ;
    int *partition = new int[num_partitions] ;
    for(int i = 0; i < num_partitions; ++i) 
      infile >> partition[i] ;
    for(int i = 0; i < num_partitions; ++i) {
      entitySet parti ;
      for(int j = 0; j < partition[i]; ++j) {
        infile >> part ;
        parti += part ;
      }
      ptn.push_back(parti);
    }
    return ptn ;
  }    

  //This routine writes out a generalized partition of entities. It
  //calls the generalized non scalable partitioning routine and writes 
  //out p partitions . 
  void write_partition(const char *fname, const vector<entitySet> &ptn) {
    if(MPI_rank == 0) {
      int num_partitions = ptn.size() ;
      ostringstream oss ;
      oss << fname << "." << num_partitions ;
      string filename = oss.str() ;
      ofstream ofile ;
      entitySet::const_iterator tei ;
      ofile.open(filename.c_str(), ios::out) ;
      ofile << num_partitions << endl ;
      for(size_t i = 0; i < ptn.size(); ++i) 
        ofile << ptn[i].size() << endl ;
      for(int i = 0; i < num_partitions; ++i) {
        entitySet temp = ptn[i];
        for(tei = temp.begin(); tei != temp.end(); ++tei)
          ofile << *tei << endl ;
      }
    }
  }
  
  //This  routine is almost similar to the above one but for the
  //introduction of the chopped partitioning and the remapping of
  //entities needed for scalable I/O
  vector<entitySet> generate_scalable_distribution(fact_db &facts, rule_db &rdb, int num_partitions) {
    if(num_partitions == 0)
      num_partitions = MPI_processes ;
    vector<entitySet> ptn ;
    if(num_partitions == 1)
      return ptn ;
    
    debugout << "Synchronising before metis_facts" << endl ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    double start = MPI_Wtime() ;
    
    metis_facts(facts,ptn,num_partitions) ;

    double end_time  = MPI_Wtime() ;
    debugout << "Time taken for metis_facts = " << end_time -start << endl ;
    std::vector<entitySet> partition(Loci::MPI_processes) ;
    dMap remap ;
    for(int i = 0 ; i < Loci::MPI_processes; ++i) {
      entitySet tmp_set = ptn[i] ;
      if(Loci::MPI_rank == i) {
	FORALL(tmp_set, ei) {
	  remap[ei] = ei ;
	} ENDFORALL ;
      }
    }

    fact_db::distribute_infoP df = new fact_db::distribute_info  ;
    //Remapping of entities keeps the partitioning contiguous . 
    df->remap = remap ;
    facts.put_init_ptn(ptn) ;
    facts.put_distribute_info(df) ;
    return ptn ;
  }
  
  void print_global(entitySet e, fact_db &facts) {
    if(facts.isDistributed()) {  
      MPI_Status *status ; 
      MPI_Request *recv_request ;
      int MAX = 100 ;
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      d = facts.get_distribute_info() ;
      l2g = facts.get_variable("l2g") ;
      if(d->myid == 0) {
	entitySet re ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	for(ti = e.begin(); ti != e.end(); ++ti)
	  re += l2g[*ti] ;
	recv_size = new int[MPI_processes-1] ;
	recv_buffer = new int*[MPI_processes-1] ;
	for(int i = 0; i < MPI_processes-1; ++i)
	  recv_buffer[i] = new int[MAX] ;
	recv_request = (MPI_Request *) malloc((MPI_processes-1) * sizeof(MPI_Request) ) ;
	status = (MPI_Status *) malloc((MPI_processes-1) * sizeof(MPI_Status) ) ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0],MAX,MPI_INT, k+1,1, MPI_COMM_WORLD, &recv_request[k] );  

#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(k = 0; k < MPI_processes-1; ++k)
	  MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
	
	for(k = 0; k < MPI_processes-1; ++k) {      
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    re += recv_buffer[k][i] ;
	  }
	}
	cout << "   " << re << endl ; 
	delete [] recv_size ;
	delete [] recv_buffer ;
	
      }
      else {
	int *send_buffer;
	int send_size ;
	
	entitySet temp;
	send_size = e.size() ;
	send_buffer = new int[send_size] ;
	
	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[j] = *ti ;
	  ++j ;
	}
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 1, MPI_COMM_WORLD) ;
	
	delete [] send_buffer ;
      } 
    }
  }
    
  Map distribute_whole_map(Map &m) {
    if(Loci::MPI_processes == 1)
      return m ;
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    dom = Loci::all_collect_entitySet(dom) ;

    nm.allocate(dom) ;
    int sz = 0 ;
    sz = dom.size() ;
    int *recv_ptr = new int[sz] ;
    if(Loci::MPI_rank == 0) {
      int tmp = 0 ;
      for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) 
	recv_ptr[tmp++] = m[*ei] ;
    }
    int tmp = 0 ;
    MPI_Bcast(recv_ptr, sz, MPI_INT, 0, MPI_COMM_WORLD) ;
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) 
      nm[*ei] = recv_ptr[tmp++] ;
    delete [] recv_ptr ;
    return nm ;
  }
  
  Map distribute_global_map(Map &m, const std::vector<entitySet> &init_ptn) {
    if(Loci::MPI_processes == 1)
      return m ;
    
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    entitySet image ;
    if(Loci::MPI_rank == 0)
      image = MapRepP(m.Rep())->image(dom) ;
    image = all_collect_entitySet(image) ;

    entitySet my ;
    if(Loci::MPI_rank == 0) {
      std::vector<entitySet> ent(MPI_processes) ;
      MPI_Status *store_status ;
      MPI_Request *store_request ;
      int sz = 0 ;
      int *s_size = new int[MPI_processes] ;
      
      for(int i = 0; i < MPI_processes; ++i) 
	ent[i] = init_ptn[i] & image ;
      for(int i = 0; i < MPI_processes-1; ++i) { 
	s_size[i] = 2 * ent[i+1].size() ;
	sz += s_size[i] ;
      }
      std::vector<entitySet> final_ent(MPI_processes) ;
      for(int i = 0; i < MPI_processes; ++i) {
	dom = MapRepP(m.Rep())->preimage(ent[i]).first ;
	final_ent[i] = dom ;
      }
      /*
	for(int i = 0; i < MPI_processes; ++i) {   
	for(int j = i+1 ; j < MPI_processes; ++j) {
	entitySet  tmp = final_ent[i] & final_ent[j] ;
	if(tmp != EMPTY)
	cerr << " ERROR: Something is screwed up ???  " << endl ;
	final_ent[j] -= tmp ;
	}
	}
      */
      int **send_ptr = new int*[MPI_processes] ;
      store_request = new MPI_Request[MPI_processes] ;
      store_status = new MPI_Status[MPI_processes] ;
      int tmp = 0 ;
      send_ptr[0] = new int[sz*10] ;
      for(int i = 1; i < MPI_processes-1; i++)
	send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
      for(int i = 0; i < MPI_processes-1; i++) {
	tmp = 0 ;
	for(entitySet::const_iterator ei = final_ent[i+1].begin(); ei != final_ent[i+1].end(); ++ei) {
	  send_ptr[i][tmp++] = *ei ;
	  send_ptr[i][tmp++] = m[*ei] ;
	  
	}
	MPI_Isend(send_ptr[i], s_size[i], MPI_INT, i+1, 3,
		  MPI_COMM_WORLD, &store_request[i]) ;
      }
#ifdef DEBUG
      int err =
#endif
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
      FATAL(err != MPI_SUCCESS) ;
      if(Loci::MPI_processes > 1)
	my = final_ent[0] ;
      else
	my = dom ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei)
	nm[*ei] = m[*ei] ;
      delete [] store_request ;
      delete [] store_status ;
      delete [] send_ptr ;
    }
    else {
      MPI_Status stat ;
      int sz = 2*(init_ptn[Loci::MPI_rank] & image).size();
      int *recv_ptr = new int[sz*10] ;
      MPI_Recv(recv_ptr, sz, MPI_INT, 0, 3, MPI_COMM_WORLD, &stat) ;
      int tmp = 1 ;
      dom = EMPTY ;
      for(int i = 0; i < sz; i += 2)
	dom += recv_ptr[i] ;
      nm.allocate(dom) ;
      for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) {
	nm[*ei] = recv_ptr[tmp] ;
	tmp += 2 ;
      }
      delete [] recv_ptr ;
    }
    return nm ;
  }
  
  Map distribute_gmap(Map &m, const std::vector<entitySet> &init_ptn) {
    if(Loci::MPI_processes == 1)
      return m ;
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    dom = Loci::all_collect_entitySet(dom) ;

    entitySet my ;
    if(Loci::MPI_rank == 0) {
      std::vector<entitySet> ent(MPI_processes-1) ;
      MPI_Status *store_status ;
      MPI_Request *store_request ;
      int sz = 0 ;
      int *s_size = new int[MPI_processes-1] ;
      for(int i = 0; i < MPI_processes-1; ++i) { 
	ent[i] = init_ptn[i+1] & dom ;
	s_size[i] = ent[i].size() ;
	sz += s_size[i] ;
      }
      int **send_ptr = new int*[MPI_processes-1] ;
      store_request = new MPI_Request[MPI_processes-1] ;
      store_status = new MPI_Status[MPI_processes-1] ;
      int tmp = 0 ;
      send_ptr[0] = new int[sz] ;
      for(int i = 1; i < MPI_processes-1; i++)
	send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
      for(int i = 0; i < MPI_processes-1; i++) {
	tmp = 0 ;
	for(entitySet::const_iterator ei = ent[i].begin(); ei != ent[i].end(); ++ei) 
	  send_ptr[i][tmp++] = m[*ei] ;
	MPI_Isend(send_ptr[i], s_size[i], MPI_INT, i+1, 3,
		  MPI_COMM_WORLD, &store_request[i]) ;
      }
#ifdef DEBUG
      int err =
#endif
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
      FATAL(err != MPI_SUCCESS) ;
      my = init_ptn[0] & dom ; ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei)
	nm[*ei] = m[*ei] ;
      delete [] store_request ;
      delete [] store_status ;
      delete [] send_ptr ;
    }
    else {
      MPI_Status stat ;
      my = init_ptn[Loci::MPI_rank] & dom ;
      int sz = my.size();
      int *recv_ptr = new int[sz] ;
      MPI_Recv(recv_ptr, sz, MPI_INT, 0, 3, MPI_COMM_WORLD, &stat) ;
      int tmp = 0 ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei) 
	nm[*ei] = recv_ptr[tmp++] ;
      delete [] recv_ptr ;
    }
    return nm ;
  }
  
  std::pair< storeRepP, Map > send_clone( storeRepP& sp, Map &m,  entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 0 ;
      FORALL(send_dom[i], si) {
	entitySet tmp = interval(si, si) ;
	send_count[i] +=  sp->pack_size(tmp) ;
      } ENDFORALL ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    Map m_int ;
    int tot = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      tot += entitySet(recv_dom[i]).size() ; //+= entitySet(recv_dom[i]) ;
    
    tmp_dom = interval(0, tot-1) ;
    
    storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] = 0 ; 
      FORALL(entitySet(recv_dom[i]), ri) {
	entitySet tmp = interval(ri, ri) ; 
	recv_count[i] +=  tmp_sp->pack_size(tmp) ;
      } ENDFORALL ;
      size_recv += recv_count[i] ;
    }
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    dmultiMap inverse ;
    entitySet m_dom = m.domain() ;
    entitySet range =  MapRepP(m.Rep())->image(m_dom) ;
    inverseMap(inverse, m, range, m_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      FORALL(send_dom[i], si) {
	entitySet tmp = interval(inverse[si][0], inverse[si][0]) ;
	sp->pack(send_store, loc_pack, size_send, tmp) ;
      } ENDFORALL ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    m_int.allocate(tmp_dom) ;
    tot = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(sequence::const_iterator si = recv_dom[i].begin() ; si != recv_dom[i].end(); ++si) {
	m_int[tot] = *si ;
	++tot ;
      }
    }
    FORALL(tmp_dom, ti) {
      entitySet tmp = interval(ti, ti) ;
      sequence tmp_seq = sequence(tmp) ;
      tmp_sp->unpack(recv_store, loc_pack, size_recv, tmp_seq) ; 
    } ENDFORALL ;
    
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return std::make_pair(tmp_sp, m_int)  ;
  }
  
  std::vector< storeRepP> fill_global_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec( MPI_processes) ;
    std::vector< storeRepP> tmp_sp( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      recv_count[i] =  sp->pack_size(entitySet(recv_dom[i])) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp[i]->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }
  
  entitySet broadcastEntitySet(entitySet s, int root) {
    int num_ivals = s.num_intervals() ;

    MPI_Bcast(&num_ivals,1,MPI_INT,root,MPI_COMM_WORLD) ;
    int *buf = new int[num_ivals*2] ;
    if(MPI_rank == root) {
      for(int i=0;i<num_ivals;++i) {
        buf[2*i] = s[i].first ;
        buf[2*i+1] = s[i].second ;
      }
    }
    MPI_Bcast(&buf,num_ivals*2,MPI_INT,root,MPI_COMM_WORLD) ;
    entitySet ret_val ;
    for(int i=0;i<num_ivals;++i)
      ret_val += interval(buf[2*i],buf[2*i+1]) ;
    delete[] buf ;
    return ret_val ;
  }

}

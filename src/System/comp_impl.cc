#include "comp_tools.h"
#include <distribute.h>
namespace Loci {
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts)  {
    do_run = true ;
    if(seq.num_intervals() == 0)
      do_run = false ;
    
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p)
  {
    do_run = true ;
    if(seq.num_intervals() == 0)
      do_run = false ;
    
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  void execute_rule::execute(fact_db &facts) {
    if(do_run)
      rp->compute(exec_seq) ;
  }
  
  void execute_rule::Print(ostream &s) const {
    if(MPI_processes > 1)
      s <<  "processor : " << MPI_rank  << "\t" ;
    s << rule_tag << " over sequence " << exec_seq << endl ;
  }
  
  execute_postcomm::execute_postcomm(std::list<comm_info> clist, fact_db &facts) {
    comm_list = clist ;
  }
   
  execute_precomm::execute_precomm(std::list<comm_info> plist, fact_db &facts) {
    comm_list = plist ;
  }
  
  void execute_precomm::execute(fact_db &facts) {
    
    MPI_Status *status ;
    MPI_Request *request ;
    std::list<comm_info>::const_iterator li ;
    std::vector<proc_details>::const_iterator vi ;
    int *p_size, *r_size ,*loc_pack, *loc_unpack ;
    unsigned char **send_ptr, **recv_ptr ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    int n, recv_count = 0, send_count = 0, recv_flag = 0 ;
    entitySet procs = d->send_neighbour | d->recv_neighbour ;
    entitySet::const_iterator ei ;
    n = procs.size() ;
    p_size = new int[n] ;
    r_size = new int[n] ;
    send_ptr = new unsigned char*[n] ;
    recv_ptr = new unsigned char*[n] ;
    int *send_procs = new int[n] ;
    int *recv_procs = new int[n] ;
    loc_pack = new int[n] ;
    loc_unpack = new int[n] ;
    Map id ;
    id.allocate(procs) ;
    int k = 0 ;
    for(ei = procs.begin(); ei != procs.end(); ++ei) {
      id[*ei] = k ;
      k++; 
    }
    for(int i = 0; i < n; ++i) {
      p_size[i] = 0 ;
      r_size[i] = 0 ;
      loc_pack[i] = 0 ;
      loc_unpack[i] = 0 ;
    }
    storeRepP sp ;
    entitySet var ;
    std::vector<std::vector<variable> > vvr(n) ;
    std::vector<std::vector<variable> > vvs(n) ;
    std::vector<variable>::const_iterator vvi ;
    std::vector<sequence> vseq(n) ; 
    std::vector<entitySet> vset(n) ;
    entitySet recv, send ;
     
    for(li = comm_list.begin(); li != comm_list.end(); ++li) {
      sp = facts.get_variable(li->v) ;
      //cout << "variable = " << li->v <<  endl ;
      for(vi = li->recv_info.begin() ; vi != li->recv_info.end(); ++vi) {
	if(vi->recv_set != EMPTY) {
	  r_size[id[vi->processor]] += sp->pack_size(entitySet(vi->recv_set)) ;
	  vseq[id[vi->processor]] += vi->recv_set ;
	  recv += vi->processor ;
	  vvr[id[vi->processor]].push_back(li->v) ; 
	  // cout << d->myid << "variable  = " << li->v << "    receiving  " << vi->recv_set << "  from   " << vi->processor << endl ;
	}
      }
      for(vi = li->send_info.begin(); vi != li->send_info.end(); ++vi) {
	if(vi->send_set != EMPTY) {
	  p_size[id[vi->processor]] += sp->pack_size(vi->send_set) ;
	  vset[id[vi->processor]] += vi->send_set ; 
	  send += vi->processor ;
	  vvs[id[vi->processor]].push_back(li->v) ;
	  //cout << d->myid <<" variable =  " << li->v <<  "    sending  " << vi->send_set << "  to   " << vi->processor << endl ;
	}
      }
    }
    
    recv_count = 0 ;
    for(ei = recv.begin(); ei != recv.end(); ++ei) {
      recv_ptr[id[*ei]] = new unsigned char[r_size[id[*ei]]] ;
      recv_procs[recv_count] = *ei ;
      recv_count++ ;
    } 
    send_count = 0 ;
    for(ei = send.begin(); ei != send.end(); ++ei) {
      send_ptr[id[*ei]]=new unsigned char[p_size[id[*ei]]] ;
      for(vvi = vvs[id[*ei]].begin(); vvi != vvs[id[*ei]].end(); ++vvi) {
	//debugout[MPI_rank] << "precomm  processor   " <<  d->myid << "variable  =  " << *vvi << " { " ;
        sp = facts.get_variable(*vvi) ;
	sp->pack(send_ptr[id[*ei]], loc_pack[id[*ei]],
		 p_size[id[*ei]], vset[id[*ei]]) ;
	
	//debugout[MPI_rank] << " } " << endl ;
	
      } 
      send_procs[send_count] = *ei ;
      send_count++ ;
    }
    request =  new MPI_Request[recv_count] ;
    status =  new MPI_Status[recv_count] ;
    
    for(int i = 0; i < recv_count; i++) {
      MPI_Irecv(recv_ptr[id[recv_procs[i]]], r_size[id[recv_procs[i]]], MPI_PACKED,
      	recv_procs[i], 1, MPI_COMM_WORLD, &request[i]) ;  
      recv_flag = 1 ;
    }
    
    for(int i = 0; i < send_count; ++i) 
      MPI_Send(send_ptr[id[send_procs[i]]], p_size[id[send_procs[i]]], MPI_PACKED, send_procs[i], 1, MPI_COMM_WORLD) ;
      
    if((recv_flag))
      MPI_Waitall(recv_count, request, status) ;
    
    if(recv_flag) {
      for(ei = recv.begin(); ei != recv.end(); ++ei) {
	for(vvi = vvr[id[*ei]].begin(); vvi != vvr[id[*ei]].end(); ++vvi) {
	  sp = facts.get_variable(*vvi) ;
	  //debugout[MPI_rank] << "precomm   processor   " <<  d->myid << "   variable  =  " << *vvi <<  " { " ; 
	  sp->unpack(recv_ptr[id[*ei]], loc_unpack[id[*ei]], 
		     r_size[id[*ei]], vseq[id[*ei]]) ;
	  //debugout[MPI_rank] << " } " << endl ;
	  
	}
      }
    }
    delete [] status ;
    delete [] request ;
    delete [] p_size ;
    delete [] r_size ;
    delete [] send_ptr ;
    delete [] recv_ptr ;
    delete [] loc_pack ;
    delete [] loc_unpack ;
    delete [] send_procs ;
    delete [] recv_procs ;
    
  }
  void execute_postcomm::execute(fact_db &facts) {
    
    MPI_Status *status ;
    MPI_Request *request ;
    std::list<comm_info>::const_iterator li ;
    std::vector<proc_details>::const_iterator vi ;
    int *p_size, *r_size ,*loc_pack, *loc_unpack ;
    unsigned char **send_ptr, **recv_ptr ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    int n, recv_count = 0, send_count = 0, recv_flag = 0 ;
    entitySet procs = d->send_neighbour | d->recv_neighbour ;
    entitySet::const_iterator ei ;
    n = procs.size() ;
    p_size = new int[n] ;
    r_size = new int[n] ;
    send_ptr = new unsigned char*[n] ;
    recv_ptr = new unsigned char*[n] ;
    int *send_procs = new int[n] ;
    int *recv_procs = new int[n] ;
    loc_pack = new int[n] ;
    loc_unpack = new int[n] ;
    Map id ;
    id.allocate(procs) ;
    int k = 0 ;
    for(ei = procs.begin(); ei != procs.end(); ++ei) {
      id[*ei] = k ;
      k++; 
    }
    for(int i = 0; i < n; ++i) {
      p_size[i] = 0 ;
      r_size[i] = 0 ;
      loc_pack[i] = 0 ;
      loc_unpack[i] = 0 ;
    }
    storeRepP sp ;
    entitySet var ;
    std::vector<std::vector<variable> > vvr(n) ;
    std::vector<std::vector<variable> > vvs(n) ;
    std::vector<variable>::const_iterator vvi ;
    std::vector<sequence> vseq(n) ; 
    std::vector<entitySet> vset(n) ;
    entitySet recv, send ;
    
    for(li = comm_list.begin(); li != comm_list.end(); ++li) {
      sp = facts.get_variable(li->v) ;
      //cout << "variable = " << li->v <<  endl ;
      for(vi = li->recv_info.begin() ; vi != li->recv_info.end(); ++vi) {
	if(vi->recv_set != EMPTY) {
	  r_size[id[vi->processor]] += sp->pack_size(entitySet(vi->recv_set)) ;
	  vseq[id[vi->processor]] += vi->recv_set ;
	  recv += vi->processor ;
	  vvr[id[vi->processor]].push_back(li->v) ; 
	  // cout << d->myid << "variable  = " << li->v << "    receiving  " << vi->recv_set << "  from   " << vi->processor << endl ;
	}
      }
      for(vi = li->send_info.begin(); vi != li->send_info.end(); ++vi) {
	if(vi->send_set != EMPTY) {
	  p_size[id[vi->processor]] += sp->pack_size(vi->send_set) ;
	  vset[id[vi->processor]] += vi->send_set ; 
	  send += vi->processor ;
	  vvs[id[vi->processor]].push_back(li->v) ;
	  //cout << d->myid <<" variable =  " << li->v <<  "    sending  " << vi->send_set << "  to   " << vi->processor << endl ;
	}
      }
    }
    
    recv_count = 0 ;
    for(ei = recv.begin(); ei != recv.end(); ++ei) {
      recv_ptr[id[*ei]] = new unsigned char[r_size[id[*ei]]] ;
      recv_procs[recv_count] = *ei ;
      recv_count++ ;
    } 
    send_count = 0 ;
    for(ei = send.begin(); ei != send.end(); ++ei) {
      send_ptr[id[*ei]]=new unsigned char[p_size[id[*ei]]] ;
      for(vvi = vvs[id[*ei]].begin(); vvi != vvs[id[*ei]].end(); ++vvi) {
	//debugout[MPI_rank] << "postcomm  processor   " <<  d->myid << "variable  =  " << *vvi << " { " ;
        sp = facts.get_variable(*vvi) ;
	sp->pack(send_ptr[id[*ei]], loc_pack[id[*ei]],
		 p_size[id[*ei]], vset[id[*ei]]) ;
	//debugout[MPI_rank] << " } " << endl ;
      }
      send_procs[send_count] = *ei ;
      send_count++ ;
    }
    
    request =  new MPI_Request[recv_count] ;
    status =  new MPI_Status[recv_count] ;
    
    for(int i = 0; i < recv_count; i++) {
      MPI_Irecv(recv_ptr[id[recv_procs[i]]], r_size[id[recv_procs[i]]], MPI_PACKED,
		recv_procs[i], 1, MPI_COMM_WORLD, &request[i]) ;  
      recv_flag = 1 ;
    }
    
    for(int i = 0; i < send_count; ++i) 
      MPI_Send(send_ptr[id[send_procs[i]]], p_size[id[send_procs[i]]], MPI_PACKED, send_procs[i], 1, MPI_COMM_WORLD) ;
      
    if((recv_flag))
      MPI_Waitall(recv_count, request, status) ;
    
    if(recv_flag) {
      for(ei = recv.begin(); ei != recv.end(); ++ei) {
	for(vvi = vvr[id[*ei]].begin(); vvi != vvr[id[*ei]].end(); ++vvi) {
	  sp = facts.get_variable(*vvi) ;
	  //debugout[MPI_rank] << "postcomm   processor   " <<  d->myid << "   variable  =  " << *vvi <<  " { " ; 
	  sp->unpack(recv_ptr[id[*ei]], loc_unpack[id[*ei]], 
		     r_size[id[*ei]], vseq[id[*ei]]) ;
	  //debugout[MPI_rank] << " } " << endl ;
	}
	 
      }
    }
    delete [] status ;
    delete [] request ;
    delete [] p_size ;
    delete [] r_size ;
    delete [] send_ptr ;
    delete [] recv_ptr ;
    delete [] loc_pack ;
    delete [] loc_unpack ;
    delete [] send_procs ;
    delete [] recv_procs ;
    
   }
  
  
  void execute_postcomm::Print(ostream &s) const {
    //cerr << MPI_rank <<"   sending the values to their respective processors " << endl ;
  }
  void execute_precomm::Print(ostream &s) const {
    
  }
     
  void impl_compiler::set_var_existence(fact_db &facts) {
    existential_rule_analysis(impl,facts) ;
  }
  
  void impl_compiler::process_var_requests(fact_db &facts) {
    exec_seq = process_rule_requests(impl,facts) ;
  }
  
  executeP impl_compiler::create_execution_schedule(fact_db &facts) {
    
#ifndef DEBUG
    if(exec_seq.size() == 0)
      return executeP(0) ;
#endif
    variableSet targets = impl.targets() ;
    WARN(targets.size() == 0) ;
    if(num_threads > 1 &&
       exec_seq.size() > num_threads*30 &&
       !impl.get_info().output_is_parameter &&
       (impl.get_info().rule_impl->get_rule_class()==rule_impl::POINTWISE||
        impl.get_info().rule_impl->get_rule_class()==rule_impl::UNIT)&&
       impl.get_info().rule_impl->thread_rule() &&
       (targets.begin()->get_info()).name != "OUTPUT") {
      execute_par *ep = new execute_par ;
      parallel_schedule(ep,exec_seq,impl,facts) ;
      return ep ;
    }
    if((targets.begin()->get_info()).name == "OUTPUT") {
      CPTR<execute_list> el = new execute_list ;
      el->append_list(new execute_rule(impl,sequence(exec_seq),facts)) ;
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
      return executeP(el) ;
    }
    return new execute_rule(impl,sequence(exec_seq),facts) ;
  }

}
 
 

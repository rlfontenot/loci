#include "comp_tools.h"
#include "distribute.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif

using std::hash_map ;
using std::list ; 

//#define VERBOSE
namespace Loci {
  class joiner_oper : public execute_modules {
    variable joiner_var ;
    storeRepP joiner_store ;
    vector<entitySet> partition ;
    vector<storeRepP> var_vec ;
    vector<CPTR<joiner> > join_ops ;
    CPTR<joiner> joiner_op ;
  public:
    joiner_oper(variable jv, storeRepP &js, vector<entitySet> &ptn,
                vector<storeRepP> &vv, CPTR<joiner> &jo) :
      joiner_var(jv), joiner_store(js), partition(ptn),var_vec(vv),
      joiner_op(jo) {
      for(int i=0;i<var_vec.size();++i) {
        join_ops.push_back(jo->clone()) ;
        join_ops[i]->SetArgs(joiner_store,var_vec[i]) ;
      }
    }
    virtual void execute(fact_db &facts) ;
    virtual void Print(ostream &s) const ;
  } ;

  void joiner_oper::execute(fact_db &facts) {
    for(int i=0;i<var_vec.size();++i) 
      join_ops[i]->Join(sequence(partition[i])) ;
    //joiner_op->Join(joiner_store,var_vec[i],sequence(partition[i])) ;
  }

  void joiner_oper::Print(ostream &s) const {
    s << "reducing thread results for variable " << joiner_var << endl ;
    s << "reducing partitions = " << endl ;
    for(int i=0;i<var_vec.size();++i)
      s << "p["<<i<< "]="<<partition[i]<<endl ;
  }
  
  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      if(*vi != reduce_var)
        sources &= facts.variable_existence(*vi) ;
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working &= facts.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }
  
  void apply_compiler::set_var_existence(fact_db &facts) {
    if(facts.isDistributed()) {

      // Compute the shadow entities produced by using this apply rules.
      // Any shadow entities that we don't own we will need to exchange
      // the partial results with other processors.
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      entitySet sources = d->my_entities ;
      entitySet constraints = d->my_entities ;
      
      warn(apply.targets().size() != 1) ;
      variable reduce_var = *apply.targets().begin() ;
      
      
      const rule_impl::info &rinfo = apply.get_info().desc ;
      
      bool outputmap = false ;
      set<vmap_info>::const_iterator si ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	if(si->mapping.size() != 0)
	  outputmap = true ;
      }
      // If there is no mapping in the output, then there will be no
      // shadow cast from this rule application.
      if(!outputmap)
        return ;
      
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        sources &= vmap_source_exist_apply(*si,facts,reduce_var) ;
      } 
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts) ;

      sources &= constraints ;
      
      entitySet context = sources & constraints ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        entitySet targets = vmap_target_exist(*si,facts,context) ;
        const variableSet &tvars = si->var ;
        variableSet::const_iterator vi ;
	for(vi=tvars.begin();vi!=tvars.end();++vi) {
#ifdef VERBOSE
	debugout[MPI_rank] << "shadow is " << targets << endl ;
	debugout[MPI_rank] << "shadow not owned is "
			   << targets - d->my_entities << endl
			   << "variable is " << *vi << endl ;
#endif
	facts.variable_shadow(*vi,targets) ;
      }
    }
  }
} 
  
  
  void apply_compiler::process_var_requests(fact_db &facts) {
    
    vdefmap tvarmap ;
    variableSet targets = apply.targets() ;
    variableSet sources = apply.sources() ;
    
    fatal(targets.size() != 1) ;
    variable tvar = *(targets.begin()) ;
    
    if(facts.get_variable(tvar)->RepType() == Loci::PARAMETER) 
      tvarmap[tvar] = facts.variable_existence(tvar) ;
    else
      tvarmap[tvar] = facts.get_variable_request(unit_tag,tvar) ;

    const rule_impl::info &rinfo = apply.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    entitySet compute ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      compute |= vmap_target_requests(*si,tvarmap,facts) ;
    }
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      compute &= d->my_entities ;
    }
    output_mapping = false ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end(); ++si) {
      variableSet::const_iterator vi ;
      entitySet comp = compute ;
      vector<variableSet>::const_iterator mi ;
      for(mi=si->mapping.begin();mi!=si->mapping.end();++mi) {
        output_mapping = true ;
        entitySet working ;
        for(vi=mi->begin();vi!=mi->end();++vi) {
          FATAL(!facts.is_a_Map(*vi)) ;
          working |= facts.image(*vi,comp) ;
        }
        comp = working ;
      }
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        if((comp - facts.variable_existence(*vi)) != EMPTY) {
          cerr << "ERROR: Apply rule " << apply <<  endl
               << " output mapping forces application to entities where unit does not exist." << endl ;
          cerr << "error occurs for entities " <<
            entitySet(comp-facts.variable_existence(*vi)) << endl ;
          cerr << "error occurs when applying to variable " << *vi << endl;
          cerr << "error is not recoverable, terminating scheduling process"
               << endl ;
          exit(-1) ;
        }
        facts.variable_request(*vi,comp) ;
      }
    }
    
    entitySet srcs = ~EMPTY ;
    entitySet cnstrnts = srcs ;
    entitySet my_entities = srcs ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      srcs = d->my_entities ;
      cnstrnts = d->my_entities ;
      my_entities = d->my_entities ;
    }
    
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si)
      srcs &= vmap_source_exist(*si,facts) ;
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      cnstrnts &= vmap_source_exist(*si,facts) ;
    if(rinfo.constraints.begin() != rinfo.constraints.end())
      if((srcs & cnstrnts) != cnstrnts) {
        cerr << "Warning, reduction rule:" << apply
             << "cannot supply all entities of constraint" << endl ;
        cerr << "constraints = " <<cnstrnts << endl ;
        entitySet sac = srcs & cnstrnts ;
        cerr << "srcs & constraints = " << sac << endl ;
        //      exit(-1) ;
        for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
          entitySet sources = vmap_source_exist(*si,facts) ;
          sources &= my_entities ;
          if((sources & cnstrnts) != cnstrnts) {
            cerr << "sources & constraints != constraints for input"
                 << endl
                 << sources  << " -- " << *si << endl ;
            
            if(si->mapping.size() > 0) {
              entitySet working = cnstrnts ;
              for(int i=0;i<si->mapping.size();++i) {
                entitySet images ;
                variableSet::const_iterator vi ;
                for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                  images |= facts.image(*vi,working) ;
                working = images ; 
              }
              variableSet::const_iterator vi ;
              for(vi=si->var.begin();vi!=si->var.end();++vi) {
                entitySet exist = facts.variable_existence(*vi) ;
                entitySet fails = working & ~exist ;
                if(fails != EMPTY) {
                  cerr << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                }
              }
            }
          }
        }
      }
    srcs &= cnstrnts ;
    
    // now trim compute to what can be computed.
    compute &= srcs ;
    exec_seq = compute ;
    /*
      if(facts.isDistributed()) {
      entitySet re = collect_entitySet(exec_seq, facts) ;
      if(Loci::MPI_rank == 0) { 
	cout << "processor  " << Loci::MPI_rank << "  exec_sequence =  " << re << endl ;
	}
	}
	else {
	
	cout  << "exec_sequence = " << exec_seq << endl ; 
	}
     */
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      entitySet requests = vmap_source_requests(*si,facts,compute) ;
      variableSet::const_iterator vi ;
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        variable v = *vi ;
        facts.variable_request(v,requests) ; 
#ifdef VERBOSE
	debugout[MPI_rank] << "rule " << apply << " requesting variable "
			   << v << " for entities " << requests << endl ;
#endif
      }
    }
        
#ifdef VERBOSE
    debugout[MPI_rank] << "rule " << apply << " computes over " << compute << endl ;
#endif
  }
  
  executeP apply_compiler::create_execution_schedule(fact_db &facts) {
#ifndef DEBUG
    if(exec_seq.size() == 0)
      return executeP(0) ;
#endif
    CPTR<execute_list> el = new execute_list ;
    if(num_threads == 1 || !apply.get_info().rule_impl->thread_rule() ||
       exec_seq.size() < num_threads*30 ) {
      el->append_list(new execute_rule(apply,sequence(exec_seq),facts)) ;
    } else if(!apply.get_info().output_is_parameter &&!output_mapping) {
      execute_par *ep = new execute_par ;
      parallel_schedule(ep,exec_seq,apply,facts) ;
      el->append_list(ep)  ;
    } else if(apply.get_info().output_is_parameter) {
      variableSet target = apply.targets() ;
      fatal(target.size() != 1) ;
      variable v = *(target.begin()) ;
      storeRepP sp = facts.get_variable(v) ;
      vector<entitySet> partition = partition_set(exec_seq,num_threads) ;
      vector<storeRepP> var_vec ;

      execute_par *ep = new execute_par ;
      for(int i=0;i<partition.size();++i) {
        storeRepP rp = sp->new_store(partition[i]) ;
        var_vec.push_back(rp) ;
        execute_sequence *es = new execute_sequence ;
        es->append_list(new execute_rule(unit_tag,sequence(partition[i]),
                                         facts,v,rp)) ;
        es->append_list(new execute_rule(apply,sequence(partition[i]),
                                         facts,v,rp)) ;
        ep->append_list(es) ;
      }
      el->append_list(ep) ;
      el->append_list(new execute_thread_sync) ;
      rule_implP arule = (apply.get_info().rule_impl);
      fatal(arule == 0) ;
      CPTR<joiner> j_op = (arule)->get_joiner() ;
      fatal(j_op == 0) ;
      el->append_list(new joiner_oper(v,sp,partition,var_vec,j_op)) ;
    } else {
      variableSet target = apply.targets() ;
      fatal(target.size() != 1) ;
      variable v = *(target.begin()) ;
      storeRepP sp = facts.get_variable(v) ;
      vector<entitySet> partition = partition_set(exec_seq,num_threads) ;

      const rule_impl::info &rinfo = apply.get_info().desc ;
      execute_par *ep = new execute_par ;
      entitySet apply_domain,all_contexts ;
      vector<entitySet> shards, shard_domains ;
      for(int i=0;i<partition.size();++i) {
        fatal(rinfo.targets.size() != 1) ;
        entitySet context = partition[i] ;
        entitySet pdom = vmap_target_exist(*rinfo.targets.begin(),facts,context) ;
      
        entitySet rem = pdom & apply_domain ;
        if(rem != EMPTY) {
          entitySet compute = rem ;
          const vmap_info &vmi = *rinfo.targets.begin() ;
          vector<variableSet>::const_reverse_iterator mi ;
          for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
            entitySet working = EMPTY ;
            variableSet::const_iterator vi;
            for(vi=mi->begin();vi!=mi->end();++vi) {
              FATAL(!facts.is_a_Map(*vi)) ;
              working |= facts.preimage(*vi,compute).second ;
            }
            compute = working ;
          }
          compute &= partition[i] ;
          shards.push_back(compute) ;
          entitySet sdom = vmap_target_exist(vmi,facts,compute) ;
          shard_domains.push_back(sdom) ;
          context &= ~compute ;
        }
        apply_domain |= pdom ;
        all_contexts |= partition[i] ;
        ep->append_list(new execute_rule(apply,sequence(context),facts)) ;
      }
      if(shards.size() == 0) {
        el->append_list(ep) ;
        el->append_list(new execute_thread_sync) ;
      } else {
        ep->append_list(new execute_sequence) ;
        bool disjoint = true ;
        entitySet dom_tot ;
        for(int i=0;i<shards.size();++i) {
          if((shard_domains[i] & dom_tot) != EMPTY)
            disjoint = false ;
          dom_tot += shard_domains[i] ;
        }
        variableSet target = apply.targets() ;
        fatal(target.size() != 1) ;
        variable v = *(target.begin()) ;
        storeRepP sp = facts.get_variable(v) ;
        fatal(sp == 0) ;
      
        vector<storeRepP> var_vec ;
      
        for(int i=0;i<shards.size();++i) {
          storeRepP rp = sp->new_store(shard_domains[i]) ;
          var_vec.push_back(rp) ;
          execute_sequence *es = new execute_sequence ;
          es->append_list(new execute_rule(unit_tag,sequence(shard_domains[i]),
                                           facts,v,rp)) ;
          es->append_list(new execute_rule(apply,sequence(shards[i]),
                                           facts,v,rp)) ;
          ep->append_list(es) ;
        }
      
        el->append_list(ep) ;
        el->append_list(new execute_thread_sync) ;
        //-----Now join the partial results

        rule_implP arule = (apply.get_info().rule_impl);
        fatal(arule == 0) ;

        if(disjoint) {
          execute_par *epj = new execute_par ;
          epj->append_list(new execute_sequence) ;
          for(int i=0;i<shard_domains.size();++i) {
            vector<entitySet> ve ;
            vector<storeRepP> vv ;
            ve.push_back(shard_domains[i]) ;
            vv.push_back(var_vec[i]) ;
            CPTR<joiner> j_op = (arule)->get_joiner() ;
            fatal(j_op == 0) ;
            epj->append_list(new joiner_oper(v,sp,ve,vv,j_op)) ;
          }
          el->append_list(epj) ;
          el->append_list(new execute_thread_sync) ;
        } else {
          for(int i=0;i<shard_domains.size();++i) {
            execute_par *epj = new execute_par ;
            vector<entitySet> decompose = partition_set(shard_domains[i],num_threads) ;
            vector<storeRepP> vv ;
            vv.push_back(var_vec[i]) ;
            for(int j=0;j<decompose.size();++j) {
              vector<entitySet> ve ;
              ve.push_back(decompose[j]) ;
              CPTR<joiner> j_op = (arule)->get_joiner() ;
              fatal(j_op == 0) ;
              epj->append_list(new joiner_oper(v,sp,ve,vv,j_op)) ;
            }
            el->append_list(epj) ;
            el->append_list(new execute_thread_sync) ;
          }
        }
      }
    }
    el->append_list(new execute_thread_sync) ;
    return executeP(el) ;
  }
  
  execute_param_red::execute_param_red(variable red, rule unit, CPTR<joiner> j_op) {
    reduce_var = red ;
    unit_rule = unit ;
    join_op = j_op ;
  }


  CPTR<joiner> global_join_op ;
  
  
  void create_user_function(void *send_ptr, void *result_ptr, int *size, MPI_Datatype* dptr) {
    storeRepP sp, tp ;
    int loc_send = 0, loc_result = 0 ;
    entitySet e ;
    sequence seq ;
    e += 1 ;
    seq += 1 ;
    sp = global_join_op->getTargetRep() ;
    tp = global_join_op->getTargetRep() ;
    tp->allocate(e) ;
    sp->allocate(e) ;
    sp->unpack(send_ptr, loc_send, *size, seq) ;
    tp->unpack(result_ptr, loc_result, *size, seq) ;
    global_join_op->SetArgs(tp, sp) ;
    global_join_op->Join(seq) ;
    loc_result = 0 ;
    loc_send = 0 ;
    tp->pack(result_ptr, loc_result, *size, e) ;
  } 
  
  void execute_param_red::execute(fact_db &facts) {
    void *send_ptr; 
    void *result_ptr ;
    int size ;
    int loc = 0 , loc_result = 0;
    storeRepP sp, tp ;
    entitySet e ;
    sequence seq ;
    MPI_Op create_join_op ;
    sp = facts.get_variable(reduce_var) ;
    size = sp->pack_size(e) ;
    send_ptr = new unsigned char[size] ;
    result_ptr = new unsigned char[size] ;
    sp->pack(send_ptr, loc, size, e) ;
    global_join_op = join_op ;
    MPI_Op_create(&create_user_function, 0, &create_join_op) ;
    MPI_Allreduce(send_ptr, result_ptr, size, MPI_PACKED, create_join_op, MPI_COMM_WORLD) ;
    sp->unpack(result_ptr, loc_result, size, seq) ;
  }
  
  void execute_param_red::Print(ostream &s) const {
    s << "param reduction on " << reduce_var << endl ;
  }
  void reduce_param_compiler::set_var_existence(fact_db &facts)  {
    
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      entitySet targets ;
      targets = facts.get_existential_info(reduce_var, unit_rule) ;
      targets += send_entitySet(targets, facts) ;
      targets &= d->my_entities ;
      targets += fill_entitySet(targets, facts) ;
      facts.set_existential_info(reduce_var,unit_rule,targets) ;
    }
  }
  
  void reduce_param_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      entitySet requests = facts.get_variable_requests(reduce_var) ;
      requests += send_entitySet(requests, facts) ;
      facts.variable_request(reduce_var,requests) ;
    }
  } 
  
  executeP reduce_param_compiler::create_execution_schedule(fact_db &facts) {
    ostringstream oss ;
    oss << "reduce param " << reduce_var ;
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ;
      el->append_list(new execute_thread_sync) ;
      el->append_list(new execute_param_red(reduce_var, unit_rule, join_op)) ; 
      return executeP(el) ;
    }
    return executeP(new execute_msg(oss.str())) ;
  }
  
  void reduce_store_compiler::set_var_existence(fact_db &facts)  {
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      entitySet targets = facts.get_existential_info(reduce_var, unit_rule) ;
      targets += send_entitySet(targets, facts) ;
      targets &= d->my_entities ;
      targets += fill_entitySet(targets, facts) ;
      facts.set_existential_info(reduce_var,unit_rule,targets) ;
    }
  }
  
  void swap_send_recv(list<comm_info> &cl) {
    list<comm_info>::iterator li ;
    for(li=cl.begin();li!=cl.end();++li) {
      entitySet tmp  = li->send_set ;
      li->send_set = entitySet(li->recv_set) ;
      li->recv_set = sequence(tmp) ;
    }
  }
  
  void reduce_store_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      variableSet vars ;
      vars += reduce_var ;
      list<comm_info> request_comm = barrier_process_rule_requests(vars,facts) ;
      
      entitySet requests = facts.get_variable_requests(reduce_var) ;
      entitySet shadow = facts.get_variable_shadow(reduce_var) ;
      shadow &= requests ;
      
      list<comm_info> slist ;
      entitySet response = send_requests(shadow, reduce_var,facts,slist) ;
      swap_send_recv(slist) ;
    
      rlist = sort_comm(slist,facts) ;
      clist = sort_comm(request_comm,facts) ;
      
#ifdef VERBOSE
    if(shadow != EMPTY) {
      debugout[MPI_rank] << "shadow = " << shadow << endl ;
      shadow -= d->my_entities ;
      debugout[MPI_rank] << "shadow/my_entites = " << shadow << endl ;
    }
#endif
    }
  }
  
class execute_comm_reduce : public execute_modules {
  vector<pair<int,vector<send_var_info> > > send_info ;
  vector<pair<int,vector<recv_var_info> > > recv_info ;
  CPTR<joiner> join_op ;
public:
  execute_comm_reduce(list<comm_info> &plist, fact_db &facts,
		      CPTR<joiner> jop) ;
  virtual void execute(fact_db &facts) ;
  virtual void Print(std::ostream &s) const ;
} ; 

execute_comm_reduce::execute_comm_reduce(list<comm_info> &plist,
					 fact_db &facts,
					 CPTR<joiner> jop) {
  join_op = jop ;
  hash_map<int,vector<send_var_info> > send_data ;
  hash_map<int,vector<recv_var_info> > recv_data ;
  list<comm_info>::const_iterator cli ;
  intervalSet send_procs, recv_procs ;
  for(cli=plist.begin();cli!=plist.end();++cli) {
    variable v = cli->v ;
    if(cli->send_set.size() > 0) {
      int send_proc = cli->processor ;
      send_procs += send_proc ;
      entitySet send_set = cli->send_set ;
      send_data[send_proc].push_back(send_var_info(v,send_set)) ;
    }
    if(cli->recv_set.size() > 0) {
      int recv_proc = cli->processor ;
      sequence recv_seq = cli->recv_set ;
      recv_procs += recv_proc ;
      recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
    }
  }
    
    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
    }
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
    }
    
  }

  void execute_comm_reduce::execute(fact_db  &facts) {
    const int nrecv = recv_info.size() ;
    int *r_size = new int[nrecv] ;
    int total_size = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = 0 ;
      for(int j=0;j<recv_info[i].second.size();++j) {
        storeRepP sp = facts.get_variable(recv_info[i].second[j].v) ;
        r_size[i] += sp->pack_size(entitySet(recv_info[i].second[j].seq)) ;
      }
      total_size += r_size[i] ;
    }
    unsigned char **recv_ptr = new unsigned char*[nrecv] ;
    recv_ptr[0] = new unsigned char[total_size] ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1]+r_size[i-1] ;
    
    const int nsend = send_info.size() ;
    int *s_size = new int[nsend] ;
    total_size = 0 ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(int j=0;j<send_info[i].second.size();++j) {
        storeRepP sp = facts.get_variable(send_info[i].second[j].v) ;
        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
      }
      total_size += s_size[i] ;
    }
    unsigned char **send_ptr = new unsigned char*[nsend] ;
    send_ptr[0] = new unsigned char[total_size] ;
    for(int i=1;i<nsend;++i)
      send_ptr[i] = send_ptr[i-1]+s_size[i-1] ;
    
    MPI_Request *request =  new MPI_Request[nrecv] ;
    MPI_Status *status =  new MPI_Status[nrecv] ;
    
    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
    }
    
    // Pack the buffer for sending 
    for(int i=0;i<nsend;++i) {
      /*
	#ifdef VERBOSE
	debugout[MPI_rank] << "sending to processor " << send_info[i].first
	<< endl ;
	#endif
      */
      int loc_pack = 0 ;
      for(int j=0;j<send_info[i].second.size();++j) {
        storeRepP sp = facts.get_variable(send_info[i].second[j].v) ;
	/*
	  #ifdef VERBOSE
	  debugout[MPI_rank] << "packing variable " << send_info[i].second[j].v
	  << endl ;
	  #endif
	*/

	/*
	  if((send_info[i].second[j].set).inSet(0)) {
	  Loci::debugout[Loci::MPI_rank] << "sending to processor " << send_info[i].first
	  << endl ;	
	  Loci::debugout[Loci::MPI_rank] << "packing variable " << send_info[i].second[j].v
	  << endl ;
	  }
	*/
	sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
      }
      
      warn(loc_pack != s_size[i]) ;
    }
    
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_PACKED,proc,1,MPI_COMM_WORLD) ;
    }
    
    if(nrecv > 0) { 
      int err = MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      /*
	#ifdef VERBOSE
	debugout[MPI_rank] << "unpacking from processor " <<recv_info[i].first
	<< endl ;
	#endif
      */
      for(int j=0;j<recv_info[i].second.size();++j) {
        storeRepP sp = facts.get_variable(recv_info[i].second[j].v) ;
        storeRepP sr = sp->new_store(entitySet(recv_info[i].second[j].seq)) ;
        /*
	  #ifdef VERBOSE
	  debugout[MPI_rank] << "unpacking variable " << recv_info[i].second[j].v
	  << endl ;
	  #endif
	*/
	
	if((entitySet(recv_info[i].second[j].seq)).inSet(0)) {
	  Loci::debugout[Loci::MPI_rank] << "  unpacking from processor " <<recv_info[i].first
					 << endl ;
	  Loci::debugout[Loci::MPI_rank] << "  unpacking variable " << recv_info[i].second[j].v
					 << endl ;
	  Loci::debugout[Loci::MPI_rank]  << "sp before join = " << endl ;
	  sp->Print( Loci::debugout[Loci::MPI_rank]) ; 
	}
	
        sr->unpack(recv_ptr[i], loc_unpack, r_size[i],
                   recv_info[i].second[j].seq) ;
	
	if((entitySet(recv_info[i].second[j].seq)).inSet(0)) {
	  Loci::debugout[Loci::MPI_rank]  << " sr   = " << endl ;
	  sr->Print( Loci::debugout[Loci::MPI_rank]) ;
	}
	
        CPTR<joiner> op = join_op->clone() ;
        op->SetArgs(sp,sr) ;
        op->Join(recv_info[i].second[j].seq) ;
	if((entitySet(recv_info[i].second[j].seq)).inSet(0)) {
	  Loci::debugout[Loci::MPI_rank] << " sp after join  = " << endl ;
	  sp->Print( Loci::debugout[Loci::MPI_rank]) ;
	}
      }
      warn(loc_unpack != r_size[i]) ;
    }
    
    delete [] status ;
    delete [] request ;
    delete [] send_ptr[0] ;
    delete [] send_ptr ;
    delete [] s_size ;
    delete [] recv_ptr[0] ;
    delete [] recv_ptr ;
    delete [] r_size ;
  }
  
  void execute_comm_reduce::Print(ostream &s) const {
    if(send_info.size()+recv_info.size() > 0) {
      s << "reduction block {" << endl ;
      if(send_info.size() > 0) {
        s << "Send:" << endl ;
        for(int i=0;i<send_info.size();++i) {
          for(int j=0;j<send_info[i].second.size();++j)
            s << "(" << send_info[i].second[j].v << "," << send_info[i].second[j].set << ") " ;
          s << " to " << send_info[i].first << endl ;
        }
      }
      if(recv_info.size() > 0) {
        s << "Recv:" << endl ;
        for(int i=0;i<recv_info.size();++i) {
          for(int j=0;j<recv_info[i].second.size();++j)
            s << "(" << recv_info[i].second[j].v << "," << recv_info[i].second[j].seq << ") " ;
          s << " from " << recv_info[i].first << endl ;
        }
      }
      s << "}" << endl ;
    }
  }
  
  executeP reduce_store_compiler::create_execution_schedule(fact_db &facts) {
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ;
      
      el->append_list(new execute_comm_reduce(rlist, facts, join_op)) ;
      el->append_list(new execute_comm(clist, facts)) ;
      ostringstream oss ;
      oss << "reduce store " << reduce_var ;
      el->append_list(new execute_msg(oss.str())) ;

      return executeP(el) ;
    }

    ostringstream oss ;
    oss << "reduce store " << reduce_var ;
    return executeP(new execute_msg(oss.str())) ;
  }
  
}

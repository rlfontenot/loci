#include "comp_tools.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;


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
  
  
  void apply_compiler::set_var_existence(fact_db &facts) {
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
    
    if(facts.isDistributed()) {
      constraint my_entities ;
      my_entities = facts.get_variable("my_entities") ;
      const rule_impl::info &finfo = apply.get_info().desc ;
      set<vmap_info>::const_iterator si ;
      entitySet compute ;
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
	compute |= vmap_target_requests(*si,tvarmap,facts) ;
      }
      output_mapping = false ;
      for(si=finfo.targets.begin();si!=finfo.targets.end(); ++si) {
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
      entitySet srcs = my_entities ;
      entitySet cnstrnts = my_entities ;
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si)
	srcs &= vmap_source_exist(*si,facts) ;
      for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
	cnstrnts &= vmap_source_exist(*si,facts) ;
      if(finfo.constraints.begin() != finfo.constraints.end())
	if((srcs & cnstrnts) != cnstrnts) {
	  cerr << "Warning, reduction rule:" << apply
	       << "cannot supply all entities of constraint" << endl ;
	  cerr << "constraints = " <<cnstrnts << endl ;
	  entitySet sac = srcs & cnstrnts ;
	  cerr << "srcs & constraints = " << sac << endl ;
	  //      exit(-1) ;
	}
      srcs &= cnstrnts ;
      
      // now trim compute to what can be computed.
      compute &= srcs ;
      exec_seq = compute ;
      
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
	entitySet requests = vmap_source_requests(*si,facts,compute) ;
	variableSet::const_iterator vi ;
	for(vi=si->var.begin();vi!=si->var.end();++vi)
	  facts.variable_request(*vi,requests) ;
      }
    }
    else {
      const rule_impl::info &finfo = apply.get_info().desc ;
      set<vmap_info>::const_iterator si ;
      entitySet compute ;
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
	compute |= vmap_target_requests(*si,tvarmap,facts) ;
      }
      
      output_mapping = false ;
      for(si=finfo.targets.begin();si!=finfo.targets.end(); ++si) {
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
      entitySet cnstrnts = ~EMPTY ;
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si)
	srcs &= vmap_source_exist(*si,facts) ;
      for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
	cnstrnts &= vmap_source_exist(*si,facts) ;
      if(finfo.constraints.begin() != finfo.constraints.end())
	if((srcs & cnstrnts) != cnstrnts) {
	  cerr << "Warning, reduction rule:" << apply
	       << "cannot supply all entities of constraint" << endl ;
	  cerr << "constraints = " <<cnstrnts << endl ;
	  entitySet sac = srcs & cnstrnts ;
	  cerr << "srcs & constraints = " << sac << endl ;
	  //      exit(-1) ;
	}
      srcs &= cnstrnts ;
      
      // now trim compute to what can be computed.
      compute &= srcs ;
      
      exec_seq = compute ;
      
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
	entitySet requests = vmap_source_requests(*si,facts,compute) ;
	variableSet::const_iterator vi ;
	for(vi=si->var.begin();vi!=si->var.end();++vi)
	  facts.variable_request(*vi,requests) ;
      }
    }
    
#ifdef VERBOSE
      cout << "rule " << apply << " computes over " << compute << endl ;
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

      const rule_impl::info &finfo = apply.get_info().desc ;
      execute_par *ep = new execute_par ;
      entitySet apply_domain,all_contexts ;
      vector<entitySet> shards, shard_domains ;
      for(int i=0;i<partition.size();++i) {
        fatal(finfo.targets.size() != 1) ;
        entitySet context = partition[i] ;
        entitySet pdom = vmap_target_exist(*finfo.targets.begin(),facts,context) ;
      
        entitySet rem = pdom & apply_domain ;
        if(rem != EMPTY) {
          entitySet compute = rem ;
          const vmap_info &vmi = *finfo.targets.begin() ;
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
    sp = global_join_op->getTargetRep() ;
    tp = global_join_op->getTargetRep() ;
    sp->unpack(send_ptr, loc_send, *size, seq) ;
    tp->unpack(result_ptr, loc_result, *size, seq) ;
    global_join_op->SetArgs(tp, sp) ;
    global_join_op->Join(seq) ;
    loc_result = 0 ;
    loc_send = 0 ;
    tp->pack(result_ptr, loc_result, *size, e) ;
    sp->pack(send_ptr, loc_send, *size, e) ;
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
    //cout << "performing param reduction " << endl ;
  }
}

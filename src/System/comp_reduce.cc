#include "comp_tools.h"
#include "distribute.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;

using std::ostream ;
using std::endl ;
using std::ostringstream ;

#include <Tools/hash_map.h>

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
      for(unsigned int i=0;i<var_vec.size();++i) {
        join_ops.push_back(jo->clone()) ;
        join_ops[i]->SetArgs(joiner_store,var_vec[i]) ;
      }
    }
    virtual void execute(fact_db &facts) ;
    virtual void Print(ostream &s) const ;
  } ;

  void joiner_oper::execute(fact_db &facts) {
    for(unsigned int i=0;i<var_vec.size();++i) 
      join_ops[i]->Join(sequence(partition[i])) ;
    //joiner_op->Join(joiner_store,var_vec[i],sequence(partition[i])) ;
  }

  void joiner_oper::Print(ostream &s) const {
    s << "reducing thread results for variable " << joiner_var << endl ;
    s << "reducing partitions = " << endl ;
    for(unsigned int i=0;i<var_vec.size();++i)
      s << "p["<<i<< "]="<<partition[i]<<endl ;
  }
  
  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var, sched_db &scheds) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      if(*vi != reduce_var)
        sources &= scheds.variable_existence(*vi) ;
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
        working &= scheds.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }
  
  void apply_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
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
        sources &= vmap_source_exist_apply(*si,facts,reduce_var, scheds) ;
      } 
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      sources &= constraints ;
      
      entitySet context = sources & constraints ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        entitySet targets = vmap_target_exist(*si,facts,context, scheds) ;
        const variableSet &tvars = si->var ;
        variableSet::const_iterator vi ;
	for(vi=tvars.begin();vi!=tvars.end();++vi) {
#ifdef VERBOSE
	debugout << "shadow is " << targets << endl ;
	debugout << "shadow not owned is "
			   << targets - d->my_entities << endl
			   << "variable is " << *vi << endl ;
#endif
	scheds.variable_shadow(*vi,targets) ;
      }
    }
  }
} 
  
  
  void apply_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    
#ifdef VERBOSE
    debugout << "in process_var_requests" << endl ;
#endif 
    vdefmap tvarmap ;
    variableSet targets = apply.targets() ;
    variableSet sources = apply.sources() ;
    
    fatal(targets.size() != 1) ;
    variable tvar = *(targets.begin()) ;
    
    if(facts.get_variable(tvar)->RepType() == Loci::PARAMETER) 
      tvarmap[tvar] = scheds.variable_existence(tvar) ;
    else
      tvarmap[tvar] = scheds.get_variable_request(unit_tag,tvar) ;
    
    const rule_impl::info &rinfo = apply.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    entitySet compute ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      compute |= vmap_target_requests(*si,tvarmap,facts, scheds) ;
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
          FATAL(!scheds.is_a_Map(*vi)) ;
          working |= scheds.image(*vi,comp) ;
        }
        comp = working ;
      }
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        if((comp - scheds.variable_existence(*vi)) != EMPTY) {
          cerr << "ERROR: Apply rule " << apply <<  endl
               << " output mapping forces application to entities where unit does not exist." << endl ;
          cerr << "error occurs for entities " <<
            entitySet(comp-scheds.variable_existence(*vi)) << endl ;
          cerr << "error occurs when applying to variable " << *vi << endl;
          cerr << "error is not recoverable, terminating scheduling process"
               << endl ;
          exit(-1) ;
        }
        scheds.variable_request(*vi,comp) ;
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
      srcs &= vmap_source_exist(*si,facts, scheds) ;
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      cnstrnts &= vmap_source_exist(*si,facts, scheds) ;
    if(rinfo.constraints.begin() != rinfo.constraints.end())
      if((srcs & cnstrnts) != cnstrnts) {
        cerr << "Warning, reduction rule:" << apply
             << "cannot supply all entities of constraint" << endl ;
        cerr << "constraints = " <<cnstrnts << endl ;
        entitySet sac = srcs & cnstrnts ;
        cerr << "srcs & constraints = " << sac << endl ;
        //      exit(-1) ;
        for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
          entitySet sources = vmap_source_exist(*si,facts, scheds) ;
          sources &= my_entities ;
          if((sources & cnstrnts) != cnstrnts) {
            cerr << "sources & constraints != constraints for input"
                 << endl
                 << sources  << " -- " << *si << endl ;
            
            if(si->mapping.size() > 0) {
              entitySet working = cnstrnts ;
              for(unsigned int i=0;i<si->mapping.size();++i) {
                entitySet images ;
                variableSet::const_iterator vi ;
                for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                  images |= scheds.image(*vi,working) ;
                working = images ; 
              }
              variableSet::const_iterator vi ;
              for(vi=si->var.begin();vi!=si->var.end();++vi) {
                entitySet exist = scheds.variable_existence(*vi) ;
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
    
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      entitySet requests = vmap_source_requests(*si,facts,compute, scheds) ;
      variableSet::const_iterator vi ;
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        variable v = *vi ;
        scheds.variable_request(v,requests) ; 
#ifdef VERBOSE
	debugout << "rule " << apply << " requesting variable "
		 << v << " for entities " << requests << endl ;
#endif
      }
    }
    
#ifdef VERBOSE
    debugout << "rule " << apply << " computes over " << compute << endl ;
#endif
  }
  
  executeP apply_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifndef DEBUG
    if(exec_seq.size() == 0)
      return executeP(0) ;
#endif
    CPTR<execute_list> el = new execute_list ;
    if(num_threads == 1 || !apply.get_info().rule_impl->thread_rule() ||
       exec_seq.size() < num_threads*30 ) {
      el->append_list(new execute_rule(apply,sequence(exec_seq),facts, scheds)) ;
    } else if(!apply.get_info().output_is_parameter &&!output_mapping) {
      execute_par *ep = new execute_par ;
      parallel_schedule(ep,exec_seq,apply,facts, scheds) ;
      el->append_list(ep)  ;
    } else if(apply.get_info().output_is_parameter) {
      variableSet target = apply.targets() ;
      fatal(target.size() != 1) ;
      variable v = *(target.begin()) ;
      storeRepP sp = facts.get_variable(v) ;
      vector<entitySet> partition = partition_set(exec_seq,num_threads) ;
      vector<storeRepP> var_vec ;

      execute_par *ep = new execute_par ;
      for(unsigned int i=0;i<partition.size();++i) {
        storeRepP rp = sp->new_store(EMPTY) ;
        rp->allocate(partition[i]) ;
        var_vec.push_back(rp) ;
        execute_sequence *es = new execute_sequence ;
        es->append_list(new execute_rule(unit_tag,sequence(partition[i]),
                                         facts,v,rp, scheds)) ;
        es->append_list(new execute_rule(apply,sequence(partition[i]),
                                         facts,v,rp, scheds)) ;
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
      for(unsigned int i=0;i<partition.size();++i) {
        fatal(rinfo.targets.size() != 1) ;
        entitySet context = partition[i] ;
        entitySet pdom = vmap_target_exist(*rinfo.targets.begin(),facts,context, scheds) ;
        entitySet rem = pdom & apply_domain ;
        if(rem != EMPTY) {
          entitySet compute = rem ;
          const vmap_info &vmi = *rinfo.targets.begin() ;
          vector<variableSet>::const_reverse_iterator mi ;
          for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
            entitySet working = EMPTY ;
            variableSet::const_iterator vi;
            for(vi=mi->begin();vi!=mi->end();++vi) {
              FATAL(!scheds.is_a_Map(*vi)) ;
              working |= scheds.preimage(*vi,compute).second ;
            }
            compute = working ;
          }
          compute &= partition[i] ;
          shards.push_back(compute) ;
          entitySet sdom = vmap_target_exist(vmi,facts,compute, scheds) ;
          shard_domains.push_back(sdom) ;
          context &= ~compute ;
        }
        apply_domain |= pdom ;
        all_contexts |= partition[i] ;
        ep->append_list(new execute_rule(apply,sequence(context),facts, scheds)) ;
      }
      if(shards.size() == 0) {
        el->append_list(ep) ;
        el->append_list(new execute_thread_sync) ;
      } else {
        ep->append_list(new execute_sequence) ;
        bool disjoint = true ;
        entitySet dom_tot ;
        for(unsigned int i=0;i<shards.size();++i) {
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
      
        for(unsigned int i=0;i<shards.size();++i) {
          storeRepP rp = sp->new_store(EMPTY) ;
          rp->allocate(shard_domains[i]) ;

          var_vec.push_back(rp) ;
          execute_sequence *es = new execute_sequence ;
          es->append_list(new execute_rule(unit_tag,sequence(shard_domains[i]),
                                           facts,v,rp, scheds)) ;
          es->append_list(new execute_rule(apply,sequence(shards[i]),
                                           facts,v,rp, scheds)) ;
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
          for(unsigned int i=0;i<shard_domains.size();++i) {
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
          for(unsigned int i=0;i<shard_domains.size();++i) {
            execute_par *epj = new execute_par ;
            vector<entitySet> decompose = partition_set(shard_domains[i],num_threads) ;
            vector<storeRepP> vv ;
            vv.push_back(var_vec[i]) ;
            for(unsigned int j=0;j<decompose.size();++j) {
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
  
  execute_param_red::execute_param_red(variable red, rule unit, CPTR<joiner> j_op) {
    reduce_var = red ;
    unit_rule = unit ;
    join_op = j_op ;
    global_join_op = join_op ;
    MPI_Op_create(&create_user_function, 0, &create_join_op) ;
  }
  execute_param_red::~execute_param_red() {
    MPI_Op_free(&create_join_op) ;
  }
  void execute_param_red::execute(fact_db &facts) {
    void *send_ptr; 
    void *result_ptr ;
    int size ;
    int loc = 0 , loc_result = 0;
    storeRepP sp, tp ;
    entitySet e ;
    sequence seq ;
    sp = facts.get_variable(reduce_var) ;
    size = sp->pack_size(e) ;
    send_ptr = new unsigned char[size] ;
    result_ptr = new unsigned char[size] ;
    sp->pack(send_ptr, loc, size, e) ;
    MPI_Allreduce(send_ptr, result_ptr, size, MPI_PACKED, create_join_op, MPI_COMM_WORLD) ;
    sp->unpack(result_ptr, loc_result, size, seq) ;
    delete [] send_ptr ;
    delete [] result_ptr ;
  }
  
  void execute_param_red::Print(ostream &s) const {
    s << "param reduction on " << reduce_var << endl ;
  }
  void reduce_param_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
    
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      entitySet targets ;
      targets = scheds.get_existential_info(reduce_var, unit_rule) ;
      targets += send_entitySet(targets, facts) ;
      targets &= d->my_entities ;
      targets += fill_entitySet(targets, facts) ;
      scheds.set_existential_info(reduce_var,unit_rule,targets) ;
    }
  }
  
  void reduce_param_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      entitySet requests = scheds.get_variable_requests(reduce_var) ;
      requests += send_entitySet(requests, facts) ;
      scheds.variable_request(reduce_var,requests) ;
    }
  } 
  
  executeP reduce_param_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
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
  
  void reduce_store_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      entitySet targets = scheds.get_existential_info(reduce_var, unit_rule) ;
      targets += send_entitySet(targets, facts) ;
      targets &= d->my_entities ;
      targets += fill_entitySet(targets, facts) ;
      scheds.set_existential_info(reduce_var,unit_rule,targets) ;
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
  
  void reduce_store_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      variableSet vars ;
      vars += reduce_var ;
      list<comm_info> request_comm = barrier_process_rule_requests(vars,facts, scheds) ;
      entitySet requests = scheds.get_variable_requests(reduce_var) ;
      entitySet shadow = scheds.get_variable_shadow(reduce_var) ;
      shadow &= requests ;
      list<comm_info> slist ;
      entitySet response = send_requests(shadow, reduce_var,facts,slist) ;
      swap_send_recv(slist) ;
      rlist = sort_comm(slist,facts) ;
      clist = sort_comm(request_comm,facts) ;
      
#ifdef VERBOSE
      if(shadow != EMPTY) {
	debugout << "shadow = " << shadow << endl ;
	shadow -= d->my_entities ;
	debugout << "shadow/my_entites = " << shadow << endl ;
      }
#endif
    }
  }
  
  class execute_comm_reduce : public execute_modules {
    vector<pair<int,vector<send_var_info> > > send_info ;
    vector<pair<int,vector<recv_var_info> > > recv_info ;
    std::vector<std::vector<storeRepP> > send_vars ; 
    std::vector<std::vector<storeRepP> > recv_vars ; 
    CPTR<joiner> join_op ;
    int *maxr_size, *maxs_size, *r_size, *s_size, *recv_sizes ;
    unsigned char **recv_ptr , **send_ptr ;
    MPI_Request *request;
    MPI_Status *status ;
  public:
    execute_comm_reduce(list<comm_info> &plist, fact_db &facts,
			CPTR<joiner> jop) ;
    ~execute_comm_reduce() ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ; 

execute_comm_reduce::execute_comm_reduce(list<comm_info> &plist,
					 fact_db &facts,
					 CPTR<joiner> jop) {
  join_op = jop ;
  HASH_MAP(int,vector<send_var_info> ) send_data ;
  HASH_MAP(int,vector<recv_var_info> ) recv_data ;
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
      send_vars.push_back(std::vector<storeRepP>()) ; 
      for(unsigned int i=0;i<send_data[*ii].size();++i) 
        send_vars.back().push_back(facts.get_variable(send_data[*ii][i].v)) ; 
    }
    
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
      recv_vars.push_back(std::vector<storeRepP>()) ; 
      for(unsigned int i=0;i<recv_data[*ii].size();++i) 
        recv_vars.back().push_back(facts.get_variable(recv_data[*ii][i].v)) ; 
    }
    
    int nsend = send_info.size() ;
    int nrecv = recv_info.size() ;
    r_size = new int[nrecv] ;
    recv_sizes = new int[nrecv] ; 
    maxr_size = new int[nrecv] ; 
    maxs_size = new int[nsend] ;
    s_size = new int[nsend] ;
    for(int i = 0; i < nrecv; ++i) {
      r_size[i] = 0 ;
      maxr_size[i] = sizeof(int) ;
    }
    for(int i = 0; i < nsend; ++i) {
      maxs_size[i] = sizeof(int) ;
      s_size[i] = 0 ;
    }
    
    recv_ptr = new unsigned char*[nrecv] ;
    send_ptr = new unsigned char*[nsend] ;
    request =  new MPI_Request[nrecv] ;
    status =  new MPI_Status[nrecv] ;
}
  execute_comm_reduce::~execute_comm_reduce() {
    delete [] maxr_size ;
    delete [] maxs_size ;
    delete [] recv_sizes ; 
    delete [] r_size ;
    delete [] s_size ;
    delete [] recv_ptr ;
    delete [] send_ptr ;
    delete [] request ;
    delete [] status ;
  }
  static unsigned char *recv_ptr_buf = 0; 
  static int recv_ptr_buf_size = 0; 
  static unsigned char *send_ptr_buf = 0 ; 
  static int send_ptr_buf_size = 0 ; 
  
  
  void execute_comm_reduce::execute(fact_db  &facts) {
    const int nrecv = recv_info.size() ;
    int resend_size = 0, rerecv_size = 0 ;
    std::vector<int> send_index ;
    std::vector<int> recv_index ; 
    int total_size = 0 ;
    MPI_Request *re_request ;
    MPI_Status *re_status ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = maxr_size[i] ;
      total_size += maxr_size[i] ;
    }
    /*
      #ifdef DEBUG
      entitySet rem = entitySet((recv_info[i].second[j].seq)) - sp->domain() ;
      if(rem != EMPTY)
      debugout << "variable " << recv_info[i].second[j].v << " not allocated, but recving for reduce entities " << rem << endl ;
      #endif
    */
    if(recv_ptr_buf_size < total_size) { 
      if(recv_ptr_buf) 
        delete[] recv_ptr_buf ; 
      recv_ptr_buf = new unsigned char[total_size] ; 
      recv_ptr_buf_size = total_size ; 
    } 
    recv_ptr[0] = recv_ptr_buf ; 
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
    
    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
    }
    total_size = 0 ;
    const int nsend = send_info.size() ;
    entitySet resend_procs, rerecv_procs ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(unsigned int j=0;j<send_info[i].second.size();++j) {
        //facts.get_variable(send_info[i].second[j].v) ;
	storeRepP sp = send_vars[i][j] ;

        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
      }
      if((s_size[i] > maxs_size[i]) || ( s_size[i] == sizeof(int))) {
	if(s_size[i] > maxs_size[i])
	  maxs_size[i] = s_size[i] ;
	int proc = send_info[i].first ;
	s_size[i] = sizeof(int) ;
	resend_procs += proc ;
	send_index.push_back(i) ;
      }
      total_size += maxs_size[i] ;
    }
    if(send_ptr_buf_size < total_size) { 
      if(send_ptr_buf) 
        delete[] send_ptr_buf ; 
      send_ptr_buf = new unsigned char[total_size] ; 
      send_ptr_buf_size = total_size ; 
    } 
    send_ptr[0] = send_ptr_buf ; 
    for(int i = 1; i < nsend; i++)
      send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
    // Pack the buffer for sending 
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      if(!resend_procs.inSet(send_info[i].first)) {
	for(unsigned int j=0;j<send_info[i].second.size();++j) {
	  storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[i].second[j].v) ;
	  sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
	}
      }
      else
	MPI_Pack(&maxs_size[i], sizeof(int), MPI_BYTE, send_ptr[i], s_size[i], &loc_pack, MPI_COMM_WORLD) ; 
      
    }
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_PACKED,proc,1,MPI_COMM_WORLD) ;
    }
    if(nrecv > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int i = 0 ; i < nrecv; i++) {
        int rcv_sizes ;
	MPI_Get_count(&status[i], MPI_BYTE, &rcv_sizes) ;  
	if(rcv_sizes == sizeof(int)) {
	  rerecv_procs += recv_info[i].first ;
	  recv_index.push_back(i) ;
	}
      }
    }
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      if(rerecv_procs.inSet(recv_info[i].first)) {
	int temp ;
	MPI_Unpack(recv_ptr[i], r_size[i], &loc_unpack, &temp, sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
	if(temp > maxr_size[i])
	  maxr_size[i] = temp ;
      }
      else
	for(unsigned int j=0;j<recv_info[i].second.size();++j) {
	  storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[i].second[j].v) ;
	  storeRepP sr = sp->new_store(EMPTY) ;
	  sr->allocate(entitySet(recv_info[i].second[j].seq)) ;
          
	  sr->unpack(recv_ptr[i], loc_unpack, r_size[i],
		     recv_info[i].second[j].seq) ;
	  CPTR<joiner> op = join_op->clone() ;
	  op->SetArgs(sp,sr) ;
	  op->Join(recv_info[i].second[j].seq) ;
	}
    }
    rerecv_size = rerecv_procs.size() ;
    resend_size = resend_procs.size() ;
    
    if(rerecv_size > 0) {
      re_request =  new MPI_Request[rerecv_size] ;
      re_status =  new MPI_Status[rerecv_size] ;
    }
    for(int i = 0; i < rerecv_size; i++) {
      int proc = recv_info[recv_index[i]].first ;
      recv_ptr[recv_index[i]] = new unsigned char[maxr_size[recv_index[i]]] ;
      MPI_Irecv(recv_ptr[recv_index[i]], maxr_size[recv_index[i]], MPI_PACKED, proc, 2, MPI_COMM_WORLD, &re_request[i]) ;
    }
    
    for(int i=0;i<resend_size;++i) {
      int loc_pack = 0 ;
      send_ptr[send_index[i]] = new unsigned char[maxs_size[send_index[i]]] ;
      for(unsigned int j=0;j<send_info[send_index[i]].second.size();++j) {
	storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[send_index[i]].second[j].v) ;
	sp->pack(send_ptr[send_index[i]], loc_pack,maxs_size[send_index[i]],send_info[send_index[i]].second[j].set);
      }
    }
    
    // Send Buffer
    for(int i=0;i<resend_size;++i) {
      int proc = send_info[send_index[i]].first ;
      MPI_Send(send_ptr[send_index[i]],maxs_size[send_index[i]],MPI_PACKED,proc,2,MPI_COMM_WORLD) ;
      delete [] send_ptr[send_index[i]] ;
    }
    if(rerecv_size > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(rerecv_size, re_request, re_status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    for(int i=0;i<rerecv_size;++i) {
      int loc_unpack = 0;
      for(unsigned int j=0;j<recv_info[recv_index[i]].second.size();++j) {
	storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[recv_index[i]].second[j].v) ;
	storeRepP sr = sp->new_store(EMPTY) ;
	sr->allocate(entitySet(recv_info[recv_index[i]].second[j].seq)) ;
	
	sr->unpack(recv_ptr[recv_index[i]], loc_unpack, maxr_size[recv_index[i]],
		   recv_info[recv_index[i]].second[j].seq) ;
	
	CPTR<joiner> op = join_op->clone() ;
	op->SetArgs(sp,sr) ;
	op->Join(recv_info[recv_index[i]].second[j].seq) ;
      }
      delete [] recv_ptr[recv_index[i]] ;
    }
    if(rerecv_size > 0) {
      delete [] re_status ;
      delete [] re_request ;
    }
  }
  
  void execute_comm_reduce::Print(ostream &s) const {
    int sz = 0 ;
    if(send_info.size()+recv_info.size() > 0) {
      s << "reduction block {" << endl ;
      if(send_info.size() > 0) {
        s << "Send:" << endl ;
        for(unsigned int i=0;i<send_info.size();++i) {
          for(unsigned int j=0;j<send_info[i].second.size();++j) {
            s << send_info[i].second[j].v << "  " ;
	    sz += (send_info[i].second[j].set).size() ;
	  }
	  s << " to " << send_info[i].first << endl ;
        }
	s << " Total entities sent = " << sz << endl ;	
      }
      sz = 0 ;
      if(recv_info.size() > 0) {
        s << "Recv:" << endl ;
        for(unsigned int i=0;i<recv_info.size();++i) {
          for(unsigned int j=0;j<recv_info[i].second.size();++j) {
            s << recv_info[i].second[j].v << "  " ;
	    sz += (recv_info[i].second[j].seq).size() ;
	  }
          s << " from " << recv_info[i].first << endl ;
        }
	s << " Total entities recieved = " << sz << endl ;	
      }
      s << "}" << endl ;
    }
  }
  
  executeP reduce_store_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
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

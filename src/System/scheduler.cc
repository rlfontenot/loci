#include "sched_tools.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

//#define VERBOSE
namespace Loci {

  // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet> schedule_dag(const digraph &g,
                                          digraph::vertexSet start_vertices = EMPTY,
                                          digraph::vertexSet only_vertices =
                                          interval(UNIVERSE_MIN,UNIVERSE_MAX)) {
    digraph gt = g.transpose() ;

    vector<digraph::vertexSet> schedule ; 
    // First schedule any vertices that have no edges leading into them and have
    // not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    if(working != EMPTY)
      schedule.push_back(working) ;
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = start_vertices + working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni)
        new_vertices += g[*ni] ;
      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      // We only schedule vertices that are also in the only_vertices set
      working = new_vertices & only_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all vertices leading
      // to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni) 
        if((gt[*ni] & visited_vertices) == gt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // and these new vertices to the schedule
      if(new_vertices != EMPTY)
        schedule.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
    return schedule ;
  }


  void existential_rule_analysis(rule f, fact_db &facts) {

    FATAL(f.type() == rule::INTERNAL) ;

    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet targets = EMPTY ;
    const rule_impl::info &finfo = f.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
      sources &= vmap_source_exist(*si,facts) ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
      constraints &= vmap_source_exist(*si,facts) ;
    if(finfo.constraints.begin() != finfo.constraints.end())
      if((sources & constraints) != constraints) {
        cerr << "Warning, rule " << f <<
          " cannot supply all entities of constraint" << endl ;
        cerr << "constraints = " << constraints ;
        cerr << "sources & constraints = " << (sources & constraints) << endl ;
        //      exit(-1) ;
      }
    sources &= constraints ;
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      targets = vmap_target_exist(*si,facts,sources) ;
      const variableSet &tvars = si->var ;
      variableSet::const_iterator vi ;
      for(vi=tvars.begin();vi!=tvars.end();++vi)
        facts.set_existential_info(*vi,f,targets) ;
    }
#ifdef VERBOSE
    cout << "rule " << f << " generating " << targets << endl ;
#endif
  }
  typedef map<variable,entitySet> vdefmap ;

  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts) {
    variableSet::const_iterator vi ;
    entitySet targets ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      FATAL(tvarmap.find(*vi) == tvarmap.end()) ;
      targets |= tvarmap.find(*vi)->second ;
    }
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,targets) ;
  
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.preimage(*vi,targets).second ;
      }
      targets = working ;
    }
    return targets ;
  }

  void vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                            entitySet compute) {
    vector<variableSet>::const_iterator mi ;
    variableSet::const_iterator vi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      entitySet working ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.image(*vi,compute) ;
      }
      compute = working ;
    }
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,compute) ;
  }

  entitySet process_rule_requests(rule f, fact_db &facts) {

    FATAL(f.type() == rule::INTERNAL) ;

    variableSet targets = f.targets() ;
  
    variableSet::const_iterator vi ;
    vdefmap tvarmap ;
    for(vi=targets.begin();vi!=targets.end();++vi) {
      if(vi->get_info().name == string("OUTPUT")) 
        facts.variable_request(*vi,facts.variable_existence(*vi)) ;
      tvarmap[*vi] = facts.get_variable_request(f,*vi) ;
    }

    const rule_impl::info &finfo = f.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    entitySet compute,isect = ~EMPTY ;
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
      compute |= tmp ;
      isect &= tmp ;
    }

    if(isect != compute) {
      entitySet working = compute ;
      vector<variableSet>::const_reverse_iterator mi ;
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
        for(mi=si->mapping.rbegin();mi!=si->mapping.rend();++mi) {
          entitySet tmp ;
          for(vi=mi->begin();vi!=mi->end();++vi)
            tmp |= facts.image(*vi,working) ;
          working = tmp ;
        }
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          facts.variable_request(*vi,working) ;
          //#define UNDERSTAND_THIS
#ifdef UNDERSTAND_THIS
          WARN(facts.get_variable_request(f,*vi) != working) ;
          if(facts.get_variable_request(f,*vi) != working) {
            cerr << "f = " << f << endl ;
            cerr << "*vi = " << *vi << endl ;
            cerr << "facts.get_variable_request(f,*vi) = "
                 << facts.get_variable_request(f,*vi)
                 << ", working = " << working << endl ;
          }
#endif
        }
      }
    }

    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
      vmap_source_requests(*si,facts,compute) ;

    for(vi=finfo.conditionals.begin();vi!=finfo.conditionals.end();++vi)
      facts.variable_request(*vi,compute) ;

#ifdef VERBOSE
    cout << "rule " << f << " computes over " << compute << endl ;
#endif
    return compute ;
  }



  dag_compiler::dag_compiler(decompose_graph *gr, digraph gin) {
    graph_ref = gr ;
    dag = gin ;
    dag_sched = schedule_dag(dag) ;
    extract_rule_sequence(rule_schedule,dag_sched) ;
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(int i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
    if(allvertices != dag.get_all_vertices()) {
      digraph::vertexSet leftout = allvertices ^ dag.get_all_vertices() ;
      cerr << "leftout rules= " << extract_rules(leftout) << endl ;
      cerr << "leftout vars = " << extract_vars(leftout) << endl ;
    }
#endif
  }


  class execute_rule : public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    bool do_run ;
  public:
    execute_rule(rule fi, sequence seq, fact_db &facts) ;
    execute_rule(rule fi, sequence seq, fact_db &facts, variable v, const storeRepP &p) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

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
    s << rule_tag << " over sequence " << exec_seq << endl ;
  }

  class execute_loop : public execute_modules {
    executeP collapse, advance ;
    variable cvar ;
    time_ident tlevel ;
  public:
    execute_loop(const variable &cv,
                 const executeP &col, const executeP &adv,
                 const time_ident &tl) :
      cvar(cv),collapse(col),advance(adv),tlevel(tl)
    { warn(col==0 || advance==0) ; control_thread = true ;}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  void execute_loop::execute(fact_db &facts) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;
    fact_db::time_infoP tinfo = facts.get_time_info(tlevel) ;
    facts.initialize_time(tinfo) ;
    for(;;) {
      collapse->execute(facts) ;
      if(*test) {
        facts.close_time(tinfo) ;
        return ;
      }
      advance->execute(facts) ;
      facts.advance_time(tinfo) ;
    }
  }

  void execute_loop::Print(ostream &s) const {
    s << "Perform loop for time level "<< tlevel << endl ;
    s << "--compute collapse rule, conditional on " << cvar << endl ;
    s << "--collapse iteration {" << tlevel << "}"<< endl ;
    collapse->Print(s) ;
    s << "--advance iteration {" << tlevel << "}" << endl ;
    advance->Print(s) ;
    s << "end of loop for time level " << tlevel << endl ;
  }


  
  class execute_conditional : public execute_modules {
    executeP conditional ;
    variable cvar ;
  public:
    execute_conditional(const executeP &cond, const variable &cv) :
      conditional(cond),cvar(cv)
    { warn(cond==0) ; control_thread = true; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  void execute_conditional::execute(fact_db &facts) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;

    if(*test) {
      conditional->execute(facts) ;
    }
  }

  void execute_conditional::Print(ostream &s) const {
    s << "--compute rule if conditional " << cvar << " true." << endl ;
    conditional->Print(s) ;
    s << "--end conditional" << endl ;
  }

  vector<entitySet> partition_set(const entitySet &s,int nthreads) {
    const int min = s.Min() ;
    const int max = s.Max() ;
    const int psize = s.size() ;

    const int div = psize/nthreads ;
    int rem = psize%nthreads ;
    vector<entitySet> partition ;
    
    for(int i=min;i<=max;) {
      int inval = div + ((rem>0)?1:0) ;
      rem-- ;
      entitySet sp = s & interval(i,i+inval-1) ;
      i+=inval ;
      while(sp.size() < inval && i<=max) {
        entitySet remain = s & interval(i,max) ;
        i=remain.Min() ;
        int remain_ival = inval - sp.size() ;
        sp += s & interval(i,i+remain_ival-1) ;
        i+=remain_ival ;
      }
      WARN(sp.size() > inval) ;
      partition.push_back(sp) ;
    }
    return partition ;
  }
  void parallel_schedule(execute_par *ep,const entitySet &exec_set,
                         const rule &impl, fact_db &facts) {
    vector<entitySet> par_set = partition_set(exec_set,num_threads) ;

    for(vector<entitySet>::const_iterator
          i=par_set.begin();i!=par_set.end();++i) {
      executeP execrule = new execute_rule(impl,sequence(*i),facts) ;
      ep->append_list(execrule) ;
    }
  }


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
    //      joiner_op->Join(joiner_store,var_vec[i],sequence(partition[i])) ;
  }

  void joiner_oper::Print(ostream &s) const {
    s << "reducing thread results for variable " << joiner_var << endl ;
    s << "reducing partitions = " << endl ;
    for(int i=0;i<var_vec.size();++i)
      s << "p["<<i<< "]="<<partition[i]<<endl ;
  }
                    
  
  void allocate_all_vars::execute(fact_db &facts) {
    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi,vii ;
    set<time_ident> time_set ;
    for(vi=vars.begin();vi!=vars.end();++vi)
      time_set.insert(vi->time()) ;
    set<time_ident>::const_iterator si ;
    for(si=time_set.begin();si!=time_set.end();++si) {
      fact_db::time_infoP tip = facts.get_time_info(*si) ;
      if(!tip->rotate_lists.empty()) {
        for(list<list<variable> >::const_iterator llv = tip->rotate_lists.begin();
            llv != tip->rotate_lists.end();++llv) {
          entitySet time_space ;
          for(list<variable>::const_iterator lv=llv->begin();lv!=llv->end();++lv) {
            variableSet aliases = facts.get_aliases(*lv) ;
            for(vii=aliases.begin();vii!=aliases.end();++vii)
              time_space += facts.get_variable_requests(*vii) ;
          }
          for(list<variable>::const_iterator lv=llv->begin();lv!=llv->end();++lv) {
            facts.variable_request(*lv,time_space) ;
          }
        }
      }
    }
  
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->domain() == EMPTY) {
        variableSet aliases = facts.get_aliases(*vi) ;
        entitySet all_requests ;
        for(vii=aliases.begin();vii!=aliases.end();++vii)
          all_requests += facts.get_variable_requests(*vii) ;
        srp->allocate(all_requests) ;
      }
    }
  }

  void allocate_all_vars::Print(ostream &s) const {
    s << "allocate all variables" << endl ;
  }


  loop_compiler::loop_compiler(decompose_graph *gr, digraph gin) {
    graph_ref = gr ;
    dag = gin ;
    ruleSet loopset =
      extract_rules(gin.get_all_vertices() & gr->looping_rules) ;
    if(loopset.size() != 1) {
      cerr << "internal consistency error, loopset size != 1" << endl ;
      exit(-1) ;
    }

    tlevel = loopset.begin()->source_time() ;
    dag.remove_vertex((*loopset.begin()).ident()) ;
    digraph dagt = dag.transpose() ;
    ruleSet collapse_rules, all_rules ;
    variableSet collapse_vars ;
  
    all_rules = extract_rules(dag.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) 
      if(ri->target_time().before(ri->source_time()))
        collapse_rules += *ri ;
  
    if(collapse_rules.size() != 1 ||
       collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
      cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
           << endl << "error is not recoverable" << endl ;
      cerr << "rules = " << collapse_rules << endl ;
      exit(-1) ;
    }

    variableSet all_vars = extract_vars(dag.get_all_vertices()) ;
    for(variableSet::const_iterator vi=all_vars.begin();vi!=all_vars.end();++vi)
      if(vi->get_info().offset == 1)
        advance_vars += *vi ;
  
    collapse_vars = collapse_rules.begin()->targets() ;
    cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

    // create Output Variable
    output = variable(variable("OUTPUT"),tlevel) ;
    variableSet outputSet ;
    outputSet = interval(output.ident(),output.ident()) ;

    // Schedule part of graph that leads to collapse
    collapse_sched = schedule_dag(dag, EMPTY,visit_vertices(dagt,collapse_vars)) ;
    extract_rule_sequence(rule_schedule,collapse_sched) ;
    extract_rule_sequence(collapse,collapse_sched) ;

    // Schedule advance part of loop.  First try to schedule any output, then
    // schedule the advance
    digraph::vertexSet visited ;
    for(int i = 0;i<collapse_sched.size();++i)
      visited += collapse_sched[i] ;
    vector<digraph::vertexSet>
      dag_sched = schedule_dag(dag,visited, visit_vertices(dagt,outputSet)) ;
    if(dag_sched.size() == 0)
      output_present = false ;
    else 
      output_present = true ;

    for(int i=0;i<dag_sched.size();++i)
      advance_sched.push_back(dag_sched[i]) ;
  
    extract_rule_sequence(rule_schedule,dag_sched) ;
    extract_rule_sequence(advance,dag_sched) ;
    for(int i = 0;i<dag_sched.size();++i)
      visited += dag_sched[i] ;
    // now schedule everything that hasn't been scheduled
    dag_sched = schedule_dag(dag,visited) ;
    extract_rule_sequence(rule_schedule,dag_sched) ;
    extract_rule_sequence(advance,dag_sched) ;
    for(int i=0;i<dag_sched.size();++i)
      advance_sched.push_back(dag_sched[i]) ;
  
  
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices = visited ;
    for(int i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
    if(allvertices != dag.get_all_vertices()) {
      cerr << " rules NOT scheduled = " << endl
           << extract_rules(dag.get_all_vertices() - allvertices) << endl ;
      cerr << " variables NOT scheduled = "
           << extract_vars(dag.get_all_vertices() - allvertices) << endl ;
      vector<digraph::vertexSet> components =
        component_sort(dag).get_components() ;
      for(int i=0;i<components.size();++i) {
        if(components[i].size() > 1) {
          cerr << "reason: graph not a dag, strongly connected component found"
               << endl ;
          cerr << "component rules = " << endl ;
          cerr << extract_rules(components[i]) << endl ;
        }
      }
      exit(-1) ;
    }
#endif
  }

  void loop_compiler::set_var_existence(fact_db &facts) {
    for(int i=0;i<rule_schedule.size();++i) 
      calc(rule_schedule[i])->set_var_existence(facts) ;
  }

  void loop_compiler::process_var_requests(fact_db &facts) {
    variableSet var_requests = advance_vars ;
    variableSet::const_iterator vi ;
    if(output_present)
      var_requests += output ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = facts.variable_existence(*vi) ;
      facts.variable_request(*vi,vexist) ;
    }
    vector<rule>::reverse_iterator ri ;
    for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
      calc(*ri)->process_var_requests(facts) ;
  }

  executeP loop_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> col = new execute_list ;
    for(int i=0;i!=collapse_sched.size();++i) {
      ruleSet rules = extract_rules(collapse_sched[i]) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        col->append_list(calc(*ri)->create_execution_schedule(facts)) ;
      if(rules.size() > 0 && num_threads > 1)
        col->append_list(new execute_thread_sync) ;
    }
    CPTR<execute_list> adv = new execute_list ;
    for(int i=0;i!=advance_sched.size();++i) {
      ruleSet rules = extract_rules(advance_sched[i]) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        adv->append_list(calc(*ri)->create_execution_schedule(facts)) ;
      if(rules.size() > 0 && num_threads > 1)
        adv->append_list(new execute_thread_sync) ;
    }
    return new execute_loop(cond_var,executeP(col),executeP(adv),tlevel) ;
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
    
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
      vmap_source_requests(*si,facts,compute) ;
    
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

  void impl_recurse_compiler::set_var_existence(fact_db &facts) {
    variableSet::const_iterator vi ;
    variableSet tvars ;
    tvars = impl.targets() ;
    fcontrol &fctrl = control_set ;
    const rule_impl::info &finfo = impl.get_info().desc ;
    warn(impl.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;

    entitySet constraints = ~EMPTY ;
    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
      if((si->var & tvars) == EMPTY)
        sources &= vmap_source_exist(*si,facts) ;
      else {
        const int num_maps = si->mapping.size() ;
        if(num_maps != 0) {
          for(vi=si->mapping[0].begin();vi!=si->mapping[0].end();++vi) {
            sources &= facts.variable_existence(*vi) ;
          }
        } 
        
        tvars += si->var ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;
        
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!facts.is_a_Map(*(miv[j]))) ;
              MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;  
              minfo.mapvec.push_back(m) ;
              minfo.mapvar.push_back(*(miv[j])) ;
            }
            minfo.v = *vi ;
            
            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.recursion_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
    } 
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      int num_maps = si->mapping.size() ;
      vector<variableSet::const_iterator> miv(num_maps) ;
      for(int i=0;i<num_maps;++i)
        miv[i] = si->mapping[i].begin() ;
      
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        int i = 0 ;
        do {
          fcontrol::mapping_info minfo ;
          for(int j=0;j!=num_maps;++j) {
            FATAL(!facts.is_a_Map(*(miv[j]))) ;
            MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
            minfo.mapvec.push_back(m) ;
            minfo.mapvar.push_back(*(miv[j])) ;
          }
          minfo.v = *vi ;
          
          for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
            miv[i] = si->mapping[i].begin() ;
          }
          fctrl.target_maps.push_back(minfo) ;
        } while(i!=num_maps) ;
      }
    }
    fctrl.use_constraints = false ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
      warn((si->var & tvars) != EMPTY) ;
      constraints &= vmap_source_exist(*si,facts) ;
      fctrl.use_constraints = true ;
    }

    fctrl.nr_sources = sources ;

    fctrl.constraints = constraints ;
    warn(fctrl.recursion_maps.size() == 0) ;

    fatal(tvars.size() != 1 ||
          control_set.recursion_maps.size() != 1 ||
          control_set.target_maps.size() != 1 ||
          fctrl.use_constraints) ;
    
    fcontrol::mapping_info &rmap = fctrl.recursion_maps[0] ;
    fcontrol::mapping_info &tmap = fctrl.target_maps[0] ;
    vector<multiMap> read_maps,read_map_inv ;

    for(int j=0;j<rmap.mapvec.size();++j) {
      read_maps.push_back(rmap.mapvec[j]->get_map()) ;
    }
    for(int j=0;j<read_maps.size();++j) {
      read_map_inv.push_back(multiMap()) ;
    }
    variable rvar = *(tvars.begin()) ;
    entitySet sdelta = facts.variable_existence(rvar) ;
    entitySet domain = fctrl.nr_sources + sdelta ;

    for(int j=read_maps.size()-1;j>=0;--j) {
      entitySet newdomain = facts.preimage(rmap.mapvar[j],domain).first ;
#ifdef VERBOSE
      cout << "j = " << j << ", domain = " << domain << ", newdomain = "
           << newdomain << endl ;
#endif
      if(domain == ~EMPTY) {
        cerr << "problem in recurse compiler for rule = "<< impl << endl ;
        domain = read_maps[j].domain() ;
      }
      inverseMap(read_map_inv[j],read_maps[j],domain,newdomain) ;
      domain = newdomain ;
    }
  
    for(int j=rmap.mapvec.size()-1;j>=0;--j)
      sdelta = rmap.mapvec[j]->preimage(sdelta).first ;
    sdelta &= fctrl.nr_sources ;
    entitySet tdelta = sdelta ;
    for(int j=0;j<tmap.mapvec.size();++j)
      tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
    cout << "sdelta_init = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif

    if(num_threads > 1)
      par_schedule.push_back(sdelta) ;
    else
      fastseq += sequence(sdelta) ;
    
    entitySet generated = tdelta ;
    const entitySet nr_sources = fctrl.nr_sources ;
    store<bool> exists ;
    exists.allocate(nr_sources) ;
    for(entitySet::const_iterator
          ei=nr_sources.begin();ei!=nr_sources.end();++ei)
      exists[*ei] = false ;
    for(entitySet::const_iterator
          ei=tdelta.begin();ei!=tdelta.end();++ei)
      exists[*ei] = true ;
    vector<const int *> mi(read_maps.size()) ;
    vector<const int *> me(read_maps.size()) ;

    do {
      for(int j=read_map_inv.size()-1;j>=0;--j) {
        entitySet candidates ;
        const int *mi ;
        for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
#ifdef DEBUG
          if(!read_map_inv[j].domain().inSet(*di)) {
            cerr << " di = " << *di << ", j = " << j << endl ;
          } else
#endif
            //        fatal(!read_map_inv[j].domain().inSet(*di)) ;
            for(mi=read_map_inv[j].begin(*di); mi!=read_map_inv[j].end(*di);++mi) 
              candidates += *mi ;
        }
        sdelta = candidates ;
      }
      // we now have candidate sdeltas, check them to see if they
      // are satisfied.
      entitySet satisfied ;
      for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
        int c = *di ;
        mi[0] = read_maps[0].begin(c) ;
        me[0] = read_maps[0].end(c) ;
        int j=0 ;
        const int je=read_maps.size();
        bool chk = true ;
        // Check that all paths from candidate go to generated cells
        while(chk && j>=0) {
          c = *mi[j] ;
          j++ ;
          if(j==je) {
            chk = chk && exists[c] && nr_sources.inSet(c) ;
            do {
              j-- ;
            } while(j>=0 && (++mi[j]) == me[j]) ;
          } else {
            mi[j] = read_maps[j].begin(c) ;
            me[j] = read_maps[j].end(c) ;
          }
        
        }
        if(chk)
          satisfied += *di ;
      }
      sdelta = satisfied ;
      entitySet tdelta = sdelta ;
      for(int j=0;j<tmap.mapvec.size();++j)
        tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
      cout << "sdelta = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif
      if(num_threads>1)
        par_schedule.push_back(sdelta) ;
      else
        fastseq += sequence(sdelta) ;

      for(entitySet::const_iterator
            ei=tdelta.begin();ei!=tdelta.end();++ei) 
        exists[*ei] = true ;
    
    } while(sdelta != EMPTY) ;

    for(entitySet::const_iterator
          ei=nr_sources.begin();ei!=nr_sources.end();++ei)
      if(exists[*ei]) {
        const int start = *ei ;
        const int end = nr_sources.Max() ;
        int finish = start ;
        for(;*ei!=end && exists[*ei];++ei)
          finish = *ei ;
        if(*ei == end && exists[end])
          finish = end ;
        generated += interval(start,finish) ;
      }
  
  
    fctrl.generated[rvar] = generated ;
#ifdef VERBOSE
    cout << "recursive rule " << impl << " generating " << generated << endl ;
#endif
    
    for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
        mi!=fctrl.generated.end();++mi)
      facts.set_existential_info(mi->first,impl,mi->second) ;
    
  }
  
  void impl_recurse_compiler::process_var_requests(fact_db &facts) {
    process_rule_requests(impl,facts) ;
  }

  executeP impl_recurse_compiler::create_execution_schedule(fact_db &facts) {

    if(num_threads > 1) {
      CPTR<execute_list> el = new execute_list ;
      sequence seq ;
      for(int i=0;i<par_schedule.size();++i) {
        if(par_schedule[i].size() < num_threads*4) {
          seq += sequence(par_schedule[i]) ;
        } else {
          if(seq.size() > 1) {
            el->append_list(new execute_rule(impl,seq,facts)) ;
            el->append_list(new execute_thread_sync) ;
            seq = sequence() ;
          }
          execute_par *ep = new execute_par ;
          parallel_schedule(ep,par_schedule[i],impl,facts) ;
          el->append_list(ep) ;
          el->append_list(new execute_thread_sync) ;
        }
      }
      if(seq.size() > 1) {
        el->append_list(new execute_rule(impl,seq,facts)) ;
      }
      return executeP(el) ;
    }

    return new execute_rule(impl,fastseq,facts) ;

  }

  void recurse_compiler::set_var_existence(fact_db &facts) {
    control_set.clear() ;
  
    ruleSet::const_iterator fi ;
    variableSet::const_iterator vi ;
    variableSet tvars ;
    const ruleSet &fset = recurse_rules ;
    for(fi=fset.begin();fi!=fset.end();++fi) 
      tvars += fi->targets() ;
    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      const rule_impl::info &finfo = fi->get_info().desc ;
      warn(fi->type() == rule::INTERNAL) ;
      entitySet sources = ~EMPTY ;
      entitySet constraints = ~EMPTY ;
      set<vmap_info>::const_iterator si ;
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
        if((si->var & tvars) == EMPTY)
          sources &= vmap_source_exist(*si,facts) ;
        else {
          tvars += si->var ;
          int num_maps = si->mapping.size() ;
          vector<variableSet::const_iterator> miv(num_maps) ;
          for(int i=0;i<num_maps;++i)
            miv[i] = si->mapping[i].begin() ;
        
          for(vi=si->var.begin();vi!=si->var.end();++vi) {
            int i = 0 ;
            do {
              fcontrol::mapping_info minfo ;
              for(int j=0;j!=num_maps;++j) {
                FATAL(!facts.is_a_Map(*(miv[j]))) ;
                minfo.mapvec.push_back(
                                       MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
              }
              minfo.v = *vi ;
              
              for(i=0;(i!=num_maps)&&(++(miv[i]) == si->mapping[i].end());++i)
                miv[i] = si->mapping[i].begin() ;

              fctrl.recursion_maps.push_back(minfo) ;
            } while(i!=num_maps) ;
          }
        }
      } 
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
        int num_maps = si->mapping.size() ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;

        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!facts.is_a_Map(*(miv[j]))) ;
              minfo.mapvec.push_back(
                                     MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
            }
            minfo.v = *vi ;
          
            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.target_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
      fctrl.use_constraints = false ;
      for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
        warn((si->var & tvars) != EMPTY) ;
        constraints &= vmap_source_exist(*si,facts) ;
        fctrl.use_constraints = true ;
      }
      fctrl.nr_sources = sources ;
      fctrl.constraints = constraints ;
      warn(fctrl.recursion_maps.size() == 0) ;
    }
  

    map<variable,entitySet> tvar_computed,tvar_update ;
    for(vi=tvars.begin();vi!=tvars.end();++vi) {
      tvar_computed[*vi] = facts.variable_existence(*vi) ;
      tvar_update[*vi] = tvar_computed[*vi] ;
    }
    map<variable,entitySet>::iterator mi ;
    for(mi=tvar_update.begin();mi!=tvar_update.end();++mi)
      mi->second = EMPTY ;

    bool finished = false;

    while(!finished) {
      for(fi=fset.begin();fi!=fset.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        // partial vmap_souce_exist
        entitySet srcs = fctrl.nr_sources ;
        for(int i=0;i<fctrl.recursion_maps.size();++i) {
          entitySet sdelta = tvar_computed[fctrl.recursion_maps[i].v] ;
          for(int j=fctrl.recursion_maps[i].mapvec.size()-1;j>=0;--j)
            sdelta = fctrl.recursion_maps[i].mapvec[j]->preimage(sdelta).first ;
          fctrl.recursion_maps[i].total += sdelta ;
          srcs &= fctrl.recursion_maps[i].total ;
        }
        if(fctrl.use_constraints) {
          if((srcs & fctrl.constraints) != fctrl.constraints) {
            cerr << "recursive rule: " << *fi
                 << " cannot supply all entitites in constraint" << endl ;
            cerr << "constraints = " << fctrl.constraints << endl ;
            entitySet sac = srcs & fctrl.constraints ;
            cerr << "srcs & constraints = " << sac << endl ;
            //exit(-1) ;
          }
          srcs &= fctrl.constraints ;
        }

        fctrl.control_list.push_back(srcs) ;

        for(int i=0;i<fctrl.target_maps.size();++i) {
          entitySet trgts = srcs ;
          for(int j=0;j<fctrl.target_maps[i].mapvec.size();++j)
            trgts = fctrl.target_maps[i].mapvec[j]->image(trgts) ;
          tvar_update[fctrl.target_maps[i].v] += trgts ;
          fctrl.generated[fctrl.target_maps[i].v] += trgts ;
        }
#ifdef VERBOSE
        cout << "recursive rule " << *fi << " generating " << srcs << endl ;
#endif
      }
      finished = true ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
        entitySet tvar_delta = tvar_update[*vi] - tvar_computed[*vi] ;
        if(tvar_delta!=EMPTY)
          finished = false ;
        tvar_computed[*vi] = tvar_update[*vi] ;
      }
    }

    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
          mi!=fctrl.generated.end();++mi)
        facts.set_existential_info(mi->first,*fi,mi->second) ;
    }
  }

  void recurse_compiler::process_var_requests(fact_db &facts) {
    ruleSet::const_iterator fi ;
    for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      entitySet control = process_rule_requests(*fi,facts) ;
      list<entitySet>::iterator ci ;
      entitySet total ;
      for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
        *ci -= total ;
        total += *ci ;
      }
      do {
        if(fctrl.control_list.back() == EMPTY)
          fctrl.control_list.pop_back() ;
        if(!fctrl.control_list.empty())
          fctrl.control_list.back() &= control ;
      } while(!fctrl.control_list.empty() && fctrl.control_list.back() == EMPTY) ;
    }
  }

  executeP recurse_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> el = new execute_list ;
    
    map<rule, list<entitySet>::const_iterator> rpos ;
    ruleSet::const_iterator ri ;
    
    for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
      rpos[*ri] = control_set[*ri].control_list.begin() ;

    bool finished ;
    do {
      finished = true ;
      
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
        const fcontrol &fctrl = control_set[*ri] ;
        list<entitySet>::const_iterator &li = rpos[*ri] ;
        if(li != fctrl.control_list.end()) {
          finished = false ;
          if(li->size() != 0) {
            const entitySet &exec_seq = *li ;

            if(num_threads > 1 && exec_seq.size() > 1 &&
               (*ri).get_info().rule_impl->thread_rule()) {
              execute_par *ep = new execute_par ;
              parallel_schedule(ep,*li,*ri,facts) ;
              el->append_list(ep) ;
            } else {
              el->append_list(new execute_rule(*ri,sequence(*li),facts)) ;
            }
          }
          li++ ;
        }
      }
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
    } while(!finished) ;
    
    if(el->size() == 0)
      return 0 ;
    else
      return executeP(el) ;
  }

  
  void dag_compiler::set_var_existence(fact_db &facts) {
    for(int i=0;i<rule_schedule.size();++i) 
      calc(rule_schedule[i])->set_var_existence(facts) ;
  }

  void dag_compiler::process_var_requests(fact_db &facts) {
    vector<rule>::reverse_iterator ri ;
    for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
      calc(*ri)->process_var_requests(facts) ;
  }

  executeP dag_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> elp = new execute_list ;

    vector<digraph::vertexSet>::const_iterator i ;
    for(i=dag_sched.begin();i!=dag_sched.end();++i) {
      ruleSet rules = extract_rules(*i) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        elp->append_list(calc(*ri)->create_execution_schedule(facts)) ;
      if(rules.size() > 0 && num_threads > 1)
        elp->append_list(new execute_thread_sync) ;
    }
    return executeP(elp) ;
  }

  conditional_compiler::conditional_compiler(decompose_graph *gr,
                                                 digraph gin,
                                                 variable conditional) {
    graph_ref = gr ;
    dag = gin ;
    cond_var = conditional ;
    dag_sched = schedule_dag(dag) ;
    extract_rule_sequence(rule_schedule,dag_sched) ;
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(int i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
#endif
  }

    
  void conditional_compiler::set_var_existence(fact_db &facts) {
    for(int i=0;i<rule_schedule.size();++i) 
      calc(rule_schedule[i])->set_var_existence(facts) ;
  }

  void conditional_compiler::process_var_requests(fact_db &facts) {
    vector<rule>::reverse_iterator ri ;
    for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
      calc(*ri)->process_var_requests(facts) ;

  }

  executeP conditional_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> elp = new execute_list ;

    vector<digraph::vertexSet>::const_iterator i ;
    for(i=dag_sched.begin();i!=dag_sched.end();++i) {
      ruleSet rules = extract_rules(*i) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        elp->append_list(calc(*ri)->create_execution_schedule(facts)) ;
      if(rules.size() > 0 && num_threads > 1)
        elp->append_list(new execute_thread_sync) ;
    }
    return new execute_conditional(executeP(elp),cond_var) ;
  }


  void promote_compiler::set_var_existence(fact_db &facts) {
  }

  void promote_compiler::process_var_requests(fact_db &facts) {
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }

  void generalize_compiler::set_var_existence(fact_db &facts) {
  }

  void generalize_compiler::process_var_requests(fact_db &facts) {
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }

  void priority_compiler::set_var_existence(fact_db &facts) {
  }

  void priority_compiler::process_var_requests(fact_db &facts) {
  }

  executeP priority_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }


  double get_timer() {
#ifdef PROFILE_CODE
    clock_t tc ;
    static double to = 0;
    double tn,t ;
    tc = clock() ;
    tn = tc/1000000.0 ;
    t = tn - to ;
    to = tn ;
    return t ;
#else
    return -1.0 ;
#endif
  }


  int num_threads = 1 ;
  
  executeP create_execution_schedule(rule_db &rdb,
                                     fact_db &facts,
                                     std::string target_string,
                                     int nth) {
    num_threads = min(nth,max_threads) ;

    double timer = get_timer() ;

    variableSet given = facts.get_typed_variables() ;
    variableSet target(expression::create(target_string)) ;

    cout << "generating dependency graph..." << endl ;
    digraph gr = dependency_graph(rdb,given,target).get_graph() ;

    // If graph is empty, return a null schedule 
    if(gr.get_target_vertices() == EMPTY)
      return executeP(0) ;

    cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr) ;
    
    cout << "decomposing graph..." << endl ;
    decompose_graph decomp(gr,given,target) ;

    timer = get_timer() ;
    cout << "Graph Processing Time: "<<timer << " seconds" << endl ;
    
    cout << "existential analysis..." << endl ;
    decomp.existential_analysis(facts) ;
    
    cout << "creating execution schedule..." << endl;
    executeP sched =  decomp.execution_schedule(facts,num_threads) ;
    
    timer = get_timer() ;
    cout << "Schedule Generation Time: " << timer << " seconds" << endl ;
    return sched ;
  }
}

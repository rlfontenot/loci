#include "visitor.h"
#include "visit_tools.h"
#include <sys/stat.h>
#include "comp_tools.h"

#include <Tools/stream.h>
#include "dist_tools.h"

using std::vector ;
using std::map ;
using std::cerr;
using std::endl;

namespace Loci {

  ///////////////////////////////////////////////////////////////////
  // orderVisitor
  ///////////////////////////////////////////////////////////////////
  // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet>
  orderVisitor::order_dag(const digraph &g,
                          digraph::vertexSet start_vertices,
                          digraph::vertexSet only_vertices) {
    digraph gt = g.transpose() ;
    
    vector<digraph::vertexSet> schedule ; 
    // First schedule any vertices that have no edges leading into them
    //and have not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    working &= only_vertices ;
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

  void orderVisitor::visit(loop_compiler& lc) {
    
    lc.collapse_sched = order_dag(lc.collapse_gr) ;
    lc.advance_sched = order_dag(lc.advance_gr) ;
  }

  void orderVisitor::visit(dag_compiler& dc) {
    dc.dag_sched = order_dag(dc.dag_gr) ;
  }

  void orderVisitor::visit(conditional_compiler& cc) {
    cc.dag_sched = order_dag(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // assembleVisitor
  /////////////////////////////////////////////////////////////////
  namespace {
    class check_dump_on_startup {
      bool checked_dump ;
      bool do_dump ;
    public :
      check_dump_on_startup() {
        do_dump = false ;
        checked_dump = false ;
      }
      bool ok() {
        if(!checked_dump) {
          do_dump=true ;
          struct stat statbuf ;
          if(stat("dump_vars",&statbuf))
            do_dump = false ;
          else if(!S_ISDIR(statbuf.st_mode)) 
            do_dump = false ;
          checked_dump = true ;
        } 
        return do_dump ;
      }
      
    } ;
    
    check_dump_on_startup check_dump_vars ;
    
    class execute_dump_var : public execute_modules {
      variableSet dump_vars ;
    public:
      execute_dump_var(variableSet vars) : dump_vars(vars) {}
      virtual void execute(fact_db &facts) ;
      virtual void Print(std::ostream &s) const ;
    } ;
    
    map<variable, int> dump_var_lookup ;
    
    void execute_dump_var::execute(fact_db &facts) {
      
      for(variableSet::const_iterator vi=dump_vars.begin();
          vi!=dump_vars.end();++vi) {
        variable v = *vi ;
        debugout << "dumping variable " << v << endl ;
        storeRepP st = facts.get_variable(v) ;
        storeRepP sc = st ;
        if(facts.isDistributed() && st->RepType() == STORE) {
          sc = collect_store(st,facts) ;
        }
        if(MPI_rank == 0) {
          ostringstream oss ;
          oss << "dump_vars/"<<v ;
          if(dump_var_lookup.find(v) ==dump_var_lookup.end())
            dump_var_lookup[v] = 0 ;
          if(dump_var_lookup[v] != 0)
            oss << "."<<dump_var_lookup[v] ;
          dump_var_lookup[v]++ ;
          
          string filename = oss.str() ;
          ofstream ofile(filename.c_str(),ios::out) ;
          ofile.precision(6) ;
          sc->Print(ofile) ;
        }
      }
    }
    
    void execute_dump_var::Print(std::ostream &s) const {
      s << "dumping variables " << dump_vars << endl ;
    }
    
    
    class dump_vars_compiler : public rule_compiler {
      variableSet dump_vars ;
    public:
      dump_vars_compiler(variableSet &vars) : dump_vars(vars) {}
      virtual void compile() {}
      virtual void accept(visitor& v) {}
      virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
      virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
      virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) {
        return executeP(new execute_dump_var(dump_vars)) ;
      }
    } ;
    
  } // end of namespace
  
  void assembleVisitor::compile_dag_sched
  (std::vector<rule_compilerP> &dag_comp,
   const std::vector<digraph::vertexSet> &dag_sched,
   const rulecomp_map &rcm,
   const digraph &dag) {
    digraph dagt = dag.transpose() ;
    if(dag_sched.size() == 0)
      return ;
   

    for(unsigned int i=0;i<dag_sched.size();++i) {
      //Loci::debugout << " in comp_dag.cc dag_sched[i] = " << dag_sched[i] << endl ;
      variableSet vars = extract_vars(dag_sched[i]) ;
      //Loci::debugout << " in comp_dag.cc vars = " << vars  << endl ;
      ruleSet rules = extract_rules(dag_sched[i]) ;
      if(rules == EMPTY && i+1<dag_sched.size()) {
        ++i ;
        vars += extract_vars(dag_sched[i]) ;
        rules = extract_rules(dag_sched[i]) ;
      }

      variableSet barrier_vars, reduce_vars,singleton_vars,all_vars ;
      variableSet::const_iterator vi ;
      
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(dagt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        bool recursive = false ;
        bool unit_rule_exists = false ;
        for(ri=var_rules.begin();ri!=var_rules.end();++ri) {
          //if(ri->get_info().rule_class != rule::INTERNAL) {
          if(!is_virtual_rule(*ri)) {
            use_rules += *ri ;

            if(ri->get_info().rule_class == rule::INTERNAL)
              if(ri->get_info().qualifier() == "CHOMP") {
                // we need to actually look into the chomping
                // graph to find out the rules that generate
                // this variable and dertermin the types of
                // rules
                rulecomp_map::const_iterator rmi ;
                rmi = rcm.find(*ri) ;
                FATAL(rmi == rcm.end()) ;
                rule_compilerP rcp = rmi->second ;
                chomp_compiler* chc =
                  dynamic_cast<chomp_compiler*>(&(*(rcp))) ;
                FATAL(chc == 0) ;
                digraph cgrt = (chc->chomp_graph).transpose() ;
                ruleSet chomp_var_rules = extract_rules(cgrt[vi->ident()]) ;
                for(ruleSet::const_iterator cri=chomp_var_rules.begin() ;
                    cri!=chomp_var_rules.end();++cri) {
                  
                  rule_implP crimp = cri->get_rule_implP() ;
                  if(crimp->get_rule_class() == rule_impl::POINTWISE)
                    pointwise = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::UNIT ||
                     crimp->get_rule_class() == rule_impl::APPLY)
                    reduction = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::UNIT)
                    unit_rule_exists = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::SINGLETON)
                    singleton = true ;
                }
                continue ;
              }
            
            rule_implP rimp = ri->get_rule_implP() ;
            if(rimp->get_rule_class() == rule_impl::POINTWISE)
              pointwise = true ;
            
            if(rimp->get_rule_class() == rule_impl::UNIT ||
               rimp->get_rule_class() == rule_impl::APPLY)
              reduction = true ;

            if(rimp->get_rule_class() == rule_impl::UNIT)
              unit_rule_exists = true ;
            
            if(rimp->get_rule_class() == rule_impl::SINGLETON)
              singleton = true ;
          } else {
            if((ri->sources() & ri->targets()) != EMPTY)
              recursive = true ;
          }
        }
        
        WARN(reduction && pointwise || pointwise && singleton ||
             reduction && singleton) ;

        if((use_rules != EMPTY)) {
          if(pointwise && !recursive && (vi->get_info().name != "OUTPUT")){
            barrier_vars += *vi ;
          }
          if(pointwise && recursive && (vi->get_info().name != "OUTPUT")) {
            all_vars += *vi ;
          }

          /*
          if(all_reduce_vars.inSet(*vi))
          reduce_vars += *vi ;
          */
          if(reduction && unit_rule_exists)
            reduce_vars += *vi ;
          
          if(singleton) {
            singleton_vars += *vi ;
          }
        }

      }

      all_vars += barrier_vars ;
      dag_comp.push_back(new barrier_compiler(barrier_vars)) ;
      
      all_vars += singleton_vars ;

      if(singleton_vars != EMPTY)
        dag_comp.push_back(new singleton_var_compiler(singleton_vars)) ;

      all_vars += reduce_vars;

      vector<CPTR<joiner> > join_op_vector ;
      vector<rule> unit_rule_vector ;
      vector<variable> reduce_var_vector ;
      
      for(variableSet::const_iterator rdvi=reduce_vars.begin();
          rdvi!=reduce_vars.end();++rdvi) {
        map<variable,pair<rule,CPTR<joiner> > >::const_iterator xi ;
        xi = reduceInfo.find(*rdvi) ;
        FATAL(xi == reduceInfo.end()) ;
        rule unit_rule = (xi->second).first ;
        CPTR<joiner> join_op = (xi->second).second ;

        if(join_op != 0) {
          storeRepP sp = join_op->getTargetRep() ;
          /* old code
             if(sp->RepType() == PARAMETER) {
             dag_comp.push_back(new reduce_param_compiler(xi->first,unit_rule,
             join_op)) ;
             }
          */
          if(sp->RepType()== PARAMETER) {
            reduce_var_vector.push_back(xi->first) ;
            unit_rule_vector.push_back(unit_rule) ;
            join_op_vector.push_back(join_op) ;
          }
          else if (sp->RepType() == BLACKBOX) {
            cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
          }else {
            dag_comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
        }
      }
      if(reduce_var_vector.size() != 0) 
        dag_comp.push_back(new reduce_param_compiler(reduce_var_vector, unit_rule_vector, join_op_vector));
      
      
      if(check_dump_vars.ok()) {
        if(all_vars != EMPTY) 
          dag_comp.push_back(new dump_vars_compiler(all_vars)) ;
      }
      
      if(rules != EMPTY) {
        ruleSet::const_iterator ri ;
        for(ri=rules.begin();ri!=rules.end();++ri) {
          rulecomp_map::const_iterator rmi ;
          rmi = rcm.find(*ri) ;
          FATAL(rmi == rcm.end()) ;
          dag_comp.push_back(rmi->second) ;
        }
      }
    }
  }

  void assembleVisitor::visit(loop_compiler& lc) {
    compile_dag_sched(lc.collapse_comp,lc.collapse_sched,
                      lc.rule_compiler_map,lc.loop_gr) ;
    compile_dag_sched(lc.advance_comp,lc.advance_sched,
                      lc.rule_compiler_map,lc.loop_gr) ;
  }
  
  void assembleVisitor::visit(dag_compiler& dc) {
    compile_dag_sched(dc.dag_comp,dc.dag_sched,
                      dc.rule_compiler_map,dc.dag_gr) ;
  }
  
  void assembleVisitor::visit(conditional_compiler& cc) {
    compile_dag_sched(cc.dag_comp,cc.dag_sched,
                      cc.rule_compiler_map,cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // simLazyAllocSchedVisitor
  /////////////////////////////////////////////////////////////////
  vector<digraph::vertexSet>
  simLazyAllocSchedVisitor::schedule(const digraph &gr) {
    // first we get all vertices except all the allocate rules
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;

    digraph::vertexSet alloc_rules_vertices ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(ri->get_info().qualifier() == "ALLOCATE")
        alloc_rules_vertices += ri->ident() ;

    vector<digraph::vertexSet> first_sched ;
    
    first_sched =
      orderVisitor::order_dag(gr,alloc_rules_vertices,
                              gr.get_all_vertices()-alloc_rules_vertices) ;
    
    // then we schedule those alloc rules
    digraph grt = gr.transpose() ;
    vector<digraph::vertexSet> final_sched ;

    digraph::vertexSet last ;
    digraph::vertexSet sched_alloc ;
    for(vector<digraph::vertexSet>::const_iterator vi=first_sched.begin();
        vi!=first_sched.end();++vi) {
      digraph::vertexSet step = *vi ;
      ruleSet step_rules = extract_rules(step) ;
      
      digraph::vertexSet pre_vertices ;
      for(ruleSet::const_iterator ri=step_rules.begin();
          ri!=step_rules.end();++ri)
        pre_vertices += grt[ri->ident()] ;

      ruleSet pre_rules = extract_rules(pre_vertices) ;
      digraph::vertexSet needed_alloc ;
      for(ruleSet::const_iterator ri=pre_rules.begin();
          ri!=pre_rules.end();++ri)
        if(ri->get_info().qualifier() == "ALLOCATE")
          needed_alloc += ri->ident() ;

      needed_alloc -= sched_alloc ;
      sched_alloc += needed_alloc ;

      digraph::vertexSet sum = needed_alloc + last ;
      if(sum != EMPTY)
        final_sched.push_back(sum) ;

      last = step ;
    }
    if(last != EMPTY)
      final_sched.push_back(last) ;

    return final_sched ;
  }
  
  void simLazyAllocSchedVisitor::visit(loop_compiler& lc) {    
    lc.collapse_sched = schedule(lc.collapse_gr) ;
    lc.advance_sched = schedule(lc.advance_gr) ;
  }

  void simLazyAllocSchedVisitor::visit(dag_compiler& dc) {
    dc.dag_sched = schedule(dc.dag_gr) ;
  }

  void simLazyAllocSchedVisitor::visit(conditional_compiler& cc) {
    cc.dag_sched = schedule(cc.cond_gr) ;
  }



  
} // end of namespace Loci

#include "visitor.h"
#include <sys/stat.h>
#include "comp_tools.h"

#include <Tools/stream.h>
#include <distribute.h>

using std::vector ;
using std::map ;


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
      
      std::map<variable,ruleSet> reduce_info ;
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(dagt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        bool recursive = false ;
        bool unit_rule_exists = false ;
        for(ri=var_rules.begin();ri!=var_rules.end();++ri)
          if(ri->get_info().rule_class != rule::INTERNAL) {
            use_rules += *ri ;
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
        
        WARN(reduction && pointwise || pointwise && singleton ||
             reduction && singleton) ;

        if((use_rules != EMPTY)) {
          if(pointwise && !recursive && (vi->get_info().name != "OUTPUT")) {
            barrier_vars += *vi ;
          }
          if(pointwise && recursive && (vi->get_info().name != "OUTPUT")) {
            all_vars += *vi ;
          }
          if(reduction) {
            reduce_info[*vi] = use_rules ;
            if(reduction && unit_rule_exists)
              reduce_vars += *vi ;
          } if(singleton) {
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
      
      if(reduce_vars != EMPTY) {
        std::map<variable,ruleSet>::const_iterator xi ;
        variableSet vars ;
        for(xi=reduce_info.begin();xi!=reduce_info.end();++xi) {
          vars += xi->first ;
          rule unit_rule ;
          CPTR<joiner> join_op = CPTR<joiner>(0) ;
          ruleSet::const_iterator ri ;
          for(ri=xi->second.begin();ri!=xi->second.end();++ri) {
            if(ri->get_rule_implP()->get_rule_class() == rule_impl::UNIT)
              unit_rule = *ri ;
            else if(ri->get_rule_implP()->get_rule_class() == rule_impl::APPLY){
              //              if(join_op == 0)
                join_op = ri->get_rule_implP()->get_joiner() ;
              //#define TYPEINFO_CHECK
#ifdef TYPEINFO_CHECK
              else
                if(typeid(*join_op) !=
                   typeid(*(ri->get_rule_implP()->get_joiner()))) {
                  cerr << "Warning:  Not all apply rules for variable " << xi->first << " have identical join operations!" << endl ;
                }
#endif
            } else {
              cerr << "Warning: reduction variable " << xi->first
                   << " has a non-reduction rule contributing to its computation,"
                   << endl << "offending rule is " << *ri << endl ;
            }
          }
          if(join_op == 0) {
            cerr << "unable to find any apply rules to complete the reduction defined by rule"
                 << endl
                 << unit_rule << endl ;
          }
          FATAL(join_op == 0) ;
          storeRepP sp = join_op->getTargetRep() ;
          if(sp->RepType() == PARAMETER) {
            dag_comp.push_back(new reduce_param_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
          else if (sp->RepType() == BLACKBOX) {
	    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
          }else {
            dag_comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
        }
      }

      
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

} // end of namespace Loci

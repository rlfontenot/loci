#include "comp_tools.h"
#include <vector>
using std::vector ;

namespace Loci {

  void barrier_compiler::set_var_existence(fact_db &facts) {
  }

  void barrier_compiler::process_var_requests(fact_db &facts) {
  }

  executeP barrier_compiler::create_execution_schedule(fact_db &facts) {
    //    if(num_threads > 1)
    ostringstream oss ;
    oss << barrier_vars ;
    return new execute_thread_sync(oss.str()) ;
  }

  
  void compile_dag_sched(std::vector<rule_compilerP> &dag_comp,
                         const std::vector<digraph::vertexSet> &dag_sched,
                         const rulecomp_map &rcm) {
    for(int i=0;i<dag_sched.size();++i) {
      variableSet vars = extract_vars(dag_sched[i]) ;
      ruleSet rules = extract_rules(dag_sched[i]) ;
      if(vars != EMPTY) {
        dag_comp.push_back(new barrier_compiler(vars)) ;
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

  
  dag_compiler::dag_compiler(rulecomp_map &rule_process, digraph dag) {

    std::vector<digraph::vertexSet> dag_sched = schedule_dag(dag) ;
    compile_dag_sched(dag_comp,dag_sched,rule_process) ;
    
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

  void dag_compiler::set_var_existence(fact_db &facts) {
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts) ;
  }

  void dag_compiler::process_var_requests(fact_db &facts) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts) ;
  }

  executeP dag_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> elp = new execute_list ;

    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts)) ;
    }

    return executeP(elp) ;
  }
}


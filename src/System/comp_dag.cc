#include "comp_tools.h"

namespace Loci {
  dag_compiler::dag_compiler(rulecomp_map &rp, digraph gin) : rule_process(rp) {
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

  void dag_compiler::set_var_existence(fact_db &facts) {
    for(int i=0;i<rule_schedule.size();++i) {
      if(rule_process.find(rule_schedule[i]) == rule_process.end())
        cerr << "didn't find compiler for " << rule_schedule[i]
             << endl ;
      calc(rule_schedule[i])->set_var_existence(facts) ;
    }
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
}


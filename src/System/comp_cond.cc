#include "comp_tools.h"

namespace Loci {
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

  conditional_compiler::conditional_compiler(rulecomp_map &rp,
                                                 digraph gin,
                                             variable conditional)
  : rule_process(rp) {
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




}

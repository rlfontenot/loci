#include "comp_tools.h"
#include <vector>
using std::vector ;

using std::ostream ;
using std::endl ;

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
  
  conditional_compiler::conditional_compiler(rulecomp_map &rule_process,
					     digraph dag,
                                             variable conditional)
  {
    cond_var = conditional ;
    std::vector<digraph::vertexSet> dag_sched = schedule_dag(dag) ;
    compile_dag_sched(dag_comp,dag_sched,rule_process,dag) ;
    
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(size_t i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
#endif
  }
  
  
  void conditional_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void conditional_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP conditional_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }
    
    return new execute_conditional(executeP(elp),cond_var) ;
  }




}

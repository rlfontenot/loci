#include "comp_tools.h"
#include <vector>
using std::vector ;

using std::ostream ;
using std::endl ;

#include "visitorabs.h"

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
                                             variable conditional,
                                             int id):cid(id)
  {
    cond_var = conditional ;
    
    cond_gr = dag ;
    ruleSet allrules = extract_rules(cond_gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }    
  }
  
  /* compile() method has been removed from the rule_compiler class
     hierarchy, this one exists here just for the DEBUG code reference.
     we may later move it into the visitor class.
  void conditional_compiler::compile()
  {
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(size_t i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != cond_gr.get_all_vertices()) ;
#endif
  }
  */
  
  void conditional_compiler::accept(visitor& v) {
    v.visit(*this) ;
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

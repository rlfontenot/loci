#include "comp_tools.h"
#include <vector>

#include <Tools/stream.h>

using std::vector ;

#include <sys/stat.h>

#include <distribute.h>
#include <parameter.h>

#include <map>
using std::map ;

#include "visitor.h"

namespace Loci {
  
  dag_compiler::dag_compiler(rulecomp_map &rule_process, digraph dag, int id):cid(id) {

    ////////////////////
    // store the dag structure and the relevant rulecompiler map
    dag_gr = dag ;
    ruleSet allrules = extract_rules(dag_gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }
    ////////////////////
  }

  /* compile() method has been removed from the rule_compiler class
     hierarchy, this one exists here just for the DEBUG code reference.
     we may later move it into the visitor class.
     
  void dag_compiler::compile() {
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(size_t i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag_gr.get_all_vertices()) ;
    if(allvertices != dag_gr.get_all_vertices()) {
      digraph::vertexSet leftout = allvertices ^ dag_gr.get_all_vertices() ;
      cerr << "leftout rules= " << extract_rules(leftout) << endl ;
      cerr << "leftout vars = " << extract_vars(leftout) << endl ;
    }
#endif    
  }
  */
  
  void dag_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void dag_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void dag_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP dag_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    return executeP(elp) ;
  }
}


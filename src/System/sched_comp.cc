#include "sched_tools.h"
#include "comp_tools.h"
#include "distribute.h"

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <set>
using std::set ;

//#define HACK ;

namespace Loci {
  class error_compiler : public rule_compiler {
  public:
    error_compiler() {}
    virtual void set_var_existence(fact_db &facts)
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual void process_var_requests(fact_db &facts) 
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual executeP create_execution_schedule(fact_db &facts)
    { cerr << "Internal consistency error" << endl ; exit(-1);
    return executeP(0);}
  } ;
 
  graph_compiler::graph_compiler(decomposed_graph &deco,variableSet
				 initial_vars) {
    fact_db_comm = new barrier_compiler(initial_vars) ;
    multiLevelGraph &mlg = deco.mlg ;
    vector<int> levels ;
    baserule = rule(mlg.toplevel) ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }
    for(vector<int>::reverse_iterator ii=levels.rbegin();
        ii!=levels.rend();
        ++ii) {
      multiLevelGraph::subGraph *p = mlg.find(*ii) ;
      ruleSet level_rules = extract_rules(p->graph_v-mlg.subgraphs) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;

      digraph grt = gr.transpose() ;


      // Collect information on unit rules on this level
      map<rule,rule> apply_to_unit ;
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
          variable unit_var = *(ri->targets().begin()) ;
          ruleSet apply_rules = extract_rules(grt[unit_var.ident()]) ;
          for(ruleSet::const_iterator rii=apply_rules.begin();
              rii!=apply_rules.end();
              ++rii) {
            apply_to_unit[*rii] = *ri ;
          }
        }
      }

      // Assign compilers to rules on this level
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->type() == rule::INTERNAL) {
          if(ri->get_info().qualifier() == "promote")
            rule_process[*ri] = new promote_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "generalize")
            rule_process[*ri] = new generalize_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "priority")
            rule_process[*ri] = new priority_compiler(*ri) ;
          else
            rule_process[*ri] = new error_compiler ;
        } else {
          if(ri->get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
            rule_process[*ri] = new impl_compiler(*ri) ;
          else
            rule_process[*ri] = new apply_compiler(*ri,apply_to_unit[*ri]) ;
        }
      }

      // We have created compilers for every rule on this level, now
      // we create a compiler for this supernode

      rule snrule = rule(*ii) ;

      if(deco.loops.inSet(*ii)) {
        // Looping supernode
        rule_process[snrule] = new loop_compiler(rule_process,gr) ;
      } else if(deco.recursive.inSet(*ii)) {
        // recursive supernode
        ruleSet recurse_rules = extract_rules(p->graph_v) ;
        if(recurse_rules.size() == 1 &&
           recurse_rules.begin()->get_info().desc.constraints.size() == 0
           && (MPI_processes ==1 || (*recurse_rules.begin()).get_rule_implP()->is_relaxed()))
          // Single rule recursion
          rule_process[snrule] =
            new impl_recurse_compiler(*(recurse_rules.begin())) ;
        else
          // Multi-rule recursion
          rule_process[snrule] =
            new recurse_compiler(rule_process, recurse_rules) ;
      } else if(deco.conditional.inSet(*ii)) {
        // conditional supernode
        variable cond_var = *(snrule.get_info().desc.conditionals.begin()) ;
        rule_process[snrule] =
          new conditional_compiler(rule_process,gr,cond_var) ;
      } else {
        // DAG supernode
        rule_process[snrule] = new dag_compiler(rule_process,gr) ;
      }
    }
  }

  void graph_compiler::existential_analysis(fact_db &facts) {
    fact_db_comm->set_var_existence(facts) ;
    (rule_process[baserule])->set_var_existence(facts) ;
    variableSet var_requests = baserule.targets() ;
    variableSet::const_iterator vi ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = facts.variable_existence(*vi) ;
      facts.variable_request(*vi,vexist) ;
    }
    (rule_process[baserule])->process_var_requests(facts) ;
    fact_db_comm->process_var_requests(facts) ;
  }
  
  class allocate_all_vars : public execute_modules {
  public:
    allocate_all_vars() { control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(ostream &s) const ;
  } ;


  fact_db *exec_current_fact_db = 0 ;
  
  void allocate_all_vars::execute(fact_db &facts) {
    exec_current_fact_db = &facts ;
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
        for(vii=aliases.begin();vii!=aliases.end();++vii) {
#ifdef HACK
	  all_requests += facts.variable_existence(*vii) ;
#else
	  all_requests += facts.get_variable_requests(*vii) ;
#endif	  
	}
#ifdef DEBUG
        //debugout << "allocating " << *vi << " for entities " << all_requests
	//       << endl ;
#endif
        srp->allocate(all_requests) ;
      }
    }
  }

  void allocate_all_vars::Print(ostream &s) const {
    s << "allocate all variables" << endl ;
  }


  executeP graph_compiler::execution_schedule(fact_db &facts, int nth) {

    CPTR<execute_list> schedule = new execute_list ;
    schedule->append_list(fact_db_comm->create_execution_schedule(facts));
    schedule->append_list(new allocate_all_vars) ;
    schedule->append_list(new execute_create_threads(nth)) ;
    executeP top_level_schedule = (rule_process[baserule])->
      create_execution_schedule(facts) ;
    if(top_level_schedule == 0) 
      return executeP(0) ;
    schedule->append_list(top_level_schedule) ;
    schedule->append_list(new execute_destroy_threads) ;
    return executeP(schedule) ;
  }
}

#include "sched_tools.h"
#include "comp_tools.h"
#include "distribute.h"

#include <Tools/stream.h>

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
  extern double total_memory_usage ;

  class error_compiler : public rule_compiler {
  public:
    error_compiler() {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds)
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) 
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds)
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
  
  void graph_compiler::existential_analysis(fact_db &facts, sched_db &scheds) {
    fact_db_comm->set_var_existence(facts, scheds) ;
    (rule_process[baserule])->set_var_existence(facts, scheds) ;
    variableSet var_requests = baserule.targets() ;
    variableSet::const_iterator vi ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi) ;
      scheds.variable_request(*vi,vexist) ;
    }
    (rule_process[baserule])->process_var_requests(facts, scheds) ;
    fact_db_comm->process_var_requests(facts, scheds) ;
  }
  
  class allocate_all_vars : public execute_modules {
    std::map<variable,entitySet> v_requests, v_existence ;
  public:
    allocate_all_vars() { control_thread = true ; }
    allocate_all_vars(fact_db &facts, sched_db &scheds) ;
    void fill_in_requests(fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  
  fact_db *exec_current_fact_db = 0 ;
  
  allocate_all_vars::allocate_all_vars(fact_db &facts, sched_db &scheds) {
    control_thread = true ;
    fill_in_requests(facts, scheds) ;
  }
  void allocate_all_vars::fill_in_requests(fact_db &facts, sched_db &scheds) {
    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi,vii ;
    
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      variableSet aliases = scheds.get_aliases(*vi) ;
      entitySet requests, existence ;
      for(vii=aliases.begin();vii!=aliases.end();++vii) {
	existence += scheds.variable_existence(*vii) ;
	requests += scheds.get_variable_requests(*vii) ;
      }
      v_requests[*vi] = requests ;
      v_existence[*vi] = existence ;
    }
  }
  
  void allocate_all_vars::execute(fact_db &facts) {
    exec_current_fact_db = &facts ;
    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi ;  
    double total_size = 0 ;
    entitySet dom, total, unused ;
    double total_wasted = 0 ;

    //#define DIAGNOSTICS
#ifdef DIAGNOSTICS
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->RepType() == Loci::STORE) {
	dom = v_requests[*vi] ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += srp->pack_size(dom) ; 
	total_wasted += srp->pack_size(unused) ;
      }
    }
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->RepType() == Loci::STORE) {
	dom = v_requests[*vi] ;
	double size, wasted_space ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	size = srp->pack_size(dom) ;
	wasted_space = srp->pack_size(unused) ;
	Loci::debugout << " ****************************************************" << endl ;
	Loci::debugout << " Total_size = " << total_size << endl ;
	Loci::debugout << "Variable = "  << *vi << endl ;
	Loci::debugout << "Domain = " << dom << endl ;
	Loci::debugout << "Size allocated = " << size << endl ; 
	if(facts.isDistributed() )  {
	  Loci::fact_db::distribute_infoP d ;
	  d   = facts.get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = srp->pack_size(clone) ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	  Loci::debugout << " My_entities = " << my_entities << endl ;
	  Loci::debugout << " Clone entities = " << clone << endl ;
	  Loci::debugout << "Memory required for the  clone region  = " << clone_size << endl ;
	  Loci::debugout << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  Loci::debugout << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	}
	Loci::debugout << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	Loci::debugout << "----------------------------------------------------" << endl;
	Loci::debugout << "Total wasted size = " << total_wasted << endl ;
	Loci::debugout << "Unused entities = " << unused << endl ;
	Loci::debugout << "Wasted space = " << wasted_space << endl ;
	Loci::debugout << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	Loci::debugout << " ***************************************************" << endl << endl << endl ;
      }
    }
#endif

    //#define HACK
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
#ifdef HACK
      entitySet alloc_dom = v_existence[*vi]  + srp->domain() ;
#else
      entitySet alloc_dom = v_requests[*vi] + srp->domain() ;
#endif
      if(srp->domain() == EMPTY) {
	srp->allocate(alloc_dom) ;
      }
      else {
	if(srp->RepType() == Loci::STORE) {
	  entitySet tmp = interval(alloc_dom.Min(), alloc_dom.Max()) ;
	  if(tmp.size() >= 2*srp->domain().size())
	    Loci::debugout << "Variable = " << *vi << "  more than twice the space allocated :  allocated over " << alloc_dom << " size = " << tmp.size()  << "  while domain is only  " << srp->domain() << " size = " << srp->domain().size() << endl ;
	  if(alloc_dom != srp->domain()) {
	    Loci::debugout << "reallocating " << *vi << "  over  " << alloc_dom << " initially it was over  " << srp->domain() << endl ;
	    srp->allocate(alloc_dom) ;
	  }
	} 
      }
    }
    total_memory_usage = total_size + total_wasted ;
  }
  
  void allocate_all_vars::Print(std::ostream &s) const {
    s << "allocate all variables" << endl ;
  }


  executeP graph_compiler::execution_schedule(fact_db &facts, sched_db &scheds, int nth) {

    CPTR<execute_list> schedule = new execute_list ;
    schedule->append_list(new allocate_all_vars(facts, scheds)) ;
    schedule->append_list(new execute_create_threads(nth)) ;
    schedule->append_list(fact_db_comm->create_execution_schedule(facts, scheds));
    executeP top_level_schedule = (rule_process[baserule])->
      create_execution_schedule(facts, scheds) ;
    if(top_level_schedule == 0) 
      return executeP(0) ;
    schedule->append_list(top_level_schedule) ;
    schedule->append_list(new execute_destroy_threads) ;
    return executeP(schedule) ;
  }
}

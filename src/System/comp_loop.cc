//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include "comp_tools.h"
#include <vector>
using std::vector ;
using std::list ;
using std::map ;

using std::ostream ;
using std::endl ;


#include "visitorabs.h"
#include "loci_globs.h"
namespace Loci {
  int printLevel = 0 ;
  
  class execute_loop : public execute_modules {
    executeP collapse, advance ;
    variable cvar ;
    variable tvar ;
    time_ident tlevel ;
    list<list<variable> > rotate_lists ;
  public:
    execute_loop(const variable &cv,
                 const executeP &col, const executeP &adv,
                 const time_ident &tl, 
                 list<list<variable> > &rl) :
      collapse(col),advance(adv),
      cvar(cv),
      tlevel(tl),rotate_lists(rl) {
      warn(col==0 || advance==0) ; tvar = variable(tlevel) ;
      control_thread = true ;}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
	virtual string getName() { return "execute_loop";};
	} ;
  
  void execute_loop::execute(fact_db &facts) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;
    // initialize conditional variables to true
    //    *test = true ;
    
    param<int> time_var ;
    time_var = facts.get_variable(tvar) ;
    // Start iteration by setting iteration variable to zero 


    *time_var = 0 ;

    for(;;) { // Begin looping
      //First evaluate the collapse condition
      collapse->execute(facts) ;

#ifdef VERBOSE
      cerr << cvar << " = " << *test << endl ;
      cerr << tvar << " = " << *time_var << endl ;
#endif

      if(*test) {// If collapse condition satisfied, were finished
        return ;
      }
      // Advance to the next iterate
      advance->execute(facts) ;

      // We are finished with this iteration, so rotate variables and
      // add one to the iterate
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
#ifdef VERBOSE
        cerr << "rotating [ " ;
        for(list<variable>::const_iterator rlii=rli->begin();rlii!=rli->end();
            ++rlii) {
          cerr << *rlii << ' ' ;
        }
        cerr << "]" << endl ;
#endif
        facts.rotate_vars(*rli) ;
      }
      *time_var += 1 ;
    }
  }
  
  void execute_loop::Print(ostream &s) const {
    printIndent(s) ;
    s << "Iteration Loop{"<< tlevel << "} {" << endl ;
    printLevel+=2 ;
    //    printIndent(s) ;
    //    s << "--compute collapse of {" <<tlevel << "} , conditional on " << cvar << endl ;
    //    printIndent(s) ;
    //    s << "--collapse iteration {" << tlevel << "}"<< endl ;
    collapse->Print(s) ;
    s<< endl ;
    printIndent(s) ;
    s << "-------------- Exit of Loop{" << tlevel << "}"<< endl ;
    printIndent(s) ;
    s << "if(" << cvar << ") break ;" << endl ;
    s << endl ;
    //    printIndent(s) ;
    //    s << "--advance iteration {" << tlevel << "}" << endl ;
    advance->Print(s) ;
    printLevel-=2;
    printIndent(s) ;
    s << "} // {" << tlevel  << "}" << endl ;
  }
  
  inline bool offset_sort(const variable &v1, const variable &v2)
  { return v1.get_info().offset > v2.get_info().offset ; }
  
  execute_modules_decorator_factory* loop_compiler::decoratorFactory = NULL;
  
  loop_compiler::loop_compiler(rulecomp_map &rule_process, digraph dag, int id):cid(id) {
    ////////////////////
    // store the graph structure and the relevant rulecompiler map
    loop_gr = dag ;
    ruleSet allrules = extract_rules(loop_gr.get_all_vertices()) ;
    
    for(ruleSet::const_iterator ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }
    ////////////////////
    ruleSet collapse_rules, all_rules, loopset ;

    all_rules = extract_rules(loop_gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        loopset += ri->ident() ;
      if(ri->target_time().before(ri->source_time()))
        collapse_rules += *ri ;
    }

    if(loopset.size() != 1) {
      cerr << "internal consistency error, loopset size != 1" << endl ;
      cerr << "looping rules = " << loopset << endl ;
      Loci::Abort() ;
    }

    tlevel = loopset.begin()->source_time() ;
    loop_gr.remove_vertex((*loopset.begin()).ident()) ;
    loop_gr.remove_dangling_vertices() ;

    
    variable OUTPUT(variable("OUTPUT"),tlevel) ;

    digraph loop_grt = loop_gr.transpose() ;
  
    if(collapse_rules.size() != 1 ||
       collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
      cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
           << endl << "error is not recoverable" << endl ;
      cerr << "rules = " << collapse_rules << endl ;
      Loci::Abort() ;
    }

    variableSet input_vars ;
    variableSet all_vars = extract_vars(loop_gr.get_all_vertices()) ;
    for(variableSet::const_iterator vi=all_vars.begin();
        vi!=all_vars.end();++vi ) {

        if(vi->get_info().offset == 1)
          advance_vars += *vi ;
        if(loop_grt[vi->ident()] == EMPTY)
          input_vars += *vi ;
    }
        
  
    collapse_vars = collapse_rules.begin()->targets() ;
    cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

    variableSet outcond ;
    ruleSet outrules = extract_rules(loop_grt[OUTPUT.ident()]) ;
    for(ruleSet::const_iterator ri=outrules.begin();ri!=outrules.end();++ri) 
      outcond += ri->get_info().desc.conditionals ;

    digraph::vertexSet collapse_search = collapse_vars ;

    collapse_search += input_vars ;
    collapse_search += outcond ;
    
    digraph::vertexSet collapse_part = visit_vertices(loop_grt,collapse_search) ;
    collapse_part += collapse_search ;
    // for the collapse part, we also need to add all the rule targets
    // into the graph. Because the collapse part is formed by searching
    // back from the collapse_search, it is possible that we won't include
    // a rule's targets if the rule's targets do not connect to other parts
    // of the collapse graph
    ruleSet collapse_rule_set = extract_rules(collapse_part) ;
    for(ruleSet::const_iterator ri=collapse_rule_set.begin();
        ri!=collapse_rule_set.end();++ri) {
      collapse_part += extract_vars(loop_gr[ri->ident()]) ;
    }

#ifdef MOVE_OUTPUT_TO_COLLAPSE
    digraph::vertexSet output_set ;

    for(ruleSet::const_iterator ri=outrules.begin();ri!=outrules.end();++ri)
      if((ri->sources() - collapse_part) == EMPTY) {
          collapse_search += ri->ident() ;
      } else {
          // Later make this verbose
          //          cerr << *ri << "put in advance because of " << extract_vars(ri->sources()-collapse_part) << endl ;
          digraph::vertexSet tmp_search  ;
          tmp_search += ri->ident() ;
          output_set += visit_vertices(loop_grt,tmp_search)+tmp_search ;
      }
    //    output_set -= visit_vertices(loop_grt,collapse_search)+collapse_search ;

    //    cerr << "collapse group = " << extract_rules(output_set) << endl ;
#endif
    
    collapse_gr = loop_gr.subgraph(collapse_part) ;
    collapse_gr.remove_dangling_vertices() ;
    
    
    digraph::vertexSet collapse_rulesV = extract_rules(collapse_gr.get_all_vertices()) ;
    digraph::vertexSet advance_subset = loop_gr.get_all_vertices() - collapse_rulesV ;
    advance_gr = loop_gr.subgraph(advance_subset) ;
    advance_gr.remove_dangling_vertices() ;

    all_loop_vars = all_vars ;

  }

  /* compile() method has been removed from the rule_compiler class
     hierarchy, this one exists here just for the DEBUG code reference.
     we may later move it into the visitor class.
  void loop_compiler::compile()  {
    collapse_sched = schedule_dag(collapse_gr) ;
    advance_sched = schedule_dag(advance_gr) ;

    compile_dag_sched(collapse_comp,collapse_sched,
                      rule_compiler_map,loop_gr) ;
    compile_dag_sched(advance_comp,advance_sched,
                      rule_compiler_map,loop_gr) ;

    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices = visited ;
    for(size_t i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
    if(allvertices != dag.get_all_vertices()) {

      cerr << "Loci INTERNAL Consistency error!" << endl ;
      cerr << " rules NOT scheduled = " << endl
           << extract_rules(dag.get_all_vertices() - allvertices) << endl ;
      cerr << " variables NOT scheduled = "
           << extract_vars(dag.get_all_vertices() - allvertices) << endl ;
      vector<digraph::vertexSet> components =
        component_sort(dag).get_components() ;
      for(size_t i=0;i<components.size();++i) {
        if(components[i].size() > 1) {
          cerr << "reason: graph not a dag, strongly connected component found"
               << endl ;
          ruleSet rs = extract_rules(components[i]) ;
          variableSet vs = extract_vars(components[i]) ;
          cerr << "component rules = " << endl ;
          cerr << rs << endl ;
          cerr << "component variables = " << vs << endl ;

          for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) {
            rule r = *ri ;
            cerr << r << ":" ;
            digraph::vertexSet edges = dag[r.ident()] & components[i] ;
            cerr <<extract_vars(edges) << ' ' << extract_rules(edges) << endl ;
          }
          for(variableSet::const_iterator vi=vs.begin();vi!=vs.end();++vi) {
            variable v = *vi ;
            cerr << v << ":" ;
            digraph::vertexSet edges = dag[v.ident()] & components[i] ;
            cerr << extract_vars(edges) << ' ' << extract_rules(edges) << endl ;
          }
        }
      }
      Loci::Abort() ;
    }
  } 
  */

  void loop_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void loop_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    if(duplicate_work) {
      for(variableSet::const_iterator vi = advance_vars.begin();
	  vi != advance_vars.end(); vi++) {
	scheds.add_policy(*vi, sched_db::NEVER);
	variable tmp_var = vi->new_offset(vi->get_info().offset - 1);
	scheds.add_policy(tmp_var, sched_db::NEVER);
      }
    
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
        variableSet rvar ;
	list<variable>::const_iterator li ;
	for(li=rli->begin();li!=rli->end();++li) {
          rvar += *li ;
	  scheds.add_policy(*li, sched_db::NEVER);
	}
        scheds.set_variable_rotations(rvar) ;
      }
    }

    list<list<variable> >::const_iterator rli ;
    for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
      variableSet rvar ;
      list<variable>::const_iterator li ;
      for(li=rli->begin();li!=rli->end();++li) {
        rvar += *li ;
      }
      scheds.set_variable_rotations(rvar) ;
    }

    std::vector<rule_compilerP>::iterator i ;
    for(i=collapse_comp.begin();i!=collapse_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
    for(i=advance_comp.begin();i!=advance_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
    
  }
  
  void loop_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    variableSet var_requests = advance_vars ;
    variableSet::const_iterator vi ;
    
    Loci::fact_db::distribute_infoP d ;
    if(facts.isDistributed()) {
      d = facts.get_distribute_info() ;
    }
    
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi) ;
      if(facts.isDistributed()) 
	vexist &= d->my_entities;
      scheds.variable_request(*vi,vexist) ;
    }
    
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=advance_comp.rbegin();ri!=advance_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;

    for(ri=collapse_comp.rbegin();ri!=collapse_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;

    // spread requests over all time iterations
    list<list<variable> >::const_iterator rli ;
    for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
      list<variable>::const_iterator li ;
      entitySet tot_request ;
      for(li=rli->begin();li!=rli->end();++li)
        tot_request += scheds.get_variable_requests(*li) ;
      for(li=rli->begin();li!=rli->end();++li)
        scheds.variable_request(*li,tot_request) ;
    }


    // this block of code is apparently not used
//     for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
//       entitySet tot_request = scheds.get_variable_requests(*vi);
//       variable tmp_var = vi->new_offset(vi->get_info().offset - 1);
//       tot_request += scheds.get_variable_requests(tmp_var);
//       //      scheds.variable_request(*vi, tot_request);
//     }

    // then we need to remove all the non-store variables from
    // var_requests to prevent any parallel communications
    variableSet non_stores ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      storeRepP sp = facts.get_variable(*vi) ;
      if(sp->RepType() != STORE)
        non_stores += *vi ;
    }
    var_requests -= non_stores ;
    
    if(facts.isDistributed()) {
      advance_variables_barrier = barrier_process_rule_requests(var_requests, facts, scheds);
      advance_variables_barrier = sort_comm(advance_variables_barrier, facts);
    }
  }
  
  executeP loop_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    
    CPTR<execute_list> col = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=collapse_comp.begin();i!=collapse_comp.end();++i) {
      col->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    CPTR<execute_list> adv = new execute_list ;
    
    for(i=advance_comp.begin();i!=advance_comp.end();++i) {
      adv->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    if(facts.isDistributed()) {
	  executeP exec_comm = new execute_comm(advance_variables_barrier, facts);
	  if(decoratorFactory != NULL)
            exec_comm = decoratorFactory->decorate(exec_comm);
      adv->append_list(exec_comm);
    }
	
	executeP execute = new execute_loop(cond_var,executeP(col),executeP(adv),tlevel,rotate_lists) ;
	if(decoratorFactory != NULL)
            execute = decoratorFactory->decorate(execute);
    return execute;
  }

}

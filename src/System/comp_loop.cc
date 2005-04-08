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
  } ;


  
  void execute_loop::execute(fact_db &facts) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;
    // initialize conditional variables to true
    *test = true ;
    
    param<int> time_var ;
    time_var = facts.get_variable(tvar) ;
    // Start iteration by setting iteration variable to zero 


    *time_var = 0 ;

    for(;;) { // Begin looping
      //First evaluate the collapse condition
      collapse->execute(facts) ;
      if(*test) {// If collapse condition satisfied, were finished
        return ;
      }
      // Advance to the next iterate
      advance->execute(facts) ;

      // We are finished with this iteration, so rotate variables and
      // add one to the iterate
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli)
        facts.rotate_vars(*rli) ;
      *time_var += 1 ;
    }
  }
  
  void execute_loop::Print(ostream &s) const {
    s << "Perform loop for time level "<< tlevel << endl ;
    s << "--compute collapse rule, conditional on " << cvar << endl ;
    s << "--collapse iteration {" << tlevel << "}"<< endl ;
    collapse->Print(s) ;
    s << "--advance iteration {" << tlevel << "}" << endl ;
    advance->Print(s) ;
    s << "end of loop for time level " << tlevel << endl ;
  }
  
  digraph::vertexSet visit_vertices(digraph dg,digraph::vertexSet begin) {
    
    digraph::vertexSet visit = begin ;
    digraph::vertexSet visited ;
    digraph::vertexSet::const_iterator ni ;
    
    // keep visiting vertices until no new vertices are found
    while(visit != EMPTY) {
      digraph::vertexSet newvertices ;
      // visit all the vertices that this graph leads to
      for(ni=visit.begin();ni!=visit.end();++ni)
        newvertices += dg[*ni] ;
      // update visit, but don't re-visit a vertex
      visit = newvertices - visited ;
      visited = visited + newvertices ;
    }
    return visited ;
  }
  
  inline bool offset_sort(const variable &v1, const variable &v2)
  { return v1.get_info().offset > v2.get_info().offset ; }
  
  loop_compiler::loop_compiler(rulecomp_map &rule_process, digraph dag, int id):cid(id)
  {
    ////////////////////
    // store the graph structure and the relevant rulecompiler map
    loop_gr = dag ;
    ruleSet allrules = extract_rules(loop_gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
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
      exit(-1) ;
    }
  
    tlevel = loopset.begin()->source_time() ;
    loop_gr.remove_vertex((*loopset.begin()).ident()) ;
    loop_gr.remove_dangling_vertices() ;

    digraph loop_grt = loop_gr.transpose() ;
  
    if(collapse_rules.size() != 1 ||
       collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
      cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
           << endl << "error is not recoverable" << endl ;
      cerr << "rules = " << collapse_rules << endl ;
      exit(-1) ;
    }

    variableSet all_vars = extract_vars(loop_gr.get_all_vertices()) ;
    for(variableSet::const_iterator vi=all_vars.begin();vi!=all_vars.end();++vi)
      if(vi->get_info().offset == 1)
        advance_vars += *vi ;
  
    collapse_vars = collapse_rules.begin()->targets() ;
    cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

    collapse_gr = loop_gr.subgraph(visit_vertices(loop_grt,collapse_vars)+collapse_vars) ;
    digraph::vertexSet collapse_rulesV = collapse_gr.get_all_vertices() & interval(UNIVERSE_MIN,-1) ;
    digraph::vertexSet advance_subset = loop_gr.get_all_vertices() - collapse_rulesV ;
    advance_gr = loop_gr.subgraph(advance_subset) ;
    advance_gr.remove_dangling_vertices() ;

    output_present = false ;

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
      exit(-1) ;
    }
  } 
  */

  void loop_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void loop_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    if(duplicate_work) {
      variableSet my_vars = advance_vars ;
      if(output_present)
	my_vars += output ;

      for(variableSet::const_iterator vi = my_vars.begin();
	  vi != my_vars.end(); vi++) {
	scheds.add_policy(*vi, sched_db::NEVER);
	variable tmp_var(*vi, vi->time());
	scheds.add_policy(tmp_var, sched_db::NEVER);
      }
    
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
	list<variable>::const_iterator li ;
	for(li=rli->begin();li!=rli->end();++li) {
	  scheds.add_policy(*li, sched_db::NEVER);
	}
      }
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
    if(output_present)
      var_requests += output ;
    
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

    return new execute_loop(cond_var,executeP(col),executeP(adv),tlevel,rotate_lists) ;
  }


}

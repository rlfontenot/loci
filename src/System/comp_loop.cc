#include "comp_tools.h"
#include <vector>
using std::vector ;
using std::list ;
using std::map ;

using std::ostream ;
using std::endl ;

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
      cvar(cv),collapse(col),advance(adv),tlevel(tl),rotate_lists(rl) {
      warn(col==0 || advance==0) ; tvar = variable(tlevel) ;
      control_thread = true ;}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;


  
  void execute_loop::execute(fact_db &facts) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;
    
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
  
  
  loop_compiler::loop_compiler(rulecomp_map &rule_process, digraph dag)  {

    ruleSet collapse_rules, all_rules, loopset ;
    variableSet collapse_vars ;

    std::vector<digraph::vertexSet> collapse_sched ;
    std::vector<digraph::vertexSet> advance_sched ;
    
    all_rules = extract_rules(dag.get_all_vertices()) ;
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
    dag.remove_vertex((*loopset.begin()).ident()) ;
    digraph dagt = dag.transpose() ;
  
    if(collapse_rules.size() != 1 ||
       collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
      cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
           << endl << "error is not recoverable" << endl ;
      cerr << "rules = " << collapse_rules << endl ;
      exit(-1) ;
    }

    variableSet all_vars = extract_vars(dag.get_all_vertices()) ;
    for(variableSet::const_iterator vi=all_vars.begin();vi!=all_vars.end();++vi)
      if(vi->get_info().offset == 1)
        advance_vars += *vi ;
  
    collapse_vars = collapse_rules.begin()->targets() ;
    cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

    // create Output Variable
    //    output = variable(variable("OUTPUT"),tlevel) ;
    //    variableSet outputSet ;
    //    outputSet = interval(output.ident(),output.ident()) ;

    // Schedule part of graph that leads to collapse
    collapse_sched = schedule_dag(dag, EMPTY,visit_vertices(dagt,collapse_vars)) ;
    compile_dag_sched(collapse_comp,collapse_sched,rule_process,dag) ;
    
    // Schedule advance part of loop.  First try to schedule any output, then
    // schedule the advance
    digraph::vertexSet visited ;
    for(int i = 0;i<collapse_sched.size();++i)
      visited += collapse_sched[i] ;
    //    vector<digraph::vertexSet>
    //      dag_sched = schedule_dag(dag,visited, visit_vertices(dagt,outputSet)) ;
    vector<digraph::vertexSet> dag_sched ;
    compile_dag_sched(advance_comp,dag_sched,rule_process,dag) ;
    
    if(dag_sched.size() == 0)
      output_present = false ;
    else 
      output_present = true ;

    for(int i=0;i<dag_sched.size();++i)
      advance_sched.push_back(dag_sched[i]) ;
  
    for(int i = 0;i<dag_sched.size();++i)
      visited += dag_sched[i] ;
    // now schedule everything that hasn't been scheduled
    dag_sched = schedule_dag(dag,visited) ;
    compile_dag_sched(advance_comp,dag_sched,rule_process,dag) ;
    
    for(int i=0;i<dag_sched.size();++i)
      advance_sched.push_back(dag_sched[i]) ;


    all_loop_vars = all_vars ;
    
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices = visited ;
    for(int i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
    if(allvertices != dag.get_all_vertices()) {
      cerr << " rules NOT scheduled = " << endl
           << extract_rules(dag.get_all_vertices() - allvertices) << endl ;
      cerr << " variables NOT scheduled = "
           << extract_vars(dag.get_all_vertices() - allvertices) << endl ;
      vector<digraph::vertexSet> components =
        component_sort(dag).get_components() ;
      for(int i=0;i<components.size();++i) {
        if(components[i].size() > 1) {
          cerr << "reason: graph not a dag, strongly connected component found"
               << endl ;
          cerr << "component rules = " << endl ;
          cerr << extract_rules(components[i]) << endl ;
        }
      }
      exit(-1) ;
    }
#endif
  }

  void loop_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
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
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi) ;
      scheds.variable_request(*vi,vexist) ;
    }
    
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=advance_comp.rbegin();ri!=advance_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
    for(ri=collapse_comp.rbegin();ri!=collapse_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
    
  }
  
  executeP loop_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {

    std::list<std::list<variable> > rotate_lists ;

    map<variable,list<variable> > vlist ;
    for(variableSet::const_iterator vi=all_loop_vars.begin();
        vi!=all_loop_vars.end();++vi) {
      if(vi->time() == tlevel && !vi->assign) {
        variable var_stationary(*vi,time_ident()) ;
        list<variable> s ;
        s.push_back(*vi) ;
        vlist[var_stationary].merge(s,offset_sort) ;
      }
    }
    map<variable,list<variable> >::const_iterator ii ;
    for(ii=vlist.begin();ii!=vlist.end();++ii) {
      if(ii->second.size() < 2) 
        continue ;
      std::list<variable>::const_iterator jj ;
      bool overlap = false ;
      variableSet vtouch ;
      for(jj=ii->second.begin();jj!=ii->second.end();++jj) {
        variableSet aliases = scheds.get_aliases(*jj) ;
        variableSet as = aliases ;
        variableSet::const_iterator vi ;
        
        for(vi=aliases.begin();vi!=aliases.end();++vi) 
          as += scheds.get_synonyms(*vi) ;
        if((as & vtouch) != EMPTY)
          overlap = true ;
        vtouch += as ;
      }
      if(!overlap) {
        rotate_lists.push_back(ii->second) ;
      } else {
        if(ii->second.size() !=2) {
          cerr << "unable to have history on variables aliased in time"
               << endl
               << "error occured on variable " << ii->first
               << "{" << tlevel << "}"
               << endl ;
          exit(-1) ;
        }
      }
    }
      
      
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

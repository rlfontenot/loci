#include "comp_tools.h"
#include <vector>
using std::vector ;

namespace Loci {

  void compile_dag_sched(std::vector<rule_compilerP> &dag_comp,
                         const std::vector<digraph::vertexSet> &dag_sched,
                         const rulecomp_map &rcm,
                         const digraph &dag) {
    digraph dagt = dag.transpose() ;
    for(int i=0;i<dag_sched.size();++i) {
      variableSet vars = extract_vars(dag_sched[i]) ;
      ruleSet rules = extract_rules(dag_sched[i]) ;
      if(rules == EMPTY && i+1<dag_sched.size()) {
        ++i ;
        vars += extract_vars(dag_sched[i]) ;
        rules = extract_rules(dag_sched[i]) ;
      }

      variableSet barrier_vars ;
      variableSet::const_iterator vi ;

      std::map<variable,ruleSet> barrier_info,reduce_info,singleton_info ;
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(dagt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        for(ri=var_rules.begin();ri!=var_rules.end();++ri)
          if(ri->get_info().rule_class != rule::INTERNAL) {
            use_rules += *ri ;
            rule_implP rimp = ri->get_rule_implP() ;
            if(rimp->get_rule_class() == rule_impl::POINTWISE)
              pointwise = true ;
            if(rimp->get_rule_class() == rule_impl::UNIT ||
               rimp->get_rule_class() == rule_impl::APPLY)
              reduction = true ;
            if(rimp->get_rule_class() == rule_impl::SINGLETON)
              singleton = true ;
          }
        WARN(reduction && pointwise || pointwise && singleton ||
             reduction && singleton) ;
        
        if((use_rules != EMPTY)) {
          if(pointwise)
            barrier_info[*vi] = use_rules ;
          if(reduction)
            reduce_info[*vi] = use_rules ;
          if(singleton)
            singleton_info[*vi] = use_rules ;
        }

      }

      dag_comp.push_back(new barrier_compiler(barrier_info)) ;

      if(singleton_info.begin() != singleton_info.end())
        dag_comp.push_back(new singleton_var_compiler(singleton_info)) ;
                           
      if(reduce_info.begin() != reduce_info.end()) {
        std::map<variable,ruleSet>::const_iterator xi ;
        variableSet vars ;
        for(xi=reduce_info.begin();xi!=reduce_info.end();++xi) {
          vars += xi->first ;
          rule unit_rule ;
          CPTR<joiner> join_op = CPTR<joiner>(0) ;
          storeRepP sp ;
          ruleSet::const_iterator ri ;
          for(ri=xi->second.begin();ri!=xi->second.end();++ri) {
            if(ri->get_rule_implP()->get_rule_class() == rule_impl::UNIT)
              unit_rule = *ri ;
            else if(ri->get_rule_implP()->get_rule_class() == rule_impl::APPLY){
              if(join_op == 0)
                join_op = ri->get_rule_implP()->get_joiner() ;
              else
                if(typeid(*join_op) !=
                   typeid(*(ri->get_rule_implP()->get_joiner()))) {
                  cerr << "Warning:  Not all apply rules for variable " << xi->first << " have identical join operations!" << endl ;
                }
            } else {
              cerr << "Warning: reduction variable " << xi->first
                   << " has a non-reduction rule contributing to its computation,"
                   << endl << "offending rule is " << *ri << endl ;
            }
          }
          FATAL(join_op == 0) ;
          sp = join_op->getTargetRep() ;
          if(sp->RepType() == PARAMETER)
            dag_comp.push_back(new reduce_param_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          else
            dag_comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                         join_op)) ;
        }
      }

      if(rules != EMPTY) {
        ruleSet::const_iterator ri ;
        for(ri=rules.begin();ri!=rules.end();++ri) {
          rulecomp_map::const_iterator rmi ;
          rmi = rcm.find(*ri) ;
          FATAL(rmi == rcm.end()) ;
          dag_comp.push_back(rmi->second) ;
        }
      }
    }
  }

  
  dag_compiler::dag_compiler(rulecomp_map &rule_process, digraph dag) {

    std::vector<digraph::vertexSet> dag_sched = schedule_dag(dag) ;
    compile_dag_sched(dag_comp,dag_sched,rule_process,dag) ;
    
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
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts) ;
  }

  void dag_compiler::process_var_requests(fact_db &facts) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts) ;
  }

  executeP dag_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> elp = new execute_list ;

    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts)) ;
    }

    return executeP(elp) ;
  }
}


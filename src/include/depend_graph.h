#ifndef DEPEND_GRAPH_H
#define DEPEND_GRAPH_H

#include "digraph.h"
#include "rule.h"

namespace Loci {

  inline ruleSet extract_rules(digraph::nodeSet ns) {
    return ruleSet(ns & interval(UNIVERSE_MIN,-1)) ;
  }

  inline variableSet extract_vars(digraph::nodeSet ns) {
    return variableSet(ns & interval(0,UNIVERSE_MAX)) ;
  }


  class dependency_graph {
    digraph gr ;
    variableSet invoke_rule(rule f) ;
    void promote_variable(variable v1, variable v2) ;
    void generalize_variable(variable v1,variable v2) ;
    void compose_graph(rule_db &rdb, variableSet given) ;
    void remove_incorrect_time_promotions() ;
    void create_looping_rules() ;
    void clean_graph(variableSet given, variableSet target) ;
  public:
    dependency_graph(rule_db &rdb, variableSet given, variableSet target) ;
    digraph get_graph() const { return gr; } 
  } ;

}  


#endif

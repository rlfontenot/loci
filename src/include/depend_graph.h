#ifndef DEPEND_GRAPH_H
#define DEPEND_GRAPH_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#include <Tools/digraph.h>
#include <Tools/cptr.h>
#include <variable.h>
#include <rule.h>

namespace Loci {

  inline ruleSet extract_rules(digraph::vertexSet ns) {
    return ruleSet(ns & interval(UNIVERSE_MIN,-1)) ;
  }

  inline variableSet extract_vars(digraph::vertexSet ns) {
    return variableSet(ns & interval(0,UNIVERSE_MAX)) ;
  }

  class dependency_graph {
    digraph gr ;
    void create_looping_rules() ;
    void clean_graph(variableSet given, variableSet target) ;
  public:
    dependency_graph(rule_db &rdb, variableSet given, variableSet target) ;
    digraph get_graph() const { return gr; } 
  };

}

#endif

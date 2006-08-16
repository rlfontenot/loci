#ifndef DEPEND_GRAPH_H
#define DEPEND_GRAPH_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

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

  // experimental version
  class dependency_graph2 {
    digraph gr ;
  public:
    dependency_graph2(const rule_db& rdb,
                      const variableSet& given,
                      const variableSet& target) ;
    
    digraph get_graph() const { return gr; } 
  };

  digraph partition_iteration(digraph gr) ;

  void clean_graph(digraph &gr,
                   const variableSet& given,
                   const variableSet& target) ;
}

#endif

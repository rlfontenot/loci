#ifndef SCHED_MLG_H
#define SCHED_MLG_H

#include "sched_tools.h"

// MultiLevel Graph Class

namespace Loci {
  struct multiLevelGraph {
    struct subGraph {
      digraph::vertexSet graph_v ;
      digraph::vertexSet incoming_v ;
      digraph::vertexSet outgoing_v ;
      digraph gr ;
      subGraph() {}
      subGraph(digraph &gin, digraph::vertexSet &grvtx) 
      { create(gin,grvtx) ; }
      void create(digraph &gin, digraph::vertexSet &grvtx) ;
    } ;
    digraph::vertexSet subgraphs ;
    std::map<int,subGraph>  subgrmap ;
    int toplevel ;
    int super_node_number ;

    multiLevelGraph() {toplevel = 0; super_node_number = 0 ;}
    multiLevelGraph(digraph gr, digraph::vertexSet grvtx) ;
    rule make_super_rule(variableSet sources, variableSet targets,
                         variable cond = variable()) ;
    void insert(int id,subGraph g) 
    { warn(find(id) != 0) ; subgrmap[id] = g ; subgraphs += id ; }
    subGraph *find(int id) ;
    int mksnode(int gr_id, digraph::vertexSet grvtx, variable cond_var = variable()) ;
  } ;

}

#endif

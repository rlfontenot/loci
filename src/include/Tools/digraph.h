#ifndef DIGRAPH_H
#define DIGRAPH_H

#include <Tools/intervalSet.h>
#include <Tools/Handle.h>
#include <map>
#include <vector>

namespace Loci {


  // Store a representation of a directed graph.  The directed graph is
  // represented by a vector of nodSets.  Each element in the vector corresponds
  // to a node, and the nodeSet is a set of edges that lead out of that
  // node.
  class digraph {
  public:
    typedef intervalSet nodeSet ;
  private:
    class digraphRep {
    public:
      typedef std::map<int,nodeSet> graph_matrix ;
      graph_matrix graph ;
      nodeSet source_nodes ;  
      int max_node ;
    
      digraphRep() ;
      ~digraphRep() ;
      void add_edge(int i, int j)
        { graph[i] += j ;  source_nodes += i ;
        max_node = max(max_node,max(i,j));}

      void add_edges(int i, const nodeSet &ns) {
        nodeSet &gedges = graph[i] ;
        nodeSet sum = gedges | ns ;
        if(sum != gedges) 
          gedges = sum ;
        source_nodes += i ;
        max_node = max(max_node,max(i,ns.Max())) ;
      }

      void add_edges(const nodeSet &ns, int j) ;

      void remove_edge(int i, int j)
        { graph[i] -= j ; if(graph[i] == EMPTY) source_nodes -= i ; }

      void remove_edges(int i, const nodeSet &ns)
        { graph[i] -= ns ; if(graph[i] == EMPTY) source_nodes -= i ; }
      void remove_node(int i) ;
      void remove_nodes(const nodeSet &ns) { subgraph(~ns) ; }

      void add_graph(const digraphRep &gr) ;
      void subtract_graph(const digraphRep &gr) ;
      void subgraph(const nodeSet &ns) ;
    
      const nodeSet &get_edges(int i) const {
        graph_matrix::const_iterator ii = graph.find(i) ;
        if(ii != graph.end())
          return ii->second ;
        else
          return EMPTY ;
      }
      nodeSet get_source_nodes() const { return source_nodes; }
      nodeSet get_target_nodes() const {
        nodeSet target_nodes = EMPTY ;
        graph_matrix::const_iterator ii ;
        for(ii=graph.begin();ii!=graph.end();++ii)
          target_nodes += ii->second ;
        return target_nodes ;
      }
    } ;
    Handle<digraphRep> Rep ;   // graph representation
    Handle<digraphRep> RepT ;  // representation transpose
    friend class digraphRep ;
  public:
    digraph() ;
    ~digraph() ;
    void add_edge(int i, int j) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->add_edge(i,j) ;
      RepT->add_edge(j,i) ;
    }
    void add_edges(const nodeSet &ns,int j) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->add_edges(ns,j) ;
      RepT->add_edges(j,ns) ;
    }
    void add_edges(int i, const nodeSet &ns) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->add_edges(i,ns) ;
      RepT->add_edges(ns,i) ;
    }
    void remove_edge(int i,int j) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->remove_edge(i,j) ;
      RepT->remove_edge(j,i) ;
    }
    void remove_node(int i) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->remove_node(i) ;
      RepT->remove_node(i) ;
    }
    void remove_nodes(const nodeSet &ns) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->remove_nodes(ns) ;
      RepT->remove_nodes(ns) ;
    }
    void add_graph(const digraph &gr) {
      Rep.MakeUnique();
      RepT.MakeUnique() ;
      Rep->add_graph(*(gr.Rep)) ;
      RepT->add_graph(*(gr.RepT)) ;
    }
    void subtract_graph(const digraph &gr) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->subtract_graph(*(gr.Rep)) ;
      RepT->subtract_graph(*(gr.RepT)) ;
    }
    const nodeSet &operator[](int i) const { return get_edges(i); }
    const nodeSet &get_edges(int i) const { return Rep->get_edges(i) ; }
    bool is_edge(int i, int j) const { return get_edges(i).inSet(j) ; }
    nodeSet get_source_nodes() const { return Rep->get_source_nodes(); }
    nodeSet get_target_nodes() const { return RepT->get_source_nodes(); }
    nodeSet get_all_nodes() const
      { return Rep->get_source_nodes() + RepT->get_source_nodes() ; }
    int max_node() const { return Rep->max_node ; }
    digraph transpose() const {
      digraph dg ;
      dg.Rep = RepT ;
      dg.RepT = Rep ;
      return dg ;
    }
    digraph subgraph(const nodeSet &ns) const {
      digraph dg = *this ;
      dg.Rep.MakeUnique() ;
      dg.RepT.MakeUnique() ;
      dg.Rep->subgraph(ns) ;
      dg.RepT->subgraph(ns) ;
      return dg ;
    }
  
  } ;

  inline digraph &operator+=(digraph &g1,const digraph &g2)
    { g1.add_graph(g2) ; return g1 ; }

  inline digraph &operator-=(digraph &g1,const digraph &g2)
    { g1.subtract_graph(g2) ; return g1 ; }


  // This class creates a topological sort of the nodes in a directed graph
  // while also creating a 
  class component_sort {
    sequence node_list ;
    std::vector<digraph::nodeSet> components ;
  public:
    component_sort(const digraph &dg) ;
    sequence node_order() const { return node_list ; }
    const std::vector<digraph::nodeSet> &get_components() { return components; }
  } ;


}

#endif

#ifndef DIGRAPH_H
#define DIGRAPH_H

#include <Tools/intervalSet.h>
#include <Tools/Handle.h>
//#include <map>
#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif

#include <vector>

namespace Loci {


  // Store a representation of a directed graph.  The directed graph is
  // represented by a vector of nodSets.  Each element in the vector corresponds
  // to a node, and the nodeSet is a set of edges that lead out of that
  // node.
  class digraph {
  public:
    typedef intervalSet vertexSet ;
    typedef vertexSet nodeSet ;
  private:
    class digraphRep {
    public:
      typedef std::hash_map<int,vertexSet> graph_matrix ;
      graph_matrix graph ;
      vertexSet source_vertices ;  
      int max_vertex ;
    
      digraphRep() ;
      ~digraphRep() ;
      void add_edge(int i, int j)
        { graph[i] += j ;  source_vertices += i ;
        max_vertex = max(max_vertex,max(i,j));}

      void add_edges(int i, const vertexSet &ns) {
        vertexSet &gedges = graph[i] ;
        vertexSet sum = gedges | ns ;
        if(sum != gedges) 
          gedges = sum ;
        source_vertices += i ;
        max_vertex = max(max_vertex,max(i,ns.Max())) ;
      }

      void add_edges(const vertexSet &ns, int j) ;

      void remove_edge(int i, int j)
        { graph[i] -= j ; if(graph[i] == EMPTY) source_vertices -= i ; }

      void remove_edges(int i, const vertexSet &ns)
        { graph[i] -= ns ; if(graph[i] == EMPTY) source_vertices -= i ; }
      void remove_vertex(int i) ;
      void remove_vertices(const vertexSet &ns) { subgraph(~ns) ; }

      void add_graph(const digraphRep &gr) ;
      void subtract_graph(const digraphRep &gr) ;
      void subgraph(const vertexSet &ns) ;
    
      const vertexSet &get_edges(int i) const {
        graph_matrix::const_iterator ii = graph.find(i) ;
        if(ii != graph.end())
          return ii->second ;
        else
          return EMPTY ;
      }
      vertexSet get_source_vertices() const { return source_vertices; }
      vertexSet get_target_vertices() const {
        vertexSet target_vertices = EMPTY ;
        graph_matrix::const_iterator ii ;
        for(ii=graph.begin();ii!=graph.end();++ii)
          target_vertices += ii->second ;
        return target_vertices ;
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
    void add_edges(const vertexSet &ns,int j) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->add_edges(ns,j) ;
      RepT->add_edges(j,ns) ;
    }
    void add_edges(int i, const vertexSet &ns) {
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
    void remove_vertex(int i) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->remove_vertex(i) ;
      RepT->remove_vertex(i) ;
    }
    void remove_vertices(const vertexSet &ns) {
      Rep.MakeUnique() ;
      RepT.MakeUnique() ;
      Rep->remove_vertices(ns) ;
      RepT->remove_vertices(ns) ;
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
    const vertexSet &operator[](int i) const { return get_edges(i); }
    const vertexSet &get_edges(int i) const { return Rep->get_edges(i) ; }
    bool is_edge(int i, int j) const { return get_edges(i).inSet(j) ; }
    vertexSet get_source_vertices() const { return Rep->get_source_vertices(); }
    vertexSet get_target_vertices() const { return RepT->get_source_vertices(); }
    vertexSet get_all_vertices() const
      { return Rep->get_source_vertices() + RepT->get_source_vertices() ; }
    int max_vertex() const { return Rep->max_vertex ; }
    digraph transpose() const {
      digraph dg ;
      dg.Rep = RepT ;
      dg.RepT = Rep ;
      return dg ;
    }
    digraph subgraph(const vertexSet &ns) const {
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


  // This class creates a topological sort of the vertices in a directed graph
  // while also creating a 
  class component_sort {
    sequence vertex_list ;
    std::vector<digraph::vertexSet> components ;
  public:
    component_sort(const digraph &dg) ;
    sequence vertex_order() const { return vertex_list ; }
    const std::vector<digraph::vertexSet> &get_components() { return components; }
  } ;


}

#endif

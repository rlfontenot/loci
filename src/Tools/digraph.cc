
#include <Tools/digraph.h>
#include <Tools/stream.h>
#include <Tools/debug.h>

using std::map ;
using std::vector ;

namespace Loci {
  
  digraph::digraph() {}
  digraph::~digraph() {}

  digraph::digraphRep::digraphRep() {max_vertex = 0 ;}

  digraph::digraphRep::~digraphRep() {}

  void digraph::digraphRep::add_edges(const vertexSet &ns, int j) 
    {
      max_vertex = max(max_vertex,max(ns.Max(),j)) ;
      vertexSet::const_iterator i ;
      for(i=ns.begin();i!=ns.end();++i)
        graph[*i] += j ;
      source_vertices += ns ;
    }


  void digraph::digraphRep::remove_vertex(int i)
    {
      if(source_vertices.inSet(i)) {
        source_vertices -= i ;
        graph[i] = EMPTY ;
      }

      graph_matrix::iterator ii ;
      for(ii=graph.begin();ii!=graph.end();++ii) {
        ii->second -= i ;

        if(ii->second == EMPTY)
          source_vertices -= ii->first ;
      }
    }

  void digraph::digraphRep::add_graph(const digraphRep &gr) {
    graph_matrix::const_iterator ii ;
    for(ii=gr.graph.begin();ii!=gr.graph.end();++ii)
      add_edges(ii->first,ii->second) ;
  }

  void digraph::digraphRep::subtract_graph(const digraphRep &gr) {
    graph_matrix::const_iterator ii ;
    for(ii=gr.graph.begin();ii!=gr.graph.end();++ii) 
      remove_edges(ii->first,ii->second) ;
  }

  void digraph::digraphRep::subgraph(const vertexSet &ns) {
    graph_matrix::iterator ii ;
    source_vertices = source_vertices & ns ;
    vertexSet empty_sources ;
    for(ii=graph.begin();ii!=graph.end();++ii) 
      if(!ns.inSet(ii->first))
        ii->second = EMPTY ;
      else {
        vertexSet isect = ii->second & ns ;
        if(ii->second != isect)
          ii->second = isect ;
        if(isect == EMPTY)
          empty_sources += ii->first ;
      }
    source_vertices -= empty_sources ;
  }

  // Perform topological sort on the name_tags based on the graph described by
  // registered function dependencies.  This algorithm is described in
  // "Intoduction to Algorithms" by T Cormen, C Leiserson, and R Rivest.

  namespace {
    enum vertex_color { WHITE, GRAY, BLACK } ;
    typedef std::map<int,vertex_color> vertex_color_map ;
    typedef map<int,sequence> priority_graph ;
    typedef digraph::vertexSet vertexSet ;
  }



  vertexSet dfs_visit(int vertex, const digraph &pg,
                    sequence &vertex_list,vertex_color_map &vertex_color)
    {
      vertex_color[vertex] = GRAY ;
      const vertexSet &edges = pg[vertex] ;

      vertexSet visited ;
      visited += vertex ;
      for(vertexSet::const_iterator ii=edges.begin();ii!=edges.end();++ii) {
        switch(vertex_color[*ii]) {
        case GRAY:
          // Backedge, indicates recursion
          //            back_edges.add_edge(vertex,*ii) ;
          //            back_edges_exist = true ;
          break ;
        case BLACK:
          // Already Searched so skip
          break ;
        case WHITE:
          // Not visited, continue depth first search
          visited += dfs_visit(*ii,pg,vertex_list,vertex_color) ;
          break ;
        default:
#ifdef DEBUG
          cerr << "Default case in dfs_visit should not occur" << endl ;
#endif
          break ;
        }
      }
      vertex_color[vertex] = BLACK ;
      vertex_list += vertex ;
      return visited ;
    }


  component_sort::component_sort(const digraph &dg)
    {
      if(dg.get_source_vertices() == EMPTY)
        return ;
      vertex_color_map vertex_color ;
      vertexSet all_vertices = dg.get_source_vertices() + dg.get_target_vertices() ;
      vertexSet::const_iterator ii ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii)
        vertex_color[*ii] = WHITE ;

      sequence order ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii) {
        switch(vertex_color[*ii]) {
        case WHITE:
          dfs_visit(*ii,dg,order,vertex_color) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid vertex color in switch" << endl ;
        }
      }
      vertex_color.clear() ;

      sequence rorder = order.Reverse() ;
      vector<int> omap ;  // order map
      map<int,int> imap ; // inverse
      for(sequence::const_iterator is=rorder.begin();is!=rorder.end();++is) {
        imap[*is] = omap.size() ;
        omap.push_back(*is) ;
      }

      digraph dgt ;
      const vertexSet &sources = dg.get_source_vertices() ;
      for(ii=sources.begin();ii!=sources.end();++ii) {
        const vertexSet &edges = dg[*ii] ;
        vertexSet::const_iterator jj ;
        warn(imap.find(*ii) == imap.end()) ;
        for(jj=edges.begin();jj!=edges.end();++jj) {
          warn(imap.find(*jj) == imap.end()) ;
          dgt.add_edge(imap[*jj],imap[*ii]) ;
        }
      }
      for(int i=0;i<omap.size();++i)
        vertex_color[i] = WHITE ;

      all_vertices = interval(0,omap.size()-1) ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii)
        vertex_color[*ii] = WHITE ;

      order = EMPTY ;
      vector<vertexSet> comp_vec ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii) {
        vertexSet n ;
        switch(vertex_color[*ii]) {
        case WHITE:
          comp_vec.push_back(dfs_visit(*ii,dgt,order,vertex_color)) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid vertex color in switch" << endl ;
        }
      }
      vector<vertexSet>::const_iterator ci ;
      for(ci=comp_vec.begin();ci!=comp_vec.end();++ci) {
        vertexSet n ;
        for(ii=ci->begin();ii!=ci->end();++ii) {
          if(*ii < omap.size())
            n += omap[*ii] ;
          else 
            cerr << "confused , *ii = " << *ii << endl ;
        }
        components.push_back(n) ;
      }

      for(sequence::const_iterator is=order.begin();is!=order.end();++is) {
        vertex_list += omap[*is] ;
      }  
    }
}

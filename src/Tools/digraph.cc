#include "digraph.h"
#include "stream.h"

#include <debug.h>

using std::map ;
using std::vector ;

namespace Loci {
  
  digraph::digraph() {}
  digraph::~digraph() {}

  digraph::digraphRep::digraphRep() {max_node = 0 ;}

  digraph::digraphRep::~digraphRep() {}

  void digraph::digraphRep::add_edges(const nodeSet &ns, int j) 
    {
      max_node = max(max_node,max(ns.Max(),j)) ;
      nodeSet::const_iterator i ;
      for(i=ns.begin();i!=ns.end();++i)
        graph[*i] += j ;
      source_nodes += ns ;
    }


  void digraph::digraphRep::remove_node(int i)
    {
      if(source_nodes.inSet(i)) {
        source_nodes -= i ;
        graph[i] = EMPTY ;
      }

      graph_matrix::iterator ii ;
      for(ii=graph.begin();ii!=graph.end();++ii) {
        ii->second -= i ;

        if(ii->second == EMPTY)
          source_nodes -= ii->first ;
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

  void digraph::digraphRep::subgraph(const nodeSet &ns) {
    graph_matrix::iterator ii ;
    source_nodes = source_nodes & ns ;
    nodeSet empty_sources ;
    for(ii=graph.begin();ii!=graph.end();++ii) 
      if(!ns.inSet(ii->first))
        ii->second = EMPTY ;
      else {
        nodeSet isect = ii->second & ns ;
        if(ii->second != isect)
          ii->second = isect ;
        if(isect == EMPTY)
          empty_sources += ii->first ;
      }
    source_nodes -= empty_sources ;
  }

  // Perform topological sort on the name_tags based on the graph described by
  // registered function dependencies.  This algorithm is described in
  // "Intoduction to Algorithms" by T Cormen, C Leiserson, and R Rivest.

  namespace {
    enum node_color { WHITE, GRAY, BLACK } ;
    typedef std::map<int,node_color> node_color_map ;
    typedef map<int,sequence> priority_graph ;
    typedef digraph::nodeSet nodeSet ;
  }



  nodeSet dfs_visit(int node, const digraph &pg,
                    sequence &node_list,node_color_map &node_color)
    {
      node_color[node] = GRAY ;
      const nodeSet &edges = pg[node] ;

      nodeSet visited ;
      visited += node ;
      for(nodeSet::const_iterator ii=edges.begin();ii!=edges.end();++ii) {
        switch(node_color[*ii]) {
        case GRAY:
          // Backedge, indicates recursion
          //            back_edges.add_edge(node,*ii) ;
          //            back_edges_exist = true ;
          break ;
        case BLACK:
          // Already Searched so skip
          break ;
        case WHITE:
          // Not visited, continue depth first search
          visited += dfs_visit(*ii,pg,node_list,node_color) ;
          break ;
        default:
#ifdef DEBUG
          cerr << "Default case in dfs_visit should not occur" << endl ;
#endif
          break ;
        }
      }
      node_color[node] = BLACK ;
      node_list += node ;
      return visited ;
    }


  component_sort::component_sort(const digraph &dg)
    {
      if(dg.get_source_nodes() == EMPTY)
        return ;
      node_color_map node_color ;
      nodeSet all_nodes = dg.get_source_nodes() + dg.get_target_nodes() ;
      nodeSet::const_iterator ii ;
      for(ii = all_nodes.begin() ; ii != all_nodes.end() ; ++ii)
        node_color[*ii] = WHITE ;

      sequence order ;
      for(ii = all_nodes.begin() ; ii != all_nodes.end() ; ++ii) {
        switch(node_color[*ii]) {
        case WHITE:
          dfs_visit(*ii,dg,order,node_color) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid node color in switch" << endl ;
        }
      }
      node_color.clear() ;

      sequence rorder = order.Reverse() ;
      vector<int> omap ;  // order map
      map<int,int> imap ; // inverse
      for(sequence::const_iterator is=rorder.begin();is!=rorder.end();++is) {
        imap[*is] = omap.size() ;
        omap.push_back(*is) ;
      }

      digraph dgt ;
      const nodeSet &sources = dg.get_source_nodes() ;
      for(ii=sources.begin();ii!=sources.end();++ii) {
        const nodeSet &edges = dg[*ii] ;
        nodeSet::const_iterator jj ;
        warn(imap.find(*ii) == imap.end()) ;
        for(jj=edges.begin();jj!=edges.end();++jj) {
          warn(imap.find(*jj) == imap.end()) ;
          dgt.add_edge(imap[*jj],imap[*ii]) ;
        }
      }
      for(int i=0;i<omap.size();++i)
        node_color[i] = WHITE ;

      all_nodes = interval(0,omap.size()-1) ;
      for(ii = all_nodes.begin() ; ii != all_nodes.end() ; ++ii)
        node_color[*ii] = WHITE ;

      order = EMPTY ;
      vector<nodeSet> comp_vec ;
      for(ii = all_nodes.begin() ; ii != all_nodes.end() ; ++ii) {
        nodeSet n ;
        switch(node_color[*ii]) {
        case WHITE:
          comp_vec.push_back(dfs_visit(*ii,dgt,order,node_color)) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid node color in switch" << endl ;
        }
      }
      vector<nodeSet>::const_iterator ci ;
      for(ci=comp_vec.begin();ci!=comp_vec.end();++ci) {
        nodeSet n ;
        for(ii=ci->begin();ii!=ci->end();++ii) {
          if(*ii < omap.size())
            n += omap[*ii] ;
          else 
            cerr << "confused , *ii = " << *ii << endl ;
        }
        components.push_back(n) ;
      }

      for(sequence::const_iterator is=order.begin();is!=order.end();++is) {
        node_list += omap[*is] ;
      }  
    }
}

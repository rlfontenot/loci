#include "visitor.h"
#include "visit_tools.h"
#include "comp_tools.h"
#include <Tools/stream.h>

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;

namespace Loci {

  namespace {
    
    void create_digraph_dot_file(const digraph &dg, const char* fname)
    {
      digraph dgt = dg.transpose() ;
      digraph::vertexSet allvertices = dg.get_all_vertices() ;
      digraph::vertexSet::const_iterator ri ;
      ofstream outf(fname) ;
      
      outf<<"digraph G {\n" ;
      //outf<<"size = \"8.5,11\";\n" ;

      // variable's ident() >= 0 rule's ident() < 0
      for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
        digraph::vertexSet outvertices = dg[*ri] ;
        digraph::vertexSet incomevertices = dgt[*ri] ;

        if(*ri < 0) {
          rule r(*ri) ;
          if(r.type() == rule::INTERNAL) {
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
          }
          else {
            outf<<"\""<<pretty_sig(r)<<"\""<<"[shape=box,style=filled,color=gold];\n" ;
          }
        }
        else {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
        }
        if(incomevertices == EMPTY) {
          if(*ri >= 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=red];\n" ;
          }
        }
        if(outvertices == EMPTY) {
          if(*ri >= 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=blue];\n" ;
          }
          continue ;
        }
        
        if(*ri < 0) {
          rule r(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r2(*ii) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""<<";\n" ;
            }else {
              variable v(*ii) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v<<"\""<<";\n" ;
            }
          }
        }else {
          variable v(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r)<<"\""<<";\n" ;
            }else {
              variable v2(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""<<";\n" ;
            }
          }        
        }
      }
      outf<<"}\n" ;
      outf.close() ;
    }

  } // end of unnamed namespace

  /////////////////////////////////////////////////////////////
  // graphVisualizeVisitor
  /////////////////////////////////////////////////////////////
  void graphVisualizeVisitor::visit(loop_compiler& lc) {
    if(Loci::MPI_rank==0) {
      cerr << "the looping part graph (collapse and advance), id: "
           << lc.cid << endl ;
      create_digraph_dot_file(lc.loop_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;

      cerr << "\tthe collapse part graph" << endl ;
      create_digraph_dot_file(lc.collapse_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
      
      cerr << "\tthe advance part graph" << endl ;
      create_digraph_dot_file(lc.advance_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

  void graphVisualizeVisitor::visit(dag_compiler& dc) {
    if(Loci::MPI_rank==0) {
      cerr << "the dag graph, id: " << dc.cid << endl ;
      create_digraph_dot_file(dc.dag_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

  void graphVisualizeVisitor::visit(conditional_compiler& cc) {
    if(Loci::MPI_rank==0) {
      cerr << "the conditional graph, id: " << cc.cid << endl ;
      create_digraph_dot_file(cc.cond_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

  ///////////////////////////////////////////////////////////////
  // dagCheckVisitor
  ///////////////////////////////////////////////////////////////
  ostream& dagCheckVisitor::visualize(ostream& s) const {
    if(Loci::MPI_rank == 0) {
      s << "visualizing detected cycle..." << endl ;
      create_digraph_dot_file(cycle,"cycle_in_dag.dot") ;
      system("dotty cycle_in_dag.dot") ;
      system("rm -fr cycle_in_dag.dot") ;
    }
    return s ;
  }

  ///////////////////////////////////////////////////////////////
  // chompRuleVisitor
  ///////////////////////////////////////////////////////////////
  ostream& chompRuleVisitor::visualize(ostream& s) const {
    if(Loci::MPI_rank == 0) {
      if(all_chains.empty()) {
        s << "NO chomping chains found!" << endl ;
        return s ;
      }
      map<int,list<chomp_chain> >::const_iterator mi ;
      list<chomp_chain>::const_iterator li ;
      int total_chomp_vars = 0 ;
      for(mi=all_chains.begin();mi!=all_chains.end();++mi) {
        list<chomp_chain> cclist = mi->second ;
        s << "There are " << cclist.size()
          << " chomping rule chain(s) in super node SN"
          << mi->first << ": " << endl ;
        int i ;
        for(li=cclist.begin(),i=1;li!=cclist.end();++li,++i) {
          int total = li->second.size() ;
          total_chomp_vars += total ;
          s << "Visualizing the NO" << i << " chomping chain..." << endl ;
          s << "variables that can be chomped in this chain: "
            << li->second << " total: " << total << endl ;
          
          s << "source variables to the chomp node: " << extract_vars(li->first.get_source_vertices() - li->first.get_target_vertices()) << endl ;
          s << "target variables of the chomp node: " << extract_vars(li->first.get_target_vertices() - li->first.get_source_vertices()) << endl ;
          
          create_digraph_dot_file(li->first,"chomping_rules.dot") ;
          system("dotty chomping_rules.dot") ;
          system("rm -fr chomping_rules.dot") ;
        }
      }
      s << "total variables can be chomped in the program: "
        << total_chomp_vars << endl ;
    }
    return s ;
  }

} // end of namespace Loci

#include "sched_tools.h"
#include "sched_mlg.h"

#define DEVELOP
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

namespace Loci {

  inline multiLevelGraph create_mlg(digraph dg,
                                    digraph::vertexSet sources,
                                    digraph::vertexSet targets) {
    digraph::vertexSet grvtx = dg.get_all_vertices() ;
    grvtx -= sources ;
    grvtx -= targets ;
    return multiLevelGraph(dg,grvtx) ;
  }

  
  void prepare_digraph(digraph &dg) {
    digraph::vertexSet grvtx = dg.get_all_vertices() ;
    ruleSet all_rules = extract_rules(grvtx) ;
    ruleSet::const_iterator ri ;
    variableSet reduce_vars ;
    for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
      // collect information about unit rules
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
        if(ri->targets().size() != 1) {
          cerr << "unit rule must have only one variable as target: "
               << endl << *ri << endl ;
          exit(-1) ;
        }
        variable unit_var = *(ri->targets().begin()) ;
        if(reduce_vars.inSet(unit_var)) {
          cerr << "only one unit rule may be defined for reduction variable "
               << unit_var << endl ;
          exit(-1) ;
        }
        reduce_vars += unit_var ;

        // Make unit rule loop on itself so that it is included in
        // components as we decompose the graph (the unit rule
        // needs to stay with its apply rules.
        dg.add_edge(unit_var.ident(),(*ri).ident()) ;

        // Now we make all apply rules depend explicitly on the unit rule
        // This will make sure that when we begin to schedule that we
        // always call the unit rule before any of the apply rules.
        // We also do a few sanity checks here as well.
        ruleSet apply_rules = extract_rules(dg.transpose()[unit_var.ident()]) ;
        apply_rules -= *ri ;
        for(ruleSet::const_iterator rii=apply_rules.begin();
            rii!=apply_rules.end();
            ++rii) {
          if(rii->get_info().rule_impl->get_rule_class() != rule_impl::APPLY) {
            cerr << "pointwise rule is producing reduction variable: "
                 << unit_var << endl ;
            cerr << "This rule should be an apply rule. (offending rule is):"
                 << endl << *rii << endl ;
            exit(-1) ;
          }
          if(rii->targets().size() != 1) {
            cerr << "apply rule should apply to only one reduction variable"
                 << endl << "offending rule is:" << endl << *rii << endl ;
            exit(-1) ;
          }
          // Here we add the dependency
          dg.add_edge(unit_var.ident(),(*rii).ident()) ;
        }
      }

    }
    
  }
  
  digraph::vertexSet sort_time_hierarchy(multiLevelGraph &mlg) {
    digraph::vertexSet new_vertices ;
    
    int toplevel = mlg.toplevel ;
    new_vertices += toplevel ;
    multiLevelGraph::subGraph &sg = *mlg.find(toplevel) ;
    digraph &dg = sg.gr ;
    variableSet vars = extract_vars(sg.graph_v) ;

    map<time_ident,digraph::vertexSet> time_sort_vertices ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
      time_sort_vertices[vi->time()] += (*vi).ident() ;
    ruleSet rules = extract_rules(dg.get_all_vertices()) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      // sort rules based on time level
      const time_ident rule_time =
        ri->type()==rule::COLLAPSE?ri->source_time():ri->target_time() ;
      time_sort_vertices[rule_time] += (*ri).ident() ;
    }

    // compute a depth first ordering of time levels
    vector<time_ident> time_stack,time_sequence ;
    time_stack.push_back(time_ident()) ;
    for(;!time_stack.empty();) {
      time_ident ctime = time_stack.back() ;
      time_stack.pop_back() ;
      time_sequence.push_back(ctime) ;
      const vector<time_ident> &children = ctime.children()  ;
      for(vector<time_ident>::const_iterator
            i=children.begin();i!=children.end();++i)
        time_stack.push_back(*i) ;
    }

    // proceed from the bottom of the time hierarchy tree, extracting each
    // time level and constructing a super node based on that time level.
    digraph main_graph = dg ;
    for(vector<time_ident>::reverse_iterator ti = time_sequence.rbegin();
        ti!=time_sequence.rend();
        ++ti) {
      if(time_sort_vertices[*ti] != EMPTY) {
        if(*ti != time_ident()) {
          int new_node = mlg.mksnode(toplevel,time_sort_vertices[*ti]) ;
          new_vertices += new_node;
          // The new node that we create will be in the level
          // above our time iteration, so insert it in that
          // iteration set.
          time_sort_vertices[*(ti+1)] += new_node ;
        }
      }
    }
    return new_vertices ;
  }

  digraph::vertexSet partition_loops(multiLevelGraph &mlg, int supernode,
                                     digraph::vertexSet &loops) {
    multiLevelGraph::subGraph &sg = *mlg.find(supernode) ;

    // First search for looping and collapse rules in this subgraph
    digraph::vertexSet looping, collapse, unit, apply, new_nodes ;

    ruleSet all_rules = extract_rules(sg.graph_v) ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        looping += ri->ident() ;
      if(ri->type() == rule::COLLAPSE) 
        collapse += ri->ident() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)
        unit += ri->ident() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
        apply += ri->ident() ;
    }

    WARN(looping == EMPTY && collapse != EMPTY) ;
    WARN(looping.size() > 1) ;

    if(looping.size() > 1) {
      cerr << extract_rules(looping) << endl ;
    }


    // Now find any recursive rules and form supernodes for these components
    // We first have to temporarily remove the looping rules so that they
    // don't appear to be big collections of recursive rules.  
    digraph tempg = sg.gr ;
    if(looping != EMPTY)
      tempg.remove_vertices(looping) ;


    // We use a component sort to identify recursive rules.
    vector<digraph::vertexSet> recurse_components =
      component_sort(tempg).get_components() ;
    vector<digraph::vertexSet>::const_iterator ci ;
   

    // Here we check all of the components (larger than on vertex).  We
    // eliminate the reduction rules (which can't be part of a recursive block)
    for(ci=recurse_components.begin();ci!=recurse_components.end();++ci)
      if(ci->size() > 1) {
        ruleSet rr = extract_rules(*ci) ;
        digraph::vertexSet sn = rr ;
        if((sn & unit) == EMPTY || (sn & apply) == EMPTY) {
          // Found a recursive rule, so make a supernode for it.
          int newnode = mlg.mksnode(supernode,sn) ;
          new_nodes += newnode ;
        } else {
          sn -= unit ;
          sn -= apply ;
          if(sn != EMPTY) {
            cerr << "reductions cannot be part of a recursive block!" << endl ;
            cerr << "offending rules = " << endl ;
            cerr << rr << endl ;
          }
        }
      }

    if(looping == EMPTY)
      return EMPTY ;

    // Add these edges so that the collapse rules will be
    // included in the loop component.
    sg.gr.add_edges(collapse,*(looping.begin())) ;
    
    vector<digraph::vertexSet> components =
      component_sort(sg.gr).get_components() ;

    for(ci=components.begin();ci!=components.end();++ci)
      if(ci->size() > 1) {
        if((*ci & looping) != EMPTY) {
          int newnode = mlg.mksnode(supernode,*ci) ; 
          new_nodes += newnode ;
          loops += newnode ;
          break ;
        }
      }
    return new_nodes ;
  }
  
  void test_decompose(digraph dg,
                      digraph::vertexSet sources,
                      digraph::vertexSet targets) {
    prepare_digraph(dg) ;
    multiLevelGraph mlg = create_mlg(dg,sources,targets) ;

    digraph::vertexSet new_vertices = sort_time_hierarchy(mlg) ;
    digraph::vertexSet::const_iterator vi ;
    digraph::vertexSet looping_components ;
    for(vi=new_vertices.begin();vi!=new_vertices.end();++vi) {
      digraph::vertexSet new_loops ;
      looping_components += partition_loops(mlg,*vi,new_loops) ;
    }
  }
    
    

}

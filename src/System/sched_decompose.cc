#include "sched_tools.h"
#include "distribute.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::pair ;
using std::make_pair ;

namespace Loci {
  void cleanup_component(const digraph &dg, digraph::vertexSet &component) {
    // Remove variables from component if they have references outside
    // the component
    variableSet cvars = extract_vars(component) ;
    digraph::vertexSet remove_vars ;
    digraph gr = dg ;
    digraph grt = gr.transpose() ;
    for(variableSet::const_iterator vi = cvars.begin();
        vi!=cvars.end();
        ++vi) {
      const int id = (*vi).ident() ;
      if((grt[id]-component) != EMPTY
         || (gr[id]-component) != EMPTY)
        remove_vars += id ;
    }

    component -= remove_vars ;
#ifdef VERBOSE
    variableSet remove_vset = extract_vars(remove_vars) ;
    debugout << "remove_vars = " << remove_vset ;
    debugout <<endl ;
#endif
  }

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

  void post_process_mlg(multiLevelGraph::subGraph &sg) {
    digraph &gr = sg.gr ;
    ruleSet gr_rules = extract_rules(sg.graph_v) ;
    ruleSet::const_iterator ri ;
    for(ri=gr_rules.begin();ri!=gr_rules.end();++ri) {
      if((ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)) {
        variable unit_var = *(ri->targets().begin()) ;
        gr.remove_edge(unit_var.ident(),ri->ident()) ;
        ruleSet reduce_rules = extract_rules(gr.transpose()[unit_var.ident()]) ;
        reduce_rules -= *ri ;
        ruleSet::const_iterator rri ;
        for(rri = reduce_rules.begin();rri != reduce_rules.end();++rri) {
          gr.remove_edge(unit_var.ident(),rri->ident()) ;
          gr.add_edge(ri->ident(),rri->ident())  ;
        }
      }
    }
    // remove extra dangling nodes
    gr.remove_dangling_vertices() ;
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
      cleanup_component(dg,time_sort_vertices[*ti]) ;
      if(time_sort_vertices[*ti] != EMPTY) {
        if(*ti != time_ident()) {
          int new_node = mlg.mksnode(toplevel,time_sort_vertices[*ti]) ;
          new_vertices += new_node;
          // The new node that we create will be in the level
          // above our time iteration, so insert it in that
          // iteration set.
          time_sort_vertices[ti->parent()] += new_node ;

          ruleSet temp = extract_rules(time_sort_vertices[*ti]) ;
        }
      }
    }
    return new_vertices ;
  }

  digraph::vertexSet partition_loops(multiLevelGraph &mlg, int supernode,
                                     digraph::vertexSet &loops,
                                     digraph::vertexSet &recursive) {
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
    if(looping == EMPTY && collapse != EMPTY) {
      ruleSet collapse_rules = extract_rules(collapse) ;
      cerr << "A collapse rule exists without a well formed loop" << endl ;
      cerr << "collapse = " << collapse_rules << endl ;
    }
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
        if((sn & unit) == EMPTY && (sn & apply) == EMPTY) {
          // Found a recursive rule, so make a supernode for it.
          int newnode = mlg.mksnode(supernode,sn) ;
          new_nodes += newnode ;
          recursive += newnode ;
          // Edit graph so that the resulting supernode is no longer recursive
          digraph::vertexSet recursive_bits = sg.gr[newnode] & sg.gr.transpose()[newnode] ;

          // remove the recursive dependency and replace with a dependency
          // that will make sure that the recursive rule block is executed last.
          digraph::vertexSet::const_iterator cvi ;
          digraph::vertexSet recurse_depend ;
          for(cvi=recursive_bits.begin();cvi!=recursive_bits.end();++cvi) {
            recurse_depend += sg.gr.transpose()[*cvi] ;
            sg.gr.remove_edge(*cvi,newnode) ;
          }
          recurse_depend -= newnode ;
            
          sg.gr.add_edges(recurse_depend,newnode);
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
          digraph::vertexSet component_parts = *ci ;
          //          cleanup_component(sg.gr,component_parts) ;
          int newnode = mlg.mksnode(supernode,component_parts) ; 
          new_nodes += newnode ;
          loops += newnode ;
          break ;
        }
      }
    return new_nodes ;
  }


  digraph::vertexSet decompose_conditional_rules(multiLevelGraph &mlg,
                                                int supernode)
  {
    multiLevelGraph::subGraph &sg = *mlg.find(supernode) ;


    map<variable,ruleSet> cond_rules ;
    ruleSet all_rules = extract_rules(sg.graph_v) ;

    digraph::vertexSet looping ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        looping += ri->ident() ;
      if(ri->get_info().desc.conditionals != EMPTY) {
        variableSet conds = ri->get_info().desc.conditionals ;
        if(conds.size() != 1) {
          cerr << "improper rule: " << *ri << endl
               << "Rule Improperly specifies more than one conditional" << endl ;
          cerr << "Error not recoverable." << endl ;
          Loci::Abort() ;
        }
        variable cond = *(conds.begin()) ;
        cond_rules[cond] += *ri ;
      }
    }
    // If no conditional rules return 
    if(cond_rules.begin() == cond_rules.end()) {
      return EMPTY ;
    }

    digraph::vertexSet new_rules ;

    digraph::vertexSet::const_iterator vi ;

    digraph tmpgr = sg.gr ;

    tmpgr.remove_vertices(looping) ;

    // Remove rule-rule dependencies.  These are not real variable
    // dependencies, but rather created by rename rules.
    ruleSet graph_rules = extract_rules(tmpgr.get_all_vertices()) ;

    for(ruleSet::const_iterator ri = graph_rules.begin();
        ri!=graph_rules.end();
        ++ri) {
      int id = (*ri).ident() ;
      digraph::vertexSet rm_v = tmpgr[id] & interval(UNIVERSE_MIN,-1) ;
      if(rm_v != EMPTY) {
        for(digraph::vertexSet::const_iterator ii = rm_v.begin() ;
            ii != rm_v.end() ;
            ++ii) {
#ifdef VERBOSE
          debugout << "removing edge " << *ri << "," << rule(*ii) << endl  ;
#endif
          tmpgr.remove_edge(id,*ii) ;
        }
      }
    }
    
    // We need to make sure that apply rules work as a group.  We use
    // a component sort to group them together and remove them
    // from the graph.  We will replace them after we've decided
    // what will be in the conditional component.
    vector<digraph::vertexSet> cs = component_sort(tmpgr).get_components() ;
    map<int,digraph::vertexSet> cm ;
    digraph::vertexSet incm ;

    for(size_t i=0;i<cs.size();++i)
      if(cs[i].size() > 1) {
        int new_vertex = tmpgr.max_vertex()+1 ;
        digraph::vertexSet in ;
        digraph::vertexSet out ;
        cm[new_vertex] = cs[i] ;
        cs[i] = extract_rules(cs[i]) ;
        for(vi=cs[i].begin();vi!=cs[i].end();++vi) {
          out += tmpgr[*vi] ;
          in += tmpgr.transpose()[*vi] ;
        }
        incm += new_vertex ;

        out -= cs[i] ;
        in -= cs[i] ;

        in -= out ;// Disable recursive loops

        tmpgr.remove_vertices(cs[i]) ;
        tmpgr.add_edges(new_vertex,out) ;
        tmpgr.add_edges(in,new_vertex) ;
      }
  
  
      
    // Loop over each conditional variable and find the part of the graph
    // that exclusively connects to the conditional components.  All of
    // these rules can be treated together when evaluating the conditional
    map<variable,ruleSet>::const_iterator mi ;
    for(mi=cond_rules.begin();mi!= cond_rules.end();++mi) {

      digraph tr = tmpgr.transpose() ;
      digraph::vertexSet component ;
      // We start with the conditional rules
      digraph::vertexSet working = mi->second ;

      while(working != EMPTY) {
        // Repeatedly search for new vertices in the graph that are exclusively
        // associated with the conditional rules.  As we find them add them
        // to the current set of component rules
        component += working ;

        // First find candidates by following all edges that lead to the
        // current working set.
        digraph::vertexSet candidates ;
        for(vi=working.begin();vi!=working.end();++vi)
          candidates += tr[*vi] ;

        // Create a new working set by searching through the candidates to
        // find those verticies that only refer to vertices already found
        // in the component.

        working = EMPTY ;
        for(vi=candidates.begin();vi!=candidates.end();++vi) {
          if((tmpgr[*vi] & component) == tmpgr[*vi]) {
            working += *vi ;
          }
        }
        // We make sure that we recursively visit the same components by
        // removing them from our next working set.
        working -= component ;
        // We also make sure that we don't follow the conditional variable.
        // This needs to be computed external to the conditional block
        working -= mi->first.ident() ;

      }

      // If any apply rules were included, add them back into the current
      // set of components

      digraph::vertexSet compcm = component & incm ;
      component -= incm ;
      for(vi=compcm.begin();vi!=compcm.end();++vi)
        component += cm[*vi] ;


      // Make sure that we don't go outside of the current subgraph
      component &= sg.graph_v ;

      // Remove variables from component if they have references outside
      // the component
      cleanup_component(sg.gr,component) ;

      int new_node =  mlg.mksnode(supernode,component,mi->first) ;

      
      new_rules += new_node ;
      
    }

    return new_rules ;
  }
  
  decomposed_graph::decomposed_graph(digraph dg,
                                     digraph::vertexSet sources,
                                     digraph::vertexSet targets) {
    prepare_digraph(dg) ;
    mlg = create_mlg(dg,sources,targets) ;

#ifdef DEBUG
    digraph::vertexSet all_vertices = mlg.find(mlg.toplevel)->graph_v ;
#endif

    digraph::vertexSet new_vertices = sort_time_hierarchy(mlg) ;
    digraph::vertexSet::const_iterator vi ;
    digraph::vertexSet looping_components ;
    
    for(vi=new_vertices.begin();vi!=new_vertices.end();++vi) {
      looping_components += partition_loops(mlg,*vi,loops,recursive) ;
    }
    digraph::vertexSet working = mlg.subgraphs ;
    digraph::vertexSet conditional_components;
    for(vi=working.begin();vi!=working.end();++vi) {
      conditional += decompose_conditional_rules(mlg,*vi) ;
    }
    
    for(vi=mlg.subgraphs.begin();vi!=mlg.subgraphs.end();++vi)
      post_process_mlg(*mlg.find(*vi)) ;
  
#ifdef DEBUG
    digraph::vertexSet visit ;
    for(vi=mlg.subgraphs.begin();vi!=mlg.subgraphs.end();++vi) {
      digraph::vertexSet vl = mlg.find(*vi)->graph_v ;
      if((vl & visit) != EMPTY) {
        cerr << "overlapping graphs in mlg" << endl ;
        cerr << "overlapping variables = " << extract_vars(vl&visit)<< endl ;
        cerr << "overlapping rules = " << extract_rules(vl&visit)<<endl ;
      }
      visit += vl ;
    }

    if((all_vertices - visit) != EMPTY) {
      cerr << "some verticies lost in mlg generation" << endl ;
      cerr << "lost variables = " << extract_vars(all_vertices-visit) << endl ;
      cerr << "lost rules = " << extract_rules(all_vertices-visit) << endl ;
    }
#endif
    
  }

}

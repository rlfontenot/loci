#include "sched_tools.h"
#include "sched_mlg.h"

//#define DEVELOP
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
    
    
  decompose_graph::decompose_graph(digraph dg,
                                   digraph::vertexSet sources,
                                   digraph::vertexSet targets) {
#ifdef DEVELOP
    test_decompose(dg,sources,targets) ;
#endif
    // Add vertices to the graph so that the source variables and target variables
    // aren't left dangling.
    variableSet vars = extract_vars(dg.get_all_vertices()) ;
    start = variable("__START__") ;
    finish = variable("__FINISH__") ;
    dg.add_edges(start.ident(),sources) ;
    dg.add_edges(targets,finish.ident()) ;

    // In the first step we decompose the graph by sorting the graph
    // into iteraton levels.
    map<time_ident,digraph::vertexSet> time_sort_vertices ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
      time_sort_vertices[vi->time()] += (*vi).ident() ;
    ruleSet rules = extract_rules(dg.get_all_vertices()) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      // sort rules based on time level
      const time_ident rule_time =
        ri->type()==rule::COLLAPSE?ri->source_time():ri->target_time() ;
      time_sort_vertices[rule_time] += (*ri).ident() ;
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
        reduce_info &rinfo = reduce_map[unit_var] ;
        rinfo.unit_rule = *ri ;
        dg.add_edge(unit_var.ident(),(*ri).ident()) ;
        rinfo.apply_rules = extract_rules(dg.transpose()[unit_var.ident()]) ;
        rinfo.apply_rules -= rinfo.unit_rule ;
        for(ruleSet::const_iterator
              rii=rinfo.apply_rules.begin();rii!=rinfo.apply_rules.end();++rii) {
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
          dg.add_edge(unit_var.ident(),(*rii).ident()) ;
        }
      }
      // collect information about looping rules
      if(ri->get_info().qualifier() == "looping")
        looping_rules += *ri ;
      // collect information about rules that are conditionally executed.
      if(ri->get_info().desc.conditionals != EMPTY) {
        variableSet conds = ri->get_info().desc.conditionals ;
        if(conds.size() != 1) {
          cerr << "improper rule: " << *ri << endl
               << "Rule Improperly specifies more than one conditional" << endl ;
          cerr << "Error not recoverable." << endl ;
          exit(-1) ;
        }
        variable cond = *(conds.begin()) ;

        if(ri->type() != rule::COLLAPSE &&
           (ri->targets().size() != 1 ||
            (ri->targets().begin())->get_info().name != string("OUTPUT"))) {
          cerr << "conditional rules must either be collapse rules or " << endl
               << "have a OUTPUT as the sole argument" << endl ;
          cerr << "offending rule: " <<*ri << endl ; 
          exit(-1) ;
        }

        conditional_rules += *ri ;
        conditional_vars += cond ;
        conditional_map[cond] += *ri ;
      }
    }

    // Create apply rule to unit rule mapping
    for(variableSet::const_iterator
          vi = reduce_vars.begin(); vi != reduce_vars.end(); ++vi) {

      reduce_info &rinfo = reduce_map[*vi] ;

      for(ruleSet::const_iterator
            it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) {
        apply_to_unit[*it] = rinfo.unit_rule ;
      }
    }

#ifdef VERBOSE
    cout << "reduction variables = " << reduce_vars << endl;
    cout << "looping rules = " << endl << looping_rules << endl ;
    cout << "conditional rules = " << endl << conditional_rules << endl ;
#endif
  
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
    for(vector<time_ident>::reverse_iterator
          ti = time_sequence.rbegin();ti!=time_sequence.rend();++ti) 
      if(time_sort_vertices[*ti] != EMPTY) {
        rule svertex = create_supernode(main_graph,time_sort_vertices[*ti]) ;
        if(*ti != time_ident()) {
          const time_ident rule_time =svertex.target_time() ;
          time_sort_vertices[rule_time] += svertex.ident() ;
        } else
          top_level_rule = svertex ; 

        if((svertex.sources() & svertex.targets()) != EMPTY) {
          cerr << "warning, time level " << svertex << " is recursive." << endl ;
          cerr << "build and collapse rules not allowed to form a "
               << "recursive block" << endl ;
          exit(-1) ;
        }
      }

    // Go over each time level and decompose it into strongly connected components
    // At this level, loops remain as a single component.
    ruleSet time_levels = supernodes ;
    for(ruleSet::const_iterator
          ri=time_levels.begin();ri!=time_levels.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      if(*ri != top_level_rule) {
        ruleSet loop_rules = ruleSet(sni.graph.get_all_vertices() & looping_rules) ;
      
        if(loop_rules.size() != 1) {
          time_ident ti ;
          variableSet vars = extract_vars(sni.graph.get_all_vertices()) ;
          for(variableSet::const_iterator vi= vars.begin();vi!=vars.end();++vi) 
            if(ti.before(vi->time()))
              ti = vi->time() ;
        
          cerr << "incorrect specification of iteration on time level " << ti
               << endl << "no loop specification was produced." << endl
               << "Perhaps no advance rules could be applied at this level?"
               << endl ;
          exit(-1) ;
        }
        // Insure that collapse rules are included in loop supernode
        sni.graph.add_edges(sni.targets,(*(loop_rules.begin())).ident()) ;
      }
      vector<digraph::vertexSet> components =
        component_sort(sni.graph).get_components() ;
      vector<digraph::vertexSet>::const_iterator ci ;
      for(ci=components.begin();ci!=components.end();++ci)
        if(ci->size() > 1) {
          variableSet vars = variableSet(*ci & reduce_vars) ;
          ruleSet loop_rules = ruleSet(sni.graph.get_all_vertices() &
                                       looping_rules) ;
          if(vars== EMPTY || loop_rules != EMPTY)
            create_supernode(sni.graph,*ci) ;
          else {
            if(vars.size() > 1) {
              cerr << "reductions can't mix with recursions." << endl ;
              cerr << "error occured when parsing variables " << vars << endl ;
              exit(-1) ;
            }
            // remove the self-referential loops and make all pointwise rules
            // depend on the unit rule.
            variable rvar = *(vars.begin()) ;
            reduce_info &rinfo = reduce_map[rvar] ;

            for(ruleSet::const_iterator
                  it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
              sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
            sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
            sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
          }
        }
    }

    ruleSet new_rules = ruleSet(supernodes - time_levels) ;
    // Now go over remaining supernodes to decompose time loops
    for(ruleSet::const_iterator
          ri=new_rules.begin();ri!=new_rules.end();++ri) {
      super_node_info &sni = supermap[*ri] ;

      ruleSet loop_rules = ruleSet(sni.graph.get_all_vertices() & looping_rules) ;

      digraph gtemp = sni.graph ;

      warn(loop_rules.size() > 1) ;

      if(loop_rules.size() == 1)
        gtemp.remove_vertex((*loop_rules.begin()).ident()) ;
      else continue ;
    
      vector<digraph::vertexSet> components =
        component_sort(gtemp).get_components() ;
      vector<digraph::vertexSet>::const_iterator ci ;
      for(ci=components.begin();ci!=components.end();++ci)
        if(ci->size() > 1) {
          variableSet vars = variableSet(*ci & reduce_vars) ;
          if(vars== EMPTY)
            create_supernode(sni.graph,*ci) ;
          else {
            if(vars.size() > 1) {
              cerr << "reductions can't mix with recursions." << endl ;
              cerr << "error occured when parsing variables " << vars << endl ;
              exit(-1) ;
            }
            // remove the self-referential loops and make all pointwise rules
            // depend on the unit rule.
            variable rvar = *(vars.begin()) ;
            reduce_info &rinfo = reduce_map[rvar] ;

            for(ruleSet::const_iterator
                  it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
              sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
            sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
            sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
          }
        }
    }
    // Create conditional execution super vertices.
    // A conditional rule may exclusively depend on other rules for execution
    // In this case, we only need to execute these other rules when the condition
    // is satisified.  This is accomplished by creating a conditional supernode
    // that encapsulates this grouped behavior.
    new_rules = supernodes ;
    ruleSet cond_supernodes  ;
    for(ruleSet::const_iterator
          ri = new_rules.begin();ri!=new_rules.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      digraph::vertexSet graph_vertices = sni.graph.get_all_vertices() &
        ~sni.sources & ~sni.targets ;
    
      variableSet cond_vars = variableSet(graph_vertices & conditional_vars) ;
      // If there are conditional variables in this graph, then decompose them
      if(cond_vars != EMPTY) {
        // loop over conditional variables and extract the functions that
        // are contingent on them
        variableSet::const_iterator vi ;
        for(vi = cond_vars.begin();vi!=cond_vars.end();++vi) {

          ruleSet cond_rules = conditional_map[*vi] ;
          // Note, all of the rules that are conditional of this variable
          // should be in this graph, if not report this error and exit.
          if((graph_vertices & cond_rules) != cond_rules) {
            cerr << "problem with conditional functions of variable " << *vi
                 << endl ;
            cerr << "rule(s) not in same level of graph:" << endl ;
            ruleSet except = ruleSet(cond_rules - (graph_vertices & cond_rules)) ;
            cerr << except << endl ;
            cerr << "error not recoverable." << endl ;
            cerr << "error occured while processing " << *ri << endl ;
            exit(-1) ;
          }
          // Find the set of vertices that are exclusive in generating the
          // rules conditional on this variable (*vi)
          digraph::vertexSet cond_vertices =
            visit_vertices_exclusive(sni.graph.transpose(), cond_rules) ;
          // remove the rules that are exclusive to the conditional evaluation
          // itself
          cond_vertices = cond_vertices -
            visit_vertices_exclusive(sni.graph.transpose(),
                                     interval((*vi).ident(),(*vi).ident())) ;

          // create a super vertex for this conditional set of rules
          cond_supernodes += create_supernode(sni.graph,cond_vertices,*vi) ;
        }
      }
    }
    // collect info on supernodes
    // remove recursive dependencies from rules
    ruleSet all_rules = supernodes ;
    for(ruleSet::const_iterator
          ri = supernodes.begin();ri!=supernodes.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      digraph::vertexSet graph_vertices = sni.graph.get_all_vertices() &
        ~sni.sources & ~sni.targets ;
      all_rules += extract_rules(graph_vertices) ;
      ruleSet loop_rules = ruleSet(sni.graph.get_all_vertices() & looping_rules) ;
      if(loop_rules != EMPTY)
        looping_supernodes += *ri ;
      if(loop_rules == EMPTY && (sni.sources & sni.targets) != EMPTY)
        recursive_supernodes += *ri ;
      else {
        // fix graph so that it is a directed acyclic graph.
        // At this stage the graph should have all cycles except
        // cycles of length two removed from the graph.  We now
        // remove these cycles and add the appropriate dependencies
        // to ensure that these cycles are scheduled last.
        digraph g = sni.graph ;
        digraph gt = sni.graph.transpose() ;
        ruleSet grules = extract_rules(graph_vertices - looping_rules) ;
        ruleSet cycle_rule ;
        ruleSet non_cycle_rule ;
        for(ruleSet::const_iterator
              rii = grules.begin(); rii != grules.end(); ++rii) {
          int id = (*rii).ident() ;
          if((g[id] & gt[id]) != EMPTY) {  // found a cycle
            if(!supernodes.inSet(*rii)) {
              cerr << "internal consistency error, offending rule = "
                   << *rii << endl ;
              digraph::vertexSet lp = g[id] & gt[id] ;
              cout << "loop = RULES:"<<endl << extract_rules(lp) <<endl
                   << "VARS: "
                   << extract_vars(lp) << endl ;
              //            exit(-1) ;
            }
          
            variableSet cvars = extract_vars(g[id] & gt[id]) ;
            cycle_rule += *rii ;
            // get other rules that generate cycle variables
            ruleSet other_rules ;
            for(variableSet::const_iterator
                  vi = cvars.begin();vi!=cvars.end();++vi)
              other_rules += extract_rules(gt[(*vi).ident()]) ;
            other_rules -= *rii ;
            // remove edges causing loop
            for(variableSet::const_iterator
                  vi = cvars.begin();vi!=cvars.end();++vi)
              sni.graph.remove_edge((*vi).ident(),id) ;
            // replace these edges with rule dependencies
            sni.graph.add_edges(other_rules,id) ;
          }
        }
        if((cycle_rule & non_cycle_rule) != EMPTY) {
          cerr << "internal consistency error, rule cycles collide for rules:"
               << endl ;
          cerr << extract_rules(cycle_rule & non_cycle_rule) ;
          exit(-1) ;
        }

      }
    }

    // Create schedule compilers for all graph vertices
    ruleSet working_set = ruleSet(all_rules - supernodes) ;
    for(ruleSet::const_iterator
          ri = working_set.begin();ri != working_set.end(); ++ri) {
      if(ri->type() != rule::INTERNAL) {
        if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
          rule_process[*ri] = new apply_compiler(this,*ri,apply_to_unit[*ri]) ;
        else
          rule_process[*ri] = new impl_compiler(this,*ri) ;
      } else if(ri->get_info().qualifier() == "promote")
        rule_process[*ri] = new promote_compiler(this,*ri) ;
      else if(ri->get_info().qualifier() == "generalize")
        rule_process[*ri] = new generalize_compiler(this,*ri) ;
      else if(ri->get_info().qualifier() == "priority")
        rule_process[*ri] = new priority_compiler(this,*ri) ;
      else
        rule_process[*ri] = new error_compiler ;
    }
    ruleSet reduction_corrections ;
    // now do same for supernodes
    for(ruleSet::const_iterator
          ri = recursive_supernodes.begin();ri!=recursive_supernodes.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      if((sni.targets & reduce_vars) != EMPTY) {
        // this is really a reduction block, so remove recursion dependencies
        // and save to schedule as a dag.
        warn((sni.targets & reduce_vars) != sni.targets) ;
        // remove the self-referential loops and make all pointwise rules
        // depend on the unit rule.
        warn(sni.targets.size() != 1) ;
        variable rvar = *(sni.targets.begin()) ;
        reduce_info &rinfo = reduce_map[rvar] ;
      
        for(ruleSet::const_iterator
              it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
          sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
        sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
        sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
        reduction_corrections += *ri ;
      
      } else {
        // recursive rule block
        ruleSet recurse_rules = extract_rules(sni.graph.get_all_vertices()) ;
        if(recurse_rules.size() == 1 &&
           recurse_rules.begin()->get_info().desc.constraints.size() == 0) {
          rule_process[*ri] =
            new impl_recurse_compiler(this,*(recurse_rules.begin())) ;
        } else {
          rule_process[*ri] = new recurse_compiler(this,recurse_rules) ;
        }
      }
    }
    recursive_supernodes -= reduction_corrections ;
  
    for(ruleSet::const_iterator
          ri = looping_supernodes.begin();ri!=looping_supernodes.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      rule_process[*ri] = new loop_compiler(this,sni.graph) ;
    }
    for(ruleSet::const_iterator
          ri = cond_supernodes.begin();ri!=cond_supernodes.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      variable cond_var = *(ri->get_info().desc.conditionals.begin()) ;
      rule_process[*ri] = new conditional_compiler(this,sni.graph,cond_var) ;
    }

    ruleSet dag_rules = ruleSet(supernodes -
                                (recursive_supernodes +
                                 looping_supernodes +
                                 cond_supernodes)) ;

    for(ruleSet ::const_iterator
          ri = dag_rules.begin();ri!=dag_rules.end();++ri) {
      super_node_info &sni = supermap[*ri] ;
      rule_process[*ri] = new dag_compiler(this,sni.graph) ;
    }

  }

  rule decompose_graph::create_supernode(digraph &g,
                                         const digraph::vertexSet vertices,
                                         variable cond_var) {
    // extract rules in this supernode 
    ruleSet srules = extract_rules(vertices) ;
    warn(srules == EMPTY) ;

    // include all variables that are sources or sinks to any rule in this
    // supernode in the supernode graph.
    digraph::vertexSet ns = vertices ;
    variableSet sources,targets,locals,local_sources,local_targets ;

    // Calculate all sources and targets for rules in this supernode
    digraph gt = g.transpose() ;
    for(ruleSet::const_iterator ii =srules.begin();ii!=srules.end();++ii) {
      int id = (*ii).ident() ;
      sources += extract_vars(gt[id]) ;
      targets += extract_vars(g[id]) ;
    }

    // Add these new vertices to the supernode graph
    ns += sources + targets ;

    variableSet all_vars = extract_vars(ns) ;
    // find all variables that are referenced exclusively from within
    // the supernode, that is find all local variables
    digraph::vertexSet not_ns = ~ns ;
    for(variableSet::const_iterator
          vi=all_vars.begin();vi!=all_vars.end();++vi) {
      int id = (*vi).ident() ;
      if((gt[id] & not_ns) == EMPTY)
        local_sources += *vi ;
      if((g[id] & not_ns) == EMPTY)
        local_targets += *vi ;
    }

    locals = local_sources & local_targets ;
    // remove local variables from the super vertex rule source and target lists
    sources -= local_sources ;
    targets -= local_targets ;

    super_node_info sninfo ;
  
    // create supernode info
    sninfo.graph = g.subgraph(ns) ;
    sninfo.graph.add_edges(start.ident(),sources) ;
    sninfo.graph.add_edges(targets,finish.ident()) ;
  
    sninfo.sources = sources ;
    sninfo.targets = targets ;
    // remove components of graph that comprise supernode from parent graph
    ns -= (sources + targets) ;
    g.remove_vertices( ns ) ;
    // create a rule for the supernode
    rule vertex_rule = make_super_rule(sources,targets,cond_var) ;
    // Add supernode rule to parent graph
    g.add_edges(sources,vertex_rule.ident()) ;
    g.add_edges(vertex_rule.ident(),targets) ;
  
    // keep a map of supernodes
    supermap[vertex_rule] = sninfo ;
    supernodes += vertex_rule ;

#ifdef VERBOSE
    cout << "vertex_rule = " << vertex_rule << endl ;
    cout << "locals  = " << locals << endl ;
    cout << "rules in block = " << endl << srules << endl ;
    if((sources & targets) != EMPTY)
      cout << "sources & targets = " << variableSet(sources & targets) << endl ;
#endif

    return vertex_rule ;
  }

  void decompose_graph::existential_analysis(fact_db &facts) {

    (rule_process[top_level_rule])->set_var_existence(facts) ;
    variableSet var_requests = top_level_rule.targets() ;
    variableSet::const_iterator vi ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = facts.variable_existence(*vi) ;
      facts.variable_request(*vi,vexist) ;
    }
    (rule_process[top_level_rule])->process_var_requests(facts) ;
  }

  executeP decompose_graph::execution_schedule(fact_db &facts, int nth) {

    CPTR<execute_list> schedule = new execute_list ;
    schedule->append_list(new allocate_all_vars) ;
    schedule->append_list(new execute_create_threads(nth)) ;
    executeP top_level_schedule = (rule_process[top_level_rule])->
      create_execution_schedule(facts) ;
    if(top_level_schedule == 0) 
      return executeP(0) ;
    schedule->append_list(top_level_schedule) ;
    schedule->append_list(new execute_destroy_threads) ;
    return executeP(schedule) ;
  }
  





}

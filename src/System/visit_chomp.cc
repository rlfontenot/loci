#include "visitor.h"
#include "visit_tools.h"
#include "comp_tools.h"
#include <Tools/stream.h>

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include<map>
using std::map ;
#include <list>
using std::list ;
#include <utility>
using std::pair ;
using std::make_pair ;
using std::min ;

// all the searching and graph
// editing algorithm for chomping 
namespace Loci {
  
  //////////////////////////////////////////////////////////////
  // some helper functions 
  //////////////////////////////////////////////////////////////
  namespace {
    // given a contrete rule, return all the sources
    // that have maps
    inline variableSet map_sources(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet sources ;
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        if(!si->mapping.empty())
          sources += si->var ;
      }
      return sources ;
    }
    
    // given a contrete rule, return all the targetss
    // that have maps
    inline variableSet map_targets(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(!si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }
    
    // given a contrete rule, return all the targets
    // that do not have maps
    inline variableSet non_map_targets(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }
    
    inline variableSet non_map_targets(const ruleSet::const_iterator& ri) {
      const rule_impl::info& rinfo = ri->get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }

    // return the vertexSet of the rule's sources, targets
    // and the rule itself
    inline digraph::vertexSet get_rule_vertexSet(const rule& r) {
      digraph::vertexSet ret ;
      variableSet vars ;
      vars += r.sources() ;
      vars += r.targets() ;
      ret = get_vertexSet(vars) ;
      ret += r.ident() ;
      return ret ;
    }

    digraph::vertexSet get_ruleSet_vertexSet(const ruleSet& rules) {
      digraph::vertexSet ret ;
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
        ret += get_rule_vertexSet(*ri) ;
      }
      return ret ;
    }
    
    inline bool thread_rule(const rule& r) {
      return r.get_rule_implP()->thread_rule() ;
    }

    // given a rule, return if there are any OUTPUT
    // in the targets
    bool has_output_in_targets(const rule& r) {
      variableSet targets = r.targets() ;
      variable output("OUTPUT") ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi) {
        variable var_stationary(*vi,time_ident()) ;
        if(var_stationary == output)
          return true ;
      }

      return false ;
    }

    // given a digraph and all the chomping candidate
    // variables inside it, this function will return
    // all chomping chains by merging variables together
    void get_chomp_chains(list<chomp_chain>& result,
                          const variableSet& chomp_vars,
                          const digraph& gr) {
      // we first get a topological order
      // of the chomp_vars 
      vector<digraph::vertexSet> gr_order =
        component_sort(gr).get_components() ;
      list<variable> chomp_vars_order ;
      for(vector<digraph::vertexSet>::const_iterator vi=gr_order.begin();
          vi!=gr_order.end();++vi) {
        FATAL(vi->size() != 1) ;
        int_type fv = *(vi->begin()) ;
        if(fv >= 0){
          variable v(fv) ;
          if(chomp_vars.inSet(v))
            chomp_vars_order.push_back(v) ;
        }
      }

      digraph grt = gr.transpose() ;
      variableSet good_chomp_vars = chomp_vars ;
      
      while(!chomp_vars_order.empty()) {
        // here is the basis, we start by selecting
        // a single chompable variable and form the
        // smallest chomp chain
        digraph::vertexSet all_vertices ;
        variable begin = chomp_vars_order.front() ;
        ruleSet all_rules = extract_rules(gr[begin.ident()] +
                                          grt[begin.ident()]) ;
        variableSet sub_chomp_vars ;
        sub_chomp_vars += begin ;
        all_vertices += get_ruleSet_vertexSet(all_rules) ;
        // remove the first variable from the list
        chomp_vars_order.pop_front() ;

        bool merge = true ;
        while(merge) {
          merge = false ;
          vector<list<variable>::iterator> remove ;
          for(list<variable>::iterator lvi=chomp_vars_order.begin();
              lvi!=chomp_vars_order.end();++lvi){
            ruleSet rvs_rules = extract_rules(grt[lvi->ident()]) ;
            ruleSet fwd_rules = extract_rules(gr[lvi->ident()]) ;
            ruleSet reachable_rules = ruleSet(rvs_rules + fwd_rules) ;
            
            if( (all_rules & reachable_rules) != EMPTY) {
              // test for self-cycle
              digraph::vertexSet cp_all_vertices = all_vertices ;
              cp_all_vertices += get_ruleSet_vertexSet(reachable_rules) ;
              digraph test_gr = gr.subgraph(cp_all_vertices) ;
              digraph::vertexSet sources = test_gr.get_source_vertices()
                - test_gr.get_target_vertices() ;
              digraph::vertexSet targets = test_gr.get_target_vertices()
                - test_gr.get_source_vertices() ;

              digraph::vertexSet internal = test_gr.get_source_vertices()
                & test_gr.get_target_vertices() ;

              FATAL(internal !=
                    (test_gr.get_all_vertices() - sources - targets)) ;
              
              variableSet problem_vars =
                variableSet(extract_vars(internal) - good_chomp_vars) ;
              
              if(!has_path(gr,targets,sources) &&
                 (problem_vars == EMPTY)
                 ) {
                sub_chomp_vars += *lvi ;
                all_rules += reachable_rules ;
                all_vertices = cp_all_vertices ;
                remove.push_back(lvi) ;
                merge = true ;
              }
            }
          } // end-of-for(chomp_vars_order)
          // remove those merged variables from the list
          for(vector<list<variable>::iterator>::size_type i=0;
              i!=remove.size();++i)
            chomp_vars_order.erase(remove[i]) ;
        }//end-of-while(merge)

        // put it into the result
        digraph sub_gr = gr.subgraph(all_vertices) ;
        // if the chain just contains one chompable variable
        // we need to test for self-cycle, if has, it is
        // then discarded
        if(sub_chomp_vars.size() == 1) {
          digraph::vertexSet sources = sub_gr.get_source_vertices()
            - sub_gr.get_target_vertices() ;
          digraph::vertexSet targets = sub_gr.get_target_vertices()
            - sub_gr.get_source_vertices() ;

          digraph::vertexSet internal = sub_gr.get_source_vertices()
            & sub_gr.get_target_vertices() ;
          
          FATAL(internal != (sub_gr.get_all_vertices() - sources - targets)) ;

          variableSet problem_vars =
            variableSet(extract_vars(internal) - good_chomp_vars) ;
          
          if(has_path(gr,targets,sources) ||
             (problem_vars != EMPTY)
             ){
            good_chomp_vars -= *(sub_chomp_vars.begin()) ;
            continue ;
          }
        }
        // we need to remove any chompable variables
        // that are the sources or targets of this
        // sub_gr. They are the cause of the self-cycles.
        // They cannot be merged into other chains because
        // they are now visible outside of this chain
        digraph::vertexSet sources = sub_gr.get_source_vertices()
          - sub_gr.get_target_vertices() ;
        digraph::vertexSet targets = sub_gr.get_target_vertices()
          - sub_gr.get_source_vertices() ;
        variableSet stv = extract_vars(sources+targets) ;
        variableSet intersection = variableSet(stv & chomp_vars) ;
        if(intersection != EMPTY) {
          for(variableSet::const_iterator rmvi=intersection.begin();
              rmvi!=intersection.end();++rmvi)
            chomp_vars_order.remove(*rmvi) ;
          good_chomp_vars -= intersection ;
        }

        result.push_back(make_pair(sub_gr,sub_chomp_vars)) ;      
      } // end-of-while(!chomp_vars_order.empty())
    }

  } // end of unamed namespace

  ///////////////////////////////////////////////////////////////
  // chompPPVisitor
  ///////////////////////////////////////////////////////////////
  //#define DISABLE_APPLY
  
  chompPPVisitor::chompPPVisitor(fact_db& fd,
                                 const map<int,variableSet>& rot_vt,
                                 const map<int,variableSet>& lsharedt,
                                 const variableSet& rv)
    :facts(fd),rename_vars(rv) {
    map<int,variableSet>::const_iterator mi ;
    for(mi=rot_vt.begin();mi!=rot_vt.end();++mi)
      rotate_vars += mi->second ;
    for(mi=lsharedt.begin();mi!=lsharedt.end();++mi)
      loop_shared_vars += mi->second ;
  }

  void chompPPVisitor::discover(const digraph& gr) {
    digraph grt = gr.transpose() ;
    // first we need to get all the contrete rules in the graph
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;
    ruleSet remove ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(is_internal_rule(*ri))
        remove += *ri ;
    rules -= remove ;

    // first get all possible variables that can be chomped
    variableSet allvars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      allvars += ri->sources() ;
      allvars += ri->targets() ;
    }

    // then we collect bad_vars...

    // first any variable that involves maps
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      bad_vars += map_sources(*ri) ;
      bad_vars += map_targets(*ri) ;
    }

    // then we get all the variables that are not suitable for
    // chomping, i.e. they are not STORE, they are in the loop
    // rotate list, loop shared variables, and variables that
    // don't have outgoing edges, and variables that connect
    // to any internal rules or rename variables
    // or variables that are generated in more that one levels
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      variableSet stvars ;
      stvars += ri->sources() ;
      stvars += ri->targets() ;
      for(variableSet::const_iterator vi=stvars.begin();
          vi!=stvars.end();++vi) {
        
        if(seen_vars.inSet(*vi)) {
          bad_vars += *vi ;
          continue ;
        }
        
        storeRepP srp = facts.get_variable(*vi) ;
        if(srp->RepType() != Loci::STORE) {
          bad_vars += *vi ;
          continue ;
        }
        
        if(rotate_vars.inSet(*vi)) {
          bad_vars += *vi ;
          continue ;
        }
        if(loop_shared_vars.inSet(*vi)) {
          bad_vars += *vi ;
          continue ;
        }

        if(rename_vars.inSet(*vi)) {
          bad_vars += *vi ;
          continue ;
        }
        
        digraph::vertexSet next_vertices = gr[vi->ident()] ;
        if(next_vertices == EMPTY) {
          bad_vars += *vi ;
          continue ;
        }

        digraph::vertexSet next_vertices_t = grt[vi->ident()] ;
        if(next_vertices_t == EMPTY) {
          bad_vars += *vi ;
          continue ;
        }

        ruleSet tmp = extract_rules(next_vertices) ;
        for(ruleSet::const_iterator rii=tmp.begin();
            rii!=tmp.end();++rii)
          if(is_internal_rule(*rii) || !thread_rule(*rii) ||
             has_output_in_targets(*rii)
#ifdef DISABLE_APPLY             
             || rii->get_info().rule_impl->get_rule_class() == rule_impl::APPLY) {
#else
            ){
#endif          
            bad_vars += *vi ;
            break ;
          }

        tmp = extract_rules(next_vertices_t) ;
        for(ruleSet::const_iterator rii=tmp.begin();
            rii!=tmp.end();++rii)
          if(is_internal_rule(*rii) || !thread_rule(*rii) ||
             has_output_in_targets(*rii)) {
            bad_vars += *vi ;
            break ;
          }
      }
    }
    // fill up...
    good_vars += allvars ;
    good_vars -= bad_vars ;
    seen_vars += allvars ;
  }

  void chompPPVisitor::visit(loop_compiler& lc) {
    discover(lc.collapse_gr) ;
    discover(lc.advance_gr) ;
  }

  void chompPPVisitor::visit(dag_compiler& dc) {
    discover(dc.dag_gr) ;
  }

  void chompPPVisitor::visit(conditional_compiler& cc) {
    discover(cc.cond_gr) ;
  }

  ///////////////////////////////////////////////////////////////
  // chompRuleVisitor
  ///////////////////////////////////////////////////////////////
  list<chomp_chain> chompRuleVisitor::find_chain(const digraph& gr) {
    variableSet allvars = extract_vars(gr.get_all_vertices()) ;
    // now we start the searching, we are first very greedy,
    // we group any potential variables that can be chomped,
    // and then later, we perform checkings and necessary
    // adjustments.

    variableSet cand_vars = variableSet(allvars & good_vars) ;

    list<chomp_chain> chomp_chain_list ;
    get_chomp_chains(chomp_chain_list,cand_vars,gr) ;

    return chomp_chain_list ;    
  } // end-of-find_chain function
  
    // edit the graph to have the chomp node,
  void chompRuleVisitor::edit_gr(digraph& gr,const list<chomp_chain>& cc,
                                 rulecomp_map& rcm) {
    if(cc.empty())
      return ;

    for(list<chomp_chain>::const_iterator li=cc.begin();li!=cc.end();++li) {
      digraph chomp_graph = li->first ;
      variableSet chomp_vars = li->second ;
      digraph::vertexSet chomp_vars_vertices = get_vertexSet(chomp_vars) ;
      
      digraph::vertexSet all_vertices = chomp_graph.get_all_vertices() ;
      ruleSet all_rules = extract_rules(all_vertices) ;
      digraph::vertexSet rules_vertices = get_vertexSet(all_rules) ;
      
      // nodes that lead to the constructed super node
      // and nodes that leave the super node
      digraph::vertexSet
        source_vars_vertices = chomp_graph.get_source_vertices() -
        chomp_graph.get_target_vertices() ;
      digraph::vertexSet
        target_vars_vertices = chomp_graph.get_target_vertices() -
        chomp_graph.get_source_vertices() ;
      
      // make a rule for the chomp_graph
      rule chomp_rule = create_rule(extract_vars(source_vars_vertices),
                                    extract_vars(target_vars_vertices),
                                    "CHOMP") ;

      rcm[chomp_rule] = new chomp_compiler(chomp_graph,
                                           chomp_vars,apply2unit) ;
      
      // the vertices to be taken out
      digraph::vertexSet takeout_vertices =
        chomp_vars_vertices + rules_vertices ;

      if(takeout_vertices != (all_vertices -
                              source_vars_vertices - target_vars_vertices)) {
        cerr << "WARNING: inconsistency in chomping graph editing!" << endl ;
        exit(-1) ;
      }

      // get other possible nodes (outside of the chomp_graph)
      // that lead to any internal vertices of the chomp graph

      digraph grt = gr.transpose() ;
      for(digraph::vertexSet::const_iterator vi=takeout_vertices.begin();
          vi!=takeout_vertices.end();++vi) {
        source_vars_vertices += grt[*vi] - all_vertices ;
      }

      // graph editing
      // take the chomp graph out from the original graph
      gr = gr.subgraph(gr.get_all_vertices() - takeout_vertices) ;

      // edit the graph again
      gr.add_edges(source_vars_vertices, chomp_rule.ident()) ;
      gr.add_edges(chomp_rule.ident(), target_vars_vertices) ;
    }
  }
  
  void chompRuleVisitor::visit(loop_compiler& lc) {
    list<chomp_chain> c ;

    map<rule,rule_compilerP> tmp ;
    
    c = find_chain(lc.collapse_gr) ;
    if(!c.empty()) {
      all_chains[-lc.cid] = c ;
      edit_gr(lc.collapse_gr,c,lc.rule_compiler_map) ;
      edit_gr(lc.loop_gr,c,tmp) ;
    }
    
    c = find_chain(lc.advance_gr) ;
    if(!c.empty()) {
      all_chains[lc.cid] = c ;
      edit_gr(lc.advance_gr,c,lc.rule_compiler_map) ;
      edit_gr(lc.loop_gr,c,tmp) ;
    }
  }

  void chompRuleVisitor::visit(dag_compiler& dc) {
    list<chomp_chain> c ;
    
    c = find_chain(dc.dag_gr) ;
    if(!c.empty()) {
      all_chains[dc.cid] = c ;
      edit_gr(dc.dag_gr,c,dc.rule_compiler_map) ;
    }
  }

  void chompRuleVisitor::visit(conditional_compiler& cc) {
    list<chomp_chain> c ;
    
    c = find_chain(cc.cond_gr) ;
    if(!c.empty()) {
      all_chains[cc.cid] = c ;
      edit_gr(cc.cond_gr,c,cc.rule_compiler_map) ;
    }
  }

  ///////////////////////////////////////////////////////////////
  // compChompVisitor
  ///////////////////////////////////////////////////////////////
  void compChompVisitor::schedule_chomp(chomp_compiler& chc) {
    chc.chomp_sched = orderVisitor::order_dag(chc.chomp_graph) ;
  }
  
  void compChompVisitor::compile_chomp(chomp_compiler& chc) {
    for(vector<digraph::vertexSet>::const_iterator vi=chc.chomp_sched.begin();
        vi!=chc.chomp_sched.end();++vi) {
      ruleSet rules = extract_rules(*vi) ;
      variableSet vars = extract_vars(*vi) ;
      if(rules == EMPTY && vi+1<chc.chomp_sched.end()) {
        ++vi ;
        vars += extract_vars(*vi) ;
        rules = extract_rules(*vi) ;
      }
      vars &= chc.chomp_vars ;
      
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri){
        // this barrier_compiler bc is just a fake compiler
        // to be used in the pair 
        rule_compilerP bc = new barrier_compiler(vars) ;
        chc.chomp_comp.push_back(make_pair(*ri,bc)) ;
      }
    }
  }

  void compChompVisitor::process_rcm(rulecomp_map& rcm) {
    rulecomp_map::iterator ri ;
    for(ri=rcm.begin();ri!=rcm.end();++ri) {
      chomp_compiler* chc =
        dynamic_cast<chomp_compiler*>(&(*(ri->second))) ;
      if(chc != 0) {
        schedule_chomp(*chc) ;
        compile_chomp(*chc) ;
      }
    }
  }

  void compChompVisitor::visit(loop_compiler& lc) {
    process_rcm(lc.rule_compiler_map) ;
  }
  
  void compChompVisitor::visit(dag_compiler& dc) {
    process_rcm(dc.rule_compiler_map) ;
  }
  
  void compChompVisitor::visit(conditional_compiler& cc) {
    process_rcm(cc.rule_compiler_map) ;
  }

} // end of namespace


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

    // find any chains that contains v and then
    // remove those chains from the list and
    // patch them back to the search list
    void notify_existchains(list<chomp_chain>& result,
                            list<variable>& chomp_vars_order,
                            variableSet& valid_chomp_vars,
                            const variable& v) {
      vector<list<chomp_chain>::iterator> remove ;
      variableSet patch_back ;
      for(list<chomp_chain>::iterator li=result.begin() ;
          li!=result.end();++li) {
        digraph& gr = li->first ;
        variableSet internal_vars =
          extract_vars(gr.get_source_vertices() & gr.get_target_vertices()) ;
        if(internal_vars.inSet(v)) {
          remove.push_back(li) ;
          patch_back += li->second ;
        }
      }
      for(vector<list<chomp_chain>::iterator>::size_type i=0;
          i!=remove.size();++i)
        result.erase(remove[i]) ;

      valid_chomp_vars += patch_back ;
      for(variableSet::const_iterator vi=patch_back.begin();
          vi!=patch_back.end();++vi)
        chomp_vars_order.push_back(*vi) ;
    }

    // get all chomping candidates in sources
    // and targets of a graph
    inline variableSet get_chompCandInST(const digraph& gr,
                                         const variableSet& chompCand) {
      digraph::vertexSet sources = gr.get_source_vertices()
        - gr.get_target_vertices() ;
      digraph::vertexSet targets = gr.get_target_vertices()
        - gr.get_source_vertices() ;
      variableSet stv = extract_vars(sources+targets) ;

      return variableSet(stv & chompCand) ;
    }

    // test for self-cycle and non-chompable
    // internal variables in a graph
    // test_gr is a subgraph of gr
    inline bool has_problem(const digraph& gr,const digraph& test_gr,
                            const variableSet& validChompVars) {
      digraph::vertexSet sources = test_gr.get_source_vertices()
        - test_gr.get_target_vertices() ;
      digraph::vertexSet targets = test_gr.get_target_vertices()
        - test_gr.get_source_vertices() ;
      
      digraph::vertexSet internal = test_gr.get_source_vertices()
        & test_gr.get_target_vertices() ;
      
      FATAL(internal !=
            (test_gr.get_all_vertices() - sources - targets)) ;
      
      variableSet internal_vars = extract_vars(internal) ;
      variableSet problem_vars =
        variableSet(internal_vars - validChompVars) ;
      
      return (has_path(gr,targets,sources) || (problem_vars != EMPTY)) ;
    }

#define CHOMP_OPT
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
      variableSet valid_chomp_vars = chomp_vars ;
      list<chomp_chain> temp_result ;
      
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

              if(!has_problem(gr,test_gr,valid_chomp_vars)){
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

        // put it into the temp_result
        digraph sub_gr = gr.subgraph(all_vertices) ;
        // if the chain just contains one chompable variable
        // we need to test for self-cycle, if has, it is
        // then discarded
        if(sub_chomp_vars.size() == 1) {
          if(has_problem(gr,sub_gr,valid_chomp_vars)) {
            variable v = *(sub_chomp_vars.begin()) ;
            valid_chomp_vars -= v ;
            chomp_vars_order.remove(v) ;
            // we need to notify existed chains of this change
            // this variable might be silently included
            // into other chains because it may be target of
            // other chompable variables
            notify_existchains(temp_result,chomp_vars_order,
                               valid_chomp_vars,v) ;
            continue ;
          }
        }
        // we need to remove any chompable variables
        // that are the sources or targets of this
        // sub_gr. They are the cause of the self-cycles.
        // They cannot be merged into other chains because
        // they are now visible outside of this chain
        variableSet intersection = get_chompCandInST(sub_gr,chomp_vars) ;
        if(intersection != EMPTY) {
          for(variableSet::const_iterator rmvi=intersection.begin();
              rmvi!=intersection.end();++rmvi)
            chomp_vars_order.remove(*rmvi) ;
          valid_chomp_vars -= intersection ;
        }

        valid_chomp_vars -= sub_chomp_vars ;
#ifdef CHOMP_OPT
        temp_result.push_back(make_pair(sub_gr,sub_chomp_vars)) ;
#else
        result.push_back(make_pair(sub_gr,sub_chomp_vars)) ;
#endif
      } // end-of-while(!chomp_vars_order.empty())
#ifdef CHOMP_OPT      
      // here we do an optimization phase
      // we try to merge the formed chains
      // and hope we can get back some of the chomping
      // variables that are excluded in the above merging
      list<chomp_chain> temp_result2 ;
      while(!temp_result.empty()) {
        chomp_chain begin = temp_result.front() ;
        digraph new_gr = begin.first ;

        temp_result.pop_front() ;

        bool merge = true ;
        while(merge) {
          merge = false ;
          vector<list<chomp_chain>::iterator> remove ;
          for(list<chomp_chain>::iterator li=temp_result.begin();
              li!=temp_result.end();++li) {
            digraph this_gr = li->first ;
            variableSet shared =
              extract_vars(new_gr.get_all_vertices() &
                           this_gr.get_all_vertices()) ;
            if( (shared & chomp_vars) != EMPTY) {
              digraph test_gr = new_gr ;
              test_gr += this_gr ;

              if(!has_problem(gr,test_gr,chomp_vars)) {
                new_gr = test_gr ;
                remove.push_back(li) ;
                merge = true ;
              }
            }
          } // end for
          // remove those merged chains from the list
          for(vector<list<variable>::iterator>::size_type i=0;
              i!=remove.size();++i)
            temp_result.erase(remove[i]) ;
        } // end while(merge)

        variableSet chomp_vars_in_chain =
          extract_vars(new_gr.get_source_vertices() &
                       new_gr.get_target_vertices()) ;
        temp_result2.push_back(make_pair(new_gr,chomp_vars_in_chain)) ;
      } // end of while(!temp_result.empty())
      
      for(list<chomp_chain>::iterator li=temp_result2.begin();
          li!=temp_result2.end();++li) {
        // it is important to use reference(&) here,
        // because we want to dynamically update temp_result2
        // so that in the following "gcio iter" we can have
        // correct result
        digraph& new_gr = li->first ;
        variableSet& chomp_vars_in_chain = li->second ;
        
        variableSet chompCandInThis =
          get_chompCandInST(li->first,chomp_vars) ;

        variableSet chompCandInOthers ;

        // gcio iter
        for(list<chomp_chain>::iterator li2=temp_result2.begin();
            li2!=temp_result2.end();++li2)
          if(li2 != li)
            chompCandInOthers += get_chompCandInST(li->first, chomp_vars) ;

        for(variableSet::const_iterator vi=chompCandInThis.begin();
            vi!=chompCandInThis.end();++vi) {
          if(!chompCandInOthers.inSet(*vi)) {
            ruleSet rvs_rules = extract_rules(grt[vi->ident()]) ;
            ruleSet fwd_rules = extract_rules(gr[vi->ident()]) ;
            ruleSet reachable_rules = ruleSet(rvs_rules + fwd_rules) ;
            digraph::vertexSet new_vertices =
              get_ruleSet_vertexSet(reachable_rules) ;
            
            digraph test_gr = gr.subgraph(new_vertices +
                                          new_gr.get_all_vertices()) ;
            variableSet new_chomp_vars_in_chain = chomp_vars_in_chain ;
            new_chomp_vars_in_chain += *vi ;
            if(!has_problem(gr,test_gr,chomp_vars)) {
              new_gr = test_gr ;
              chomp_vars_in_chain = new_chomp_vars_in_chain ;
            }
          }
        } 
        // finally we accept the chain
        result.push_back(*li) ;
      }// end of for(temp_result2)
#endif
    }// end of the get_chomp_chains function

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
    // variables that are sources or targets of internal rules
    variableSet intstv ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(is_internal_rule(*ri)) {
        remove += *ri ;
        intstv += ri->sources() ;
        intstv += ri->targets() ;
      }

    rules -= remove ;

    // first get all possible variables that can be chomped
    variableSet allvars, remaining_vars, apply_targets ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      allvars += ri->sources() ;
      allvars += ri->targets() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
        apply_targets += ri->targets() ;
    }
    remaining_vars = allvars ;
    // intstv cannot be chomped for sure, we remove them
    remaining_vars -= intstv ;
    // targets of apply rule cannot be chomped
    remaining_vars -= apply_targets ;

    // then we collect bad_vars...

    // first any variable that involves maps
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      bad_vars += map_sources(*ri) ;
      bad_vars += map_targets(*ri) ;
    }
    remaining_vars -= bad_vars ;

    // then we get all the variables that are not suitable for
    // chomping, i.e. they are not STORE, they are in the loop
    // rotate list, loop shared variables, and variables that
    // don't have outgoing edges, and variables that connect
    // to any internal rules or rename variables
    // or variables that are generated in more that one levels
    for(variableSet::const_iterator vi=remaining_vars.begin();
        vi!=remaining_vars.end();++vi) {
      
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
    remaining_vars -= bad_vars ;

    // fill up...
    good_vars += remaining_vars ;
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

    // collect all chomped variables
    for(list<chomp_chain>::const_iterator li=chomp_chain_list.begin();
        li!=chomp_chain_list.end();++li)
      all_chomped_vars += li->second ;

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

      // check to see if there's any rules inside them, if it is,
      // then they are removed and their sources or targets are added
      ruleSet ruleinsource = extract_rules(source_vars_vertices) ;
      ruleSet ruleintarget = extract_rules(target_vars_vertices) ;
      if(ruleinsource != EMPTY) {
        digraph grt = gr.transpose() ;
        digraph::vertexSet add ;
        for(ruleSet::const_iterator ri=ruleinsource.begin();
            ri!=ruleinsource.end();++ri)
          add += grt[ri->ident()] ;
        source_vars_vertices -= get_vertexSet(ruleinsource) ;
        source_vars_vertices += add ;
      }
      if(ruleintarget != EMPTY) {
        digraph::vertexSet add ;
        for(ruleSet::const_iterator ri=ruleintarget.begin();
            ri!=ruleintarget.end();++ri)
          add += gr[ri->ident()] ;
        target_vars_vertices -= get_vertexSet(ruleintarget) ;
        target_vars_vertices += add ;
      }
      
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
        digraph::vertexSet ast ;
        ast = all_vertices - source_vars_vertices - target_vars_vertices ;
        digraph::vertexSet diff ;
        diff = takeout_vertices - ast ;
        if(diff == EMPTY)
          diff = ast - takeout_vertices ;
        variableSet diff_vars = extract_vars(diff) ;
        ruleSet diff_rules = extract_rules(diff) ;
        if(diff_vars != EMPTY)
          cerr << "These vars should be taken out: " << diff_vars << endl ;
        if(diff_rules != EMPTY)
          cerr << "These rules should be taken out: " << diff_rules << endl ;
        cerr << "The chomped vars in this chain are: " << chomp_vars << endl ;
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
  namespace {
    // scheduling utilities for chomping graph
    chompingPrio cpf ;
    graphSchedulerVisitor cgsv(cpf) ;
  }
  
  void compChompVisitor::schedule_chomp(chomp_compiler& chc) {
    //chc.chomp_sched = orderVisitor::order_dag(chc.chomp_graph) ;
    chc.chomp_sched = cgsv.schedule(chc.chomp_graph) ;
  }
  
  void compChompVisitor::compile_chomp(chomp_compiler& chc,
                                       const rulecomp_map& rcm) {
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
        rulecomp_map::const_iterator rmi ;
        rmi = rcm.find(*ri) ;
        FATAL(rmi == rcm.end()) ;
        //rule_compilerP bc = new barrier_compiler(vars) ;
        chc.chomp_comp.push_back(make_pair(*ri,rmi->second)) ;
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
        compile_chomp(*chc,rcm) ;
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


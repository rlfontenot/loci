#include <depend_graph.h>
#include <map>
using std::map ;
#include <vector>
using std::vector ;
#include <set>
using std::set ;

#ifndef DEBUG
#define PRUNE_GRAPH
#else
// Comment this line out if we want stricter debug options
#define PRUNE_GRAPH
#endif
//#define VERBOSE

namespace Loci {

  namespace {
    rule create_rule(variable sv, variable tv, string qualifier) {
      ostringstream oss ;
      oss << "source(" << sv << ')' ;
      oss << ",target(" << tv << ')' ;
      oss << ",qualifier(" << qualifier << ')' ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
    
    rule create_rule(variableSet source, variableSet target, string qualifier) {
      ostringstream oss ;
      oss << "source(" << source << ')' ;
      oss << ",target(" << target << ')' ;
      oss << ",qualifier(" << qualifier << ")" ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }

    void print_graph_from(variableSet given, const digraph &rule_graph) {
      variableSet working = given ;
      variableSet processed ;
      ruleSet processed_rules ;
      cout << "GRAPH Traversal results:" << endl ;
      while(working != EMPTY) {
        variableSet new_vars ;
        ruleSet new_rules ;
        for(variableSet::const_iterator vi = working.begin();
            vi!=working.end();
            ++vi) {
          ruleSet var_rules = extract_rules(rule_graph.get_edges(vi->ident())) ;
          var_rules -= processed_rules ;
          new_rules += var_rules ;
          for(ruleSet::const_iterator ri = var_rules.begin();
              ri!=var_rules.end();
              ++ri) {
            new_vars += extract_vars(rule_graph.get_edges(ri->ident())) ;
          }
        }
        cout << "------------------------------------------------------------"
             << endl ;
        cout << "VARS="<<working << endl ;
        cout << "RULES=" << endl << new_rules ;
        processed_rules += new_rules ;
        processed += working ;
        new_vars -= processed ;
        working = new_vars ;
      }
      cout << "------------------------------------------------------------"
           << endl ;
    }
    
    inline void invoke_rule(rule f, digraph &gr) {
      gr.add_edges(f.sources(),f.ident()) ;
      gr.add_edges(f.ident(),f.targets()) ;
    }

    void invoke_rule_wp(rule f, digraph &gr) {
      variableSet targets = f.targets() ;
      variableSet sources = f.sources() ;
      invoke_rule(f,gr) ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();
          ++vi) {
        if(vi->get_info().priority.size() != 0) {
          variable v = *vi ;
          variable vold = v ;
          while(v.get_info().priority.size() != 0) {
            v = v.drop_priority() ;
            rule priority_rule = create_rule(vold,v,"priority") ;
            invoke_rule(priority_rule,gr) ;
            vold = v ;
          }
        }
      }
    }

    ruleSet fill_graph(variableSet start, const digraph &rule_graph,
                       digraph &gr,variableSet known) {
#ifdef VERBOSE
      cout << "fillgraph(start="<<start<<",known="<<known<<")"<<endl ;
#endif
      digraph rgt = rule_graph.transpose() ;
      variableSet working = start ;
      variableSet processed ;
      ruleSet processed_rules ;
      while(working != EMPTY) {
#ifdef VERBOSE
        cout << "fill_graph: working = " << working << endl ;
#endif
        variableSet new_vars ;
        ruleSet working_rules ;
        processed += working ;
        variableSet known_processed = processed ;
        known_processed += known ;
        for(variableSet::const_iterator vi = working.begin();
            vi!=working.end();
            ++vi) {
          ruleSet var_rules =
            extract_rules(rule_graph.get_edges(vi->ident())+
                          rgt.get_edges(vi->ident())) ;
          var_rules -= processed_rules ;
          ruleSet reject_rules ;
          for(ruleSet::const_iterator ri = var_rules.begin() ;
              ri!=var_rules.end();
              ++ri) {
            if(ri->type() == rule::INTERNAL &&
               ri->qualifier()=="iterating_rule") {
              if(vi->time() == ri->source_time()) {
                invoke_rule(*ri,gr) ;
                new_vars += extract_vars(rule_graph.get_edges(ri->ident())) ;
              } else {
                cerr << "rejecting iterating rule " << *ri << " for variable " << *vi << endl ;
                cerr << "Please report when this happens" << endl ;
                reject_rules += *ri ;
              }
            } else {
              variableSet rule_depend = ri->constraints() ;
#ifdef PRUNE_GRAPH
              // Adding this line will prune dependency graph to only those
              // rules that can execute.
              rule_depend += ri->sources() ;
#endif
              if(rule_depend == EMPTY)
                rule_depend = ri->sources() ;
                  
              rule_depend -= rule_graph.get_edges(ri->ident()) ;
              rule_depend -= known_processed ;
              variableSet time_vars ;
              for(variableSet::const_iterator vi = rule_depend.begin();
                  vi != rule_depend.end();
                  ++vi) {
                if(vi->get_info().tvar)
                  time_vars += *vi ;
              }
              rule_depend -= time_vars ;
              if(rule_depend == EMPTY) {
                invoke_rule_wp(*ri,gr) ;
                new_vars += extract_vars(rule_graph.get_edges(ri->ident())) ;
              } else {
                reject_rules += *ri ;
              }
            }
          }
          var_rules -= reject_rules ;
          processed_rules += var_rules ;
          working_rules += var_rules ;

#ifdef VERBOSE
          cout << "rules involved = " << var_rules << endl ;
#endif
          
          time_ident vtime = vi->time() ;
          if(vtime != time_ident()) {
            variable stationary_var(*vi,time_ident()) ;
            ruleSet promote_rules =
              extract_rules(rule_graph.get_edges(stationary_var.ident())+
                            rgt.get_edges(stationary_var.ident())) ;

            for(ruleSet::const_iterator ri = promote_rules.begin() ;
                ri!=promote_rules.end();
                ++ri) {
              if(ri->type() != rule::TIME_SPECIFIC &&
                 !(ri->type() == rule::INTERNAL &&
                   ri->qualifier()=="iterating_rule")) {
                rule pr(*ri,vtime) ;
                if(!processed_rules.inSet(pr)) {
                  variableSet rule_depend = pr.get_info().constraints() ;
                  // Hack
                  rule_depend += pr.sources() ;
                  if(rule_depend == EMPTY) 
                    rule_depend = pr.sources() ;
                  rule_depend -= known_processed ;
                  if(rule_depend == EMPTY) {
#ifdef VERBOSE
                    cout << "promote rule = " << pr << endl ;
#endif
                    processed_rules += pr ;
                    working_rules += pr ;
                    invoke_rule(pr,gr) ;
                    new_vars += pr.targets() ;
                  }
                }
              }
            }
          }
        }
        new_vars -= processed ;
        working = new_vars ;
      }
#ifdef VERBOSE
      cout << "return from fill_graph" << endl ;
#endif
      return ruleSet(EMPTY) ;
    }

    struct graph_info ;
    
    struct iteration_info {
      bool active ;
      digraph iteration_graph ;
      time_ident iteration_time ;
      rule iteration_rule ;
      ruleSet build,advance,collapse ;
      variableSet known_vars ;
      void build_graph(const digraph &rule_graph,variableSet input_vars) ;
    } ;

    struct graph_info {
    map<time_ident,iteration_info> iteration_rules ;
    map<rule,time_ident> iteration_time_ident ;
    } ;

    
    void iteration_info::build_graph(const digraph &rule_graph,
                                     variableSet input_vars) {
      active = false ;
      
    
      ruleSet rules_that_pass ;
      for(ruleSet::const_iterator ri = build.begin();ri!=build.end();++ri) {
        variableSet bs = ri->sources() ;
        bs -= input_vars ;
        if(bs == EMPTY)
          rules_that_pass += *ri ;
        else {
          //          cerr << "build rule failed because of " << bs << endl ;
        }
      }
      if(rules_that_pass == EMPTY)
        return ;   // No build rules can execute, we are finished

      
      rules_that_pass = EMPTY ;
      for(ruleSet::const_iterator ri = collapse.begin();
          ri!=collapse.end();
          ++ri) {
        variableSet cs = ri->sources() ;
        variableSet csp ;
        for(variableSet::const_iterator vi=cs.begin();vi!=cs.end();++vi) {
          if(vi->time() != iteration_time)
            csp += *vi ;
        }
        csp -= input_vars ;
        if(csp == EMPTY)
          rules_that_pass += *ri ;
        else {
          //          cerr << "collapse rule failed because of " << csp << endl ;
        }
      }
      if(rules_that_pass == EMPTY)
        return ;    // No collapse rules can execute, we are finished


      // This iteration can be scheduled so fill out the graph
      active = true ;

      ruleSet all_iteration_rules ;
      all_iteration_rules += build ;
      all_iteration_rules += advance  ;
      all_iteration_rules += collapse ;
      variableSet build_vars ;
      known_vars = input_vars ;

#ifdef VERBOSE
      cout << "build_rules = " << build << endl ;
      cout << "advance_rules = " << advance << endl ;
      cout << "collapse_rules = " << collapse << endl ;
#endif

      map<variable,intervalSet> build_offsets ;
      
      for(ruleSet::const_iterator ri = build.begin();ri!=build.end();++ri) {
        variableSet btarget = ri->targets() ;
        for(variableSet::const_iterator vi = btarget.begin();
            vi!=btarget.end();
            ++vi) {
          if(!vi->assign) {
            cerr << "incorrect specification of build rule " << *ri
                 << endl
                 << "A correct build rule should have an assign in the target"
                 << " variable.  For example, instead of " << *vi
                 << " use " << vi->name << "{"<< vi->time_id<<"=0}"
                 << endl ;
          }
          warn(!vi->assign) ;

          build_vars += *vi ;
          variable vt = vi->drop_assign() ;
          variable vbase = vt.new_offset(0) ;
          build_offsets[vbase] += vi->offset ;
          rule gv =  create_rule(*vi,vt,"generalize") ;
          invoke_rule(gv,iteration_graph) ;
          if(vi->offset == 0)
            build_vars += vt ;
        }
      }

      variableSet looping_input ;
      variableSet looping_output ;
      map<variable,intervalSet>::iterator mvi ;
      for(mvi = build_offsets.begin();mvi!=build_offsets.end();++mvi) {
        variable bv = mvi->first ;
        for(intervalSet::const_iterator ii=mvi->second.begin();
            ii!=mvi->second.end();
            ++ii) {
          looping_output += bv.new_offset(*ii) ;
        }
        looping_input += bv.new_offset(mvi->second.Max()+1) ;
      }
      variable ov("OUTPUT") ;
      looping_input += variable(ov,iteration_time) ;
      looping_output+= variable(ov,iteration_time) ;
      
      looping_output += variable(iteration_time) ;
      ostringstream oss ;
      oss << "source("<< looping_input
          << "),target(" << looping_output
          << "),qualifier(looping)" ;
      rule floop(oss.str()) ;
      known_vars += looping_output ;
      
      for(variableSet::const_iterator vi =input_vars.begin();
          vi!=input_vars.end();
          ++vi) {
        variable nv(*vi,iteration_time) ;
        known_vars += nv ;
      }
      for(ruleSet::const_iterator ri = all_iteration_rules.begin();
          ri != all_iteration_rules.end();
          ++ri) {
        invoke_rule_wp(*ri,iteration_graph) ;
      }
      build_vars += variable(iteration_time) ;
      //      known_vars += variable(iteration_time) ;
      
      fill_graph(build_vars,rule_graph,iteration_graph,known_vars) ;
      known_vars += extract_vars(iteration_graph.get_target_vertices()) ;

      invoke_rule(floop,iteration_graph) ;
    }

    variableSet convert_stationary(const variableSet &v) {
      variableSet result ;
      variableSet::const_iterator vi ;
      time_ident stationary_time ;
      for(vi=v.begin();vi!=v.end();++vi) {
        if(vi->time() != stationary_time)
          result += variable(*vi,time_ident()) ;
        else
          result += *vi ;
      }
      return result ;
    }

    void add_rename_dependencies(digraph &gr) {
      variableSet all_vars = extract_vars(gr.get_all_vertices()) ;
      ruleSet     all_rules = extract_rules(gr.get_all_vertices()) ;
      
      // extract the qualified rules, these are rules that are
      // automatically generated by the system.  Since these rules
      // cannot provide type information directly, they are
      // singled out so that they can be handled as a separate case.
      ruleSet qualified_rules,rename_rules ;
      for(ruleSet::const_iterator ri=all_rules.begin();
          ri!=all_rules.end();
          ++ri) {
        if(ri->type() == rule::INTERNAL)
          qualified_rules += *ri ;
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
          if(vmsi->assign.size() != 0) 
            rename_rules += *ri ;
      }
      // We need the transpose of the graph in order to find the rules that
      // generate a particular variable

      ruleSet check_rules = all_rules ;
      check_rules -= qualified_rules ;
      for(ruleSet::const_iterator ri=check_rules.begin();
          ri != check_rules.end();
          ++ri) {
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
          if(vmsi->assign.size() != 0) 
            for(int i=0;i<vmsi->assign.size();++i) {
              variable orig_name = vmsi->assign[i].second ;
              //              digraph grt = gr.transpose() ;
              ruleSet depend_rules = extract_rules(gr[orig_name.ident()]) ;
              depend_rules -= rename_rules ;
              // We ignore renaming dependencies in collapse rules, however
              // in reality we may need to follow the dependencies back to
              // the build rules.  However, since we know that the collapse
              // must follow the build, there is no need to add the rename
              // dependencies to collapse rules.
              if(ri->type() != rule::COLLAPSE)
                gr.add_edges(depend_rules,ri->ident()) ;
#ifdef VERBOSE
              cout << "adding edges from " <<depend_rules<<
                " to " <<*ri << endl ;
#endif
            }
      }
      
    }
  
  }
  
  dependency_graph::dependency_graph(rule_db &rdb, variableSet given,
                                             variableSet target) {


#ifdef VERBOSE
    cout << "in dependency_graph:" << endl ;
    cout << "given = " << given << endl ;
    cout << "target = " << target << endl ;
#endif
    
    ruleSet all_rules = rdb.all_rules() ;

    
    ruleSet::const_iterator ri ;

    graph_info gi ;
    
    ruleSet iteration_set ;
    
    ruleSet working_rules ;
    for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
      const rule::rule_type type = ri->type() ;
      if(type==rule::BUILD) {
        gi.iteration_rules[ri->target_time()].build += *ri ;
      } else if(type == rule::COLLAPSE) {
        gi.iteration_rules[ri->source_time()].collapse += *ri ;
      } else if(!ri->time_advance) {
        working_rules += *ri ;
      }
      if(ri->time_advance) {
        if(type != rule::COLLAPSE)
          gi.iteration_rules[ri->target_time()].advance += *ri ;
      }
    }


    map<time_ident,iteration_info>::iterator mi ;
    
    for(mi=gi.iteration_rules.begin();mi!=gi.iteration_rules.end();++mi) {
      mi->second.iteration_time = mi->first ;
      ruleSet build = mi->second.build ;
      ruleSet collapse = mi->second.collapse;
      if(build == EMPTY)  {
        cerr << "malformed iteration: " << mi->first <<
          ", no build rules" << endl ;
        abort() ;
      }
      if(collapse == EMPTY)  {
        cerr << "malformed iteration: " << mi->first <<
          ", no collapse rules" << endl ;
        abort() ;
      }
      ruleSet iteration_vars ;
      variableSet source,target ;
      for(ri=build.begin();ri!=build.end();++ri) {
        const rule::rule_type type = ri->type() ;
        source += ri->sources() ;
      }
      for(ri=collapse.begin();ri!=collapse.end();++ri) {
        const rule::rule_type type = ri->type() ;
        target += ri->targets() ;
      }
      rule i_rule = create_rule(source,target,"iterating_rule") ;
      gi.iteration_time_ident[i_rule] = mi->first ;
      mi->second.iteration_rule  = i_rule ;
      iteration_set += i_rule ;
      if(i_rule.get_info().time_advance) {
        gi.iteration_rules[i_rule.target_time()].advance += i_rule ;
      } else {
        working_rules += i_rule ;
      }
    }

    ruleSet dump_iterations ;
    // Check iterations for validity
    for(mi=gi.iteration_rules.begin();mi!=gi.iteration_rules.end();++mi) {
      variableSet build_targets,advance_targets ;
      ruleSet build = mi->second.build ;
      ruleSet advance = mi->second.advance;
      if(advance == EMPTY) {
        cerr << "malformed iteration: " << mi->first <<
          ", no advance rules" << endl ;
        abort() ;
      }
    }
    
    digraph rule_graph ;

    working_rules -= dump_iterations ;

    
    for(ruleSet::const_iterator ri = working_rules.begin();
        ri != working_rules.end();
        ++ri) {
      invoke_rule_wp(*ri,rule_graph) ;
    }
  
    // Fill graph with rules that will compute target.
    fill_graph(given,rule_graph,gr,given) ;


    ruleSet scheduled_iteration_rules =
      extract_rules(gr.get_all_vertices()&iteration_set) ;
    ruleSet visited_iteration_rules ;
    variableSet baseline_vars = extract_vars(gr.get_all_vertices()) ;
    baseline_vars += given ;
    while(scheduled_iteration_rules != EMPTY) {
      ruleSet new_iteration_rules ;
      for(ruleSet::const_iterator ri = scheduled_iteration_rules.begin();
          ri != scheduled_iteration_rules.end();
          ++ri) {
        map<rule,time_ident>::const_iterator tp =
          gi.iteration_time_ident.find(*ri) ;
        if(tp == gi.iteration_time_ident.end()) {
          cerr << "something is wrong"<<endl;
          abort() ;
        }
        map<time_ident,iteration_info>::iterator ip =
          gi.iteration_rules.find(tp->second) ;
        if(ip == gi.iteration_rules.end()) {
          cerr << " rule doesn't exist, confused!"<< endl ;
          abort() ;
        }
        time_ident iteration_time = tp->second ;
        time_ident parent_time = iteration_time.parent() ;

        variableSet known_vars ;
        if(parent_time == time_ident()) {
          known_vars = baseline_vars ;
        } else {
          map<time_ident,iteration_info>::iterator ipp =
            gi.iteration_rules.find(parent_time) ;
          if(ipp == gi.iteration_rules.end()) {
            cerr << "parent iteration rule doesn't exit for time level"
                 << parent_time << endl ;
            known_vars = baseline_vars ;
          } else {
            known_vars = ipp->second.known_vars ;
          }
        }
        ip->second.build_graph(rule_graph,known_vars) ;
        new_iteration_rules +=
          extract_rules(ip->second.iteration_graph.get_all_vertices() &
                        iteration_set) ;
      }
      visited_iteration_rules += scheduled_iteration_rules ;
#ifdef VERBOSE
      cout << "new_iteration_rules = " << new_iteration_rules << endl ;
      cout << "scheduled_iteration_rules = " << scheduled_iteration_rules << endl ;
#endif
      
      warn((new_iteration_rules & visited_iteration_rules) != EMPTY) ;

#ifdef DEBUG
      ruleSet tmp = new_iteration_rules ;
      tmp &= visited_iteration_rules ;
      if(tmp != EMPTY) {
        cerr << "new_iteration_rules = " << new_iteration_rules << endl ;
        cerr << "new_iteration_rules & visited_iteration_rules = " <<
          tmp << endl ;
      }
#endif
        
      
      scheduled_iteration_rules = new_iteration_rules ;
      scheduled_iteration_rules -= visited_iteration_rules ;
      
    }
    map<time_ident, iteration_info>::iterator ip ;
    // Assemble All of the Iteration Graphs
    for(ip = gi.iteration_rules.begin();ip!=gi.iteration_rules.end();++ip) {
      if(ip->second.active) {
        gr.add_graph(ip->second.iteration_graph) ;
      }
    }
    // Add Promotion Rules!
    for(ip = gi.iteration_rules.begin();ip!=gi.iteration_rules.end();++ip) {
      if(ip->second.active) {
        digraph grt = gr.transpose() ;
        digraph &ig = ip->second.iteration_graph ;
        digraph igt = ig.transpose() ;
        variableSet ivars = extract_vars(ig.get_all_vertices()) ;
#ifdef VERBOSE
        cout << "ivars = " << ivars << endl ;
#endif

        for(variableSet::const_iterator vi=ivars.begin();vi!=ivars.end();++vi) {
          if(grt.get_edges(vi->ident()) == EMPTY) {
            const variable::info &kvi = vi->get_info() ;
            if(!kvi.assign && !kvi.tvar && kvi.offset == 0
               && kvi.time() != time_ident()) {
              //               && kvi.time() == ip->first) {
              time_ident parent = kvi.time().parent() ;
              variable vs = *vi ;
              while(parent != time_ident()) {
                variable vp(*vi,parent) ;
                invoke_rule(create_rule(vp,vs,"promote"),gr) ;
#ifdef VERBOSE
                cout << "adding promote for iteration " << ip->first << endl ;
                cout << vp << " to " << vs << endl ;
#endif
                if(grt.get_edges(vp.ident()) != EMPTY) {
                  break ;
                }
                parent = parent.parent() ;
                vs = vp ;
              }
              if(parent == time_ident()) {
                variable vp(*vi,parent) ;
                invoke_rule(create_rule(vp,vs,"promote"),gr) ;
#ifdef VERBOSE
                cout << "adding promote for iteration " << ip->first << endl ;
                cout << vp << " to " << vs << endl ;
#endif
              }
                
            }
          }
        }

      }
    }

    gr.remove_vertices(visited_iteration_rules) ;
    create_looping_rules() ;

#ifdef VERBOSE
    print_graph_from(given,gr) ;
#endif

#ifdef PRUNE_GRAPH
    clean_graph(given,target) ;
#endif

    add_rename_dependencies(gr) ;
  }

  void dependency_graph::create_looping_rules() {
    //#define OLDE
#ifdef OLDE
    digraph::vertexSet all_vertices = gr.get_all_vertices() ;
    variableSet all_vars = extract_vars(all_vertices) ;
    variableSet::const_iterator vi ;
    map<time_ident,variableSet> tvarmap ;
    for(vi=all_vars.begin();vi!=all_vars.end();++vi) 
      if(vi->get_info().offset == 1 || vi->get_info().name==string("OUTPUT")) 
        tvarmap[vi->time()] += *vi ;

    map<time_ident,variableSet>::const_iterator ii ;
    for(ii=tvarmap.begin();ii!=tvarmap.end();++ii) {
      variableSet source = ii->second ;
      variableSet target ;
      target += variable(ii->first) ;
      bool advance_vars = false ;
      for(vi=source.begin();vi!=source.end();++vi)
        if(vi->get_info().offset == 1)
          advance_vars = true ;

      if(advance_vars) {
        for(vi=source.begin();vi!=source.end();++vi)
          target += vi->new_offset(0) ;
        ostringstream oss ;
        oss << "source(" << source
            << "),target(" << target
            << "),qualifier(looping)" ;
        rule f(oss.str()) ;
        gr.add_edges(source,f.ident()) ;
        gr.add_edges(f.ident(),target) ;
      }
    }
#endif
  }

  void dependency_graph::clean_graph(variableSet given, variableSet target) {
    // Remove unnecessary vertices from graph.
    int virtual_vertex = gr.max_vertex() + 1 ;
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    variableSet allvars = extract_vars(allvertices) ;
    variableSet::const_iterator vi ;
    //  target += variable(expression::create("OUTPUT")) ;

    gr.add_edges(virtual_vertex, given) ;
    gr.add_edges(target,virtual_vertex) ;

    const vector<digraph::vertexSet> components =
      component_sort(gr).get_components() ;

    digraph::vertexSet subset = EMPTY ;
    
    for(int i=0;i<components.size();++i) 
      if(components[i].inSet(virtual_vertex)) {
        subset = components[i] ;
        break ;
      }
    subset -= virtual_vertex ;
    ruleSet rules = extract_rules(subset) ;
    ruleSet::const_iterator fi ;
    for(fi=rules.begin();fi!=rules.end();++fi)
      subset += fi->targets() ;

    // Check for looping rules here, don't clean looping rule if it is
    // in the subset.
    digraph grt = gr.transpose() ;
    digraph::vertexSet  cleanout ;
    for(fi=rules.begin();fi!=rules.end();++fi)
      if((subset & fi->sources()) != fi->sources()) {
        cleanout += fi->ident() ;
      }


    subset -= cleanout ;
    
    WARN(subset == EMPTY) ;
#ifdef VERBOSE
    cout << "cleaning out rules: " << endl ;
    cout << extract_rules(gr.get_source_vertices()-subset) ;
#endif
    gr = gr.subgraph(subset) ;
  }


}// End namespace Loci

#include <depend_graph.h>

#include <map>
#include <vector>


using std::map ;
using std::vector ;

//#define VERBOSE

namespace Loci {

  dependency_graph::dependency_graph(rule_db &rdb, variableSet given,
                                     variableSet target) {
    compose_graph(rdb,given) ;
    remove_incorrect_time_promotions() ;
    create_looping_rules() ;
    clean_graph(given,target) ;
  }

  variableSet dependency_graph::invoke_rule(rule f) {
#ifdef VERBOSE
    cout << "invoking rule: " << f << endl ;
#endif
    gr.add_edges(f.sources(),f.ident()) ;
    variableSet targets = f.targets() ;
    gr.add_edges(f.ident(),targets) ;
    variableSet new_targets = targets ;
    for(variableSet::const_iterator vii=targets.begin();
        vii != targets.end(); ++vii) {
      variable v = *vii ;
      if(v.get_info().priority.size() != 0) {
        variable vold = v ;
        while(v.get_info().priority.size() != 0) {
          v = v.drop_priority() ;
          ostringstream oss ;
          oss << "source(" << vold << "),target("<<v<<"),qualifier(priority)" ;
          rule fpri(oss.str()) ;
          gr.add_edge(vold.ident(),fpri.ident()) ;
          gr.add_edge(fpri.ident(),v.ident()) ;
          vold  = v ;
        }
        new_targets -= *vii ;
        new_targets += v ;
      }
    }
    return new_targets ;
  }



  void dependency_graph::promote_variable(variable v1, variable v2)
    {
      if(v1.get_info().name != string("OUTPUT")) {
        ostringstream oss ;
        oss << "source(" << v1 << "),target("<<v2<<"),qualifier(promote)" ;
        rule f(oss.str()) ;
        gr.add_edge(v1.ident(),f.ident()) ;
        gr.add_edge(f.ident(),v2.ident()) ;
      }
    }

  void dependency_graph::generalize_variable(variable v1, variable v2)
    {
      ostringstream oss ;
      oss << "source(" << v1 << "),target("<<v2<<"),qualifier(generalize)" ;
      rule f(oss.str()) ;
      gr.add_edge(v1.ident(),f.ident()) ;
      gr.add_edge(f.ident(),v2.ident()) ;
    }


  void dependency_graph::compose_graph(rule_db &rdb, variableSet given)
    {
      ruleSet promote_rule ;
      ruleSet::const_iterator fi ;
      //  digraph gr ;
      const ruleSet &known_rules = rdb.all_rules() ;
      for(fi=known_rules.begin();fi!=known_rules.end();++fi)
        if(fi->type() == rule::GENERIC)
          promote_rule += *fi ;

      variableSet known = given ;
      variableSet working = given ;
      variableSet processed = given ;
      variableSet newtime ;
      variableSet newvars ;
      intervalSet time_levels ;
      ruleSet rset ; // set of all invoked rules
      time_levels += 0 ;
      for(variableSet::const_iterator i=newvars.begin();i!=newvars.end();++i) 
        if(!time_levels.inSet(i->time().ident())) 
          newtime += *i ;
      working -= newtime ;
      while(working != EMPTY) {
#ifdef VERBOSE    
        cout << "WORKING = " << working << endl ;
#endif            
        variableSet::const_iterator i,j ;
        for(i=working.begin();i!=working.end();++i) {
          const variable::info &vi = (*i).get_info();
            
          if(vi.assign) {
            variable v2 = i->drop_assign() ; ;

#ifdef VERBOSE
            cout << "converting variable: " << *i << " to " << v2 << endl ;
#endif
            //        gr.add_edge((*i).ident(),v2.ident()) ;
            generalize_variable(*i,v2) ;

            if(!processed.inSet(v2))
              newvars += v2 ;
          } else if(vi.priority.begin() != vi.priority.end()) {

            variable v2 = i->drop_priority() ;

#ifdef VERBOSE
            cout << "converting  variable: " << *i << " to " << v2 << endl ;
#endif
            gr.add_edge((*i).ident(),v2.ident()) ;

            if(!processed.inSet(v2))
              newvars += v2 ;
          }

          const ruleSet &fncs = rdb.rules_by_source(*i) ;
          for(fi=fncs.begin();fi!=fncs.end();++fi) {
            if(!rset.inSet(*fi) && (fi->sources()-known)==EMPTY) {
              rset += *fi ;
              variableSet targets = invoke_rule(*fi) ;
              newvars += targets - processed ;
              if(fi->type() == rule::BUILD) {
                variable tvar(fi->target_time()) ;
                if(!processed.inSet(tvar)) {
#ifdef VERBOSE
                  cout << "tvar = " << tvar << endl ;
#endif
                  // add time variable dependency to graph
                  //              gr.add_edge((*fi).ident(),tvar.ident()) ;
                  newvars += tvar ;
                }
              }
            }
          }
        }
        variableSet pset ;
        for(j=newvars.begin();j!=newvars.end();++j) 
          for(intervalSet::const_iterator ti=time_levels.begin();
              ti!=time_levels.end();++ti) {
            if(*ti == 0)
              continue ;
            const variable::info &kvi = (*j).get_info() ;
            time_ident parent = time_ident(*ti).parent() ;
            if(!kvi.assign && kvi.offset == 0 && kvi.time_id == parent){
              variable v2(*j,time_ident(*ti)) ;
              if(!known.inSet(v2)) {
#ifdef VERBOSE
                cout << "promoting " << *j << " to " << v2 << endl ;
#endif
            
                //            gr.add_edge((*j).ident(),v2.ident()) ;
                promote_variable(*j,v2) ;
                if(!processed.inSet(v2))
                  pset += v2 ;
              }
            }
          }
        newvars += pset ;

        processed += working ;

        for(i=newvars.begin();i!=newvars.end();++i) 
          if(!time_levels.inSet(i->time().ident())) 
            newtime += *i ;
        newvars -= newtime ;

        if(newvars == EMPTY && newtime != EMPTY) {
          i = newtime.begin() ;
          time_ident time_id = i->time() ;
          if(!time_levels.inSet(time_id.ident())) {
#ifdef VERBOSE
            cout << "new time level found with variable " << *i << endl ;
#endif
            time_ident parent = time_id.parent() ;
            variableSet pset ;
            for(j=known.begin();j!=known.end();++j) {
              const variable::info &kvi = (*j).get_info() ;
              if(!kvi.assign && kvi.offset == 0 && kvi.time_id == parent){
                variable v2(*j,time_id) ;
                if(!known.inSet(v2)) {
#ifdef VERBOSE
                  cout << "promoting " << *j << " to " << v2 << endl ;
#endif
                  //              gr.add_edge((*j).ident(),v2.ident()) ;
                  promote_variable(*j,v2) ;
                  pset += v2 ;
                }
              }
            }
            known += pset ;
            //                newvars += pset ;
            for(fi=promote_rule.begin();fi!=promote_rule.end();++fi) {
              rdb.add_rule(rule(*fi,time_id)) ;
            }
            time_levels += time_id.ident() ;
          }
          newvars += *i ;
          newtime -= newvars ;
        }
    
        working = newvars ;
        known += newvars ;
        newvars = EMPTY ;
      }

      for(variableSet::const_iterator vi=processed.begin();
          vi != processed.end();++vi) {
        const ruleSet &fncs = rdb.rules_by_target(*vi) ;
        ruleSet::const_iterator fi ;
        for(fi=fncs.begin();fi!=fncs.end();++fi) 
          if(!rset.inSet(*fi) && (fi->sources()-known)==EMPTY) {
            rset += *fi ;
            variableSet new_targets = invoke_rule(*fi) ;
            WARN((new_targets - processed) != EMPTY) ;
#ifdef DEBUG
            if((new_targets-processed) != EMPTY) {
              cerr << "offending rule is " << *fi << endl ;
              cerr << "new_targets = " << new_targets << endl;
              cerr << "processed = " << processed << endl ;
            }
#endif
          }
      }

#ifdef VERBOSE
      cout << "processed_variables = " << processed << endl ;
#endif
    }

  void dependency_graph::remove_incorrect_time_promotions() {
    digraph grt = gr.transpose() ;
    digraph::nodeSet all_nodes = gr.get_all_nodes() ;
    ruleSet all_rules = extract_rules(all_nodes) ;
    variableSet all_vars  = extract_vars(all_nodes) ;
    ruleSet promote_rules,unit_rules ;
    ruleSet::const_iterator fi ;
    for(fi=all_rules.begin();fi!=all_rules.end();++fi) {
      if(fi->type() == rule::INTERNAL && fi->qualifier() == "promote")
        promote_rules += *fi ;
      if(fi->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) 
        unit_rules += *fi ;
    }
  
    variableSet::const_iterator vi ;
    ruleSet remove_rules ;
    for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
      ruleSet var_rules = extract_rules(grt[(*vi).ident()]) ;
      ruleSet var_promote = var_rules ;
      var_promote &= promote_rules ;
      if(var_promote != EMPTY && var_rules != var_promote)
        remove_rules += var_promote ;
    }

#ifdef VERBOSE
    cout << "removing rules: " << endl << remove_rules ;
#endif
  
    // looping unit rules to themselves, this is so that
    // all reduction rules remain together.
    for(fi = unit_rules.begin();fi!=unit_rules.end();++fi)
      gr.add_edges(fi->targets(),(*fi).ident()) ;

    gr.remove_nodes(remove_rules) ;
  }

  void dependency_graph::create_looping_rules() {
    digraph::nodeSet all_nodes = gr.get_all_nodes() ;
    variableSet all_vars = extract_vars(all_nodes) ;
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
  }

  void dependency_graph::clean_graph(variableSet given, variableSet target) {
    // Remove unnecessary nodes from graph.
    int virtual_node = gr.max_node() + 1 ;
    digraph::nodeSet allnodes = gr.get_all_nodes() ;
    variableSet allvars = extract_vars(allnodes) ;
    variableSet::const_iterator vi ;
    //  target += variable(expression::create("OUTPUT")) ;

    gr.add_edges(virtual_node, given) ;
    gr.add_edges(target,virtual_node) ;

    const vector<digraph::nodeSet> components =
      component_sort(gr).get_components() ;

    digraph::nodeSet subset = EMPTY ;
    
    for(int i=0;i<components.size();++i) 
      if(components[i].inSet(virtual_node)) {
        subset = components[i] ;
        break ;
      }
    subset -= virtual_node ;
    ruleSet rules = extract_rules(subset) ;
    ruleSet::const_iterator fi ;
    for(fi=rules.begin();fi!=rules.end();++fi)
      subset += fi->targets() ;
  
    WARN(subset == EMPTY) ;
#ifdef VERBOSE
    cout << "cleaning out rules: " << endl ;
    cout << extract_rules(gr.get_source_nodes()-subset) ;
#endif
    gr = gr.subgraph(subset) ;

     
  }

}

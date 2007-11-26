#include <param_rule.h>

using std::string ;
using std::ostringstream ;
#include <distribute.h>

//#define VERBOSE

namespace Loci {
  using std::pair ;
  using std::string ;
  using std::map ;
  using std::vector ;
  using std::cout ;
  
  typedef pair<string,int> param_rule_key ;
  typedef map<param_rule_key, ruleSet> param_rule_db ;

  string var2key(variable v) {
    if(v.get_info().namespac.size() == 0) {
      return v.get_info().name ;
    }
    vector<string> ns = v.get_info().namespac ;
    string n ;
    for(size_t i = 0;i<ns.size();++i)
      n+= ns[i] + '@' ;
    n += v.get_info().name ;
    return n ;
  }
  
  // Return set of parametic variables that were accessed by the given ruleSet
  variableSet scanRulesForParametrics(ruleSet rs) {
    ruleSet::const_iterator ri ;
    variableSet pvars ;
    for(ri=rs.begin();ri!=rs.end();++ri) {
      variableSet sources = ri->sources() ;
      variableSet::const_iterator vi ;
      for(vi=sources.begin();vi!=sources.end();++vi)
        if(variable(*vi).get_arg_list().size() > 0) {
          variable v = *vi ;
          while(v.time() != time_ident())
            v=v.parent() ;
          while(v.get_info().priority.size() != 0)
            v = v.drop_priority() ;

          pvars += v ;
        }
    }
    return pvars ;
  }

  variable applySubstitution(map<string,variable> &transform_map, variable v) {
    variable::info vinfo = v.get_info() ;
    // Substitute name with underscore separation
    string name = vinfo.name ;
    vector<string> vname ;
    vname.push_back(string()) ;
    for(size_t i=0;i<name.size();++i) {
      if(name[i] == '_')
        vname.push_back(string()) ;
      else
        vname.back() += name[i] ;
    }

    vector<variable> subst ;
    for(size_t i=0;i<vname.size();++i) {
      map<string,variable>::const_iterator ii = transform_map.find(vname[i]) ;
      if(ii != transform_map.end()) {
        subst.push_back(ii->second) ;
        vname[i] = ii->second.get_info().name ;
      }
    }
    string sname =vname[0] ;
    for(size_t i=1;i<vname.size();++i)
      sname += '_' + vname[i] ;
    vinfo.name = sname ;

    // Now substitute arguments
    vector<int> args = vinfo.v_ids ;
    for(size_t i=0;i<subst.size();++i) {
      if(subst[i].get_info().v_ids.size() != 0) {
        if(args.size() != 0) {
          cerr << "problem with substitution argument conflict when processing variable " << v << endl ;
          cerr << "parametric rule may have irregular substitution" << endl;
        }
        args = subst[i].get_info().v_ids ;
      }
    }
    for(size_t i=0;i<args.size();++i) {
      variable tmp(args[i]) ;
      tmp = applySubstitution(transform_map,tmp) ;
      args[i] = tmp.ident() ;
    }
    vinfo.v_ids = args ;

    return variable(vinfo) ;
  }
  
  rule instantiateRule(rule r, variable v) {
    variableSet target = r.targets() ;
    variable vt ;
    rule_implP rp = r.get_rule_implP() ;
    FATAL(rp == 0) ;
    if(rp->is_parametric_provided()) {
      vt = rp->get_parametric_variable() ;
    } else {
      if(target.size() != 1) {
        cerr << "can't handle parametric rules with multiple outputs!" << endl ;
        cerr << "rule = " << r << "variable " << v << endl ;
      }
      vt = *target.begin() ;
    }
    // Find substution rules
    map<string,variable> transform_map ;
    std::vector<int> alt = vt.get_arg_list() ;
    std::vector<int> aln = v.get_arg_list() ;

    if(alt.size() != aln.size()) {
      cerr << "arg lists for parametric rule " << r << " don't match variable " << v << endl ;
      cerr << "this should not happen!" << endl ;
      Loci::Abort() ;
    }
#ifdef VERBOSE
    cout << "transform_map =" ;
#endif
    for(size_t i=0;i<alt.size();++i) {
      variable tv(alt[i]) ;
      variable sv(aln[i]) ;
      transform_map[tv.get_info().name] = sv ;
#ifdef VERBOSE
      cout << " ("<< tv.get_info().name << "," << sv << ")" << endl ;
#endif
    }


    // variable renaming map
    map<variable,variable> vm ;
    // Now loop over variables and apply substitution rules
    variableSet all_vars = target ;
    all_vars += r.sources() ;
    variableSet:: const_iterator vi ;

    for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
      variable sv = applySubstitution(transform_map,*vi) ;
      if(sv != *vi) {
        vm[*vi] = sv ;
      }
    }

    rp->rename_vars(vm) ;
    return(rule(rp)) ;
  }

  ruleSet instantiateParametrics(param_rule_db &prule_db, variableSet pvars) {
    variableSet::const_iterator vi ;
    ruleSet newrules ;
    for(vi=pvars.begin();vi!=pvars.end();++vi) {
      variable v = *vi ;
      int psize = v.get_arg_list().size() ;
      param_rule_key rk(var2key(v),psize) ;
      param_rule_db::const_iterator pri = prule_db.find(rk) ;
      if(pri==prule_db.end()) {
        //        if(MPI_rank == 0)
        //          cerr << "WARNING: Unable to intantiate parametric rules for variable "
        //               << v << endl ;
        debugout << "Unable to intantiate parametric rules for variable "
                 << v << endl ;
        continue ;
      }
      ruleSet::const_iterator ri ;
      for(ri=pri->second.begin();ri!=pri->second.end();++ri) {
        newrules += instantiateRule(*ri,*vi) ;
      }
    }
    return newrules ;
  }
  
  rule_db parametric_rdb(const rule_db &rdb,variableSet query_vars) {
    Loci::rule_db par_rdb ;
    param_rule_db prule_db ; // Parametric rule db ;

    // put parametric rules in prule_db, all else in par_rdb
    ruleSet rset = rdb.all_rules() ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end();++rsi) {
      variableSet target = rsi->targets() ;
      bool param=false ;
      rule_implP rp = rsi->get_rule_implP() ;
      if(rp == 0)
        continue ;
      if(!rp->is_specialized()) {
        if(rp->is_parametric_provided()) {
          variable pv = rp->get_parametric_variable() ;
          int psize = pv.get_arg_list().size() ;
          param_rule_key rk(var2key(pv),psize) ;
          prule_db[rk] += *rsi ;
          param = true ;
        } else {
          for(variableSet::const_iterator vsi = target.begin(); vsi !=
                target.end(); ++vsi) { 
            int psize = variable(*vsi).get_arg_list().size() ;
            if(psize > 0) {
              param_rule_key rk(var2key(*vsi),psize) ;
              prule_db[rk] += *rsi ;
              param = true ;
            }
          }
        }
      }
      if(!param)
        par_rdb.add_rule(*rsi) ;
    }

    variableSet parvars = scanRulesForParametrics(par_rdb.all_rules()) ;
    for(variableSet::const_iterator vi=query_vars.begin();vi!=query_vars.end();++vi) {
      if(variable(*vi).get_arg_list().size() > 0) {
        parvars += *vi ;
      }
    }
    
#ifdef VERBOSE
    std::cout << "parvars = " << parvars << endl ;
#endif
    
    ruleSet processed ;
    int cnt = 0 ;
    while(parvars != EMPTY) {
#ifdef VERBOSE
      cout << "prcessing variables " << parvars << endl ;
#endif
      processed += parvars ;
      ruleSet newrules = instantiateParametrics(prule_db,parvars) ;

#ifdef VERBOSE
      cout << "generating rules = " << endl << newrules << endl ;
#endif
      parvars = scanRulesForParametrics(newrules) ;
      ruleSet::const_iterator ri ;
      for(ri=newrules.begin();ri!=newrules.end();++ri)
        par_rdb.add_rule(*ri) ;
      parvars -= processed ;
      cnt++ ;
      if(cnt == 100) {
        if(MPI_rank == 0) 
          cerr << "Warning, parametric rule instatiation 100 levels deep!" << endl
               << "Probably this is caused by a recursive parametric rule" << endl
               << "currently working on the parametric variables:" << endl
               << parvars << endl ;
      }
    }
    
    return par_rdb ;
  }
}

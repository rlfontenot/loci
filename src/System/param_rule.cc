#include <param_rule.h>

using std::string ;
using std::ostringstream ;
namespace Loci {

  rule prepend_time_level(Loci::rule r, string prepend) {
    std::map<variable, variable> vm ;
    rule_implP rp = r.get_rule_implP() ;
    variableSet vset = rp->get_var_list() ;
    
    for(variableSet::const_iterator vsi=vset.begin(); vsi!=vset.end(); ++vsi) {
      variable v = variable(*vsi) ;
      time_ident vt = v.time() ;
      time_ident vtnew(prepend,vt) ;
      variable vnew = v.change_time(vtnew) ;
      vm[v] = vnew ; 
    }
    rp->rename_vars(vm) ;
    return rule(rp) ;
  }
  
  std::string get_complete_name_base(variable v) {
    std::string name ;
    std::vector<std::string> ns_vec = v.get_info().get_namespace() ;
    for(size_t i = 0; i < ns_vec.size(); ++i) {
      name.append(ns_vec[i]) ;
      name.append("@") ;
    }
    name.append(v.get_info().name) ;
    return name ;
  }

  std::string get_complete_name(variable v) {
    std::string name ;
    std::vector<std::string> prior_vec = v.get_info().priority  ;
    for(size_t i = 0; i < prior_vec.size(); ++i) {
      name.append(prior_vec[i]) ;
      name.append("::") ;
    }
    std::vector<std::string> ns_vec = v.get_info().get_namespace() ;
    for(size_t i = 0; i < ns_vec.size(); ++i) {
      name.append(ns_vec[i]) ;
      name.append("@") ;
    }
    name.append(v.get_info().name) ;
    return name ;
  }
 
  std::string get_complete_name_wot(variable v) { 
    std::string name ;
    std::vector<std::string> prior_vec = v.get_info().priority  ;
    for(size_t i = 0; i < prior_vec.size(); ++i) {
      name.append(prior_vec[i]) ;
      name.append("::") ;
    }
    std::vector<std::string> ns_vec = v.get_info().get_namespace() ;
    for(size_t i = 0; i < ns_vec.size(); ++i) {
      name.append(ns_vec[i]) ;
      name.append("@") ;
    }
    name.append(v.get_info().name) ;
    std::vector<int> v_ids = v.get_arg_list() ; 
    ostringstream oss ;
    if(v_ids.size() != 0) {
      oss << "(";
      std::vector<int>::const_iterator vi = v_ids.begin() ;
      variable temp = variable(*vi) ;
      temp.get_info().Print(oss) ;
      ++vi ;
      for(;vi != v_ids.end(); ++vi) {
	oss << ",";
	temp = variable(*vi) ;
	temp.get_info().Print(oss) ;
      }
      oss << ")" ;
    }
    name.append(oss.str()) ;
    return name ; 
  }
  rule rename_rule(rule r, std::map<variable, variable> &vm) {

    std::vector<string> str_vec ;//This is to deal with the renaming
    //of other variables eg. W_f, f_W etc
    std::set<string> ren_tars ; // Create a set  of names of
    // renamed target variables  
    variableSet source = r.sources() ;
    for(variableSet::const_iterator tvsi = source.begin(); tvsi != source.end(); ++tvsi ) { 
      string name = get_complete_name(variable(*tvsi)) ;
      if(((name.find("_")) != string::npos) && (!variable(*tvsi).get_arg_list().size())) {
        for(std::map<variable, variable>::const_iterator mi = vm.begin(); mi != vm.end(); ++mi) {
          std::string sub_name = get_complete_name(mi->first) ;
	  std::string::size_type i ;
          std::string tmp_name = name ;
          do {
            i = tmp_name.find(sub_name) ;
            if(i != string::npos) {
	      if(i > 0) {
                if(tmp_name[i-1] == '_') {
		  str_vec.push_back(name) ;
                  tmp_name.erase(i-1, i+sub_name.size()+1) ;
                }
                else if((i+sub_name.size()) < tmp_name.size())
                  if(tmp_name[i+sub_name.size()] == '_') {
                    str_vec.push_back(name) ;
                    tmp_name.erase(i, i+sub_name.size()+1) ;
                  }
              }
              else {
                if((i+sub_name.size()) < tmp_name.size())
                  if(tmp_name[i+sub_name.size()] == '_') 
                    str_vec.push_back(name) ;
                tmp_name.erase(i, i+sub_name.size()+1) ;
              }
            }
          }while(i != std::string::npos) ;
        }
      }
    }
    variableSet target = r.targets() ;
    for(variableSet::const_iterator tvsi = target.begin(); tvsi
          != target.end(); ++tvsi ) { 
      variable var = variable(*tvsi) ;
      std::vector<int> tvint = var.get_arg_list() ; 
      for(std::vector<int>::const_iterator tvi = tvint.begin(); tvi != tvint.end(); ++tvi) 
        ren_tars.insert(get_complete_name(variable(*tvi))) ;
    }
    for(std::vector<string>::const_iterator svs =
          str_vec.begin(); svs != str_vec.end(); ++svs) {
      string name = *svs ;
      string orig = name ;
      string replace ;
      for(std::set<string>::const_iterator ssi =
            ren_tars.begin(); ssi != ren_tars.end(); ++ssi)
        {
          std::string::size_type i ;
          variable s = variable(*ssi) ; 
	  std::string tmp_name = name ;
          std::string final_name ;
          std::string::size_type pos = 0 ;
          do {
            i = tmp_name.find(*ssi) ;
            if(i != string::npos)  {
              std::string ts = tmp_name.substr(0, i);
              final_name.insert(pos, ts) ;
              pos += i ;
              if(i > 0) {
                if(tmp_name[i-1] == '_') 
                  if((i+(*ssi).size()) < tmp_name.size()) {
                    if(tmp_name[i+(*ssi).size()] == '_') {
                      replace = get_complete_name(vm[s]); //vm[s].get_info().name ;
		      name = tmp_name.substr(i+(*ssi).size(), tmp_name.size()) ;
                      tmp_name = name ;
                      final_name.insert(pos, replace) ;
                      pos += replace.size() ;
                    }
                    else {
                      final_name.insert(pos, (*ssi)) ;
                      pos += (*ssi).size() ;
                      name = tmp_name.substr(i+(*ssi).size(), tmp_name.size()) ;
                      tmp_name = name ;
                    }
                  }
              }
              else {
                if((i+(*ssi).size()) < tmp_name.size())
                  if(tmp_name[i+(*ssi).size()] == '_') {
                    replace = get_complete_name(vm[s]) ; //.get_info().name ;
		    name = tmp_name.substr(i+(*ssi).size(), tmp_name.size()) ;
                    tmp_name = name ;
                    final_name.insert(pos, replace) ;
                    pos += replace.size() ;
                  }
                  else {
                    final_name.insert(pos, *ssi) ;
                    pos += (*ssi).size() ;
                    name = tmp_name.substr(i+(*ssi).size(), tmp_name.size()) ;
                    tmp_name = name ;
                  }
              }
            }
          } while(i != string::npos) ;
          final_name.append(tmp_name) ;
          name = final_name ;
          std::vector<int> rep = vm[s].get_arg_list() ;
          if(rep.size()) {
            string par ;
            int p_open = 0 ;
            par.append("(") ;
            p_open++ ;
            std::vector<int> trep = rep ;
            while(trep.size()) {
              trep.clear() ;
              for(std::vector<int>::const_iterator pi = rep.begin() ; pi != rep.end(); ++pi) {
                trep = variable(*pi).get_arg_list() ;
                par.append("(") ;
                par.append(get_complete_name(variable(*pi))) ;
                p_open++ ;
              }
              rep = trep ;
            }
	  
            while(p_open) {
              par.append(")") ;
              p_open-- ;
            }
            name.append(par) ;
          }
        }
      if(name != orig) {

        vm[variable(orig)] = variable(name) ;
	//std::cout << "vm[ " << variable(orig) << " ] = " << variable(name) << endl ;
      }
    }
    rule_implP rp = r.get_rule_implP() ;
    rp->rename_vars(vm) ;
    return(rule(rp)) ;
  }

  variable recursive_rename(variable &v, std::map<variable, variable> &vm) {
    if(!v.get_arg_list().size()) {
      if(vm.find(v) != vm.end())
        return(vm[v]) ;
      else
        return v ;
    }
    std::vector<int> vint = v.get_arg_list();
    std::vector<int> tmp_vint ;
    for(size_t i = 0 ; i < vint.size(); ++i) {
      variable tmp = variable(vint[i]) ;
      variable tmp_var = recursive_rename(tmp, vm) ;
      tmp_vint.push_back(tmp_var.ident()) ;
    }
    return(variable(v, tmp_vint)) ;
  }
  
  rule_db parametric_rdb(const rule_db &rdb,variableSet query_vars) {
    Loci::rule_db par_rdb ;
    ruleSet param_target, added_rules ;
    variableSet target ;
    std::map<variable, ruleSet> mruleset ; // This data structure is
    // used to map the parametric variable with the corresponding
    // parametric rule in the rule database.
  
    std::map<std::string, variableSet> mvarset ; //This data structure
    // is used to store the parametric variables corresponding to a
    // particular name. eg. grad->grad(X), grad(y), grad(T), grad(X,y)
    // etc.   
    std::vector<int> vint ;
    ruleSet rset = rdb.all_rules() ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end();++rsi) {
      target = rsi->targets() ;
      std::vector<std::vector<int> > target_args ;
      for(variableSet::const_iterator vsi = target.begin(); vsi !=
            target.end(); ++vsi) { 
        vint = variable(*vsi).get_arg_list() ;
        if(vint.size()) {
          target_args.push_back(vint) ;
          variable tmp = variable(*vsi) ;
          mvarset[get_complete_name_base(tmp)] += *vsi ;
          mruleset[tmp] += *rsi ;
        }
      }
      if(target_args.size()) {
        std::vector<std::vector<int> >::const_iterator vi =
          target_args.begin(); 
        vint = *vi ;
	//#ifdef DEBUG
        //for( ; vi != target_args.end(); ++vi) 
	//FATAL(vint != *vi) ;
	//#endif      
        // Make a ruleSet of all the rules having parametric variables
        // as the target 
        if(vint.size() && (!param_target.inSet(*rsi))) 
          param_target += *rsi ;
      }
      
    }
  
    //Remove the rules having parametric variables in the head a rule
    //from the rule database. 
    rset -= param_target ;
  
    //Loop over the sources of the remaining rules and find out the
    //parametric variables(add them to param_vars). Add the
    //corresponding rule to use_param_rule ruleSet.   
  
    variableSet all_source_vars = query_vars ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end(); ++rsi) {
      all_source_vars += rsi->sources() ;
    }
  
    variableSet param_vars ;
    for(variableSet::const_iterator vsi = all_source_vars.begin();
        vsi != all_source_vars.end(); ++vsi) {
      if(variable(*vsi).get_arg_list().size() != 0) {
        param_vars += *vsi ;
      }
    }

    //This part handles the recursive parametric rules. 
    ruleSet wrule, nrule ;
    variableSet newvars ;
    if(param_vars != EMPTY) {
      for(variableSet::const_iterator vsi = param_vars.begin(); vsi !=
            param_vars.end(); ++vsi) {
        ruleSet wrule ;
        if(mruleset.find(variable(*vsi)) != mruleset.end())
          wrule = mruleset.find(variable(*vsi))->second ;
        nrule = wrule ;
        while(nrule != EMPTY) {
          nrule = EMPTY ;
          for(ruleSet::const_iterator rsi = wrule.begin(); rsi != wrule.end();
              ++rsi) {
            variableSet source = rsi->sources() ;
            for(variableSet::const_iterator vci = source.begin(); vci !=
                  source.end(); ++vci) {
              if(variable(*vci).get_arg_list().size()) {
                variableSet vset ;
                if(mvarset.find(get_complete_name(variable(*vci))) != mvarset.end())
                  vset = mvarset.find(get_complete_name(variable(*vci)))->second ;
                for(variableSet::const_iterator ivci = vset.begin(); ivci !=
                      vset.end(); ++ivci) {
                  if(variable(*vci) == variable(*ivci))  
                    if(mruleset.find(variable(*ivci)) != mruleset.end()) {
                      nrule = mruleset.find(variable(*ivci))->second ;
                      for(ruleSet::const_iterator irsi = nrule.begin(); irsi 
                            != nrule.end(); ++irsi) {
                        for(variableSet::const_iterator tvsi =
                              irsi->sources().begin(); tvsi != irsi
                              ->sources().end(); ++tvsi)
                          if(variable(*tvsi).get_arg_list().size()) {
                            if(mvarset.find(get_complete_name(variable(*tvsi))) != mvarset.end())
                              if(!mvarset.find(get_complete_name(variable(*tvsi)))->second.inSet(*tvsi)) {
                                newvars += *vci ;
                                std::vector<int> vin = variable(*tvsi).get_arg_list() ;
                                if(vin.size()) 
                                  if(!variable(*(vin.begin())).get_arg_list().size()) 
                                    newvars += *tvsi ;
                              }
                          }
                      }
                    }
                }
              }
            }
          }
          wrule = nrule ;
        }
      }
    }
    param_vars += newvars ;
  
    //Add the non parametric rules to the rule database.  
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end(); ++rsi) {
      rule_implP rp = rsi->get_rule_implP() ;
      par_rdb.add_rule(rule(rp)) ;
      added_rules += rule(rp) ;
    }
  
    variableSet working, newset ;
    variableSet processed ;
    working = param_vars;
    newset = param_vars ;
    while(newset != EMPTY) {
      newset = EMPTY ;
      for(variableSet::const_iterator vsi = working.begin(); vsi != working.end(); ++vsi) {
        std::map<variable, variable> vm ;
        variable v = variable(*vsi) ;
        variableSet vs;
        ruleSet rs ;
	if(v.time() != time_ident())
	  v = variable(get_complete_name_wot(v)) ;
	if(mvarset.find(get_complete_name(v)) != mvarset.end()) 
          vs = mvarset.find(get_complete_name(v))->second ;
        for(variableSet::const_iterator mvsi = vs.begin(); mvsi != vs.end(); ++mvsi)
          if(v.get_arg_list().size() == variable(*mvsi).get_arg_list().size()) {
            if(mruleset.find(variable(*mvsi)) != mruleset.end())
              rs = mruleset.find(variable(*mvsi))->second ;
            else
              cerr << "can't find parametric rule for variable " << variable(*mvsi) << endl ;
            for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi ) {
              target = rsi->targets() ;
              for(variableSet::const_iterator tvsi = target.begin(); tvsi
                    != target.end(); ++tvsi ) { 
                variable var = variable(*tvsi) ;
                std::vector<int> tvint = var.get_arg_list() ; 
                vint = v.get_arg_list()  ;
                variable v2 = v ;
                if(var.get_info().priority.size() != 0) {
                  variable::info vinfo = v.get_info() ;
                  vinfo.priority = var.get_info().priority ;
                  v2 = variable(vinfo) ;
                }
                vm[var] = v2 ;
                std::vector<int>::const_iterator tvi = tvint.begin() ; 
                for(std::vector<int>::const_iterator vi = vint.begin(); vi
                      != vint.end(), tvi != tvint.end(); ++vi, ++tvi) {
                  if(variable(*vi).get_arg_list().size()) 
                    newset += variable(*vi) ;
                  if(vm.find(variable(*tvi)) == vm.end()) {
                    vm[variable(*tvi)] = variable(*vi) ;
                  }
                }
              }
              variableSet source = rsi->sources() ;
              for(variableSet::const_iterator tvsi = source.begin(); tvsi != source.end(); ++tvsi ) {    
                std::vector<int> tmp_vint = variable(*tvsi).get_arg_list() ;
                variable tmp_var ;
                variable v = variable(*tvsi) ;
                tmp_var = recursive_rename(v, vm) ; 
                if(v != tmp_var) {
                  vm[v] = tmp_var ;
                  if(tmp_var.get_arg_list().size())
                    newset += tmp_var ;
                }
                else
                  if(v.get_arg_list().size())
                    newset += v ;
              }
              rule r = rename_rule(*rsi,vm) ;
              rule_implP rp = r.get_rule_implP() ;
              if(!added_rules.inSet(rule(rp))) {
                par_rdb.add_rule(rule(rp)) ;
                added_rules += rule(rp) ;
              }
            }
          }
      }
      // If a variable has already been processed, no need to process it again
      // (This is needed to make sure that things work when parametric rules
      // are recursive
      newset -= processed ;
      // Update working set to remainder
      working = newset ;
      // Next time we get here, working will have been processed
      processed += working ;
    }
    return par_rdb ;
  }
  
}

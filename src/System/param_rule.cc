#include <param_rule.h>

namespace Loci {
  rule_db parametric_rdb(rule_db& rdb) {
    Loci::rule_db par_rdb ;
    ruleSet param_target, param_source, use_param_rule, added_rules ;
    variableSet source, target ;
    std::map<variable, ruleSet> mruleset ;
    std::map<std::string, variableSet> mvarset ;
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
	  mvarset[tmp.get_info().name] += *vsi ;
	  mruleset[tmp] += *rsi ;
	}
      }
      if(target_args.size()) {
	std::vector<std::vector<int> >::const_iterator vi =
	  target_args.begin(); 
	vint = *vi ;
	for( ; vi != target_args.end(); ++vi) 
	  FATAL(vint != *vi) ;
	
	
	if(vint.size() && (!param_target.inSet(*rsi))) {
	  param_target += *rsi ;
	  //cout << "removing rule  " << *rsi << endl ;
	}
      }
      
      source = rsi->sources() ;
      for(variableSet::const_iterator vsi = source.begin(); vsi !=
	    source.end(); ++vsi) { 
	vint = variable(*vsi).get_arg_list() ;
	if(vint.size()) 
	  param_source += *rsi ;
      }
    }
    rset -= param_target ;
    variableSet param_vars ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end(); ++rsi) {
      source = rsi->sources() ;
      for(variableSet::const_iterator vsi = source.begin(); vsi != source.end(); ++vsi) {
	vint = variable(*vsi).get_arg_list() ;
	if(vint.size()) {
	  use_param_rule += *rsi ;
	  param_vars += *vsi ;
	}
      }
    }
    
    param_source -= use_param_rule ;
    for(ruleSet::const_iterator rsi = param_source.begin(); rsi != param_source.end(); ++rsi) { 
      target = rsi->targets() ;
      for(variableSet::const_iterator vsi = target.begin(); vsi !=
	    target.end(); ++vsi) { 
	vint = variable(*vsi).get_arg_list() ;
	if(vint.size()) {
	  variable tmp = variable(*vsi) ;
	  mvarset[tmp.get_info().name] += *vsi ;
	  mruleset[tmp] += *rsi ;
	}
      }
      source = rsi->sources() ;
      for(variableSet::const_iterator vsi = source.begin(); vsi !=
	    source.end(); ++vsi) 
	vint = variable(*vsi).get_arg_list() ;
    }
    
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end(); ++rsi) {
      rule_implP rp = rsi->get_rule_implP() ;
      par_rdb.add_rule(rule(rp)) ;
      added_rules += rule(rp) ;
    }
    
    
    variableSet working, newset ;
    working = param_vars;
    newset = param_vars ;
    while(newset != EMPTY) {
      std::vector<string> str_vec ;
      newset = EMPTY ;
      for(variableSet::const_iterator vsi = working.begin(); vsi != working.end(); ++vsi) {
	variable v = variable(*vsi) ;
	variableSet vs = mvarset.find(v.get_info().name)->second ;
	for(variableSet::const_iterator mvsi = vs.begin(); mvsi != vs.end(); ++mvsi)
	  if(v.get_arg_list().size() == variable(*mvsi).get_arg_list().size()) {
	    ruleSet rs = mruleset.find(variable(*mvsi))->second ;
	    std::set<string> ren_tars ;
	    for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi ) {
	      std::map<variable, variable> vm ;
	      target = rsi->targets() ;
	      for(variableSet::const_iterator tvsi = target.begin(); tvsi
		    != target.end(); ++tvsi ) { 
		variable var = variable(*tvsi) ;
		std::vector<int> tvint = var.get_arg_list() ; 
		vint = v.get_arg_list()  ;
		vm[var] = v ;
		std::vector<int>::const_iterator tvi = tvint.begin() ; 
		for(std::vector<int>::const_iterator vi = vint.begin(); vi
		      != vint.end(), tvi != tvint.end(); ++vi, ++tvi) {
		  if(variable(*vi).get_arg_list().size()) 
		    newset += variable(*vi) ;
		  vm[variable(*tvi)] = variable(*vi) ;
		  ren_tars.insert(variable(*tvi).get_info().name) ;
		}
	      }
	      source = rsi->sources() ;
	      for(variableSet::const_iterator tvsi = source.begin(); tvsi != source.end(); ++tvsi ) { 
		std::vector<int> tmp_vint = variable(*tvsi).get_arg_list() ;
		std::vector<int> tmp_vec, vec_int ;
		vec_int = tmp_vint ;
		vec_int.push_back(variable(*tvsi).ident()) ;
		while(tmp_vint.size()) {
		  tmp_vint.clear() ;
		  for(std::vector<int>::const_iterator vi = vec_int.begin(); vi != vec_int.end(); ++vi) { 
		    variable tmp = variable(*vi) ;
		    std::vector<int> pi = tmp.get_arg_list() ;
		    if(pi.size())
		      for(std::vector<int>::const_iterator vp =
			    pi.begin(); vp != pi.end(); ++vp)
			tmp_vint.push_back(variable(*vp).ident()) ;
		    else {
		      std::map<variable, variable>::const_iterator mi = vm.find(tmp) ;
		      if(mi != vm.end())
			vm[tmp] = mi->second ;
		    }
		  }
		  std::vector<int> ulti ;
		  for(std::vector<int>::const_iterator vi = tmp_vint.begin(); vi != tmp_vint.end(); ++vi) { 
		    variable tmp = variable(*vi) ;
		    tmp_vec = tmp.get_arg_list() ;
		    std::vector<int> vec ;
		    for(std::vector<int>::const_iterator ivi = tmp_vec.begin(); ivi != tmp_vec.end(); ++ivi) {
		      variable tmp_var = variable(*ivi) ;
		      std::map<variable, variable>::const_iterator mi = vm.find(tmp_var) ;
		      if(mi != vm.end()) 
			vec.push_back(mi->second.ident()) ;
		      else
			vec.push_back(tmp_var.ident()) ;
		    }
		    if(tmp_vec != vec) {
		      
		      vm[variable(tmp)] = variable(tmp, vec) ;
		    }
		    else {
		      ulti.push_back(tmp.ident()) ;
		    }
		  }
		  vec_int = ulti ;
		}
		tmp_vint = variable(*tvsi).get_arg_list() ;
		std::vector<int> vec ;
		for(std::vector<int>::const_iterator vi =
		      tmp_vint.begin(); vi != tmp_vint.end(); ++vi) {
		  variable tmp_var = variable(*vi) ;
		  std::map<variable, variable>::const_iterator mi =
		    vm.find(tmp_var) ;
		  if(mi != vm.end())
		    vec.push_back(mi->second.ident()) ;
		  else
		    vec.push_back(tmp_var.ident()) ;
		}
		if(tmp_vint != vec) { 
		  vm[variable(variable(*tvsi))] =
		    variable(variable(*tvsi), vec) ;
		  newset +=  variable(variable(*tvsi), vec) ;
		}
		string name = variable(*tvsi).get_info().name ;
		if((name.find("_")) != string::npos)
		  str_vec.push_back(name) ;
	      }
	      for(std::vector<string>::const_iterator svs =
		    str_vec.begin(); svs != str_vec.end(); ++svs) {
		string name = *svs ;
		string orig = name ;
		string replace ;
		for(std::set<string>::const_iterator ssi =
		      ren_tars.begin(); ssi != ren_tars.end(); ++ssi) {
		  int i = name.find(*ssi) ;
		  variable s ;
		  if(i != string::npos)  {
		    s = variable(*ssi) ;
		    replace = vm[s].get_info().name ;
		    name.erase(i, i+(*ssi).size()) ;
		    name.insert(i, replace) ;
		  }
		}
		if(name != orig)
		  vm[variable(orig)] = variable(name) ;
	      }
	      rule_implP rp = rsi->get_rule_implP() ;
	      rp->rename_vars(vm) ;
	      if(!added_rules.inSet(rule(rp))) {
		par_rdb.add_rule(rule(rp)) ;
		//cout << "adding rule  " << rule(rp) << endl ;
		added_rules += rule(rp) ;
	      }
	    }
	  }
      }
      working = newset ;
      
    }
    return par_rdb ;
  }
}

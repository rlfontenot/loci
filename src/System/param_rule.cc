#include <param_rule.h>

namespace Loci {

rule_db local_modify_time_vars(rule_db& rdb, const std::string &sn) {
  rule_db new_rdb ;
  ruleSet rset = rdb.all_rules() ;
  for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end(); ++rsi) {
    variableSet vset, tvarset ;
    std::map<variable, variable> vm ;
    rule_implP rp = rsi->get_rule_implP() ;
    vset = rp->get_var_list() ;
    for(variableSet::const_iterator vsi = vset.begin(); vsi !=
	  vset.end(); ++vsi ) {
      if((variable(*vsi).time()).level_name() != "*") 
	tvarset += *vsi ;
    }
    for(variableSet::const_iterator vsi = tvarset.begin(); vsi != tvarset.end(); ++vsi) {
      variable v = variable(*vsi) ;
      ostringstream oss ;
      v.get_info().Print(oss) ;
      std::string name = oss.str() ;
      std::string new_var ;
      new_var.append(v.get_info().name) ;
      new_var.append("{") ;
      new_var.append(sn) ;
      new_var.append(",");
      int i = name.find("{") ;
      i++ ;
      string sub = name.substr(i, name.size()-i) ;
      new_var.append(sub) ;
      variable tmp = variable(new_var) ;
      vm[v] = tmp ;
    }
    if(tvarset != EMPTY)
      rp->rename_vars(vm) ;
    new_rdb.add_rule(rule(rp)) ;
  }
  return new_rdb ;
}
  
rule_db parametric_rdb(rule_db& rdb) {
  Loci::rule_db par_rdb ;
  ruleSet param_target, param_source, use_param_rule, added_rules ;
  variableSet source, target ;
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
	mvarset[tmp.get_info().name] += *vsi ;
	mruleset[tmp] += *rsi ;
      }
    }
    if(target_args.size()) {
      std::vector<std::vector<int> >::const_iterator vi =
	target_args.begin(); 
      vint = *vi ;
#ifdef DEBUG
      for( ; vi != target_args.end(); ++vi) 
	FATAL(vint != *vi) ;
#endif      
      // Make a ruleSet of all the rules having parametric variables
      // as the target 
      if(vint.size() && (!param_target.inSet(*rsi))) 
	param_target += *rsi ;
    }
    
    source = rsi->sources() ;
    for(variableSet::const_iterator vsi = source.begin(); vsi !=
	  source.end(); ++vsi) { 
      vint = variable(*vsi).get_arg_list() ;
      //Make a ruleSet of all the rules having paramtric variables as
      //the source 
      if(vint.size()) 
	param_source += *rsi ;
    }
  }
  //Remove the rules having parametric variables in the head a rule
  //from the rule database. 
  rset -= param_target ;
  variableSet param_vars ;
  
  //Loop over the sources of the remaining rules and find out the
  //parametric variables(add them to param_vars). Add the
  //corresponding rule to use_param_rule ruleSet.   
  
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
  
  //Remove the rules which uses the parametric rules(use_param_rule)
  //from the param_source ruleSet. The remaining ruleSet will have the 
  //parametric rules with variables(that could be parametric) in the source.
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
  }
  variableSet newvars ;
  //This part handles the recursive parametric rules. 
  ruleSet wrule, nrule ;
  if(param_vars != EMPTY) {
    for(variableSet::const_iterator vsi = param_vars.begin(); vsi !=
	  param_vars.end(); ++vsi) {
      ruleSet wrule ;
      //Only if there is a recursive parametric rule there will be
      //a ruleSet to loop over. 
      if(mruleset.find(variable(*vsi)) != mruleset.end())
	wrule = mruleset.find(variable(*vsi))->second ;
      nrule = wrule ;
      while(nrule != EMPTY) {
	nrule = EMPTY ;
	for(ruleSet::const_iterator rsi = wrule.begin(); rsi != wrule.end();
	    ++rsi) {
	  source = rsi->sources() ;
	  for(variableSet::const_iterator vci = source.begin(); vci !=
		source.end(); ++vci) {
	    if(variable(*vci).get_arg_list().size()) {
	      variableSet vset ;
	      if(mvarset.find(variable(*vci).get_info().name) != mvarset.end())
		vset = mvarset.find(variable(*vci).get_info().name)->second ;
	      for(variableSet::const_iterator ivci = vset.begin(); ivci !=
		    vset.end(); ++ivci) {
		if(mruleset.find(variable(*ivci)) != mruleset.end()) {
		  nrule = mruleset.find(variable(*ivci))->second ;
		  for(ruleSet::const_iterator irsi = nrule.begin(); irsi 
			!= nrule.end(); ++irsi) {
		    for(variableSet::const_iterator tvsi =
			  irsi->sources().begin(); tvsi != irsi
			  ->sources().end(); ++tvsi)
		      if(variable(*tvsi).get_arg_list().size())
			if(mvarset.find(variable(*tvsi).get_info().name) != mvarset.end()) 
			  if(!mvarset.find(variable(*tvsi).get_info().name)->second.inSet(*tvsi))			
			    newvars += *vci ;
		    
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
  working = param_vars;
  newset = param_vars ;
  while(newset != EMPTY) {
    std::vector<string> str_vec ;//This is to deal with the renaming
    //of other variables eg. W_f, f_W etc
    newset = EMPTY ;
    std::map<variable, variable> vm ;
    for(variableSet::const_iterator vsi = working.begin(); vsi != working.end(); ++vsi) {
      variable v = variable(*vsi) ;
      variableSet vs = mvarset.find(v.get_info().name)->second ;
      for(variableSet::const_iterator mvsi = vs.begin(); mvsi != vs.end(); ++mvsi)
	if(v.get_arg_list().size() == variable(*mvsi).get_arg_list().size()) {
	  ruleSet rs = mruleset.find(variable(*mvsi))->second ;
	  std::set<string> ren_tars ; // Create a set  of names of
	  // renamed target variables  
	  for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi ) {
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
		if(vm.find(variable(*tvi)) == vm.end()) {
		  vm[variable(*tvi)] = variable(*vi) ;
		}
		ren_tars.insert(variable(*tvi).get_info().name) ;
	      }
	    }
	    source = rsi->sources() ;
	    for(variableSet::const_iterator tvsi = source.begin(); tvsi != source.end(); ++tvsi ) { 
	      std::vector<int> tmp_vint = variable(*tvsi).get_arg_list() ;
	      std::vector<int> tmp_vec, vec_int ;
	      vec_int = tmp_vint ;
	      while(tmp_vint.size()) {
		tmp_vint.clear() ;
		for(std::vector<int>::const_iterator vi = vec_int.begin(); vi != vec_int.end(); ++vi) { 
		  variable tmp = variable(*vi) ;
		  std::vector<int> pi = tmp.get_arg_list() ;
		  if(pi.size()) {
		    tmp_vint.push_back(tmp.ident()) ;
		    for(std::vector<int>::const_iterator vp =
			  pi.begin(); vp != pi.end(); ++vp) {
		      tmp_vint.push_back(variable(*vp).ident()) ;
		    }
		  }
		  else {
		    std::map<variable, variable>::const_iterator mi = vm.find(tmp) ;
		    if(mi != vm.end())
		      vm[tmp] = mi->second ;
		  }
		}
		std::vector<int> ulti = tmp_vint;
		tmp_vint.clear() ;
		for(std::vector<int>::const_iterator vi = ulti.begin(); vi != ulti.end(); ++vi) { 
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
		  if(tmp_vec != vec) 
		    vm[variable(tmp)] = variable(tmp, vec) ;
		  
		  else 
		    tmp_vint.push_back(tmp.ident()) ;
		}
		vec_int = tmp_vint ;
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
	      if(((name.find("_")) != string::npos) && (!variable(*tvsi).get_arg_list().size()))
		str_vec.push_back(name) ;
	    }
	    
	    for(std::vector<string>::const_iterator svs =
		  str_vec.begin(); svs != str_vec.end(); ++svs) {
	      string name = *svs ;
	      string orig = name ;
	      string replace ;
	      for(std::set<string>::const_iterator ssi =
		    ren_tars.begin(); ssi != ren_tars.end(); ++ssi)
		{
		  int i ;
		  variable s = variable(*ssi) ; 
		  do {
		    i = name.find(*ssi) ;
		    if(i != string::npos)  {
		      replace = vm[s].get_info().name ;
		      name.erase(i, i+(*ssi).size()) ;
		      name.insert(i, replace) ;
		    }
		  }while(i != string::npos) ;
		  
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
			par.append(variable(*pi).get_info().name) ;
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
	      if(name != orig)
		vm[variable(orig)] = variable(name) ;
	    }
	    
	    rule_implP rp = rsi->get_rule_implP() ;
	    rp->rename_vars(vm) ;
	    if(!added_rules.inSet(rule(rp))) {
	      par_rdb.add_rule(rule(rp)) ;
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

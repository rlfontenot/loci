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
	  //cout << " name of variable in the map =  " << tmp.get_info().name << endl ;
	  mvarset[tmp.get_info().name] += *vsi ;
	  mruleset[tmp] += *rsi ;
	  
	}
      }
      if(target_args.size()) {
	std::vector<std::vector<int> >::const_iterator vi =
	  target_args.begin(); 
	vint = *vi ;
	for( ; vi != target_args.end(); ++vi) {
	  FATAL(vint != *vi) ;
	} 
	
	if(vint.size() && (!param_target.inSet(*rsi))) {
	  //cout << "removing rule due to target.... " << *rsi << endl ;
	  param_target += *rsi ;
      }
      }
      
      source = rsi->sources() ;
      for(variableSet::const_iterator vsi = source.begin(); vsi !=
	    source.end(); ++vsi) { 
	vint = variable(*vsi).get_arg_list() ;
	if(vint.size()) {
	  param_source += *rsi ;
	}
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
	    source.end(); ++vsi) {
	vint = variable(*vsi).get_arg_list() ;
	if(vint.size()) 
	  param_vars += *vsi ;
      }
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
      newset = EMPTY ;
      for(variableSet::const_iterator vsi = working.begin(); vsi != working.end(); ++vsi) {
	variable v = variable(*vsi) ;
	variableSet vs = mvarset.find(v.get_info().name)->second ;
	for(variableSet::const_iterator mvsi = vs.begin(); mvsi != vs.end(); ++mvsi)
	  if(v.get_arg_list().size() == variable(*mvsi).get_arg_list().size()) {
	    ruleSet rs = mruleset.find(variable(*mvsi))->second ;
	    for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi ) {
	      std::map<variable, variable> vm, tmpvm ;
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
		  tmpvm[variable(*vi)] = variable(*tvi) ;
		}
	      }
	      source = rsi->sources() ;
	      for(variableSet::const_iterator tvsi = source.begin(); tvsi != source.end(); ++tvsi ) { 
		std::vector<int> vec_int = v.get_arg_list();
		for(std::vector<int>::const_iterator vi = vec_int.begin(); vi != vec_int.end(); ++vi)   
		  vm[tmpvm[variable(*vi)]] = variable(*vi) ;
		
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

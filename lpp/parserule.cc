//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
/*
  technique notes:
  //rule_impl_type {POINTWISE,SINGLETON,UNIT,APPLY,DEFAULT,OPTIONAL,CONSTRAINT_RULE,MAP_RULE,BLACKBOX_RULE,INSERTION,DELETION,ERASE,SUPER_RULE,UNKNOWN}

  Loci preprocessor is restructured so that  a rule is read in first and will be output in different versions according to devices that the program is run on.
  class ruleparse is designed to store the data of a rule. The output of different device version can be writted into one file or different files,
  which is controlled by two variables:

  vector<pair<string, ostream*> > ofs
  bool oneOutput

  currently, we assume in case of mutiple output files, inside non-cpu files, for non-loci-specific code, only include and namespace part will be written. In the future, we need
  choose what go where depend on the scheduler and compiler. If only one output file/stream, even non-cpu version, all non-loci-specific code will be written.



  currently, for cuda version:
  1. replace the loci rule type name with cuda rule type name.( in function loci_type(string rule_type))
  2. all the member variables of a rule are placed by pointers. and multiMaps are replaced by two pointers, the extra pointer is to store the  offset of each element
  3. constructor in cpu version is replaced by cuda version constructors, name_store() are replaced by bind() function.
  4. added member functions such as getDomain(), setDomain()
  5. calculate() function in cpu version is replaced by operator() function

  After a rule is read in, the preprocessor will decide if it is cpu-only according the rule type and its member variables. The current criteria are:

  1. if its rule_type is singleton, default, optional, constraint or blackbox, then it is cpu_only;
  2. if it use_prelude or is_specialized , then it is cpu-only
  3. if it has constraints,  cpu_only = true;
  4. if it is  singletonApply,  cpu_only = true;
  5. if it has inplace variables,  cpu_only = true;
  6. if it use  multiStore , storeVec, or storeMat,   cpu_only = true;

  In the future:
  1. the criteria for cpu-only need to be refined
  2. inside the operator() function, the variable access need to be modified so that only one-dimensional arrays are used.
  3. tons of tedious work need to be done to filter what is/is not in non-cpu version.


*/

#include "parserule.h"

variable convertVariable(variable v) {
  variable::info vinfo = v.get_info() ;
  vinfo.priority = std::vector<std::string>() ;
  for(size_t i=0;i<vinfo.v_ids.size();++i) {
    std::ostringstream ss ;
    ss << 'X' << i << endl ;
    variable xi = variable(ss.str()) ;
    vinfo.v_ids[i] = xi.ident() ;
  }
  return variable(vinfo) ;
}





namespace {
  inline void fill_descriptors(set<vmap_info> &v, const exprList &in) {

    for(exprList::const_iterator i = in.begin();i!=in.end();++i) {
      // This needs to be improved to use an actual variable syntax
      // certification.  This test will just get the blindingly obvious
      if((*i)->op != OP_ARROW &&
	 (*i)->op != OP_NAME &&
	 (*i)->op != OP_FUNC &&
	 (*i)->op != OP_NAME_BRACE &&
	 (*i)->op != OP_FUNC_BRACE &&
	 (*i)->op != OP_SCOPE &&
	 (*i)->op != OP_AT &&
	 (*i)->op != OP_DOLLAR) {
	cerr << "malformed descriptor: " ;
	(*i)->Print(cerr) ;
	cerr << endl ;
	throw parseError("rule signature error") ;
      }
      vmap_info di(*i) ;
      if(v.find(di) != v.end())
	cerr << "Warning, duplicate variable in var set." << endl ;
      else
	v.insert(di) ;
    }
  }
}


//replace the loci rule type name with cuda rule type name
string cuda_type(const string& rule_type){
  if(rule_type == "apply") return "ApplyRule";
  if(rule_type == "pointwise") return "PointwiseRule";
  if(rule_type == "unit") return "UnitRule";

  //rules that will never run on GPU, such as default, singleton??
  cerr << " rule_type: " << rule_type << " should not have cuda version" << endl;
  return rule_type;

}

//replace the loci container type names with cuda pointer names
//multiStore, storeVec is still kept here for future use
string var2ptr(const pair<string, string>& var){
  string ctype = var.first; //container type
  string dtype = var.second; //data type
  string result;
  if(ctype == "store" || ctype == "param" || ctype == "multiStore" || ctype == "storeVec" || ctype == "storeMat"   ){
    result = dtype.substr(1, dtype.find_last_of('>')-1) + "*";
  }else if(ctype == "Map" || ctype == "multiMap" ){
    result = "Loci::int_type*";
  }else{
    std::cerr << " ERROR: container type " <<  var.first << endl;

  }
  return result;
}

string var2name(variable v) {
  string vn = v.str() ;
  string name ;
  if(!prettyOutput)
    name += "L_" ;
  for(size_t si=0;si!=vn.size();++si) {
    if(isalpha(vn[si]) || isdigit(vn[si]) || vn[si] == '_')
      name += vn[si] ;
    if(vn[si]=='{' || vn[si] == '}')
      name += '_' ;
    if(vn[si]=='=')
      name += "_EQ_" ;
    if(vn[si]=='+')
      name += "_P_" ;
    if(vn[si]=='-')
      name += "_M_" ;
  }
  if(!prettyOutput)
    name += "_" ;
  return name ;
}

// expand mapping list into all possible map strings
std::vector<list<variable> > expand_mapping(std::vector<variableSet> vset) {
  // if we have sliced off all of the variable sets, then the list is empty
  if(vset.size() == 0) {
    return std::vector<list<variable> >() ;
  }
  // get map set for the last item in the list
  variableSet mlast = vset.back() ;
  vset.pop_back() ;

  // expand remainder of list
  std::vector<list<variable> > tmp  = expand_mapping(vset) ;

  // Now build list by enumerating all maps from this level
  std::vector<list<variable> > tmp2 ;
  int tsz = tmp.size() ;
  if(tmp.size() == 0) {
    for(variableSet::const_iterator vi=mlast.begin();vi!=mlast.end();++vi) {
      list<variable> l1 ;
      l1.push_back(*vi) ;
      tmp2.push_back(l1) ;
    }
  } else {
    for(int i=0;i<tsz;++i) {
      for(variableSet::const_iterator vi=mlast.begin();vi!=mlast.end();++vi) {
	list<variable> l1= tmp[i] ;
	l1.push_back(*vi) ;
	tmp2.push_back(l1) ;
      }
    }
  }
  return tmp2 ;
}



istream &parserule::get(std::string& filename, int& cnt, int start_line, const std::map<Loci::variable,std::pair<std::string,std::string> > &type_map, istream &is) {
  //this function also serve as initializer
  first_line = start_line;
  cpu_only = false;
  use_compute = true ;
  output_param = false ;

  parsebase::killsp(is) ;
  if(is_name(is)) {
    rule_type = get_name(is) ;
  } else
    throw parseError("syntax error") ;
  signature.get(is) ;
  lines += signature.num_lines() ;
  parsebase::killsp(is) ;
  if(rule_type == "apply") {
    if(is.peek() != '[')
      throw parseError("apply rule missing '[operator]'") ;
    apply_op.get(is) ;
    lines += apply_op.num_lines() ;
    parsebase::killsp(is) ;
  }


  use_prelude = false ;
  is_specialized = false ;
  while(is.peek() == ',') {
    is.get() ;
    parsebase::killsp(is) ;
    if(!is_name(is))
      throw parseError("syntax error") ;

    string s = get_name(is) ;
    if(s == "constraint") {
      nestedparenstuff con ;
      con.get(is) ;
      if(constraint == "")
	constraint = con.str() ;
      else
	constraint += "," + con.str() ;
      lines += con.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      if(parametric_var != "") {
	throw parseError("syntax error: canot specify more than one parametric variable") ;
      }

      parametric_var = con.str() ;
      lines += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
	throw parseError("syntax error: canot specify more than one conditional variable") ;
      }
      conditional = con.str() ;
      lines += con.num_lines() ;
      // Check variable
      variable v(conditional) ;
      v = convertVariable(v) ;
      map<variable,pair<string,string> >::const_iterator mi ;
      if((mi = type_map.find(v)) == type_map.end()) {

	v = v.new_offset(0) ;
	v = v.drop_assign() ;
	while(v.time() != time_ident())
	  v = v.parent() ;

	if((mi = type_map.find(v)) == type_map.end()) {
	  while(v.get_info().namespac.size() != 0)
	    v = v.drop_namespace() ;
	  mi = type_map.find(v) ;
	}
      }
      if(mi == type_map.end()) {
	cerr << filename << ':' << start_line+lines << ":0: warning: type of conditional variable '" << v << "' not found!"  << endl  ;
      } else {

	// clean up type string
	string val = mi->second.first + mi->second.second ;
	string val2 ;
	int valsz = val.size() ;
	for(int i=0;i<valsz;++i)
	  if(val[i] != ' ' && val[i] != '\t' && val[i] != '\r' && val[i] != '\n')
	    val2 += val[i] ;

	if(val2 != "param<bool>") {
	  throw(parseError("conditional variable must be typed as a param<bool>")) ;
	}
      }
    } else if(s == "inplace") {
      nestedparenstuff ip ;
      ip.get(is) ;
      lines += ip.num_lines() ;
      exprP p = expression::create(ip.str()) ;
      exprList l = collect_associative_op(p,OP_OR) ;
      if(l.size() != 2)
	throw parseError("inplace needs two variables with a '|' separator") ;

      exprList::const_iterator i = l.begin() ;
      variable v1(*i) ;
      ++i ;
      variable v2(*i) ;
      inplace.push_back(pair<variable,variable>(v1,v2)) ;

    } else if(s == "prelude") {
      use_prelude=true ;
      parsebase::killsp(is) ;
      continue ;
    } else if(s == "specialized") {
      is_specialized = true ;
    } else if(s == "option") {
      nestedparenstuff ip ;
      ip.get(is) ;
      lines += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      lines += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    parsebase::killsp(is) ;
  }


  sig = signature.str() ;
  head=0;
  body=0 ;
  for(size_t i=0;i<sig.size()-1;++i) {
    if(sig[i]=='<' && sig[i+1]=='-') {
      heads = sig.substr(0,i) ;
      bodys = sig.substr(i+2,sig.size()) ;
      head = expression::create(heads) ;
      body = expression::create(bodys) ;
      if(rule_type == "optional" || rule_type == "default") {
	throw parseError("'optional' or 'default' rules should not have a body (defined by '<-' operator)!") ;
      }
    }
  }
  if(head == 0) {
    heads = sig ;
    head = expression::create(heads) ;
    if(rule_type == "optional" || rule_type == "default") {
      if(constraint != "")
	throw parseError("'optional' or 'default' rules should not have a constraint!") ;
    } else {
      if(constraint == "") {
	throw parseError("rules without bodies should have a defined constraint as input!") ;
      }
    }
  }


  class_name = "file_" ;
  for(size_t i=0;i<filename.size();++i) {
    char c = filename[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      class_name += c ;
    if(c == '.')
      break ;
  }
  class_name += '0' + (cnt/100)%10 ;
  class_name += '0' + (cnt/10)%10 ;
  class_name += '0' + (cnt)%10 ;
  timeb tdata ;
  ftime(&tdata) ;

  ostringstream tss ;
  tss <<  '_' << tdata.time << 'm'<< tdata.millitm;

  class_name += tss.str() ;
  cnt++ ;
#ifdef OLD
  class_name += "_rule_" ;
  if(conditional != "")
    sig += "_" + conditional ;
  if(constraint != "")
    sig += "_" + constraint ;
  for(size_t i=0;i<sig.size();++i) {
    if(isalpha(sig[i]) || isdigit(sig[i]))
      class_name += sig[i] ;
    if(sig[i] == ',' || sig[i] == '-' || sig[i] == '>' || sig[i] == '('||
       sig[i] == ')' || sig[i] == '{' || sig[i] == '}' || sig[i] == '='||
       sig[i] == '+' || sig[i] == '_')
      class_name += '_' ;
  }
#endif

  //input, output
  if(body != 0)
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;

  set<vmap_info>::const_iterator i ;
  for(i=sources.begin();i!=sources.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    input += i->var ;
  }

  for(i=targets.begin();i!=targets.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    output += i->var ;
  }



  for(i=sources.begin();i!=sources.end();++i) {
    if(i->mapping.size() == 0) {
      variableSet::const_iterator vi ;
      for(vi=i->var.begin();vi!=i->var.end();++vi) {
	std::list<variable> vbasic ;

	vbasic.push_back(*vi) ;
	validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
	variableSet::const_iterator vi ;
	std::list<variable> mapping_list = maplist[j] ;
	validate_set.insert(mapping_list) ;
	for(vi=i->var.begin();vi!=i->var.end();++vi) {
	  std::list<variable> mapping_list2 = maplist[j] ;
	  mapping_list2.push_back(*vi) ;
	  validate_set.insert(mapping_list2) ;
	}
	mapping_list.pop_back() ;
	while(!mapping_list.empty()) {
	  validate_set.insert(mapping_list) ;
	  mapping_list.pop_back() ;
	}
      }
    }
  }

  for(i=targets.begin();i!=targets.end();++i) {
    if(i->mapping.size() == 0) {
      variableSet::const_iterator vi ;
      for(vi=i->var.begin();vi!=i->var.end();++vi) {
	std::list<variable> vbasic ;
	variable vt = *vi ;
	while(vt.get_info().priority.size() != 0)
	  vt = vt.drop_priority() ;
	vbasic.push_back(vt) ;
	validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
	variableSet::const_iterator vi ;
	std::list<variable> mapping_list = maplist[j] ;
	validate_set.insert(mapping_list) ;
	for(vi=i->var.begin();vi!=i->var.end();++vi) {
	  std::list<variable> mapping_list2 = maplist[j] ;
	  variable vt = *vi ;
	  while(vt.get_info().priority.size() != 0)
	    vt = vt.drop_priority() ;
	  mapping_list2.push_back(vt) ;
	  validate_set.insert(mapping_list2) ;
	}
	mapping_list.pop_back() ;
	while(!mapping_list.empty()) {
	  validate_set.insert(mapping_list) ;
	  mapping_list.pop_back() ;
	}
      }
    }
  }

  variableSet::const_iterator vi ;
  all_vars = input;
  all_vars += output ;

  for(vi=input.begin();vi!=input.end();++vi) {
    if(vi->get_info().priority.size() != 0) {
      ostringstream oss ;
      oss<< "improper use of priority annotation on rule input, var=" << *vi << endl ;
      throw parseError(oss.str()) ;
    }
  }

  if(rule_type != "pointwise" && rule_type != "default") {
    for(vi=output.begin();vi!=output.end();++vi) {
      if(vi->get_info().priority.size() != 0) {
	ostringstream oss ;
	oss << "only pointwise rules can use priority annotation, var="<< *vi << endl ;
	throw parseError(oss.str()) ;
      }
    }
  }

  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    vnames[*vi] = var2name(*vi) ;
    if(vi->get_info().priority.size() != 0) {
      variable v = *vi ;
      while(v.get_info().priority.size() != 0)
	v = v.drop_priority() ;
      vnames[v] = vnames[*vi] ;
    }
  }


  variableSet checkset ;//this is local
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    checkset += v ;
  }
  list<pair<variable,variable> >::const_iterator ipi ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
    variable v = ipi->first ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->first << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }
    v = ipi->second ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->second << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }

  }

  //local_type_map
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    if(v.is_time_variable()) {
      local_type_map[v] = pair<string,string>("param","<int> ") ;
      continue ;
    }
    v = convertVariable(v) ;
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = type_map.find(v)) == type_map.end()) {

      v = v.new_offset(0) ;
      v = v.drop_assign() ;
      while(v.time() != time_ident())
	v = v.parent() ;
      if((mi = type_map.find(v)) == type_map.end()) {
	while(v.get_info().namespac.size() != 0)
	  v = v.drop_namespace() ;
	if((mi = type_map.find(v)) == type_map.end()) {
	  string s ;
	  s = "unable to determine type of variable " ;
	  s += v.str() ;
	  throw parseError(s) ;
	}
      }
    }
    local_type_map[*vi] = mi->second ;
  }

  if(rule_type == "pointwise") {
    for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi) {
      if(local_type_map[*vi].first == "param" && vi->get_info().name != "OUTPUT") {
	throw(parseError("pointwise rule cannot compute param, use singleton")) ;
      }
    }
  }
  if(rule_type == "singleton") {
    for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi) {
      string t = local_type_map[*vi].first ;
      if(t == "store" || t == "storeVec" || t == "multiStore") {
	throw(parseError("singleton rule cannot compute store's, use pointwise")) ;
      }
    }
  }


  //singletonApply
  singletonApply = false ;
  if(rule_type == "apply") {
    if(output.size() != 1)
      throw parseError("apply rule should have only one output variable") ;
    variable av = *(output.begin()) ;
    pair<string,string> tinfo = local_type_map[av] ;

    if(tinfo.first == "param") {
      variableSet::const_iterator vi ;
      bool allparam = true ;
      for(vi=input.begin();vi!=input.end();++vi) {
	pair<string,string> tinfo2 = local_type_map[*vi] ;
	if(tinfo2.first != "param") {
	  allparam = false ;
	}
      }
      if(allparam)
	singletonApply = true ;
    }
  }

  //check output/input variables
  outs = output ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  ins = input ;
  ins -= outs ;

  cpu_only = false; //start check if cpu-only

  for(vi=ins.begin();vi!=ins.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if((mi->second).first == "multiStore" || (mi->second).first == "storeVec" || (mi->second).first == "storeMat") cpu_only = true ; //check container type
  }

  //check output_param
  output_param = false ;
  for(vi=outs.begin();vi!=outs.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if((mi->second).first == "multiStore" || (mi->second).first == "storeVec") cpu_only = true ; //check container type
    if(vi->get_info().name != "OUTPUT" && mi->second.first == "param") {
      output_param= true ;
    }
  }

  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  //check constraint?
  if(constraint!="") {
    // Check to see that the constraint is typed
    exprP C = expression::create(constraint) ;
    set<vmap_info> Cdigest ;
    fill_descriptors(Cdigest,collect_associative_op(C,OP_COMMA)) ;
    set<vmap_info>::const_iterator i ;
    variableSet constraint_vars ;
    for(i=Cdigest.begin();i!=Cdigest.end();++i) {
      for(size_t j=0;j<i->mapping.size();++j)
	constraint_vars += i->mapping[j] ;
      constraint_vars += i->var ;
    }

    for(variableSet::const_iterator vi=constraint_vars.begin();
	vi!=constraint_vars.end();++vi) {
      variable v = *vi ;
      v = convertVariable(v) ;
      map<variable,pair<string,string> >::const_iterator mi ;
      if((mi = type_map.find(v)) == type_map.end()) {
	v = v.new_offset(0) ;
	v = v.drop_assign() ;
	while(v.time() != time_ident())
	  v = v.parent() ;

	if((mi = type_map.find(v)) == type_map.end()) {
	  while(v.get_info().namespac.size() != 0)
	    v = v.drop_namespace() ;
	  mi = type_map.find(v) ;
	}
      }
      if(mi == type_map.end()  && v.get_info().name != "UNIVERSE" &&
	 v.get_info().name != "EMPTY") {
	cerr << filename << ':' << start_line+lines << ":0: warning: type of constraint variable '" << v << "' not found!"  << endl  ;

      }
    }
  }

  signature_line = start_line + lines;

  prelude_line = start_line+lines;
  use_compute = true ;
  if(use_prelude) {
    prelude.get(is) ;
    lines += prelude.num_lines();
    parsebase::killsp(is) ;
    if(is.peek() == ';') {
      is.get() ;
      use_compute = false ;
    }
    if(is_name(is)) {
      string s = get_name(is) ;
      if(s != "compute") {
	throw parseError("syntax error, expecting 'compute'") ;
      }
    }
    parsebase::killsp(is) ;
  }


  if(use_compute && is.peek() != '{')
    throw parseError("syntax error, expecting '{' place 1") ;

  compute_line = start_line+lines;

  if(use_compute){
    compute.get(is);
    lines += compute.num_lines();
  }

  last_line  = start_line+lines;

  bool sized_outputs = false;
  variableSet outsmi = outs ;
  outsmi -= input ;
  for(vi=outsmi.begin();vi!=outsmi.end();++vi) {
    string ot = local_type_map[*vi].first ;
    if(ot == "storeVec" || ot == "storeMat" || ot == "multiStore")
      sized_outputs = true ;
  }

  use_calculate = true;
  if(rule_type == "singleton" ||
     rule_type == "optional"  ||
     rule_type == "default" ||
     rule_type == "constraint" ||
     (output_param && rule_type != "apply" ) ) {
    if(use_prelude) {
      string error = "inappropriate prelude on " + rule_type + " rule." ;
      throw parseError(error) ;
    }
    use_calculate = false;
  } else {
    if(use_compute) {
      if(singletonApply) {
	cerr << "NOTE: parameter only apply rule on '" << output << "' now executes single instance." << endl ;
      }
    }
  }


  if(!use_prelude && sized_outputs && (rule_type != "apply"))
    throw parseError("need prelude to size output type!") ;


  if(rule_type == "singleton" || rule_type == "default" || rule_type == "optional"
     || rule_type == "constraint" || rule_type =="blackbox") cpu_only = true;
  if(use_prelude || is_specialized ) cpu_only = true;
  if(constraint != "") cpu_only = true;
  if( singletonApply) cpu_only = true;
  if(!inplace.empty())  cpu_only = true;
  //if multiStore or storeVec exist, aleady done

  return is ;
}

//write out cpu version of signature part, and constructor
void parserule::out_sig_cpu(std::string& filename, std::ostream &outputFile) {
  //cout<<" line " << first_line << " to " << last_line << " prelude " << prelude_line << " compute " << compute_line << endl;

  //rule name
  if(!prettyOutput)
    outputFile << "namespace {" ;

  outputFile << "class " << class_name << " : public Loci::" << rule_type << "_rule" ;

  //apply rule, the sig is a little bit longer
  if(rule_type == "apply") {
    variable av = *(output.begin()) ;
    pair<string,string> tinfo = local_type_map[av] ;
    outputFile << "< " << tinfo.first << tinfo.second <<","
	       << apply_op.str() ;

    if(tinfo.first == "storeVec") {
      outputFile << "<Vect" << tinfo.second <<" > " ;
    } else if(tinfo.first == "storeMat") {
      outputFile << "<Mat" << tinfo.second <<" > " ;
    } else {
      outputFile << tinfo.second ;
    }

    outputFile << "> " ;
  }

  outputFile << " {" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;

  variableSet::const_iterator vi ;

  //input variables
  for(vi=ins.begin();vi!=ins.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }


    if(!prettyOutput)
      outputFile << "    Loci::const_" << mi->second.first <<  mi->second.second ;
    else
      outputFile << "    const_" << mi->second.first <<  mi->second.second ;

    outputFile << " " << vnames[*vi] << " ; " << endl ;
    ::syncFile(signature_line,filename, outputFile) ;
  }

  //output variables
  for(vi=outs.begin();vi!=outs.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(!prettyOutput)
      outputFile << "    Loci::" << mi->second.first <<  mi->second.second ;
    else
      outputFile << "    " << mi->second.first <<  mi->second.second ;

    outputFile << " " << vnames[*vi] << " ; " << endl ;

    ::syncFile(signature_line,filename, outputFile) ;

  }



  outputFile << "public:" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;

  outputFile <<   "    " << class_name << "() {" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;


  //name_store part,
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
	       << vnames[*vi] << ") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }

  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }
  list<pair<variable,variable> >::const_iterator ipi;
  for(set<vmap_info>::const_iterator  i=targets.begin();i!=targets.end();++i) {
    outputFile <<   "       output(\"" ;
    for(size_t j=0;j<i->mapping.size();++j)
      outputFile << i->mapping[j] << "->" ;

    // Output target variables, adding inplace notation if needed
    variableSet::const_iterator vi ;
    if(i->var.size() > 1)
      outputFile << '(' ;
    for(vi=i->var.begin();vi!=i->var.end();++vi) {
      if(vi != i->var.begin())
	outputFile << ',' ;
      for( ipi=inplace.begin();ipi!=inplace.end();++ipi) {
	if((ipi->first) == *vi)
	  break ;
      }
      if(ipi!=inplace.end()) {
	if(i->mapping.size() == 0 || i->var.size() > 1)
	  outputFile << ipi->first << "=" << ipi->second ;
	else
	  outputFile << '('<<ipi->first << "=" << ipi->second <<')';
      } else
	outputFile << *vi ;
    }
    if(i->var.size() > 1)
      outputFile << ')' ;

    outputFile <<  "\") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }

  //output constraints,  parametric_var , set_specialized, conditional, options and comments
  if(constraint!="") {
    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }

  if(parametric_var != "") {
    outputFile <<   "       set_parametric_variable(\""
	       << parametric_var << "\") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }
  if(is_specialized) {
    outputFile <<   "       set_specialized() ; " << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }
  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }
  list<string>::const_iterator lsi ;
  for(lsi=options.begin();lsi!=options.end();++lsi) {
    string s = *lsi ;
    bool has_paren = false ;
    for(size_t i = 0;i<s.size();++i)
      if(s[i] == '(')
	has_paren = true ;
    outputFile <<   "       " << s ;
    if(!has_paren)
      outputFile << "()" ;
    outputFile << " ;" << endl;
    ::syncFile(signature_line,filename, outputFile) ;

  }
  for(lsi=comments.begin();lsi!=comments.end();++lsi) {
    outputFile <<   "       comments(" << *lsi << ") ;" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;

  }

  outputFile <<   "    }" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;
}


//write out cuda version of signature part, including constructors and bind(), setDomain() and getDomain() functions
void parserule::out_sig_cuda(std::string& filename, std::ostream &outputFile) {

  if(!prettyOutput)
    outputFile << "namespace {" ;

  outputFile << "class " << class_name << " : public " << cuda_type(rule_type) ;

  if(rule_type == "apply") {

    variable av = *(output.begin()) ;
    pair<string,string> tinfo = local_type_map[av] ;

    //will Cuda support apply_op?
    outputFile << "< " << tinfo.first << tinfo.second <<","
	       << apply_op.str() ;

    if(tinfo.first == "storeVec") {
      outputFile << "<Vect" << tinfo.second <<" > " ;
    } else if(tinfo.first == "storeMat") {
      outputFile << "<Mat" << tinfo.second <<" > " ;
    } else {
      outputFile << tinfo.second ;
    }

    outputFile << "> " ;
  }
  outputFile << " {" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;

  //added constructors
  outputFile << "public:" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;

  outputFile <<   "    " << class_name << "() {}" << endl ;
  outputFile <<   "    " << class_name << "(int ctxId)" <<endl;
  outputFile <<   "    " <<":" << cuda_type(rule_type)<<"(ctxId)  {}" << endl ;
  ::syncFile(signature_line,filename, outputFile) ;
  outputFile << endl << endl;
  outputFile << "    "  << "class Compute { " << endl;
  outputFile << "      " << "public: " << endl;
  outputFile << "      " << "Loci::sequence domain ;" << endl;

  //main changes here, the variables specify as pointers

  variableSet::const_iterator vi ;
  for(vi=ins.begin();vi!=ins.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }

    outputFile << "      const " << var2ptr(mi->second) ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;

    if(mi->second.first.find("multi") != mi->second.first.npos){
      outputFile << "      const Loci::int_type*" ;
      outputFile << " " << vnames[*vi] << "Offset  ; " << endl ;
    }
    ::syncFile(signature_line,filename, outputFile) ;
  }


  for(vi=outs.begin();vi!=outs.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }

    outputFile << "      " << var2ptr(mi->second) ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;

    if(mi->second.first.find("multi") != mi->second.first.npos){
      outputFile << "    Loci::int_type" ;
      outputFile << " " << vnames[*vi] << "Offset  ; " << endl ;
    }
    ::syncFile(signature_line,filename, outputFile) ;
  }


  //cuda output bind function instead of name_store
  //domain need to be set separately
  //so add setDomain() and getDomain() function

  outputFile << endl << "      void bind(StoreDB<GI, T> & db) {" << endl;
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       " << vnames[*vi] <<" =  db."<< *vi << ";" << endl ;
    ::syncFile(signature_line,filename, outputFile) ;
  }

  outputFile << "      }" << endl<<endl;
  ::syncFile(signature_line,filename, outputFile) ;
  outputFile << "      __host__ __device__" << endl;
  outputFile << "      Loci::sequence getDomain() const {" << endl;
  outputFile << "        return domain ;" << endl;
  outputFile << "      }" << endl << endl;
  outputFile << "      __host__ __device__" << endl;
  outputFile << "      void setDomain(Loci::sequence& dom)  {" << endl;
  outputFile << "       domain = dom ;" << endl; //is it copy here?
  outputFile << "      }" << endl << endl;
  ::syncFile(signature_line,filename, outputFile) ;
}

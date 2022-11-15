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



#ifndef PARSERULE_H
#define PARSERULE_H

#include "parsebase.h"
#include "variable.h"

#include <ctype.h>
#include <set>
#include <iostream>
#include <sstream>
#include <sys/timeb.h>
#include <cstring>

using std::istringstream ;
using std::ostringstream ;

using std::pair ;
using std::list ;
using std::string ;
using std::set ;
using std::map ;
using std::vector ;
using std::istream ;
using std::ifstream ;
using std::ofstream ;
using std::ostream ;
using std::ios ;
using std::endl ;
using std::cerr ;
using std::cout ;
using namespace Loci ;



//this class is added to parse a rule as an object
class parserule : public parsebase {
 public:
  string rule_type ;
  nestedparenstuff signature ;
  nestedbracketstuff apply_op ;
  string constraint, conditional ;
  string parametric_var ;
  list<string> options ;
  list<string> comments ;
  list<pair<variable,variable> > inplace ;

  string sig ;
  string heads,bodys ;
  exprP head=0,body=0 ;

  string class_name ;
  set<vmap_info> sources ;
  set<vmap_info> targets ;
  variableSet input,output ;
  variableSet ins, outs;

  set<std::list<variable> > validate_set ;
  variableSet all_vars;

  map<variable,pair<string,string> > local_type_map ;
  map<variable,string> vnames ;

  nestedbracestuff prelude;
  nestedbracestuff compute;

  bool use_prelude  ;
  bool is_specialized ;
  bool cpu_only ;
  bool use_compute;
  bool output_param  ;
  bool singletonApply;
  bool use_calculate;
  int first_line, signature_line, prelude_line, compute_line, last_line; //for syncFile

  istream &get(std::string& filename, int& cnt, int start_line, const std::map<Loci::variable,std::pair<std::string,std::string> > &type_map, istream &is) ;
  
  //write out cpu version of signature part, and constructor
  void out_sig_cpu(std::string& filename, std::ostream &outputFile);

  //write out cuda version of signature part, including constructors and bind(), setDomain() and getDomain() functions
  void out_sig_cuda(std::string& filename, std::ostream &outputFile);

  int num_lines() {
    return lines ;
  }
} ;





#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
template<class T> class templlist : public parsebase {
 public:
  list<T> flist ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '<')
      return s ;
    char c ;
    s.get(c);
    parsebase::killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      parsebase::killsp(s) ;
      if(s.peek() == '>') {
	s.get() ;
	return s ;
      }
      if(s.peek() != ',') {
	throw parseError("syntax error, expected comma") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "<" ;
      typename list<T>::const_iterator ii ;
      ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
	s+="," ;
	s+= ii->str() ;
      }
      s += "> " ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    typename list<T>::const_iterator ii ;
    for(ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;



#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
template<class T> class funclist : public parsebase {
 public:
  list<T> flist ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '(')
      return s ;
    char c = s.get();
    parsebase::killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      parsebase::killsp(s) ;
      if(s.peek() == ')') {
	c = s.get() ;
	return s ;
      }
      if(s.peek() != ',') {
	throw parseError("syntax error") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "(" ;
      typename list<T>::const_iterator ii ;
      ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
	s+="," ;
	s+= ii->str() ;
      }
      s += ")" ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    typename list<T>::const_iterator ii ;
    for(ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;

class typestuff : public parsebase {
 public:
  string name ;
  templlist<typestuff> templ_args ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(isalpha(s.peek()) || s.peek() == '_') {
      char c = s.peek() ;
      while(isalpha(c = s.peek()) || isdigit(c) || c == '_' ||  c == ':')
	name += s.get() ;
    } else if(isdigit(s.peek())) {
      while(isdigit(s.peek())) {
	name += s.get() ;
      }
    } else
      throw parseError("syntax error") ;
    templ_args.get(s) ;
    return s ;
  }
  string str() const {
    string s ;
    s+= name ;
    s+= templ_args.str() ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i+= templ_args.num_lines() ;
    return i ;
  }
} ;

variable convertVariable(variable v);


#endif

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
#include "lpp.h"
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


bool is_name(istream &s) {
  int ch = s.peek() ;
  return isalpha(ch) || ch == '_' ;
}

string get_name(istream &s) {
  if(!is_name(s))
    throw parseError("expected name ") ;

  string str ;
  while(!s.eof() && (s.peek() != EOF) &&
        (isalnum(s.peek()) || (s.peek() == '_')) )
    str += s.get() ;

  return str ;
}

bool is_string(istream &s) {
  return s.peek() == '\"' ;
}

string get_string(istream &s) {
  if(!is_string(s))
    throw parseError("expected string") ;
  string str ;
  if(s.eof())
    throw parseError("unexpected EOF") ;
  s.get() ;
  int ch = s.get() ;
  while(ch != '\"' &&!s.eof()) {
    str += ch ;
    ch = s.get() ;
  }
  if(ch!='\"')
    throw parseError("no closing \" for string") ;
  return str ;
}

//namespace definition, not "using namespace...".
//maybe also killcomments somewhere
bool is_namespace(istream &s) {
  bool result = false;
  if(is_name(s)){
    string ns_str = get_name(s);
    if(ns_str == "namespace"){
      string str;
      int ch = s.get() ;
      while(ch != '{' && ch != ';' && !s.eof()) {
        str += ch ;
        ch = s.get() ;
      }
      str += ch ;
      if(ch =='{')result = true;
      for (std::string::reverse_iterator rit=str.rbegin(); rit!=str.rend(); ++rit) s.putback(*rit);
    }
    for (std::string::reverse_iterator rit=ns_str.rbegin(); rit!=ns_str.rend(); ++rit) s.putback(*rit);
  }
  return result ;
}

//multi-ostream version, keep namespace definition inside non-cpu output
istream &killNamespOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
  if(!is_namespace(s))
    throw parseError("expected namespace") ;
  string str ;
  if(s.eof())
    throw parseError("unexpected EOF in killNamespOuts") ;
  int ch = s.get() ;
  while(ch != '{' &&!s.eof()) {
    if(ch == '\n')lines++;
    str += ch ;
    ch = s.get() ;
  }
  str += ch ;
  if(ch!='{')
    throw parseError("no beginning { for namespace") ;

  if(oneOutput){
    *(out[0].second) << str ;
  }else{
    for(size_t i = 0; i < out.size(); i++){
      *(out[i].second) << str ;
    }
  }
  return s;
}


bool is_include(istream &s) {
  bool result = false;
  if(s.peek()=='#'){
    string str;
    str += s.get();
    if(is_name(s)){
      str += get_name(s);
      if(str == "#include"){
        result = true;
      }
    }
    for (std::string::reverse_iterator rit=str.rbegin(); rit!=str.rend(); ++rit) s.putback(*rit);
  }
  return result ;
}
//multi-ostream version, keep "#include ..."  inside non-cpu output
istream &killIncludeOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
  if(!is_include(s))
    throw parseError("expected include") ;
  string str ;
  if(s.eof())
    throw parseError("unexpected EOF in killIncludeOuts") ;
  int ch = s.get() ;
  while(ch != '\n' &&!s.eof()) {
    str += ch ;
    ch = s.get() ;
  }
  str += ch;
  lines++;

  if(ch!='\n')
    throw parseError("no newline for include") ;


  if(oneOutput){
    *(out[0].second) << str ;
  }else{
    for(size_t i = 0; i < out.size(); i++){
      *(out[i].second) << str ;
    }
  }
  return s;
}



bool is_comment(istream &s) {
  if(s.peek() != '/')
    return false ;

  s.get() ;
  char c = s.peek() ;
  s.unget() ;
  if(c == '/' || c == '*')
    return true ;
  return false ;
}



istream &killComment(istream &s, int & lines) {
  s.get() ;
  char c = s.get()  ;
  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n') {
      s.get() ;
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
    }
    return s ;
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}

istream &killsp(istream &s, int &lines) {

  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      s.get();
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killComment(s,lines) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  return s ;
}

//return line count
int kill_sp(istream &s, int &lines) {
  int l = lines;
  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      s.get();
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killComment(s,lines) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;

  return lines-l ;
}

//multi-ostream version, it will be odd if a rule is cpu-only but the comments outside it
//is written out into non-cpu version, but currently we can not decide in advance
istream &killCommentOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
  s.get() ;

  if(oneOutput){
    *(out[0].second) << '/' ;
  }else{
    for(size_t i = 0; i < out.size(); i++){
      *(out[i].second) << '/' ;
    }
  }

  char c = s.get()  ;
  if(oneOutput){
    *(out[0].second) << c ;
  }else{
    for(size_t i = 0; i < out.size(); i++){
      *(out[i].second) << c ;
    }
  }

  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n') {
      char c = s.get() ;
      if(oneOutput){
        *(out[0].second) << c ;
      }else{
        for(size_t i = 0; i < out.size(); i++){
          *(out[i].second) << c ;
        }
      }
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
      if(oneOutput){
        *(out[0].second) << '\n' ;
      }else{
        for(size_t i = 0; i < out.size(); i++){
          *(out[i].second) << '\n' ;
        }
      }
    }
    return s ;
  }

  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    if(oneOutput){
      *(out[0].second) << c ;
    }else{
      for(size_t i = 0; i < out.size(); i++){
        *(out[i].second) << c ;
      }
    }
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        if(oneOutput){
          *(out[0].second) << '/' ;
        }else{
          for(size_t i = 0; i < out.size(); i++){
            *(out[i].second) << '/' ;
          }
        }
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}


istream &killCommentOut(istream &s, int & lines,ostream &out) {
  s.get() ;
  out << '/' ;
  char c = s.get()  ;
  out << c ;
  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n') {
      char c = s.get() ;
      out << c ;
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
      out << '\n' ;
    }
    return s ;
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    out << c ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        out << '/' ;
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}

istream &killspOut(istream &s, int &lines, ostream &out) {

  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      char c = s.get();
      out << c ;
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killCommentOut(s,lines,out) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  return s ;
}

// added this function, for nestedbracestuff
istream &killspOstr(istream &s, int &lines, string &str) {
  std::ostringstream out;
  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      char c = s.get();
      out << c ;
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killCommentOut(s,lines,out) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  str += out.str();
  return s ;
}





void syncFile(int line_no, const string& filename, std::ostream &outputFile) {
  if(!prettyOutput)
    outputFile << "#line " << line_no << " \"" << filename << "\"" << std::endl ;
}

int parseFile::killsp() {
  int l = line_no ;
  ::killsp(is,line_no) ;
  return l-line_no ;
}

int parseFile::killspout(std::ostream &outputFile) {
  int l = line_no ;
  ::killspOut(is,line_no,outputFile) ;
  return l-line_no ;
}

class parsebase {
public:
  int lines ;
  parsebase() {lines = 0 ; }
  istream &killsp(istream &s) {
    ::killsp(s,lines) ;
    return s ;
  }
  istream &killspOstr(istream &s, string& str) {
    ::killspOstr(s,lines, str) ;
    return s ;
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

#ifdef __GNUC__
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


class bracestuff : public parsebase {
public:
  string stuff ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() == '{') {
      char c = s.get() ;
      while(s.peek() != EOF && s.peek() != '}') {
        c = s.get() ;
        if(c == '{')
          throw parseError("syntax error") ;
        stuff += c ;
        parsebase::killsp(s) ;
      }
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      c = s.get() ;
      parsebase::killsp(s) ;
    }
    return s ;
  }

  string str() const {
    string s ;
    if(stuff == "")
      return s ;
    s += "{" ;
    s += stuff ;
    s += "}" ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    return i ;
  }
} ;


class var : public parsebase {
public:
  bool isdollar ;
  string name ;
  list<string> prio_list ;
  list<string> nspace_list ;
  funclist<var> param_args ;
  bracestuff bs ;
  var() : isdollar(false) {}

  istream &get(istream &s) {
    isdollar = false ;
    parsebase::killsp(s) ;
    if(s.peek() == '$') {
      s.get() ;
      isdollar=true ;
    }
    if(!is_name(s))
      throw parseError("syntax error: expecting name after '$'") ;
    name = get_name(s) ;
    parsebase::killsp(s) ;
    if(s.peek() == ':') {
      while(s.peek() == ':') {
        s.get() ;
        if(s.peek() != ':') {
	  string err = "syntax error, improper trailing colon, use parenthesis around variable '"+ name+"' to fix." ;
	  throw parseError(err.c_str()) ;
	}
        s.get() ;
        parsebase::killsp(s) ;
        prio_list.push_back(name);
        if(!is_name(s))
          throw parseError("syntax error near ':'") ;
        name = get_name(s) ;
        parsebase::killsp(s) ;
      }
    }
    if(s.peek() == '@') {
      while(s.peek() == '@') {
        s.get() ;
        parsebase::killsp(s) ;
        nspace_list.push_back(name);
        if(!is_name(s))
          throw parseError("syntax error near '@'") ;
        name = get_name(s) ;
        parsebase::killsp(s) ;
      }
    }

    param_args.get(s) ;
    bs.get(s) ;

    return s ;
  }
  string str() const {
    string s ;
    list<string>::const_iterator li ;
    for(li=prio_list.begin();li!=prio_list.end();++li)
      s+= *li + "::" ;
    if(isdollar)
      s+="$" ;
    for(li=nspace_list.begin();li!=nspace_list.end();++li)
      s += *li + "@" ;
    s+=name ;
    s+= param_args.str() ;
    s+= bs.str() ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i += param_args.num_lines() ;
    i += bs.num_lines() ;
    return i ;
  }
} ;

class nestedparenstuff : public parsebase {
public:
  string paren_contents ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '(')
      throw parseError("syntax error, expecting '('") ;
    s.get() ;
    int open_parens = 0 ;
    parsebase::killsp(s) ;
    while(s.peek() != ')' || open_parens != 0) {
      if(s.peek() == '"') { // grab string
	paren_contents += s.get() ;
	while(s.peek() != '"') {
	  if(s.peek() == EOF) {
	    throw parseError("unexpected EOF parsing sting") ;
	  }
	  if(s.peek() == '\n' || s.peek() == '\r') {
	    lines++ ;
	  }
	  paren_contents += s.get() ;
	}
	paren_contents += s.get() ;
	continue ;
      }
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '(')
        open_parens++ ;
      if(s.peek() == ')')
        open_parens-- ;
      if(s.peek() == '\n' || s.peek() == '\r') {
        s.get() ;
        lines++ ;
        continue ;
      }
      paren_contents += s.get() ;
      parsebase::killsp(s) ;
    }
    s.get() ;
    parsebase::killsp(s) ;
    return s ;
  }
  string str() {
    return paren_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

class nestedbracketstuff : public parsebase {
public:
  string bracket_contents ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '[')
      throw parseError("syntax error, expecting '['") ;
    s.get() ;
    parsebase::killsp(s) ;
    int open_brackets = 0 ;
    while(s.peek() != ']' || open_brackets != 0) {
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '[')
        open_brackets++ ;
      if(s.peek() == ']')
        open_brackets-- ;
      if(s.peek() == '\n' || s.peek() == '\r') {
        s.get() ;
        lines++ ;
        continue ;
      }
      bracket_contents += s.get() ;
      parsebase::killsp(s) ;
    }
    s.get() ;
    return s ;
  }
  string str() {
    return bracket_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

//read in everything between { and } , no modification
//this class is added to read in compute part of a rule
class nestedbracestuff : public parsebase {
public:
  string brace_contents ;
  istream &get(istream &s) {

    parsebase::killspOstr(s, brace_contents) ;
    if(s.peek() != '{')
      throw parseError(string("syntax error, expecting '{' place 2")) ;
    s.get() ;

    parsebase::killspOstr(s, brace_contents) ;
    int open_braces = 0 ;
    while(s.peek() != '}' || open_braces != 0) {
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '{')
        open_braces++ ;
      if(s.peek() == '}')
        open_braces-- ;
      if(s.peek() == '\n' || s.peek() == '\r') {
        brace_contents += s.get() ;
        lines++ ;
        continue;
      }
      brace_contents += s.get() ;
      parsebase::killspOstr(s, brace_contents) ;
    }
    brace_contents += s.get() ;
    return s ;
  }
  string str() {
    return brace_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;





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

/*collect the variables in type_map
  output original as comments and newlines
  no difference between cpu and non-cpu version
*/
void parseFile::setup_Type(std::vector<pair<string, std::ostream*> > &outputFile, bool oneOutput) {
  var vin ;
  vin.get(is) ;
  typestuff tin ;
  tin.get(is) ;
  while(is.peek() == ' ' || is.peek() == '\t')
    is.get() ;
  if(is.peek() != ';')
    throw parseError("syntax error, missing ';'") ;
  is.get() ;
  variable v(vin.str()) ;
  v = convertVariable(v) ;
  if(oneOutput){
    *(outputFile[0].second) << "// $type " << v << ' ' << tin.str() ;
  }else{
    for(size_t j = 0; j < outputFile.size(); j++){
      *(outputFile[j].second) << "// $type " << v << ' ' << tin.str() ;
    }
  }
  int nl = vin.num_lines()+tin.num_lines() ;
  line_no += nl ;
  for(int i=0;i<nl;++i){
    if(oneOutput){
      *(outputFile[0].second) << endl ;
    }else{
      for(size_t j = 0; j < outputFile.size(); j++){
        *(outputFile[j].second) << endl ;
      }
    }
  }
  type_map[v] = pair<string,string>(tin.name,tin.templ_args.str()) ;
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

class specialCommand{
public:
  string name;
  string body;
};


void  parseFile::process_SpecialCommand(istream  &str_is,
                                        int& lines,
                                        std::ostream &outputFile,
                                        const map<variable,string> &vnames,
                                        int &openbrace) {

  str_is.get() ; // get leading [
  string name = get_name(str_is) ;
  if(str_is.peek() != ']') {
    cerr << "expecting ']' to close special command '" << name << "'" << endl ;
    throw parseError("syntax error") ;
  }
  str_is.get() ;

  int nsz = name.size() ;
  for(int i=0;i<nsz;++i)
    if(name[i] >= 'A' || name[i] <= 'Z')
      name[i] = std::tolower(name[i]) ;

  if(name == "once") {
    ::killsp(str_is, lines) ;
    if(str_is.peek() != '{') {
      cerr << "expecting '{' after $[Once] command" << endl ;
      cerr << "found " << char(str_is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    outputFile << "if(Loci::is_leading_execution()) " ;

  } else if(name == "atomic") {
    ::killsp(str_is, lines) ;
    if(str_is.peek() != '{') {
      cerr << "expecting '{' after $[Atomic] command" << endl ;
      cerr << "found " << char(str_is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    str_is.get() ;
    openbrace++ ;
    outputFile << "{ Loci::atomic_region_helper L__ATOMIC_REGION ; " << endl ;
  } else {
    cerr << "unknown special command '[" << name << "]' !" << endl ;
    throw parseError("syntax error") ;
  }
}

//multi-ostream version, don't know if non-cpu version should have special command
//keep it for later
void  parseFile::process_SpecialCommands(istream  &str_is,
                                         int& lines,
                                         std::vector<pair<std::string, std::ostream*> > &outputFile,
                                         bool oneOutput,
                                         const map<variable,string> &vnames,
                                         int &openbrace) {

  str_is.get() ; // get leading [
  string name = get_name(str_is) ;
  if(str_is.peek() != ']') {
    cerr << "expecting ']' to close special command '" << name << "'" << endl ;
    throw parseError("syntax error") ;
  }
  str_is.get() ;

  int nsz = name.size() ;
  for(int i=0;i<nsz;++i)
    if(name[i] >= 'A' || name[i] <= 'Z')
      name[i] = std::tolower(name[i]) ;

  if(name == "once") {
    ::killsp(str_is, lines) ;
    if(str_is.peek() != '{') {
      cerr << "expecting '{' after $[Once] command" << endl ;
      cerr << "found " << char(str_is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }

    if(oneOutput){
      *(outputFile[0].second) << "if(Loci::is_leading_execution()) " ;
    }else{
      for(size_t i = 0; i < outputFile.size(); i++){
        if(outputFile[i].first == "cpu")
          *(outputFile[i].second) << "if(Loci::is_leading_execution()) " ;
      }
    }

  } else if(name == "atomic") {
    ::killsp(str_is, lines) ;
    if(str_is.peek() != '{') {
      cerr << "expecting '{' after $[Atomic] command" << endl ;
      cerr << "found " << char(str_is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    str_is.get() ;
    openbrace++ ;
    if(oneOutput){
      *(outputFile[0].second) << "{ Loci::atomic_region_helper L__ATOMIC_REGION ; " << endl ;
    }else{
      for(size_t i = 0; i < outputFile.size(); i++){
        if(outputFile[i].first == "cpu")
          *(outputFile[i].second) << "{ Loci::atomic_region_helper L__ATOMIC_REGION ; " << endl ;
      }
    }
  } else {
    cerr << "unknown special command '[" << name << "]' !" << endl ;
    throw parseError("syntax error") ;
  }
}

void parseFile::process_Prelude(std::string& in,
                                std::ostream &outputFile,
                                const map<variable,string> &vnames) {

  //modified, instead of processing from istream of parseFile, process from a string in

  istringstream str_is(in) ;//redefine istream
  int lines = 0;

  outputFile << "    virtual void prelude(const Loci::sequence &seq) { " ;

  int openbrace = 1 ;
  for(;;) {
    ::killspOut(str_is,lines,outputFile) ;

    if(str_is.peek() == EOF){
      throw parseError("unexpected EOF") ;
    }
    if(str_is.peek() == '}') {
      str_is.get() ;
      outputFile << '}' ;

      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(str_is.peek() == '{') {
      str_is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(str_is.peek() == '$') {
      string name ;
      variable v ;
      str_is.get() ;
      if(str_is.peek() == '[') {
        process_SpecialCommand(str_is, lines, outputFile,vnames,openbrace) ;
        continue ;
      }
      var vin ;
      vin.get(str_is) ;
      v = variable(vin.str()) ;

      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      outputFile << vmi->second  ;
    }

    char c = str_is.get() ;
    if(c == '\n')
      lines++ ;
    outputFile << c ;
  } ;
}

void parseFile::process_Compute(std::string& in, std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  //modified, instead of processing from istream of parseFile, process from a string in

  istringstream str_is(in) ; //redefine istream
  int lines = 0;

  outputFile << "    void compute(const Loci::sequence &seq) { " ;

  int openbrace = 1 ;
  for(;;) {
    ::killspOut(str_is, lines, outputFile) ;
    if(str_is.peek() == EOF)
      throw parseError("unexpected EOF") ;

    if(str_is.peek() == '}') {
      str_is.get() ;
      outputFile << '}' ;

      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(str_is.peek() == '{') {
      str_is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(str_is.peek() == '"') {
      str_is.get() ;
      outputFile << '"' ;
      while(str_is.peek() != '"' && str_is.peek() != EOF) {
        char c = str_is.get() ;
        outputFile << c ;
      }
      str_is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(str_is.peek() == '\'') {
      str_is.get() ;
      outputFile << '\'' ;
      while(str_is.peek() != '\'' && str_is.peek() != EOF) {
        char c = str_is.get() ;
        outputFile << c ;
      }
      str_is.get() ;
      outputFile << '\'' ;
      continue ;
    }
    if(str_is.peek() == '$') {
      variable v ;
      str_is.get() ;
      if(str_is.peek() == '[') {
        process_SpecialCommand(str_is, lines, outputFile,vnames,openbrace) ;
        continue ;
      }
      bool deref = true ;
      if(str_is.peek() == '*') {
        str_is.get() ;
        deref = false ;
      }


      var vin ;
      vin.get(str_is) ;
      v = variable(vin.str()) ;

      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " str_is unknown to thstr_is rule!" << endl ;
        throw parseError("type error") ;
      }
      map<variable,pair<string,string> >::const_iterator mi ;
      if(((mi = type_map.find(v)) != type_map.end()) &&
         (mi->second.first == "Constraint" || !deref)) {
        outputFile << vmi->second  ;
      } else {
        outputFile << "(*" << vmi->second << ')' ;
      }

    }
    char c = str_is.get() ;
    if(c == '\n')
      lines++ ;
    outputFile << c ;
  }
}



string getNumber(std::istream &is) {
  string num ;
  while(isdigit(is.peek()))
    num+= is.get();
  if(is.peek() == '.') {
    num += is.get() ;
    while(isdigit(is.peek()))
      num+= is.get();
  }
  if(is.peek() == 'e' || is.peek() == 'E') {
    num += is.get() ;
    if(is.peek() == '-' || is.peek() == '+')
      num += is.get() ;
    while(isdigit(is.peek()))
      num += is.get() ;
  }
  return num ;
}

string parseFile::process_String(string in,
                                 const map<variable,string> &vnames,
				 const set<list<variable> > &validate_set) {
  ostringstream outputFile ;
  istringstream is(in) ;

  int line_no = 0 ;

  for(;;) {
    ::killspOut(is,line_no,outputFile) ;

    if(is.peek() == EOF)
      break ;

    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      continue ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      continue ;
    }
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
      continue ;
    }
    if(is.peek() == '/') {
      is.get() ;
      outputFile << '/' ;
      if(is.peek() == '/') { // comment line
        is.get() ;
        outputFile << '/' ;
        while(is.peek() != '\n') {
          char c = is.get() ;
          outputFile << c ;
        }
        ::killspOut(is,line_no,outputFile) ;
      }
      continue ;
    }

    if(is.peek() == '#') {
      is.get() ;
      outputFile << '#' ;
      while(is.peek() != '\n') {
        char c = is.get() ;
        outputFile << c ;
      }
      ::killspOut(is,line_no,outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }

    if(is_name(is) || is.peek() == '$') {
      bool first_name = is_name(is) ;
      string name ;
      variable v ;
      string brackets ;
      if(first_name)
        name = get_name(is) ;
      else {
        is.get() ;
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;

          variable v(vin.str()) ;
          map<variable,string>::const_iterator vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }

          outputFile << vmi->second ;
          continue ;
        }

        var vin ;
        vin.get(is) ;
        v = variable(vin.str()) ;
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames,validate_set) ;
          brackets = "[" + binfo + "]" ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          ::killspOut(is,line_no,outputFile) ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            ::killspOut(is,line_no,outputFile) ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames,validate_set) ;
              brk = "[" + binfo +"]";
            }
            blist.push_back(brk) ;
          } else {
            dangling_arrow = true ;
            break ;
          }
        } else {
          is.unget() ;
          break ;
        }
      }
      if(dangling_arrow && first_name) {
        outputFile << name << " ->" ;
        continue ;
      }
      if(dangling_arrow)
        throw parseError("syntax error, near '->' operator") ;

      validate_VariableAccess(v,vlist,first_name,vnames,validate_set) ;

      if(first_name && (vlist.size() == 0)) {
        outputFile << name << ' ' ;
        continue ;
      }
      list<variable>::reverse_iterator ri ;
      for(ri=vlist.rbegin();ri!=vlist.rend();++ri) {
        map<variable,string>::const_iterator vmi = vnames.find(*ri) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << '[' ;
      }
      if(first_name) {
        outputFile << '*' << name ;
      } else {
        map<variable,string>::const_iterator vmi = vnames.find(v) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << v << " is unknown to this rule!" << endl ;
          throw parseError("type error: is this variable in the rule signature?") ;
        }
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      list<string>::const_iterator rbi ;
      for(rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

    }
    if(is.peek() != EOF) {
      char c = is.get() ;
      outputFile << c ;
    }
  }


  return outputFile.str() ;
}


void parseFile::validate_VariableAccess(variable v, const list<variable> &vlist,
					bool first_name,
					const map<variable,string> &vnames,
					const set<list<variable> > &validate_set) {

  list<variable> vlistall ;
  list<variable>::const_iterator vitmp ;
  for(vitmp=vlist.begin();vitmp!=vlist.end();++vitmp) {
    variable vt = *vitmp ;
    while(vt.get_info().priority.size() != 0)
      vt = vt.drop_priority() ;
    vlistall.push_back(vt) ;
  }
  variable vt = v ;
  while(vt.get_info().priority.size() != 0)
    vt = vt.drop_priority() ;
  vlistall.push_front(vt) ;

  if(!first_name && !vlistall.empty()
     && validate_set.find(vlistall) == validate_set.end()) {
    ostringstream msg ;
    msg << "variable access " ;
    list<variable>::const_iterator lvi ;
    for(lvi=vlistall.begin();lvi!=vlistall.end();) {
      msg << *lvi ;
      ++lvi ;
      if(lvi!=vlistall.end())
	msg << "->" ;
    }
    msg << " not consistent with rule signature!" ;
    throw parseError(msg.str()) ;
  }

  list<variable>::const_reverse_iterator ri ;
  for(ri=vlist.rbegin();ri!=vlist.rend();++ri) {
    map<variable,string>::const_iterator vmi = vnames.find(*ri) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
      throw parseError("type error") ;
    }
  }
  if(!first_name) {
    map<variable,string>::const_iterator vmi = vnames.find(v) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << v << " is unknown to this rule!" << endl ;
      throw parseError("type error: is this variable in the rule signature?") ;
    }
  }
}

void parseFile::process_Calculate(int start_line,
                                  std::string& in,
                                  std::ostream &outputFile,
                                  const map<variable,string> &vnames,
                                  const set<list<variable> > &validate_set) {
  //modified. instead of processing from istream of parseFile, process from a string in
  //start_line is for syncFile

  istringstream str_is(in) ; //redefine istream
  int lines = 0;

  ::syncFile(start_line+lines,filename, outputFile) ;

  while(str_is.peek() == ' ' || str_is.peek() == '\t')
    str_is.get() ;
  if(str_is.peek() == '\n') {
    str_is.get() ;
    lines++ ;
  }

  ::killspOut(str_is, lines, outputFile) ;
  int openbrace = 1 ;
  for(;;) {
    ::killspOut(str_is, lines, outputFile) ;
    if(str_is.peek() == EOF)
      throw parseError("unexpected EOF in process_Calculate") ;

    if(str_is.peek() == '}') {
      str_is.get() ;
      outputFile << '}' ;

      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(str_is.peek() == '{') {
      str_is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(str_is.peek() == '"') {
      str_is.get() ;
      outputFile << '"' ;
      while(str_is.peek() != '"' && str_is.peek() != EOF) {
        char c = str_is.get() ;
        outputFile << c ;
      }
      str_is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(str_is.peek() == '\'') {
      str_is.get() ;
      outputFile << '\'' ;
      while(str_is.peek() != '\'' && str_is.peek() != EOF) {
        char c = str_is.get() ;
        outputFile << c ;
      }
      str_is.get() ;
      outputFile << '\'' ;
      continue ;
    }
    if(str_is.peek() == '/') {
      str_is.get() ;
      outputFile << '/' ;
      if(str_is.peek() == '/') { // comment line
        str_is.get() ;
        outputFile << '/' ;
        while(str_is.peek() != '\n') {
          char c = str_is.get() ;
          outputFile << c ;
        }
        ::killspOut(str_is, lines, outputFile) ;
      }
      continue ;
    }

    if(str_is.peek() == '#') {
      str_is.get() ;
      outputFile << '#' ;
      while(str_is.peek() != '\n') {
        char c = str_is.get() ;
        outputFile << c ;
      }
      ::killspOut(str_is, lines, outputFile) ;
      continue ;
    }

    if(isdigit(str_is.peek())) {
      outputFile << getNumber(str_is) ;
      continue ;
    }

    if(is_name(str_is) || str_is.peek() == '$') {
      int lcount = 0 ;
      bool first_name = is_name(str_is) ;
      if(!first_name) {
        str_is.get() ;
        if(str_is.peek() == '[') { // special command
          process_SpecialCommand(str_is, lines, outputFile,vnames,openbrace) ;
          continue ;
        }
      }
      string name ;
      variable v ;
      string brackets ;
      if(first_name)
        name = get_name(str_is) ;
      else {
        if(str_is.peek() == '*') {
          str_is.get() ;
          var vin ;
          vin.get(str_is) ;
          lines += vin.num_lines() ;
          lcount += vin.num_lines() ;

          variable v(vin.str()) ;
          map<variable,string>::const_iterator vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }

          outputFile << vmi->second ;
          continue ;
        }

        var vin ;
        vin.get(str_is) ;
        lines += vin.num_lines() ;
        lcount += vin.num_lines();
        v = variable(vin.str()) ;
        ::killsp(str_is, lines) ;
        if(str_is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(str_is) ;
          string binfo = process_String(nb.str(),vnames,validate_set) ;
          brackets = "[" + binfo + "]" ;
          lines += nb.num_lines() ;
          lcount += nb.num_lines() ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        lcount += ::kill_sp(str_is, lines) ;
        if(str_is.peek() != '-')
          break ;
        char c=str_is.get() ;
        if(c== '-' && str_is.peek() == '>') {
          c=str_is.get() ;
          lcount += ::kill_sp(str_is, lines) ;
          if(str_is.peek() == '$') {
            str_is.get() ;
            var vin ;
            vin.get(str_is) ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            lcount += ::kill_sp(str_is, lines) ;
            if(str_is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(str_is) ;
              string binfo = process_String(nb.str(),vnames,validate_set) ;
              brk = "[" + binfo +"]";
              lines += nb.num_lines() ;
              lcount += nb.num_lines() ;
            }
            blist.push_back(brk) ;
          } else {
            dangling_arrow = true ;
            break ;
          }
        } else {
          str_is.unget() ;
          break ;
        }
      }
      if(dangling_arrow && first_name) {
        outputFile << name << " ->" ;
        continue ;
      }
      if(dangling_arrow)
        throw parseError("syntax error, near '->' operator") ;

      if(first_name && (vlist.empty())) {
        outputFile << name << ' ' ;
        continue ;
      }

      validate_VariableAccess(v,vlist,first_name,vnames,validate_set) ;

      list<variable>::reverse_iterator ri ;
      for(ri=vlist.rbegin();ri!=vlist.rend();++ri) {
        map<variable,string>::const_iterator vmi = vnames.find(*ri) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << '[' ;
      }
      if(first_name) {
        outputFile << '*' << name ;
      } else {
        map<variable,string>::const_iterator vmi = vnames.find(v) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << v << " is unknown to this rule!" << endl ;
          throw parseError("type error: is this variable in the rule signature?") ;
	}
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      list<string>::const_iterator rbi ;
      for(rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

      for(int i=0;i<lcount;++i)
        outputFile << endl ;
      continue ;

    }

    char c = str_is.get() ;
    if(c == '\n')
      lines++ ;
    outputFile << c ;

  }
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

  istream &get(std::string& filename, int& cnt, int start_line, const std::map<Loci::variable,std::pair<std::string,std::string> > &type_map, istream &is) {
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
  void out_sig_cpu(std::string& filename, std::ostream &outputFile) {
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
  void out_sig_cuda(std::string& filename, std::ostream &outputFile) {

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

  int num_lines() {
    return lines ;
  }
} ;


void parseFile::setup_Rule(std::vector<pair<std::string, std::ostream*> > & dev_outs, bool oneOutput) {
  //read in the rule
  parserule r;
  int start_line  = line_no;
  r.get(filename, cnt, start_line, type_map, is);
  line_no += r.lines;



  for(size_t i = 0 ; i < dev_outs.size(); i++){//cpu version
    if(dev_outs[i].first == "cpu"){
      r.out_sig_cpu(filename, *(dev_outs[i].second));
      ostream& outputFile = *(dev_outs[i].second);
      //output prelude
      if(r.use_prelude){
        string in = r.prelude.str();
        process_Prelude(in, outputFile, r.vnames) ;
      }

      //output compute/calculate
      string in = r.compute.str();
      if(!(r.use_calculate)){
        process_Compute(in,outputFile,r.vnames) ;
      }else{
        if(r.use_compute){
          if(prettyOutput)
            outputFile << "    void calculate(Loci::Entity e) { " << endl ;
          else
            outputFile << "    void calculate(Loci::Entity _e_) { " << endl ;
          process_Calculate(r.compute_line+1, in,outputFile,r.vnames, r.validate_set) ;
        }
        outputFile <<   "    void compute(const Loci::sequence &seq) { " << endl ;
        syncFile(outputFile) ;

        if(r.use_compute) {
          if(r.singletonApply) {
            cerr << "NOTE: parameter only apply rule on '" << r.output << "' now executes single instance." << endl ;
            outputFile <<   "      if(Loci::MPI_rank == 0) calculate(0) ;" << endl ;
            syncFile(outputFile) ;
          } else {
            outputFile <<   "      do_loop(seq,this) ;" << endl ;
            syncFile(outputFile) ;
          }
        }
        outputFile <<   "    }" << endl ;
        syncFile(outputFile) ;
      }
      outputFile <<   "} ;" << endl ;
      syncFile(outputFile) ;

      //register rule
      if(!prettyOutput)
        outputFile << "Loci::register_rule<"<<r.class_name<<"> register_"<<r.class_name
                   << " ;" << endl ;
      else
        outputFile << "register_rule<"<<r.class_name<<"> register_"<<r.class_name
                   << " ;" << endl ;
      syncFile(outputFile) ;

      //end namespace
      if(!prettyOutput) {
        outputFile << "}" << endl ;
        syncFile(outputFile) ;
      }



    }else if(dev_outs[i].first == "cuda"){ //cuda version .
      //  if r is cpu only, it won't have cuda version output
      if(r.cpu_only) continue;


      r.out_sig_cuda(filename, *(dev_outs[i].second));

      ostream& outputFile = *(dev_outs[i].second);

      //no prelude
      //output prelude
      // if(r.use_prelude){
      //   string in = r.prelude.str();
      //   process_Prelude(in, outputFile, r.vnames) ;
      // }

      //output compute/calculate
      string in = r.compute.str();
      //this part won't happen
      // if(!(r.use_calculate)){
      //   process_Compute(in,outputFile,r.vnames) ;
      // }else{
      if(r.use_compute){
        outputFile << "      __host__ __device__ " << endl;
        outputFile << "      void operator()(Loci::Entity e) { " << endl ;
        process_Calculate(r.compute_line+1, in,outputFile,r.vnames, r.validate_set) ;
        outputFile <<   "} ;" << endl ;
      }
      syncFile(outputFile) ;
      outputFile <<   "} ;" << endl ;
      syncFile(outputFile) ;


      //end namespace
      if(!prettyOutput) {
        outputFile << "}" << endl ;
        syncFile(outputFile) ;
      }

    }else{
      throw parseError("unknow device type") ;
    }
  }


}




void parseFile::processFile(string file,vector<pair< string, ostream*> > &outputFile, bool oneOutput) {


  int open_namespace = 0;
  int open_brace = 0;
  bool error = false ;
  filename = file ;
  line_no = 1 ;

  is.open(file.c_str(),ios::in) ;
  if(is.fail()) {
    list<string>::const_iterator li ;
    for(li=include_dirs.begin();li!=include_dirs.end();++li) {
      string s = *li + "/" + file ;
      is.clear() ;
      is.open(s.c_str(),ios::in) ;
      if(!is.fail())
        break ;
    }
    if(is.fail()) {
      string s = "can't open include file '" ;
      s += file ;
      s += "'" ;
      throw parseError(s) ;
    }
  }
  char c ;
  syncFiles(outputFile, oneOutput) ;
  do {
    while(is.peek() == ' ' || is.peek() == '\t') {
      is.get(c) ;
      if(oneOutput){
        *(outputFile[0].second) << c ;
      }else{
        for(size_t i = 0; i < outputFile.size(); i++){
          *(outputFile[i].second) << c ;
        }
      }
    }
    try {
      if(is.peek() == '$') { // Loci specific code!
        is.get(c) ; // get the $
        if(is.peek() == '[') {
          map<variable,string> vnames ;
          int openbrace = 0 ;
          //assume non-cpu version also has special commands
          process_SpecialCommands(is, line_no, outputFile, oneOutput, vnames,openbrace) ;
        } else  if(is_name(is)) {
          std::string key = get_name(is) ;
          if(key == "type") {
            setup_Type(outputFile, oneOutput) ;
          } else if(key == "rule") {
            setup_Rule(outputFile, oneOutput) ;
          } else if(key == "include") {
            killsp() ;
            if(!is_string(is))
              throw parseError("syntax error") ;
            string newfile = get_string(is) ;
            parseFile parser ;
            parser.processFile(newfile,outputFile, oneOutput) ;
            syncFiles(outputFile, oneOutput) ;
            map<variable,pair<string,string> >::const_iterator mi ;
            for(mi=parser.type_map.begin();mi!=parser.type_map.end();++mi)
              type_map[mi->first] = mi->second ;


          } else {
            throw parseError("syntax error: unknown key") ;
          }
        } else {
          throw parseError("unable to process '$' command") ;
        }
      } else {
	bool foundComment = false ;
        while(is.peek() != '\n' && is.peek() != EOF) {
	  if(is_comment(is)) {
            killCommentOuts(is,line_no,outputFile, oneOutput) ;
            foundComment = true ;
	    break ;
	  }
          if(is_namespace(is)) {
            killNamespOuts(is,line_no,outputFile, oneOutput) ;
            open_namespace++;
            open_brace++;
            foundComment = true ;
            break;
          }

          if(is_include(is)) {
            killIncludeOuts(is,line_no,outputFile, oneOutput) ;
            foundComment = true ;
            break;
          }


          is.get(c) ;
          if(c == '{')open_brace++;
          if(c == '}')open_brace--;

          if(oneOutput){
            *(outputFile[0].second) << c ;
          }else{
            for(size_t i = 0; i < outputFile.size(); i++){
              if(outputFile[i].first == "cpu")
                *(outputFile[i].second) << c ;
              else if( c == '}' && open_brace == (open_namespace-1)){
                *(outputFile[i].second) << c ;
                open_namespace--;
              }
            }
          }
        }

	if(!foundComment) {
	  is.get(c) ;
          if(oneOutput){
            *(outputFile[0].second) << endl ;
          }else{
            for(size_t i = 0; i < outputFile.size(); i++){
              *(outputFile[i].second) << endl ;
            }
          }
	  line_no++ ;
	}
      }
      if(is.peek() == EOF){
        is.get(c) ;
      }
    }
    catch(parseError pe) {
      cerr << filename << ':' << line_no << ": " << pe.error_type << endl ;
      char buf[512] ;
      is.getline(buf,512) ;
      line_no++ ;
      //      cerr << "remaining line = " << buf << endl ;
      error = true ;
    }

  } while(!is.eof()) ;

  if(error)
    throw parseError("syntax error") ;
}

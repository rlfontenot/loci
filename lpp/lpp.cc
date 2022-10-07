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

#include "lpp.h"
#include "parsestruct.h"
#include "parserule.h"
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


/*start with " */
bool is_string(istream &s) {
  return s.peek() == '\"' ;
}

/* get everything between " and " */
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


bool is_namespace(istream &s) {
  return is_token(s, "namespace");
}

//multi-ostream version, keep namespace definition inside non-cpu output
//need consider using namespace ***;
istream &killNamespOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
  // if(!is_namespace(s))
  //throw parseError("expected namespace") ;
  string str="namespace" ;
  if(!get_token(s, str)){
    std::cerr << " error in getting namespace " << endl;
  }
  if(s.eof())
    throw parseError("unexpected EOF in killNamespOuts") ;
  int ch = s.get() ;
  while(ch != '{' && ch!= ';' &&!s.eof()) {
    if(ch == '\n')lines++;
    str += ch ;
    ch = s.get() ;
  }
  str += ch ;
  // if(ch!='{')
  //   throw parseError("no beginning { for namespace") ;
  
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
  return is_token(s, "#include");
}

//multi-ostream version, keep "#include ..."  inside non-cpu output
istream &killIncludeOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
  // if(!is_include(s))
  //   throw parseError("expected include") ;
  string str ="#include";
  if(!get_token(s, str)){
    std::cerr << " error in getting #include" << endl;
  }
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


//return line count
static int kill_sp(istream &s, int &lines) {
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
static istream &killCommentOuts(istream &s, int & lines, vector<pair<string, ostream*> >&out, bool oneOutput) {
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




static istream &killspOut(istream &s, int &lines, ostream &out) {

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

//-----------------------------------------------------------------------------------------

static string getNumber(std::istream &is) {
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


//-----------------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------

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


//-----------------------------------------------------------------------------------------


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


//-----------------------------------------------------------------------------------------

class specialCommand{
public:
  string name;
  string body;
};

//-----------------------------------------------------------------------------------------

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


//-----------------------------------------------------------------------------------------

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


//-----------------------------------------------------------------------------------------

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


//-----------------------------------------------------------------------------------------

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




//-----------------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------------

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


//-----------------------------------------------------------------------------------------

void parseFile::setup_Struct(std::vector<pair<std::string, std::ostream*> > & dev_outs, bool oneOutput) {
  //read in the struct
  parsestruct s;
  int start_line  = line_no;
  s.get(filename, cnt, start_line, is);
  line_no += s.lines;
  //output struct
  for(size_t i = 0 ; i < dev_outs.size(); i++){//cpu version
    if(dev_outs[i].first == "cpu"){
      ostream& outputFile = *(dev_outs[i].second);
      outputFile << s.str();
      if(!prettyOutput) syncFile(outputFile) ;
    }
  }
  
}

//-----------------------------------------------------------------------------------------

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


	  }else if(key == "struct"){
	    setup_Struct(outputFile, oneOutput) ;
	  }else {
	    throw parseError("syntax error: unknown key "+key) ;
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
	    //if #include and namespace need special treatment in non-cpu output,
	    //maybe define another bool variable
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

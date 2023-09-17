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
#include <ctype.h>
#include <set>
#include <iostream>
#include <sstream>
//#include <sys/timeb.h>
#include <time.h>
#include <vector>

using std::istringstream ;
using std::ostringstream ;

using std::pair ;
using std::list ;
using std::string ;
using std::set ;
using std::map ;

using std::istream ;
using std::ifstream ;
using std::ofstream ;
using std::ostream ;
using std::ios ;
using std::endl ;
using std::cerr ;
using std::cout ;
using std::vector ;
using namespace Loci ;

bool is_name(istream &s) {
  int ch = s.peek() ;
  return isalpha(ch) || ch == '_' ;
}
    
string get_name(istream &s) {
  if(!is_name(s))
    throw parseError("expected name") ;
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

void parseFile::setup_Type(std::ostream &outputFile) {
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
  outputFile << "// $type " << v << ' ' << tin.str() ;
  int nl = vin.num_lines()+tin.num_lines() ;
  line_no += nl ;
  for(int i=0;i<nl;++i)
    outputFile << endl ;
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

void parseFile::process_SpecialCommand(std::ostream &outputFile,
                                       const map<variable,string> &vnames,
                                       int &openbrace) {
  is.get() ; // get leading [
  string name = get_name(is) ;
  if(is.peek() != ']') {
    cerr << "expecting ']' to close special command '" << name << "'" << endl ;
    throw parseError("syntax error") ;
  }
  is.get() ;

  int nsz = name.size() ;
  for(int i=0;i<nsz;++i)
    if(name[i] >= 'A' || name[i] <= 'Z')
      name[i] = std::tolower(name[i]) ;
  
  if(name == "once") {
    killsp() ;
    if(is.peek() != '{') {
      cerr << "expecting '{' after $[Once] command" << endl ;
      cerr << "found " << char(is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    outputFile << "if(Loci::is_leading_execution()) " ;

  } else if(name == "atomic") {
    killsp() ;
    if(is.peek() != '{') {
      cerr << "expecting '{' after $[Atomic] command" << endl ;
      cerr << "found " << char(is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    is.get() ;
    openbrace++ ;
    outputFile << "{ Loci::atomic_region_helper L__ATOMIC_REGION ; " << endl ;
  } else {
    cerr << "unknown special command '[" << name << "]' !" << endl ;
    throw parseError("syntax error") ;
  }
}

void parseFile::process_Prelude(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    virtual void prelude(const Loci::sequence &seq) { " ;
  is.get() ;
  
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '$') {
      string name ;
      variable v ;
      is.get() ;
      if(is.peek() == '[') {
        process_SpecialCommand(outputFile,vnames,openbrace) ;
        continue ;
      } 
      var vin ;
      vin.get(is) ;
      v = variable(vin.str()) ;
        
      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      outputFile << vmi->second  ;
    }
  
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } ;
}

void parseFile::process_Compute(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    void compute(const Loci::sequence &seq) { " ;
  is.get() ;
  
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
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
    if(is.peek() == '$') {
      variable v ;
      is.get() ;
      if(is.peek() == '[') {
        process_SpecialCommand(outputFile,vnames,openbrace) ;
        continue ;
      }
      bool deref = true ;
      if(is.peek() == '*') {
        is.get() ;
        deref = false ;
      }
      

      var vin ;
      vin.get(is) ;
      v = variable(vin.str()) ;

      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
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
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
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


void parseFile::process_Calculate(std::ostream &outputFile,
                                  const map<variable,string> &vnames,
                                  const set<list<variable> > &validate_set) {
  if(prettyOutput)
    outputFile << "    void calculate(Loci::Entity e) { " << endl ;
  else
    outputFile << "    void calculate(Loci::Entity _e_) { " << endl ;
  is.get() ;
  while(is.peek() == ' ' || is.peek() == '\t')
    is.get() ;
  if(is.peek() == '\n') {
    is.get() ;
    line_no++ ;
  }
  syncFile(outputFile) ;
  killspout(outputFile) ;
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF in process_Calculate") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
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
        killspout(outputFile) ;
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
      killspout(outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }
    
    if(is_name(is) || is.peek() == '$') {
      int lcount = 0 ;
      bool first_name = is_name(is) ;
      if(!first_name) {
        is.get() ;
        if(is.peek() == '[') { // special command
          process_SpecialCommand(outputFile,vnames,openbrace) ;
          continue ;
        }
      }
      string name ;
      variable v ;
      string brackets ;
      if(first_name) 
        name = get_name(is) ;
      else {
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;
          line_no += vin.num_lines() ;
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
        vin.get(is) ;
        line_no += vin.num_lines() ;
        lcount += vin.num_lines();
        v = variable(vin.str()) ;
        killsp() ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames,validate_set) ;
          brackets = "[" + binfo + "]" ;
          line_no += nb.num_lines() ;
          lcount += nb.num_lines() ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        lcount += killsp() ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          lcount += killsp() ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            lcount += killsp() ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames,validate_set) ;
              brk = "[" + binfo +"]";
              line_no += nb.num_lines() ;
              lcount += nb.num_lines() ;
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
    
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
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

// Abstract Syntax Tree
class AST_type : public CPTR_type {
public:
  typedef CPTR<AST_type> ASTP ;
  typedef std::list<ASTP> ASTList ;
  enum elementType {
		    OP_SCOPE=0x000,
		    OP_AT=0x080, // For using @ to separate namespaces
		    // Traditional C operators
		    OP_ARROW=0x100, 
		    OP_TIMES = 0x300, OP_DIVIDE, OP_MODULUS,
		    OP_PLUS  = 0x400, OP_MINUS, 
		    OP_SHIFT_RIGHT = 0x500, OP_SHIFT_LEFT,
		    OP_LT = 0x600, OP_GT, OP_GE, OP_LE,
		    OP_EQUAL = 0x700, OP_NOT_EQUAL, 
		    OP_AND=0x800, OP_EXOR=0x900, OP_OR=0xa00,
		    OP_LOGICAL_AND=0xb00, OP_LOGICAL_OR=0xc00,
		    OP_TERTIARY,
		    OP_ASSIGN=0xd00,
		    OP_TIMES_ASSIGN,
		    OP_DIVIDE_ASSIGN,
		    OP_MODULUS_ASSIGN,
		    OP_PLUS_ASSIGN,
		    OP_MINUS_ASSIGN,
		    OP_SHIFT_LEFT_ASSIGN,
		    OP_SHIFT_RIGHT_ASSIGN,
		    OP_AND_ASSIGN,
		    OP_OR_ASSIGN,
		    OP_EXOR_ASSIGN,
		    OP_COMMA=0xe00, OP_DOT,
		    OP_COLON=0xf00,
		    OP_SEMICOLON,
		    // terminal for empty statement
		    OP_NIL=0x1000,
		    // terminals for variable name, function, array or name{args}
		    OP_INCREMENT,OP_DECREMENT,
		    OP_POSTINCREMENT, OP_POSTDECREMENT,
		    OP_COMMENT,
		    OP_BRACEBLOCK,
		    OP_NAME, OP_FUNC, OP_ARRAY, OP_NAME_BRACE, OP_FUNC_BRACE,
		    // terminal for string, integer, or unspecified error condition
		    OP_STRING, OP_NUMBER, OP_ERROR,
		    // Unary operations
		    OP_UNARY_PLUS, OP_UNARY_MINUS, OP_NOT, OP_TILDE,
		    OP_AMPERSAND, OP_DOLLAR, OP_STAR,
		    OP_GROUP,OP_GROUP_ERROR,
		    OP_OPENPAREN,OP_CLOSEPAREN,OP_OPENBRACKET,OP_CLOSEBRACKET,
		    OP_OPENBRACE,OP_CLOSEBRACE,
		    OP_LOCI_DIRECTIVE,OP_LOCI_VARIABLE,OP_LOCI_CONTAINER,
		    OP_TERM, OP_SPECIAL,
		    TK_BRACEBLOCK=0x2000,
		    TK_SCOPE,
		    TK_AT, // For using @ to separate namespaces
		    // Traditional C operators
		    TK_ARROW, 
		    TK_TIMES, TK_DIVIDE, TK_MODULUS,
		    TK_PLUS, TK_MINUS, 
		    TK_SHIFT_RIGHT, TK_SHIFT_LEFT,
		    TK_LT, TK_GT, TK_GE, TK_LE,
		    TK_EQUAL, TK_NOT_EQUAL, 
		    TK_AND, TK_EXOR, TK_OR,
		    TK_LOGICAL_AND, TK_LOGICAL_OR, 
		    TK_ASSIGN,
		    TK_TIMES_ASSIGN,
		    TK_DIVIDE_ASSIGN,
		    TK_MODULUS_ASSIGN,
		    TK_PLUS_ASSIGN,
		    TK_MINUS_ASSIGN,
		    TK_SHIFT_LEFT_ASSIGN,
		    TK_SHIFT_RIGHT_ASSIGN,
		    TK_AND_ASSIGN,
		    TK_OR_ASSIGN,
		    TK_EXOR_ASSIGN,
		    TK_COMMA, TK_DOT,
		    TK_COLON,
		    TK_SEMICOLON,
		    // terminal for empty statement
		    TK_NIL,
		    // terminals for variable name, function, array or name{args}
		    TK_INCREMENT,TK_DECREMENT,TK_COMMENT,
		    TK_NAME, 
		    // terminal for string, integer, or unspecified error condition
		    TK_STRING, TK_NUMBER, TK_ERROR,
		    // Unary operations
		    TK_UNARY_PLUS, TK_UNARY_MINUS, TK_NOT, TK_TILDE,
		    TK_QUESTION, TK_AMPERSAND, TK_STAR,
		    TK_OPENPAREN,TK_CLOSEPAREN,TK_OPENBRACKET,TK_CLOSEBRACKET,
		    TK_OPENBRACE,TK_CLOSEBRACE,
		    TK_LOCI_DIRECTIVE,TK_LOCI_VARIABLE,TK_LOCI_CONTAINER,
		    // Now the keywords
		    TK_ALIGNAS, TK_ALIGNOF, TK_ASM, 
		    TK_BOOL,TK_FALSE,TK_TRUE,TK_CHAR,TK_INT,TK_LONG,
		    TK_SHORT,TK_SIGNED,TK_UNSIGNED,TK_DOUBLE,TK_FLOAT,TK_ENUM,
		    TK_MUTABLE,TK_CONST,TK_STATIC,TK_VOLATILE,TK_AUTO,
		    TK_REGISTER,TK_EXPORT,TK_EXTERN,TK_INLINE,TK_NAMESPACE,
		    TK_EXPLICIT,TK_DYNAMIC_CAST,TK_STATIC_CAST,
		    TK_REINTERPRET_CAST,
		    TK_OPERATOR,TK_PROTECTED,TK_NOEXCEPT,TK_NULLPTR,
		    TK_RETURN,TK_SIZEOF,TK_THIS,TK_TYPEID,
		    TK_SWITCH,TK_CASE,TK_BREAK,TK_DEFAULT,
		    TK_FOR,TK_DO,TK_WHILE,TK_CONTINUE,
		    TK_CLASS,TK_STRUCT,TK_PUBLIC,TK_PRIVATE,TK_FRIEND,
		    TK_UNION,TK_TYPENAME,TK_TEMPLATE,TK_TYPEDEF,
		    TK_VIRTUAL,TK_VOID,TK_TRY,TK_CATCH,TK_THROW,
		    TK_IF,TK_ELSE,TK_GOTO,TK_NEW,TK_DELETE,
		    TK_SENTINEL 
		    
  } ;
  elementType nodeType ;
  virtual void DiagPrint(ostream &s, int &line) const = 0 ;
  
} ;

class AST_Token : public AST_type {
public:
  string text ;
  int lineno ;
  void DiagPrint(ostream &s, int &line) const {
    if(line != lineno) {
      cout << endl ;
      
      if(line < 0 || line+1 != lineno) {
	cout << "#line " << lineno << endl ;
      }
      line = lineno ;
    }
    s <<text << ' ' ;
  }

} ;

class AST_SimpleStatement: public AST_type {
public:
  AST_SimpleStatement(ASTP e, ASTP t) : exp(e),Terminal(t) {}
  ASTP exp ;
  ASTP Terminal ;
  void DiagPrint(ostream  &s, int &line) const {
    exp->DiagPrint(s,line) ;
    Terminal->DiagPrint(s,line) ;
  }
} ;

class AST_Block : public AST_type {
public:
  ASTList elements ;
  void DiagPrint(ostream &s, int &lineno) const {
    for(ASTList::const_iterator ii=elements.begin();ii!=elements.end();++ii)
      (*ii)->DiagPrint(s,lineno) ;
  }
} ;

class AST_typeDecl : public AST_type {
public:
  ASTList type_decl ;
  void DiagPrint(ostream &s, int &lineno) const {
    for(ASTList::const_iterator ii=type_decl.begin();ii!=type_decl.end();++ii)
      (*ii)->DiagPrint(s,lineno) ;
  }
} ;
  
class AST_declaration : public AST_type {
public:
  ASTP type_decl ;
  ASTP decls ;
  void DiagPrint(ostream &s, int &lineno) const {
    type_decl->DiagPrint(s,lineno) ;
    decls->DiagPrint(s,lineno) ;
  }
} ;

string OPtoString(AST_type::elementType val) {
  switch(val) {
  case AST_type::OP_SCOPE:
    return string("::") ;
  case AST_type::OP_AT:
    return string("@") ;
  case AST_type::OP_ARROW:
    return string("->") ;
  case AST_type::OP_TIMES:
    return string("*") ;
  case AST_type::OP_DIVIDE:
    return string("/") ;
  case AST_type::OP_MODULUS:
    return string("%") ;
  case AST_type::OP_PLUS:
    return string("+") ;
  case AST_type::OP_MINUS:
    return string("-") ;
  case AST_type::OP_SHIFT_RIGHT:
    return string(">>") ;
  case AST_type::OP_SHIFT_LEFT:
    return string("<<") ;
  case AST_type::OP_LT:
    return string("<") ;
  case AST_type::OP_GT:
    return string(">") ;
  case AST_type::OP_GE:
    return string(">=") ;
  case AST_type::OP_LE:
    return string("<=") ;
  case AST_type::OP_EQUAL:
    return string("==") ;
  case AST_type::OP_NOT_EQUAL:
    return string("!=") ;
  case AST_type::OP_AND:
    return string("&") ;
  case AST_type::OP_EXOR:
    return string("^") ;
  case AST_type::OP_OR:
    return string("|") ;
  case AST_type::OP_LOGICAL_AND:
    return string("&&") ;
  case AST_type::OP_LOGICAL_OR:
    return string("||") ;
  case AST_type::OP_ASSIGN:
    return string("=") ;
  case AST_type::OP_TIMES_ASSIGN:
    return string("*=") ;
  case AST_type::OP_DIVIDE_ASSIGN:
    return string("/=") ;
  case AST_type::OP_MODULUS_ASSIGN:
    return string("%=") ;
  case AST_type::OP_PLUS_ASSIGN:
    return string("+=") ;
  case AST_type::OP_MINUS_ASSIGN:
    return string("-=") ;
  case AST_type::OP_SHIFT_LEFT_ASSIGN:
    return string("<<=") ;
  case AST_type::OP_SHIFT_RIGHT_ASSIGN:
    return string(">>=") ;
  case AST_type::OP_AND_ASSIGN:
    return string("&=") ;
  case AST_type::OP_OR_ASSIGN:
    return string("|=") ;
  case AST_type::OP_EXOR_ASSIGN:
    return string("^=") ;
  case AST_type::OP_COMMA:
    return string(",") ;
  case AST_type::OP_DOT:
    return string(".") ;
  case AST_type::OP_COLON:
    return string(":") ;
  case AST_type::OP_SEMICOLON:
    return string(";") ;
  case AST_type::OP_INCREMENT:
    return string(" ++") ;
  case AST_type::OP_DECREMENT:
    return string(" --") ;
  case AST_type::OP_POSTINCREMENT:
    return string("++ ") ;
  case AST_type::OP_POSTDECREMENT:
    return string("-- ") ;
  case AST_type::OP_UNARY_PLUS:
    return string("+") ;
  case AST_type::OP_UNARY_MINUS:
    return string("-") ;
  case AST_type::OP_NOT:
    return string("!") ;
  case AST_type::OP_TILDE:
    return string("~") ;
  case AST_type::OP_AMPERSAND:
    return string("&") ;
  case AST_type::OP_TERTIARY:
    return string("?") ;
  case AST_type::OP_DOLLAR:
    return string("$") ;
  case AST_type::OP_STAR:
    return string("*") ;
  default:
    return string("/*error*/") ;
  }
  return string("/*error*/") ;
}

class AST_exprOper : public AST_type {
public:
  ASTList terms ;
  void DiagPrint(ostream &s, int &lineno) const {
    switch (nodeType) {
    case OP_GROUP:
      s << '(' ;
      for(ASTList::const_iterator ii=terms.begin();ii!=terms.end();++ii)
	(*ii)->DiagPrint(s,lineno) ;
      s << ')' ;
      break ;
    case OP_FUNC:
      {
	ASTList::const_iterator ii=terms.begin() ;
	(*ii)->DiagPrint(s,lineno) ;
	++ii ;
	s << '(' ;
	(*ii)->DiagPrint(s,lineno) ;
	s << ')' ;
	++ii ;
	if(ii!=terms.end()) {
	  cerr << "syntax error" ;
	  (*ii)->DiagPrint(cerr,lineno) ;
	}
      }
      break ;
    case OP_ARRAY:
      {
	ASTList::const_iterator ii=terms.begin() ;
	(*ii)->DiagPrint(s,lineno) ;
	++ii ;
	s << '[' ;
	(*ii)->DiagPrint(s,lineno) ;
	s << ']' ;
	++ii ;
	if(ii!=terms.end()) {
	  cerr << "syntax error" ;
	  (*ii)->DiagPrint(cerr,lineno) ;
	}
      }
      break ;
    case OP_TERTIARY:
      {
	ASTList::const_iterator ii=terms.begin() ;
	(*ii)->DiagPrint(s,lineno) ;
	++ii ;
	s << '?' ;
	(*ii)->DiagPrint(s,lineno) ;
	s << ':' ;
	if(ii==terms.end()) {
	  cerr << "syntax error on tertiary operator" << endl ;
	} else
	  ++ii ;
	(*ii)->DiagPrint(s,lineno) ;

      }
      break ;
      
    case OP_UNARY_PLUS:
    case OP_UNARY_MINUS:
    case OP_NOT:
    case OP_AMPERSAND:
    case OP_STAR:
    case OP_INCREMENT:
    case OP_DECREMENT:
      {
	string op = OPtoString(nodeType) ;
	s << op ;
	for(ASTList::const_iterator ii=terms.begin();ii!=terms.end();++ii)
	  (*ii)->DiagPrint(s,lineno) ;
      }
      break ;
    case OP_POSTINCREMENT:
    case OP_POSTDECREMENT:
      {
	string op = OPtoString(nodeType) ;
	for(ASTList::const_iterator ii=terms.begin();ii!=terms.end();++ii)
	  (*ii)->DiagPrint(s,lineno) ;
	s << op ;
      }
      break ;
    default:
      {
	//	s << "[" ;
	string op = OPtoString(nodeType) ;
	for(ASTList::const_iterator ii=terms.begin();ii!=terms.end();) {
	  (*ii)->DiagPrint(s,lineno) ;
	  ++ii ;
	  if(ii != terms.end())
	    s << op ;
	}
	//	s << "]" ;
      }
      break ;
    }
  }
} ;

class AST_term : public AST_type {
public:
  enum TermTypes {
		 TERM_NUMBER,
		 TERM_STRING,
		 TERM_BOOLEAN,
		 TERM_VARIABLE,
		 TERM_LOCIVARIABLE,
		 TERM_LOCICONTAINER,
		 TERM_INVALID
		 
  } ;
  TermTypes TermType ;
  ASTP term ;
  void DiagPrint(ostream&s, int &lineno) const {
    term->DiagPrint(s,lineno) ;
  }
} ;

class AST_ifStatement : public AST_type {
public:
  ASTP iftok ;
  ASTP conditional ;
  ASTP ifblock ;
  ASTP elseblock ;
  AST_ifStatement(ASTP tok, ASTP C, ASTP IF, ASTP ELSE):
    iftok(tok),conditional(C),ifblock(IF),elseblock(ELSE) {} 
  void DiagPrint(ostream&s, int &lineno) const {
    iftok->DiagPrint(s,lineno) ;
    s << "(" ;
    conditional->DiagPrint(s,lineno) ;
    s << ")" ;
    ifblock->DiagPrint(s,lineno) ;
    if(elseblock != 0) {
      s << " else " ;
      elseblock->DiagPrint(s,lineno) ;
    }
  }
} ;
class AST_loopStatement : public AST_type {
public:
  ASTP loop ;
  ASTP initializer ;
  ASTP conditional ;
  ASTP advance ;
  ASTP body ;
  AST_loopStatement(ASTP L, ASTP I, ASTP C, ASTP A, ASTP B):
    loop(L),initializer(I),conditional(C),advance(A),body(B) {}
  AST_loopStatement(ASTP L, ASTP C, ASTP B):
    loop(L),initializer(0),conditional(C),advance(0),body(B) {}
  void DiagPrint(ostream&s, int &lineno) const {
    loop->DiagPrint(s,lineno) ;
    if(loop->nodeType == AST_type::TK_FOR) {
      s << "(" ;
      initializer->DiagPrint(s,lineno) ;
      s << ";" ;
      conditional->DiagPrint(s,lineno) ;
      s << ";" ;
      advance->DiagPrint(s,lineno) ;
      s << ")" ;
      body->DiagPrint(s,lineno) ;
    } else if(loop->nodeType == AST_type::TK_WHILE) {
      s << "(" ;
      conditional->DiagPrint(s,lineno) ;
      s << ")" ;
      body->DiagPrint(s,lineno) ;
    } else {
      body->DiagPrint(s,lineno) ;
      s << "while(" ;
      conditional->DiagPrint(s,lineno) ;
      s << ") ;" ;
    }
  }
} ;
  
class AST_switchStatement : public AST_type {
public:
  ASTP statement ;
  ASTP conditional ;
  ASTList body ;
  AST_switchStatement() {}
  void DiagPrint(ostream&s, int &lineno) const {
    statement->DiagPrint(s,lineno) ;
    s << "(" ;
    conditional->DiagPrint(s,lineno) ;
    s << ") {" ;
    for(ASTList::const_iterator ii=body.begin();ii!=body.end();++ii)
      (*ii)->DiagPrint(s,lineno) ;
    s << "}" ;
    
  }
} ;

vector<CPTR<AST_Token> >  tokenStack ;

inline void pushToken(CPTR<AST_Token> &pt) {
  tokenStack.push_back(pt) ;
}

extern CPTR<AST_Token> getToken(std::istream &is, int &linecount) ;
bool isTerm(AST_type::elementType e) {
  return ((e == AST_type::TK_STRING) ||
	  (e == AST_type::TK_NAME) ||
	  (e == AST_type::TK_NUMBER) ||
	  (e == AST_type::TK_TRUE) ||
	  (e == AST_type::TK_FALSE) ||
	  (e == AST_type::TK_LOCI_VARIABLE) ||
	  (e == AST_type::TK_LOCI_CONTAINER) ) ;
}

AST_type::ASTP parseTerm(std::istream &is, int &linecount) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
  CPTR<AST_term> info = new AST_term ;
  info->nodeType = AST_type::OP_TERM ;
  info->term = AST_type::ASTP(token) ;
  info->TermType = AST_term::TERM_INVALID ;
  switch(token->nodeType) {
  case AST_type::TK_NAME:
    info->TermType = AST_term::TERM_VARIABLE ;
    return AST_type::ASTP(info) ;
  case AST_type::TK_NUMBER:
    info->TermType = AST_term::TERM_NUMBER ;
    return AST_type::ASTP(info) ;
  case AST_type::TK_STRING:
    info->TermType = AST_term::TERM_STRING ;
    return AST_type::ASTP(info) ;
  case AST_type::TK_TRUE:
  case AST_type::TK_FALSE:
    info->TermType = AST_term::TERM_BOOLEAN ;
    return AST_type::ASTP(info) ;
  case AST_type::TK_LOCI_VARIABLE:
    info->TermType = AST_term::TERM_LOCIVARIABLE ;
    return AST_type::ASTP(info) ;
  case AST_type::TK_LOCI_CONTAINER:
    info->TermType = AST_term::TERM_LOCICONTAINER ;
    return AST_type::ASTP(info) ;
  default:
    return AST_type::ASTP(0) ;
  }
  return AST_type::ASTP(0) ;
}
  
CPTR<AST_Token> getNumberToken(std::istream &is, int &linecount) {
  // Process elements that start with 0-9, '.'
  CPTR<AST_Token> AST_data = new AST_Token() ;
  AST_data->lineno = linecount ;
  string numberdata ;
  bool hasPoint = false ;
  if(is.peek() == '.')  { // This might be a floating point number
    numberdata += is.get() ;
    if(is.peek() <'0' || is.peek() > '9') {
      // Not a floating point number
      AST_data->nodeType = AST_type::TK_DOT ;
      AST_data->text = numberdata ;
      return AST_data ;
    }
    hasPoint = true ;
  }
  if((is.peek() >= '0' && is.peek() <='9')) {
    AST_data->nodeType = AST_type::TK_NUMBER ;
    string numberdata ;
    if(is.peek() == '0' && !hasPoint) {
      // This path is either binary, octal, or hexidecimal
      numberdata += is.get() ;
      if(is.peek() == 'x' || is.peek() == 'X') { //hex number
	numberdata += is.get() ;
	while((is.peek() >='0' && is.peek() <= '9') ||
	      (is.peek() >='a' && is.peek() <= 'f') ||
	      (is.peek() >='A' && is.peek() <= 'F')) {
	  numberdata += is.get() ;
	}
	if(numberdata.size() == 2)
	  AST_data->nodeType = AST_type::TK_ERROR ;
	while(is.peek() == 'l' || is.peek() =='L' ||
	      is.peek() == 'u' || is.peek() =='U')
	  numberdata += is.get() ;
	AST_data->text = numberdata ;
	return AST_data ;
      }
      if(is.peek() == 'b' || is.peek() == 'B') { // binary number
	numberdata += is.get() ;
	while(is.peek() >= '0' && is.peek() <= '1') {
	  numberdata += is.get() ;
	}
	if(numberdata.size() == 2 || (is.peek() >= '2' && is.peek() <='9'))
	  AST_data->nodeType = AST_type::TK_ERROR ;
	while(is.peek() == 'l' || is.peek() =='L' ||
	      is.peek() == 'u' || is.peek() =='U')
	  numberdata += is.get() ;
	AST_data->text = numberdata ;
	return AST_data ;
      }
      // octal number
      while(is.peek() >= '0' && is.peek() <= '7') {
	numberdata += is.get() ;
      }
      if(is.peek() >= '8' && is.peek() <='9')
	AST_data->nodeType = AST_type::TK_ERROR ;
      while(is.peek() == 'l' || is.peek() =='L' ||
	    is.peek() == 'u' || is.peek() =='U')
	numberdata += is.get() ;
      AST_data->text = numberdata ;
      return AST_data ;
    }
    while(is.peek() >= '0' && is.peek() <= '9') {
      numberdata += is.get() ;
    }
    if(!hasPoint && is.peek() == '.') {
      numberdata += is.get() ;
      hasPoint = true ;
      while(is.peek() >= '0' && is.peek() <= '9') {
	numberdata += is.get() ;
      }
    }
    bool isFloat = hasPoint ;
    if(is.peek() == 'e' || is.peek() == 'E') {
      isFloat = true ;
      numberdata += is.get() ;
      // get exponent
      if(is.peek()=='+' || is.peek() == '-')
	numberdata += is.get() ;
      while(is.peek() >= '0' && is.peek() <= '9') {
	numberdata += is.get() ;
      }
    }
    if(isFloat) {
      if(is.peek()=='f' || is.peek()=='F' || is.peek()=='l' || is.peek()=='L')
	numberdata += is.get() ;
    } else {
      while(is.peek()=='l' || is.peek() == 'L' || is.peek() == 'u' || is.peek()=='U')
	numberdata += is.get() ;
    }
    AST_data->text=numberdata ;
    return AST_data ;
  }
  AST_data->text=numberdata ;
  AST_data->nodeType = AST_type::TK_ERROR ;
  return AST_data ;
}

struct keywords {
  const char *keyword ;
  AST_type::elementType nodeType;
} ;

keywords keywordDictionary[] = {
  {"alignas", AST_type::TK_ALIGNAS},
  {"alignof", AST_type::TK_ALIGNOF},
  {"and", AST_type::TK_LOGICAL_AND},
  {"and_eq", AST_type::TK_AND_ASSIGN},
  {"asm", AST_type::TK_ASM},
  {"auto", AST_type::TK_AUTO},
  {"bitand", AST_type::TK_AND},
  {"bitor", AST_type::TK_OR},
  {"bool", AST_type::TK_BOOL},
  {"break", AST_type::TK_BREAK},
  {"case", AST_type::TK_CASE},
  {"catch", AST_type::TK_CATCH},
  {"char", AST_type::TK_CHAR},
  {"class", AST_type::TK_CLASS},
  {"const", AST_type::TK_CONST},
  {"continue", AST_type::TK_CONTINUE},
  {"default", AST_type::TK_DEFAULT},
  {"delete", AST_type::TK_DELETE},
  {"do", AST_type::TK_DO},
  {"double", AST_type::TK_DOUBLE},
  {"dynamic_cast", AST_type::TK_DYNAMIC_CAST},
  {"else", AST_type::TK_ELSE},
  {"enum", AST_type::TK_ENUM},
  {"explicit", AST_type::TK_EXPLICIT},
  {"export", AST_type::TK_EXPORT},
  {"extern", AST_type::TK_EXTERN},
  {"false", AST_type::TK_FALSE},
  {"float", AST_type::TK_FLOAT},
  {"for", AST_type::TK_FOR},
  {"friend", AST_type::TK_FRIEND},
  {"goto", AST_type::TK_GOTO},
  {"if", AST_type::TK_IF},
  {"inline", AST_type::TK_INLINE},
  {"int", AST_type::TK_INT},
  {"long", AST_type::TK_LONG},
  {"mutable", AST_type::TK_MUTABLE},
  {"namespace", AST_type::TK_NAMESPACE},
  {"new", AST_type::TK_NEW},
  {"noexcept", AST_type::TK_NOEXCEPT},
  {"not", AST_type::TK_NOT},
  {"not_eq", AST_type::TK_EQUAL},
  {"nullptr", AST_type::TK_NULLPTR},
  {"operator", AST_type::TK_OPERATOR},
  {"or", AST_type::TK_LOGICAL_OR},
  {"or_eq", AST_type::TK_OR_ASSIGN},
  {"private", AST_type::TK_PRIVATE},
  {"protected",AST_type::TK_PROTECTED},
  {"public", AST_type::TK_PUBLIC},
  {"register", AST_type::TK_REGISTER},
  {"reinterpret_cast", AST_type::TK_REINTERPRET_CAST},
  {"return", AST_type::TK_RETURN},
  {"short", AST_type::TK_SHORT},
  {"signed", AST_type::TK_SIGNED},
  {"sizeof", AST_type::TK_SIZEOF},
  {"static", AST_type::TK_STATIC},
  {"static_cast", AST_type::TK_STATIC_CAST},
  {"struct", AST_type::TK_STRUCT},
  {"switch", AST_type::TK_SWITCH},
  {"template", AST_type::TK_TEMPLATE},
  {"this", AST_type::TK_THIS},
  {"throw", AST_type::TK_THROW},
  {"true", AST_type::TK_TRUE},
  {"try", AST_type::TK_TRY},
  {"typedef", AST_type::TK_TYPEDEF},
  {"typeid", AST_type::TK_TYPEID},
  {"typename", AST_type::TK_TYPENAME},
  {"union", AST_type::TK_UNION},
  {"unsigned", AST_type::TK_UNSIGNED},
  {"virtual", AST_type::TK_VIRTUAL},
  {"void", AST_type::TK_VOID},
  {"volatile", AST_type::TK_VOLATILE},
  {"while", AST_type::TK_WHILE},
  {"xor",  AST_type::TK_EXOR},
  {"xor_eq", AST_type::TK_EXOR_ASSIGN}
} ;

inline AST_type::elementType findToken(const string &name) {
  const int ntok = sizeof(keywordDictionary)/sizeof(keywords) ;
  int start = ntok/2 ;
  if(name == keywordDictionary[start].keyword) 
    return keywordDictionary[start].nodeType ;
  int min = 0 ;
  int max = ntok-1 ;
  if(name > keywordDictionary[start].keyword) 
    min = start ;
  else
    max = start ;
  while (max-min > 1) {
    start = (min+max)/2 ;
    if(name == keywordDictionary[start].keyword)
      return keywordDictionary[start].nodeType ;
    if(name > keywordDictionary[start].keyword) 
      min = start ;
    else
      max = start ;
  }
  if(name == keywordDictionary[min].keyword)
    return keywordDictionary[min].nodeType ;
  if(name == keywordDictionary[max].keyword)
    return keywordDictionary[max].nodeType ;
  return AST_type::TK_NAME ;
}
    
     
    
CPTR<AST_Token> getToken(std::istream &is, int &linecount) {
  if(tokenStack.size() != 0) {
    CPTR<AST_Token> tok = tokenStack.back() ;
    tokenStack.pop_back() ;
    return tok ;
  }
  
  killsp(is,linecount) ;
  if(is.fail() || is.eof()) {
    CPTR<AST_Token> AST_data = new AST_Token() ;
    AST_data->lineno = linecount ;
    AST_data->nodeType = AST_type::TK_ERROR ;
    return AST_data ;
  }
  if(is_name(is)) {
    CPTR<AST_Token> AST_data = new AST_Token() ;
    AST_data->text = get_name(is) ;
    AST_data->lineno = linecount ;
    AST_data->nodeType = findToken(AST_data->text) ;
    return AST_data ;
  }
  if(is.peek() == '.' || (is.peek() >= '0' && is.peek() <='9')) {
    return getNumberToken(is,linecount) ;
  }
  CPTR<AST_Token> AST_data = new AST_Token() ;
  AST_data->lineno = linecount ;
  AST_data->text += is.get() ;
  AST_data->nodeType = AST_type::TK_ERROR ;
  switch(AST_data->text[0]) {
  case '+':
    if(is.peek()=='+') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_INCREMENT ;
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_PLUS_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_PLUS ;
    return AST_data ;
  case '-':
    if(is.peek()=='-') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_DECREMENT ;
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_MINUS_ASSIGN ;
      return AST_data ;
    }
    if(is.peek()=='>') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_ARROW ;
      return AST_data ;
    }      
    AST_data->nodeType = AST_type::TK_MINUS ;
    return AST_data ;
  case '*':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_TIMES_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_TIMES ;
    return AST_data ;
  case '/':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_DIVIDE_ASSIGN ;
      return AST_data ;
    }
    if(is.peek()=='/') {
      AST_data->text += is.get() ;
      while(is.peek() != '\n' && !is.eof() && !is.fail()) {
	AST_data->text += is.get() ;
      }
      AST_data->nodeType = AST_type::TK_COMMENT ;
      return AST_data ;
    }
    if(is.peek()=='*') {
      AST_data->text += is.get() ;
      for(;;) {
	if(is.eof() || is.fail())
	  break ;
	if(is.peek() != '*') {
	  AST_data->text += is.get() ;
	  if(is.peek() == '/') {
	    AST_data->text = is.get() ;
	    AST_data->nodeType=AST_type::TK_COMMENT ;
	    return AST_data ;
	  }
	} else {
	  if(is.peek() == '\n')
	    linecount++ ;
	  AST_data->text += is.get() ;
	}
      }
      AST_data->nodeType=AST_type::TK_ERROR ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_DIVIDE ;
    return AST_data ;
  case '%':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_MODULUS_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_MODULUS ;
    return AST_data ;
  case ',':
    AST_data->nodeType = AST_type::TK_COMMA ;
    return AST_data ;
  case '@':
    AST_data->nodeType = AST_type::TK_AT ;
    return AST_data ;
  case '&':
    if(is.peek() == '&' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LOGICAL_AND ;
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_AND_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_AND ;
    return AST_data ;
  case '|':
    if(is.peek() != '|') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LOGICAL_OR ;
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_OR_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_OR ;
    return AST_data ;
  case '^':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_EXOR_ASSIGN ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_EXOR ;
    return AST_data ;
  case '=':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_EQUAL ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_ASSIGN ;
    return AST_data ;
  case '!':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_NOT_EQUAL ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_NOT ;
    return AST_data ;
  case ':':
    if(is.peek() == ':') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_SCOPE ;
      return AST_data ;
    }

    AST_data->nodeType = AST_type::TK_COLON ;
    return AST_data ;
  case '<':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LE ;
      return AST_data ;
    }
    if(is.peek() == '<') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = AST_type::TK_SHIFT_LEFT_ASSIGN ;
	return AST_data ;
      }
      AST_data->nodeType = AST_type::TK_SHIFT_LEFT ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_LT ;
    return AST_data ;
  case '>':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_GE ;
      return AST_data ;
    }
    if(is.peek() == '>') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = AST_type::TK_SHIFT_RIGHT_ASSIGN ;
	return AST_data ;
      }
      AST_data->nodeType = AST_type::TK_SHIFT_RIGHT ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_GT ;
    return AST_data ;
  case '~':
    AST_data->nodeType = AST_type::TK_TILDE ;
    return AST_data ;
  case '?':
    AST_data->nodeType = AST_type::TK_QUESTION ;
    return AST_data ;
  case '"':
    while(is.peek() != '"' && !is.eof() && !is.fail()) {
      AST_data->text += is.get() ;
    }
    if(is.peek() == '"') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_STRING ;
      return AST_data ;
    }
    break ;
  case '\'':
    while(is.peek() != '\'' && !is.eof() && !is.fail()) {
      AST_data->text += is.get() ;
    }
    if(is.peek() == '\'') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_STRING ;
      return AST_data ;
    }
    break ;
  case ';':
    AST_data->nodeType = AST_type::TK_SEMICOLON ;
    return AST_data ;
  case '(':
    AST_data->nodeType = AST_type::TK_OPENPAREN ;
    return AST_data ;
  case ')':
    AST_data->nodeType = AST_type::TK_CLOSEPAREN ;
    return AST_data ;
  case '[':
    AST_data->nodeType = AST_type::TK_OPENBRACKET ;
    return AST_data ;
  case ']':
    AST_data->nodeType = AST_type::TK_CLOSEBRACKET ;
    return AST_data ;
  case '{':
    AST_data->nodeType = AST_type::TK_OPENBRACE ;
    return AST_data ;
  case '}':
    AST_data->nodeType = AST_type::TK_CLOSEBRACE ;
    return AST_data ;
  case '$':
    // Now we have a Loci variable or a Loci command
    if(is.peek() == '[') { // This is a Loci command
      AST_data->text += is.get() ;
      while(!is.fail() && !is.eof() && is.peek() != ']') {
	if(is.peek() == '[') {
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '[')
	      count++ ;
	    if(is.peek() == ']')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else 
	  AST_data->text += is.get() ;
      }
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LOCI_DIRECTIVE ;
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_LOCI_VARIABLE ;
    if(is.peek() == '*') {
      AST_data->nodeType = AST_type::TK_LOCI_CONTAINER ;
      AST_data->text += is.peek() ;
    }
    if(!is.fail() && !is.eof() &&
       ((is.peek() >='a' && is.peek() <='z') ||
	(is.peek() >='A' && is.peek() <='Z') ||
	is.peek() == '_' || is.peek() == '@')) {
      // extract variable name
      while(!is.fail() && !is.eof() &&
	    ((is.peek() >='a' && is.peek() <='z') ||
	     (is.peek() >='A' && is.peek() <='Z') ||
	     is.peek() == '_' || is.peek() == '@' ||
	     is.peek() == '(' || is.peek() == '{')) {
	if(is.peek() == '(') {
	  AST_data->text +=is.get() ;
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '(')
	      count++ ;
	    if(is.peek() == ')')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else if(is.peek() == '{') {
	  AST_data->text +=is.get() ;
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '{')
	      count++ ;
	    if(is.peek() == '}')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else
	  AST_data->text += is.get() ;
      }

    }
    return AST_data ;
  } 
     
      
  AST_data->nodeType = AST_type::TK_ERROR ;
  return AST_data ;

}

AST_type::ASTP parseBlock(std::istream &is, int &linecount) ;

AST_type::ASTP parseOperator(std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  switch(openToken->nodeType) {
  case AST_type::TK_SCOPE:
    openToken->nodeType = AST_type::OP_SCOPE ;
    return AST_type::ASTP(openToken) ;
    //-----------------------------------
    // For using @ to separate namespaces
  case AST_type::TK_AT:
    openToken->nodeType = AST_type::OP_AT ;
    return AST_type::ASTP(openToken) ;
    //-----------------------------------
    // Traditional C operators
  case AST_type::TK_ARROW:
    openToken->nodeType = AST_type::OP_ARROW ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_TIMES:
    openToken->nodeType = AST_type::OP_TIMES ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DIVIDE:
    openToken->nodeType = AST_type::OP_DIVIDE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MODULUS:
    openToken->nodeType = AST_type::OP_MODULUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_PLUS:
    openToken->nodeType = AST_type::OP_PLUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MINUS:
    openToken->nodeType = AST_type::OP_MINUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_RIGHT:
    openToken->nodeType = AST_type::OP_SHIFT_RIGHT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_LEFT:
    openToken->nodeType = AST_type::OP_SHIFT_LEFT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LT:
    openToken->nodeType = AST_type::OP_LT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_GT:
    openToken->nodeType = AST_type::OP_GT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_GE:
    openToken->nodeType = AST_type::OP_GE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LE:
    openToken->nodeType = AST_type::OP_LE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EQUAL:
    openToken->nodeType = AST_type::OP_EQUAL ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_NOT_EQUAL:
    openToken->nodeType = AST_type::OP_NOT_EQUAL ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_AND:
    openToken->nodeType = AST_type::OP_AND ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EXOR:
    openToken->nodeType = AST_type::OP_EXOR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_OR:
    openToken->nodeType = AST_type::OP_OR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LOGICAL_AND:
    openToken->nodeType = AST_type::OP_LOGICAL_AND;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LOGICAL_OR:
    openToken->nodeType = AST_type::OP_LOGICAL_OR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_ASSIGN:
    openToken->nodeType = AST_type::OP_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_TIMES_ASSIGN:
    openToken->nodeType = AST_type::OP_TIMES_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DIVIDE_ASSIGN:
    openToken->nodeType = AST_type::OP_DIVIDE_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MODULUS_ASSIGN:
    openToken->nodeType = AST_type::OP_MODULUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_PLUS_ASSIGN:
    openToken->nodeType = AST_type::OP_PLUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MINUS_ASSIGN:
    openToken->nodeType = AST_type::OP_MINUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_LEFT_ASSIGN:
    openToken->nodeType = AST_type::OP_SHIFT_LEFT_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_RIGHT_ASSIGN:
    openToken->nodeType = AST_type::OP_SHIFT_RIGHT_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_AND_ASSIGN:
    openToken->nodeType = AST_type::OP_AND_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_OR_ASSIGN:
    openToken->nodeType = AST_type::OP_OR_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EXOR_ASSIGN:
    openToken->nodeType = AST_type::OP_EXOR_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_COMMA:
    openToken->nodeType = AST_type::OP_COMMA ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_QUESTION:
    openToken->nodeType = AST_type::OP_TERTIARY ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DOT:
    openToken->nodeType = AST_type::OP_DOT ;
    return AST_type::ASTP(openToken) ;
  default:
    pushToken(openToken) ;
    return AST_type::ASTP(0) ;
  }
}

bool checkUnaryToken(AST_type::elementType e) {
  return (e == AST_type::TK_PLUS ||
	  e == AST_type::TK_MINUS ||
	  e == AST_type::TK_NOT ||
	  e == AST_type::TK_AND ||
	  e == AST_type::TK_TIMES ||
	  e == AST_type::TK_INCREMENT ||
	  e == AST_type::TK_DECREMENT) ;
}

AST_type::elementType unaryOperator(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_PLUS:
    return AST_type::OP_UNARY_PLUS ;
  case AST_type::TK_MINUS:
    return AST_type::OP_UNARY_MINUS ;
  case AST_type::TK_NOT:
    return AST_type::OP_NOT ;
  case AST_type::TK_AND:
    return AST_type::OP_AMPERSAND ;
  case AST_type::TK_TIMES:
    return AST_type::OP_STAR ;
  case AST_type::TK_INCREMENT:
    return AST_type::OP_INCREMENT ;
  case AST_type::TK_DECREMENT:
    return AST_type::OP_DECREMENT ;
  default:
    return AST_type::OP_ERROR ;
  }
}

bool checkPostFixToken(AST_type::elementType e) {
  return (e == AST_type::TK_INCREMENT ||
	  e == AST_type::TK_DECREMENT) ;
}
  
AST_type::elementType postFixOperator(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_INCREMENT:
    return AST_type::OP_POSTINCREMENT ;
  case AST_type::TK_DECREMENT:
    return AST_type::OP_POSTDECREMENT ;
  default:
    return AST_type::OP_ERROR ;
  }
  return AST_type::OP_ERROR ;
}
extern AST_type::ASTP parseExpression(std::istream &is, int &linecount) ;

AST_type::ASTP applyPostFixOperator(AST_type::ASTP expr,
				    std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  if(checkPostFixToken(openToken->nodeType)) {
    CPTR<AST_exprOper> post = new AST_exprOper ;
    post->nodeType = postFixOperator(openToken->nodeType) ;
    post->terms.push_back(expr) ;
    return AST_type::ASTP(post) ;
  }
  if(openToken->nodeType == AST_type::TK_OPENBRACKET) {
    AST_type::ASTP index = parseExpression(is,linecount) ;
    openToken = getToken(is,linecount) ;
    if(openToken->nodeType != AST_type::TK_CLOSEBRACKET) {
      pushToken(openToken) ;
      cerr << "syntax error line " << linecount << endl ;
    }
    CPTR<AST_exprOper> array = new AST_exprOper ;
    array->nodeType = AST_type::OP_ARRAY ;
    array->terms.push_back(expr) ;
    array->terms.push_back(index) ;
    return AST_type::ASTP(array) ;
  }
  pushToken(openToken) ;
    
  return expr ;
}


AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  if(openToken->nodeType == AST_type::TK_OPENPAREN) {
    AST_type::ASTP exp = parseExpression(is,linecount) ;
    CPTR<AST_Token> closeToken = getToken(is,linecount) ;
    CPTR<AST_exprOper> group = new AST_exprOper ;
    group->nodeType = AST_type::OP_GROUP ;
    group->terms.push_back(AST_type::ASTP(exp)) ;
    if(closeToken->nodeType != AST_type::TK_CLOSEPAREN) {
      pushToken(closeToken) ;
      group->nodeType=AST_type::OP_GROUP_ERROR ;
    }
    return applyPostFixOperator(AST_type::ASTP(group),is,linecount) ;
  }
  if(checkUnaryToken(openToken->nodeType)) {
    //check for unary operators
    AST_type::ASTP expr = parseExpressionPartial(is,linecount) ;
    if(expr!= 0) {
      CPTR<AST_exprOper> unary = new AST_exprOper ;
      unary->nodeType = unaryOperator(openToken->nodeType) ;
      unary->terms.push_back(expr) ;
      return AST_type::ASTP(unary) ;
    }
  }
  if(isTerm(openToken->nodeType)) {
    pushToken(openToken) ;
    AST_type::ASTP exp = parseTerm(is,linecount);

    openToken = getToken(is,linecount) ;
    //    if(checkPostFixToken(openToken->nodeType)) {
    //      CPTR<AST_exprOper> post = new AST_exprOper ;
    //      post->nodeType = postFixOperator(openToken->nodeType) ;
    //      post->terms.push_back(AST_type::ASTP(exp)) ;
    //      return AST_type::ASTP(post) ;
    //    }
    if(openToken->nodeType == AST_type::TK_OPENPAREN) {
      // Function
      AST_type::ASTP args = parseExpression(is,linecount) ;
      CPTR<AST_Token> closeToken = getToken(is,linecount) ;
      CPTR<AST_exprOper> func = new AST_exprOper ;
      func->nodeType = AST_type::OP_FUNC ;
      func->terms.push_back(AST_type::ASTP(exp)) ;
      func->terms.push_back(args) ;
      if(closeToken->nodeType != AST_type::TK_CLOSEPAREN) {
	CPTR<AST_Token> err = new AST_Token ;
	*err = *closeToken ;
	err->nodeType = AST_type::OP_ERROR ;
	func->terms.push_back(AST_type::ASTP(err)) ;
	pushToken(closeToken) ;
      }
      return applyPostFixOperator(AST_type::ASTP(func),is,linecount) ;
    }
    pushToken(openToken) ;
    return applyPostFixOperator(exp,is,linecount) ;
    //    return exp ;
  }
  pushToken(openToken) ;
  return 0 ;
}
  


AST_type::ASTP parseExpression(std::istream &is, int &linecount) {
  AST_type::ASTP expr = parseExpressionPartial(is,linecount) ;
  if(expr == 0) // If no valid expression then return null
    return expr ;
  vector<CPTR<AST_exprOper> > exprStack ;
  {
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = AST_type::OP_NIL ;
    tmp->terms.push_back(expr) ;
    exprStack.push_back(tmp) ;
  }
  const unsigned int mask = ~0x7f ;
  // After getting the first term we are in a loop of searching for operators
  do {
    // binary operator check
    AST_type::ASTP op = parseOperator(is, linecount) ;
    if(op == 0)
      break ;
    
    expr = parseExpressionPartial(is,linecount) ;
    
    if(expr == 0) {
      cerr << "expecting expression after binary operator" << endl ;
    }
    
    if(exprStack.back()->nodeType == AST_type::OP_NIL) {
      // If no operator is parsed yet, we just get started and initilize
      // the left branch
      exprStack.back()->nodeType = op->nodeType ;
      exprStack.back()->terms.push_back(expr) ;
    } else {
      // Now we reorder the tree based on operator precedence
      while(exprStack.size() >1 &&
	    ((op->nodeType&mask) >= (mask&exprStack.back()->nodeType))) {
	exprStack.pop_back() ;
      }
      if(op->nodeType == exprStack.back()->nodeType) {
	// If operator is the same, just chain the terms
	exprStack.back()->terms.push_back(expr) ;
      } else if(((op->nodeType)&mask) < ((exprStack.back()->nodeType)&mask)) {
        // if operator is lower precedence
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
	np->terms.push_back(exprStack.back()->terms.back()) ;
	np->terms.push_back(expr) ;
	exprStack.back()->terms.back() = AST_type::ASTP(np) ;
	if(op->nodeType == AST_type::OP_TERTIARY ) {
	  CPTR<AST_Token> op2 = getToken(is,linecount) ;
	  if(op2->nodeType == AST_type::TK_COLON) {
	    expr = parseExpressionPartial(is,linecount) ;
	    np->terms.push_back(expr) ;
	  } else {
	    pushToken(op2) ;
	    cerr << "syntax error parsing tertiary operator" << endl ;
	  }
	}
	exprStack.push_back(np) ;
      } else {
	if(op->nodeType == AST_type::OP_TERTIARY) {
	  cerr << "unexpected TERTIARY operator!" << endl ;
	} 
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
	np->terms.push_back(AST_type::ASTP(exprStack.back())) ;
	np->terms.push_back(expr) ;
	exprStack.back() = np ;
      }
    }
  } while(true) ;

  return AST_type::ASTP(exprStack.front()) ;
}

AST_type::ASTP parseDeclaration(std::istream &is, int &linecount) ;
AST_type::ASTP parseStatement(std::istream &is, int &linecount) ;
AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount) ;
AST_type::ASTP parseIfStatement(std::istream &is, int &linecount) ;

AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount)  ;

AST_type::ASTP parseCaseStatement(std::istream &is, int &linecount) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_CASE &&
     token->nodeType != AST_type::TK_DEFAULT) {
    cerr << "internal error parsing switch statement on line " << linecount
	 << endl ;
  }
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::TK_CASE ;
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  if(token->nodeType == AST_type::TK_CASE ) {
    AST_type::ASTP expr = parseExpression(is,linecount) ;
    AST_data->elements.push_back(expr) ;
  }
  
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_COLON) {
    cerr << "syntax error on case statement, line = " << linecount << endl ;
  }
  
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
  CPTR<AST_switchStatement> sw = new AST_switchStatement ;
  sw->nodeType = AST_type::TK_SWITCH ;
  sw->statement = AST_type::ASTP(token) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENPAREN) {
    cerr << "syntax error, switch statement should have parenthesis around input" << endl ; ;
  }
  AST_type::ASTP conditional = parseExpression(is,linecount) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_CLOSEPAREN) {
    cerr << "syntax error, switch statement should have parenthesis around input" << endl ; ;
  }
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENBRACE) {
    cerr << "syntax error, switch statement missing open brace, line = " << linecount << endl ;
  }
  sw->conditional = conditional ;
  token = getToken(is,linecount) ;
  while(token->nodeType != AST_type::TK_CLOSEBRACE &&
	token->nodeType != AST_type::TK_ERROR) {
    if(token->nodeType == AST_type::TK_CASE ||
       token->nodeType == AST_type::TK_DEFAULT) {
      pushToken(token) ;
      sw->body.push_back(parseCaseStatement(is,linecount)) ;
    } else {
      pushToken(token) ;
      sw->body.push_back(parseStatement(is,linecount)) ;
    }
    token = getToken(is,linecount) ;
  }
  
  
  return AST_type::ASTP(sw) ;

}

AST_type::ASTP parseIfStatement(std::istream &is, int &linecount)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_IF) {
    cerr << "unexpected error in parseIfStatement" << endl ;
    return AST_type::ASTP(token) ;
  }
  AST_type::ASTP iftok = AST_type::ASTP(token) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENPAREN) {
    cerr << "syntax error, line=" << linecount << " parsing if statement" << endl ;
    return AST_type::ASTP(token) ;
  }

  AST_type::ASTP conditional = parseExpression(is,linecount) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_CLOSEPAREN) {
    cerr << "syntax error, line=" << linecount << " missing closing paren" << endl ;
    return AST_type::ASTP(token) ;
  }
  
  AST_type::ASTP body = parseStatement(is,linecount) ;

  AST_type::ASTP ebody = 0 ;
  token = getToken(is,linecount) ;
  
  if(token->nodeType == AST_type::TK_ELSE) {
    ebody = parseStatement(is,linecount) ;
  } else {
    pushToken(token) ;
  }
  return AST_type::ASTP(new AST_ifStatement(iftok, conditional,body,ebody)) ;
}			 

bool isTypeDecl(CPTR<AST_Token> p) {
  switch(p->nodeType) {
  case AST_type::TK_CHAR:
  case AST_type::TK_FLOAT:
  case AST_type::TK_DOUBLE:
  case AST_type::TK_INT:
  case AST_type::TK_BOOL:
  case AST_type::TK_LONG:
  case AST_type::TK_SIGNED:
  case AST_type::TK_UNSIGNED:
  case AST_type::TK_CONST:
    return true ;
  default:
    return false ;
  }
}


AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
  AST_type::ASTP loop = AST_type::ASTP(token) ;
  switch(token->nodeType) {
  case AST_type::TK_FOR:
    {
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	cerr << "syntax error, expecting '(' in for loop on line" << linecount << endl ;
	return loop ;
      }
      token = getToken(is,linecount) ;
      AST_type::ASTP initializer = 0 ;
      if(isTypeDecl(token)) {
	pushToken(token) ;
	initializer = parseDeclaration(is,linecount) ;
      } else {
	pushToken(token) ;
	initializer = parseExpression(is,linecount) ;
      }
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
	cerr << "syntax error after initializer in for loop, line " << linecount << endl ;
      }
      AST_type::ASTP conditional = parseExpression(is,linecount) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
	cerr << "syntax error after conditional in for loop, line " << linecount << endl ;
	return loop ;
      }
      AST_type::ASTP advance = parseExpression(is,linecount) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_CLOSEPAREN) {
	cerr << "syntax error findng close paren in for loop, line " << linecount << endl ;
	return loop ;
      }
      AST_type::ASTP body = parseStatement(is,linecount) ;
      return AST_type::ASTP(new AST_loopStatement(loop,initializer,conditional, advance, body)) ;
    }
    break ;
  case AST_type::TK_WHILE:
    {
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	cerr << "syntax error, expecting '(' in for loop on line" << linecount << endl ;
	return loop ;
      }
      AST_type::ASTP conditional = parseExpression(is,linecount) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_CLOSEPAREN) {
	cerr << "syntax error findng close paren in for loop, line " << linecount << endl ;
	return loop ;
      }
      AST_type::ASTP body = parseStatement(is,linecount) ;
      return AST_type::ASTP(new AST_loopStatement(loop,conditional, body)) ;
    }
    break ;
  case AST_type::TK_DO:
    {
      AST_type::ASTP body = parseStatement(is,linecount) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_WHILE) {
	cerr << "syntax error in do loop, expecting while, line " << linecount << endl ;
	return loop ;
      }
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	cerr << "syntax error in do loop, expecting '(', line " << linecount << endl ;
	return loop ;
      }
      AST_type::ASTP conditional = parseExpression(is,linecount) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_CLOSEPAREN) {
	cerr << "syntax error in do loop, expecting ')', line " << linecount << endl ;
      }
      token = getToken(is,linecount) ;
      
      if(token->nodeType != AST_type::TK_SEMICOLON) {
	cerr << "syntax error in do loop, expecting ';', line " << linecount << endl ;
      }
      
      return AST_type::ASTP(new AST_loopStatement(loop,conditional, body)) ;
    }
    break ;
  default:
    return loop ;
  }
  
  return AST_type::ASTP(getToken(is,linecount)) ;
}

AST_type::ASTP parseType(std::istream &is, int &linecount) {
  //  cout << "parsing type" << endl ;
  
  CPTR<AST_typeDecl> AST_data = new AST_typeDecl ;

  CPTR<AST_Token> token = getToken(is,linecount) ;
  while(isTypeDecl(token)) {
    AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
    token = getToken(is,linecount) ;
  }
  pushToken(token) ;
  
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseDeclaration(std::istream &is, int &linecount) {

  //  cout << "parsing declaration" << endl ;
  
  CPTR<AST_declaration> AST_data = new AST_declaration ;

  AST_data->type_decl = parseType(is,linecount) ;
  AST_data->decls = parseExpression(is,linecount) ;
  
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSpecialControlStatement(std::istream &is, int &linecount) {
  
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::OP_SPECIAL ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_SEMICOLON) {
    if(AST_data->elements.back()->nodeType == AST_type::TK_RETURN) {
      pushToken(token) ;
      AST_type::ASTP exp = parseExpression(is,linecount) ;
      AST_data->elements.push_back(exp) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON)
	cerr << "expected semicolon on return statement, line=" << linecount
	     << endl ;
    } else
      cerr << "syntax error near line " << linecount << endl ;
  }
  AST_data->elements.push_back(AST_type::ASTP(token)) ;

  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseStatement(std::istream &is, int &linecount) {
  CPTR<AST_Token> firstToken = getToken(is,linecount) ;
  
  switch(firstToken->nodeType) {
  case AST_type::TK_OPENBRACE:
  case AST_type::TK_OPENPAREN:
    pushToken(firstToken) ;
    return parseBlock(is,linecount) ;
  case AST_type::TK_CHAR:
  case AST_type::TK_FLOAT:
  case AST_type::TK_DOUBLE:
  case AST_type::TK_INT:
  case AST_type::TK_BOOL:
  case AST_type::TK_LONG:
  case AST_type::TK_SIGNED:
  case AST_type::TK_UNSIGNED:
  case AST_type::TK_CONST:
    pushToken(firstToken) ;
    return parseDeclaration(is,linecount) ;
  case AST_type::TK_FOR:
  case AST_type::TK_WHILE:
  case AST_type::TK_DO:
    pushToken(firstToken) ;
    return parseLoopStatement(is,linecount) ;
  case AST_type::TK_IF:
    pushToken(firstToken) ;
    return parseIfStatement(is,linecount) ;
  case AST_type::TK_SWITCH:
    pushToken(firstToken) ;
    return parseSwitchStatement(is,linecount) ;
  case AST_type::TK_BREAK:
  case AST_type::TK_CONTINUE:
  case AST_type::TK_RETURN:
    pushToken(firstToken) ;
    return parseSpecialControlStatement(is,linecount) ;
  case AST_type::TK_SEMICOLON:
    return AST_type::ASTP(firstToken) ;
  case AST_type::TK_NAME:
    {
      pushToken(firstToken) ;
      AST_type::ASTP exp = parseExpression(is,linecount) ;
      AST_type::ASTP term = AST_type::ASTP(getToken(is,linecount)) ;
      if(term->nodeType != AST_type::TK_SEMICOLON) {
	cerr << "syntax error, linecount = " << linecount << endl ;
      }
      AST_type::ASTP stat = new AST_SimpleStatement(exp,term) ;
      return stat ;
    }
  
  default:
    break ;
  }
  
  return AST_type::ASTP(firstToken) ;
}

AST_type::ASTP parseBlock(std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;

  AST_type::elementType closeType = AST_type::TK_CLOSEBRACE ;
  switch(openToken->nodeType) {
  case AST_type::TK_OPENBRACE:
    closeType = AST_type::TK_CLOSEBRACE ;
    break ;
  case AST_type::TK_OPENBRACKET:
    closeType = AST_type::TK_CLOSEBRACKET ;
    break ;
  case AST_type::TK_OPENPAREN:
    closeType = AST_type::TK_CLOSEPAREN ;
    break ;
  default:
    return AST_type::ASTP(openToken) ;
  }

  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::TK_BRACEBLOCK ;
  AST_data->elements.push_back(AST_type::ASTP(openToken)) ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  while(token->nodeType != closeType) {
    pushToken(token) ;
    CPTR<AST_type> statement = parseStatement(is,linecount) ;
    AST_data->elements.push_back(statement) ;
    token = getToken(is,linecount) ;
    if(is.fail() || is.eof()) 
      break ;
  }
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(AST_data) ;
}

void parseFile::setup_Test(std::ostream &outputFile) {
  CPTR<AST_Token> token = getToken(is,line_no) ;
  if(token->nodeType == AST_type::TK_OPENBRACE) {
    pushToken(token) ;
    CPTR<AST_type> ap = parseBlock(is,line_no) ;
    //    outputFile << "Parsed TEST:" << endl ;
    int lineno=-1 ;
    ap->DiagPrint(outputFile,lineno) ;
  }
}

void parseFile::setup_Rule(std::ostream &outputFile) {
  killsp() ;
  string rule_type ;
  if(is_name(is)) {
    rule_type = get_name(is) ;
  } else 
    throw parseError("syntax error") ;
  nestedparenstuff signature ;

  signature.get(is) ;
  line_no += signature.num_lines() ;
  nestedbracketstuff apply_op ;
  killsp() ;
  if(rule_type == "apply") {
    if(is.peek() != '[') 
      throw parseError("apply rule missing '[operator]'") ;
    apply_op.get(is) ;
    line_no += apply_op.num_lines() ;
    killsp() ;
  }
  

  string constraint, conditional ;
  string parametric_var ;
  list<string> options ;
  list<string> comments ;
  list<pair<variable,variable> > inplace ;
  
  bool use_prelude = false ;
  bool is_specialized = false ;
  while(is.peek() == ',') {
    is.get() ;
    killsp() ;
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
      line_no += con.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      if(parametric_var != "") {
        throw parseError("syntax error: canot specify more than one parametric variable") ;
      }
        
      parametric_var = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
        throw parseError("syntax error: canot specify more than one conditional variable") ;
      }
      conditional = con.str() ;
      line_no += con.num_lines() ;
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
        cerr << filename << ':' << line_no << ":0: warning: type of conditional variable '" << v << "' not found!"  << endl  ;
      } else {
        //        cout << "mi->first=" << mi->first << endl ;
        //        cout << "mi->second.first=" << mi->second.first
        //             << "mi->second.second='" << mi->second.second <<"'"<< endl ;
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
      line_no += ip.num_lines() ;
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
      killsp() ;
      continue ;
    } else if(s == "specialized") {
      is_specialized = true ;
    } else if(s == "option") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    killsp() ;
  }

  string sig = signature.str() ;
  string heads,bodys ;
  exprP head=0,body=0 ;
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
  
  string class_name = "file_" ;
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

  //  timeb tdata ;
  //  ftime(&tdata) ;
  timespec tdata ;
  clock_gettime(CLOCK_MONOTONIC,&tdata) ;
  
  
  ostringstream tss ;
  //  tss <<  '_' << tdata.time << 'm'<< tdata.millitm;
  tss << '_' << tdata.tv_sec << 'm' << tdata.tv_nsec/1000000 ;
  
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
  
  set<vmap_info> sources ;
  set<vmap_info> targets ;
  if(body != 0)
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;

  set<vmap_info>::const_iterator i ;
  variableSet input,output ;
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

  set<std::list<variable> > validate_set ;
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

  

  map<variable,string> vnames ;
  variableSet::const_iterator vi ;
  variableSet all_vars = input;
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


  variableSet checkset ;
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
  
  map<variable,pair<string,string> > local_type_map ;
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

  if(!prettyOutput)
    outputFile << "namespace {" ;
  outputFile << "class " << class_name << " : public Loci::" << rule_type << "_rule" ;
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
    
              
  bool singletonApply = false ;
  if(rule_type == "apply") {
    if(output.size() != 1) 
      throw parseError("apply rule should have only one output variable") ;
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
    outputFile << "> " ;
  }
  outputFile << " {" << endl ;
  syncFile(outputFile) ;

  variableSet outs = output ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  variableSet ins = input ;
  ins -= outs ;
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
    syncFile(outputFile) ;
  }
  bool output_param = false ;
  for(vi=outs.begin();vi!=outs.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(vi->get_info().name != "OUTPUT" && mi->second.first == "param") {
      output_param= true ;
    }
    if(!prettyOutput)
      outputFile << "    Loci::" << mi->second.first <<  mi->second.second ;
    else
      outputFile << "    " << mi->second.first <<  mi->second.second ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "public:" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    " << class_name << "() {" << endl ;
  syncFile(outputFile) ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
               << vnames[*vi] << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  for(i=targets.begin();i!=targets.end();++i) {
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
      for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
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
    syncFile(outputFile) ;
  }
  //  outputFile <<   "       output(\"" << heads << "\") ;" << endl ;
  //  syncFile(outputFile) ;

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

        cerr << filename << ':' << line_no << ":0: warning: type of constraint variable '" << v << "' not found!"  << endl  ;

      } 
    }
    
    

    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  if(parametric_var != "") {
    outputFile <<   "       set_parametric_variable(\""
               << parametric_var << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(is_specialized) {
    outputFile <<   "       set_specialized() ; " << endl ;
    syncFile(outputFile) ;
  }
  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    syncFile(outputFile) ;
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
    syncFile(outputFile) ;
  }
  for(lsi=comments.begin();lsi!=comments.end();++lsi) {
    outputFile <<   "       comments(" << *lsi << ") ;" << endl ;
    syncFile(outputFile) ;
  }

  outputFile <<   "    }" << endl ;
  syncFile(outputFile) ;

  bool use_compute = true ;

  if(use_prelude) {
    process_Prelude(outputFile,vnames) ;
    killsp() ;
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
    killsp() ;
  }

  
  if(use_compute && is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;

  bool sized_outputs = false;
  variableSet outsmi = outs ;
  outsmi -= input ;
  for(vi=outsmi.begin();vi!=outsmi.end();++vi) {
    string ot = local_type_map[*vi].first ;
    if(ot == "storeVec" || ot == "storeMat" || ot == "multiStore")
      sized_outputs = true ;
  }

    
    

  if(rule_type == "singleton" ||
     rule_type == "optional"  ||
     rule_type == "default" ||
     rule_type == "constraint" ||
     (output_param && rule_type != "apply" ) ) {
    if(use_prelude) {
      string error = "inappropriate prelude on " + rule_type + " rule." ;
      throw parseError(error) ;
    }
    process_Compute(outputFile,vnames) ;
  } else {
    if(use_compute)
      process_Calculate(outputFile,vnames,validate_set) ;

    outputFile <<   "    void compute(const Loci::sequence &seq) { " << endl ;
    syncFile(outputFile) ;
//     if(use_prelude) {
//       outputFile <<   "      prelude(seq) ;" << endl ;
//       syncFile(outputFile) ;
//     }
    if(use_compute) {
      if(singletonApply) {
        cerr << "NOTE: parameter only apply rule on '" << output << "' now executes single instance." << endl ;
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

  if(!prettyOutput)
    outputFile << "Loci::register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  else
    outputFile << "register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  syncFile(outputFile) ;

  if(!prettyOutput) {
    outputFile << "}" << endl ;
    syncFile(outputFile) ;
  }

  if(!use_prelude && sized_outputs && (rule_type != "apply")) 
    throw parseError("need prelude to size output type!") ;
}

void parseFile::processFile(string file, ostream &outputFile) {
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
  syncFile(outputFile) ;
  do {
    while(is.peek() == ' ' || is.peek() == '\t') {
      is.get(c) ;
      outputFile << c ;
    }
    try {
      if(is.peek() == '$') { // Loci specific code!
        is.get(c) ; // get the $
        if(is.peek() == '[') {
          map<variable,string> vnames ;
          int openbrace = 0 ;
          process_SpecialCommand(outputFile,vnames,openbrace) ;
        } else  if(is_name(is)) {
          std::string key = get_name(is) ;
          if(key == "type") {
            setup_Type(outputFile) ;
          } else if(key == "rule") {
            setup_Rule(outputFile) ;
	  } else if(key == "test") {
	    setup_Test(outputFile) ;
          } else if(key == "include") {
            killsp() ;
            if(!is_string(is))
              throw parseError("syntax error") ;
            string newfile = get_string(is) ;
            
            parseFile parser ;
            parser.processFile(newfile,outputFile) ;
            syncFile(outputFile) ;
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
	    killCommentOut(is,line_no,outputFile) ;
	    foundComment = true ;
	    break ;
	  }

          is.get(c) ;
          outputFile << c ;
        }
	if(!foundComment) {
	  is.get(c) ;
	  outputFile << endl ;
	  line_no++ ;
	}
      }
      if(is.peek() == EOF)
        is.get(c) ;
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

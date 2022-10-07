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
#ifndef PARSEBASE_H
#define PARSEBASE_H

#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>

using std::istream;
using std::string;
using std::ostream;


istream &killsp(istream &s, int &lines);
istream &killspOstr(istream &s, int &lines, string &str);

bool is_token(istream &s, const string &token);
bool get_token(istream &s, const string &token);

bool is_name(istream &s);
string get_name(istream &s);

bool is_comment(istream &s);
istream &killComment(istream &s, int & lines);
istream &killCommentOut(istream &s, int & lines,ostream &out);

extern bool prettyOutput ;

inline void syncFile(int line_no, const string& filename, std::ostream &outputFile) {
  if(!prettyOutput)
    outputFile << "#line " << line_no << " \"" << filename << "\"" << std::endl ;
}

struct parseError {
  std::string error_type ;
parseError(std::string errs) : error_type(errs) {}
} ;


class parsebase {
 public:
  int lines ;
  parsebase() {lines = 0 ; }
  
  istream &killsp(istream &s) {
    ::killsp(s,lines) ;
    return s ;
  }

  // added this function, for nestedbracestuff
  istream &killspOstr(istream &s, string& str) {
    ::killspOstr(s,lines, str) ;
    return s ;
  }
} ;

//read in everything between ( and )
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
	    throw parseError("unexpected EOF parsing string") ;
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

//read in everything between [ and ] 
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

//read in everything between { and } , 
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


#endif

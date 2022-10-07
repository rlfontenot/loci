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
#include "parsebase.h"
#include <iostream>
#include <sstream>

using std::ostream;

bool  is_token(istream &s, const string &token) {
  const int sz = token.size() ;
  for(int i=0;i<sz;++i) {
    if(s.peek() != token[i]) {
      for(--i;i>=0;--i)
	s.putback(token[i]) ;
      return false ;
    }
    s.get() ;
  }
  for(int i=token.size()-1;i>=0;--i) 
    s.putback(token[i]) ;
  return true ;
}
    
bool get_token(istream &s, const string &token) {
  const int sz = token.size() ;
  for(int i=0;i<sz;++i) {
    if(s.peek() != token[i]) {
      for(--i;i>=0;--i)
	s.putback(token[i]) ;
      return false ;
    }
    s.get() ;
  }
  return true ;
}
      
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

/* remove ' '\t\r\n and comments */
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

/*read in comments and send to out*/
istream &killCommentOut(istream &s, int & lines, ostream &out) {
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



/* added this function, for nestedbracestuff
   read in ' '\t\r\n comments, send them to str
*/
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


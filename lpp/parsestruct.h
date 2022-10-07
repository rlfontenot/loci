

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
#ifndef PARSESTRUCT_H
#define PARSESTRUCT_H

#include "parsebase.h"
#include <string>
#include <vector>
#include <iostream>

using std::list;
using std::string;
using std::istream;
using std::vector;
using std::cerr;
using std::endl;

/*read in everything between ( and ) ,
  compared with class nestedparenstuff,
  this class won't remove white space.
  used in struct for funcion signature 
*/
class nestedparenstuff_func : public parsebase {
public:
  string paren_contents ;
  istream &get(istream &s);
 
  string str() {
    return paren_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;


//----------------------------------------------------------------

/*this class is added to parse a struct as an object
  namespace :: not considered yet
  decorators such as const, & * not considered yet
*/
class parsestruct : public parsebase {
private:
  //all types allowed in $struct definition, use vector here instead of set so that is_type() can return which type it is
  static  vector<string> type_list;

  int is_type(istream &s) { //
    for(size_t i = 0; i< type_list.size(); i++){
      if(is_token(s, type_list[i])){
	return i;
      }
    }
    return -1;
  }

  void add_type(const string& t){
    for(size_t i= 0; i< type_list.size(); i++){
      if(type_list[i]==t){
	std::cerr << " ERROR: data type " <<  t << " already exist"<< endl;
	return ;
      }
    }
    type_list.push_back(t);
  }
  
  struct funct{
    string ret_type; //return type
    string name; //function name
    string op; //operator
    nestedparenstuff_func sig; //argument list
    nestedbracestuff body; //function body
  };

  struct constr{
    string init;//initializer
    nestedparenstuff_func sig; //argument list
    nestedbracestuff body; //function body
  };

 
public:
 
  string struct_name ;
  list<string> data_type;
  list<string> data_name;
  list<funct> funct_list;
  list<constr> constr_list;
  int first_line, data_line, last_line;;

  //read in
  istream &get(std::string& filename, int& cnt, int start_line,  istream &is); 

  //output
  string str();
  
};


  
#endif

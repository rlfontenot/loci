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
  to define a struct in .loci file
  1. start with: $struct, followed by: struct_name {
  2. end with };
  3. inside {}, data type followed by data name to defined data members;
  4. constructors
  5. member functions
  6. overloading operators as member functions

  current restrctions:
  1. data type is the build-in types such as double, float, int
  and  previous use-defined types using $struct
   
  2. no friend functions

  3. inside {}, c++ keyword such as const, static, inline ... not considered

  4. namespace not considered yet

  5. the decorater after a type, such as &, * not considered yet

  things done:
  parse the struct definition into an object

  things to do:
  modify str() function, identify special token such as types, and replace it with other tokens  

*/
#include "parsestruct.h"
using std::cout;

//what about void here? can only be return type, 
vector<string> parsestruct :: type_list ({"double", "int", "float", "void"});
/*
  split a string into a list of tokens,
  replace some tokens,
  and convert back to a string
*/
std::string process_string(std::string& s){ 
  list<string> tokens;
  string delim = " ";
  auto start = 0U;
  auto end = s.find(delim);
  while (end != std::string::npos)
    {
      std::cout << s.substr(start, end - start) << std::endl;
      tokens.push_back(s.substr(start, end - start));
      start = end + delim.length();
      end = s.find(delim, start);
    }
  
  std::cout << s.substr(start, end);
  tokens.push_back(s.substr(start, end - start));
  cout << "done split in" << endl;
  for(auto x:tokens)cout << x << " ";
  cout << endl;

  string out;
  return out;
    

}

  
istream &nestedparenstuff_func::get(istream &s) {

  parsebase::killspOstr(s, paren_contents) ;
  if(s.peek() != '(')
    throw parseError(string("syntax error, expecting '(' for function dec")) ;
  s.get() ;

  parsebase::killspOstr(s, paren_contents) ;
  int open_parens = 0 ;
  while(s.peek() != ')' || open_parens != 0) {
    if(s.peek() == EOF)
      throw parseError("unexpected EOF") ;
    if(s.peek() == '(')
      open_parens++ ;
    if(s.peek() == ')')
      open_parens-- ;
    if(s.peek() == '\n' || s.peek() == '\r') {
      paren_contents += s.get() ;
      lines++ ;
      continue;
    }
    paren_contents += s.get() ;
    parsebase::killspOstr(s, paren_contents) ;
  }
  paren_contents += s.get() ;
  return s ;
}


//read in struct
istream &parsestruct::get(std::string& filename, int& cnt, int start_line,  istream &is) {
  //this function also serve as initializer
  first_line = start_line;
  parsebase::killsp(is) ; //after $struct is white space
  if(is_name(is)) { //should be followed by the name of struct
    struct_name = get_name(is) ;
  } else
    throw parseError("syntax error, struct name expected") ;
  parsebase::killsp(is) ; //space after struct name
  if(is.peek() != '{')
    throw parseError("syntax error, expecting '{' after struct name") ;
  add_type(struct_name);
  int ch =  is.get(); //get '{'
  int openbrace = 1 ;
    
  parsebase::killsp(is) ; //space after '{'
  int type_index = -1;
  while(true ){ //parse the member data and member function

    //first, it should start with a type,
    //might be decorator such as 'const', 'inline', consider this later
    type_index = is_type(is);
    if(type_index <0 || type_index > type_list.size()) break; //if not start with a type, break

    //get type
    if(!get_token(is, type_list[type_index])) throw parseError("syntax error") ;

    parsebase::killsp(is) ;
     
    //after type, it can be a name(function name or data_name)
    if(is_name(is)){
      string s = get_name(is) ;
      parsebase::killsp(is) ;
      string op;
      if(s == "operator"){ //operator definition
	if(is.peek()=='('){ //function call operator
	  op += is.get();
	  parsebase::killsp(is) ;
	  if(is.peek() != ')')throw parseError("syntax error: ')' expected") ;
	  op += is.get();
	}else{ //other operator
	  while(is.good() && is.peek() != ' ' && is.peek() != '(') op += is.get();
	}
	//if(!is_operator(op)) throw parseError("syntax error, expecting operator") ;
	parsebase::killsp(is) ;
      }
	
      if( is.peek() != '('){ //member data declaration
	data_type.push_back(type_list[type_index]);
	data_name.push_back(s);
	while(is.peek() == ',') {
	  is.get() ;
	  parsebase::killsp(is) ;
	  if(!is_name(is))
	    throw parseError("syntax error") ;
	    
	  string s = get_name(is) ;
	  data_type.push_back(type_list[type_index]);
	  data_name.push_back(s);
	}
	parsebase::killsp(is) ;
	if(is.peek() != ';')
	  throw parseError("syntax error, expecting ';' ") ;
	is.get();
	parsebase::killsp(is) ;
      }else{//this is function definition
	funct afunct;
	afunct.name = s;
	afunct.op = op;
	afunct.ret_type = type_list[type_index];
	afunct.sig.get(is) ;
	lines += afunct.sig.num_lines() ; 
	parsebase::killsp(is) ;
	//if no '{', then the function definition is outside, should end with ';' 
	if(is.peek() == '{'){
	  afunct.body.get(is);
	  lines += afunct.body.num_lines() ;
	  funct_list.push_back(afunct);
	}else if(is.peek()==';'){
	  funct_list.push_back(afunct);
	}else{
	  throw parseError("syntax error, expecting ';' or '{' ") ; 
	}
	//this is constructor definition
	parsebase::killsp(is) ;
	
      }
    }else if(is.peek()=='('){ //after type, it can also be '(' for constructors

      constr aconstr;

      aconstr.sig.get(is) ;
      lines += aconstr.sig.num_lines() ;
      parsebase::killsp(is) ;
      //if no '{', then the function definition is outside, should end with ';' 
      if(is.peek() == '{'){
	aconstr.body.get(is);
	lines += aconstr.body.num_lines() ;
	constr_list.push_back(aconstr);
      }else if(is.peek()==';'){
	constr_list.push_back(aconstr);
      }else if(is.peek()==':'){
	is.get();
	parsebase::killsp(is) ;
	while(is.peek() != '{'&& !is.eof()){
	  aconstr.init += is.get();
	}
	parsebase::killsp(is) ;
	if(is.peek() != '{')throw parseError("syntax error, expecting '{' ") ;
	aconstr.body.get(is);
	lines += aconstr.body.num_lines() ;
	constr_list.push_back(aconstr);
	  
      }else{
	cerr << " actual :  " << is.peek() << endl;
	throw parseError("syntax error, expecting ';' or '{', and it is ") ; 
      }
      //this is constructor definition
      parsebase::killsp(is) ;
    }
  }
  parsebase::killsp(is) ;
  // $struct should end with "};"
  if(is.peek() != '}')throw parseError("syntax error, expecting '}' ") ;
  ch = is.get();
  openbrace--;
  parsebase::killsp(is) ;
    
  if(is.peek() != ';')throw parseError("syntax error, expecting '}' ") ;
  ch = is.get();
  parsebase::killsp(is) ;
  last_line = lines;
  return is;
}


/* output struct */
string parsestruct::str() {
  string s ;
  //declaration of struct
  s += "struct " + struct_name + " { \n" ;

  //output data members
  auto itr1 = data_type.begin();
  auto itr2 = data_name.begin();
  for ( ; itr1 != data_type.end() && itr2 != data_name.end(); itr1++, itr2++){
    s += *itr1 + " " + *itr2 +";\n";
  }

  //output constructors
  for(auto itr = constr_list.begin(); itr != constr_list.end(); itr++){
    string s1 = '(' + (itr->sig).str();
    string s2 = '{' + (itr->body).str();
    s += struct_name + s1;
    if(!(itr->init).empty()) s += ':' + itr->init;
    s +=  s2 +"\n";
  }

  //output functions( include operator overloading functions)
  for(auto itr = funct_list.begin(); itr != funct_list.end(); itr++){
    string s1 = itr->ret_type +' ' + itr->name + ' ' + itr->op + '(' + (itr->sig).str();
    string s2 = '{' + (itr->body).str();
    s += s1;
    s +=  s2 +"\n";
  }
    
  s += "}; \n" ;
  return s ;
}

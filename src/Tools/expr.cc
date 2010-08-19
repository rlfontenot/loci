//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#include <Tools/expr.h>
#include <Tools/parse.h>
#include <Tools/stream.h>
#include <Tools/except.h>
#include <Tools/tools.h>
#include <map> 
using std::map ;
#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <string>
using std::string ;

// Todo list:

// 0) add double constants (done (except number starting with a '.')
// 1) expression substitution (done)
// 2) Compiling expressions
// 3) update exception handling
// 4) add vector handling
// 5) add gradient and divergence operator evaluation


namespace Loci {
  using namespace expr ;

  // Compare two exprsesions.  Return 0 if equal, return -1 if e1 is cannonically before e2 and 1 if otherwise.
  int compare_expressions(exprP e1, exprP e2) {
    if(e1->op!=e2->op)
      return (e1->op < e2->op)?-1:1 ;
    switch(e1->op) {
    case OP_NAME:
    case OP_STRING:
      if(e1->name == e2->name)
	return 0 ;
      return (e1->name < e2->name)?-1:1 ;
    case OP_INT:
      if(e1->int_val == e2->int_val)
	return 0 ;
      return (e1->int_val < e2->int_val)?-1:1 ;
    case OP_DOUBLE:
      if(e1->real_val == e2->real_val)
	return 0 ;
      return (e1->real_val < e2->real_val)?-1:1 ;
    case OP_FUNC:
    case OP_ARRAY:
    case OP_NAME_BRACE:
    case OP_FUNC_BRACE:
      if(e1->name != e2->name) 
	return (e1->name < e2->name)?-1:1 ;
      break ;
    default:
      break ;
    }
    // Now all thats left is to check the expr_list

    exprList::const_iterator li1,li2 ;
    for(li1=e1->expr_list.begin(),li2=e2->expr_list.begin();
	li1!=e1->expr_list.end()&&li2!=e2->expr_list.end();++li1,++li2) {
      int cmp = compare_expressions(*li1,*li2) ;
      if(cmp != 0)
	return cmp ;
    }
    if(li1!=e1->expr_list.end())
      return 1 ; 
    if(li2!=e2->expr_list.end())
      return -1 ;
    return 0 ;
  }


  using std::cout ;    
  void expression::PrintOperation(ostream &s, std::string oper,
				  char poChar, char pcChar) const
  {
    s << poChar ;
    //    if(expr_list.size() == 1) 
    //      s << "??" <<oper<<"??" ;
    
    if(!expr_list.empty()) {
      exprList::const_iterator i = expr_list.begin() ;
      (*i)->Print(s) ;
      ++i ;
      for(;i!=expr_list.end();++i) {
	s<<oper ;
	(*i)->Print(s) ;
      }
    }
    s << pcChar ;
  }

  void expression::Print(ostream &s) const
  {
    switch(op) {
    case OP_AT:
      PrintOperation(s,"@") ;
      break;
		
    case OP_NAME:
      s << name ;
      break ;
    case OP_FUNC:
      s << name ;
      PrintOperation(s,",",'(',')') ;
      break ;
      
    case OP_ARRAY:
      s << name ;
      PrintOperation(s,",",'[',']') ;
      break ;
      
    case OP_NAME_BRACE:
      s << name ;
      PrintOperation(s,",",'{','}') ;
      break ;
      
    case OP_FUNC_BRACE:
      { 
	s << name ;
	warn(expr_list.size() != 2) ;
	s << '(';
	exprList::const_iterator func = expr_list.begin() ;
	(*func)->Print(s) ;
	s << ')' ;
	
	exprList::const_iterator brace = ++func ;
	s << "{" ;
	(*brace)->Print(s) ;
	s << "}" ;
      }
      break ;
      
    case OP_STRING:
      s << '"' << name <<'"' ;
      break ;
      
    case OP_INT:
      s << int_val ;
      break ;
    case OP_DOUBLE:
      s.precision(14) ;
      s << real_val ;
      break ;
      
    case OP_PLUS:
      PrintOperation(s,"+") ;
      break ;
    case OP_MINUS:
      PrintOperation(s,"-") ;
      break ;
    case OP_TIMES:
      PrintOperation(s,"*") ;
      break ;
    case OP_DIVIDE:
      PrintOperation(s,"/") ;
      break ;
    case OP_AND:
      PrintOperation(s,"&") ;
      break ;
    case OP_OR:
      PrintOperation(s,"|") ;
      break ;
    case OP_MODULUS:
      PrintOperation(s,"%") ;
      break ;
    case OP_ASSIGN:
      PrintOperation(s,"=") ;
      break ;
    case OP_EXOR:
      PrintOperation(s,"^") ;
      break ;
    case OP_LOGICAL_AND:
      PrintOperation(s,"&&") ;
      break ;
    case OP_LOGICAL_OR:
      PrintOperation(s,"||") ;
      break ;
    case OP_COMMA:
      PrintOperation(s,",") ;
      break ;
    case OP_COLON:
      PrintOperation(s,":") ;
      break ;
    case OP_ARROW:
      PrintOperation(s,"->") ;
      break ;
    case OP_SHIFT_RIGHT:
      PrintOperation(s,">>") ;
      break ;
    case OP_SHIFT_LEFT:
      PrintOperation(s,"<<") ;
      break ;
    case OP_LE:
      PrintOperation(s,"<=") ;
      break ;
    case OP_LT:
      PrintOperation(s,"<") ;
      break ;
    case OP_GE:
      PrintOperation(s,">=") ;
      break ;
    case OP_GT:
      PrintOperation(s,">") ;
      break ;
    case OP_EQUAL:
      PrintOperation(s,"==") ;
      break ;
    case OP_NOT_EQUAL:
      PrintOperation(s,"!=") ;
      break ;

    case OP_SCOPE:
      PrintOperation(s,"::") ;
      break ;
        
    case OP_UNARY_PLUS:
      s << "+" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_UNARY_MINUS:
      s << "-" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_NOT:
      s << "!" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_TILDE:
      s << "~" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_AMPERSAND:
      s << "&" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_DOLLAR:
      s << "$" ;
      expr_list.front()->Print(s) ;
      break ;
    case OP_STAR:
      s << "*" ;
      expr_list.front()->Print(s) ;
      break ;

    case OP_ERROR:
      s << "*ERR*" ;
      break ;
        
    case OP_NIL:
      s <<"" ;
      break ;
    default:
      throw exprError("Undefined","unexpected operation in Print()",ERR_UNDEF) ;
      break ;
    }
  }
    
  OpType expression::get_oper(istream &s) {
    if(parse::get_token(s,"@")) {
      return OP_AT;
    }
    if(parse::get_token(s,"&&")) {
      return OP_LOGICAL_AND ;
    }
    if(parse::get_token(s,"||")) {
      return OP_LOGICAL_OR ;
    }
    if(parse::get_token(s,">>")) {
      return OP_SHIFT_RIGHT ;
    }
    if(parse::get_token(s,"<<")) {
      return OP_SHIFT_LEFT ;
    }
    if(parse::get_token(s,"->")) {
      return OP_ARROW ;
    }
    if(parse::get_token(s,"<=")) {
      return OP_LE ;
    }
    if(parse::get_token(s,">=")) {
      return OP_GE ;
    }
    if(parse::get_token(s,"==")) {
      return OP_EQUAL ;
    }
    if(parse::get_token(s,"!=")) {
      return OP_NOT_EQUAL ;
    }
    if(parse::get_token(s,"::")) {
      return OP_SCOPE ;
    }
    if(parse::get_token(s,"<")) {
      return OP_LT ;
    }
    if(parse::get_token(s,">")) {
      return OP_GT ;
    }
    
    if(parse::get_token(s,"%")) {
      return OP_MODULUS ;
    }
    if(parse::get_token(s,"+")) {
      return OP_PLUS ;
    } 
    
    if(parse::get_token(s,"-")) {
      return OP_MINUS ;
    }
    if(parse::get_token(s,"*")) {
      return OP_TIMES ;
    }
    if(parse::get_token(s,"/")) {
      return OP_DIVIDE ;
    }
    if(parse::get_token(s,"&")) {
      return OP_AND ;
    }
    if(parse::get_token(s,"|")) {
      return OP_OR ;
    }
    if(parse::get_token(s,"^")) {
      return OP_EXOR ;
    }
      
    if(parse::get_token(s,",")) {
      return OP_COMMA ;
    }
    if(parse::get_token(s,"=")) {
      return OP_ASSIGN ;
    }
    if(parse::get_token(s,":")) {
      return OP_COLON ;
    }
    return OP_ERROR ;
  }

  OpType expression::get_unary_oper(istream &s)
  {
    if(parse::get_token(s,"+"))
      return OP_UNARY_PLUS ;
    if(parse::get_token(s,"-"))
      return OP_UNARY_MINUS ;
    if(parse::get_token(s,"!"))
      return OP_NOT ;
    if(parse::get_token(s,"~"))
      return OP_TILDE ;
    if(parse::get_token(s,"&"))
      return OP_AMPERSAND ;
    if(parse::get_token(s,"$"))
      return OP_DOLLAR ;
    if(parse::get_token(s,"*"))
      return OP_STAR ;
    return OP_ERROR ;
  }

  exprP expression::create(istream &s)
  {
    return expression::create(s,';') ;
  }

  exprP expression::create(const string &s)
  {
    istringstream is(s) ;
    return expression::create(is,';') ;
  }

  exprP expression::expand_oper(istream &s, exprP &p)
  {
    exprList estack ;
    estack.push_back(p) ;
    OpType ot = expression::get_oper(s) ;
    const unsigned int mask = ~0x7f ;
    while(ot != OP_ERROR) {
      exprP p2 = expression::get_term(s) ;
      while(estack.size()>1 && ((ot&mask) >= (mask&(estack.back()->op)))) {
	estack.pop_back() ;
      }
      if(estack.back()->op == ot) {
	estack.back()->expr_list_priv.push_back(p2) ;
      } else if((mask&ot) < (mask&(estack.back()->op))) {
	exprP np = new expression() ;
	np->op_priv = ot ;
	np->expr_list_priv.push_back(estack.back()->expr_list.back()) ;
	np->expr_list_priv.push_back(p2) ;
	estack.back()->expr_list_priv.back() = np ;
	estack.push_back(np) ;
      } else {
	exprP np = new expression() ;
	np->op_priv = ot ;
	np->expr_list_priv.push_back(estack.back()) ;
	np->expr_list_priv.push_back(p2) ;
	estack.back() = np ;
      }
      ot = expression::get_oper(s) ;
    }
    return estack.front() ;
  }


  exprP expression::create(istream &s, char closing) {
    exprP p = new expression ;
    parse::kill_white_space(s) ;
    if(s.eof() || s.peek() == EOF) {
      p->op_priv = OP_NIL ;
      return p ;
    }
    if(s.peek() == EOF || s.peek() == closing) {
      s.get() ;
      p->op_priv = OP_NIL ;
      return p ;
    }
  
    p = expression::get_term(s) ;
    parse::kill_white_space(s) ;
    if(s.eof() || s.peek() == EOF) {
      return p ;
    }
  
    if(s.peek() == closing) {
      s.get() ;
      return p ;
    }
    exprP p1 = p ;
    p = new expression ;
    p->op_priv = expression::get_oper(s) ;
    if(p->op == OP_ERROR) {
      ostringstream oss ;
      string near ;
      s >> near ;
      oss << "Unable to determine operator, near text '" << near << "'" ;
      throw exprError("Syntax Error",oss.str(),ERR_SYNTAX) ;
      return p ;
    }
    p->expr_list_priv.push_back(p1) ;
    p->expr_list_priv.push_back(expression::get_term(s)) ;
    p = expression::expand_oper(s,p) ;
    parse::kill_white_space(s) ;
    if((s.eof() || s.peek() == EOF) && closing == ';')
      return p ;
    if(s.peek() == closing) {
      s.get() ;
      return p ;
    }
    p = new expression ;
    return p ;
  }

  exprP expression::get_name(istream &s) {
    if(parse::get_token(s,"(")) {
      return expression::create(s,')') ;
    }
    if(parse::is_name(s)) {
      exprP name = new expression ;
      name->op_priv = OP_NAME ;
      name->name_priv = parse::get_name(s) ;
      char closing = 0 ;
      if(parse::get_token(s,"(")) {
	name->op_priv = OP_FUNC ;
	closing = ')' ;
      } else if(parse::get_token(s,"[")) {
	name->op_priv = OP_ARRAY ;
	closing = ']' ;
      } else if(parse::get_token(s,"{")) {
	name->op_priv = OP_NAME_BRACE ;
	closing = '}' ;
      }
      else
	return name ;
      exprP args = expression::create(s,closing) ;
    
      if((closing == ')' )&&(parse::get_token(s,"{"))) {
	name->op_priv = OP_FUNC_BRACE;

	name->expr_list_priv.push_back(args) ;
	closing = '}' ;
	exprP brace_args = expression::create(s, closing) ;
	name->expr_list_priv.push_back(brace_args) ;
      }
      else {
	if(args->op == OP_COMMA)
	  name->expr_list_priv = args->expr_list ;
	else
	  name->expr_list_priv.push_back(args) ;
      }
      return name ;
    }
    exprP p = new expression ;
    return p ;
  }

  exprP expression::get_term(istream &s)
  {
    if(parse::get_token(s,"(")) {
      exprP temp = (expression::create(s,')')) ;
      return temp ;
    }
    if(s.peek() == '-' || s.peek() == '+') {
      char ch = s.get() ;
      while(isspace(s.peek()))
	s.get() ;
      char ch2 = s.peek() ;
      if(ch2 == '.') {
	s.putback('0') ;
	ch2 = '0' ;
      }
      s.putback(ch) ;
      if(isdigit(ch2)) {
	exprP ival = new expression ;
	ival->int_val_priv = parse::get_int(s) ;
	ival->op_priv = OP_INT ;
	if(s.peek() == '.' || s.peek() == 'e' || s.peek() =='E') {
	  double val1 = ival->int_val_priv ;
	  double val2=0 ;
	  if(s.peek() == '.') {
	    s.get() ;
	    double digit = 0.1 ;
	    while(isdigit(s.peek())) {
	      char d = s.get() ;
	      val2 += double(d-'0')*digit ;
	      digit*=0.1 ;
	    }
	  }
	  double exp = 1;
	  if(s.peek() == 'e' || s.peek() == 'E') {
	    s.get() ;
	    int e = parse::get_int(s) ;
#ifdef NO_CMATH
	    exp = ::pow(10.0,double(e)) ;
#else
	    exp = std::pow(10.0,double(e)) ;
#endif
	  }
	  double val = (val1+val2)*exp ;
	  ival->real_val_priv = val ;
	  ival->op_priv = OP_DOUBLE ;
	}
	return ival ;
      }
    }
  
    OpType ot = get_unary_oper(s) ;
    if(ot != OP_ERROR) {
      exprP p = new expression ;
      p->op_priv = ot ;
      p->expr_list_priv.push_back(expression::get_term(s)) ;
      return p ;

    }
    if(parse::is_name(s)) {
      exprP temp  = (expression::get_name(s)) ;
      return temp ;
    }
    if (parse::is_int(s) || s.peek() == '.') {
      if(s.peek() == '.') 
	s.putback('0') ;
      exprP ival = new expression ;
      ival->int_val_priv = parse::get_int(s) ;
      ival->op_priv = OP_INT ;
      if(s.peek() == '.' || s.peek() == 'e' || s.peek() == 'E') {
	double val1 = ival->int_val_priv ;
	double val2=0 ;
	if(s.peek() == '.') {
	  s.get() ;
	  double digit = 0.1 ;
	  while(isdigit(s.peek())) {
	    char d = s.get() ;
	    val2 += double(d-'0')*digit ;
	    digit*=0.1 ;
	  }
	}
	double exp = 1;
	if(s.peek() == 'e' || s.peek() == 'E') {
	  s.get() ;
	  int e = parse::get_int(s) ;
#ifdef NO_CMATH
	  exp = ::pow(10.0,double(e)) ;
#else
	  exp = std::pow(10.0,double(e)) ;
#endif
	}
	double val = (val1+val2)*exp ;
	ival->real_val_priv = val ;
	ival->op_priv = OP_DOUBLE ;
      }
      return ival ;
    }
    if (parse::is_string(s)) {
      exprP sval = new expression ;
      sval->name_priv = parse::get_string(s) ;
      sval->op_priv = OP_STRING ;
      return sval ;
    }
    
    exprP error = new expression ;
    return error ;
    
  }

  exprList collect_associative_op(const exprP &e, const OpType op) {
  
    if(op == e->op) {
      exprList l = e->expr_list, o ;
      while(l.begin() != l.end()) {
	exprP s = l.front() ;
	l.pop_front() ;
	if(op == s->op) {
	  exprList::const_reverse_iterator ii ;
	  for(ii=s->expr_list.rbegin();ii!=s->expr_list.rend();++ii)
	    l.push_front(*ii) ;
	} else {
	  o.push_back(s) ;
	}
      }
      return o ;
	
    } else {
      exprList l ;
      l.push_front(e) ;
      return l ;
    }
  }


  double expression::evaluate(const std::map<std::string, double> &varmap) const {
    double tmp ;
    exprList::const_iterator li ;
    std::map<std::string,double>::const_iterator mi ;
    switch(op) {
    case OP_INT:
      return double(int_val) ;
    case OP_DOUBLE:
      return (real_val) ;
    case OP_PLUS: 
      tmp = 0 ;
      for(li=expr_list.begin();li!=expr_list.end();++li) {
	tmp += (*li)->evaluate(varmap) ;
      }
      return tmp ;
    case OP_MINUS:
      tmp = 0 ;
      li = expr_list.begin() ;
      if(li!=expr_list.end())
	tmp = (*li)->evaluate(varmap) ;
      for(++li;li!=expr_list.end();++li) {
	tmp -= (*li)->evaluate(varmap) ;
      }
      return tmp ;
    case OP_TIMES:
      tmp = 1 ;
      for(li=expr_list.begin();li!=expr_list.end();++li) {
	tmp *= (*li)->evaluate(varmap) ;
      }
      return tmp ;
    case OP_DIVIDE:
      li = expr_list.begin() ;
      tmp = 1 ;
      if(li!=expr_list.end()) 
	tmp =  (*li)->evaluate(varmap) ;
      for(++li;li!=expr_list.end();++li) {
	tmp /= (*li)->evaluate(varmap) ;
      }
      return tmp ;
    case OP_FUNC:
      li = expr_list.begin() ;
      tmp = 0 ;
      if(li!=expr_list.end()) 
	tmp = (*li)->evaluate(varmap) ;
      ++li ;
      if(name == "pow") {
	double tmp2 = 1 ;
	if(li!=expr_list.end()) 
	  tmp2 = (*li)->evaluate(varmap) ;
#ifdef NO_CMATH
        return ::pow(tmp,tmp2) ;
#else
	return std::pow(tmp,tmp2) ;
#endif
      }
#ifdef NO_CMATH
      if(name == "sin")
        return ::sin(tmp) ;
      if(name == "cos")
        return ::cos(tmp) ;
      if(name == "tan")
	return ::tan(tmp) ;
      if(name == "asin") 
	return ::asin(tmp) ;
      if(name == "acos")
	return ::acos(tmp) ;
      if(name == "atan") 
	return ::atan(tmp) ;
      if(name == "sinh") 
	return ::sinh(tmp) ;
      if(name == "cosh")
	return ::cosh(tmp) ;
      if(name == "tanh")
	return ::tanh(tmp) ;
      if(name == "exp")
	return ::exp(tmp) ;
      if(name == "sqrt")
	return ::sqrt(tmp) ;
      if(name == "ln") 
	return ::log(tmp) ;
      if(name == "log") 
	return ::log(tmp) ;
      if(name == "log10")
	return ::log10(tmp) ;
#else
      if(name == "sin")
        return std::sin(tmp) ;
      if(name == "cos")
        return std::cos(tmp) ;
      if(name == "tan")
	return std::tan(tmp) ;
      if(name == "asin") 
	return std::asin(tmp) ;
      if(name == "acos")
	return std::acos(tmp) ;
      if(name == "atan") 
	return std::atan(tmp) ;
      if(name == "sinh") 
	return std::sinh(tmp) ;
      if(name == "cosh")
	return std::cosh(tmp) ;
      if(name == "tanh")
	return std::tanh(tmp) ;
      if(name == "exp")
	return std::exp(tmp) ;
      if(name == "sqrt")
	return std::sqrt(tmp) ;
      if(name == "ln") 
	return std::log(tmp) ;
      if(name == "log") 
	return std::log(tmp) ;
      if(name == "log10")
	return std::log10(tmp) ;
#endif
      {
        string msg = "in expression evaluation, function " + name
          + " has no definition";
        throw exprError("Undefined",msg,ERR_UNDEF) ;
      }
    case OP_NAME:
      mi = varmap.find(name) ;
      if(mi != varmap.end()) 
	return mi->second ;
      if(name == "pi")
	return M_PI ;

      {
        string msg = "in expression evaluation, variable " + name
          + " has no definition";
        throw exprError("Undefined",msg,ERR_UNDEF) ;
      }
	
    default:
      throw exprError("Undefined","operation not defined in evaluate()",ERR_UNDEF) ;
      return 0 ;
    }
    return 0 ;
  }

  
  exprP const_group(exprP p) {
    exprList lgroup,lmgroup ;
    exprList::const_iterator li,lis ;
    
    if(p==0)
      return 0 ;

    switch(p->op) {
    case OP_PLUS:
      {
	lgroup = collect_associative_op(p, OP_PLUS) ;
	if(lgroup.size() == 1)
	  return lgroup.front() ;
	int sum = 0 ;
	for(li=lgroup.begin();li!=lgroup.end();++li) 
	  if((*li)->op==OP_INT)
	    sum += (*li)->int_val ;
	if(sum != 0)
	  lmgroup.push_back(e_int(sum)) ;
	for(li=lgroup.begin();li!=lgroup.end();++li) 
	  if((*li)->op!=OP_INT)
	    lmgroup.push_back(const_group(*li)) ;
	
	if(lmgroup.empty())
	  return exprP(e_int(0)) ;
	if(lmgroup.size() > 1)
	  return exprP(new expression(OP_PLUS,"",lmgroup,0)) ;
	return lmgroup.front() ;
      }
    case OP_TIMES:
      {
	lgroup = collect_associative_op(p, OP_TIMES) ;
	if(lgroup.size() == 1)
	  return lgroup.front() ;

	int prod = 1 ;
	for(li=lgroup.begin();li!=lgroup.end();++li) 
	  if((*li)->op==OP_INT)
	    prod *= (*li)->int_val ;
	if(prod != 1)
	  lmgroup.push_back(e_int(prod)) ;
	if(prod == 0)
	  return e_int(0) ;
	for(li=lgroup.begin();li!=lgroup.end();++li) 
	  if((*li)->op!=OP_INT)
	    lmgroup.push_back(const_group(*li)) ;
	
	if(lmgroup.empty()) 
	  return exprP(e_int(1)) ;
	if(lmgroup.size() > 1)
	  return exprP(new expression(OP_TIMES,"",lmgroup,0)) ;
	return lmgroup.front() ;
      }
    case OP_DIVIDE:
      {
	lgroup = p->expr_list ;
	if(lgroup.size() == 1)
	  return lgroup.front() ;
	int prod = 1 ;
	li = lgroup.begin() ;
	li++ ;
	for(;li!=lgroup.end();++li) 
	  if((*li)->op==OP_INT)
	    prod *= (*li)->int_val ;
	if(prod != 1)
	  lmgroup.push_back(e_int(prod)) ;

	li = lgroup.begin() ;
	li++ ;
	for(;li!=lgroup.end();++li) 
	  if((*li)->op!=OP_INT)
	    lmgroup.push_back(const_group(*li)) ;
	
	if(lmgroup.empty()) 
	  return exprP(const_group(lgroup.front())) ;
	exprList frac ;
	frac.push_back(const_group(lgroup.front())) ;
	if(lmgroup.size() > 1)
	  frac.push_back(new expression(OP_TIMES,"",lmgroup,1)) ;
	else
	  frac.push_back(lmgroup.front()) ;
	return exprP(new expression(OP_DIVIDE,"",frac,1)) ;
      }	
    default:
      break ;
    }      

    if(p->expr_list.empty())
      return p ;
    else {
      for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
	lgroup.push_back(const_group((*li))) ;
      return(new expression(p->op,p->name,lgroup,p->int_val,p->real_val)) ;
    }
  }

  exprP remove_minus(exprP p) {
    switch(p->op) {
    case OP_MINUS:
      { 
	exprList l ;
	exprList::const_iterator li ;
	for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
	  if(li == p->expr_list.begin())
	     l.push_back(remove_minus(*li)) ;
	  else {
	    if( (*li)->op == OP_TIMES) {
	      exprList lp ;
	      lp.push_back(e_int(-1)) ;
	      exprList::const_iterator lpi ;
	      for(lpi=(*li)->expr_list.begin();
		  lpi!=(*li)->expr_list.end();++lpi) {
		lp.push_back(remove_minus(*lpi)) ;
	      }
	      l.push_back(new expression(OP_TIMES,"",lp,1)) ;
	    } else {
	      l.push_back(-1*(remove_minus(*li))) ;
	    }
	  }
	return new expression(OP_PLUS,"",l,0) ;
      }
    case OP_UNARY_PLUS:
      return p->expr_list.front() ;
    case OP_UNARY_MINUS:
	p = p->expr_list.front() ;
	if(p->op == OP_TIMES) {
	  exprList l ;
	  exprList::const_iterator li ;
	  l.push_back(e_int(-1)) ;
	  for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
	    l.push_back(remove_minus(*li)) ;
	  return new expression(p->op,p->name,l,p->int_val,p->real_val) ;
	} else {
	  return -1*remove_minus(p) ;
	}
    default:
      if(p->expr_list.empty())
        return p ;
      else { 
	exprList l ;
	exprList::const_iterator li ;
	for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
	  l.push_back(remove_minus(*li)) ;
	return new expression(p->op,p->name,l,p->int_val,p->real_val) ;
      }
    }      

  }

  exprP remove_divide(exprP p) {
    exprList l ;
    exprList::const_iterator li ;
    for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
      l.push_back(remove_divide(*li)) ;
  
    if(p->op == OP_DIVIDE) {
      exprList l2 ;
      for(li=l.begin();li!=l.end();++li)
	if(li == l.begin()) 
	  l2.push_back(*li) ;
	else if((*li)->op == OP_TIMES) {
	  exprList::const_iterator li2 ;
	  for(li2=(*li)->expr_list.begin();li2!=(*li)->expr_list.end();++li2) {
	    l2.push_back(pow(*li2,-1)) ;
	  }
	} else {
	  l2.push_back(pow(*li,-1)) ;
	}
      return exprP(new expression(OP_TIMES,"",l2,1)) ;
    }
    return exprP(new expression(p->op,p->name,l,p->int_val,p->real_val)) ;
  }

  exprP add_divide(exprP p) {
    exprList l ;
    exprList::const_iterator li ;
    for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
      l.push_back(add_divide(*li)) ;
  
    if(p->op == OP_TIMES) {
      exprList num,denom ;
      
      for(li=l.begin();li!=l.end();++li) {
	exprP e = *li ;
	if(e->op == OP_FUNC && e->name == "pow" && 
	   e->expr_list.back()->op==OP_INT &&
	   e->expr_list.back()->int_val < 0) {
	  int val = e->expr_list.back()->int_val ;
	  exprP m = e->expr_list.front() ;
	  val = -val ;
	  if(val == 1)
	    denom.push_back(m) ;
	  else {
	    denom.push_back(pow(m,val)) ;
	  }
	} else
	  num.push_back(e) ;
      }
      if(denom.empty()) {
	if(num.empty()) 
	  return e_int(1) ;
	else
	  return exprP(new expression(OP_TIMES,"",num,1)) ;
      }
	
      exprList divList ;
      if(num.empty())
	divList.push_back(e_int(1)) ;
      else if(num.size() == 1)
	divList.push_back(num.front()) ;
      else
	divList.push_back(new expression(OP_TIMES,"",num,1)) ;

      if(denom.size() == 1)
	divList.push_back(denom.front()) ;
      else
	divList.push_back(new expression(OP_TIMES,"",denom,1)) ;
      return exprP(new expression(OP_DIVIDE,"",divList,1)) ;
    }
    return exprP(new expression(p->op,p->name,l,p->int_val,p->real_val)) ;
  }

  exprP simplify_expr(exprP e) {
    switch(e->op) {
    case OP_INT:
    case OP_DOUBLE:
    case OP_STRING:
    case OP_ERROR:
    case OP_NAME:
      return e ;
      break ;
    case OP_PLUS:
      {
	// Condense ints and doubles
	exprList::const_iterator li ;
	exprList lstart ;
	int ival = 0;
	double dval = 0 ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprP p = simplify_expr(*li) ;
	  if(p->op == OP_INT)
	    ival += p->int_val ;
	  else if (p->op == OP_DOUBLE)
	    dval += p->real_val ;
	  else
	    lstart.push_back(p) ;
	}
	if(dval == 0 && ival != 0) {
	  lstart.push_front(e_int(ival)) ;
	} else if(dval != 0) {
	  lstart.push_front(new expression(OP_DOUBLE,"",exprList(),ival,
					   double(ival)+dval)) ;
	}

        if(lstart.empty())
          return e_int(0) ;
        if(lstart.size() == 1)
          return lstart.front() ;

	exprList l ;
	map<exprP,int> exp_map ;
	map<exprP,int>::iterator mi ;
	for(li=lstart.begin();li!=lstart.end();++li) {
	  exprP p = (*li) ;
	  if(p->op == OP_TIMES && p->expr_list.back()->op==OP_INT) {
	    exprList l = p->expr_list ;
	    int val = l.back()->int_val ;
	    l.pop_back() ;
	    if(l.size()==1)
	      p = l.front() ;
	    else
	      p = exprP(new expression(OP_TIMES,"",l,1)) ;
	    if((mi=exp_map.find(p)) == exp_map.end())
	      exp_map[p] = val ;
	    else {
	      mi->second += val ;
	    }
	    
	  } else {
	    if((mi=exp_map.find(p)) == exp_map.end())
	      exp_map[p] = 1 ;
	    else {
	      mi->second += 1 ;
	    }
	  }
	}
	exprList l2 ;
	for(mi=exp_map.begin();mi!=exp_map.end();++mi) {
	  if(mi->second == 1) {
	    l2.push_back(mi->first) ;
	  } else {
	    exprP c = new expression(OP_INT,"",exprList(),mi->second) ;
	    exprList m ;
	    m.push_back(c) ;
	    m.push_back(mi->first) ;
	    l2.push_back(exprP(new expression(OP_TIMES,"",m,1))) ;
	  }
	}
	return exprP(new expression(OP_PLUS,"",l2,0)) ;
      }

    case OP_TIMES:
      {
	// Condense ints and doubles
	exprList::const_iterator li ;
	exprList lstart ;
	int ival = 1;
	double dval = 1 ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprP p = simplify_expr(*li) ;
	  if(p->op == OP_INT)
	    ival *= p->int_val ;
	  else if (p->op == OP_DOUBLE)
	    dval *= p->real_val ;
	  else
	    lstart.push_back(p) ;
	}

	if(dval == 1 && ival != 1) {
	  lstart.push_front(e_int(ival)) ;
	} else if(dval != 1) {
	  if(dval == 0)
	    lstart.push_front(e_int(0)) ;
	  else 
	    lstart.push_front(new expression(OP_DOUBLE,"",exprList(),ival,
					     double(ival)*dval)) ;
        }

        if(lstart.empty())
          return e_int(1) ;
        if(lstart.size() == 1)
          return lstart.front() ;
	
	// Combine pow() functions
	exprList l ;
	map<exprP,int> exp_map ;
	map<exprP,int>::iterator mi ;
	for(li=lstart.begin();li!=lstart.end();++li) {
	  exprP p = (*li) ;
	  if(p->op == OP_FUNC && p->name == "pow" && p->expr_list.back()->op==OP_INT) {
	    int val = p->expr_list.back()->int_val ;
	    p = p->expr_list.front() ;
	    if((mi=exp_map.find(p)) == exp_map.end()) 
	      exp_map[p] = val ; 
	    else {
	      mi->second += val ; 
	    }
	  } else {
	    if((mi=exp_map.find(p)) == exp_map.end()) 
	      exp_map[p] = 1 ;
	    else {
	      mi->second += 1 ;
	    }
	  }
	}
	exprList l2 ;
	for(mi=exp_map.begin();mi!=exp_map.end();++mi) {
	  if(mi->second == 1) {
	    l2.push_back(mi->first) ;
	  } else {
	    if(mi->second == 0)
	      l2.push_back(e_int(1)) ;
	    else {
	      l2.push_back(pow(mi->first,mi->second)) ;
	    }
	  }
	}
	return exprP(new expression(OP_TIMES,"",l2,0)) ;
      }
    case OP_FUNC:
      if(e->name == "pow" ) { // Check for power of power
	if(e->expr_list.front()->op == OP_FUNC &&
	   e->expr_list.front()->name == "pow") {
	  exprP p1 = const_group(simplify_expr(e->expr_list.back())) ;
	  exprP p2 = const_group(simplify_expr(e->expr_list.front()->expr_list.back())) ;
	  
	  exprP arg1 = simplify_expr(e->expr_list.front()->expr_list.front()) ;
	  if(p1->op == OP_INT && p2->op == OP_INT) 
	    return simplify_expr(pow(arg1,e_int(p1->int_val*p2->int_val))) ;
	  return pow(pow(arg1,p2),p1) ;
	    
	} else {
	  exprP arg1 = simplify_expr(e->expr_list.front()) ;
	  exprP arg2 = simplify_expr(e->expr_list.back()) ;
	  if(arg1->op == OP_INT) {
	    if(arg1->int_val == 0) 
	      return e_int(0) ;
	    if(arg1->int_val == 1)
	      return e_int(1) ;
	    if(arg2->op == OP_INT) {
	      int exp = arg2->int_val ;
	      int p = 1 ;
              if(exp>=0) {
                for(int i=0;i<exp;++i)
                  p *= arg1->int_val ;
                return e_int(p) ;
              }
              exp*=-1 ;
              for(int i=0;i<exp;++i)
                p *= arg1->int_val ;
              return pow(e_int(p),e_int(-1)) ;
	    }
	    if(arg2->op == OP_DOUBLE) {
	      double m = arg1->int_val ;
	      double e = arg2->real_val ;
	      double v = ::pow(m,e) ;
	      return new expression(OP_DOUBLE,"",exprList(),0,v) ;
	    }
	  }
	  if(arg1->op == OP_DOUBLE && arg2->op == OP_DOUBLE) {
	    double m = arg1->real_val ;
	    double e = arg2->real_val ;
	    double v = ::pow(m,e) ;
	    return new expression(OP_DOUBLE,"",exprList(),0,v) ;
	  }
	  if(arg1->op == OP_DOUBLE && arg2->op == OP_INT) {
	    double m = arg1->real_val ;
	    double e = arg2->int_val ;
	    double v = ::pow(m,e) ;
	    return new expression(OP_DOUBLE,"",exprList(),0,v) ;
	  }
	  // pow(x,1) == x
	  if(arg2->op == OP_INT && arg2->int_val == 1)
	    return arg1 ;
	  return pow(arg1,const_group(arg2)) ;
	}
      } else {
	exprList l ;
	exprList::const_iterator li ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
	  l.push_back(simplify_expr(*li)) ;
	return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
      }
    default:
      {
	exprList l ;
	exprList::const_iterator li ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
	  l.push_back(simplify_expr(*li)) ;
	return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
      }
      break ;
    }
    return e ;
  }


  exprP derivative(exprP e,const std::string &var) {
    exprList p ;
    switch(e->op) {
    case OP_INT:
    case OP_DOUBLE:
      return e_int(0) ;
    case OP_NAME:
      if(e->name == var) 
	return e_int(1) ;
      else
	return e_int(0) ;
    case OP_COMMA:
      { exprList nlist ;
	exprList::const_iterator li ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprP tmp = derivative(*li,var) ;
	  nlist.push_back(tmp) ;
	}
	return exprP(new expression(OP_COMMA,e->name,nlist,e->int_val)) ;
      }
    case OP_PLUS:
      { exprList nlist ;
	exprList::const_iterator li ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprP tmp = derivative(*li,var) ;
	  if(!(tmp->op == OP_INT && tmp->int_val == 0))
	    nlist.push_back(tmp) ;
	}
	if(nlist.size() == 0)
	  return e_int(0) ;
	  
	return exprP(new expression(OP_PLUS,e->name,nlist,e->int_val)) ;
      }
    case OP_MINUS:
      { exprList nlist ;
	exprList::const_iterator li ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprP tmp = derivative(*li,var) ;
	  if(!(tmp->op == OP_INT && tmp->int_val == 0))
	    nlist.push_back(tmp) ;
	}
	if(nlist.size() == 0)
	  return e_int(0) ;
	  
	return exprP(new expression(OP_MINUS,e->name,nlist,e->int_val)) ;
      }
    case OP_TIMES:
      { exprList nlist ; // Terms from the product rule
	exprList::const_iterator li ;
	exprList::const_iterator lip ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprList plist ; // list of products ;
	  bool iszero = false ;
	  for(lip=e->expr_list.begin();lip!=e->expr_list.end();++lip) {
	    if(lip == li) {
	      exprP tmp = (*lip)->derivative(var) ;
	      if(tmp->op == OP_INT && tmp->int_val == 0)
		iszero = true ;
	      if(!(tmp->op == OP_INT && tmp->int_val == 1))
		plist.push_back(tmp) ;
	    }else
	      plist.push_back(*lip) ;
	  }

	  if(!iszero) {
	    nlist.push_back(exprP(new expression(OP_TIMES,"",plist,0))) ;
	  }
	}
	if(nlist.empty()) 
	  return e_int(0) ;
	else
	  return exprP(new expression(OP_PLUS,"",nlist,0)) ;
      }
    case OP_DIVIDE:
      { exprList terms ; // Terms from the product rule
	exprList::const_iterator li ;
	exprList::const_iterator lip ;
	for(li=e->expr_list.begin();li!=e->expr_list.end();++li) {
	  exprList nlist ; // list of numerators ;
	  exprList dlist ; // list of denominators ;
	  bool iszero = false ;
	  for(lip=e->expr_list.begin();lip!=e->expr_list.end();++lip) {
	    if(lip == li) {
	      if(lip == e->expr_list.begin()) { // First term is derivative
		exprP dfx = (*lip)->derivative(var) ;
		if(dfx->op == OP_INT && dfx->int_val == 0)
		  iszero = true ;
		nlist.push_back(dfx) ; 
	      } else { // all other terms are 1/f(x)
		exprP dfx = (*lip)->derivative(var) ; 
		if(dfx->op == OP_INT && dfx->int_val == 0)
		  iszero = true ;
		nlist.push_back(e_int(-1)) ;
		nlist.push_back(dfx) ;
		dlist.push_back(pow(*lip,2)) ;
	      }
	    }else
	      if(lip == e->expr_list.begin())
		nlist.push_back(*lip) ;
	      else
		dlist.push_back(*lip) ;
	  }

	  if(!iszero) {
	    exprList pl ;
	    pl.push_back(new expression(OP_TIMES,"",nlist,1)) ;
	    pl.push_back(new expression(OP_TIMES,"",dlist,1)) ;
	    terms.push_back(new expression(OP_DIVIDE,"",pl,1)) ;
	  }
	}
	if(terms.empty()) 
	  return e_int(0) ;

	return exprP(new expression(OP_PLUS,"",terms,0)) ;
      }
    case OP_FUNC:
      {
	if(e->name == "pow") {
	  exprP arg1 = const_group(e->expr_list.front()) ;
	  exprP arg2 = const_group(e->expr_list.back()) ;
	  exprP darg1 = Loci::derivative(arg1,var) ;
	  exprP darg2 = Loci::derivative(arg2,var) ;

	  if(darg2->op == OP_INT && darg2->int_val == 0) {
	    if(darg1->op == OP_INT && darg1->int_val == 0) {
	      return darg2 ;
	    }
	    exprList l ;
	    exprList plist ;
	    plist.push_back(e->expr_list.front()) ;
	    l.push_back(arg1) ;
	   
	    if(arg2->op == OP_INT) {
	      if(arg2->int_val == 1) 
		return darg1 ;
	      if(arg2->int_val == 2) {
		exprList l2 ;
		l2.push_back(arg2) ;
		l2.push_back(arg1) ;
		l2.push_back(darg1) ;
		return exprP(new expression(OP_TIMES,"",l2,1)) ;
	      }
	      l.push_back(e_int(arg2->int_val-1)) ;
	    } else {
	      l.push_back(arg2-1) ;
	    }	      
	      
	    exprList l2 ;
	    // pow derivative push
	    l2.push_back(darg1) ; 
	    l2.push_back(arg2) ;
	    l2.push_back(exprP(new expression(OP_FUNC,"pow",l,0)));
	    return exprP(new expression(OP_TIMES,"",l2,1)) ;
	  } else {
	    // exponent is a function of var
	    if(darg1->op == OP_INT && darg1->int_val == 0) {
	      // base is not function of var, then the derivative of
	      // a^f(x) == ln(a) a^f(x) f'(x)
	      exprList l ;
	      l.push_back(ln(arg1)) ;
	      l.push_back(e) ;
	      l.push_back(darg2) ;
	      return exprP(new expression(OP_TIMES,"",l,1)) ;
	    } else {
	      // most general case f(x)^g(x)
	      // = f(x)^(g(x)-1)*f'(x) g(x)+f(x)^g(x) g'(x) ln(f(x))
	      
	      exprList l1 ;
	      // f(x)^(g(x)-1)
	      l1.push_back(pow(arg1,arg2-1)) ;
	      l1.push_back(darg1) ;// f'(x)
	      l1.push_back(arg2) ; // g(x) ;

	      exprList ls ;
	      ls.push_back(exprP(new expression(OP_TIMES,"",l1,1))) ;

	      exprList l2 ;
	      l2.push_back(ln(arg1)) ;
	      l2.push_back(e) ;
	      l2.push_back(darg2) ;
	      ls.push_back(exprP(new expression(OP_TIMES,"",l2,1))) ;
	      return exprP(new expression(OP_PLUS,"",ls,0)) ;
	    }
	    
	  }
	}
	      
	if(e->expr_list.size() != 1) 
	  return 0 ;
	exprP darg = Loci::derivative(e->expr_list.front(),var) ;
	if(darg->op == OP_INT && darg->int_val == 0)
	  return darg ;
	exprP df ;
	if(e->name == "sin")
	  df = cos(e->expr_list.front()) ;
	else if (e->name == "cos") {
	  df = -1*sin(e->expr_list.front()) ;
	} else if(e->name == "tan") {
	  df = darg/pow(cos(e->expr_list.front()),e_int(2)) ;
	  return df ;
	} else if(e->name == "asin") {
	  df = darg/sqrt(e_int(1)-pow(e->expr_list.front(),e_int(2))) ;
	  return df ;
	} else if(e->name == "acos") {
	  df = e_int(-1)*darg/sqrt(e_int(1)-pow(e->expr_list.front(),e_int(2))) ;
	  return df ;
	} else if(e->name == "atan") {
	  df = darg/(e_int(1)+pow(e->expr_list.front(),e_int(2))) ;
	  return df ;
	} else if(e->name == "ln" || e->name == "log") {
	  df = darg/e->expr_list.front() ;
	  return df ;
	}else if(e->name == "log10") {
	  df = darg/(ln(e_int(10))*e->expr_list.front()) ;
	  return df ;
	} else if(e->name == "exp") {
	  df = exp(e->expr_list.front()) ;
	} else if(e->name == "exp") {
	  df = exp(e->expr_list.front()) ;
	} else if(e->name == "sinh") {
	  df = cosh(e->expr_list.front()) ;
	} else if(e->name == "cosh") {
	  df = sinh(e->expr_list.front()) ;
	} else if(e->name == "tanh") {
	  df = e_int(1)-pow(sinh(e->expr_list.front()),e_int(2))/pow(cosh(e->expr_list.front()),e_int(2)) ;
	} else if(e->name == "sqrt") {
	  df = darg/(2*sqrt(e->expr_list.front())) ;
	  return df ;
	} else 
	  return 0 ;
	if(darg->op == OP_INT && darg->int_val == 1)
	  return df ;
	return darg*df ;
      }
    default:
      {
        ostringstream oss ;
        oss << "derivative of " << e << endl ;
        
        throw exprError("Not Supported",oss.str(),ERR_UNDEF) ;
      }
      return exprP(e_int(0)) ;
    }
    
    return 0 ;
  }

  exprP substitute_expr(exprP p, exprP s, exprP e) {
    if(compare_expressions(p,s) == 0)
      return e ;
    exprList l ;
    exprList::const_iterator li ;
    for(li=p->expr_list.begin();li!=p->expr_list.end();++li)
      l.push_back(substitute_expr(*li,s,e)) ;
    return new expression(p->op,p->name,l,p->int_val,p->real_val) ;
  }

  exprP expression::substitute(exprP s, exprP e) const {
    return substitute_expr(new expression(op,name,expr_list,int_val,real_val),
			   s,e) ;
  }

  void getVarNames(exprP e, set<string> &namelist) {
    if(e->op == OP_NAME)
      namelist.insert(e->name) ;
    exprList::const_iterator li ;
    for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
      getVarNames(*li,namelist) ;
  }

  exprP substitutionEngine(exprP target, exprP list) {
    map<string,exprP> sub_map ;
    exprList l ;
    if(list->op == OP_COMMA) 
      l = list->expr_list ;
    else
      l.push_back(list) ;
    exprList::const_iterator li ;

    for(li=l.begin();li!=l.end();++li) {
      if((*li)->op == OP_ASSIGN) {
	if((*li)->expr_list.size() != 2)
          throw exprError("Syntax Error","unable to interpret substitution",ERR_BADFORM) ;
	if((*li)->expr_list.front()->op != OP_NAME)
          throw exprError("Syntax Error","substitution rhs should be name",ERR_BADFORM) ;
        sub_map[(*li)->expr_list.front()->name] = (*li)->expr_list.back() ;
      }
    }

    exprP work = target ;

    const int max_depth = 200 ;
    
    for(int i=0;i<max_depth;++i) {
      set<string> namelist ;
      getVarNames(work,namelist) ;
      set<string>::const_iterator  si ;
      bool sub_found = false ;
      for(si=namelist.begin();si!=namelist.end();++si) {
	string name = *si ;
	map<string,exprP>::const_iterator mi ;
	if((mi = sub_map.find(name)) != sub_map.end()) {
	  exprP namep = new expression(OP_NAME,name,exprList(),0,0.0) ;
	  work = substitute_expr(work,namep,mi->second) ;
	  sub_found = true ;
	}
      }
      if(!sub_found)
	return work ;
    }

    throw exprError("Limit Exceeded","Recursive Depth Exceeded in Substitution",ERR_LIMIT) ;
    return target ;
  }

  
  exprP expression::derivative(std::string var) const {
    exprP p = remove_minus(new expression(op,name,expr_list,int_val,real_val)) ;
    p = const_group(p) ;
    return const_group(Loci::derivative(p,var)) ;
  }

  exprP expression::constant_grouping() const {
    return const_group(exprP(new expression(op,name,expr_list,int_val,real_val))) ;
  }

  

  exprP makeCannon(exprP e) {
    exprList l ;
    exprList::const_iterator li ;
    for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
      l.push_back(makeCannon(*li)) ;
    switch(e->op) {
    case OP_PLUS:
    case OP_TIMES:
      l.sort() ;
      break ;
    default:
      break ;
    }

    switch(e->op) {
    case OP_FUNC:
    case OP_ARRAY:
    case OP_NAME_BRACE:
    case OP_FUNC_BRACE:
    case OP_UNARY_PLUS:
    case OP_UNARY_MINUS:
    case OP_NOT:
    case OP_TILDE:
    case OP_AMPERSAND:
    case OP_DOLLAR:
    case OP_STAR:
      return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
      
    default:
      if(l.size() == 1) {
	return l.front() ;
      }else
	return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
    }

  }

  exprP factor_expression(exprP e) {

    switch(e->op) {
    case OP_INT:
    case OP_DOUBLE:
    case OP_STRING:
    case OP_ERROR:
    case OP_NAME:
      return e ;
    default:
      break ;
    }
      
    exprList l ;
    exprList::const_iterator li ;
    for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
      l.push_back(factor_expression(*li)) ;


    map<exprP,vector<int> > exp_map ;
    map<exprP,vector<int> >::iterator mi ;
    vector<exprList> vlist ;
    switch(e->op) {
    case OP_PLUS:
      {
	for(li=l.begin();li!=l.end();++li) {
	  if((*li)->op == OP_TIMES) {
	    vlist.push_back((*li)->expr_list) ;
	    exprList::const_iterator vli ;
	    for(vli=vlist.back().begin();vli!=vlist.back().end();++vli) {
	      if((mi=exp_map.find(*vli)) != exp_map.end()) {
		if(mi->second.back() != int(vlist.size())-1)
		  mi->second.push_back(int(vlist.size())-1) ;
	      } else
		exp_map[*vli].push_back(int(vlist.size())-1) ;
	    }
	  } else {
	    exprList l ;
	    l.push_back(*li) ;
	    vlist.push_back(l) ;
	    // Add term to map
	    if((mi=exp_map.find(*li)) != exp_map.end()) {
	      if(mi->second.back() != int(vlist.size())-1)
		mi->second.push_back(int(vlist.size())-1) ;
	    } else
	      exp_map[*li].push_back(int(vlist.size())-1) ;
	  }
	}
	int mx = 0 ;
	exprP mp ;
	vector<int> terms ;
	for(mi=exp_map.begin();mi!=exp_map.end();++mi) {
	  if(int(mi->second.size()) > mx) {
	    mx = mi->second.size() ;
	    terms = mi->second ;
	    mp = mi->first ;
	  }
	}
	// Now found most common sub-expression, factor out.
	if(mx > 1) {
	  exprList factor ;
	  exprList nonfactor ;
	  for(size_t i=0;i<vlist.size();++i) {
	    bool found = false ;
	    for(size_t j=0;j<terms.size();++j)
	      if(terms[j]==int(i)) 
		found = true ;
	    if(found) {
	      // This is a term of the factor
	      exprList parts ;
	      found = false ;
	      for(li=vlist[i].begin();li!=vlist[i].end();++li) {
		if(compare_expressions(*li,mp) != 0) // if not equal
		  parts.push_back(*li) ;
		else {
		  if(found)
		    parts.push_back(*li) ; // always push after found
		  found = true ; // only find the factor once
		}
	      }
	      if(parts.empty()) 
		factor.push_back(e_int(1)) ;
	      else if(parts.size() == 1)
		factor.push_back(parts.front()) ;
	      else
		factor.push_back(new expression(OP_TIMES,"",parts,1)) ;
	    } else {
	      if(vlist[i].size() == 1) {
		nonfactor.push_back(vlist[i].front()) ;
	      } else
		nonfactor.push_back(new expression(OP_TIMES,"",vlist[i],1)) ;
	    }
	  }
	  exprP termp = new expression(OP_PLUS,"",factor,0) ;
	  exprP factorp = termp*mp ;
	  if(nonfactor.empty()) {
	    return factorp ;
	  }
	  nonfactor.push_back(factorp) ;
	  exprP total = new expression(OP_PLUS,"",nonfactor,0) ;
	  return total ;
	}
      }
      // Now loop over factors and find if one appears more than once
      // factor out greatest common one.
      break ;
    default:
      break ;
    }
    return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
  }


  exprP symbolic_evaluate(exprP e) {

    switch(e->op) {
    case OP_INT:
    case OP_DOUBLE:
    case OP_STRING:
    case OP_ERROR:
    case OP_NAME:
      return e ;
    default:
      break ;
    }

    exprList l ;
    exprList::const_iterator li ;
    for(li=e->expr_list.begin();li!=e->expr_list.end();++li)
      l.push_back(symbolic_evaluate(*li)) ;

    if(e->op == OP_FUNC && e->name == "del") {
      if(l.size() != 2) {
        ostringstream oss ;
        oss << "del operator malformed, needs two arguments e= " << e << endl ;
        throw exprError("Syntax Error",oss.str(),ERR_BADFORM) ;
      }
      if(l.back()->op != OP_NAME) {
        ostringstream oss ;
        oss << "malformed del operator = " << e << endl ;
        throw exprError("Syntax Error",oss.str(),ERR_BADFORM) ;
      }
      exprP p = Loci::derivative(l.front(),l.back()->name) ;
      return p ;
    }
    
    return exprP(new expression(e->op,e->name,l,e->int_val,e->real_val)) ;
  }

    
  exprP expression::symbolic_eval() const {
    return symbolic_evaluate(exprP(new expression(op,name,expr_list,int_val,real_val))) ;
  }

  exprP expression::simplify() const {
    exprP p = remove_minus(exprP(new expression(op,name,expr_list,int_val,real_val))) ;
    p = remove_divide(p) ;
    p = const_group(p) ;

    exprP ps = p ;
    do {
      p = ps ;
      ps = makeCannon(p) ; // Make cannonical form
      ps = simplify_expr(ps) ;
      ps = const_group(ps) ; 
      ps = const_group(ps) ; // Two const groups should contract any constant
      //factors added in simplify expression
      ps = factor_expression(ps) ;// factor out common terms
      // Repeat until expression remains the same
    } while(compare_expressions(ps,p) != 0) ;

    p = add_divide(p) ;
    return p ;
  }
}

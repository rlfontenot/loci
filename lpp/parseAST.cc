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
#include "parseAST.h"
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

void AST_Token::DiagPrint(ostream &s, int &line) const {
  if(line != lineno) {
    cout << endl ;
    
    if(line < 0 || line+1 != lineno) {
      cout << "#line " << lineno << endl ;
    }
    line = lineno ;
  }
  s <<text << ' ' ;
}

void AST_SimpleStatement::DiagPrint(ostream  &s, int &line) const {
  exp->DiagPrint(s,line) ;
  Terminal->DiagPrint(s,line) ;
}

void AST_Block::DiagPrint(ostream &s, int &lineno) const {
  for(ASTList::const_iterator ii=elements.begin();ii!=elements.end();++ii)
    (*ii)->DiagPrint(s,lineno) ;
}

void AST_typeDecl::DiagPrint(ostream &s, int &lineno) const {
  for(ASTList::const_iterator ii=type_decl.begin();ii!=type_decl.end();++ii)
    (*ii)->DiagPrint(s,lineno) ;
}

void AST_declaration::DiagPrint(ostream &s, int &lineno) const {
  type_decl->DiagPrint(s,lineno) ;
  decls->DiagPrint(s,lineno) ;
}

void AST_exprOper::DiagPrint(ostream &s, int &lineno) const {
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

void AST_term::DiagPrint(ostream&s, int &lineno) const {
  term->DiagPrint(s,lineno) ;
}

void AST_ifStatement::DiagPrint(ostream&s, int &lineno) const {
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


void AST_loopStatement::DiagPrint(ostream&s, int &lineno) const {
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

void AST_switchStatement::DiagPrint(ostream&s, int &lineno) const {
  statement->DiagPrint(s,lineno) ;
  s << "(" ;
  conditional->DiagPrint(s,lineno) ;
  s << ") {" ;
  for(ASTList::const_iterator ii=body.begin();ii!=body.end();++ii)
    (*ii)->DiagPrint(s,lineno) ;
  s << "}" ;
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

#ifndef PARSEAST_H
#define PARSEAST_H

#include <Tools/cptr.h>
using Loci::CPTR ;
using Loci::CPTR_type ;
#include <list>
using std::list ;
#include <vector>
using std::vector ;
#include <string>
using std::string ;

#include <iostream>
using std::istream ;
using std::ifstream ;
using std::ofstream ;
using std::ostream ;
using std::ios ;
using std::endl ;
using std::cerr ;
using std::cout ;

// Setup for facilties for parsing and creating abstract syntax trees (AST)


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

extern std::string OPtoString(AST_type::elementType val) ;

  class AST_Token : public AST_type {
public:
  string text ;
  int lineno ;
  void DiagPrint(ostream &s, int &line) const ;
} ;

class AST_SimpleStatement: public AST_type {
public:
  AST_SimpleStatement(ASTP e, ASTP t) : exp(e),Terminal(t) {}
  ASTP exp ;
  ASTP Terminal ;
  void DiagPrint(ostream  &s, int &line) const ;
} ;

class AST_Block : public AST_type {
public:
  ASTList elements ;
  void DiagPrint(ostream &s, int &lineno) const ;
} ;

class AST_typeDecl : public AST_type {
public:
  ASTList type_decl ;
  void DiagPrint(ostream &s, int &lineno) const ;
} ;
  
class AST_declaration : public AST_type {
public:
  ASTP type_decl ;
  ASTP decls ;
  void DiagPrint(ostream &s, int &lineno) const ;
} ;

class AST_exprOper : public AST_type {
public:
  ASTList terms ;
  void DiagPrint(ostream &s, int &lineno) const ;
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
  void DiagPrint(ostream&s, int &lineno) const ;
} ;

class AST_ifStatement : public AST_type {
public:
  ASTP iftok ;
  ASTP conditional ;
  ASTP ifblock ;
  ASTP elseblock ;
 AST_ifStatement(ASTP tok, ASTP C, ASTP IF, ASTP ELSE):
  iftok(tok),conditional(C),ifblock(IF),elseblock(ELSE) {} 
  void DiagPrint(ostream&s, int &lineno) const ;
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

  void DiagPrint(ostream&s, int &lineno) const ;
} ;
  
class AST_switchStatement : public AST_type {
public:
  ASTP statement ;
  ASTP conditional ;
  ASTList body ;
  AST_switchStatement() {}

  void DiagPrint(ostream&s, int &lineno) const ;
} ;

extern AST_type::ASTP parseExpression(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseDeclaration(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseStatement(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseIfStatement(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount)  ;
extern AST_type::ASTP parseBlock(std::istream &is, int &linecount) ;
extern AST_type::ASTP parseTerm(std::istream &is, int &linecount) ;
extern CPTR<AST_Token> getToken(std::istream &is, int &linecount) ;
extern void pushToken(CPTR<AST_Token> &pt) ;
extern bool isTerm(AST_type::elementType e) ;

extern istream &killsp(istream &s, int &lines) ;
extern bool is_name(istream &s) ;
extern string get_name(istream &s) ;




#endif

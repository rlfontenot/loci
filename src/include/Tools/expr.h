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
#ifndef EXPR_H
#define EXPR_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <iostream>
#include <string>

#include <Tools/cptr.h>
#include <Tools/except.h>

#include <list>
#include <string>
#include <map>
#include <set>

namespace Loci {

  enum exprErrorType {
    ERR_SYNTAX=1,
    ERR_UNDEF=2,
    ERR_BADFORM=3,
    ERR_LIMIT=4
  } ;
    
  struct exprError : public BasicException {
  public:
    std::string err,msg ;
    exprError() {}
    exprError(std::string e, std::string m, exprErrorType c) :err(e),msg(m) {code = c;}
    virtual ~exprError() {}
    virtual std::ostream &Print(std::ostream &s) const {
      s << err << ": " << msg << std::endl ;
      return s ;
    }
  } ;
  
  enum OpType {
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
    OP_ASSIGN=0xd00,
    OP_COMMA=0xe00, OP_COLON=0xf00, 
    // terminal for empty statement
    OP_NIL=0x1000,
    // terminals for variable name, function, array or name{args}
    OP_NAME, OP_FUNC, OP_ARRAY, OP_NAME_BRACE, OP_FUNC_BRACE,
    // terminal for string, integer, or unspecified error condition
    OP_STRING, OP_INT, OP_DOUBLE, OP_ERROR,
    // Unary operations
    OP_UNARY_PLUS, OP_UNARY_MINUS, OP_NOT, OP_TILDE,
    OP_AMPERSAND, OP_DOLLAR, OP_STAR
  } ;

  class expression : public CPTR_type {
  public:
    typedef CPTR<expression> exprP ;
    typedef std::list<exprP> exprList ;
  private:
    OpType              op_priv ;
    exprList            expr_list_priv ;
    std::string         name_priv ;
    int                 int_val_priv ;
    double              real_val_priv ;
      
  expression() : op(op_priv), expr_list(expr_list_priv), name(name_priv),int_val(int_val_priv),real_val(real_val_priv) { op_priv = OP_ERROR ; }

    static exprP  create( std::istream &s, char closing ) ;
    static exprP  get_term( std::istream &s ) ;
    static exprP  get_name( std::istream &s ) ;
    static OpType get_oper( std::istream &s ) ;
    static OpType get_unary_oper( std::istream &s ) ;
    static exprP  expand_oper( std::istream &s, exprP &p ) ;

    void PrintOperation(std::ostream &s, std::string oper,
                        char poChar = '(', char pcChar = ')' ) const ;
    exprP constant_grouping() const ;
  public:
    const OpType              &op ;
    const exprList            &expr_list ;
    const std::string         &name ;
    const int                 &int_val ;
    const double              &real_val ;
  expression(OpType opin, const std::string nm, 
	     const exprList &elist, int ival = 0, double fval = 0) :
    op(op_priv), expr_list(expr_list_priv), name(name_priv),
      int_val(int_val_priv), real_val(real_val_priv)
    { op_priv = opin; expr_list_priv = elist ; name_priv = nm; 
      int_val_priv = ival ; real_val_priv = fval ;}
    void Print(std::ostream &s) const ;
    double evaluate(const std::map<std::string,double> &varmap) const ;
    exprP simplify() const ; // Simplify the expression
    exprP substitute(exprP s, exprP r) const ; // Substitute all s for r
    exprP derivative(std::string var) const ; // symbolic differentation
    exprP symbolic_eval() const ; // Return evaluation of symbolic opers
    static exprP create(std::istream &s) ;
    static exprP create(const std::string &s) ;

  } ;

  typedef expression::exprP     exprP ;
  typedef expression::exprList  exprList ;

  inline std::ostream &operator<<(std::ostream &s, const exprP &exp) {
    exp->Print(s) ;
    return s;
  }

  inline std::istream &operator>>(std::istream &s, exprP &exp) {
    exp = expression::create(s) ;
    return s ;
  }


  exprList collect_associative_op(const exprP &e, const OpType op) ;

  int compare_expressions(exprP e1, exprP e2) ;

  inline exprP e_int(int val) {
    return exprP(new expression(OP_INT,"",exprList(),val)) ;
  }


  inline bool operator<(exprP e1, exprP e2) {
    return compare_expressions(e1,e2) < 0 ;
  }
  inline bool operator<=(exprP e1, exprP e2) {
    return compare_expressions(e1,e2) <= 0 ; 
  }
  inline bool operator==(exprP e1, exprP e2) {
    return compare_expressions(e1,e2) == 0 ; 
  }
  inline bool operator!=(exprP e1, exprP e2) {
    return compare_expressions(e1,e2) != 0 ; 
  }


  void getVarNames(exprP e, std::set<std::string> &namelist) ;
  exprP substitutionEngine(exprP target, exprP list) ;

  namespace expr {
    inline exprP operator+(exprP p1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(p2) ;
      return new expression(OP_PLUS,"",lgroup,0) ;
    }

    inline exprP operator+(exprP p1, int i2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(e_int(i2)) ;
      return new expression(OP_PLUS,"",lgroup,0) ;
    }

    inline exprP operator+(int i1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(e_int(i1)) ;
      lgroup.push_back(p2) ;
      return new expression(OP_PLUS,"",lgroup,0) ;
    }

    inline exprP operator-(exprP p1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(p2) ;
      return new expression(OP_MINUS,"",lgroup,0) ;
    }

    inline exprP operator-(exprP p1, int i2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(e_int(i2)) ;
      return new expression(OP_MINUS,"",lgroup,0) ;
    }

    inline exprP operator-(int i1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(e_int(i1)) ;
      lgroup.push_back(p2) ;
      return new expression(OP_MINUS,"",lgroup,0) ;
    }

    inline exprP operator*(exprP p1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(p2) ;
      return new expression(OP_TIMES,"",lgroup,1) ;
    }

    inline exprP operator*(exprP p1, int i2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(e_int(i2)) ;
      return new expression(OP_TIMES,"",lgroup,1) ;
    }

    inline exprP operator*(int i1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(e_int(i1)) ;
      lgroup.push_back(p2) ;
      return new expression(OP_TIMES,"",lgroup,1) ;
    }

    inline exprP operator/(exprP p1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(p2) ;
      return new expression(OP_DIVIDE,"",lgroup,1) ;
    }

    inline exprP operator/(exprP p1, int i2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(e_int(i2)) ;
      return new expression(OP_DIVIDE,"",lgroup,1) ;
    }

    inline exprP operator/(int i1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(e_int(i1)) ;
      lgroup.push_back(p2) ;
      return new expression(OP_DIVIDE,"",lgroup,1) ;
    }

    inline exprP pow(exprP p1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(p2) ;
      return new expression(OP_FUNC,"pow",lgroup,0) ;
    }
    inline exprP pow(exprP p1, int i2) {
      exprList lgroup ;
      lgroup.push_back(p1) ;
      lgroup.push_back(e_int(i2)) ;
      return new expression(OP_FUNC,"pow",lgroup,0) ;
    }
    inline exprP pow(int i1, exprP p2) {
      exprList lgroup ;
      lgroup.push_back(e_int(i1)) ;
      lgroup.push_back(p2) ;
      return new expression(OP_FUNC,"pow",lgroup,0) ;
    }
  
    inline exprP sin(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"sin",lgroup,0) ;
    }
    inline exprP cos(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"cos",lgroup,0) ;
    }
    inline exprP tan(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"tan",lgroup,0) ;
    }
    inline exprP sinh(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"sinh",lgroup,0) ;
    }
    inline exprP cosh(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"cosh",lgroup,0) ;
    }
    inline exprP tanh(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"tanh",lgroup,0) ;
    }
    inline exprP exp(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"exp",lgroup,0) ;
    }
    inline exprP ln(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"ln",lgroup,0) ;
    }
    inline exprP log(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"log",lgroup,0) ;
    }
    inline exprP log10(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"log10",lgroup,0) ;
    }
    inline exprP sqrt(exprP p) {
      exprList lgroup ;
      lgroup.push_back(p) ;
      return new expression(OP_FUNC,"sqrt",lgroup,0) ;
    }
  }

}
    
#endif

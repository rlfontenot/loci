#ifndef EXPR_H
#define EXPR_H

#include <Tools/stream.h>
#include <Tools/cptr.h>

#include <list>
#include <string>

namespace Loci {
enum OpType {
    OP_SCOPE=0x000,
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
    OP_STRING, OP_INT, OP_ERROR,
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

    expression() : op(op_priv), expr_list(expr_list_priv), name(name_priv),int_val(int_val_priv) { op_priv = OP_ERROR ; }

    static exprP  create( std::istream &s, char closing ) ;
    static exprP  get_term( std::istream &s ) ;
    static exprP  get_name( std::istream &s ) ;
    static OpType get_oper( std::istream &s ) ;
    static OpType get_unary_oper( std::istream &s ) ;
    static exprP  expand_oper( std::istream &s, exprP &p ) ;

    void PrintOperation(std::ostream &s, std::string oper,
                        char poChar = '(', char pcChar = ')' ) const ;
  public:
    const OpType              &op ;
    const exprList            &expr_list ;
    const std::string         &name ;
    const int                 &int_val ;
    expression(OpType op, const string nm, 
	       const exprList &elist, int ival = 0) :
      op(op_priv), expr_list(expr_list_priv), name(name_priv),
      int_val(int_val_priv) 
      { op_priv = op; expr_list_priv = elist ; name_priv = nm; 
      int_val_priv = ival ;}
    void Print(std::ostream &s) const ;
    static exprP create(std::istream &s) ;
    static exprP create(const std::string &s) ;

} ;

typedef expression::exprP     exprP ;
typedef expression::exprList  exprList ;

inline ostream &operator<<(ostream &s, const exprP &exp) {
    exp->Print(s) ;
    return s;
}

inline istream &operator>>(istream &s, exprP &exp) {
    exp = expression::create(s) ;
    return s ;
}

exprList collect_associative_op(const exprP &e, const OpType op) ;

}
    
#endif

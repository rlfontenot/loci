#include <Tools/expr.h>
#include <Tools/parse.h>
#include <Tools/stream.h>

namespace Loci {
    
void expression::PrintOperation(ostream &s, std::string oper,
                                char poChar, char pcChar) const
{
    s << poChar ;
    if(expr_list.size() != 0) {
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

    case OP_STRING:
        s << '"' << name <<'"' ;
        break ;

    case OP_INT:
        s << int_val ;
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
        cout << "*ERR*" ;
        break ;
        
    default:
        cerr << "unexpected operation in void expression::Print(ostream &s)"
             << endl ;
        break ;
    }
}
    
OpType expression::get_oper(istream &s) {
    if(parse::get_token(s,"&&"))
      return OP_LOGICAL_AND ;
    if(parse::get_token(s,"||"))
      return OP_LOGICAL_OR ;
    if(parse::get_token(s,">>"))
      return OP_SHIFT_RIGHT ;
    if(parse::get_token(s,"<<"))
      return OP_SHIFT_LEFT ;
    if(parse::get_token(s,"->"))
      return OP_ARROW ;
    if(parse::get_token(s,"<="))
      return OP_LE ;
    if(parse::get_token(s,">="))
      return OP_GE ;
    if(parse::get_token(s,"=="))
      return OP_EQUAL ;
    if(parse::get_token(s,"!="))
      return OP_NOT_EQUAL ;
    if(parse::get_token(s,"::"))
      return OP_SCOPE ;
    if(parse::get_token(s,"<"))
      return OP_LT ;
    if(parse::get_token(s,">"))
      return OP_GT ;
    
    if(parse::get_token(s,"%"))
      return OP_MODULUS ;
    if(parse::get_token(s,"+"))
      return OP_PLUS ;
    if(parse::get_token(s,"-"))
      return OP_MINUS ;
    if(parse::get_token(s,"*"))
      return OP_TIMES ;
    if(parse::get_token(s,"/"))
      return OP_DIVIDE ;
    if(parse::get_token(s,"&"))
      return OP_AND ;
    if(parse::get_token(s,"|"))
      return OP_OR ;
    if(parse::get_token(s,"^"))
      return OP_EXOR ;
    if(parse::get_token(s,","))
      return OP_COMMA ;
    if(parse::get_token(s,"="))
      return OP_ASSIGN ;
    if(parse::get_token(s,":"))
      return OP_COLON ;
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
    const unsigned int mask = ~0xff ;
    while(ot != OP_ERROR) {
        exprP p2 = expression::get_term(s) ;
        while(estack.size()>1 && ((ot&mask) >= (mask&(estack.back()->op))))
          estack.pop_back() ;
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


exprP expression::create(istream &s, char closing)
{
    exprP p = new expression ;
    parse::kill_white_space(s) ;
    if(s.eof() || s.peek() == EOF) {
        //        warn(closing!=';') ;
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
        //        warn(closing != ';') ;
        return p ;
    }
    if(s.peek() == closing) {
        s.get() ;
        return p ;
    }
    exprP p1 = p ;
    p = new expression ;
    p->op_priv = expression::get_oper(s) ;
    if(p->op == OP_ERROR) 
      return p ;
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

exprP expression::get_name(istream &s)
{
    if(parse::get_token(s,"("))
      return expression::create(s,')') ;
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
        } else
          return name ;
        
        exprP args = expression::create(s,closing) ;
        if(args->op == OP_COMMA)
          name->expr_list_priv = args->expr_list ;
        else
          name->expr_list_priv.push_back(args) ;
        return name ;
    }
    exprP p = new expression ;
    return p ;
}

exprP expression::get_term(istream &s)
{
    if(parse::get_token(s,"("))
      return expression::create(s,')') ;
    if(s.peek() == '-' || s.peek() == '+') {
        char ch = s.get() ;
        while(isspace(s.peek()))
          s.get() ;
        char ch2 = s.peek() ;
        s.putback(ch) ;
        if(isdigit(ch2)) {
            exprP ival = new expression ;
            ival->int_val_priv = parse::get_int(s) ;
            ival->op_priv = OP_INT ;
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
        return expression::get_name(s) ;
    }
    if (parse::is_int(s)) {
        exprP ival = new expression ;
        ival->int_val_priv = parse::get_int(s) ;
        ival->op_priv = OP_INT ;
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
    

}

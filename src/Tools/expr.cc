#include <Tools/expr.h>
#include <Tools/parse.h>
#include <Tools/stream.h>

namespace Loci {
  using std::cout ;    
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
        
  case OP_NIL:
    s <<"" ;
    break ;
    
    default:
        cerr << "unexpected operation in void expression::Print(ostream &s)"
             << endl ;
        break ;
    }
}
    
OpType expression::get_oper(istream &s) {
  if(parse::get_token(s,"&&")) {
    //cout << "in get oper && = " << OP_LOGICAL_AND << endl ; 
    return OP_LOGICAL_AND ;
  }
  if(parse::get_token(s,"||")) {
    //cout << "in get oper || = " << OP_LOGICAL_OR << endl ; 
    return OP_LOGICAL_OR ;
  }
  if(parse::get_token(s,">>")) {
    //cout << "in get oper >> = " << OP_SHIFT_RIGHT << endl ; 
    return OP_SHIFT_RIGHT ;
  }
  if(parse::get_token(s,"<<")) {
    //cout << "in get oper << = " << OP_SHIFT_LEFT << endl ; 
    return OP_SHIFT_LEFT ;
  }
  if(parse::get_token(s,"->")) {
    //cout << "in get oper ->  = " << OP_ARROW << endl ; 
    return OP_ARROW ;
  }
  if(parse::get_token(s,"<=")) {
    //cout << "in get oper  less than or equal to = " << OP_LE << endl ; 
    return OP_LE ;
  }
  if(parse::get_token(s,">=")) {
    //cout << "in get oper  greater than or equal to = " << OP_GE << endl ; 
    return OP_GE ;
  }
    if(parse::get_token(s,"==")) {
      //cout << "in get oper equals =  " << OP_EQUAL << endl ; 
      return OP_EQUAL ;
    }
    if(parse::get_token(s,"!=")) {
      //cout << "in get oper not_equals =  " << OP_NOT_EQUAL << endl ; 
      return OP_NOT_EQUAL ;
    }
    if(parse::get_token(s,"::")) {
      //cout << "in get oper scope =  " << OP_SCOPE << endl ; 
      return OP_SCOPE ;
    }
    if(parse::get_token(s,"<")) {
      //cout << "in get oper  less than = " << OP_LT << endl ; 
      return OP_LT ;
    }
    if(parse::get_token(s,">")) {
      //cout << "in get oper  greater than = " << OP_GT << endl ; 
      return OP_GT ;
    }
    
    if(parse::get_token(s,"%")) {
      //cout << "in get oper  modulus =  " << OP_MODULUS << endl ; 
      return OP_MODULUS ;
    }
    if(parse::get_token(s,"+")) {
      //cout << "in get oper  plus =  " << OP_PLUS << endl ; 
      return OP_PLUS ;
	} 
    
    if(parse::get_token(s,"-")) {
      //cout << "in get oper minus = " << OP_MINUS << endl ; 
      return OP_MINUS ;
    }
    if(parse::get_token(s,"*")) {
      //cout << "in get oper star =  " << OP_TIMES << endl ; 
      return OP_TIMES ;
    }
    if(parse::get_token(s,"/")) {
      //cout << "in get oper slash =  " << OP_DIVIDE << endl ; 
      return OP_DIVIDE ;
    }
    if(parse::get_token(s,"&")) {
      //cout << "in get oper & =  " << OP_AND << endl ; 
      return OP_AND ;
    }
    if(parse::get_token(s,"|")) {
      //cout << "in get oper | =  " << OP_OR << endl ; 
      return OP_OR ;
    }
    if(parse::get_token(s,"^")) {
      //cout << "in get oper ^ =  " << OP_EXOR << endl ; 
      return OP_EXOR ;
    }
      
    if(parse::get_token(s,",")) {
      //cout << "in get oper comma  = " << OP_COMMA << endl ; 
      return OP_COMMA ;
    }
    if(parse::get_token(s,"=")) {
      //cout << "in get oper assign  =  " << OP_ASSIGN << endl ; 
      return OP_ASSIGN ;
    }
    if(parse::get_token(s,":")) {
      //cout << "in get oper colon =" << OP_COLON << endl ; 
      return OP_COLON ;
    }
    //cout << "in get oper OP_ERROR  = " << OP_ERROR << endl ; 
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
  //cout << "In expand oper " << endl ;
  exprList estack ;
  estack.push_back(p) ;
  OpType ot = expression::get_oper(s) ;
  const unsigned int mask = ~0xff ;
  while(ot != OP_ERROR) {
    exprP p2 = expression::get_term(s) ;
    //cout << " in expand_oper p2 = " ;
    //p2->Print(cout) ;
    //cout << endl ;
    while(estack.size()>1 && ((ot&mask) >= (mask&(estack.back()->op)))) {
      //cout << " ot&mask = " << (ot&mask) << endl ;
      //cout << "mask&(estack.back()->op)  =  " << (mask&(estack.back()->op)) << endl ;
      //cout << " in expand oper popping out " ;
      //estack.back()->Print(cout)  ;
      //cout << endl ;
      estack.pop_back() ;
    }
    if(estack.back()->op == ot) {
      //cout << " ot&mask = " << (ot&mask) << endl ;
      //cout << "mask&(estack.back()->op)  =  " << (mask&(estack.back()->op)) << endl ;
      estack.back()->expr_list_priv.push_back(p2) ;
    } else if((mask&ot) < (mask&(estack.back()->op))) {
      //cout << "mask&ot less than mask&estack.back()->op" << endl ;
      //cout << " ot&mask = " << (ot&mask) << endl ;
      //cout << "mask&(estack.back()->op)  =  " << (mask&(estack.back()->op)) << endl ;
      exprP np = new expression() ;
      np->op_priv = ot ;
      np->expr_list_priv.push_back(estack.back()->expr_list.back()) ;
      //cout << " in expand oper after first push_back " ;
      //np->Print(cout) ;
      //cout << endl ;
      np->expr_list_priv.push_back(p2) ;
      //cout << " in expand oper after second push_back " ;
      //np->Print(cout) ;
      //cout << endl ;
      estack.back()->expr_list_priv.back() = np ;
      estack.push_back(np) ;
    } else {
      //cout << "mask&ot greater than mask&estack.back()->op" << endl ;
      //cout << " ot&mask = " << (ot&mask) << endl ;
      //cout << "mask&(estack.back()->op)  =  " << (mask&(estack.back()->op)) << endl ;
      exprP np = new expression() ;
      np->op_priv = ot ;
      np->expr_list_priv.push_back(estack.back()) ;
      //cout << " in expand oper after first push_back " ;
      //np->Print(cout) ;
      //cout << endl ;
      np->expr_list_priv.push_back(p2) ;
      //cout << " in expand oper after second push_back " ;
      //np->Print(cout) ;
      //cout << endl ;
      estack.back() = np ;
    }
    ot = expression::get_oper(s) ;
  }
  //cout << " in expand oper returning  " ;
  //estack.front()->Print(cout) ;
  //cout << endl ;
  return estack.front() ;
}


exprP expression::create(istream &s, char closing) {
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
  //cout << " in  create  p1 = "  ;
  //p1->Print(cout) ;
  //cout << endl ;
  p = new expression ;
  p->op_priv = expression::get_oper(s) ;
  if(p->op == OP_ERROR) 
    return p ;
  p->expr_list_priv.push_back(p1) ;
  //cout << " in create after pushing back p1  =  "  ;
  //p->Print(cout) ;
  //cout << endl ; 
  p->expr_list_priv.push_back(expression::get_term(s)) ;
  //cout << " in create after pushing back expression::get_term  =  "  ;
  //p->Print(cout) ;
  //cout << endl ; 
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
    //cout << "in get_name calling expression::create " << endl ;
    return expression::create(s,')') ;
  }
  if(parse::is_name(s)) {
    exprP name = new expression ;
    name->op_priv = OP_NAME ;
    name->name_priv = parse::get_name(s) ;
    char closing = 0 ;
    if(parse::get_token(s,"(")) {
      name->op_priv = OP_FUNC ;
      //cout << " in get_name oper =  OP_FUNC" << endl ; 
      closing = ')' ;
    } else if(parse::get_token(s,"[")) {
      name->op_priv = OP_ARRAY ;
      //cout << " in get_name oper =  OP_ARRAY" << endl ; 
      closing = ']' ;
    } else if(parse::get_token(s,"{")) {
      name->op_priv = OP_NAME_BRACE ;
      //cout << " in get_name oper =  OP_NAME_BRACE" << endl ; 
      closing = '}' ;
    }
    else
      return name ;
    exprP args = expression::create(s,closing) ;
    
    if((closing == ')' )&&(parse::get_token(s,"{"))) {
      name->op_priv = OP_FUNC_BRACE;
      //cout << "setting the op_priv as op_func_brace instead of op_func" << endl ;
      
      name->expr_list_priv.push_back(args) ;
      closing = '}' ;
      exprP brace_args = expression::create(s, closing) ;
      //cerr << "brace_args = " ;
      //brace_args->Print(cerr) ;
      //cerr << endl ;
      
      name->expr_list_priv.push_back(brace_args) ;
    }
    else {
      if(args->op == OP_COMMA)
	name->expr_list_priv = args->expr_list ;
      else
	name->expr_list_priv.push_back(args) ;
    }
    //cout << "in get_name returning  " ;
    //name->Print(cout) ;
    //cout << endl ;
    return name ;
  }
  exprP p = new expression ;
  return p ;
}

exprP expression::get_term(istream &s)
{
  if(parse::get_token(s,"(")) {
    //cout << "in get_term returning expression create " ;
    exprP temp = (expression::create(s,')')) ;
    //temp->Print(cout) ;
    //cout << endl ;
    return temp ;
    //return expression::create(s,')') ;
  }
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
      //cout << "in get_term returning ival after + or minus =   ";
      //ival->Print(cout) ;
      //cout << endl ;
      return ival ;
    }
  }
  
  OpType ot = get_unary_oper(s) ;
  if(ot != OP_ERROR) {
    exprP p = new expression ;
    p->op_priv = ot ;
    p->expr_list_priv.push_back(expression::get_term(s)) ;
    //p->Print(cout) ;
    //cout << endl ;
    return p ;
  }
  if(parse::is_name(s)) {
    //cout << "in get_term returning after calling is_name  =   ";
    exprP temp  = (expression::get_name(s)) ;
    //temp->Print(cout) ;
    //cout << endl ;
    return temp ;
    //return expression::get_name(s) ;
    }
    if (parse::is_int(s)) {
        exprP ival = new expression ;
        ival->int_val_priv = parse::get_int(s) ;
        ival->op_priv = OP_INT ;
	//cout << "in get_term returning after calling is_int  =   ";
	//ival->Print(cout) ;
	//cout << endl ;
        return ival ;
    }
    if (parse::is_string(s)) {
      exprP sval = new expression ;
      sval->name_priv = parse::get_string(s) ;
      sval->op_priv = OP_STRING ;
      //cout << "in get_term returning sval ater calling is_string =   ";
      //sval->Print(cout) ;
      //cout << endl ;
      return sval ;
    }
    
    exprP error = new expression ;
    //cout << "returning error  " ;
    //error->Print(cout) ;
    //cout << endl ;
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

#include <variable.h>

#include <Tools/stream.h>

using std::vector ;
using std::string ;
using std::make_pair ;

namespace Loci {
time_ident::time_hierarchy *time_ident::thp = 0 ;

int time_ident::time_hierarchy::add_level(const string &lname,int level) {
    vector<int>::const_iterator vi ;
    for(vi=time_db[level].children.begin();
        vi!=time_db[level].children.end();++vi) {
        if(time_db[*vi].level_name == lname)
          break ;
    }
    if(vi==time_db[level].children.end()) {
        time_db.push_back(time_info(lname,level)) ;
        int n = time_db.size() - 1 ;
        time_db[level].children.push_back(n) ;
        return n ;
    }
    else
      return *vi ;
}


time_ident::time_ident(const exprP &exp) {
    create_thp() ;
    id = 0 ;
    if(exp->op == OP_NAME) 
      return ;
    if(exp->op != OP_NAME_BRACE) {
        cerr << "syntax error interpreting time list in expression "
             << exp << endl ;
        return ;
    }
        
    exprList::const_iterator bgn,end,ii ;
    bgn = exp->expr_list.begin() ;
    end = exp->expr_list.end() ;
    --end ;
    for(ii=bgn;ii!=end;++ii) {
        if(OP_NAME != (*ii)->op) {
            cerr << "syntax error interpreting time list in expression "
                 << exp << endl ;
            cerr << "expecting OP_NAME, got " << *ii << endl ;
            break ;
        } else
          id = thp->add_level((*ii)->name,id) ;
    }

    if(OP_NAME == (*end)->op) {
        id = thp->add_level((*ii)->name,id) ;
        return ;
    }
    
    if(OP_PLUS == (*end)->op ||
       OP_MINUS == (*end)->op ||
       OP_ASSIGN == (*end)->op) {
        if((*end)->expr_list.size() != 2) {
            cerr << "syntax error in time label of name expression "
                 << exp << endl ;
            cerr << "error found while parsing " << *end << endl ;
            id = 0 ;
            return ;
        }
        exprP tname = (*end)->expr_list.front() ;
        if(OP_NAME != tname->op) {
            cerr << "syntax error in time label of name expression "
                 << exp << endl ;
            cerr << "time level should be of form NAME <op> int" << endl;
            id = 0 ;
            return ;
        } else {
            id = thp->add_level(tname->name,id) ;
            return ;
        }
    }

    cerr << "syntax error in time label of name expression "
         << exp << endl ;
    id = 0 ;

}

ostream &time_ident::Print(ostream &s) const {
    if(id==0)
      return s ;
    vector<int> v ;
    int l = id ;
    while(l != 0) {
        v.push_back(l) ;
        l = thp->parent(l) ;
    }
    s << thp->name(v.back()) ;
    v.pop_back() ;
    while(v.size() != 0) {
        s << "," << thp->name(v.back()) ;
        v.pop_back() ;
    }
    return s ;
}

bool time_ident::before(const time_ident &t) const {
  int i = t.id ;
  while(i != 0) {
    i = thp->parent(i) ;
    if(i == id)
      return true ;
  }
  return false ;
}

vector<time_ident> time_ident::children() {
  const vector<int> &vc = thp->get_children(id) ;
  vector<time_ident> v ;
  for(int i=0;i<vc.size();++i)
    v.push_back(time_ident(vc[i])) ;
  return v ;
}


variable::variable_db *variable::vdb = 0 ;

bool variable::info::operator<(const info &v) const {
    if(tvar != v.tvar)
      return tvar ;
    else if(assign != v.assign)
      return assign ;
    else
      return  name < v.name                               ||
          (name    == v.name    && ( time_id < v.time_id  ||
          (time_id == v.time_id && (  offset < v.offset   ||
          ( offset == v.offset  &&  priority < v.priority))))) ; 
}

bool variable::info::operator==(const info &v) const {
    return
           tvar == v.tvar    &&
         assign == v.assign  &&
           name == v.name    &&
        time_id == v.time_id &&
         offset == v.offset ;
}

variable::variable(const exprP &p) {
    create_vdb() ;
    info v ;
    exprP e = p ;
    

    if(OP_SCOPE == e->op) {
        exprList l = collect_associative_op(e,OP_SCOPE) ;
        while(l.begin() != l.end()) {
            exprP s = l.front() ;
            l.pop_front() ;
            if(l.begin() == l.end()) 
              e = s ;
            else 
              switch(s->op) {
              case OP_NAME:
                  v.priority.push_back(s->name) ;
                  break ;
              default:
                  cerr << "unable to interpret priority list in expression "
                       << s << endl
                       << "error occured while parsing " << e << endl ;
                    
              }
        }
    }
    
    if(e->op == OP_DOLLAR) {
        v.tvar = true ;
        e = e->expr_list.front() ;
    }
    exprP end ;
    switch(e->op) {
    case OP_NAME:
        v.name = e->name ;
        break ;
    case OP_NAME_BRACE:
        v.name = e->name ;
        v.time_id = time_ident(e) ;
        end = e->expr_list.back() ;
        if(OP_PLUS == end->op ||
           OP_MINUS == end->op ||
           OP_ASSIGN == end->op) {
            if(end->expr_list.size() != 2) {
                cerr << "syntax error in time label of name expression "
                     << p << endl ;
                cerr << "error found while parsing " << end << endl ;
                break ;
            }
            exprP tname = end->expr_list.front() ;
            int time_offset = 0 ;
            if(OP_NAME != tname->op) {
                cerr << "syntax error in time label of name expression "
                     << p << endl ;
                cerr << "time level should be of form NAME <op> int" << endl;
                break ;
            } else
                 if(OP_INT == end->expr_list.back()->op) {
                time_offset = end->expr_list.back()->int_val ;
            } else if((OP_UNARY_PLUS  == end->expr_list.back()->op ||
                       OP_UNARY_MINUS == end->expr_list.back()->op) &&
                      OP_INT == end->expr_list.back()->expr_list.front()->op){
                time_offset = end->expr_list.back()->
                              expr_list.front()->int_val ;
                if(OP_UNARY_MINUS == end->expr_list.back()->op)
                  time_offset = -time_offset ;
            } else {
                cerr << "time level should be of form NAME <op> int in expr "
                     << p << endl ;
                cerr << "got " << end << endl;
                break ;
            }
            if(OP_MINUS == end->op)
              time_offset = -time_offset ;
            
            v.offset = time_offset ;
            
            v.assign = false ;
            if(OP_ASSIGN == end->op) 
              v.assign = true ;
        }
        break ;
    default:
        cerr << "unable to interpret expression " << e << endl ;
        break ;
    }
    id = vdb->vars.get_id(v) ; 
}

variable::variable(const time_ident &t) {
    info v ;
    v.tvar = true ;
    v.time_id = t ;
    v.name = t.level_name() ;
    create_vdb() ;
    id = vdb->vars.get_id(v) ;
}

variable::variable(const variable &v, const time_ident &t) {
    create_vdb() ;
    info v2 = v.get_info() ;
    v2.time_id = t ;
    id = vdb->vars.get_id(v2) ;
}

variable variable::info::parent() const {
  info vi = *this ;
  vi.time_id = time_id.parent() ;
  return variable(vi) ;
}

variable variable::info::drop_assign() const {
  info vi = *this ;
  vi.assign = false ;
  return variable(vi) ;
}

variable variable::info::drop_priority() const {
  info vi = *this ;
  for(int i = 0;i<vi.priority.size()-1;++i)
    vi.priority[i] = vi.priority[i+1] ;
  vi.priority.pop_back() ;
  return variable(vi) ;
}

variable variable::info::new_offset(int o) const {
  info vi = *this ;
  vi.offset = o ;
  return variable(vi) ;
}

ostream &variable::info::Print(ostream &s ) const {
    if(tvar)
      s << "$" ;
    if(priority.begin() != priority.end()) {
        for(vector<string>::const_iterator i =priority.begin();
            i!=priority.end();++i)
          s<< *i <<"::" ;
    }

    s<< name ;
    if(time_id != time_ident()) {
        s << "{" << time_id ;
        if(assign)
          s << "=" << offset ;
        else if(offset != 0) {
            if(offset>0) 
              s << "+" << offset ;
            else
              s << offset ;
        }
        s << "}" ;
    }
    return s;
}

variableSet::variableSet(const exprP &e) {
    exprList l = collect_associative_op(e,OP_COMMA) ;
    
    for(exprList::const_iterator i=l.begin();i!=l.end();++i)
      *this += variable(*i) ; 
}

ostream &variableSet::Print(ostream &s) const
{
    variableSet::const_iterator i = begin() ;
    int sz = size() ;
    if(sz == 0)
      s << "()" ;
    else if(sz == 1)
        s << *i ;
    else {
      s << "(" << *i ;
      ++i ;
      for(;i!=end();++i) 
        s << "," << *i ;
      s << ")" ;
    }
    return s;
}


        
vmap_info::vmap_info(const exprP &e) {
    exprList l = collect_associative_op(e,OP_ARROW) ;
    exprList v = collect_associative_op(l.back(),OP_COMMA) ;
    for(exprList::const_iterator i=v.begin();i!=v.end();++i) {
        if((*i)->op == OP_ASSIGN) {
            variable varid((*i)->expr_list.front()) ;
            variable asnid((*i)->expr_list.back()) ;
                                     
            var += varid ;
            assign.push_back(make_pair(varid,asnid)) ;
        } else {
            var += variable(*i) ;
        }
    }

    l.pop_back() ;
    for(exprList::const_iterator j = l.begin();j!=l.end();++j) 
      mapping.push_back(variableSet(*j)) ;
}

ostream &vmap_info::Print(ostream &s) const {
    for(int j=0;j<mapping.size();++j) 
      s << mapping[j] << "->" ;
    s << var ;
    return s ;
}

//ostream &operator<<(ostream &s, const set<vmap_info> &v) {
//    set<vmap_info>::const_iterator i ;
//    for(i = v.begin();i!=v.end();) {
//        s << (*i) ;
//        ++i ;
//        if(i!=v.end())
//          s << "," ;
//    }
//    return s ;
//}



}

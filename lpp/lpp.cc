#include "lpp.h"
#include <ctype.h>
#include <set>

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
using namespace Loci ;

bool is_name(istream &s) {
  int ch = s.peek() ;
  return isalpha(ch) || ch == '_' ;
}
    
string get_name(istream &s) {
  if(!is_name(s))
    throw parseError("expected name") ;
  string str ;
  while(!s.eof() && (s.peek() != EOF) &&
        (isalnum(s.peek()) || (s.peek() == '_')) )
    str += s.get() ;
  
  return str ;
}

bool is_string(istream &s) {
  return s.peek() == '\"' ;
}
    
string get_string(istream &s) {
  if(!is_string(s))
    throw parseError("expected string") ;
  string str ;
  if(s.eof())
    throw parseError("unexpected EOF") ;
  s.get() ;
  int ch = s.get() ;
  while(ch != '\"' &&!s.eof()) {
    str += ch ;
    ch = s.get() ;
  }
  if(ch!='\"')
    throw parseError("no closing \" for string") ;
  return str ;
}

class parsebase {
public:
  int lines ;
  parsebase() {lines = 0 ; }
  istream &killsp(istream &s) {
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      s.get();
    }
    return s ;
  }
} ;
  
template<class T> class funclist : public parsebase {
public:
  list<T> flist ;
  istream &get(istream &s) {
    killsp(s) ;
    if(s.peek() != '(')
      return s ;
    char c = s.get();
    killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      killsp(s) ;
      if(s.peek() == ')') {
        c = s.get() ;
        return s ;
      }
      if(s.peek() != ',') {
        throw parseError("syntax error") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "(" ;
      typename list<T>::const_iterator ii ;
      ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
        s+="," ;
        s+= ii->str() ;
      }
      s += ")" ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    typename list<T>::const_iterator ii ;
    for(ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;
  

template<class T> class templlist : public parsebase {
public:
  list<T> flist ;
  istream &get(istream &s) {
    killsp(s) ;
    if(s.peek() != '<')
      return s ;
    char c ;
    s.get(c);
    killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      killsp(s) ;
      if(s.peek() == '>') {
        s.get() ;
        return s ;
      }
      if(s.peek() != ',') {
        throw parseError("syntax error, expected comma") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "<" ;
      typename list<T>::const_iterator ii ;
      ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
        s+="," ;
        s+= ii->str() ;
      }
      s += "> " ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    typename list<T>::const_iterator ii ;
    for(ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;

class typestuff : public parsebase {
public:
  string name ;
  templlist<typestuff> templ_args ;
  istream &get(istream &s) {
    killsp(s) ;
    if(is_name(s)) {
      name = get_name(s) ;
    } else if(isdigit(s.peek())) {
      while(isdigit(s.peek())) {
        name += s.get() ;
      }
    } else
      throw parseError("syntax error") ;
    templ_args.get(s) ;
    return s ;
  }
  string str() const {
    string s ;
    s+= name ;
    s+= templ_args.str() ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i+= templ_args.num_lines() ;
    return i ;
  }
} ;
class bracestuff : public parsebase {
public:
  string stuff ;
  istream &get(istream &s) {
    killsp(s) ;
    if(s.peek() == '{') {
      char c = s.get() ;
      while(s.peek() != EOF && s.peek() != '}') {
        c = s.get() ;
        if(c == '{')
          throw parseError("syntax error") ;
        stuff += c ;
      }
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      c = s.get() ;
    }
    return s ;
  }
    
  string str() const {
    string s ;
    if(stuff == "")
      return s ;
    s += "{" ;
    s += stuff ;
    s += "}" ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    return i ;
  }
} ;
  

class var : public parsebase {
public:
  bool isdollar ;
  string name ;
  list<string> prio_list ;
  funclist<var> param_args ;
  bracestuff bs ;
  var() : isdollar(false) {}
  
  istream &get(istream &s) {
    isdollar = false ;
    killsp(s) ;
    if(s.peek() == '$') {
      s.get() ;
      isdollar=true ;
    }
    if(!is_name(s))
      throw parseError("syntax error") ;
    name = get_name(s) ;
    killsp(s) ;
    if(s.peek() == ':') {
      while(s.peek() == ':') {
        s.get() ;
        if(s.peek() != ':')
          throw parseError("syntax error") ;
        s.get() ;
        killsp(s) ;
        prio_list.push_back(name);
        if(!is_name(s)) 
          throw parseError("syntax error") ;
        name = get_name(s) ;
        killsp(s) ;
      }
    }
          
    param_args.get(s) ;
    bs.get(s) ;

    return s ;
  }
  string str() const {
    string s ;
    list<string>::const_iterator li ;
    for(li=prio_list.begin();li!=prio_list.end();++li)
      s+= *li + "::" ;
    if(isdollar)
      s+="$" ;
    s+=name ;
    s+= param_args.str() ;
    s+= bs.str() ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i += param_args.num_lines() ;
    i += bs.num_lines() ;
    return i ;
  }
} ;

class nestedparenstuff : public parsebase {
public:
  string paren_contents ;
  istream &get(istream &s) {
    killsp(s) ;
    if(s.peek() != '(')
      throw parseError("syntax error, expecting '('") ;
    s.get() ;
    int open_parens = 0 ;
    while(s.peek() != ')' || open_parens != 0) {
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '(')
        open_parens++ ;
      if(s.peek() == ')')
        open_parens-- ;
      if(s.peek() == '\n' || s.peek() == '\r') {
        s.get() ;
        lines++ ;
        continue ;
      }
      paren_contents += s.get() ;
    }
    s.get() ;
    return s ;
  }
  string str() {
    return paren_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

class nestedbracketstuff : public parsebase {
public:
  string bracket_contents ;
  istream &get(istream &s) {
    killsp(s) ;
    if(s.peek() != '[')
      throw parseError("syntax error, expecting '['") ;
    s.get() ;
    int open_brackets = 0 ;
    while(s.peek() != ']' || open_brackets != 0) {
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '[')
        open_brackets++ ;
      if(s.peek() == ']')
        open_brackets-- ;
      if(s.peek() == '\n' || s.peek() == '\r') {
        s.get() ;
        lines++ ;
        continue ;
      }
      bracket_contents += s.get() ;
    }
    s.get() ;
    return s ;
  }
  string str() {
    return bracket_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

void parseFile::setup_Type(std::ostream &outputFile) {
  var vin ;
  vin.get(is) ;
  typestuff tin ;
  tin.get(is) ;
  while(is.peek() == ' ' || is.peek() == '\t') 
    is.get() ;
  if(is.peek() != ';')
    throw parseError("syntax error, missing ';'") ;
  is.get() ;
  variable v(vin.str()) ;
  outputFile << "// $type " << v << ' ' << tin.str() ;
  int nl = vin.num_lines()+tin.num_lines() ;
  line_no += nl ;
  for(int i=0;i<nl;++i)
    outputFile << endl ;
  type_map[v] = pair<string,string>(tin.name,tin.templ_args.str()) ;
}

namespace {
  inline void fill_descriptors(set<vmap_info> &v, const exprList &in) {
    
    for(exprList::const_iterator i = in.begin();i!=in.end();++i) {
      vmap_info di(*i) ;
      if(v.find(di) != v.end())
        cerr << "Warning, duplicate variable in var set." << endl ;
      else
        v.insert(di) ;
    }
  }
}


void parseFile::process_Compute(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    void compute(const sequence &seq) { " ;
  is.get() ;
  
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '$') {
      string name ;
      variable v ;
      is.get() ;
      var vin ;
      vin.get(is) ;
      v = variable(vin.str()) ;

      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      outputFile << "(*" << vmi->second << ')' ;
    }
    char c = is.get() ;
    outputFile << c ;
  } ;
}

void parseFile::process_Calculate(std::ostream &outputFile,
                                  const map<variable,string> &vnames) {
  outputFile << "    void calculate(Entity _e_) { " << endl ;
  is.get() ;
  while(is.peek() == ' ' || is.peek() == '\t')
    is.get() ;
  if(is.peek() == '\n') {
    is.get() ;
  }
  syncFile(outputFile) ;
  killspout(outputFile) ;
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF in process_Calculate") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '/') {
      is.get() ;
      outputFile << '/' ;
      if(is.peek() == '/') { // comment line
        is.get() ;
        outputFile << '/' ;
        while(is.peek() != '\n') {
          char c = is.get() ;
          outputFile << c ;
        }
        killspout(outputFile) ;
      }
      continue ;
    }
          
    
    if(is_name(is) || is.peek() == '$') {
      bool first_name = is_name(is) ;
      string name ;
      variable v ;
      if(first_name) 
        name = get_name(is) ;
      else {
        is.get() ;
        var vin ;
        vin.get(is) ;
        v = variable(vin.str()) ;
      }
      list<variable> vlist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        killsp() ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          killsp() ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
          } else {
            dangling_arrow = true ;
            break ;
          }
        } else {
          is.unget() ;
          break ;
        }
      }
      if(dangling_arrow && first_name) {
        outputFile << name << " ->" ;
        continue ;
      }
      if(dangling_arrow)
        throw parseError("syntax error, near '->' operator") ;

      if(first_name && (vlist.size() == 0)) {
        outputFile << name << ' ';
        continue ;
      }
      list<variable>::reverse_iterator ri ;
      for(ri=vlist.rbegin();ri!=vlist.rend();++ri) {
        map<variable,string>::const_iterator vmi = vnames.find(*ri) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << '[' ;
      }
      if(first_name) {
        outputFile << '*' << name ;
      } else {
        map<variable,string>::const_iterator vmi = vnames.find(v) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << v << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << "[_e_]" ;
      }
      for(size_t i=0;i<vlist.size();++i)
        outputFile << ']' ;
    }
    char c = is.get() ;
    outputFile << c ;
  } ;
}

void parseFile::setup_Rule(std::ostream &outputFile) {
  killsp() ;
  string rule_type ;
  if(is_name(is)) {
    rule_type = get_name(is) ;
  } else 
    throw parseError("syntax error") ;
  nestedparenstuff signature ;

  signature.get(is) ;
  line_no += signature.num_lines() ;
  nestedbracketstuff apply_op ;
  killsp() ;
  if(rule_type == "apply") {
    if(is.peek() != '[') 
      throw parseError("apply rule missing '[operator]'") ;
    apply_op.get(is) ;
    line_no += apply_op.num_lines() ;
    killsp() ;
  }
  

  string constraint, conditional ;
  
  while(is.peek() == ',') {
    is.get() ;
    killsp() ;
    if(!is_name(is))
      throw parseError("syntax error") ;

    string s = get_name(is) ;
    if(s == "constraint") {
      nestedparenstuff con ;
      con.get(is) ;
      constraint = con.str() ;
      line_no += con.num_lines() ;
    }
    if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      conditional = con.str() ;
      line_no += con.num_lines() ;
    }
    killsp() ;
  }
  
  if(is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;

  string sig = signature.str() ;
  string heads,bodys ;
  exprP head=0,body=0 ;
  for(size_t i=0;i<sig.size()-1;++i) {
    if(sig[i]=='<' && sig[i+1]=='-') {
      heads = sig.substr(0,i) ;
      bodys = sig.substr(i+2,sig.size()) ;
      head = expression::create(heads) ;
      body = expression::create(bodys) ;
    }
  }
  if(head == 0) {
    heads = sig ;
    head = expression::create(heads) ;
  }
  
  string class_name = "file_" ;
  for(size_t i=0;i<filename.size();++i) {
    char c = filename[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      class_name += c ;
    if(c == '.')
      break ;
  }
  class_name += "_rule_" ;
  if(conditional != "")
    sig += "_" + conditional ;
  if(constraint != "")
    sig += "_" + constraint ;
  for(size_t i=0;i<sig.size();++i) {
    if(isalpha(sig[i]) || isdigit(sig[i]))
      class_name += sig[i] ;
    if(sig[i] == ',' || sig[i] == '-' || sig[i] == '>' || sig[i] == '('||
       sig[i] == ')' || sig[i] == '{' || sig[i] == '}' || sig[i] == '='||
       sig[i] == '+' || sig[i] == '_')
      class_name += '_' ;
  }

  set<vmap_info> sources ;
  set<vmap_info> targets ;
  if(body != 0)
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;

  set<vmap_info>::const_iterator i ;
  variableSet input,output ;
  for(i=sources.begin();i!=sources.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    input += i->var ;
  }

  for(i=targets.begin();i!=targets.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    output += i->var ;
  }

  map<variable,string> vnames ;
  variableSet::const_iterator vi ;
  variableSet all_vars = input;
  all_vars += output ;
  
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    string vn = (*vi).str() ;
    string name ;
    for(size_t si=0;si!=vn.size();++si) {
      if(isalpha(vn[si]) || isdigit(vn[si]) || vn[si] == '_')
        name += vn[si] ;
      if(vn[si]=='{' || vn[si] == '}')
        name += '_' ;
      if(vn[si]=='=')
        name += "_EQ_" ;
      if(vn[si]=='+')
        name += "_P_" ;
      if(vn[si]=='-')
        name += "_M_" ;
    }
    vnames[*vi] = name ;
  }
  map<variable,pair<string,string> > local_type_map ;
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    if(v.is_time_variable()) {
      local_type_map[v] = pair<string,string>("param","<int> ") ;
      continue ;
    }
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;

    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = type_map.find(v)) == type_map.end()) {
      v = v.new_offset(0) ;
      v = v.drop_assign() ;
      while(v.time() != time_ident())
        v = v.parent() ;
      
      if((mi = type_map.find(v)) == type_map.end()) {
        while(v.get_info().namespac.size() != 0)
          v = v.drop_namespace() ;
        if((mi = type_map.find(v)) == type_map.end()) {
          cerr << "unable to determine type of variable " << *vi << endl ;
        }
      }
    }
    local_type_map[*vi] = mi->second ;
  }

  outputFile << "class " << class_name << " : public " << rule_type << "_rule" ;
  if(rule_type == "apply") {
    if(output.size() != 1) 
      throw parseError("apply rule should have only one output variable") ;
    variable av = *(output.begin()) ;
    pair<string,string> tinfo = local_type_map[av] ;
    outputFile << "< " << tinfo.first << tinfo.second <<","
               << apply_op.str() ;
    if(tinfo.first == "storeVec") {
      outputFile << "Vect<" << tinfo.second <<" > " ;
    } else if(tinfo.first == "storeMat") {
      outputFile << "Mat<" << tinfo.second <<" > " ;
    } else {
      outputFile << tinfo.second ;
    }
    outputFile << "> " ;
  }
  outputFile << " {" << endl ;
  syncFile(outputFile) ;
  
  input -= output ;
  variableSet ins = input ;
  for(vi=ins.begin();vi!=ins.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    outputFile << "    const_" << mi->second.first <<  mi->second.second ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  bool output_param = false ;
  for(vi=output.begin();vi!=output.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(mi->second.first == "param") {
      output_param= true ;
    }
    outputFile << "    " << mi->second.first <<  mi->second.second ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "public:" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    " << class_name << "() {" << endl ;
  syncFile(outputFile) ;
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
               << vnames[*vi] << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  outputFile <<   "       output(\"" << heads << "\") ;" << endl ;
  syncFile(outputFile) ;
  if(constraint!="") {
    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    syncFile(outputFile) ;
  }  
  outputFile <<   "    }" << endl ;
  syncFile(outputFile) ;

  if(rule_type == "singleton" ||
     rule_type == "optional"  ||
     rule_type == "default" ||
     (output_param && rule_type != "apply") ) {
    process_Compute(outputFile,vnames) ;
  } else {
    process_Calculate(outputFile,vnames) ;

    outputFile <<   "    void compute(const sequence &seq) { " << endl ;
    syncFile(outputFile) ;
    outputFile <<   "      do_loop(seq,this) ;" << endl ;
    syncFile(outputFile) ;
    outputFile <<   "    }" << endl ;
    syncFile(outputFile) ;
  }
  outputFile <<   "} ;" << endl ;
    syncFile(outputFile) ;
  outputFile << "register_rule<"<<class_name<<"> register_"<<class_name
             << " ;" << endl ;
  syncFile(outputFile) ;
}

void parseFile::processFile(string file, ostream &outputFile) {
  filename = file ;
  line_no = 1 ;
  is.open(file.c_str(),ios::in) ;
  if(is.fail()) {
    throw parseError("can't open") ;
  }
  char c ;
  syncFile(outputFile) ;
  do {
    while(is.peek() == ' ' || is.peek() == '\t') {
      is.get(c) ;
      outputFile << c ;
    }
    try {
      if(is.peek() == '$') { // Loci specific code!
        is.get(c) ; // get the $
        if(is_name(is)) {
          std::string key = get_name(is) ;
          if(key == "type") {
            setup_Type(outputFile) ;
          } else if(key == "rule") {
            setup_Rule(outputFile) ;
          } else if(key == "include") {
            killsp() ;
            if(!is_string(is))
              throw parseError("syntax error") ;
            string newfile = get_string(is) ;
            
            parseFile parser ;
            parser.processFile(newfile,outputFile) ;
            syncFile(outputFile) ;
            map<variable,pair<string,string> >::const_iterator mi ;
            for(mi=parser.type_map.begin();mi!=parser.type_map.end();++mi)
              type_map[mi->first] = mi->second ;
            
            
          } else {
            
            cerr << filename << ':'<<line_no << " Loci Preprocessor key " <<
              key << " unknown!" << endl ;
          }
        } else {
          cerr << filename << ':' << line_no << " syntax error" << endl ;
          string s ;
          is >> s ;
        }
      } else {
        while(is.peek() != '\n' && is.peek() != EOF) {
          is.get(c) ;
          outputFile << c ;
        }
        is.get(c) ;
        outputFile << endl ;
        line_no++ ;
      }
      if(is.peek() == EOF)
        is.get(c) ;
    }
    catch(parseError pe) {
      cerr << filename << ':' << line_no << ": " << pe.error_type << endl ;
      char buf[512] ;
      is.getline(buf,512) ;
      line_no++ ;
      cerr << "remaining line = " << buf << endl ;
    }
  } while(!is.eof()) ;
    
}

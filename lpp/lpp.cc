#include "lpp.h"
#include <ctype.h>
#include <set>
#include <iostream>
#include <sstream>

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

bool is_comment(istream &s) {
  if(s.peek() != '/')
    return false ;
  s.get() ;
  char c = s.peek() ;
  s.unget() ;
  if(c == '/' || c == '*')
    return true ;
  return false ;
}

istream &killComment(istream &s, int & lines) {
  s.get() ;
  char c = s.get()  ;
  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n')
      s.get() ;
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
    }
    return s ;
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}
    
istream &killsp(istream &s, int &lines) {

  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      s.get();
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killComment(s,lines) ;
      foundstuff = true ;    }
  } while(foundstuff) ;
  return s ;
}

istream &killCommentOut(istream &s, int & lines,ostream &out) {
  s.get() ;
  out << '/' ;
  char c = s.get()  ;
  out << c ;
  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n') {
      char c = s.get() ;
      out << c ;
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
      out << '\n' ;
    }
    return s ;
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    out << c ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        out << '/' ;
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}
    
istream &killspOut(istream &s, int &lines, ostream &out) {

  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      char c = s.get();
      out << c ;
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killCommentOut(s,lines,out) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  return s ;
}

int parseFile::killsp() {
  int l = line_no ;
  ::killsp(is,line_no) ;
  return l-line_no ;
}

int parseFile::killspout(std::ostream &outputFile) {
  int l = line_no ;
  ::killspOut(is,line_no,outputFile) ;
  return l-line_no ;
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
    if(isalpha(s.peek()) || s.peek() == '_') {
      char c = s.peek() ;
      while(isalpha(c = s.peek()) || isdigit(c) || c == '_' ||  c == ':')
        name += s.get() ;
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

variable convertVariable(variable v) {
  variable::info vinfo = v.get_info() ;
  vinfo.priority = std::vector<std::string>() ;
  for(size_t i=0;i<vinfo.v_ids.size();++i) {
    std::ostringstream ss ;
    ss << 'X' << i << endl ;
    variable xi = variable(ss.str()) ;
    vinfo.v_ids[i] = xi.ident() ;
  }
  return variable(vinfo) ;
}

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
  v = convertVariable(v) ;
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


void parseFile::process_Prelude(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    void prelude(const Loci::sequence &seq) { " ;
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
      outputFile << vmi->second  ;
    }
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } ;
}

void parseFile::process_Compute(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    void compute(const Loci::sequence &seq) { " ;
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
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
      continue ;
    }      
    if(is.peek() == '$') {
      variable v ;
      is.get() ;
      bool deref = true ;
      if(is.peek() == '*') {
        is.get() ;
        deref = false ;
      }
      

      var vin ;
      vin.get(is) ;
      v = variable(vin.str()) ;

      map<variable,string>::const_iterator vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      map<variable,pair<string,string> >::const_iterator mi ;
      if(((mi = type_map.find(v)) != type_map.end()) &&
         (mi->second.first == "Constraint" || !deref)) {
        outputFile << vmi->second  ;
      } else {
        outputFile << "(*" << vmi->second << ')' ;
      }
    }
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } ;
}

string getNumber(std::istream &is) {
  string num ;
  while(isdigit(is.peek()))
    num+= is.get();
  if(is.peek() == '.') {
    num += is.get() ;
    while(isdigit(is.peek()))
      num+= is.get();
  }
  if(is.peek() == 'e' || is.peek() == 'E') {
    num += is.get() ;
    if(is.peek() == '-' || is.peek() == '+')
      num += is.get() ;
    while(isdigit(is.peek()))
      num += is.get() ;
  }
  return num ;
}

string parseFile::process_String(string in,
                                 const map<variable,string> &vnames) {
  ostringstream outputFile ;
  istringstream is(in) ;

  int line_no = 0 ;

  for(;;) {
    ::killspOut(is,line_no,outputFile) ;

    if(is.peek() == EOF)
      break ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      continue ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      continue ;
    }
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
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
        ::killspOut(is,line_no,outputFile) ;
      }
      continue ;
    }
          
    if(is.peek() == '#') {
      is.get() ;
      outputFile << '#' ;
      while(is.peek() != '\n') {
        char c = is.get() ;
        outputFile << c ;
      }
      ::killspOut(is,line_no,outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }

    if(is_name(is) || is.peek() == '$') {
      bool first_name = is_name(is) ;
      string name ;
      variable v ;
      string brackets ;
      if(first_name) 
        name = get_name(is) ;
      else {
        is.get() ;
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;
          
          variable v(vin.str()) ;
          map<variable,string>::const_iterator vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }
          
          outputFile << vmi->second ;
          continue ;
        }
        
        var vin ;
        vin.get(is) ;
        v = variable(vin.str()) ;
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames) ;
          brackets = "[" + binfo + "]" ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          ::killspOut(is,line_no,outputFile) ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            ::killspOut(is,line_no,outputFile) ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames) ;
              brk = "[" + binfo +"]";
            }
            blist.push_back(brk) ;
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
        outputFile << name << ' ' ;
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
          throw parseError("type error: is this variable in the rule signature?") ;
        }
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      list<string>::const_iterator rbi ;
      for(rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

    }
    if(is.peek() != EOF) {
      char c = is.get() ;
      outputFile << c ;
    }
  } 

  
  return outputFile.str() ;
}
                                 
void parseFile::process_Calculate(std::ostream &outputFile,
                                  const map<variable,string> &vnames) {
  if(prettyOutput)
    outputFile << "    void calculate(Entity e) { " << endl ;
  else
    outputFile << "    void calculate(Entity _e_) { " << endl ;
  is.get() ;
  while(is.peek() == ' ' || is.peek() == '\t')
    is.get() ;
  if(is.peek() == '\n') {
    is.get() ;
    line_no++ ;
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
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
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
          
    if(is.peek() == '#') {
      is.get() ;
      outputFile << '#' ;
      while(is.peek() != '\n') {
        char c = is.get() ;
        outputFile << c ;
      }
      killspout(outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }

    if(is_name(is) || is.peek() == '$') {
      int lcount = 0 ;
      bool first_name = is_name(is) ;
      string name ;
      variable v ;
      string brackets ;
      if(first_name) 
        name = get_name(is) ;
      else {
        is.get() ;
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;
          line_no += vin.num_lines() ;
          lcount += vin.num_lines() ;
          
          variable v(vin.str()) ;
          map<variable,string>::const_iterator vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }
          
          outputFile << vmi->second ;
          continue ;
        }
        
        var vin ;
        vin.get(is) ;
        line_no += vin.num_lines() ;
        lcount += vin.num_lines();
        v = variable(vin.str()) ;
        killsp() ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames) ;
          brackets = "[" + binfo + "]" ;
          line_no += nb.num_lines() ;
          lcount += nb.num_lines() ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        lcount += killsp() ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          lcount += killsp() ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            lcount += killsp() ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames) ;
              brk = "[" + binfo +"]";
              line_no += nb.num_lines() ;
              lcount += nb.num_lines() ;
            }
            blist.push_back(brk) ;
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
        outputFile << name << ' ' ;
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
          throw parseError("type error: is this variable in the rule signature?") ;
        }
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      list<string>::const_iterator rbi ;
      for(rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

      for(int i=0;i<lcount;++i)
        outputFile << endl ;
    }
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } 
}


string var2name(variable v) {
  string vn = v.str() ;
  string name ;
  if(!prettyOutput)
    name += "L_" ;
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
  if(!prettyOutput)
    name += "_" ;
  return name ;
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
  string parametric_var ;
  list<string> options ;
  list<string> comments ;
  list<pair<variable,variable> > inplace ;
  
  bool use_prelude = false ;
  bool is_specialized = false ;
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
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      parametric_var = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      conditional = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "inplace") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      exprP p = expression::create(ip.str()) ;
      exprList l = collect_associative_op(p,OP_OR) ;
      if(l.size() != 2) 
        throw parseError("inplace needs two variables with a '|' separator") ;
        
      exprList::const_iterator i = l.begin() ;
      variable v1(*i) ;
      ++i ;
      variable v2(*i) ;
      inplace.push_back(pair<variable,variable>(v1,v2)) ;
    
    } else if(s == "prelude") {
      use_prelude=true ;
      killsp() ;
      continue ;
    } else if(s == "specialized") {
      is_specialized = true ;
    } else if(s == "option") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    killsp() ;
  }

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
    exprList tmp ;
    tmp = collect_associative_op(head,OP_COMMA) ; 
    if(tmp.size() > 1)
      throw parseError("rule signature missing <- operator!") ;
  }
  
  string class_name = "file_" ;
  for(size_t i=0;i<filename.size();++i) {
    char c = filename[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      class_name += c ;
    if(c == '.')
      break ;
  }
  class_name += '0' + (cnt/100)%10 ;
  class_name += '0' + (cnt/10)%10 ;
  class_name += '0' + (cnt)%10 ;
  cnt++ ;
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
    vnames[*vi] = var2name(*vi) ;
    if(vi->get_info().priority.size() != 0) {
      variable v = *vi ;
      while(v.get_info().priority.size() != 0)
        v = v.drop_priority() ;
      vnames[v] = vnames[*vi] ;
    }
  }

  list<pair<variable,variable> >::const_iterator ipi ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
  }
  
  map<variable,pair<string,string> > local_type_map ;
  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    if(v.is_time_variable()) {
      local_type_map[v] = pair<string,string>("param","<int> ") ;
      continue ;
    }
    v = convertVariable(v) ;
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = type_map.find(v)) == type_map.end()) {

      v = v.new_offset(0) ;
      v = v.drop_assign() ;
      while(v.time() != time_ident())
        v = v.parent() ;
      debugout << "dropped v = " << v << endl ;
      if((mi = type_map.find(v)) == type_map.end()) {
        while(v.get_info().namespac.size() != 0)
          v = v.drop_namespace() ;
        if((mi = type_map.find(v)) == type_map.end()) {
          string s ;
          s = "unable to determine type of variable " ;
          s += v.str() ;
          throw parseError(s) ;
        }
      }
    }

    local_type_map[*vi] = mi->second ;
  }

  if(!prettyOutput)
    outputFile << "namespace {" ;
  outputFile << "class " << class_name << " : public Loci::" << rule_type << "_rule" ;
  if(rule_type == "pointwise") {
    for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi) {
      if(local_type_map[*vi].first == "param" && vi->get_info().name != "OUTPUT") {
        throw(parseError("pointwise rule cannot compute param, use singleton")) ;
      }
    }
  }
  if(rule_type == "singleton") {
    for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi) {
      string t = local_type_map[*vi].first ;
      if(t == "store" || t == "storeVec" || t == "multiStore") {
        throw(parseError("singleton rule cannot compute store's, use pointwise")) ;
      }
    }
  }
    
              

  if(rule_type == "apply") {
    if(output.size() != 1) 
      throw parseError("apply rule should have only one output variable") ;
    variable av = *(output.begin()) ;
    pair<string,string> tinfo = local_type_map[av] ;
    outputFile << "< " << tinfo.first << tinfo.second <<","
               << apply_op.str() ;
    if(tinfo.first == "storeVec") {
      outputFile << "<Vect" << tinfo.second <<" > " ;
    } else if(tinfo.first == "storeMat") {
      outputFile << "<Mat" << tinfo.second <<" > " ;
    } else {
      outputFile << tinfo.second ;
    }
    outputFile << "> " ;
  }
  outputFile << " {" << endl ;
  syncFile(outputFile) ;

  variableSet outs = output ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  variableSet ins = input ;
  ins -= outs ;
  for(vi=ins.begin();vi!=ins.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(!prettyOutput) 
      outputFile << "    Loci::const_" << mi->second.first <<  mi->second.second ;
    else 
      outputFile << "    const_" << mi->second.first <<  mi->second.second ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  bool output_param = false ;
  for(vi=outs.begin();vi!=outs.end();++vi) {
    map<variable,pair<string,string> >::const_iterator mi ;
    if((mi = local_type_map.find(*vi)) == local_type_map.end()) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(vi->get_info().name != "OUTPUT" && mi->second.first == "param") {
      output_param= true ;
    }
    if(!prettyOutput)
      outputFile << "    Loci::" << mi->second.first <<  mi->second.second ;
    else
      outputFile << "    " << mi->second.first <<  mi->second.second ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "public:" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    " << class_name << "() {" << endl ;
  syncFile(outputFile) ;
  for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
               << vnames[*vi] << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  for(i=targets.begin();i!=targets.end();++i) {
    outputFile <<   "       output(\"" ;
    for(size_t j=0;j<i->mapping.size();++j)
      outputFile << i->mapping[j] << "->" ;

    // Output target variables, adding inplace notation if needed
    variableSet::const_iterator vi ;
    if(i->var.size() > 1)
      outputFile << '(' ;
    for(vi=i->var.begin();vi!=i->var.end();++vi) {
      if(vi != i->var.begin())
        outputFile << ',' ;
      for(ipi=inplace.begin();ipi!=inplace.end();++ipi) {
        if((ipi->first) == *vi)
          break ;
      }
      if(ipi!=inplace.end()) {
        if(i->mapping.size() == 0 || i->var.size() > 1)
          outputFile << ipi->first << "=" << ipi->second ;
        else
          outputFile << '('<<ipi->first << "=" << ipi->second <<')';
      } else
        outputFile << *vi ;
    }
    if(i->var.size() > 1)
      outputFile << ')' ;

    outputFile <<  "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  //  outputFile <<   "       output(\"" << heads << "\") ;" << endl ;
  //  syncFile(outputFile) ;

  if(constraint!="") {
    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  if(parametric_var != "") {
    outputFile <<   "       set_parametric_variable(\""
               << parametric_var << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(is_specialized) {
    outputFile <<   "       set_specialized() ; " << endl ;
    syncFile(outputFile) ;
  }
  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  list<string>::const_iterator lsi ;
  for(lsi=options.begin();lsi!=options.end();++lsi) {
    string s = *lsi ;
    bool has_paren = false ;
    for(size_t i = 0;i<s.size();++i)
      if(s[i] == '(')
        has_paren = true ;
    outputFile <<   "       " << s ;
    if(!has_paren)
      outputFile << "()" ;
    outputFile << " ;" << endl;
    syncFile(outputFile) ;
  }
  for(lsi=comments.begin();lsi!=comments.end();++lsi) {
    outputFile <<   "       comments(" << *lsi << ") ;" << endl ;
    syncFile(outputFile) ;
  }

  outputFile <<   "    }" << endl ;
  syncFile(outputFile) ;

  bool use_compute = true ;

  if(use_prelude) {
    process_Prelude(outputFile,vnames) ;
    killsp() ;
    if(is.peek() == ';') {
      is.get() ;
      use_compute = false ;
    }
    if(is_name(is)) {
      string s = get_name(is) ;
      if(s != "compute") {
        throw parseError("syntax error, expecting 'compute'") ;
      }
    }
    killsp() ;
  }

  
  if(use_compute && is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;

  bool sized_outputs = false;
  variableSet outsmi = outs ;
  outsmi -= input ;
  for(vi=outsmi.begin();vi!=outsmi.end();++vi) {
    string ot = local_type_map[*vi].first ;
    if(ot == "storeVec" || ot == "storeMat" || ot == "multiStore")
      sized_outputs = true ;
  }

    
    

  if(rule_type == "singleton" ||
     rule_type == "optional"  ||
     rule_type == "default" ||
     rule_type == "constraint" ||
     (output_param && rule_type != "apply" ) ) {
    if(use_prelude) 
      throw parseError("inappropriate prelude") ;
    process_Compute(outputFile,vnames) ;
  } else {
    if(use_compute)
      process_Calculate(outputFile,vnames) ;

    outputFile <<   "    void compute(const Loci::sequence &seq) { " << endl ;
    syncFile(outputFile) ;
    if(use_prelude) {
      outputFile <<   "      prelude(seq) ;" << endl ;
      syncFile(outputFile) ;
    }
    if(use_compute) {
      outputFile <<   "      do_loop(seq,this) ;" << endl ;
      syncFile(outputFile) ;
    }
    outputFile <<   "    }" << endl ;
    syncFile(outputFile) ;
  }
  outputFile <<   "} ;" << endl ;
  syncFile(outputFile) ;

  if(!prettyOutput)
    outputFile << "Loci::register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  else
    outputFile << "register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  syncFile(outputFile) ;

  if(!prettyOutput) {
    outputFile << "}" << endl ;
    syncFile(outputFile) ;
  }

  if(!use_prelude && sized_outputs && (rule_type != "apply")) 
    throw parseError("need prelude to size output type!") ;
}

void parseFile::processFile(string file, ostream &outputFile) {
  bool error = false ;
  filename = file ;
  line_no = 1 ;
  
  is.open(file.c_str(),ios::in) ;
  if(is.fail()) {
    list<string>::const_iterator li ;
    for(li=include_dirs.begin();li!=include_dirs.end();++li) {
      string s = *li + "/" + file ;
      is.clear() ;
      is.open(s.c_str(),ios::in) ;
      if(!is.fail())
        break ;
    }
    if(is.fail()) {
      string s = "can't open include file '" ;
      s += file ;
      s += "'" ;
      throw parseError(s) ;
    }
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
            throw parseError("syntax error: unknown key") ;
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
      //      cerr << "remaining line = " << buf << endl ;
      error = true ;
    }

  } while(!is.eof()) ;
  if(error)
    throw parseError("syntax error") ;
}

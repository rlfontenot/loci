#include <Tools/options_list.h>
#include <Tools/parse.h>
#include <Tools/debug.h>

using std::string ;
using std::ostream ;
using std::istream ;
using std::cerr ;
using std::endl ;


using std::make_pair ;

namespace Loci {

  struct parseError {
    string error_string ;
    parseError() {error_string = "parse error"; }
    parseError(const string &s) { error_string = s ; }
  } ;
  
  options_list::options_list(const string &s) {
    int sz = s.length() ;
    restrict_set = true ;
    if(sz == 0) {
      restrict_set = false ;
      return ;
    }
    string option ;
    for(int i=0;i<sz;++i) {
      if(s[i] == ':') {
        if(option.size() != 0)
          set_of_options.insert(option) ;
        option = "" ;
      } else 
        option += s[i] ;
    }
    if(option.size() != 0)
      set_of_options.insert(option) ;
  }

  option_value_type
    options_list::getOptionValueType(const string &option) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end())
        return NOT_ASSIGNED ;
      return (*tmp).second.value_type ;
    }

  option_values
    options_list::getOption(const string &option) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end())
        return option_values() ;
      return (*tmp).second ;
    }  
  
  void options_list::getOption(const string &option, bool &value) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve BOOLEAN type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != BOOLEAN) ;
    if((*tmp).second.value_type == BOOLEAN)
      value = (*tmp).second.boolean_value ;
  }
  
  void options_list::getOptionUnits(const string &option,
                                    const string &units,
                                    double &value) const {
    if(optionExists(option)) {
      if(getOptionValueType(option) == Loci::REAL) {
        getOption(option,value) ;
      } else if(getOptionValueType(option) == Loci::UNIT_VALUE) {
        Loci::UNIT_type Tu ;
        getOption(option,Tu) ;
        if(!Tu.is_compatible(units)) {
          cerr << "wrong type of unit for " << option <<": " << Tu << endl ;
          cerr << "should have units compatible with " << units << endl ;
        } else {
          value = Tu.get_value_in(units) ;
        }
      } else {
        cerr << "incorrect type for "<< option << endl ;
      }
    }

  }

  void options_list::getOption(const string &option, double &value) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve REAL type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != REAL) ;
    if((*tmp).second.value_type == REAL)
      value = (*tmp).second.real_value ;
  }

  void options_list::getOption(const string &option, UNIT_type &uvalue) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve UNIT_VALUE type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != UNIT_VALUE) ;
    if((*tmp).second.value_type == UNIT_VALUE)
      uvalue = (*tmp).second.units_value ;
  }

  void options_list::getOption(const string &option, string &name) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve NAME type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != NAME) ;
    name = (*tmp).second.name ;
  }

  void options_list::getOption(const string &option, string &name,
                               arg_list &value_list) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve FUNCTION type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != FUNCTION) ;
    
    name = (*tmp).second.name ;
    value_list = (*tmp).second.value_list ;
  }

  void options_list::getOption(const string &option, arg_list &value_list) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      cerr << "WARNING:attempt to retrieve LIST type option " << option
           << " failed." << endl ;
      return ;
    }
    warn((*tmp).second.value_type != LIST) ;
    
    value_list = (*tmp).second.value_list ;
  }

  void options_list::setOption(const string &option, bool value) {
    option_map::iterator tmp ;

    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = BOOLEAN ;
    (*tmp).second.boolean_value = value ;
  }
  
  void options_list::setOption(const string &option, double value) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = REAL ;
    (*tmp).second.real_value = value ;
  }

  void options_list::setOption(const string &option, UNIT_type uvalue) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = UNIT_VALUE ;
    (*tmp).second.units_value = uvalue ;
  }

  void options_list::setOption(const string &option, const string &name) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = NAME ;
    (*tmp).second.name = name ;
    arg_list l ;
    (*tmp).second.value_list = l ;
  }

  void options_list::setOption(const string &option, const arg_list &value_list) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = LIST ;
    (*tmp).second.value_list = value_list ;
  }

  void options_list::setOption(const string &option, const string &name,
                               const arg_list &value_list) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = FUNCTION ;
    (*tmp).second.name = name ;
    (*tmp).second.value_list = value_list ;
  }

  bool
    options_list::checkOption(const string &option, const string &name) const
    {
      option_map::const_iterator tmp ;
      if((tmp = options_db.find(option)) == options_db.end()) {
        return false ;
      }
      warn((*tmp).second.value_type != NAME) ;
      return name == (*tmp).second.name ;
    }

  bool options_list::optionExists(const string &option) const {
    option_map::const_iterator tmp ;
    if((tmp = options_db.find(option)) == options_db.end()) {
      return false ;
    }
    return (*tmp).second.value_type != NOT_ASSIGNED ;
  }

  ostream &option_values::Print(ostream &s) const {
    switch(value_type) {
    case BOOLEAN:
      s << boolean_value?"$true":"$false" ;
      break ;
    case REAL:
      s << real_value ;
      break ;
    case UNIT_VALUE:
      s << units_value ;
      break ;
    case NAME:
      s << name ;
      break ;
    case NAME_ASSIGN:
      s << name << '=' ;
      {
        value_list_type::const_iterator i=value_list.begin() ;
        i->Print(s) ;
        ++i ;
        if(i!=value_list.end())
          cerr << "confused setup in NAME_ASSIGN" << endl;
      }
      break ;
    case STRING:
      s << "\"" << name << "\"" ;
      break ;
    case FUNCTION:
      s << name << "(" ;
      for(value_list_type::const_iterator i=value_list.begin();;) {
        s << *i ;
        ++i ;
        if(i==value_list.end())
          break ;
        s << "," ;
      }
      s << ")" ;
      break ;
    case LIST:
      s << "[" ;
      for(value_list_type::const_iterator i=value_list.begin();;) {
        i->Print(s) ;
        ++i ;
        if(i==value_list.end())
          break ;
        s << "," ;
      }
      s << "]" ;
      break ;
    default:
      break ;
    }
          
    return s ;
  }

  istream &option_values::Input(istream &s) {
    value_list_type l ;
    value_list = l ;
    parse::kill_white_space(s) ;
    if(parse::is_real(s)) {
      real_value = parse::get_real(s) ;
      parse::kill_white_space(s) ;
      if(parse::is_name(s)) {
        string units ;
        int ch = EOF ;
        int opens = 0 ;
        do {
          while(!s.eof() &&((ch=s.peek()) != EOF) &&
                (isalnum(ch) || ch=='*' || ch == '/' ))
            units += s.get() ;
          if(ch == '(') {
            units += s.get() ;
            opens++ ;
          }
          if(opens != 0 && ch == ')') {
            units += s.get() ;
            opens-- ;
          }
        } while(opens!=0) ;

        units_value = UNIT_type(UNIT_type::MKS,"general",real_value,units) ;
        value_type = UNIT_VALUE ;
      } else
        value_type = REAL ;
    } else if(parse::is_name(s)) {
      name = parse::get_name(s) ;
      parse::kill_white_space(s) ;
      value_type = NAME ;
      if(s.peek() == '(') {
        value_type = FUNCTION ;
        s.get() ;
        parse::kill_white_space(s) ;
        while(!s.eof() && s.peek()!=')') {
          option_values ov ;
          s >> ov ;
          value_list.push_back(ov) ;
          parse::kill_white_space(s) ;
          if(s.peek() == ',') {
            s.get() ;
            parse::kill_white_space(s) ;
          }
        }
        if(s.peek() == ')')
          s.get() ;
      }
      if(s.peek() == '=') {
        s.get() ;
        parse::kill_white_space(s) ;
        value_type = NAME_ASSIGN ;
        option_values ov ;
        s >> ov ;
        value_list.push_back(ov) ;
        parse::kill_white_space(s) ;
      }
    } else if(s.peek() == '[') {
      s.get() ;
      value_type = LIST ;
      parse::kill_white_space(s) ;
      while(!s.eof() && s.peek()!=']') {
        option_values ov ;
        s >> ov ;
        parse::kill_white_space(s) ;
        if(s.peek() == '=') {
          s.get() ;
          parse::kill_white_space(s) ;
          
          if(ov.value_type != NAME && ov.value_type != STRING) {
            cerr << "improper assignement in list" << endl ;
            ov.value_type = NOT_ASSIGNED ;
          } else
            ov.value_type = NAME_ASSIGN ;
          option_values ov2 ;
          s >> ov2 ;
          if(ov2.value_type == LIST) {
            ov.value_list = ov2.value_list ;
          } else
            ov.value_list.push_back(ov2) ;
        }
        value_list.push_back(ov) ;
        if(s.peek() == ',') {
          s.get() ;
          parse::kill_white_space(s) ;
        }
      }
      if(s.peek() == ']')
        s.get() ;
    } else if(s.peek() == '"') {
      s.get() ;
      value_type = STRING ;
      char ch = s.get() ;
      while(ch != '"') {
        name += ch ;
        ch = s.get() ;
      }
    } else if(s.peek() == '$') {
      bool v = true ;
      if(s.peek() == '$') {
        s.get() ;
        string bvs = parse::get_name(s) ;
        if(bvs == "false")
          v = false ;
        else if(bvs != "true") {
          cerr << "option_list warning:" << endl ;
          cerr << "boolean value can only be \"$true\"  or \"$false\""
               << endl;
        }
      }
      value_type = BOOLEAN ;
      boolean_value = v ;
    } else {
      parseError err("error reading option_values") ;
      throw err ;
    }
    return s ;
  }

  ostream &options_list::Print(ostream &s) const {
    s << "<" ;
    option_map::const_iterator tmp = options_db.begin() ;
    if(tmp != options_db.end())
      s<<(*tmp).first<<"="<<(*tmp).second ;
    for(++tmp;tmp!=options_db.end();++tmp)
      s<<","<<(*tmp).first<<"="<<(*tmp).second ;
    s << ">" ;
    return s ;
  }

  istream &options_list::Input(istream &s) {
    try {
      parse::kill_white_space(s) ;
      if(s.peek() != '<') {
        cerr << "format error in options_list::Input" << endl ;
        cerr << "expected '<', got '" << s.peek() << "'" << endl ;
        return s ;
      }

      s.get() ;
    
      for(;;) {
        parse::kill_white_space(s) ;
        if(!parse::is_name(s)) {
          cerr << "format error while reading option in option_list::Input"
               << endl ;
          return s ;
        }
        string option = parse::get_name(s) ;
        if(set_of_options.find(option) == set_of_options.end()) {
          if(restrict_set) {
            cerr << "Invalid option name " << option << " in read options"
                 << endl ;
          } else {
            set_of_options.insert(option) ;
          }
        }
        parse::kill_white_space(s) ;
        option_values v ;

        if(s.peek() == '=') {
          s.get() ;
          parse::kill_white_space(s) ;

          s >> v ;
        } else {
          v.value_type = BOOLEAN ;
          v.boolean_value = true ;
        }

        option_map::iterator tmp ; 
        if((tmp = options_db.find(option)) == options_db.end()) {
          tmp = options_db.insert(make_pair(option,v)).first ;
        } else
          (*tmp).second = v ;
        parse::kill_white_space(s) ;
        if(s.peek()=='>') {
          s.get() ;
          break ;
        }
        if(s.peek()!=',')
          cerr << "error reading option " << option << endl ;
        else
          s.get() ;
      }
    }
    catch (parseError err) {
      cerr << err.error_string << endl ;
    }
    return s ;
  }

  void options_list::Input(const arg_list &l) {
    options_db.clear() ;

    for(arg_list::const_iterator li=l.begin();li!=l.end();++li) {
      option_values v  ;
      if(li->type_of() == NAME_ASSIGN) {
        arg_list lv ;
        li->get_value(lv) ;
        warn(lv.size() != 1) ;
        v = *lv.begin() ;
      } else if(li->type_of() == NAME) {
        v.value_type = BOOLEAN ;
        v.boolean_value = true ;
      } else {
        cerr << "unable to parse input value " << *li << endl ;
        continue ;
      }
      string option  ;
      li->get_value(option) ;
      if(set_of_options.find(option) == set_of_options.end()) {
        if(restrict_set) {
          continue ;
        } else {
          set_of_options.insert(option) ;
        }
      }

      option_map::iterator tmp ; 
      if((tmp = options_db.find(option)) == options_db.end()) {
        tmp = options_db.insert(make_pair(option,v)).first ;
      } else
        (*tmp).second = v ;
    }
  }
}

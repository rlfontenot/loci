#include <Tools/options_list.h>
#include <Tools/parse.h>
#include <Tools/debug.h>
#include <Tools/stream.h>

using std::make_pair ;

namespace Loci {

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

  void options_list::setOption(const string &option, double value) {
    option_map::iterator tmp ;
    
    if((tmp = options_db.find(option)) == options_db.end()) {
      option_values v ;
      tmp = options_db.insert(make_pair(option,v)).first ;
    }
    (*tmp).second.value_type = REAL ;
    (*tmp).second.real_value = value ;
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
    case REAL:
      s << real_value ;
      break ;
    case NAME:
      s << name ;
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
        s << *i ;
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
    } else if(s.peek() == '[') {
      s.get() ;
      value_type = LIST ;
      parse::kill_white_space(s) ;
      while(!s.eof() && s.peek()!=']') {
        option_values ov ;
        s >> ov ;
        value_list.push_back(ov) ;
        parse::kill_white_space(s) ;
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
    } else {
      cerr << "error reading option_values" << endl ;
      value_type = NOT_ASSIGNED ;
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
      if(s.peek() != '=')
        cerr << "error reading option " << option << endl ;
      else
        s.get() ;
      parse::kill_white_space(s) ;

      option_map::iterator tmp ;

      option_values v ;
      s >> v ;
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
    
    return s ;
  }
}

#ifndef OPTIONS_LIST_H
#define OPTIONS_LIST_H

#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include <map>

namespace Loci {
  class options_list ;

  enum option_value_type {NOT_ASSIGNED,REAL,NAME,FUNCTION,LIST,STRING} ;

  class option_values {
  public:
    typedef std::vector<option_values> value_list_type ;
  private:
    option_value_type value_type ;
    value_list_type value_list ;
    std::basic_string<char> name ;
    double real_value ;
    friend class options_list ;
  public:
    option_values() { value_type = NOT_ASSIGNED ; }

    option_value_type type_of() { return value_type ; }

    void get_value(double &r) { r = real_value ; }
    void get_value(value_list_type &l) { l = value_list ; }
    void get_value(std::basic_string<char> &n) { n = name ; }

    std::ostream & Print(std::ostream &s) const ;
    std::istream & Input(std::istream &s) ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const option_values &ov)
    { return ov.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, option_values &ov)
    { return ov.Input(s) ; }    

  class options_list {
  public:
    typedef option_values::value_list_type arg_list ;
    typedef std::map<std::string,option_values> option_map ;
    typedef std::set<std::string> option_set ;

  private:
    option_set set_of_options ;
    option_map options_db ;

  protected:
    options_list() {} ;


  public:
    options_list(const std::string &s) ;
        
    option_value_type getOptionValueType(const std::string &option) const ;
    option_values getOption(const std::string &option) const ;
    void getOption(const std::string &option, double &value) const ;
    void getOption(const std::string &option, std::string &name) const ;
    void getOption(const std::string &option, arg_list &value_list) const ;
    void getOption(const std::string &option, std::string &name,
                   arg_list &value_list) const ;

    void setOption(const std::string &option, double value) ;
    void setOption(const std::string &option, const std::string &name) ;
    void setOption(const std::string &option, const arg_list &value_list) ;
    void setOption(const std::string &option, const std::string &name,
                   const arg_list &value_list) ;

    bool checkOption(const std::string &option, const std::string &name) const ;
    bool optionExists(const std::string &option) const ;

    std::ostream & Print(std::ostream &s) const ;
    std::istream & Input(std::istream &s) ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const options_list &ol)
    { return ol.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, options_list &ol)
    { return ol.Input(s) ; }
}
#endif

#ifndef OPTIONS_LIST_H
#define OPTIONS_LIST_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/unit_type.h>

#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include <map>

namespace Loci {
  class options_list ;

  enum option_value_type {NOT_ASSIGNED,REAL,NAME,FUNCTION,LIST,STRING,BOOLEAN,
                          UNIT_VALUE,NAME_ASSIGN} ;

  class option_values {
  public:
    typedef std::vector<option_values> value_list_type ;
  private:
    option_value_type value_type ;
    value_list_type value_list ;
    std::string name ;
    double real_value ;
    bool boolean_value ;
    UNIT_type units_value ;
    friend class options_list ;
  public:
    option_values() { value_type = NOT_ASSIGNED ; real_value = 0 ; boolean_value = false ; }

    option_value_type type_of() const { return value_type ; }

    void get_value(bool &b) const { b = boolean_value; }
    void get_value(double &r) const { r = real_value ; }
    void get_value(value_list_type &l) const { l = value_list ; }
    void get_value(std::string &n) const { n = name ; }
    void get_value(UNIT_type &ut) const { ut = units_value ; }
    
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
    typedef std::list<std::string> option_namelist ;

  private:
    option_set set_of_options ;
    bool restrict_set ;
    option_map options_db ;

  protected:

  public:
    options_list(const std::string &s) ;
    options_list() {restrict_set=false; } ;


    option_namelist getOptionNameList() const {
      option_namelist l ;
      for(option_map::const_iterator mi=options_db.begin();
          mi!=options_db.end();
          mi++) {
        l.push_back(mi->first) ;
      }
      return l ;
    }
    option_value_type getOptionValueType(const std::string &option) const ;
    option_values getOption(const std::string &option) const ;
    void getOption(const std::string &option, bool &value) const ;
    void getOption(const std::string &option, double &value) const ;
    void getOption(const std::string &option, UNIT_type &uvalue) const ;
    void getOption(const std::string &option, std::string &name) const ;
    void getOption(const std::string &option, arg_list &value_list) const ;
    void getOption(const std::string &option, std::string &name,
                   arg_list &value_list) const ;
    void getOptionUnits(const std::string &option, const std::string &units,
                        double &value) const ;

    void setOption(const std::string &option, bool value) ;
    void setOption(const std::string &option, double value) ;
    void setOption(const std::string &option, UNIT_type uvalue) ;
    void setOption(const std::string &option, const std::string &name) ;
    void setOption(const std::string &option, const arg_list &value_list) ;
    void setOption(const std::string &option, const std::string &name,
                   const arg_list &value_list) ;


    bool checkOption(const std::string &option, const std::string &name) const ;
    bool optionExists(const std::string &option) const ;

    std::ostream & Print(std::ostream &s) const ;
    std::istream & Input(std::istream &s) ;
    void Input(const arg_list &l) ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const options_list &ol)
    { return ol.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, options_list &ol)
    { return ol.Input(s) ; }
}
#endif

#ifndef UNIT_TYPE_H
#define UNIT_TYPE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/parse.h>
#include <Tools/expr.h>
#include <stack>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

namespace Loci {

  class UNIT_type{

  public:
    std::string unit_kind;//the kind of unit ,such as pressure, time and so on
    std::string input_unit;// the unit which you input
    double input_value;// the value which you input
    enum unit_mode{MKS, CGS, check_available};
    // MKS: use MKS system
    // CGS: use CGS system
    //check_available: check the system if the unit available
    unit_mode mode;

    double value;//temp container of value calculation
    std::map<std::string,int> unit_num_map,unit_den_map;//containers of numerator and dnominator
    enum basic_unit_type {Length,Mass,Time,Temperature,Electric_current,
			  Amount_of_substance,Luminous_intensity, Angle, NoDim};
    double conversion_factor;

    //three tables of unit type - basic, composite, reference types----//
    struct basic_units{
      const char* name;
      basic_unit_type unit_type;
      double convert_factor;};  
    static basic_units basic_unit_table[] ;//MKS
    static basic_units cgs_basic_unit_table[];//CGS
   
    struct composite_units{
      const char* name;
      const char* derived_unit;
      double convert_factor;};  
    static composite_units composite_unit_table[] ;//MKS
    static composite_units cgs_composite_unit_table[] ;//CGS
  
    struct reference_units{
      const char* name;
      const char* refer_unit;
      double convert_factor;};  
    static reference_units reference_unit_table[] ;

    struct default_units{//If no unit input, check the default table 
      const char* default_type;//according to the unit_kind
      const char* default_name;
    };
    static default_units default_unit_table[] ;
  
  private:
    //check there is single temperature or temperature internal
    int is_single_temperature(const exprP input_expr);
    // if single temperature, do the conversion
    void calculate_temperature(exprP &input_expr, double &value);
    // when you need to convert to other temperature, reserve do it.
    void reverse_calculate_temperature(exprP &input_expr,double &value);

  public:
    exprP input(std::istream &in);// get the input unit
    void output(exprP &in_exp);// ouput the converted basic unit
    bool is_in_db(const std::string &str);// check is the unit available in db
    bool is_compatible(const std::string unit_str);//check two units compatible
    //private_is_compatible//
    // check with unit_kind for the input unit, for example: unit_kind=pressure    //input_unit=second, then it is not compatible. But you can set the 
    //unit_kind to "general" then it will not check again. Or you may not need 
    //this function for checking
    bool private_is_compatible();
    //get the value in converted unit
    double get_value_in(const std::string unit_str);

    UNIT_type(unit_mode in_mode, std::string in_kind, double in_value, std::string in_unit) {mode=in_mode,unit_kind=in_kind,value=in_value,input_unit=in_unit;
    input_value=value;
    exprP exp;
    exp=expression::create(input_unit);
    output(exp);
    }

    //    UNIT_type(const UNIT_type &ut) {unit_num_map=ut.unit_num_map,unit_den_map=ut.unit_den_map,value=ut.value,input_value=ut.input_value,conversion_factor=ut.conversion_factor;}
    UNIT_type() { mode=MKS, unit_kind="Pressure";}

  private:
    bool is_reference_unit(std::string str);
    int where_reference_unit(std::string str);
    bool is_composite_unit(std::string str);
    int where_composite_unit(std::string str);
    bool is_basic_unit(std::string str);
    int where_basic_unit(std::string str);

    void printing(std::ostream &s,exprList exlist);//print undecomposed unit

    void build_lists(exprList &numerator, exprList &denominator, exprP input, bool isnum=true);
    void count_dim(std::map<std::string,int> &c_map, const exprList c_list);
    void rem_dup(std::map<std::string,int> &num_map,std::map<std::string,int> &den_map);
    std::map<std::string,int> combine_units(std::map<std::string,int> num_map,std::map<std::string,int> den_map);
    void seperate_unit(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,exprP p);

    void change_to_basic_unit(std::map<std::string,int>initial_map,std::map<std::string,int>&num_map,std::map<std::string,int>&den_map,double &conversion_factor);
    void get_conversion(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,double &conversion_factor);
    bool check_unit(std::istream &in, double &value);//check input unit and get the value
    exprP set_default_unit(exprP &in_exp);
    int in_unit_kind();// check the unit_kind is available in db

  };

  inline std::ostream &operator<<(std::ostream &s, const UNIT_type &o_unit){
    std::map<std::string,int>::iterator mi,mj;
    std::map<std::string,int> n_map=o_unit.unit_num_map,d_map=o_unit.unit_den_map;
    /*cout<<"Numerator: "<<std::endl;
      for(mi= n_map.begin();mi!=n_map.end();++mi)
	cout<<mi->first<<"  "<<mi->second<<std::endl;
      cout<<"Denominator: "<<std::endl;
      for(mj= d_map.begin();mj!=d_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<std::endl;*/

    s<<std::endl;
    s<<o_unit.conversion_factor*o_unit.input_value;
    s<<" ";
    for(mi= n_map.begin();mi!=n_map.end();++mi){
      if(mi->second==1)
	s<<mi->first;
      else
	s<<"("<<mi->first<<"^"<<mi->second<<")";
      if(++mi!=n_map.end())
	s<<"*";
      --mi;
    }
    mj= d_map.begin();
    if(mj!=d_map.end()){
      s<<"/"<<"(";
      for(;mj!=d_map.end();++mj){
	if(mj->second==1)
	  s<<mj->first;
	else
	  s<<"("<<mj->first<<"^"<<mj->second<<")";   
      if(++mj!=d_map.end())
	s<<"*";
      --mj;
      }
      s<<")";
    }
    return s;
    }

  inline std::istream &operator>>(std::istream &s, UNIT_type &unit){
    exprP exp1;
    exp1=unit.input(s);//get the expression for unit
    unit.output(exp1); // get the unit exressed by basic units
    return s;
  }

  inline void unit_error(const int e_no, const std::string &err)
    {
      std::cerr << "error:" << err << std::endl;
    }
}
#endif

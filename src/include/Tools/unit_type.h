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
    double value;
    std::map<std::string,int> unit_num_map,unit_den_map;
    enum basic_unit_type {Length,Mass,Time,Temperature,Electric_current,
			  Amount_of_substance,Luminous_intensity};  
    double conversion_factor;

    //three tables of unit type - basic, composite, reference types----//
    struct basic_units{
      char* name;
      basic_unit_type unit_type;
      double convert_factor;};  
    static basic_units basic_unit_table[] ;
    
    struct composite_units{
      char* name;
      std::string derived_unit;
      double convert_factor;};  
    static composite_units composite_unit_table[] ;
  
    struct reference_units{
      char* name;
      std::string refer_unit;
      double convert_factor;};  
    static reference_units reference_unit_table[] ;
  
  private:
    int is_single_temperature(const exprP input_expr);
    void calculate_temperature(exprP &input_expr, double &value);

  public:
    exprP input(std::istream &in);
    void output(exprP &in_exp);
    double get_value();

  private:
    bool is_reference_unit(std::string str);
    int where_reference_unit(std::string str);
    bool is_composite_unit(std::string str);
    int where_composite_unit(std::string str);
    bool is_basic_unit(std::string str);
    int where_basic_unit(std::string str);

    void printing(std::ostream &s,exprList exlist);

    void build_lists(exprList &numerator, exprList &denominator, exprP input, bool isnum=true);
    void count_dim(std::map<std::string,int> &c_map, const exprList c_list);
    void rem_dup(std::map<std::string,int> &num_map,std::map<std::string,int> &den_map);
    std::map<std::string,int> combine_units(std::map<std::string,int> num_map,std::map<std::string,int> den_map);
    void seperate_unit(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,exprP p);

    void change_to_basic_unit(std::map<std::string,int>initial_map,std::map<std::string,int>&num_map,std::map<std::string,int>&den_map,double &conversion_factor);
    void get_conversion(std::map<std::string,int> &num_map, std::map<std::string,int> &den_map,double &conversion_factor);
    bool check_unit(istream &in, double &value);

  };

  inline std::ostream &operator<<(std::ostream &s, UNIT_type &unit){
    std::map<std::string,int>::iterator mi,mj;
    std::map<std::string,int> n_map=unit.unit_num_map,d_map=unit.unit_den_map;
      cout<<"Numerator: "<<endl;
      for(mi= n_map.begin();mi!=n_map.end();++mi)
	cout<<mi->first<<"  "<<mi->second<<endl;
      cout<<"Denominator: "<<endl;
      for(mj= d_map.begin();mj!=d_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<endl;

      cout<<"conversion factor is: "<<unit.conversion_factor<<endl;
      cout<<"input value is: "<<unit.value<<endl;


    return s;
    }

  inline std::istream &operator>>(std::istream &s, UNIT_type &unit){
    exprP exp1;
    exp1=unit.input(s);//get the expression for unit
    unit.output(exp1); // get the unit exressed by basic units
    return s;
  }
}

//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef UNIT_TYPE_H
#define UNIT_TYPE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/except.h>

#include <Tools/parse.h>
#include <Tools/expr.h>
#include <Tools/autodiff.h>

#include <stack>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

namespace Loci
{

class UNIT_type
{
  /** **************************************************************************
   * @brief
   * @param str
   * @return true
   * @return false
   ****************************************************************************/
   bool is_reference_unit(std::string str);

  /** **************************************************************************
   * @brief
   * @param str
   * @return true
   * @return false
   ****************************************************************************/
  bool is_composite_unit(std::string str);

  /** **************************************************************************
   * @brief
   * @param str
   * @return true
   * @return false
   ****************************************************************************/
  bool is_basic_unit(std::string str);

  /** **************************************************************************
   * @brief  check input unit and get the value
   * @param  in
   * @param  value
   * @return true
   * @return false
   ****************************************************************************/
  bool check_unit(std::istream &in, double &value);

  /** **************************************************************************
   * @brief check input unit and get the value
   * @param in
   * @param value
   * @return true
   * @return false
   ****************************************************************************/
  bool check_unit(std::istream &in, MFADd &value);

  /** **************************************************************************
   * @brief check input unit and get the value
   * @param in
   * @param value
   * @return true
   * @return false
   ****************************************************************************/
  bool check_unit(std::istream &in, FAD2d &value);

  /** **************************************************************************
   * @brief
   * @param in
   * @param value
   * @return true
   * @return false
   ****************************************************************************/
  bool check_unit(std::istream &in, FADd &value)
  {
    FAD2d v;
    bool  ret = check_unit(in, v);
    value     = FADd(v.value, v.grad);
    return ret;
  }

  /** **************************************************************************
   * @brief
   * @param str
   * @return int
   ****************************************************************************/
  int where_reference_unit(std::string str);

  /** **************************************************************************
   * @brief
   * @param str
   * @return int
   ****************************************************************************/
  int where_composite_unit(std::string str);

  /** **************************************************************************
   * @brief
   * @param str
   * @return int
   ****************************************************************************/
  int where_basic_unit(std::string str);

  /** **************************************************************************
   * @brief check there is single temperature or temperature internal
   * @param input_expr
   * @return int
   ****************************************************************************/
  int is_single_temperature(const exprP input_expr);

  /** **************************************************************************
   * @brief  check the unit_kind is available in db
   * @return int
   ****************************************************************************/
  int in_unit_kind();

  /** **************************************************************************
   * @brief Set the default unit object
   * @param in_exp
   * @return exprP
   ****************************************************************************/
  exprP set_default_unit(exprP &in_exp);

  /** **************************************************************************
   * @brief print undecomposed unit
   * @param s
   * @param exlist
   ****************************************************************************/
  void printing(std::ostream &s, exprList exlist);

  /** **************************************************************************
   * @brief
   * @param numerator
   * @param denominator
   * @param input
   * @param isnum
   ****************************************************************************/
  void build_lists(exprList &numerator, exprList &denominator, exprP input,
                   bool isnum=true);

  /** **************************************************************************
   * @brief
   * @param c_map
   * @param c_list
   ****************************************************************************/
  void count_dim(std::map<std::string, int> &c_map, const exprList c_list);

  /** **************************************************************************
   * @brief
   * @param num_map
   * @param den_map
   ****************************************************************************/
  void rem_dup(std::map<std::string, int> &num_map,
               std::map<std::string, int> &den_map);

  /** **************************************************************************
   * @brief
   * @param num_map
   * @param den_map
   * @param p
   ****************************************************************************/
  void seperate_unit(std::map<std::string, int> &num_map,
                     std::map<std::string, int> &den_map, exprP p);

  /** **************************************************************************
   * @brief
   * @param initial_map
   * @param num_map
   * @param den_map
   * @param conversion_factor
   ****************************************************************************/
  void change_to_basic_unit(std::map<std::string, int> initial_map,
                            std::map<std::string, int> &num_map,
                            std::map<std::string, int> &den_map,
                            double                     &conversion_factor);

  /** **************************************************************************
   * @brief Get the conversion object
   * @param num_map
   * @param den_map
   * @param conversion_factor
   ****************************************************************************/
  void get_conversion(std::map<std::string, int> &num_map,
                      std::map<std::string, int> &den_map,
                      double                     &conversion_factor);

  /** **************************************************************************
   * @brief if single temperature, do the conversion
   * @param input_expr
   * @param value
   ****************************************************************************/
  void calculate_temperature(exprP &input_expr, MFADd &value);

  /** **************************************************************************
   * @brief
   * @param input_expr
   * @param value
   ****************************************************************************/
  void calculate_temperature(exprP &input_expr, FAD2d &value);

  /** **************************************************************************
   * @brief
   * @param input_expr
   * @param value
   ****************************************************************************/
  void calculate_temperature(exprP &input_expr, FADd &value)
  {
    FAD2d v;
    calculate_temperature(input_expr, v);
    value = FADd(v.value, v.grad);
  }

  /** **************************************************************************
   * @brief when you need to convert to other temperature, reserve do it.
   * @param input_expr
   * @param value
   ****************************************************************************/
  void reverse_calculate_temperature(exprP &input_expr, MFADd &value);

  /** **************************************************************************
   * @brief
   * @param input_expr
   * @param value
   ****************************************************************************/
  void reverse_calculate_temperature(exprP &input_expr, FAD2d &value);

  /** **************************************************************************
   * @brief
   * @param input_expr
   * @param value
   ****************************************************************************/
  void reverse_calculate_temperature(exprP &input_expr, FADd &value)
  {
    FAD2d v;
    reverse_calculate_temperature(input_expr, v);
    value = FADd(v.value, v.grad);
  }

  /** **************************************************************************
   * @brief
   * @param num_map
   * @param den_map
   * @return std::map<std::string, int>
   ****************************************************************************/
  std::map<std::string, int> combine_units(std::map<std::string, int> num_map,
                                           std::map<std::string, int> den_map);

public:

  /// @brief
  enum eUnit_Mode
  {
    MKS,            ///< MKS: use MKS system
    CGS,            ///< CGS: use CGS system
    check_available ///< check the system if the unit available
  };

  /// @brief
  enum eBasic_Unit_Type
  {
    Length,              ///<
    Mass,                ///<
    Time,                ///<
    Temperature,         ///<
    Electric_current,    ///<
    Amount_of_substance, ///<
    Luminous_intensity,  ///<
    Angle,               ///<
    NoDim                ///<
  };

  /// @brief Tables of unit type basic
  struct basic_units
  {
    const char*      name;           ///<
    eBasic_Unit_Type unit_type;      ///<
    double           convert_factor; ///<
  };

  /// @brief Tables of unit type composite
  struct composite_units
  {
    const char* name;           ///<
    const char* derived_unit;   ///<
    double      convert_factor; ///<
  };

  /// @brief Tables of unit type reference
  struct reference_units
  {
    const char* name;           ///<
    const char* refer_unit;     ///<
    double      convert_factor; ///<
  };

  /// @brief If no unit input, check the default table according to the unit_kind
  struct default_units
  {
    const char* default_type; ///<
    const char* default_name; ///<
  };

  double      conversion_factor; ///<
  std::string unit_kind;         ///< The kind of unit , such as pressure, time and so on
  std::string input_unit;        ///< The unit which you input
  FAD2d       input_value;       ///< The value which you input
  MFADd       input_value_mfad;  ///< The value which you input
  eUnit_Mode  mode;              ///<
  FAD2d       value;             ///< temp container of value calculation
  MFADd       value_mfad;        ///< temp container of value calculation

  std::map<std::string, int> unit_num_map; ///< containers of numerator
  std::map<std::string, int> unit_den_map; ///< containers of denominator

  static basic_units     basic_unit_table[];         ///< MKS
  static basic_units     cgs_basic_unit_table[];     ///< CGS
  static composite_units composite_unit_table[];     ///< MKS
  static composite_units cgs_composite_unit_table[]; ///< CGS
  static reference_units reference_unit_table[];     ///<
  static default_units   default_unit_table[];       ///<

  /** **************************************************************************
   * @brief Construct a new unit type object
   ****************************************************************************/
  UNIT_type()
  {
    mode              = MKS;
    unit_kind         = "";
    value             = 0;
    conversion_factor = 1;
    input_value       = 0;
    input_value_mfad  = 0;
  }

  /** **************************************************************************
   * @brief Construct a new unit type object
   * @param in_mode
   * @param in_kind
   * @param in_value
   * @param in_unit
   ****************************************************************************/
  UNIT_type(eUnit_Mode in_mode, std::string in_kind, FADd in_value, std::string in_unit)
  {
    mode                   = in_mode;
    unit_kind              = in_kind;
    value                  = FAD2d(in_value.value, in_value.grad, 0.0);
    input_unit             = in_unit;
    input_value            = value;
    input_value_mfad.value = input_value.value;

    for(size_t i=0; i<input_value_mfad.maxN; ++i)
      input_value_mfad.grad[i] = 0;

    exprP exp;
    exp = expression::create(input_unit);
    output(exp);
  }

  /** **************************************************************************
   * @brief Construct a new unit type object
   * @param in_mode
   * @param in_kind
   * @param in_value
   * @param in_unit
   ****************************************************************************/
  UNIT_type(eUnit_Mode in_mode, std::string in_kind, FAD2d in_value, std::string in_unit)
  {
    mode                   = in_mode;
    unit_kind              = in_kind;
    value                  = in_value;
    input_unit             = in_unit;
    input_value            = value;
    input_value_mfad.value = input_value.value;

    for(size_t i=0; i<input_value_mfad.maxN; ++i)
      input_value_mfad.grad[i] = 0;

    exprP exp;
    exp = expression::create(input_unit);
    output(exp);
  }

  /** **************************************************************************
   * @brief Construct a new unit type object
   * @param in_mode
   * @param in_kind
   * @param in_value
   * @param in_unit
   ****************************************************************************/
  UNIT_type(eUnit_Mode in_mode, std::string in_kind, MFADd in_value, std::string in_unit)
  {
    mode             = in_mode;
    unit_kind        = in_kind;
    value_mfad       = in_value;
    input_unit       = in_unit;
    input_value      = value_mfad.value;
    input_value_mfad = value_mfad;
    exprP exp;
    exp = expression::create(input_unit);
    output(exp);
  }

  /** **************************************************************************
   * @brief Construct a new unit type object
   * @param in_mode
   * @param in_kind
   * @param in_value
   * @param in_unit
   ****************************************************************************/
  UNIT_type(eUnit_Mode in_mode, std::string in_kind, double in_value, std::string in_unit)
  {
    mode                   = in_mode;
    unit_kind              = in_kind;
    value                  = in_value;
    input_unit             = in_unit;
    input_value            = value;
    input_value_mfad.value = value.value;

    for(size_t i=0; i<input_value_mfad.maxN; ++i)
      input_value_mfad.grad[i] = 0;

    exprP exp;
    exp = expression::create(input_unit);
    output(exp);
  }

  /** **************************************************************************
   * @brief  get the input unit
   * @param in
   * @return exprP
   ****************************************************************************/
  exprP input(std::istream &in);

  /** **************************************************************************
   * @brief  ouput the converted basic unit
   * @param in_exp
   ****************************************************************************/
  void output(exprP &in_exp);

  /** **************************************************************************
   * @brief check is the unit available in db
   * @param str
   * @return true
   * @return false
   ****************************************************************************/
  bool is_in_db(const std::string &str);

  /** **************************************************************************
   * @brief check two units compatible
   * @param unit_str
   * @return true
   * @return false
   ****************************************************************************/
  bool is_compatible(const std::string unit_str);

  /** **************************************************************************
   * @brief Check with unit_kind for the input unit, for example:
   * `unit_kind=pressure` and  `input_unit=second`, then it is not compatible.
   * But you can set the `unit_kind=general` then it will not check again. Or
   * you may not need this function for checking.
   ****************************************************************************/
  bool private_is_compatible();

  /** **************************************************************************
   * @brief get the value in converted unit
   * @param unit_str
   * @return double
   ****************************************************************************/
  double get_value_in(const std::string unit_str);

  /** **************************************************************************
   * @brief Get the value inD object
   * @param unit_str
   * @return FAD2d
   ****************************************************************************/
  FAD2d get_value_inD(const std::string unit_str);

  /** **************************************************************************
   * @brief Get the value inM object
   * @param unit_str
   * @return MFADd
   ****************************************************************************/
  MFADd get_value_inM(const std::string unit_str);

}; // End Class UNIT_type

/** ****************************************************************************
 * @brief
 * @param s
 * @param o_unit
 * @return std::ostream&
 ******************************************************************************/
inline std::ostream &operator<<(std::ostream &s, const UNIT_type &o_unit)
{
  s << o_unit.input_value << ' ' << o_unit.input_unit;
  return s;
}

/** ****************************************************************************
 * @brief
 * @param s
 * @param unit
 * @return std::istream&
 ******************************************************************************/
inline std::istream &operator>>(std::istream &s, UNIT_type &unit)
{
  exprP exp1;
  exp1 = unit.input(s); // get the expression for unit
  unit.output(exp1);    // get the unit exressed by basic units
  return s;
}

/** ****************************************************************************
 * @brief
 * @param e_no
 * @param err
 ******************************************************************************/
inline void unit_error(const int e_no, const std::string &err)
{
  throw StringError(err);
}
} // End Namespace Loci
#endif

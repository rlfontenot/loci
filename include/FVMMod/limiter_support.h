/** ****************************************************************************
 * @file      limiter_support.h
 * @authors   Ed Luke (MSU)
 *            Raymond Fontenot (CFDRC)
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Limiter support classes and limiter function defintions.
 * @details   This file is a part of the Loci Framework, a free software.
 * You can redistribute it and/or modify it under the terms of the Lesser
 * GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * The Loci Framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
 ******************************************************************************/
#ifndef _LIMITER_SUPPORT_H
#define _LIMITER_SUPPORT_H
#include <Loci.h>

using std::cerr;
using std::endl;

namespace Loci
{
typedef real_t           real;
typedef vector3d<real_t> vect3d;

/// @brief Class for storing max/min/norm for scalars
class stenMaxMinNorm
{
public:
  real max;   //!< [-] max value for stenMaxMinNorm
  real min;   //!< [-] min value for stenMaxMinNorm
  real norm;  //!< [-] norm for stenMaxMinNorm
};

/// @brief Class for storing max/min/norm for vect3d
class stenMaxMinNormv3d
{
public:
  vect3d max;   //!< [-] max value for stenMaxMinNormv3d
  vect3d min;   //!< [-] min value for stenMaxMinNormv3d
  real   norm;  //!< [-] norm for stenMaxMinNormv3d
};


/** ****************************************************************************
 * @brief     Template operation for Max and Min of a scalar
 * @tparam T  stenMaxMinNorm
 ******************************************************************************/
template <class T> struct MaxMin
{
  /**
   * @brief
   * @param r in/out T
   * @param s in T for comparison
   */
  void operator()(T &r, const T &s)
  {
    r.max = max(r.max,s.max);
    r.min = min(r.min,s.min);
  }
};

/** ****************************************************************************
 * @brief  Template operation for summing Max and Min of a scalar
 * @tparam T  stenMaxMinNorm
 ******************************************************************************/
template <class T> struct SumMaxMin
{
  /**
   * @brief
   * @param r    in/out T
   * @param s    in T for comparison
   */
  void operator()(T &r, const T &s)
  {
    r.max += s.max;
    r.min += s.min;
  }
};


/** ****************************************************************************
 * @brief Template operation for Max and Min of a vect3d
 * @tparam T   stenMaxMinNormv3d
 ******************************************************************************/
template <class T> struct MaxMinV3D
{
  /**
   * @brief
   * @param r in/out T
   * @param s in T for comparison
   */
  void operator()(T &r, const T &s)
  {
    r.max.x = max(r.max.x,s.max.x);
    r.max.y = max(r.max.y,s.max.y);
    r.max.z = max(r.max.z,s.max.z);
    r.min.x = min(r.min.x,s.min.x);
    r.min.y = min(r.min.y,s.min.y);
    r.min.z = min(r.min.z,s.min.z);
  }
};

/** ****************************************************************************
 * @brief Template operation for summing Max and Min of a vector
 * @tparam T stenMaxMinNorm
 ******************************************************************************/
template <class T> struct SumMaxMinV3D
{
  /**
   * @brief
   * @param r in/out T
   * @param s in T for comparison
   */
  void operator()(T &r, const T &s)
  {
    r.max.x += s.max.x;
    r.max.y += s.max.y;
    r.max.z += s.max.z;
    r.min.x += s.min.x;
    r.min.y += s.min.y;
    r.min.z += s.min.z;
  }
};

/** ****************************************************************************
 * @brief Template operation for Max and Min of a vector
 * @tparam T stenMaxMinNorm
 ******************************************************************************/
template <class T> struct MaxMinv
{
  /**
   * @brief do nothing, as max/min for this is handled explicitly
   * @param r in/out T
   * @param s in T for comparison
   */
  void operator()(T &r, const T &s)
  {}
};


/** ****************************************************************************
 * @brief Template operation for summing Max and Min of a vect3d
 * @tparam T stenMaxMinNormv3d
 ******************************************************************************/
template <class T> struct SumMaxMinv
{
  /**
   * @brief do nothing, as max/min for this is handled explicitly
   * @param r in/out T
   * @param s in T for comparison
   */
  void operator()(T &r, const T &s)
  {}
};

/** ****************************************************************************
 * @brief Venkatakrishnan Limiter Function
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param eps2     [-] cell volume function
 * @return real    [-] limiter
 ******************************************************************************/
real vlimit(real Xcc, real qmin, real qmax, real qdif, real eps2);


/** ****************************************************************************
 * @brief Barth Limiter function
 * @param Xcc       [-] Cell center value
 * @param qdif      [-] face - cell center value
 * @param qmax      [-] maximum value in local stencil
 * @param qmin      [-] minimum value in local stencil
 * @return real     [-] limiter value
 ******************************************************************************/
real barth_limit(real Xcc, real qdif, real qmax, real qmin);

/** ****************************************************************************
 * @brief Nishikawa limiter function. See () for reference
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param epsp     [-] cell volume function
 * @param nisPow   [-] order of the limiter
 * @return real    [-] limiter
 ******************************************************************************/
real nis_limit(real Xcc, real qmin, real qmax, real qdif, real epsp, real nisPow);


/** ****************************************************************************
 * @brief  Register for stenMaxMinNorm class
 * @tparam stenMaxMinNorm
 ******************************************************************************/
template<> struct data_schema_traits<stenMaxMinNorm>
{
  typedef IDENTITY_CONVERTER Schema_Converter;
  static DatatypeP get_type()
  {
    CompoundDatatypeP ct = CompoundFactory(stenMaxMinNorm());
    LOCI_INSERT_TYPE(ct,stenMaxMinNorm,max);
    LOCI_INSERT_TYPE(ct,stenMaxMinNorm,min);
    LOCI_INSERT_TYPE(ct,stenMaxMinNorm,norm);
    return DatatypeP(ct);
  }
};

/** ****************************************************************************
 * @brief Register for stenMaxMinNormv3d class
 * @tparam stenMaxMinNormv3d
 ******************************************************************************/
template<> struct data_schema_traits<stenMaxMinNormv3d>
{
  typedef IDENTITY_CONVERTER Schema_Converter;
  static DatatypeP get_type()
  {
    CompoundDatatypeP ct = CompoundFactory(stenMaxMinNormv3d());
    LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,max);
    LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,min);
    LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,norm);
    return DatatypeP(ct);
  }
};
} // End Namespace Loci

#endif
/** ****************************************************************************
 * @file      limiter_support.cpp
 * @authors   Ed Luke (MS State)
 *            Raymond Fontenot (CFDRC)
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Support rules for limiter rules
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
#include "FVMMod/limiter_support.h"

namespace Loci
{

/** ****************************************************************************
 * @brief Venkatakrishnan Limiter Function,
 * written in standard form, see @cite nishikawa2022new for reference
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param eps2     [-] cell volume function
 * @return real    [-] limiter
 ******************************************************************************/
real vlimit(real Xcc, real qmin, real qmax, real qdif, real eps2)
{
  const real delp = (qdif>=0.0) ? qmax-Xcc   : qmin-Xcc;             // delta +
  const real delm = (qdif > 0)  ? qdif+1e-30 : qdif-1e-30;           // delta -
  const real num  = ((delp*delp+eps2)*delm+ 2.0*delm*delm*delp);     // numerator of limiter
  const real den  = (delm*(delp*delp+2.0*delm*delm+delm*delp+eps2)); // denominator of limiter
  real       lim  = num/den;   // make the limiting case of 0/0 work as expected
  return max(0.0,min(1.0,lim));
}

/** ****************************************************************************
 * @brief Barth Limiter function
 * @param Xcc       [-] Cell center value
 * @param qdif      [-] face - cell center value
 * @param qmax      [-] maximum value in local stencil
 * @param qmin      [-] minimum value in local stencil
 * @return real     [-] limiter value
 ******************************************************************************/
real barth_limit(real Xcc, real qdif, real qmax, real qmin)
{
  real lim;

  if(qdif > 0) lim = (qmax-Xcc)/(qdif+1e-30);
  else         lim = (qmin-Xcc)/(qdif-1e-30);

  return max(0.0,min(1.0,lim));
}

/** ****************************************************************************
 * @brief Nishikawa limiter function. See @cite nishikawa2022new for reference
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param epsp     [-] cell volume function
 * @param nisPow   [-] order for Nishikawa Limiter
 * @return real    [-] limiter
 ******************************************************************************/
real nis_limit(real Xcc, real qmin, real qmax, real qdif, real epsp, real nisPow)
{
  const real delp = (qdif>=0.0) ? qmax-Xcc : qmin-Xcc;  // delta +
  const real delm = qdif;                               // delta -
  const real a    = abs(delp);
  const real b    = abs(delm) + 1e-30;

  // save on computation
  if(a > 2.*b)
    return 1.0;

  const real fun1 = pow(a,nisPow) + epsp; // numerator of limiter
  real       Sp;

  if(nisPow == 3)      Sp = 4.0*b*b;
  else if(nisPow == 4) Sp = 2.0*b*(  a*a - 2.0*b*(a-b));
  else if(nisPow == 5) Sp = 8.0*b*b*(a*a - 2.0*b*(a-b));

  const real num = fun1 + a*Sp;
  const real den = fun1 + b*(pow(delp,nisPow-1.0) + Sp);  // denominator of limiter
  real       lim = num/den;
  return max(0.0,lim);
}
}
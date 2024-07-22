/** ****************************************************************************
 * @file      slope_rescale.cpp
 * @author    Ed Luke (MS State)
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Scale vectors, function definition
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

#include "FVMMod/slope_rescale.h"

namespace Loci {
/** ****************************************************************************
 * @brief Void function to rescale slope in muscl reconstruction for vectors
 * @param[out] slopes  real, rescaled slope for reconstuction
 * @param[in]  vs      int, size of vector
 ******************************************************************************/
void slope_rescale(real *slopes, int vs)
{
  real sum  = 0.0;
  real sumn = 0.0;
  real sump = 0.0;
  for(int i=0;i<vs;++i)
  {
    sum  += slopes[i];
    sumn += min(real(0.0),slopes[i]);
    sump += max(real(0.0),slopes[i]);
  }
  if(sum < 0.0)
  {
    for(int i=0;i<vs;++i)
    {
      slopes[i] = slopes[i]*((slopes[i]<0.0) ? (sumn-sum)/(sumn-(real)1e-30) : (real)1.0);
    }
    if(sum > 0.0)
    {
      for(int i=0;i<vs;++i)
      {
        slopes[i] = slopes[i]*((slopes[i]>0.0) ? (sump-sum)/(sump-(real)1e-30) : (real)1.0);
      }
    }
  }
}
}
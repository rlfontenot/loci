/** ****************************************************************************
 * @file      slope_rescale.h
 * @author    Ed Luke (MSU)
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Sets frozen/minimum limiter for vect3d state variables
 * @details   The software tool Loci/GGFS module consisting is being furnished
 *            to NASA Personnel under SBIR Data Rights. The SBIR DATA rights are
 *            asserted by CFD Research Corporation.
 *
 *            These SBIR data are furnished with SBIR rights under Contract No.
 *            80NSSC18P2154. For a period of 4 years, unless extended in
 *            accordance with FAR 27.409(h), after acceptance of all items to be
 *            delivered under this contract, the Government will use these data
 *            for Government purposes only, and they shall not be disclosed
 *            outside the Government (including disclosure for procurement
 *            purposes) during such period without permission of the Contractor,
 *            except that, subject to the foregoing use and disclosure
 *            prohibitions, these data may be disclosed for use by support
 *            Contractors. After the protection period, the Government has a
 *            paid-up license to use, and to authorize others to use on its
 *            behalf, these data for Government purposes, but is relieved of all
 *            disclosure prohibitions and assumes no liability for unauthorized
 *            use of these data by third parties. This notice shall be affixed
 *            to any reproductions of these data, in whole or in part.
 * @attention Distribution C: Limited to Government Employees only. A release
 *            under Distribution B and A is being considered and may be done for
 *            future releases of the code.
 ******************************************************************************/
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <Loci.h>
#include "FVMAdapt/sciTypes.h"

#ifndef _SLOPE_RESCALE_H
#define _SLOPE_RESCALE_H
namespace Loci {

typedef real_t real;
/** ****************************************************************************
 * @brief Function to rescale slope in muscl reconstruction for vectors
 * @param[out] slopes  real, rescaled slope for reconstuction
 * @param[in]  vs      int, size of vector
 ******************************************************************************/
void slope_rescale(real *slopes, int vs);
}

#endif
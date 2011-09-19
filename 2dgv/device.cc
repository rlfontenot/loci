#include "device.h"

void grDevice::DrawPoly(const positions *pnts, int num_pnts) {
  const positions *p = pnts ;
  for(int i=0;i<num_pnts-1;++i) {
    const positions *pold = p ;
    p++ ;
    DrawLine(*pold,*p) ;
  }
}

   

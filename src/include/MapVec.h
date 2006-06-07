#ifndef MAPVEC_H
#define MAPVEC_H

#include <DMapVec.h>
#include <MapVec_def.h>
#include <MapVec_impl.h>

namespace Loci {
  template <int M>
  inline void remapMap(MapVec<M> &m, const Map &remap, Map &tmp) {
    for(int i=0;i<M;++i) {
      FORALL(m.domain(),fc) {
        tmp[remap[fc]] = m[fc][i] ;
      } ENDFORALL ;
      FORALL(m.domain(),fc) {
        m[fc][i] = tmp[fc] ;
      } ENDFORALL ;
    }
  }
}
    

#endif

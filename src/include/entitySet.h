#ifndef ENTITYSET_H
#define ENTITYSET_H 1

#include <Tools/intervalSet.h>

namespace Loci {
  typedef intervalSet entitySet ;

  template<class T> inline entitySet create_entitySet(T start,T end) {
    return create_intervalSet(start,end) ;
  }

}

#endif

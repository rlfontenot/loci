#ifndef MAP_REP_H
#define MAP_REP_H

#include "store_rep.h"

namespace Loci {

  class multiMap ;
  
  class MapRep : public storeRep {
  public:
    virtual ~MapRep() ;
    virtual entitySet image(const entitySet &domain) const = 0 ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const = 0 ;
    virtual multiMap get_map() = 0 ;
    virtual store_type RepType() const ;
  } ;

  typedef NPTR<MapRep> MapRepP ;
}

#endif

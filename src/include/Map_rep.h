#ifndef MAP_REP_H
#define MAP_REP_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#include <store_rep.h>

namespace Loci {

  class multiMap ;
  class Map ;
  
  class MapRep : public storeRep {
  public:
    virtual ~MapRep() ;
    virtual entitySet image(const entitySet &domain) const = 0 ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const = 0 ;
    virtual multiMap get_map() = 0 ;
    virtual void compose(const dMap &m, const entitySet &context) = 0 ;
    virtual store_type RepType() const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) = 0 ;
    virtual storeRepP thaw() = 0 ;
  } ;

  typedef NPTR<MapRep> MapRepP ;
}

#endif

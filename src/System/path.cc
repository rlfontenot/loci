#include "path.h"

#include <stream.h>

namespace Loci {
    
pathRep::pathRep()
{
    store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ;
}

pathRep::pathRep(const sequence &p)
{
    store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ;
    path = p ;
}

pathRep::~pathRep()
{
}

void pathRep::allocate(const entitySet &p)
{
    store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ;
    dispatch_notify() ;
}

storeRep *pathRep::new_store(const entitySet &p) const
{
    return new pathRep() ;
}

store_type pathRep::RepType() const
{
    return PATH ;
}

const entitySet &pathRep::domain() const {
    return store_domain ;
}

ostream &pathRep::Print(ostream &s) const {
    s << path << endl ;
    return s ;
}


istream &pathRep::Input(istream &s) {
    sequence e ;
    s >> path ;
    store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ;
    dispatch_notify() ;
    return s ;
}

path::path() {
    setRep(new pathType);
}
    
path::path(path &var) {
    setRep(var.Rep()) ;
}

path::path(const sequence &ptn) {
    setRep(new pathType(ptn));
}

path::~path() {}

void path::notification()
{
    NPTR<pathType> p(Rep());
    if(p!=0)
      data = p->get_path() ;
    warn(p==0);
}

}

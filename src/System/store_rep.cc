#include <store_rep.h>

#include <Tools/stream.h>
#include <Map.h>

namespace Loci {
  storeRep::~storeRep() {}
  void storeRep::set_elem_size(int sz) { warn(true) ; }

  storeRepP storeRep::getRep() { return storeRepP(this) ; }
  storeRepP storeRep::getRep() const { return storeRepP(const_cast<storeRep *>(this)) ; }
  store_instance::instance_type store_instance::access() const {
    return READ_WRITE ;
  }

  store_ref::~store_ref() {}

  void store_ref::allocate(const entitySet &ptn) {
    Rep()->allocate(ptn) ;
  }

  storeRep *store_ref::new_store(const entitySet &p) const  {
    return Rep()->new_store(p) ;
  }

  store_type store_ref::RepType() const {
    return Rep()->RepType() ;
  }

  storeRepP store_ref::remap(const dMap &m) const {
    return Rep()->remap(m) ;
  }

  void store_ref::copy(storeRepP &st, const entitySet &context) {
    Rep()->copy(st,context) ;
  }

  void store_ref::gather(const dMap &m,storeRepP &st,
                         const entitySet &context) {
    Rep()->gather(m,st,context) ;
  }

  void store_ref::scatter(const dMap &m, storeRepP &st,
                          const entitySet &context) {
    Rep()->scatter(m,st,context) ;
  }  
  int store_ref::pack_size(const entitySet &e )  {
    return(Rep()->pack_size(e)) ;
  }
  void store_ref::pack(void * ptr, int &loc, int &size, const entitySet &e )  {
    Rep()->pack(ptr, loc, size, e) ;
  }
  void store_ref::unpack(void * ptr, int &loc, int &size, const sequence &seq )  {
     Rep()->unpack(ptr, loc, size, seq) ;
   }

  ostream &store_ref::Print(ostream &s) const {
    return Rep()->Print(s) ;
  }

  istream &store_ref::Input(istream &s) {
    return Rep()->Input(s) ;
  }

  entitySet store_ref::domain() const {
    return Rep()->domain() ;
  }

  storeRepP store_ref::getRep() {
    return Rep()->getRep() ;
  }

  storeRepP store_ref::getRep() const {
    return Rep()->getRep() ;
  }
    
  void store_ref::notification() {
    dispatch_notify() ;
  }
}

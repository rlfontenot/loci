#include <store_rep.h>

#include <Tools/stream.h>


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

    ostream &store_ref::Print(ostream &s) const {
        return Rep()->Print(s) ;
    }

    istream &store_ref::Input(istream &s) {
        return Rep()->Input(s) ;
    }

    const entitySet &store_ref::domain() const {
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

#ifndef STORE_REP_H
#define STORE_REP_H

#include "debug.h"
#include "entitySet.h"
#include "nptr.h"

namespace Loci {
    enum store_type { STORE, PARAMETER, MAP, CONSTRAINT, PATH, GLOBAL } ;

    class storeRep ;
    typedef NPTR<storeRep> storeRepP ;
    typedef const_NPTR<storeRep> const_storeRepP ;
    
    class storeRep : public NPTR_type {
      public:
        virtual ~storeRep() ;
        virtual void allocate(const entitySet &p) = 0 ;
        virtual void set_elem_size(int sz) ;
        virtual storeRep *new_store(const entitySet &p) const = 0 ;
        virtual store_type RepType() const = 0 ;
        virtual std::ostream &Print(std::ostream &s) const = 0 ;
        virtual std::istream &Input(std::istream &s) = 0 ;
        virtual const entitySet &domain() const = 0 ;
        virtual storeRepP getRep() ;
        virtual storeRepP getRep() const ;
    } ;


    class store_instance : public eventNotify {
        friend class store_ref ;
        storeRepP rep ;
      public:
        enum instance_type { READ_WRITE, READ_ONLY } ;
        void setRep(const storeRepP &p)
        { rep = p ; rep.set_notify(this); notification() ; }
        storeRepP Rep() { return rep->getRep(); }
        storeRepP Rep() const { return rep->getRep(); }
        virtual instance_type access() const ;
    } ;

    class store_ref : public store_instance, public storeRep {
      public:
        store_ref() {}
        store_ref(storeRepP &p) { setRep(p); }
        virtual ~store_ref() ;

        store_ref &operator=(storeRepP &p) {
            setRep(p) ;
            return *this ;
        }
        
        store_ref &operator=(const store_ref &sr) {
            setRep(sr.rep) ;
            return *this ;
        }

        virtual void allocate(const entitySet &ptn) ;
        virtual storeRep *new_store(const entitySet &p) const ;
        virtual store_type RepType() const ;
        virtual std::ostream &Print(std::ostream &s) const ;
        virtual std::istream &Input(std::istream &s) ;
        virtual const entitySet &domain() const ;
        virtual storeRepP getRep() ;
        virtual storeRepP getRep() const ;
        
        virtual void notification() ;
    } ;
    typedef NPTR<store_ref> store_refP ;
    
}

#endif

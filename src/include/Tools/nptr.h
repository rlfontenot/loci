#ifndef NPTR_H
#define NPTR_H

#include "debug.h"
#include "eventNotify.h"
#include "cptr.h"

namespace Loci {

class NPTR_type : public CPTR_type {
    mutable eventDispatcher ed ;
  public:
    void engage(eventNotify *en) const { ed.engage(en) ; }
    void disengage(eventNotify *en) const {ed.disengage(en) ; }
    void dispatch_notify() const { ed.dispatch_notify() ; }
} ;

template <class T> class const_NPTR ;

template <class T> class NPTR {
    T *ptr ;
    eventNotify * receiver ;
    void unlink_ptr()
      { disengage() ; if(ptr) ptr->NPTR_type::unlink() ; ptr = 0 ; }
    void link_ptr() { if(ptr) ptr->NPTR_type::link() ; engage() ; }
    void set_ptr(T *p) {
        if(p)
          p->NPTR_type::link() ;
        unlink_ptr() ;
        ptr = p ;
        engage() ;
    }
    
    void engage() { if(ptr && receiver) ptr->engage(receiver) ; }
    void disengage() { if(ptr && receiver) ptr->disengage(receiver) ; }
  public:
    // OCF Methods
    template <class S> friend class NPTR ;
    template <class S> friend class const_NPTR ;
    
    NPTR(T *p=0,eventNotify *r=0) { ptr = 0 ; receiver=r ; set_ptr(p) ; }
    NPTR(const NPTR<T> &p,eventNotify *r=0)
    {ptr = 0 ; receiver=r ; set_ptr(p.ptr) ; }

    template <class S> explicit NPTR(const NPTR<S> &p, eventNotify *r=0)
    { ptr = 0 ; receiver=r ;  T *pt = dynamic_cast<T *>(p.ptr) ;
    set_ptr(pt) ; warn(pt == 0) ; }
    
    ~NPTR() { unlink_ptr(); }
    NPTR<T> &operator=(const NPTR<T> &p) { set_ptr(p.ptr) ; return *this; }
    NPTR<T> &operator=(T *p) { set_ptr(p); return *this; }

    bool operator==( const NPTR<T> &p) const { return  ptr == p->ptr; }
    bool operator==( const const_NPTR<T> &p) const { return ptr == p->ptr ; }
    bool operator==( const T *p ) const { return ptr == p ; }

    bool operator!=( const NPTR<T> &p) const { return  ptr != p->ptr; }
    bool operator!=( const const_NPTR<T> &p) const { return ptr != p->ptr ; }
    bool operator!=( const T *p ) const { return ptr != p ; }

    void set_notify(eventNotify *np)
      { disengage() ; receiver = np ; engage() ; }

    // Dereference Operator
    T &operator*() {return *ptr ; }
    const T &operator* () const { return *ptr ; }

    // Delegation Operator
    T * operator-> () {return ptr ; }
    const T * operator-> () const { return ptr ; }

    //    T * get_ptr() { return ptr ; }
    //    const T *get_ptr() const { return ptr ; }
} ;

template <class T> class const_NPTR {
    const T *ptr ;
    eventNotify * receiver ;
    void unlink_ptr() 
      { disengage() ; if(ptr) ptr->NPTR_type::unlink() ; ptr = 0 ; }
    void link_ptr() { if(ptr) ptr->NPTR_type::link() ; engage() ; }
    void set_ptr(const T *p) {
        if(p)
          p->NPTR_type::link() ;
        unlink_ptr() ;
        ptr = p ;
        engage() ;
    }
    void engage() const { if(ptr && receiver) ptr->engage(receiver) ; }
    void disengage() const { if(ptr && receiver) ptr->disengage(receiver) ; }
  public:
    // OCF Methods
    template <class S> friend class NPTR ;
    template <class S> friend class const_NPTR ;
    const_NPTR(const T *p=0,eventNotify *r=0)
    { ptr = 0 ; receiver=r ; set_ptr(p) ; }
    const_NPTR(const const_NPTR<T> &p)
    {ptr = 0 ; receiver=0 ; set_ptr(p.ptr) ; }
    const_NPTR(const NPTR<T> &p)
    {ptr =0 ; receiver = 0 ; set_ptr(p.ptr) ; }

    template <class S> explicit const_NPTR(const const_NPTR<S> &p) {
      ptr = dynamic_cast<T *>(p.ptr) ;
      warn(ptr==0) ;
      link_ptr() ;
    }

    template <class S> explicit const_NPTR(const NPTR<S> &p) {
      ptr = dynamic_cast<T *>(p.ptr) ;
      warn(ptr==0) ;
      link_ptr() ;
    }
    
    ~const_NPTR() { unlink_ptr(); }
    const_NPTR<T> &operator=(const const_NPTR<T> &p)
    { set_ptr(p.ptr) ; return *this; }
    const_NPTR<T> &operator=(const NPTR<T> &p)
    { set_ptr(p.ptr) ; return *this; }
    const_NPTR<T> &operator=(const T *p)
    { set_ptr(p); return *this; }

    bool operator==( const_NPTR<T> &p) const { return  ptr == p->ptr; }
    bool operator==( NPTR<T> &p) const { return ptr == p->ptr ; }
    bool operator==( const T *p ) const { return ptr == p ; }

    void set_notify(eventNotify *np)
      { disengage() ; receiver = np ; engage() ; }

    // Dereference Operator
    const T &operator* () const { return *ptr ; }

    // Delegation Operator
    const T * operator-> () const { return ptr ; }

    //    T * get_ptr() { return ptr ; }
    //    const T *get_ptr() const { return ptr ; }
} ;
}

#endif

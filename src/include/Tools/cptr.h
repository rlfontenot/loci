#ifndef CPTR_H
#define CPTR_H

#include <Tools/debug.h>

namespace Loci {
// Counted Pointer

class CPTR_type {
    mutable int count ;
    //    CPTR_type *operator&() { warn(true) ; return 0 ; }
  public:
    CPTR_type() {count = 0 ; }
    virtual ~CPTR_type() { warn(count!=0) ; }
    void link() const  { ++count ; }
    void unlink() const { --count ; if(count == 0) delete this ; }
} ;

template <class T> class const_CPTR ;

template <class T> class CPTR {
    T *ptr ;
    void unlink_ptr() { if(ptr) ptr->CPTR_type::unlink() ; ptr = 0 ; }
    void link_ptr() { if(ptr) ptr->CPTR_type::link() ; }
    void set_ptr(T *p) { if(p) p->CPTR_type::link() ; unlink_ptr() ; ptr = p ;}
  public:
    template <class S> friend class CPTR ;
    template <class S> friend class const_CPTR ;
    
    // OCF Methods
    CPTR(T *p = 0) { ptr = p ; link_ptr(); }
    CPTR(const CPTR<T> &p) {ptr = p.ptr;  ; link_ptr(); }
    template <class S> explicit CPTR(const CPTR<S> &p) {
      ptr = dynamic_cast<T *>(p.ptr) ;
      warn(ptr==0) ;
      link_ptr() ;
    }

    ~CPTR() { unlink_ptr(); }
    CPTR<T> &operator=(const CPTR<T> &p) { set_ptr(p.ptr) ; return *this; }
    CPTR<T> &operator=(T *p) { set_ptr(p); return *this; }

    bool operator==( const CPTR<T> &p) const { return  ptr == p.ptr; }
    bool operator==( const const_CPTR<T> &p) const { return  ptr == p.ptr; }
    bool operator==( const T *p ) const { return ptr == p ; }

    bool operator!=( const CPTR<T> &p) const { return  ptr != p.ptr; }
    bool operator!=( const const_CPTR<T> &p) const { return  ptr != p.ptr; }
    bool operator!=( const T *p ) const { return ptr != p ; }

    // Dereference Operator
    T &operator*() {return *ptr ; }
    const T &operator* () const { return *ptr ; }

    // Delegation Operator
    T * operator-> () {return ptr ; }
    const T * operator-> () const { return ptr ; }
} ;

template <class T> class const_CPTR {
    const T *ptr ;
    void unlink_ptr() { if(ptr) ptr->CPTR_type::unlink() ; ptr = 0 ; }
    void link_ptr() { if(ptr) ptr->CPTR_type::link() ; }
    void set_ptr(T *p) { if(p) p->CPTR_type::link() ; unlink_ptr() ; ptr = p ; }
  public:
    template <class S> friend class CPTR ;
    template <class S> friend class const_CPTR ;
    
    // OCF Methods
    const_CPTR(const T *p = 0) { ptr = p ; link_ptr(); }
    const_CPTR(const CPTR<T> &p) {ptr = p.ptr;  ; link_ptr(); }
    const_CPTR(const const_CPTR<T> &p) {ptr = p.ptr;  ; link_ptr(); }
    template <class S> explicit const_CPTR(const CPTR<S> &p) {
      ptr = dynamic_cast<T *>(p.ptr) ;
      warn(ptr==0) ;
      link_ptr() ;
    }
    template <class S> explicit const_CPTR(const const_CPTR<S> &p) {
      ptr = dynamic_cast<T *>(p.ptr) ;
      warn(ptr==0) ;
      link_ptr() ;
    }

    ~const_CPTR() { unlink_ptr(); }
    const_CPTR<T> &operator=(const CPTR<T> &p)
    { set_ptr(p.ptr) ; return *this; }
    const_CPTR<T> &operator=(const const_CPTR<T> &p)
    { set_ptr(p.ptr) ; return *this; }
    const_CPTR<T> &operator=(const T *p)
    { set_ptr(p); return *this; }

    bool operator==( const CPTR<T> &p) const { return  ptr == p.ptr; }
    bool operator==( const const_CPTR<T> &p) const { return  ptr == p.ptr; }
    bool operator==( const T *p ) const  { return ptr == p ; }

    bool operator!=( const CPTR<T> &p) const { return  ptr != p.ptr; }
    bool operator!=( const const_CPTR<T> &p) const { return  ptr != p.ptr; }
    bool operator!=( const T *p ) const  { return ptr != p ; }

    // Dereference Operator
    const T &operator* () const { return *ptr ; }

    // Delegation Operator
    const T * operator-> () const { return ptr ; }
} ;

}

#endif




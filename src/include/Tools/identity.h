#ifndef IDENTITY_H
#define IDENTITY_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#include <typeinfo>

#include <Tools/debug.h>
#include <Tools/lmutex.h>

namespace Loci {
typedef unsigned long id_type ;

class IdentityIMPL {
  id_type id ;
  IdentityIMPL & operator=(const IdentityIMPL &i)
  { warn(true) ; return *this; }
  public:
    IdentityIMPL() { id = reinterpret_cast<id_type>(this) ; }
    IdentityIMPL(const IdentityIMPL &i) { id = reinterpret_cast<id_type>(this);} 
    ~IdentityIMPL() {} ;
    id_type ident() { return id ; }
} ;

class Identity : virtual public IdentityIMPL {} ;
}

#endif

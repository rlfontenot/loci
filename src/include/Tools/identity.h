#ifndef IDENTITY_H
#define IDENTITY_H

#include <Tools/debug.h>

namespace Loci {
typedef unsigned long id_type ;

class IdentityIMPL {
    id_type id ;
    static id_type id_ ;
    IdentityIMPL & operator=(const IdentityIMPL &i)
    { warn(true) ; return *this; }
  public:
    IdentityIMPL() { id = ++id_ ; }
    IdentityIMPL(const IdentityIMPL &i) { id = ++id_; } 
    ~IdentityIMPL() {} ;
    id_type ident() { return id ; }
} ;

class Identity : virtual public IdentityIMPL {} ;
}

#endif

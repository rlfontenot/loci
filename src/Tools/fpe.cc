#include <Tools/fpe.h>


#ifdef SPARC

extern "C" {
  extern void SPARC_set_fpe_abort() ;
}

namespace Loci {
  void set_fpe_abort() {
    SPARC_set_fpe_abort() ;
  }
}
#else

#ifdef SGI

#include <stdio.h>
#include <sigfpe.h>

namespace Loci {
  void set_fpe_abort() {
    sigfpe_[_DIVZERO].abort = 1 ;
    sigfpe_[_OVERFL].abort = 1 ;
    sigfpe_[_INVALID].abort = 1 ;
    handle_sigfpes(_DEBUG,_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
                   0,_ABORT_ON_ERROR,0) ;
  }
}
                      
#else
namespace Loci {
  
  void set_fpe_abort()
    {
    }
}
#endif
#endif


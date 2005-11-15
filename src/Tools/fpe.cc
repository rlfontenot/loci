#include <Tools/fpe.h>
#include <Tools/debugger.h>
#include <stdlib.h>


#ifdef SPARC

extern "C" {
  extern void SPARC_set_fpe_abort() ;
  void fpe_debugger_() { Loci::debugger_() ; }
}

namespace Loci {
  void set_fpe_abort() {
    SPARC_set_fpe_abort() ;
  }
}
#else

#ifdef SGI

#include <iostream>
#include <stdio.h>
#include <sigfpe.h>
extern "C" {
  void sgi_fpe_user_(unsigned int v[5], int y[2]) {
    std::cerr << "Floating Point Exception" << std::endl ;
    Loci::debugger_() ;
  }
  void sgi_fpe_abort_(unsigned int **p) {
    std::cerr << "Floating Point Exception" << std::endl ;
    Loci::debugger_() ;
  }
}

namespace Loci {
  void set_fpe_abort() {
    sigfpe_[_DIVZERO].abort = 1 ;
    sigfpe_[_OVERFL].abort = 1 ;
    sigfpe_[_INVALID].abort = 1 ;
    handle_sigfpes(_ON,_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
                   sgi_fpe_user_,_ABORT_ON_ERROR,sgi_fpe_abort_) ;
                   //                   0,_ABORT_ON_ERROR,0) ;
  }
}
                      
#else
#ifdef LINUX

#include <stdio.h>
#include <signal.h>
#include <iostream>

// Floating point exception environment
#ifndef _GNU_SOURCE
#define _GNU_SOURCE     /* needed to get non POSIX extensions. */
#endif
#include <fenv.h>

using std::cerr ;
using std::endl ;

extern "C" {
  void fpe_debugger_(int i)
  {
    std::cerr << "floating point exception" << std::endl ;
    Loci::debugger_() ;
  }
}

namespace Loci {
  
  void set_fpe_abort()
  {
#ifdef SYSTEM_ITANIUM64
    
    if(feenableexcept((FE_DIVBYZERO)) == -1) {
      std::cerr << "feenableexcept had error, floating point exceptions not caught" << std::endl ;
    } else {
      signal(SIGFPE,fpe_debugger_) ;
    }
#else
    if(feenableexcept((FE_DIVBYZERO|FE_OVERFLOW|FE_INVALID)) == -1) {
      std::cerr << "feenableexcept had error, floating point exceptions not caught" << std::endl ;
    } else {
      signal(SIGFPE,fpe_debugger_) ;
    }
#endif
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
#endif


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
#ifdef HAVE_FENWM


#include <fenvwm.h>
#include <iostream>
using std::cerr ;
using std::endl ;

namespace Loci {
  int report_overflow(struct sigcontext_struct *sig,
                      struct i387_hard_struct *fpu)
  {
    cerr << "FPU Overflow exception occured at "
         << (void *) fpu->fip << endl ;
    abort() ;
    return 0 ;
  }
  int report_divzero(struct sigcontext_struct *sig,
                      struct i387_hard_struct *fpu)
  {
    cerr << "FPU Divide by Zero exception occured at "
         << (void *) fpu->fip << endl ;
    abort() ;
    return 0 ;
  }
  int report_invalid(struct sigcontext_struct *sig,
                     struct i387_hard_struct *fpu)
  {
    cerr << "FPU INVALID operation exception occured at "
         << (void *) fpu->fip << endl ;
    abort() ;
    return 0 ;
  }

  ftrap_t overflow = (ftrap_t) &report_overflow;
  ftrap_t divzero  = (ftrap_t) &report_divzero ;
  ftrap_t invalid  = (ftrap_t)&report_invalid ;
  
  void set_fpe_abort() {
    fesetvector(&overflow,FE_OVERFLOW) ;
    feenabletraps(FE_OVERFLOW) ;
    fesetvector(&divzero,FE_DIVBYZERO) ;
    feenabletraps(FE_DIVBYZERO) ;
    fesetvector(&invalid,FE_INVALID) ;
    feenabletraps(FE_INVALID) ;
  }
}  
#else
#ifdef LINUX
#include <fpu_control.h>
#include <signal.h>
#include <iostream>

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
    fpu_control_t cw =
      _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
    _FPU_SETCW(cw);
    signal(SIGFPE,fpe_debugger_) ;
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
#endif


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
#ifdef LINUX_X86

#include <stdlib.h>
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
namespace Loci {
  
  void set_fpe_abort()
    {
    }
}
#endif
#endif
#endif


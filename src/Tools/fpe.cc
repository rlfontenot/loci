#include "fpe.h"

namespace Loci {

#ifdef NO_FPE_EXCEPTIONS

void set_fpe_abort()
{
}

#else
#include <stdio.h>
#include <sigfpe.h>
//#include <signal.h>

void set_fpe_abort()
{
    sigfpe_[_DIVZERO].abort = 1 ;
    sigfpe_[_OVERFL].abort = 1 ;
    sigfpe_[_INVALID].abort = 1 ;
    handle_sigfpes(_DEBUG,_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
                   0,_ABORT_ON_ERROR,0) ;
}

#endif

}

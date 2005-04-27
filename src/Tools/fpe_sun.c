#include <Config/conf.h>
#ifdef SPARC
#include <stdio.h>
#include <stdlib.h>
#include <floatingpoint.h>
#include <siginfo.h>
#include <ucontext.h>
#include <sunmath.h>

extern void fpe_debugger_() ;

void SPARC_ieee_abort(int sig, siginfo_t *sip, ucontext_t *uap) {
  char *label;
  
  switch (sip->si_code) {
  case FPE_FLTINV: label = "invalid operand"; break;
  case FPE_FLTRES: label = "inexact"; break;
  case FPE_FLTDIV: label = "division-by-zero"; break;
  case FPE_FLTUND: label = "underflow"; break;
  case FPE_FLTOVF: label = "overflow"; break;
  default: label = "???"; break;
  }
    
  fprintf(stderr, "FP exception %s (0x%x) occurred at address %p.\n",
          label, sip->si_code, (void *) sip->si_addr);
  fpe_debugger_() ;
  abort();
}

void SPARC_set_fpe_abort() {
  sigfpe_handler_type hdl ;
  hdl = (sigfpe_handler_type) SPARC_ieee_abort ;
  ieee_handler("set","common",hdl) ;
}
#else
void SOARC_set_fpe_abort() {
}
#endif


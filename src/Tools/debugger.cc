#include <stdio.h>
#include <signal.h>
#ifdef SUN
#include <floatingpoint.h>
#include <siginfo.h>
#endif
#include <ucontext.h>

#ifdef SUN
#include <sys/machsig.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>

#include <stdlib.h>
#include <unistd.h>

#include "Tools/debugger.h"

#define HOST_ID "localhost"

namespace Loci {
  char *debug_hostname = HOST_ID ;
  char *debug_execname = "a.out" ;
  char *debug_program = "gdb" ;
  bool debugger_setup = false ;

  void setup_debugger(char *execname, char *gdb, char *hostname) {
    debug_execname = execname ;
    debug_program = gdb ;
    debug_hostname = hostname ;
    debugger_setup = true ;
  }
  
  void debugger_()
  {
    if(!debugger_setup)
      abort() ;
    
    int pid = getpid() ;
    char buf[512] ;
    int breakout ;
    struct stat sb ;
#ifdef SUN
    char *xtermpath = "/usr/openwin/bin/xterm" ;
    char *xtermlibpath = "/usr/openwin/lib:/usr/lib" ;
#endif
#ifdef LINUX
    char *xtermpath = "/usr/X11R6/bin/xterm" ;
    char *xtermlibpath = "/usr/X11R6/lib:/usr/lib:/lib/i686:/lib" ;
#endif
    
    sprintf(buf,"export LD_LIBRARY_PATH;LD_LIBRARY_PATH=%s; %s  -display %s:0 -e %s %s %d &",
            xtermlibpath,xtermpath,
            debug_hostname,debug_program,debug_execname,pid) ;
    system(buf) ;
    breakout = 1 ;
    sleep(10) ; /* Wait for debugger to attach */
  }



  void program_trap(int sig) /*,code,scp,addr)*/
  {
    char *sigtype = "(undefined)" ;

    switch(sig) {
    case SIGBUS:
      sigtype = "a Bus Error" ;
      break ;
    case SIGSEGV:
      sigtype = "a Segmentation Violation" ;
      break ;
    case SIGILL:
      sigtype = "an Illegal Instruction Call" ;
      break ;
    case SIGSYS:
      sigtype = "an Illegal System Call" ;
      break ;
    case SIGFPE:
      sigtype = "a Floating Point Exception" ;
      break ;
    }
    fprintf(stderr,"ERROR: Program terminated due to %s\n",sigtype) ;
    debugger_() ;
    exit(-1) ;
  }

#ifdef SUN
  void FPE_program_trap(int sig, siginfo_t *sip, ucontext_t *uap) {
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
    debugger_() ;
    exit(-1) ;
  }
#endif

  void chopsigs_()
  {
#ifdef SUN
    sigfpe_handler_type hdl ;
#endif
  
    signal(SIGBUS,program_trap) ;
    signal(SIGSEGV,program_trap) ;
    signal(SIGILL,program_trap) ;
    signal(SIGSYS,program_trap) ;
#ifdef LINUX
    signal(SIGFPE,program_trap) ;
#endif
#ifdef SUN  
    hdl = (sigfpe_handler_type) FPE_program_trap ;
    sigfpe(FPE_FLTINV,FPE_program_trap) ;
    sigfpe(FPE_FLTDIV,FPE_program_trap) ;
    sigfpe(FPE_FLTOVF,FPE_program_trap) ;
#endif
  }


}  

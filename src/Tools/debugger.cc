#ifdef SPARC
#define SUN
#endif

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

#include <mpi.h>


#include <iostream>

using std::cerr ;
using std::endl ;

#define HOST_ID "localhost"

namespace Loci {
  const char *debug_hostname = HOST_ID ;
  const char *debug_execname = "a.out" ;
  const char *debug_program = "gdb" ;
  bool debugger_setup = false ;

  void setup_debugger(const char *execname, const char *gdb, const char *hostname) {
    debug_execname = execname ;
    debug_program = gdb ;
    debug_hostname = hostname ;
    debugger_setup = true ;
  }
  
  void debugger_()
  {
    int MPI_rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    cerr << "program failed on processor " << MPI_rank << endl ;

    if(!debugger_setup)
      MPI_Abort(MPI_COMM_WORLD,-1);
    
    int pid = getpid() ;
    char buf[512] ;
    int breakout ;
    struct stat sb ;
#ifdef SUN
    const char *xtermpath = "/usr/openwin/bin/xterm" ;
    const char *xtermlibpath = "/usr/openwin/lib:/usr/lib" ;
#endif
#ifdef LINUX
    const char *xtermpath = "/usr/X11R6/bin/xterm" ;
    const char *xtermlibpath = "/usr/X11R6/lib:/usr/lib:/lib/i686:/lib" ;
#endif
#ifdef SGI
    const char *xtermpath = "/usr/bin/X11/xterm" ;
    const char *xtermlibpath = "/usr/lib" ;
#endif
#ifndef SGI
    sprintf(buf,"export LD_LIBRARY_PATH;LD_LIBRARY_PATH=%s; %s  -display %s -e %s %s %d &",
            xtermlibpath,xtermpath,
            debug_hostname,debug_program,debug_execname,pid) ;
#else
    sprintf(buf,"export LD_LIBRARY_PATH;LD_LIBRARY_PATH=%s; %s  -display %s -e %s -p %d %s &",
            xtermlibpath,xtermpath,
            debug_hostname,debug_program,pid,debug_execname) ;
#endif
    system(buf) ;
    breakout = 1 ;
    sleep(100) ; /* Wait for debugger to attach */
  }



  void program_trap(int sig) /*,code,scp,addr)*/
  {
    const char *sigtype = "(undefined)" ;

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
    MPI_Abort(MPI_COMM_WORLD,-1) ;
  }


  void chopsigs_()
  {
    signal(SIGBUS,program_trap) ;
    signal(SIGSEGV,program_trap) ;
    signal(SIGILL,program_trap) ;
    signal(SIGSYS,program_trap) ;
#ifdef LINUX
    signal(SIGFPE,program_trap) ;
#endif
  }
}  

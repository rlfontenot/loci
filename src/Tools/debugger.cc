//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Config/conf.h>

#include <stdio.h>
#include <signal.h>
#ifdef SUN
#include <floatingpoint.h>
#include <siginfo.h>
#endif

#ifndef __CYGWIN__
#include <ucontext.h>
#endif

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
#ifdef SGI
  const char *debug_program = "dbx" ;
#else
  const char *debug_program = "gdb" ;
#endif
  
  bool debugger_setup = false ;

  void (*debug_callback)() = 0 ;

  void set_debug_callback(void (*dcb)()) {
    debug_callback = dcb ;
  }

  
  void setup_debugger(const char *execname, const char *gdb, const char *hostname) {
    debug_execname = execname ;
    debug_program = gdb ;
    debug_hostname = hostname ;
    debugger_setup = true ;
  }
  
  void debugger_()
  {


    void (*dbc)() = debug_callback ;
    debug_callback = 0 ;

    int MPI_rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    cerr << "program failed on processor " << MPI_rank << endl ;

    if(dbc == 0)
      MPI_Abort(MPI_COMM_WORLD,-1) ;

    if(dbc != 0)
      (*dbc)() ;

    if(!debugger_setup)
      MPI_Abort(MPI_COMM_WORLD,-1);
    
    int pid = getpid() ;
    char buf[512] ;
    const char *xtermpath = "xterm" ;
#ifdef BSD
    exit(-1) ;
#else
#ifndef SGI
    sprintf(buf,"%s  -display %s -e %s %s %d &",
            xtermpath,
            debug_hostname,debug_program,debug_execname,pid) ;
#else
    sprintf(buf,"%s  -display %s -e %s -p %d %s &",
            xtermpath,
            debug_hostname,debug_program,pid,debug_execname) ;
#endif
    cerr << buf << endl ;
    int err = system(buf) ;
    if(err != 0)
      cerr << "system call failed with on '"<< buf << "'" << endl ;


    sleep(100) ; /* Wait for debugger to attach */
#endif
  }

}

extern "C" {
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
    Loci::debugger_() ;
    MPI_Abort(MPI_COMM_WORLD,-1) ;
  }
}


namespace Loci {
  

  void chopsigs_()
  {
    signal(SIGBUS,program_trap) ;
    signal(SIGSEGV,program_trap) ;
    signal(SIGILL,program_trap) ;
    signal(SIGSYS,program_trap) ;
    signal(SIGFPE,program_trap) ;
    signal(SIGABRT,program_trap) ;
  }
}  

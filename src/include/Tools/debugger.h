#ifndef DEBUGGER_H
#define DEBUGGER_H

namespace Loci {
  void setup_debugger(const char *execname,const char *gdb, const char *hostname) ;
  void debugger_() ;
  void chopsigs_() ;
}



#endif

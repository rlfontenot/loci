#include <rule.h>
#include "dist_tools.h"
#include "loci_globs.h"
#include <Tools/debug.h>

#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <unistd.h>

#include "Tools/debugger.h"

#include <mpi.h>

#include <iostream>
#include <fstream>
using std::cout ;
using std::cerr ; 
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::ostream ;
using std::ofstream ;
using std::string ;
using std::ostringstream ;
using std::cout ;

namespace Loci {
  int MPI_processes = 1;
  int MPI_rank ;
  int num_threads = 1 ;
  int method = 3 ;
  /////////////////////////////
  // flags to turn on/off the visualization feature
  bool show_graphs = false ;
  bool show_decoration = false ;
  // flag to enable/disable the dynamic memory management
  bool use_dynamic_memory = true ;
  // flag to enable/disable output of dynamic memory
  // and multilevel graph decoration information
  bool show_dmm_verbose = false ;
  // flag to enable/disable chomping
  bool use_chomp = true ;
  // flag to visualize the chomping graph
  bool show_chomp = false ;
  // flag to turn on the summary report on chomping
  bool chomp_verbose = false ;
  // flag to enable outputing schedule to file
  bool schedule_output = false ;
  // flag to enable dynamic scheduling feature
  bool use_dynamic_scheduling = false ;
  // chomping size(unit is KB), default is 128KB
  int chomping_size = 128 ;
  // flag to enable memory profiling
  bool profile_memory_usage = false ;
  // flag to use a different scheduler
  // (more memory conservative and more synchronization
  //  points will be generated)
  bool memory_greedy_schedule = true ;
  bool use_old_dependency_graph = false ;
  /////////////////////////////
  
  ofstream debugout ;

  double total_memory_usage = 0 ;

  extern int current_rule_id ;

  void debug_print_rule() {
    if(current_rule_id != 0) {
      rule r(current_rule_id) ;
      cerr << "crash occured in rule " << r << endl ;

      if(exec_current_fact_db != 0) {
        char buf[512] ;
        sprintf(buf,"crash_dump.%d",MPI_rank) ;
        ofstream cfile(buf,ios::out) ;
        cfile << "rule: " << r << endl ;
        
        variableSet v = r.sources() ;

        cfile << "facts = {" << endl ;
        for(variableSet::const_iterator vi=v.begin();vi!=v.end();++vi) {
          cfile << *vi << ':' ;
          storeRepP p = Loci::exec_current_fact_db->get_variable(*vi);
          p->Print(cfile) ;
        }
        cfile << "}" << endl ;
      }
    }
  }
  //This is the first call to be made for any Loci program be it
  //sequential or parallel. 
  void Init(int* argc, char*** argv)  {
    char *execname = (*argv)[0] ;
    const char *hostname = "localhost" ;
    const char *debug = "gdb" ;
    //Setting up of the global variables for processor ID and the
    //total number of processes.  
    MPI_Init(argc, argv) ;
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN) ;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    //Create a debug file for each process
    ostringstream oss ;
    // if output directory doesn't exist, create one
    bool debug_is_directory = true ;
    struct stat statbuf ;
    if(GLOBAL_OR(stat("debug",&statbuf)!=0)) {
      if(MPI_rank == 0)
        mkdir("debug",0755) ;
      for(int i=0;i<1000;++i) {
        if(GLOBAL_AND(stat("debug",&statbuf)==0))
          break ;
      }
    } else {
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file 'debug' should be a directory!, rename 'output' and start again."
             << endl ;
        debug_is_directory = false ;

      }
    }

    if(debug_is_directory) {
      if(MPI_processes == 1)
        oss << "debug/debug" ;
      else
        oss << "debug/debug."<<MPI_rank ;
    } else
      oss << "debug."<< MPI_rank ;
    
    string filename  = oss.str() ;
    debugout.open(filename.c_str(),ios::out) ;
    // All the rules in an unnamed namespace are first copied into the 
    // global rule list. To add rules to the rule database we just
    // neeed to use the global_rule_list. Inititally when the rules
    // are registered using the register rule it gets pushed into the
    // register_rule_list which is a static rule list. 
    if(!register_rule_list.empty()) {
      global_rule_list.copy_rule_list(register_rule_list) ;
      register_rule_list.clear() ;
    }
    bool debug_setup = false ;
    int i = 1 ;
    while(i<*argc) {
      if(!strcmp((*argv)[i],"--display")) {
        debug_setup = true ;
        hostname = (*argv)[i+1] ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--debug")) {
        debug = (*argv)[i+1] ;
        schedule_output = true ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--threads")) {
        cerr << "warning --threads not yet implemented" << endl ;
        num_threads = atoi((*argv)[i+1]) ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--scheduleoutput")) {
        schedule_output = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--graphs")) {
        show_graphs = true ; // visualize the dependency graph &
                             // the decomposed graph & every supernode
        i++ ;
      } else if(!strcmp((*argv)[i],"--decoration")) {
        show_decoration = true ; // visualize the decorated multilevel graph
        i++ ;
      } else if(!strcmp((*argv)[i],"--dmm")) {
        use_dynamic_memory = true ; // use the dynamic memory management
        i++ ;
      } else if(!strcmp((*argv)[i],"--nodmm")) {
        use_dynamic_memory = false ; // use the dynamic memory management
        i++ ;
      } else if(!strcmp((*argv)[i],"--dmmverbose")) {
        // output some info about dmm
        show_dmm_verbose = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chomp")) {
        // use the chomping scheme
        use_chomp = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--nochomp")) {
        // use the chomping scheme
        use_chomp = false ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--showchomp")) {
        // visualize the chomp graph
        show_chomp = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chompverbose")) {
        // summary report
        chomp_verbose = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chompingsize")) {
        std::stringstream ss ;
        ss << (*argv)[i+1] ;
        ss >> chomping_size ;
        if(chomping_size > 500*1024*1024/*500MB*/)
          cerr << "WARNING: setting chompingsize too large"
               << " may crash the program" << endl ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--memprofile")) {
        // profiling memory usage at run time
        profile_memory_usage = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--memgreedy")) {
        memory_greedy_schedule = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--nomemgreedy")) {
        memory_greedy_schedule = false ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--method")) {
        method = atoi((*argv)[i+1]);
        i+=2;
      } else if(!strcmp((*argv)[i],"--balance")) {
        use_dynamic_scheduling = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--olddg")) {
        use_old_dependency_graph = true ;
        i++ ;
      } else
        break ;
    }
    if(i!=1) {
      *argc -= (i-1) ;
      for(int k=1;k<*argc;++k)
        (*argv)[k] = (*argv)[k+i-1] ;
    }
    
    set_debug_callback(debug_print_rule) ;
    if(debug_setup) {
      setup_debugger(execname,debug,hostname) ;
    }
    chopsigs_() ;
  }
  //All Loci programs must end with this call. 
  void Finalize() {
    MPI_Finalize() ;
  }

  void Abort() {
    debugger_() ;
  }

}

#ifndef LOCI_GLOBS_H
#define LOCI_GLOBS_H

#include <fstream>
#include <string>
#include <vector>

namespace Loci {
  extern int num_threads ;
  /////////////////////////////
  // flags to turn on/off the visualization feature
  extern bool verbose ;
  extern bool show_graphs ;
  extern bool show_decoration ;
  // flag to enable/disable the dynamic memory management
  extern bool use_dynamic_memory ;
  // flag to enable/disable output of dynamic memory
  // and multilevel graph decoration information
  extern bool show_dmm_verbose ;
  // flag to enable/disable chomping
  extern bool use_chomp ;
  // flag to enable outputing schedule to file
  extern bool schedule_output ;
  //flag for duplicating computations to save communication
  extern bool duplicate_work;
  extern bool multilevel_duplication;
  extern bool reduction_duplication;
  extern bool pointwise_duplication;
  extern bool extended_duplication;
  extern bool collect_timings;
  extern std::ofstream timeout;
  extern double time_duration_to_collect_data;
  extern bool use_duplicate_model;
  extern char * model_file;
  /////////////////////////////
  extern bool random_partition;

  extern bool measure_rule_timings;
  extern std::ofstream ruleTimeOut;
  extern int printLevel ;
  inline std::ostream &printIndent(std::ostream &s) {
    for(int i=0;i<printLevel;++i)
      s << "  " ;
    return s ;
  }
      
}

#endif

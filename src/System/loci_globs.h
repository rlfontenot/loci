#ifndef LOCI_GLOBS_H
#define LOCI_GLOBS_H

namespace Loci {
  extern int num_threads ;
  /////////////////////////////
  // flags to turn on/off the visualization feature
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
  /////////////////////////////
}

#endif

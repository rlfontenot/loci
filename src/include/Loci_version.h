#ifndef LOCI_VERSION_H
#define LOCI_VERSION_H
#include <string>

// Tag that this is a distributed memory version of Loci
#define LOCI_DISTRIBUTED_MEMORY

// Tag that this is a distribution with a revised fact database distribution
// interface

#define LOCI_MODIFIED_FACT_DISTRIBUTE

namespace Loci {
  std::string version() ;
}

#endif

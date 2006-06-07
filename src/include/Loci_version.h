#ifndef LOCI_VERSION_H
#define LOCI_VERSION_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <string>

// Tag that this is a distributed memory version of Loci
#define LOCI_DISTRIBUTED_MEMORY

// Tag that this is a distribution with a revised fact database distribution
// interface

#define LOCI_MODIFIED_FACT_DISTRIBUTE

#define LOCI_VERSION_3_1_RELEASE

namespace Loci {
  std::string version() ;
}

#endif

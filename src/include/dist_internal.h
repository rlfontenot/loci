#ifndef DIST_INTERNAL_H
#define DIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <hdf5_readwrite.h>
#include <vector>
namespace Loci {

  void read_multi_vector_int(hid_t group_id, const char* name, int dim,  std::vector<int>& vint) ;
  void read_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size) ;
  void write_vector_int(hid_t group_id, const char* name, std::vector<int>& vint) ;

  std::vector<int> all_collect_sizes(int size) ;

}

#endif

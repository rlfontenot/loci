#ifndef FTRN_READER_H
#define FTRN_READER_H 1
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#include <Config/conf.h>

namespace Loci {
class fortran_binary_file {
    int fd ;
    unsigned long record_size ;
    const char *error ;
  public:
    fortran_binary_file(const char *filename);
    ~fortran_binary_file() ;
    const char *error_msg() { return error; }
    int get_array(int *num_array, unsigned long  sz) ;
    int get_array(float *num_array,unsigned long  sz) ;
    int get_array(double *num_array,unsigned long  sz) ;
} ;
}

#endif

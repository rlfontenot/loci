#ifndef FTRN_READER_H

#define FTRN_READER_H 1

namespace Loci {
class fortran_binary_file {
    int fd ;
    unsigned long record_size ;
    char *error ;
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

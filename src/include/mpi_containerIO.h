#ifndef MPI_CONTAINERIO_H
#define MPI_CONTAINERIO_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <dist_internal.h>
#include <distribute_container.h>
#include <distribute.h>
#include <distribute_io.h>
#include <mpi.h>
namespace Loci {
  /*
    mpi version to read/write stores. parallel io only
    Notice: each variable is read in from /write into a separate file
    argument xfer_type is used to control collective or independent data transfer
    global variable PHDF5_MPI_Info is used.
  */
  void pmpi_redistribute_write_container(std::string& filename,
                                         Loci::storeRepP var, Loci::fact_db &facts, int xfer_type);
  void pmpi_read_container_redistribute(std::string& filename,
                                        Loci::storeRepP var, Loci::entitySet read_set,
                                        Loci::fact_db &facts, int xfer_type);
  inline void  pmpi_writeContainer(std::string& filename, Loci::storeRepP var, Loci::fact_db &facts, int xfer_type){
    pmpi_redistribute_write_container(filename, var, facts, xfer_type) ;
  }

  inline void pmpi_readContainer(std::string& filename, Loci::storeRepP var, Loci::entitySet readSet, Loci::fact_db &facts, int xfer_type) {
    pmpi_read_container_redistribute(filename, var, readSet, facts, xfer_type) ;
  }
}
#endif

#ifndef DISTRIBUTE_IO_H
#define DISTRIBUTE_IO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <store_rep.h>
#include <DMap.h>
#include <fact_db.h>

namespace Loci {

  void write_container(hid_t group_id, storeRepP qrep) ;
  void read_container(hid_t group_id, storeRepP qrep, entitySet &dom) ;

  storeRepP collect_reorder_store(storeRepP &sp, dMap &remap, fact_db &facts) ;
  void distribute_reorder_store(storeRepP &new_sp, storeRepP sp_init, dMap &remap, fact_db &facts) ;

  void redistribute_write_container(hid_t file_id, std::string vname,
                                    Loci::storeRepP var, fact_db &facts) ;
  void read_container_redistribute(hid_t file_id, std::string vname,
                                   Loci::storeRepP var, entitySet read_set,
                                   fact_db &facts) ;

  inline hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id) {
    if(Loci::MPI_rank==0)
      return H5Fcreate(name,flags,create_id,access_id) ;
    else
      return 0 ;
  }

  inline hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id) {
    if(Loci::MPI_rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
  }

  inline herr_t hdf5CloseFile(hid_t file_id) {
    if(Loci::MPI_rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }
    

  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var) {
    if(Loci::exec_current_fact_db == 0) {
      fact_db facts ;// Hack!!
      redistribute_write_container(file_id,vname,var,facts) ;
    } else
      redistribute_write_container(file_id,vname,var,
                                   *Loci::exec_current_fact_db) ;
  }
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet) {
    if(Loci::exec_current_fact_db == 0) {
      fact_db facts ; // Hack!!
      read_container_redistribute(file_id,vname,var,readSet,facts) ;
    } else
      read_container_redistribute(file_id,vname,var,readSet,
                                  *Loci::exec_current_fact_db) ;
  }

                             
}

#endif

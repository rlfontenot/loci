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

}

#endif

#ifndef DISTRIBUTE_OLD_SUPPORT_H
#define DISTRIBUTE_OLD_SUPPORT_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include "dist_tools.h"
namespace Loci {

  //from dist_tools.h //////
  std::vector<entitySet> read_partition(const char *fname,int num_partitions) ;
  void write_partition(const char *fname, const std::vector<entitySet> &ptn) ;

  std::vector<entitySet> generate_scalable_distribution(fact_db &facts, rule_db &rdb, int num_partitions = 0) ;

  void print_global(entitySet e, fact_db &facts) ;

  Map distribute_whole_map(Map &m) ;

  //from distribute.h //////
  Map distribute_global_map(Map &m, const std::vector<entitySet> &vset) ;
  Map distribute_gmap(Map &m, const std::vector<entitySet> &vset) ;

  std::pair<storeRepP, Map > send_clone(storeRepP& sp, Map &m,  entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> fill_global_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;

}

#endif

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <typeinfo>
#include <vector>
#include <store.h>
#include <Map.h>
#include <DMap.h>
#include <DMultiMap.h>
#include <constraint.h>
#include <fact_db.h>

#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>

#define LOCI_HAS_READ_WRITE_PARTITION 1
 
namespace Loci {

  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  
  void Init(int* argc, char*** argv) ;
  void Finalize() ; 
  void Abort() ;
  
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn, int num_partitions = 0 ) ;
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<dMap> send_global_map(Map &attrib_data, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  void fill_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  std::pair<storeRepP, Map > send_clone(storeRepP& sp, Map &m,  entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  storeRepP send_clone_non(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> fill_global_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  std::vector<entitySet> generate_distribution(fact_db &facts, rule_db &rdb, int num_partitions = 0) ;

  void distribute_facts(std::vector<entitySet> &ptn, fact_db &facts,
                        rule_db &rdb) ;
  
  entitySet collect_entitySet(entitySet e, fact_db &facts) ;
  storeRepP collect_store(storeRepP &sp, fact_db &facts) ;
  storeRepP collect_global_store(storeRepP &sp) ;
  Map distribute_global_map(Map &m, const std::vector<entitySet> &vset) ;
  Map distribute_gmap(Map &m, const std::vector<entitySet> &vset) ;
  storeRepP collect_reorder_store(storeRepP &sp, dMap &remap, fact_db &facts) ;
  void distribute_reorder_store(storeRepP &new_sp, storeRepP sp_init, dMap &remap, fact_db &facts) ;
  
  void write_container(hid_t group_id, storeRepP qrep) ;
  void read_container(hid_t group_id, storeRepP qrep, entitySet &dom) ;

  storeRepP distribute_store(storeRepP &sp, fact_db &facts) ;

  extern fact_db *exec_current_fact_db ;

  extern storeRepP collect_store(storeRepP &sp, fact_db &facts) ;
  
  inline storeRepP collect_store(storeRepP &sp) 
    { return collect_store(sp,*exec_current_fact_db) ;}
  
  inline entitySet collect_entitySet(entitySet e)
    { return collect_entitySet(e,*exec_current_fact_db) ; }
  
  inline storeRepP distribute_store(storeRepP &sp)
    { return distribute_store(sp,*exec_current_fact_db) ; }
  
  void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn); 
  entitySet all_collect_entitySet(const entitySet &e) ;
  std::vector<entitySet> all_collect_vectors(entitySet &e) ;
  std::vector<int> all_collect_sizes(int size) ;
  int GLOBAL_OR(int b) ;
  int GLOBAL_AND(int b) ;
  int GLOBAL_MAX(int b) ;
  int GLOBAL_MIN(int b) ;
  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) ;

  void read_multi_vector_int(hid_t group_id, const char* name, int dim,  std::vector<int>& vint) ;
  void read_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size) ;
  void write_vector_int(hid_t group_id, const char* name, std::vector<int>& vint) ;
}
#endif
 

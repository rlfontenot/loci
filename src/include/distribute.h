#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H
#include <typeinfo>
#include <vector>
#include <store.h>
#include <Map.h>
#include <DMap.h>
#include <DMultiMap.h>
#include <constraint.h>
#include <fact_db.h>

#define LOCI_HAS_READ_WRITE_PARTITION 1
 
namespace Loci {

  extern ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  extern double barrier_time ;
  
  void Init(int* argc, char*** argv) ;
  void Finalize() ;
  
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn, int num_partitions = 0 ) ;
  void categories(fact_db &facts,std::vector<interval> &pvec) ;
  
  void get_mappings(rule_db &rdb, fact_db &facts, std::set<std::vector<variableSet> > &maps) ;
  
  
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  
  
  std::vector<entitySet> generate_distribution(fact_db &facts, rule_db &rdb, int num_partitions = 0) ;

  std::vector<entitySet> read_partition(const char *fname,int num_partitions) ;
  void write_partition(const char *fname, const std::vector<entitySet> &ptn) ;

  void distribute_facts(std::vector<entitySet> &ptn, fact_db &facts,
                        rule_db &rdb) ;
  
  entitySet fill_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> fill_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  entitySet send_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> send_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  void print_global(entitySet e, fact_db &facts) ;
  entitySet collect_entitySet(entitySet e, fact_db &facts) ;
  storeRepP collect_store(storeRepP &sp, fact_db &facts) ;
  storeRepP distribute_store(storeRepP &sp, fact_db &facts) ;

  extern fact_db *exec_current_fact_db ;

  extern storeRepP collect_store(storeRepP &sp, fact_db &facts) ;

  inline storeRepP collect_store(storeRepP &sp) 
  { return collect_store(sp,*exec_current_fact_db) ;}

  inline entitySet collect_entitySet(entitySet e)
  { return collect_entitySet(e,*exec_current_fact_db) ; }

  inline storeRepP distribute_store(storeRepP &sp)
  { return distribute_store(sp,*exec_current_fact_db) ; }

   void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) ;
 
 void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) ;
  
 void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) ;
 
 void distributed_inverseMap(dmultiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts); 
 
}


#endif

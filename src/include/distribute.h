#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H
#include <typeinfo>
#include <vector>
#include <store.h>
#include <Map.h>
#include <constraint.h>
#include <fact_db.h>

 
namespace Loci {

  extern ofstream debugout ;
  
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  extern bool EXIT ;
  extern bool READ_DIST ;
  extern int NUM_PARTITIONS ;
  void Init(int* argc, char*** argv) ;
  void Finalize() ;
  
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn ) ;
  void categories(fact_db &facts,std::vector<interval> &pvec) ;
  
  void get_mappings(rule_db &rdb, fact_db &facts, std::set<std::vector<variableSet> > &maps) ;
  
  
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  
  
  std::vector<entitySet> generate_distribution(fact_db &facts, rule_db &rdb) ;
  
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
}


#endif

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H
#include <typeinfo.h>
#include <vector>
#include <store.h>
#include <Map.h>
#include <constraint.h>
#include <fact_db.h>


namespace Loci {
 
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  
  void Init(int* argc, char*** argv) ;
  
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn, store<int> &partition ) ;
  void categories(fact_db &facts,std::vector<interval> &pvec) ;
  
  void get_mappings(rule_db &rdb, std::set<std::vector<variableSet> > &maps) ;
  
  
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  
  
  void generate_distribution(fact_db &facts, rule_db &rdb,std::vector<std::vector<entitySet> > &get_entities ) ;
  
  void distribute_facts(std::vector<std::vector<entitySet> > &get_entities, fact_db &facts)   ;
  entitySet fill_entitySet(entitySet& e, fact_db &facts) ;
  entitySet send_entitySet(entitySet& e, fact_db &facts) ;
  void print_global(entitySet e, fact_db &facts) ;
}


#endif

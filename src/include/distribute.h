#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H
#include <typeinfo.h>
#include <vector>
#include <store.h>
#include <Map.h>
#include <constraint.h>
#include <fact_db.h>


namespace Loci {
 
  struct dist_facts {
    int myid ;
    store<int> isDistributed ;
    constraint my_entities ;
    Map g2l ;
    Map l2g ;
    entitySet send_neighbour ;
    entitySet recv_neighbour ;
    store<entitySet> send_entities ;
    store<entitySet> recv_entities ;
  } ;
  
  class distribute_info {
  public:
    distribute_info();
    dist_facts* get_dist_facts() ; 
    void set_dist_facts(dist_facts dist) ;
    void operator = (distribute_info&) ;
  private:
    dist_facts distributed_facts ;
    
  } ;
  
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
  
  distribute_info distribute_facts(std::vector<std::vector<entitySet> > &get_entities, fact_db &facts)   ;
  entitySet fill_entitySet(entitySet& e, distribute_info& dist, fact_db &facts) ;
  
}


#endif

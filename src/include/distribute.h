#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H
#include <typeinfo.h>
#include <vector>
//#include <ostream>
#include <store.h>
#include <Map.h>
#include <constraint.h>
#include <fact_db.h>


namespace Loci {
  
  struct dist_facts {
    store<int> isDistributed ;
    constraint my_entities ;
    Map g2l ;
    Map l2g ;
    store<entitySet> send_neighbour ;
    store<entitySet> recv_neighbour ;
  } ;
  
 
  class distribute_info {
  public:
    distribute_info() ;
    distribute_info(int);
    dist_facts& get_dist_facts() ; 
    void set_dist_facts(int, store<int>, constraint, Map, Map, store<entitySet>, store<entitySet>) ;
  private:
    dist_facts distributed_facts ;
    
    } ;
    
  void Init(int argc, char** argv, int& num_procs, int& myid) ;
 
  void metis_facts(fact_db &facts,int num_partitions,
		   std::vector<entitySet> &ptn, store<int> &partition ) ;
  void categories(fact_db &facts,std::vector<interval> &pvec) ;
  
  void get_mappings(rule_db &rdb, std::set<std::vector<variableSet> > &maps) ;
  
  
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  
  
  void generate_distribution(fact_db &facts, rule_db &rdb, int num_procs,std::vector<std::vector<entitySet> > &get_entities ) ;
  
  void distribute_facts(std::vector<std::vector<entitySet> > &get_entities, fact_db &facts, int myid)   ;
  
  
}


#endif

#include <distribute.h>

namespace Loci {
  void get_clone(fact_db &facts, const rule_db &rdb) ;
  void categories(fact_db &facts,std::vector<entitySet> &pvec) ;
  std::vector<entitySet> modified_categories(fact_db &facts, std::map<variable, entitySet> &vm, std::vector<interval> &pvec) ;
  void get_mappings(const rule_db &rdb, fact_db &facts, std::set<std::vector<variableSet> > &maps) ;
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  entitySet dist_expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  std::vector<entitySet> generate_scalable_distribution(fact_db &facts, rule_db &rdb, int num_partitions = 0) ;

  std::vector<entitySet> read_partition(const char *fname,int num_partitions) ;
  void write_partition(const char *fname, const std::vector<entitySet> &ptn) ;

  entitySet fill_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> fill_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  entitySet send_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> send_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  void print_global(entitySet e, fact_db &facts) ;
  Map distribute_whole_map(Map &m) ;
}



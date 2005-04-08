#ifndef DIST_TOOLS_H
#define DIST_TOOLS_H

#include <distribute.h>
#include <distribute_io.h>
#include <distribute_container.h>

namespace Loci {

  extern bool use_dynamic_scheduling ;
  void get_clone(fact_db &facts, const rule_db &rdb) ;
  void categories(fact_db &facts,std::vector<entitySet> &pvec) ;
  void get_mappings(const rule_db &rdb, fact_db &facts, std::set<std::vector<variableSet> > &maps, unsigned int rule_part = 0) ;
  
  entitySet expand_map(entitySet domain, fact_db &facts,
		       const std::set<std::vector<variableSet> > &maps) ;
  entitySet dist_expand_map(entitySet domain, fact_db &facts,
			    const std::set<std::vector<variableSet> > &maps) ;
  entitySet fill_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> fill_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  entitySet send_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> send_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn, int num_partitions = 0 );
  entitySet context_for_map_output(entitySet domain, fact_db &facts,
				   const std::set<std::vector<variableSet> > &maps);
  entitySet dist_reverse_expand_map(fact_db &facts,
                                    const std::set<std::vector<variableSet> > &maps);
  entitySet dist_special_expand_map(entitySet domain, fact_db &facts,
				    const std::set<std::vector<variableSet> > &maps);
  // function that restores the fact_db back to its global numbering
  void restore_global_facts(fact_db& facts) ;
  
}

#endif

#ifndef DIST_TOOLS_H
#define DIST_TOOLS_H

#include <distribute.h>
#include <distribute_io.h>
#include <distribute_container.h>
#include <vector>

namespace Loci {

  extern bool use_dynamic_scheduling ;
  void get_clone(fact_db &facts, const rule_db &rdb) ;
  void categories(fact_db &facts,std::vector<entitySet> &pvec) ;
  entitySet dist_collect_entitySet(entitySet inSet, const std::vector<entitySet> &ptn) ;
  entitySet dist_expand_entitySet(entitySet inSet, entitySet copy,
                                  const std::vector<entitySet> &ptn) ;
  void get_mappings(const rule_db &rdb, fact_db &facts, std::set<std::vector<variableSet> > &maps, unsigned int rule_part = 0) ;
  
  entitySet dist_expand_map(entitySet domain, fact_db &facts,
			    const std::set<std::vector<variableSet> > &maps) ;
  entitySet fill_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> fill_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  entitySet send_entitySet(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> send_entitySetv(const entitySet& e, fact_db &facts) ;
  std::vector<entitySet> send_entitySet(const std::vector<entitySet>& e,
                                        fact_db &facts) ;
  entitySet context_for_map_output(entitySet domain, fact_db &facts,
				   const std::set<std::vector<variableSet> > &maps);
  entitySet dist_special_expand_map(entitySet domain, fact_db &facts,
				    const std::set<std::vector<variableSet> > &maps);
  std::set<std::vector<variableSet> > classify_moderate_maps(fact_db &facts, const std::set<std::vector<variableSet> > &maps);
  // function that restores the fact_db back to its global numbering
  void restore_global_facts(fact_db& facts) ;
  entitySet findBoundingSet(entitySet in) ;
  rule_db replace_map_constraints(fact_db& facts, const rule_db& rdb) ;
}

#endif

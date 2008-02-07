//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
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

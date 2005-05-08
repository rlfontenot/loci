#include "sched_tools.h"
#include "loci_globs.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;
using std::ostringstream ;

namespace Loci {
// Create a sequential ordering of rules based on the concurrent dag schedule
  void extract_rule_sequence(vector<rule> &rule_seq,
                             const vector<digraph::vertexSet> &v) {
    vector<digraph::vertexSet>::const_iterator i ;
    for(i=v.begin();i!=v.end();++i) {
      ruleSet rules = extract_rules(*i) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        rule_seq.push_back(*ri) ;
    }
  }
  
  rule make_rename_rule(variable new_name, variable old_name) {
    ostringstream oss ;
    oss << "source(" << old_name << "),target(" << new_name
        << "),qualifier(rename)"  ;
    return rule(oss.str()) ;
  }

  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts, sched_db &scheds) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      if(duplicate_work) {
	entitySet temp = scheds.variable_existence(*vi) ;
	if(facts.isDistributed()) {
	  //Some constraints are applied on map variables
	  //For that case, the existence of a map is actually its domain on a processor
	  //If the map is actually not used in the rule other than constraints,
	  //then that map may not have expanded enough to include necessrary 
	  //domain on a processor.  We need to communicate to collect existence.
	  if((facts.get_variable(*vi))->RepType() == MAP) {
	    temp = collect_entitySet(temp);
	    temp = all_collect_entitySet(temp);
	    Map l2g = facts.get_variable("l2g");
	    temp = l2g.preimage(temp).first;
	  }
	}
	sources &= temp;
      }
      else
	sources &= scheds.variable_existence(*vi) ;
    }
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
	working &= scheds.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }


  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute, sched_db &scheds) {
    vector<variableSet>::const_iterator mi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      if(mi->size() == 1) {
        variable v = *(mi->begin()) ;
        FATAL(!scheds.is_a_Map(v)) ;
        compute = scheds.image(v,compute) ;
      } else {
        variableSet::const_iterator vi ;
        entitySet images ;
        for(vi=mi->begin();vi!=mi->end();++vi) {
          variable v = *vi ;
          FATAL(!scheds.is_a_Map(v)) ;
          images |= scheds.image(v,compute) ;
        }
        compute = images ;
      }
    }
    return compute ;
  }

}

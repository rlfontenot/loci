#include "sched_tools.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

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
        
  rule make_super_rule(variableSet sources, variableSet targets,
                       variable cond) {
    FATAL(targets == EMPTY) ;
    static int super_node_number = 0 ;
    ostringstream oss ;
    oss << "source("<<sources << "),target(" << targets << ")," ;
    if(cond != variable()) 
      oss<< "conditional(" << cond << ")," ;
    oss << "qualifier(SN" << super_node_number++ << ")" ;
   
    return rule(oss.str()) ;
  }

  rule make_rename_rule(variable new_name, variable old_name) {
    ostringstream oss ;
    oss << "source(" << old_name << "),target(" << new_name
        << "),qualifier(rename)"  ;
    return rule(oss.str()) ;
  }

  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      sources &= facts.variable_existence(*vi) ;
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working &= facts.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }


  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute) {
    vector<variableSet>::const_iterator mi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      if(mi->size() == 1) {
        variable v = *(mi->begin()) ;
        FATAL(!facts.is_a_Map(v)) ;
        compute = facts.image(v,compute) ;
      } else {
        variableSet::const_iterator vi ;
        entitySet images ;
        for(vi=mi->begin();vi!=mi->end();++vi) {
          variable v = *vi ;
          FATAL(!facts.is_a_Map(v)) ;
          images |= facts.image(v,compute) ;
        }
        compute = images ;
      }
    }
    return compute ;
  }

  digraph::vertexSet visit_vertices(digraph dg,digraph::vertexSet begin) {

    digraph::vertexSet visit = begin ;
    digraph::vertexSet visited ;
    digraph::vertexSet::const_iterator ni ;

    // keep visiting vertices until no new vertices are found
    while(visit != EMPTY) {
      digraph::vertexSet newvertices ;
      // visit all the vertices that this graph leads to
      for(ni=visit.begin();ni!=visit.end();++ni)
        newvertices += dg[*ni] ;
      // update visit, but don't re-visit a vertex
      visit = newvertices - visited ;
      visited = visited + newvertices ;
    }
    return visited ;
  }

  digraph::vertexSet visit_vertices_exclusive(digraph dg, digraph::vertexSet begin) {
    digraph dgt = dg.transpose() ;
    digraph::vertexSet visited, visit ;

    visit = begin ;
    visited = visit ;

    // keep visiting vertices until no new vertices are found
    while(visit != EMPTY) {
      digraph::vertexSet::const_iterator ni,nii,no ;
      digraph::vertexSet newvertices, not_visited = ~visited ;
      // visit all new vertices and check to see if they contain paths to
      // new candidate vertices
      for(ni=visit.begin();ni!=visit.end();++ni) 
        for(nii=dg[*ni].begin();nii!=dg[*ni].end();++nii) { // loop over edges
          // if this edge leads to a vertex that can only be reached through
          // the visited set, then it is a new vertex to be visited
          digraph::vertexSet out_bound_vertices = dgt[*nii] & not_visited ;

          // Check to see if out of bound vertices loop exclusively back to this
          // vertex, if so then add it to the list of new vertices.
          // This is a little bit of a hack since the graphs analyzed may
          // contain loops no larger than 2 edges and 2 vertexes.
          bool flg = true ;
          digraph::vertexSet local_not_visit = ~(visited + interval(*nii,*nii));
          for(no=out_bound_vertices.begin();no!=out_bound_vertices.end();++no) {
            flg = flg && (dgt[*no] & local_not_visit) == EMPTY ;
          }
          if(flg) 
            newvertices += *nii ;
        }
      // next time we will visit the vertices found in this iteration
      // update our list of visited vertices.
      visit = newvertices - visited ;
      visited += newvertices ;
    }
    return visited ;
  }



}

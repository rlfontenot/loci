#include "comp_tools.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <distribute.h>
namespace Loci {
    // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet> schedule_dag(const digraph &g,
                                          digraph::vertexSet start_vertices,
                                          digraph::vertexSet only_vertices) {

    digraph gt = g.transpose() ;

    vector<digraph::vertexSet> schedule ; 
    // First schedule any vertices that have no edges leading into them and have
    // not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    if(working != EMPTY)
      schedule.push_back(working) ;
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = start_vertices + working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni)
        new_vertices += g[*ni] ;
      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      // We only schedule vertices that are also in the only_vertices set
      working = new_vertices & only_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all vertices leading
      // to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni) 
        if((gt[*ni] & visited_vertices) == gt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // and these new vertices to the schedule
      if(new_vertices != EMPTY)
        schedule.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
    return schedule ;
  }

  
  void existential_rule_analysis(rule r, fact_db &facts) {
    
    FATAL(r.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet::const_iterator ei ;
    
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = new fact_db::distribute_info ;
      d = facts.get_distribute_info() ;
      constraint my_entities ;
      my_entities = facts.get_variable("my_entities") ;
      
      sources &= my_entities ;
      constraints &= my_entities ;
      
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;
      // First we compute the intersection of all the sources (inputs) to
      // the rule
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) 
	sources &= vmap_source_exist(*si,facts) ;
      // We compute the intersection of all the constraints
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
	constraints &= vmap_source_exist(*si,facts) ;
      
      // We do a sanity check to make sure that the rule is able to
      // compute all entities in the constraint
      if(rinfo.constraints.begin() != rinfo.constraints.end())
	if((sources & constraints) != constraints) {
	  cerr << "Warning  " << "in processor   " << d->myid << "   , rule " << r <<
	    " cannot supply all entities of constraint" << endl ;
	  cerr << "constraints = " << constraints ;
	  cerr << "sources & constraints = " << (sources & constraints) << endl ;
	  // exit(-1) ;
	}
      sources &= constraints ;
      
      // Context is the set of entities that the rule will loop over
      // to compute its values.  If we didn't have maps in the output
      // we would be pretty much finished here we could just assign
      // the context to the list of output variables in the existential
      // analysis.
      
      entitySet context = sources & constraints ;
      
      // Instead we loop over all of the targets (outputs) of the rule
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	// Here we transform the context through any mappings used in
	// the output.
	entitySet targets ;
	if((si->mapping).size() == 0) {
	  cout << *si << "  case with no mapping in output " << endl ;
	  context = fill_entitySet(context, facts) ;
	  targets = vmap_target_exist(*si,facts,context) ;
	}
	else {
	  cout << *si <<  "  case with mapping in output " << endl ;
	  targets = vmap_target_exist(*si,facts,context) ;
	  targets = fill_entitySet(targets, facts) ;
	  targets = send_entitySet(targets, facts) ;
	}
	
	// now targets contains the entities for which values will be
	// generated by applying this rule over context.
	const variableSet &tvars = si->var ;
	variableSet::const_iterator vi ;
	// Now we loop over the variables that the above mapping applies to
	// for example in a->(b,c) we will be looping over b and c.
	// For each of these values we update the existential info
	// concerning entities which can be assigned these values.
	for(vi=tvars.begin();vi!=tvars.end();++vi) {
	  facts.set_existential_info(*vi,r,targets) ;
	  //cout << "processor  =  " << d->myid << "  vi =   " << *vi << "\t targets =   " 
	  //<< targets << endl ;  
	}
      }
    }
    else {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
	sources &= vmap_source_exist(*si,facts) ;
      } 
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
	constraints &= vmap_source_exist(*si,facts) ;
      if(rinfo.constraints.begin() != rinfo.constraints.end())
	if((sources & constraints) != constraints) {
	  cerr << "Warning, rule " << r <<
	    " cannot supply all entities of constraint" << endl ;
	  cerr << "constraints = " << constraints ;
	  cerr << "sources & constraints = " << (sources & constraints) << endl ;
	  // exit(-1) ;
	}
      sources &= constraints ;
      
      entitySet context = sources & constraints ;
      
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	entitySet targets = vmap_target_exist(*si,facts,context) ;
	const variableSet &tvars = si->var ;
	variableSet::const_iterator vi ;
	for(vi=tvars.begin();vi!=tvars.end();++vi) {
	  facts.set_existential_info(*vi,r,targets) ;
	}
      }
    }
#ifdef VERBOSE
    cout << "rule " << r << " generating " << targets << endl ;
#endif
    
  }
  
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts) {
    // Here we will compute the context implied by a particular target
    // mapping
    variableSet::const_iterator vi ;
    entitySet targets ;
    // First we get the queries for all the variables that the mapping
    // is applied to.
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      // The variable should be in tvarmap, but we will report an error
      // if something fishy happens
      FATAL(tvarmap.find(*vi) == tvarmap.end()) ;
      // Get the requests for variable *vi and union it with the target set.
      targets |= tvarmap.find(*vi)->second ;
    }
    // Here we do a hack to make sure that if there are multiple
    // variable that a map applies to, then they will request their
    // union.  (e.g. for a->(b,c) we make sure that b and c both have
    // the same requests.
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,targets) ;

    // Now we are applying the mapping that is applied to the target
    // variables.  We do this by finding the preimage of each map.
    // (We use the union preimage since any value that is touched
    // will need to be computed.  The union primage is in the second
    // element of the pair that preimage returns. 
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      // working is the entityset the becomes the union of all primages
      // on this level
      entitySet working = EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.preimage(*vi,targets).second ;
      }
      // Now we have evaluated this map, we move targets to this level
      targets = working ;
    }
    // When we are finished, we have followed the maps back from the targets to
    // their root.  We now have the set of entities that will be in the context
    // of the rule that will be used to satisfy this set of requests.
    return targets ;
  }

  entitySet vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                            entitySet context) {
    // this routine computes the set of entities that a source mapping will
    // imply.  It does this by following the images of the mapping.
    // The resulting entitySet contains all entities that will be accessed
    // when a loop over context is executed.
    entitySet compute = context ;
    vector<variableSet>::const_iterator mi ;
    variableSet::const_iterator vi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      entitySet working ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.image(*vi,compute) ;
      }
      compute = working ;
    }
    return compute ;
  }

  entitySet process_rule_requests(rule r, fact_db &facts) {

    // Internal rules should be handling the appropriate rule requests via
    // their associated compiler.
    FATAL(r.type() == rule::INTERNAL) ;


    // First we get the target variables of this rule ;
    variableSet targets = r.targets() ;
    // We will be iterating over the target variables so we need an iterator
    variableSet::const_iterator vi ;  

    // The vdefmap data structure is a map from variables to entitySets.
    // We use the tvarmap to record the requests for target variables
    // Here we are filling in the requests.
    vdefmap tvarmap ;
    // Loop over target variables and get requests from fact database
    for(vi=targets.begin();vi!=targets.end();++vi) {
      // This is a hack for the special case of a rule with OUTPUT
      // as a target.  In that case we will request OUTPUT for
      // all entities that exist.  So we add a request for OUTPUT
      // to the fact database
      if(vi->get_info().name == string("OUTPUT")) 
        facts.variable_request(*vi,facts.variable_existence(*vi)) ;

      // Now fill tvarmap with the requested values for variable *vi
      tvarmap[*vi] = facts.get_variable_request(r,*vi) ;
    }

    
    const rule_impl::info &rinfo = r.get_info().desc ;

    // Here we compute the context of the rule.  This is the union of all of
    // the requests for the variables that this rule produces
    set<vmap_info>::const_iterator si ;
    entitySet context,isect = ~EMPTY ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      // Transform the variable requests using the mapping constructs
      // in *si
      entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
      // The context is the union
      context |= tmp ;
      isect &= tmp ;
    }

    // If the interstection and the union are not equal, then we are in
    // danger of not properly allocating variables for computations.  It is
    // an optimization to check this.  For the distributed memory version it
    // may be useful to always do this if there is a mapping in the
    // targets of the rule.

    if(isect != context) {
      entitySet working = context ;

      vector<variableSet>::const_reverse_iterator mi ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        for(mi=si->mapping.rbegin();mi!=si->mapping.rend();++mi) {
          entitySet tmp ;
          for(vi=mi->begin();vi!=mi->end();++vi)
            tmp |= facts.image(*vi,working) ;
          working = tmp ;
        }
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          facts.variable_request(*vi,working) ;
        }
      }
    }


    // Loop over all sources for this rule and pass on the requests.
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      // First map the context through source mappings
      entitySet requests = vmap_source_requests(*si,facts,context) ;
      // Now we have the actual requests we are making of other rules
      // so we can tell the fact database that we are now requesting
      // these values.
      for(vi=si->var.begin();vi!=si->var.end();++vi)
        facts.variable_request(*vi,requests) ;
    }

    // We also need to pass the requests on to any conditional variables
    // this rule may have.
    for(vi=rinfo.conditionals.begin();vi!=rinfo.conditionals.end();++vi)
      facts.variable_request(*vi,context) ;

#ifdef VERBOSE
    cout << "rule " << r << " computes over " << context << endl ;
#endif
    return context ;
  }



  vector<entitySet> partition_set(const entitySet &s,int nthreads) {
    const int min = s.Min() ;
    const int max = s.Max() ;
    const int psize = s.size() ;

    const int div = psize/nthreads ;
    int rem = psize%nthreads ;
    vector<entitySet> partition ;
    
    for(int i=min;i<=max;) {
      int inval = div + ((rem>0)?1:0) ;
      rem-- ;
      entitySet sp = s & interval(i,i+inval-1) ;
      i+=inval ;
      while(sp.size() < inval && i<=max) {
        entitySet remain = s & interval(i,max) ;
        i=remain.Min() ;
        int remain_ival = inval - sp.size() ;
        sp += s & interval(i,i+remain_ival-1) ;
        i+=remain_ival ;
      }
      WARN(sp.size() > inval) ;
      partition.push_back(sp) ;
    }
    return partition ;
  }

  void parallel_schedule(execute_par *ep,const entitySet &exec_set,
                         const rule &impl, fact_db &facts) {
    vector<entitySet> par_set = partition_set(exec_set,num_threads) ;

    for(vector<entitySet>::const_iterator
          i=par_set.begin();i!=par_set.end();++i) {
      executeP execrule = new execute_rule(impl,sequence(*i),facts) ;
      ep->append_list(execrule) ;
    }
  }

}

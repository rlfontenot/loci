#include "comp_tools.h"

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


  void existential_rule_analysis(rule f, fact_db &facts) {

    FATAL(f.type() == rule::INTERNAL) ;

    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet targets = EMPTY ;
    const rule_impl::info &finfo = f.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
      sources &= vmap_source_exist(*si,facts) ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
      constraints &= vmap_source_exist(*si,facts) ;
    if(finfo.constraints.begin() != finfo.constraints.end())
      if((sources & constraints) != constraints) {
        cerr << "Warning, rule " << f <<
          " cannot supply all entities of constraint" << endl ;
        cerr << "constraints = " << constraints ;
        cerr << "sources & constraints = " << (sources & constraints) << endl ;
        //      exit(-1) ;
      }
    sources &= constraints ;
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      targets = vmap_target_exist(*si,facts,sources) ;
      const variableSet &tvars = si->var ;
      variableSet::const_iterator vi ;
      for(vi=tvars.begin();vi!=tvars.end();++vi)
        facts.set_existential_info(*vi,f,targets) ;
    }
#ifdef VERBOSE
    cout << "rule " << f << " generating " << targets << endl ;
#endif
  }

  
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts) {
    variableSet::const_iterator vi ;
    entitySet targets ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      FATAL(tvarmap.find(*vi) == tvarmap.end()) ;
      targets |= tvarmap.find(*vi)->second ;
    }
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,targets) ;
  
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.preimage(*vi,targets).second ;
      }
      targets = working ;
    }
    return targets ;
  }

  void vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                            entitySet compute) {
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
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,compute) ;
  }

  entitySet process_rule_requests(rule f, fact_db &facts) {

    FATAL(f.type() == rule::INTERNAL) ;

    variableSet targets = f.targets() ;
  
    variableSet::const_iterator vi ;
    vdefmap tvarmap ;
    for(vi=targets.begin();vi!=targets.end();++vi) {
      if(vi->get_info().name == string("OUTPUT")) 
        facts.variable_request(*vi,facts.variable_existence(*vi)) ;
      tvarmap[*vi] = facts.get_variable_request(f,*vi) ;
    }

    const rule_impl::info &finfo = f.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    entitySet compute,isect = ~EMPTY ;
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
      compute |= tmp ;
      isect &= tmp ;
    }

    if(isect != compute) {
      entitySet working = compute ;
      vector<variableSet>::const_reverse_iterator mi ;
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
        for(mi=si->mapping.rbegin();mi!=si->mapping.rend();++mi) {
          entitySet tmp ;
          for(vi=mi->begin();vi!=mi->end();++vi)
            tmp |= facts.image(*vi,working) ;
          working = tmp ;
        }
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          facts.variable_request(*vi,working) ;
          //#define UNDERSTAND_THIS
#ifdef UNDERSTAND_THIS
          WARN(facts.get_variable_request(f,*vi) != working) ;
          if(facts.get_variable_request(f,*vi) != working) {
            cerr << "f = " << f << endl ;
            cerr << "*vi = " << *vi << endl ;
            cerr << "facts.get_variable_request(f,*vi) = "
                 << facts.get_variable_request(f,*vi)
                 << ", working = " << working << endl ;
          }
#endif
        }
      }
    }

    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
      vmap_source_requests(*si,facts,compute) ;

    for(vi=finfo.conditionals.begin();vi!=finfo.conditionals.end();++vi)
      facts.variable_request(*vi,compute) ;

#ifdef VERBOSE
    cout << "rule " << f << " computes over " << compute << endl ;
#endif
    return compute ;
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

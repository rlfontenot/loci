#include "comp_tools.h"
#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <distribute.h>
//#define VERBOSE1
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
      constraint my_entities ;
      my_entities = facts.get_variable("my_entities") ;
      sources &= my_entities ;
      constraints &= my_entities ;
    }
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
    // union.  e.g. for a->(b,c) we make sure that b and c both have
    // the same requests.
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      facts.variable_request(*vi,targets) ;
    
    // Now we are applying the mapping that is applied to the target
    // variables. We do this by finding the preimage of each map.
    // We use the union preimage since any value that is touched
    // will need to be computed.  The union preimage is in the second
    // element of the pair that preimage returns. 
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      // working is the entityset that becomes the union of all preimages
      // on this level
      entitySet working = EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.preimage(*vi,targets).second ;
	//cout << "mi = " << *mi << "   vi =  " << *vi  <<"    working = " << working << endl ;
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
     entitySet::const_iterator ei, ti ;
     // The vdefmap data structure is a map from variables to entitySets.
     // We use the tvarmap to record the requests for target variables
     // Here we are filling in the requests.
     vdefmap tvarmap ;
     // Loop over target variables and get requests from fact database
     
     // Here we compute the context of the rule.  This is the union of all of
     // the requests for the variables that this rule produces
     set<vmap_info>::const_iterator si ;
     entitySet context,isect = ~EMPTY ;
     if(facts.isDistributed()) {
       constraint my_entities ;
       my_entities = facts.get_variable("my_entities") ;
       
       for(vi=targets.begin();vi!=targets.end();++vi) {
	 // This is a hack for the special case of a rule with OUTPUT
	 // as a target.  In that case we will request OUTPUT for
	 // all entities that exist.  So we add a request for OUTPUT
	 // to the fact database
	
	 if(vi->get_info().name == string("OUTPUT")) 
	   facts.variable_request(*vi,facts.variable_existence(*vi)) ;
	 
	 // Now fill tvarmap with the requested values for variable *vi
	 tvarmap[*vi] = facts.get_variable_request(r,*vi) ;
	 //cout << d->myid << "    variable  =  "<< *vi << "   tvarmap  =  " <<
	 //tvarmap[*vi] << endl ;
       }
       const rule_impl::info &rinfo = r.get_info().desc ;
       
       for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
         // Transform the variable requests using the mapping constructs
         // in *si
         entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
         //The context is the union
         context |= tmp ;
         isect &= tmp ;
         //cout <<d->myid <<"      si =  " << *si <<  "   context = " << context << endl ;
       }

       // If the interstection and the union are not equal, then we are in
       // danger of not properly allocating variables for computations.  It is
       // an optimization to check this. For the distributed memory version it
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

       entitySet working = context ;      

      // Loop over all sources for this rule and pass on the requests.

       for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        // First map the context through source mappings
         entitySet requests;
         requests = vmap_source_requests(*si,facts,context) ;
         //cout <<d->myid <<  "   *si  =  "  << *si << "   requests  =  " << re
         //quests << endl ;
         entitySet var ;
         context &= my_entities ;
        
         // Now we have the actual requests we are making of other rules
         // so we can tell the fact database that we are now requesting
         // these values.
	
         for(vi=si->var.begin();vi!=si->var.end();++vi)
           facts.variable_request(*vi,requests) ;
	
         // We also need to pass the requests on to any conditional variables
         // this rule may have.
	
         for(vi=rinfo.conditionals.begin();vi!=rinfo.conditionals.end();++vi) 
           facts.variable_request(*vi,context) ;
       }
     }  else {
       for(vi=targets.begin();vi!=targets.end();++vi) {
        if(vi->get_info().name == string("OUTPUT")) 
          facts.variable_request(*vi,facts.variable_existence(*vi)) ;
        tvarmap[*vi] = facts.get_variable_request(r,*vi) ;
       }
      const rule_impl::info &rinfo = r.get_info().desc ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
        context |= tmp ;
        isect &= tmp ;
      } 
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
        for(vi=si->var.begin();vi!=si->var.end();++vi)
          facts.variable_request(*vi,requests) ;
      }
      
      // We also need to pass the requests on to any conditional variables
      // this rule may have.
      for(vi=rinfo.conditionals.begin();vi!=rinfo.conditionals.end();++vi)
        facts.variable_request(*vi,context) ;
     }
#ifdef VERBOSE
     cout << "rule " << r << " computes over " << context << endl ;
#endif
    return context ;
  }

  std::list<comm_info> put_postcomm_info(std::map<variable, ruleSet> barrier_info, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int k = 0, recv_flag = 0, recv_count = 0 ;
    entitySet::const_iterator ei, ti ;
    entitySet context,isect = ~EMPTY, temp ;
    std::list<comm_info> clist ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    constraint my_entities ;
    Map l2g ;
    my_entities = facts.get_variable("my_entities") ;
    l2g = facts.get_variable("l2g") ;
    int **send_buffer, **recv_buffer ;
    int *recv_size, *send_size ;
    std::map<variable, ruleSet>::iterator mi ;
    entitySet requests ;
    for(mi = barrier_info.begin() ; mi != barrier_info.end(); ++mi) {
      variable v = mi->first ;
      requests = facts.get_variable_requests(v) ;
      comm_info ci ;
      ci.v = v ;
      entitySet tempset ;
      std::vector<proc_details> send_pvec ;
      std::vector<proc_details> recv_pvec ;
      recv_buffer = new int*[d->send_neighbour.size()] ;
      recv_size = new int[d->send_neighbour.size()] ;
      
      send_buffer = new int*[d->recv_neighbour.size()] ;
      send_size = new int[d->recv_neighbour.size()] ;
      recv_count = 0 ;
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	recv_size[k] = (d->send_entities[*ei]).size() ;
	recv_buffer[k] = new int[recv_size[k]] ;
	recv_count++ ;
	k++ ;
      }
      recv_request = (MPI_Request *) malloc(recv_count * sizeof(MPI_Request) ) ;
      status = (MPI_Status *) malloc(recv_count * sizeof(MPI_Status) ) ;
      k = 0 ;
      
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	MPI_Irecv(&recv_buffer[k][0], recv_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD, 
		  &recv_request[k] ) ;  
	recv_flag = 1 ;
	k++ ;
      }
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	temp = EMPTY ;
	for(ti = requests.begin(); ti != requests.end(); ++ti)
	  temp += l2g[*ti] & d->recv_entities[*ei] ;
	send_size[k] = temp.size()  ;
	send_buffer[k] = new int[send_size[k]] ;
	k++ ;
      }
      Map id ;
      id.allocate(interval(0,d->send_neighbour.size()-1)) ;
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	id[k] = *ei ;
	k++ ;
      }
      k = 0 ;
      
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	proc_details pd ;
	pd.processor = *ei ;
	temp = EMPTY ;
	for(ti = requests.begin(); ti != requests.end(); ++ti)
	  temp += l2g[*ti] & d->recv_entities[*ei] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[k][j] = *ti ;
	  pd.recv_set += d->g2l[*ti] ;
	  ++j ;
	}
	//cout << d->myid << "post comm variable  = " << *vi << "    receiving from   " << *ei << "   " << pd.recv_set << endl ;
	pd.send_set = EMPTY ;
	MPI_Send(&send_buffer[k][0], send_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD) ;
	k++ ;
	recv_pvec.push_back(pd) ;
      }
      if(recv_flag)
	MPI_Waitall(recv_count, recv_request, status) ;
      
      for(k = 0; k < recv_count; ++k)
	MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
      
      for(k = 0; k < d->send_neighbour.size(); ++k) {      
	if((recv_flag == 1) && (recv_size[k] > 0)) {
	  proc_details pd ;
	  pd.processor = id[k] ;
	  pd.send_set = EMPTY ;
	  pd.recv_set = EMPTY ;
	  for(int i = 0; i < recv_size[k]; ++i) {
	    pd.send_set += d->g2l[recv_buffer[k][i]] ;
	  }
	  //cout << d->myid << "post comm variable  = " << *vi << "    sending to   " << id[k] << "   " << pd.send_set << endl ;
	  send_pvec.push_back(pd) ;
	}
      }
      
      ci.send_info = send_pvec ;
      ci.recv_info = recv_pvec ;
      clist.push_back(ci) ;
      
      delete [] send_size ;
      delete [] recv_size ;
      delete [] send_buffer ;
      delete [] recv_buffer ; 
      
    }
    return clist ;
  } 
  
  std::list<comm_info> put_precomm_info(std::map<variable, ruleSet > barrier_info, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int k = 0, recv_flag = 0, recv_count = 0 ;
    variableSet::const_iterator vi ;  
    entitySet::const_iterator ei, ti ;
    entitySet context, isect = ~EMPTY, temp ;
    std::list<comm_info> plist ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    constraint my_entities ;
    Map l2g ; 
    my_entities = facts.get_variable("my_entities") ;
    l2g = facts.get_variable("l2g") ;
    int **send_buffer, **recv_buffer ;
    int *recv_size, *send_size ;
    std::map<variable, ruleSet>::iterator mi ;
    entitySet requests ;
    for(mi = barrier_info.begin(); mi != barrier_info.end(); ++mi) {
      variable v = mi->first ;
      ruleSet rs = mi->second ;
      ruleSet::const_iterator rsi ;
      entitySet temp_requests ;
      requests = facts.get_variable_requests(v) ;
      for(rsi = rs.begin(); rsi != rs.end(); ++rsi) {
	const rule_impl::info &rinfo = rsi->get_info().desc ;
	std::set<vmap_info>::const_iterator si ;
	entitySet t_requests = requests ;
	for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	  vector<variableSet>::const_reverse_iterator rmi ;
	  variableSet::const_iterator rvi ;
	  if(si->mapping.size() != 0) {
	    for(rmi=si->mapping.rbegin();rmi!=si->mapping.rend();++rmi) {
	      for(rvi=rmi->begin();rvi!=rmi->end();++rvi) {
		t_requests = facts.preimage(*rvi,requests).first ;
		t_requests -= my_entities ;
	      }
	    }
	    vector<variableSet>::const_iterator mi ;
	    for(mi=si->mapping.begin();mi!=si->mapping.end();++mi) {
	      for(rvi=mi->begin();rvi!=mi->end();++rvi) 
		t_requests = facts.image(*rvi, t_requests) ;
	    }
	  }
	  else 
	    t_requests = EMPTY ;
	}
	temp_requests += t_requests ;
      }
      recv_count = 0 ;
      comm_info ci ;
      ci.v = v ;
      entitySet tempset ;
      std::vector<proc_details> send_pvec ;
      std::vector<proc_details> recv_pvec ;
      
      recv_buffer = new int*[d->recv_neighbour.size()] ;
      recv_size = new int[d->recv_neighbour.size()] ;
      
      send_buffer = new int*[d->send_neighbour.size()] ;
      send_size = new int[d->send_neighbour.size()] ;
      
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	recv_size[k] = (d->recv_entities[*ei]).size() ;
	recv_buffer[k] = new int[recv_size[k]] ;
	recv_count++ ;
	k++ ;
      }
      recv_request = (MPI_Request *) malloc(recv_count * sizeof(MPI_Request) ) ;
      status = (MPI_Status *) malloc(recv_count * sizeof(MPI_Status) ) ;
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	MPI_Irecv(&recv_buffer[k][0], recv_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD, 
		  &recv_request[k]) ;  
	recv_flag = 1 ;
	k++ ;
      }
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	temp = EMPTY ;
	for(ti = temp_requests.begin(); ti != temp_requests.end(); ++ti)
	  temp += l2g[*ti] & d->send_entities[*ei] ;
	send_size[k] = temp.size() ;
	send_buffer[k] = new int[send_size[k]] ;
	k++ ;
      }
      Map id ;
      id.allocate(interval(0,d->recv_neighbour.size()-1)) ;
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	id[k] = *ei ;
	k++ ;
      }
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	proc_details pd ;
	pd.processor = *ei ;
	temp = EMPTY ;
	for(ti = temp_requests.begin(); ti != temp_requests.end(); ++ti)
	  temp += l2g[*ti] & d->send_entities[*ei] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[k][j] = *ti ;
	  pd.recv_set += d->g2l[*ti] ;
	  ++j ;
	} 
	pd.send_set = EMPTY ;
	//if(pd.recv_set != EMPTY) {
	//cout << d->myid << "precomm   sending to " << *ei << endl ; 
	//cout << d->myid << " pre comm  variable =  " << *vi << "   receiving from   " << *ei << "   " << pd.recv_set << endl ;
	//}
	MPI_Send(&send_buffer[k][0], send_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD) ;
	k++ ;
	
	recv_pvec.push_back(pd) ;
      }
      if(recv_flag)
	MPI_Waitall(recv_count, recv_request, status) ;
      
      for(k = 0; k < recv_count; ++k)
	MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
      
      for(k = 0; k < d->recv_neighbour.size(); ++k) {      
	if((recv_flag == 1) && (recv_size[k] > 0)) {
	  proc_details pd ;
	  pd.processor = id[k] ;
	  pd.send_set = EMPTY ;
	  pd.recv_set = EMPTY ;
	  for(int i = 0; i < recv_size[k]; ++i) {
	    pd.send_set += d->g2l[recv_buffer[k][i]] ;
	  }
	  //cout << d->myid << "pre  comm  variable  =   " << *vi << "  sending to   " << id[k] << "   " << pd.send_set << endl ;
	  send_pvec.push_back(pd) ;
	      //cout << d->myid << "precomm  receiving from  " << id[k] << endl ; 
	}
      }
      ci.send_info = send_pvec ;
      ci.recv_info = recv_pvec ;
      plist.push_back(ci) ;
      delete [] send_size ;
      delete [] recv_size ;
      delete [] send_buffer ;
      delete [] recv_buffer ; 
    }
    return plist ;
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
  void barrier_existential_rule_analysis(std::map<variable, ruleSet > barrier_info, fact_db &facts) {
    std::map<variable, ruleSet>::iterator mi ;
    constraint my_entities ;
    my_entities = facts.get_variable("my_entities") ;
    for(mi = barrier_info.begin() ; mi != barrier_info.end(); ++mi) {
      variable v = mi->first ;
      ruleSet rs = mi->second ;
      ruleSet::const_iterator rsi ;
      for(rsi = rs.begin(); rsi != rs.end(); ++rsi) {
	entitySet targets ;
	targets = facts.get_existential_info(v, *rsi) ;
	targets = send_entitySet(targets, facts) ;
	targets &= my_entities ;
	targets = fill_entitySet(targets, facts) ;
	facts.set_existential_info(v,*rsi,targets) ;
      }
    }
  }
  
  void barrier_process_rule_requests(std::map<variable, ruleSet > barrier_info, fact_db &facts) {
    std::map<variable, ruleSet>::iterator mi ;
    for(mi = barrier_info.begin() ; mi != barrier_info.end(); ++mi) {
      variable v = mi->first ;
      entitySet requests = facts.get_variable_requests(v) ;
      requests = send_entitySet(requests, facts) ;
      facts.variable_request(v,requests) ;
    }
  }
  
  void barrier_compiler::set_var_existence(fact_db &facts) {
    if(facts.isDistributed())
      barrier_existential_rule_analysis(barrier_info, facts) ;
  }
  
  void barrier_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      barrier_process_rule_requests(barrier_info, facts) ;
      plist = put_precomm_info(barrier_info, facts) ;
      clist = put_postcomm_info(barrier_info, facts) ;
    }
  }
  
  
  executeP barrier_compiler::create_execution_schedule(fact_db &facts) {
    //    if(num_threads > 1)
    ostringstream oss ;
    std::map<variable,ruleSet>::const_iterator mi ;
    for(mi=barrier_info.begin();mi!=barrier_info.end();++mi) {
      oss << mi->first << " " ;
    }
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ;
      el->append_list(new execute_thread_sync) ;
      el->append_list(new execute_precomm(plist, facts) ) ; 
      el->append_list(new execute_postcomm(clist, facts)) ;
      return executeP(el) ;
    }
    return new execute_thread_sync(oss.str()) ;
  }
  
  class execute_msg : public execute_modules {
    std::string msg ;
  public:
    execute_msg(string m) : msg(m) {}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  void execute_msg::execute(fact_db &facts) {  }
  
  void execute_msg::Print(std::ostream &s) const { s << msg << endl ; }
  
  
  
  void singleton_var_compiler::set_var_existence(fact_db &facts)  {
    if(facts.isDistributed())
      barrier_existential_rule_analysis(barrier_info, facts) ;
  }
  
  void singleton_var_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      barrier_process_rule_requests(barrier_info, facts) ;
    }
  }
  
  executeP singleton_var_compiler::create_execution_schedule(fact_db &facts) {
    variableSet vars ;
    std::map<variable, ruleSet>::const_iterator ri ;
    for(ri=barrier_info.begin();ri!=barrier_info.end();++ri)
      vars += ri->first ;
    ostringstream oss ;
    oss << "singleton param " << vars ;
    return executeP(new execute_msg(oss.str())) ;
  }
  
  void reduce_param_compiler::set_var_existence(fact_db &facts)  {
  }

  void reduce_param_compiler::process_var_requests(fact_db &facts) {
  }

  executeP reduce_param_compiler::create_execution_schedule(fact_db &facts) {
    ostringstream oss ;
    oss << "reduce param " << reduce_var ;
    return executeP(new execute_msg(oss.str())) ;
  }
  
  void reduce_store_compiler::set_var_existence(fact_db &facts)  {
  }

  void reduce_store_compiler::process_var_requests(fact_db &facts) {
  }

  executeP reduce_store_compiler::create_execution_schedule(fact_db &facts) {
    ostringstream oss ;
    oss << "reduce store " << reduce_var ;
    return executeP(new execute_msg(oss.str())) ;
  }
  

}

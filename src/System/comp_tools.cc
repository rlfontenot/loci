#include "comp_tools.h"

#include <mpi.h> 

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <list>
using std::list ;
#include <map>
using std::map ;
#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif
using std::hash_map ;

using std::pair ;
using std::make_pair ;


#include <distribute.h>
//#define VERBOSE

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

  /*The existential information is required to generate an execution
    schedule . This routine returns a set of entities such that the
    rule can be applied over those entities. */
  void existential_rule_analysis(rule r, fact_db &facts) {
    
    FATAL(r.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet::const_iterator ei ;
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and 
      // constraints to be within my_entities.  
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      sources &= d->my_entities ;
      constraints &= d->my_entities ;
      my_entities = d->my_entities ;
    }
    const rule_impl::info &rinfo = r.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    /*The function vmap_source_exist takes into consideration the maps 
      in the body of the rule . By looping over each of the sources in 
      the rule and also the constraints we make sure that the
      attribute specified by the target is implied by the satisfaction 
      of the attributes in the body of the rule. */
     
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      sources &= vmap_source_exist(*si,facts) ;
    }
    
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      constraints &= vmap_source_exist(*si,facts) ;
    
    if(rinfo.constraints.begin() != rinfo.constraints.end())
      if((sources & constraints) != constraints) {
	cerr << "Warning, rule " << r <<
	  " cannot supply all entities of constraint" << endl ;
	cerr << "constraints = " << constraints << endl ;
	cerr << "sources & constraints = " << (sources & constraints) << endl ;

        for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
          entitySet sources = vmap_source_exist(*si,facts) ;
          sources &= my_entities ;
          if((sources & constraints) != constraints) {
            cerr << "sources & constraints != constraints for input"
                 << endl
                 << sources  << " -- " << *si << endl ;
            
            if(si->mapping.size() > 0) {
              entitySet working = constraints ;
              for(int i=0;i<si->mapping.size();++i) {
                entitySet images ;
                variableSet::const_iterator vi ;
                for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                  images |= facts.image(*vi,working) ;
                working = images ;
              }
              variableSet::const_iterator vi ;
              for(vi=si->var.begin();vi!=si->var.end();++vi) {
                entitySet exist = facts.variable_existence(*vi) ;
                entitySet fails = working & ~exist ;
                if(fails != EMPTY) {
                  cerr << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                }
              }
            }
          }
        }
      }
    sources &= constraints ;
    //The context over which the rule is applied is given by the intersection
    // of the existential information of the sources with that of the 
    //  constraints. 
    entitySet context = sources & constraints ;
    
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      entitySet targets = vmap_target_exist(*si,facts,context) ;
      const variableSet &tvars = si->var ;
      variableSet::const_iterator vi ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
	facts.set_existential_info(*vi,r,targets) ;
#ifdef VERBOSE
        debugout << "rule " << r << " generating variable " << *vi
                           << " for entities " << targets << endl ;
#endif
      }
    }
    /*Since the unit rules are used for the existential deduction for
      the case of reduction rules these need to be treated
      separately. The information need to be communicated at this
      stage because the unit rules initializes the entities. */
    if(facts.isDistributed()) {
      if(r.get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
        WARN(r.targets().size() != 1) ;
        variable v = *r.targets().begin() ;
        entitySet exist = facts.get_existential_info(v, r) ;
        exist += fill_entitySet(exist,facts) ;
        facts.set_existential_info(v,r,exist) ;
      }
    }

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

     entitySet filter = ~EMPTY ;
     if(facts.isDistributed()) {
       fact_db::distribute_infoP d = facts.get_distribute_info() ;
       filter = d->my_entities ;
       isect = d->my_entities ;
     }

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

     // Unit rules need to apply in the clone region as well, so
     // here we make an exception for unit rules.  (this is because
     // we will be reducing to the clone region and then communicating
     // partial results.
     if(r.get_info().rule_impl->get_rule_class() != rule_impl::UNIT) {
       context &= filter ;
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
        
       // Now we have the actual requests we are making of other rules
       // so we can tell the fact database that we are now requesting
       // these values.
#ifdef VERBOSE
       debugout << "rule " << r << " requesting variables "
                          << si->var << " for entities " << requests << endl ;
#endif
       for(vi=si->var.begin();vi!=si->var.end();++vi)
         facts.variable_request(*vi,requests) ;
	
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

  /* This routine, in addition to sending the entities that are not
     owned by a particular processor,  information is stored for
     performing this communication during the execution of the
     schedule . We know the entities that a particular processor is
     supposed to send (send_entities)  . But we need to inform its neighbours
     that they are supposed to receive those entities. */ 
  std::list<comm_info>
  put_precomm_info(vector<pair<variable,entitySet> > send_entities,
                   fact_db &facts) {
    
    
    std::list<comm_info> plist ;
    if(send_entities.size() == 0)
      return plist ;
    
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      const int sesz = send_entities.size() ;
      int **send_buffer, **recv_buffer ;
      int *recv_size ;
      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[d->xmit_total_size*sesz+sesz*d->xmit.size()] ;
        recv_size[0] = d->xmit[0].size*sesz + sesz ;

        for(int i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+recv_size[i-1] ;
          recv_size[i] = d->xmit[i].size*sesz+sesz ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size*sesz+sesz*d->copy.size()] ;
        for(int i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size*sesz+sesz ;
      }
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;
      
      for(int i=0;i<d->xmit.size();++i) {
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 2,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  
      }
      for(int i=0;i<d->copy.size();++i) {
        int j=sesz ;
        for(int k=0;k<sesz;++k) {
          entitySet temp = send_entities[k].second & d->copy[i].entities ;
          send_buffer[i][k] = temp.size() ;
	  
          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
            send_buffer[i][j++] = l2g[*ei] ;
	  
          if(temp != EMPTY) {
            comm_info ci ;
            ci.v = send_entities[k].first ;
            ci.processor = d->copy[i].proc ;
            ci.send_set = temp ;
            plist.push_back(ci) ;
          }
        }
        int send_size = j ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 2,MPI_COMM_WORLD) ;
      }

      if(d->xmit.size() > 0) {
	int err = MPI_Waitall(d->xmit.size(), recv_request, status) ;
        FATAL(err != MPI_SUCCESS) ;
      }
      

      for(int i=0;i<d->xmit.size();++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=sesz ;
        for(int k=0;k<sesz;++k) {
          sequence seq ;
          for(int l=0;l<recv_buffer[i][k];++l)
            seq += d->g2l[recv_buffer[i][j++]] ;
          if(seq != EMPTY) {
            comm_info ci ;
            ci.v = send_entities[k].first ;
            ci.processor = d->xmit[i].proc ;
            ci.recv_set = seq ;
            plist.push_back(ci) ;
          }
        }
        WARN(j!=recieved) ;
      }


      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;
      
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

  bool rule_has_mapping_in_output(rule r) {
    std::set<vmap_info>::const_iterator si;
    const rule_impl::info &rinfo = r.get_info().desc ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      if(si->mapping.size() != 0)
        return true ;
    }
    return false ;
  }

  ruleSet extract_rules_with_mapping_in_output(ruleSet rs) {
    ruleSet ret ;
    for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri)
      if(rule_has_mapping_in_output(*ri))
        ret += *ri ;
    return ret ;
  }
  
   /* In the case with mapping in the output we might end up computing
     values for some of the entities in the clone region. In that case 
     we need to send these values to the processor that actually owns
     them. The information as to what entities are to be send for a
     particular variable is returned by the barrier_existential_rule_analysis routine. */

  vector<pair<variable,entitySet> >
  barrier_existential_rule_analysis(variableSet vlst,
                                    fact_db &facts) {
    vector<pair<variable,entitySet> > send_entities ;
    std::map<variable, ruleSet>::iterator mi ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    vector<entitySet> exinfo ;
    vector<ruleSet> rules ;
    vector<variable> vars ;

    vector<int> exent ;
    vector<variable> send_vars ;
    vector<rule> send_rule ;

    int ent = 0 ;
    for(variableSet::const_iterator vi=vlst.begin();vi!=vlst.end();++vi) {
      variable v = *vi ;
      ruleSet r = facts.get_existential_rules(v) ;
      
      vars.push_back(v) ;
      rules.push_back(r) ;
    }
    
    for(int i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        if(rule_has_mapping_in_output(*rsi)) {
          exent.push_back(ent) ;
          send_vars.push_back(v) ;
          send_rule.push_back(*rsi) ;
        }
        exinfo.push_back(facts.get_existential_info(v, *rsi)) ;
        ent++ ;
      }
    }
    
    
    vector<entitySet> seinfo ;
    
    map<variable,entitySet> vmap ;
    for(int i=0;i<send_vars.size();++i) {
      variable v = send_vars[i] ;
      rule rsend = send_rule[i] ;
      entitySet send_ents = exinfo[exent[i]] - d->my_entities ;
      seinfo.push_back(send_ents) ;
      vmap[v] += send_ents ;
    }
    
    for(int i=0;i<send_vars.size();++i) {
      variable v = send_vars[i] ;
      send_entities.push_back(make_pair(v,vmap[v])) ;
    }



    if(seinfo.size() != 0) {
      vector<entitySet> send_sets = send_entitySet(seinfo,facts) ;
      for(int i=0;i<seinfo.size();++i) {
        exinfo[exent[i]] += send_sets[i] ;
        exinfo[exent[i]] &= d->my_entities ;
      }
    }
    int j = 0 ;
#ifdef VERBOSE
    for(int i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        debugout << "v=" << v << ",rule ="<<*rsi
                           <<"exinfo="<<exinfo[j++] << endl ;
      }
    }
#endif
    vector<entitySet> fill_sets = fill_entitySet(exinfo,facts) ;
    
    j=0;
    for(int i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
	exinfo[j] += fill_sets[j] ;
        facts.set_existential_info(v,*rsi,exinfo[j]) ;
        ++j ;
      }
    }
    return send_entities ;
  }


  /*In this routine we fill in the communication data structure needed 
    for filling in the clone region . From the "copy" data structure we know
    what all we have to receive from the neighbouring processors(clone
    region information) . But we need to inform the neighbouring 
    processors to send us those information also. The
    clist data structure is set up so that it stores what information
    need to be send or received from a particular processor so that
    the clone region is filled up . */
  
entitySet send_requests(const entitySet& e, variable v, fact_db &facts,
                          list<comm_info> &clist) {
    entitySet re ;
    if(facts.isDistributed()) {  
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      
      int **send_buffer, **recv_buffer ;
      int *recv_size ;
      
      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;
	
        recv_buffer[0] = new int[d->xmit_total_size] ;
        recv_size[0] = d->xmit[0].size ;
	
        for(int i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->xmit[i-1].size ;
          recv_size[i] = d->xmit[i].size ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size] ;
        for(int i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(int i=0;i<d->xmit.size();++i) {
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 3,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  
      }
      
      for(int i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;
	
        comm_info ci ;
        ci.v = v ;
        ci.processor = d->copy[i].proc ;
        ci.recv_set = temp ;
        clist.push_back(ci) ;
        
        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;
	
        int send_size = temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 3,MPI_COMM_WORLD) ;
      }
      
      if(d->xmit.size() > 0) {
	int err = MPI_Waitall(d->xmit.size(), recv_request, status) ;
        FATAL(err != MPI_SUCCESS) ;
      }
      
      for(int i=0;i<d->xmit.size();++i) {
        int recieved ;
        MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        entitySet temp ;
        for(int j=0;j<recieved;++j)
          temp += d->g2l[recv_buffer[i][j]] ;
        re += temp ;
        comm_info ci ;
        ci.v = v ;
        ci.processor = d->xmit[i].proc ;
        ci.send_set = temp ;
        clist.push_back(ci) ;
      }

      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;
      
    }
    return re ;
  }
  
 
  list<comm_info>
  barrier_process_rule_requests(variableSet vars, fact_db &facts) {
    list<comm_info> clist ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      variable v = *vi ;
      entitySet requests = facts.get_variable_requests(v) ;
      requests += send_requests(requests, v, facts, clist ) ;
      requests += fill_entitySet(requests, facts) ;
      facts.variable_request(v,requests) ;
    }
    return clist ;
  }
  
  
  execute_comm::execute_comm(list<comm_info> &plist, fact_db &facts) {
    hash_map<int,vector<send_var_info> > send_data ;
    hash_map<int,vector<recv_var_info> > recv_data ;
    list<comm_info>::const_iterator cli ;
    intervalSet send_procs, recv_procs ;
    for(cli=plist.begin();cli!=plist.end();++cli) {
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        int send_proc = cli->processor ;
        send_procs += send_proc ;
        entitySet send_set = cli->send_set ;
        send_data[send_proc].push_back(send_var_info(v,send_set)) ;
      }
      if(cli->recv_set.size() > 0) {
        int recv_proc = cli->processor ;
        sequence recv_seq = cli->recv_set ;
        recv_procs += recv_proc ;
        recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
      }
    }
    
    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
    }
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
    }

    /* This part sets up the memory allocation needed for the execute routine. 
       Instead of allocating and deallocating memory each time for
       receiving a message from a processor we allocate a fixed
       message size for the receive buffer. Initially the maximum
       receive size and the maximum send size is set to be the size of
       an integer. This approach also reduces the cost incurred in
       sending the sizes in advance before the actual message is
       sent. 
    */ 

    int nsend = send_info.size() ;
    int nrecv = recv_info.size() ;
    r_size = new int[nrecv] ;
    maxr_size = new int[nrecv] ; 
    maxs_size = new int[nsend] ;
    s_size = new int[nsend] ;
    for(int i = 0; i < nrecv; ++i) {
      r_size[i] = 0 ;
      maxr_size[i] = sizeof(int) ;
    }
    for(int i = 0; i < nsend; ++i) {
      maxs_size[i] = sizeof(int) ;
      s_size[i] = 0 ;
    }
    
    recv_ptr = new unsigned char*[nrecv] ;
    send_ptr = new unsigned char*[nsend] ;
    request =  new MPI_Request[nrecv] ;
    status =  new MPI_Status[nrecv] ;
  }
  
  execute_comm::~execute_comm() {
    delete [] maxr_size ;
    delete [] maxs_size ;
    delete [] r_size ;
    delete [] s_size ;
    delete [] recv_ptr ;
    delete [] send_ptr ;
    delete [] request ;
    delete [] status ;
    
  }
  void execute_comm::execute(fact_db  &facts) {
    const int nrecv = recv_info.size() ;
    int resend_size = 0, rerecv_size = 0 ;
    std::vector<int> send_index ;
    std::vector<int> recv_index ;
    int total_size = 0 ;
    unsigned char *send_alloc, *recv_alloc ; 
    MPI_Request *re_request ;
    MPI_Status *re_status ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = 0 ;
      for(int j=0;j<recv_info[i].second.size();++j) {
	storeRepP sp = facts.get_variable(recv_info[i].second[j].v) ;
	r_size[i] += sp->pack_size(entitySet((recv_info[i].second[j].seq))) ;

#ifdef DEBUG
        entitySet rem = entitySet((recv_info[i].second[j].seq)) - sp->domain() ;
        if(rem != EMPTY)
          debugout << "variable " << recv_info[i].second[j].v << " not allocated, but recving entities " << rem << endl ;
#endif
      }
      if(r_size[i] > maxr_size[i])
	maxr_size[i] = r_size[i] ;
      else
	r_size[i] = maxr_size[i] ;
      total_size += r_size[i] ;
    }
    recv_ptr[0] = new unsigned char[total_size] ;
    recv_alloc = recv_ptr[0] ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
    
    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
    }

    /*First we find out the size of the message we are trying to
      receive using the pack_size method associated with that
      container. For static containers pack_size returns the correct
      size for the messages to be received. But for containers like
      multiStore and storeVec whose sizes gets set at run time
      pack_size returns a erroneous value. To avoid allocating the
      wrong buffer size we send the sizes first if the size returned by
      pack_size is the size of an integer or if it is greater than the
      maximum send_size(to that particular processor) . In this
      approach - in the worst case we might end up sending two messages
      always . 
      
    */
    total_size = 0 ;
    const int nsend = send_info.size() ;
    entitySet resend_procs, rerecv_procs ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(int j=0;j<send_info[i].second.size();++j) {
        storeRepP sp = facts.get_variable(send_info[i].second[j].v) ;
        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
#ifdef DEBUG
        entitySet rem = send_info[i].second[j].set - sp->domain() ;
        if(rem != EMPTY)
          debugout << "variable " << send_info[i].second[j].v << " not allocated, but sending for entities " << rem << endl ;
#endif
      }
      if((s_size[i] > maxs_size[i]) || (s_size[i] == sizeof(int))) {
	maxs_size[i] = s_size[i] ;
	int proc = send_info[i].first ;
	s_size[i] = sizeof(int) ;
	resend_procs += proc ;
	send_index.push_back(i) ;
      }
      total_size += s_size[i] ;
    }
    
    send_ptr[0] = new unsigned char[total_size] ;
    send_alloc = send_ptr[0] ;
    for(int i = 1; i < nsend; i++)
      send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
    
    // Pack the buffer for sending 
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      if(!resend_procs.inSet(send_info[i].first)) {
	for(int j=0;j<send_info[i].second.size();++j) {
	  storeRepP sp = facts.get_variable(send_info[i].second[j].v) ;
	  sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
	}
      }
      else
	MPI_Pack(&maxs_size[i], sizeof(int), MPI_BYTE, send_ptr[i], s_size[i], &loc_pack, MPI_COMM_WORLD) ; 
    }
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_PACKED,proc,1,MPI_COMM_WORLD) ;
    }
    /* We receive a message from all the processors in the
    neighbourhood. Whether the message needs to be received a second
    time is determined from the size of the message received. If the
    size of the message is equal to the size of an integer it is added
    to the list to be received a second time(even if the sent value is
    a store value or the size of the message to be received. */ 
    if(nrecv > 0) { 
      int err = MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      int *recv_sizes = new int[nrecv] ;
      for(int i = 0 ; i < nrecv; i++) {
	MPI_Get_count(&status[i], MPI_BYTE, &recv_sizes[i]) ;  
	if(recv_sizes[i] == sizeof(int)) {
	  rerecv_procs += recv_info[i].first ;
	  recv_index.push_back(i) ;
	}
      }
      delete [] recv_sizes ;
    }
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      if(rerecv_procs.inSet(recv_info[i].first)) {
	int temp ;
	/*If the size of the message received is that of an integer
	then we need to check whether it is greater than the maximum
	size received so far from that processor. If it is not then
	the maximum size is set to that value. */
	MPI_Unpack(recv_ptr[i], r_size[i], &loc_unpack, &temp, sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
	if(temp > maxr_size[i])
	  maxr_size[i] = temp ;
      }
      else
	for(int j=0;j<recv_info[i].second.size();++j) {
	  storeRepP sp = facts.get_variable(recv_info[i].second[j].v) ;
	  sp->unpack(recv_ptr[i], loc_unpack, r_size[i],
		     recv_info[i].second[j].seq) ;
	}
    }
    rerecv_size = rerecv_procs.size() ;
    resend_size = resend_procs.size() ;
    if(rerecv_size > 0) {
      re_request =  new MPI_Request[rerecv_size] ;
      re_status =  new MPI_Status[rerecv_size] ;
    }
    for(int i = 0; i < rerecv_size; i++) {
      int proc = recv_info[recv_index[i]].first ;
      recv_ptr[recv_index[i]] = new unsigned char[maxr_size[recv_index[i]]] ;
      MPI_Irecv(recv_ptr[recv_index[i]], maxr_size[recv_index[i]], MPI_PACKED, proc, 2, MPI_COMM_WORLD, &re_request[i]) ;
    }
    
    for(int i=0;i<resend_size;++i) {
      int loc_pack = 0 ;
      send_ptr[send_index[i]] = new unsigned char[maxs_size[send_index[i]]] ;
      for(int j=0;j<send_info[send_index[i]].second.size();++j) {
        storeRepP sp = facts.get_variable(send_info[send_index[i]].second[j].v) ;
	sp->pack(send_ptr[send_index[i]], loc_pack,maxs_size[send_index[i]],send_info[send_index[i]].second[j].set);
      }
    }
    
    // Send Buffer
    for(int i=0;i<resend_size;++i) {
      int proc = send_info[send_index[i]].first ;
      MPI_Send(send_ptr[send_index[i]],maxs_size[send_index[i]],MPI_PACKED,proc,2,MPI_COMM_WORLD) ;
      delete [] send_ptr[send_index[i]] ;
    }
    if(rerecv_size > 0) { 
      int err = MPI_Waitall(rerecv_size, re_request, re_status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    for(int i=0;i<rerecv_size;++i) {
      int loc_unpack = 0;
      for(int j=0;j<recv_info[recv_index[i]].second.size();++j) {
	storeRepP sp = facts.get_variable(recv_info[recv_index[i]].second[j].v) ;
	sp->unpack(recv_ptr[recv_index[i]], loc_unpack, maxr_size[recv_index[i]],
		   recv_info[recv_index[i]].second[j].seq) ;
      }
      delete [] recv_ptr[recv_index[i]] ;
    }
    
    delete recv_alloc ;
    delete send_alloc ;
    if(rerecv_size > 0) {
      delete [] re_status ;
      delete [] re_request ;
    }
  }
  
  void execute_comm::Print(ostream &s) const {
    if(send_info.size()+recv_info.size() > 0) {
      s << "communication block {" << endl ;
      if(send_info.size() > 0) {
        s << "Send:" << endl ;
        for(int i=0;i<send_info.size();++i) {
          for(int j=0;j<send_info[i].second.size();++j)
            s << "(" << send_info[i].second[j].v << "," << send_info[i].second[j].set << ") " ;
          s << " to " << send_info[i].first << endl ;
        }
      }
      if(recv_info.size() > 0) {
        s << "Recv:" << endl ;
        for(int i=0;i<recv_info.size();++i) {
          for(int j=0;j<recv_info[i].second.size();++j)
            s << "(" << recv_info[i].second[j].v << "," << recv_info[i].second[j].seq << ") " ;
          s << " from " << recv_info[i].first << endl ;
        }
      }
      s << "}" << endl ;
    }
  }
  
  
  // Sort the communication list so that the receive sequence is in the
  // order corresponding to the sending entitySet
  list<comm_info> sort_comm(list<comm_info> slist, fact_db &facts) {
    vector<pair<int,vector<send_var_info> > > send_info ;
    vector<pair<int,vector<recv_var_info> > > recv_info ;
    
    
    // First collect information from slist

    hash_map<int,vector<send_var_info> > send_data ;
    hash_map<int,vector<recv_var_info> > recv_data ;
    list<comm_info>::const_iterator cli ;
    intervalSet send_procs, recv_procs ;
    for(cli=slist.begin();cli!=slist.end();++cli) {
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        int send_proc = cli->processor ;
        send_procs += send_proc ;
        entitySet send_set = cli->send_set ;
	send_data[send_proc].push_back(send_var_info(v,send_set)) ;
      }
      if(cli->recv_set.size() > 0) {
        int recv_proc = cli->processor ;
        sequence recv_seq = cli->recv_set ;
        recv_procs += recv_proc ;
	recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
      }
    }
    
    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
    }
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
    }
    
    
    // Now build sorted comm list
    
    list<comm_info> clist ;
    
    Map l2g ;
    l2g = facts.get_variable("l2g") ;
    
    const int nrecv = recv_info.size() ;
    int *r_size = new int[nrecv] ;
    int total_size = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = 0 ;
      for(int j=0;j<recv_info[i].second.size();++j) {
        r_size[i] += recv_info[i].second[j].seq.size() ;
      }
      total_size += r_size[i] ;
    }
    int  **recv_ptr = new int*[nrecv] ;
    recv_ptr[0] = new int[total_size] ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1]+r_size[i-1] ;
    
    const int nsend = send_info.size() ;
    int *s_size = new int[nsend] ;
    total_size = 0 ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(int j=0;j<send_info[i].second.size();++j) {
        s_size[i] += send_info[i].second[j].set.size() ;
      }
      total_size += s_size[i] ;
    }
    int **send_ptr = new int*[nsend] ;
    send_ptr[0] = new int[total_size] ;
    for(int i=1;i<nsend;++i)
      send_ptr[i] = send_ptr[i-1]+s_size[i-1] ;
    
    MPI_Request *request =  new MPI_Request[nrecv] ;
    MPI_Status *status =  new MPI_Status[nrecv] ;
    
    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_INT, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
     
    }
    
    // Pack the buffer for sending 
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      for(int j=0;j<send_info[i].second.size();++j) {
        comm_info ci ;
        ci.v = send_info[i].second[j].v ;
        ci.processor = send_info[i].first ;
        ci.send_set = send_info[i].second[j].set ;
        clist.push_back(ci) ;
        for(entitySet::const_iterator ei=ci.send_set.begin();
            ei!=ci.send_set.end();
            ++ei) {
          send_ptr[i][loc_pack++] = l2g[*ei] ;
	}
      }
      WARN(loc_pack != s_size[i]) ;
    } 

    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_INT,proc,1,MPI_COMM_WORLD) ;
    }
    
    if(nrecv > 0) {
      int err = MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      for(int j=0;j<recv_info[i].second.size();++j) {
        sequence seq ;
        for(int k=0;k<recv_info[i].second[j].seq.size();++k) {
	  seq += d->g2l[recv_ptr[i][loc_unpack++]] ;
	}
	comm_info ci ;
        ci.v = recv_info[i].second[j].v ;
        ci.processor = recv_info[i].first ;
        ci.recv_set = seq ;
	clist.push_back(ci) ;
      }
      WARN(loc_unpack != r_size[i]) ;
    }
    
    delete [] status ;
    delete [] request ;
    delete [] send_ptr[0] ;
    delete [] send_ptr ;
    delete [] s_size ;
    delete [] recv_ptr[0] ;
    delete [] recv_ptr ;
    delete [] r_size ;
    
    return clist ;
  }
  
  void barrier_compiler::set_var_existence(fact_db &facts) {
    if(facts.isDistributed())
      send_entities = barrier_existential_rule_analysis(barrier_vars, facts) ;
  }
  
  void barrier_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      list<comm_info> request_comm ;
      /* The list<comm_info> returned by the
	 barrier_process_rule_requests contains the communication
	 information to send and receive the entities in the clone region*/
      request_comm = barrier_process_rule_requests(barrier_vars, facts) ;
      
      vector<pair<variable,entitySet> >::const_iterator vi ;
      vector<pair<variable,entitySet> > send_requested ;
      
      for(vi=send_entities.begin();vi!=send_entities.end();++vi) {
        variable v = vi->first ;
        entitySet send_set = vi->second ;
        send_requested.push_back(make_pair(v,send_set &
                                           facts.get_variable_requests(v))) ;
      }
      /*The put_precomm_info is used in case there is a mapping in the 
	output for any of the rules. */
      plist = put_precomm_info(send_requested, facts) ;
      clist = request_comm ;
      clist = sort_comm(request_comm,facts) ;
    }
  }

  executeP barrier_compiler::create_execution_schedule(fact_db &facts) {
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ;
      el->append_list(new execute_thread_sync) ;
      el->append_list(new execute_comm(plist, facts) ) ; 
      el->append_list(new execute_comm(clist, facts)) ;
      return executeP(el) ;
    }
    ostringstream oss ;
    oss << barrier_vars << endl ;
    return new execute_thread_sync(oss.str()) ;
  }
  
  void execute_msg::execute(fact_db &facts) {  }
  
  void execute_msg::Print(std::ostream &s) const { s << msg << endl ; }
  
  
  
  void singleton_var_compiler::set_var_existence(fact_db &facts)  {
    if(facts.isDistributed())
      barrier_existential_rule_analysis(barrier_vars, facts) ;
  }
  
  void singleton_var_compiler::process_var_requests(fact_db &facts) {
    if(facts.isDistributed()) {
      barrier_process_rule_requests(barrier_vars, facts) ;
    }
  }
  
  executeP singleton_var_compiler::create_execution_schedule(fact_db &facts) {
    variableSet vars ;
    vars = barrier_vars ;
    ostringstream oss ;
    oss << "singleton param " << vars ;
    return executeP(new execute_msg(oss.str())) ;
  }
  
}

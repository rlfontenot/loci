#include <distribute.h>
#include <Tools/debug.h>

#include "metis.h"
#include "mpi.h"
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::swap ;

namespace Loci {
  int MPI_processes ;
  int MPI_rank ;
  int num_threads = 1 ;
  
  void Init(int* argc, char*** argv)  {
    cerr << "before MPI_Init" << endl ;
    MPI_Init(argc, argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    cerr << "after MPI_Init with " << MPI_processes << ", my_rank = " << MPI_rank << endl ;
  }
  
  void metis_facts(fact_db &facts, std::vector<entitySet> &ptn, store<int> &partition ) {
    int num_partitions = MPI_processes ;
    variableSet fact_vars ;
    fact_vars = facts.get_typed_variables() ;
    entitySet map_entities ;
    variableSet::const_iterator vi ;
    entitySet::const_iterator ei ;
    for(vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;
      if(vp->RepType() == MAP) {
        MapRepP mp = MapRepP(vp->getRep()) ;
        FATAL(mp == 0) ;
        entitySet dom = mp->domain() ;
        map_entities += dom ;
        map_entities += mp->image(dom) ;
      }
    }
   
    store<entitySet> dynamic_map ;
    dynamic_map.allocate(map_entities) ;
    for(vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;      if(vp->RepType() == MAP) {
        MapRepP mp = MapRepP(vp->getRep()) ;
        FATAL(mp == 0) ;
        multiMap m = mp->get_map() ;
        entitySet dom = mp->domain() ;
        for(ei=dom.begin();ei!=dom.end();++ei) {
          for(const int *i = m.begin(*ei);i != m.end(*ei); ++i) {
            // Two associations (*ei,*i), (*i,*ei)
            dynamic_map[*i] += *ei ;
            dynamic_map[*ei]+= *i ;
          }
	 
        }
      }
    }
    int size_map = map_entities.size() ;
    Map entities ;
    Map reverse ;
    store<int> size_adj ;
    int count  = 0 ;
    entitySet dom_map = interval(0, size_map-1) ;
    entities.allocate(map_entities) ;
    partition.allocate(map_entities) ;
    size_adj.allocate(dom_map) ;
    reverse.allocate(dom_map) ;
    count = 0 ;
    for(ei = map_entities.begin(); ei!=map_entities.end(); ++ei) {
      entities[*ei] = count ;
      ++count ;
    }
    count = 0 ;
    for(ei = map_entities.begin(); ei != map_entities.end(); ++ei) {
      size_adj[count] = dynamic_map[*ei].size() ;  
      ++count ;
    }
    count = 0; 
    for(ei = map_entities.begin(); ei!=map_entities.end(); ++ei) {
      reverse[count] = *ei ;
      ++count ;
    }
   
    int *xadj = new int[size_map+1] ;
    int options, numflag, edgecut, wgtflag ;
    int *part = new int[size_map] ;
    options = 0 ;
    numflag = 0 ;
    wgtflag = 0 ;
    edgecut = 0 ;
    xadj[0] = 0 ;
    for(int i = 0; i < size_map; ++i)
      xadj[i+1] = xadj[i] + size_adj[i] ;
    int *adjncy = new int[xadj[size_map]] ;
    count = 0 ;
    for(ei = map_entities.begin(); ei != map_entities.end(); ++ei) 
      for(entitySet::const_iterator di = dynamic_map[*ei].begin(); di!=dynamic_map[*ei].end(); ++di)        {
        adjncy[count] = entities[*di] ;
        count ++ ;
      }
   
   
    METIS_PartGraphKway(&size_map,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,&num_partitions,&options,&edgecut,part) ;
    cerr << " Edge cut   " <<  edgecut << endl ;
    entitySet num_parts = interval(0, num_partitions-1) ;
    store<int> number ;
    store<int> dummy_number ;
    number.allocate(num_parts) ;
    dummy_number.allocate(num_parts) ;
    for(ei = num_parts.begin(); ei!=num_parts.end(); ++ei)
      number[*ei] = 0 ;
    for(int i = 0; i < size_map; i++)
      number[part[i]] += 1 ;
   
    for(ei = num_parts.begin(); ei!=num_parts.end(); ++ei) {
      dummy_number[*ei] = number[*ei] ;
    }
    multiMap epp ;
    epp.allocate(number) ;
    for(int i = 0; i < size_map; i++)
      epp[part[i]][--dummy_number[part[i]]] = reverse[i] ;
    int p = 0 ;
    for(ei=num_parts.begin();ei!=num_parts.end();++ei) {
      entitySet parti ;
      for(const int *ii=epp.begin(*ei);ii!= epp.end(*ei);++ii) {
        parti += *ii ;
        partition[*ii] = p ;
      }
      p++ ;
      ptn.push_back(parti) ;
    }
  }
  
  
  void get_mappings(rule_db &rdb, std::set<std::vector<variableSet> > &maps){
    ruleSet rules = rdb.all_rules() ;
    maps.clear() ;
    for(ruleSet::const_iterator ri = rules.begin(); ri != rules.end(); ++ri) {
      std::set<vmap_info>::const_iterator vmsi ;
      for(vmsi = ri->get_info().desc.targets.begin();
          vmsi != ri->get_info().desc.targets.end();
          ++vmsi) {
        if(vmsi->mapping.size() != 0) {
	  for(int i = 0; i < vmsi->mapping.size(); ++i) {
            for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
                vi != vmsi->mapping[i].end();
                ++vi) {
              variableSet v ;
              v += *vi ;
	      std::vector<variableSet> vvs ;
	      vvs.push_back(v) ;
              maps.insert(vvs) ;
            }
	  }
	}
      }
    
      for(vmsi = ri->get_info().desc.sources.begin();
          vmsi != ri->get_info().desc.sources.end();
          ++vmsi) {
        if(vmsi->mapping.size() != 0) {
	  for(int i = 0; i < vmsi->mapping.size(); i++) {
	    for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
		vi != vmsi->mapping[i].end();
		++vi) {
              variableSet v ;
              v += *vi ;
              std::vector<variableSet> vvs ;
              vvs.push_back(v) ;
              maps.insert(vvs) ;
	    }
	  }
	}
      }
    
    
      for(vmsi = ri->get_info().desc.constraints.begin();
          vmsi != ri->get_info().desc.constraints.end();
          ++vmsi) {
        if(vmsi->mapping.size() != 0) {
          for(int i = 0; i < vmsi->mapping.size(); i++) {
            for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
                vi != vmsi->mapping[i].end();
                ++vi) {
              variableSet v ;
              v += *vi ;
              std::vector<variableSet> vvs ;
              vvs.push_back(v) ;
              maps.insert(vvs) ;
            }
	  
          }
        }
      }
    }
  }
  
  
  entitySet expand_map(entitySet domain, fact_db &facts,
                       const std::set<std::vector<variableSet> > &maps) {
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      entitySet image ;
      const std::vector<variableSet> &mv = *smi ;
      for(int i = 0; i < mv.size(); ++i) {
        variableSet v = mv[i] ;
        v &= vars ;
        for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
          storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() == MAP) {
            MapRepP mp = MapRepP(p->getRep()) ;
            image += mp->image(p->domain() & locdom) ;
          }
          dom += image ;
          locdom = image ;
          image = EMPTY ;
        }
      }
    }
    return dom ;
  }

  void categories(fact_db &facts,std::vector<interval> &pvec) {
    using std::set ;
    entitySet active_set ;
    
    set<entitySet> set_of_sets ;
    set_of_sets.insert(~EMPTY) ;
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP p = facts.get_variable(*vi) ;
      if(p->RepType() == MAP) {
        MapRepP mp = MapRepP(p->getRep()) ;
        active_set += p->domain() ;
        active_set += mp->image(p->domain()) ;
      }
      if(p->RepType() == STORE) {
        active_set += p->domain() ;
      }
      
      set_of_sets.insert(p->domain()) ;
    }
    std::vector<int> vals,vals2 ;
    set<entitySet>::iterator si ;
    for(si=set_of_sets.begin();si!=set_of_sets.end();++si) {
      entitySet s = *si & active_set ;
      for(int i = 0;i < s.num_intervals(); ++i) {
        vals.push_back(s[i].first-1) ;
        vals.push_back(s[i].first) ;
        vals.push_back(s[i].second) ;
        vals.push_back(s[i].second+1) ;
      }
    }
    
    std::sort(vals.begin(),vals.end()) ;
   
    vals2.push_back(vals[0]) ;
    for(int i=1;i<vals.size();++i) {
      if(vals[i-1] == vals[i])
        continue ;
      vals2.push_back(vals[i]) ;
    }
    
    FATAL(vals2.size() < 4) ;
    pvec.clear() ;
    for(int i=1;i<vals2.size()-1;) {
      int i1 = vals2[i] ;
      int i2 = vals2[i+1] ;
      if(i1+1 == i2)
        i2 = i1 ;
      else
        ++i ;
      ++i ;
      interval iv(i1,i2) ;
      pvec.push_back(iv) ;
    }
  }
  
  
  
  void generate_distribution(fact_db &facts, rule_db &rdb, std::vector<std::vector<entitySet> > &get_entities ) {
    int num_procs = MPI_processes ;
    std::vector<entitySet> ptn ;
    store<int> partition ;
    std::set<std::vector<Loci::variableSet> > maps ;
    std::vector<entitySet> copy(num_procs) ;
    std::vector<entitySet> image(num_procs) ;
    metis_facts(facts,ptn,partition) ;
    get_mappings(rdb,maps) ;
    std::set<std::vector<Loci::variableSet> >::const_iterator smi ;
    for(int pnum = 0; pnum < num_procs; pnum++) {
      image[pnum] = expand_map(ptn[pnum], facts, maps) ;
      copy[pnum] = image[pnum] - ptn[pnum] ;
      for(int i = 0; i < num_procs; ++i) {
	entitySet slice = copy[pnum] & ptn[i] ;
	get_entities[pnum].push_back(slice) ;
      }
      get_entities[pnum][pnum] = ptn[pnum] ;
    }
  }
  
  void  distribute_facts(std::vector<std::vector<entitySet> > &get_entities, fact_db &facts)  {
    int num_procs = MPI_processes ;
    int myid = MPI_rank ;
    int j = 0 ;
    int size = 0 ;
    fact_db::distribute_infoP df = new fact_db::distribute_info  ;
    Map l2g ;
    constraint my_entities ;
    int isDistributed ;
    std::vector<Loci::interval> iv ;
    entitySet::const_iterator ei ;
    std::vector<entitySet> proc_entities ;
    categories(facts,iv) ;
    entitySet e ;
    for(int i = 0; i < iv.size(); ++i) {
      e = get_entities[myid][myid] & iv[i] ; 
      if(e != EMPTY){
	proc_entities.push_back(e) ;
	size += e.size() ;
      }
    }
    
    for(int i = 0; i < num_procs; ++i) 
      if(myid != i) {
	e = get_entities[myid][i] ;
	if(e != EMPTY) {
	  proc_entities.push_back(e) ;
	  size += e.size() ;
	}
      }
    
    entitySet g ;
    e = interval(0, size - 1) ;
    l2g.allocate(e) ;
    for(int i = 0; i < proc_entities.size(); ++i) {
      g += proc_entities[i] ;
      for(ei = proc_entities[i].begin(); ei != proc_entities[i].end(); ++ei ) {
	l2g[j] = *ei ;
	++j ;
      }
    }
    df->g2l.allocate(g) ;
    j = 0 ;
    for(int i = 0; i < proc_entities.size(); ++i) {
      for(ei = proc_entities[i].begin(); ei != proc_entities[i].end(); ++ei ) {
	df->g2l[*ei] = j ;
	++j ;
      }
    }
    
    for(int i = 0 ; i < num_procs; ++i ) 
      if(myid != i )
	if(get_entities[myid][i] != EMPTY) 
	  df->recv_neighbour += i ; 
    
    for(int i = 0; i < num_procs; ++i)
      if(myid != i)
	if(get_entities[i][myid] != EMPTY)
	  df->send_neighbour += i ;
    
    df->send_entities.allocate(df->send_neighbour) ;
    df->recv_entities.allocate(df->recv_neighbour) ;
    
    for(ei = df->recv_neighbour.begin(); ei!= df->recv_neighbour.end(); ++ei)
      df->recv_entities[*ei] = get_entities[myid][*ei] ;
    for(ei = df->send_neighbour.begin(); ei!= df->send_neighbour.end(); ++ei)
      df->send_entities[*ei] = get_entities[*ei][myid] ;
    reorder_facts(facts, df->g2l) ;
    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    g = EMPTY ;
    for(ei = get_entities[myid][myid].begin(); ei != get_entities[myid][myid].end(); ++ei)
      g+= df->g2l[*ei] ;
    my_entities = g ;
    df->myid = myid ;
    facts.put_distribute_info(df) ;
    facts.create_fact("l2g", l2g) ;
    facts.create_fact("my_entities", my_entities) ;
  }

  
  entitySet fill_entitySet(entitySet& e, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int k = 0, recv_flag = 0, recv_count = 0 ;
    store<int> is ;
    Map l2g ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    constraint my_entities ;
    entitySet re, temp;
    entitySet::const_iterator ei, ti ;
    if(facts.isDistributed()) {  
      int **send_buffer, **recv_buffer ;
      int *recv_size, *send_size ;
      
      recv_buffer = new int*[d->recv_neighbour.size()] ;
      recv_size = new int[d->recv_neighbour.size()] ;
      
      send_buffer = new int*[d->send_neighbour.size()] ;
      send_size = new int[d->send_neighbour.size()] ;
      
      l2g = facts.get_variable("l2g") ;
      my_entities = facts.get_variable("my_entities") ;
      
      re = my_entities & e ;
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	recv_size[k] = (d->recv_entities[*ei]).size() ;
	recv_buffer[k] = new int[recv_size[k]] ;
	recv_count++ ;
	k++ ;
      }
      recv_request = (MPI_Request *) malloc(recv_count * sizeof(MPI_Request) ) ;
      status = (MPI_Status *) malloc((recv_count) * sizeof(MPI_Status) ) ;
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, *ei, 1, MPI_COMM_WORLD, &recv_request[k] ) ;  
	recv_flag = 1 ;
	k++ ;
      }
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	temp = EMPTY ;
	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] & d->send_entities[*ei] ;
	send_size[k] = temp.size()  ;
	send_buffer[k] = new int[send_size[k]] ;
	k++ ;
      }
      
      k = 0 ;
      for(ei = d->send_neighbour.begin(); ei != d->send_neighbour.end(); ++ei) {
	temp = EMPTY ;
	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] & d->send_entities[*ei] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[k][j] = *ti ;
	  ++j ;
	}
	MPI_Send(&send_buffer[k][0], send_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD) ;
	k++ ;
      }
      if(recv_flag)
	MPI_Waitall(recv_count, recv_request, status) ;
      
      for(k = 0; k < recv_count; ++k)
	MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
      
      for(k = 0; k < d->recv_neighbour.size(); ++k) {      
	if((recv_flag == 1) && (recv_size[k] > 0)) {
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    
	    re += d->g2l[recv_buffer[k][i]] ;
	  }
	}
      }
      
      delete [] send_size ;
      delete [] recv_size ;
      delete [] send_buffer ;
      delete [] recv_buffer ;
      
    } 
    return re ;
  }

  entitySet send_entitySet(entitySet& e, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int k = 0, recv_flag = 0, recv_count = 0 ;
    store<int> is ;
    Map l2g ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    constraint my_entities ;
    entitySet re, temp ;
    entitySet::const_iterator ei, ti ;
    if(facts.isDistributed()) {  
      int **send_buffer, **recv_buffer ;
      int *recv_size, *send_size ;
      
      recv_buffer = new int*[d->send_neighbour.size()] ;
      recv_size = new int[d->send_neighbour.size()] ;
      
      send_buffer = new int*[d->recv_neighbour.size()] ;
      send_size = new int[d->recv_neighbour.size()] ;
      
      l2g = facts.get_variable("l2g") ;
      my_entities = facts.get_variable("my_entities") ;
      
      re = my_entities & e ;
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
	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] & d->recv_entities[*ei] ;
	send_size[k] = temp.size()  ;
	send_buffer[k] = new int[send_size[k]] ;
	k++ ;
      }
      
      k = 0 ;
      for(ei = d->recv_neighbour.begin(); ei != d->recv_neighbour.end(); ++ei) {
	temp = EMPTY ;
	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] & d->recv_entities[*ei] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[k][j] = *ti ;
	  ++j ;
	}
	MPI_Send(&send_buffer[k][0], send_size[k], MPI_INT, *ei, 1, MPI_COMM_WORLD) ;
	k++ ;
      }
      if(recv_flag)
	MPI_Waitall(recv_count, recv_request, status) ;
      
      for(k = 0; k < recv_count; ++k)
	MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
      
      for(k = 0; k < d->send_neighbour.size(); ++k) {      
	if((recv_flag == 1) && (recv_size[k] > 0)) {
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    re += d->g2l[recv_buffer[k][i]] ;
	  }
	}
      }
      
      delete [] send_size ;
      delete [] recv_size ;
      delete [] send_buffer ;
      delete [] recv_buffer ;
      
    }
    return re ;
  }

  void print_global(entitySet e, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int MAX = 100 ;
    store<int> is ;
    Map l2g ;
    entitySet::const_iterator ti ;
    fact_db::distribute_infoP d = new fact_db::distribute_info ;
    d = facts.get_distribute_info() ;
    l2g = facts.get_variable("l2g") ;
	
    if(facts.isDistributed()) {  
      if(d->myid == 0) {
	entitySet re ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	for(ti = e.begin(); ti != e.end(); ++ti)
	  re += l2g[*ti] ;
	recv_size = new int[MPI_processes-1] ;
	recv_buffer = new int*[MPI_processes-1] ;
	for(int i = 0; i < MPI_processes-1; ++i)
	  recv_buffer[i] = new int[MAX] ;
	recv_request = (MPI_Request *) malloc((MPI_processes-1) * sizeof(MPI_Request) ) ;
	status = (MPI_Status *) malloc((MPI_processes-1) * sizeof(MPI_Status) ) ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0],MAX,MPI_INT, k+1,1, MPI_COMM_WORLD, &recv_request[k] );  
	
	MPI_Waitall(MPI_processes-1, recv_request, status) ;
	
	for(k = 0; k < MPI_processes-1; ++k)
	  MPI_Get_count(&status[k], MPI_INT, &recv_size[k]) ;
	
	for(k = 0; k < MPI_processes-1; ++k) {      
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    re += recv_buffer[k][i] ;
	  }
	}
	cout << "   " << re << endl ; 
	delete [] recv_size ;
	delete [] recv_buffer ;
	
      }
      else {
	int *send_buffer;
	int send_size ;

	entitySet temp;
	send_size = e.size() ;
	send_buffer = new int[send_size] ;

	for(ti = e.begin(); ti != e.end(); ++ti)
	  temp += l2g[*ti] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[j] = *ti ;
	  ++j ;
	}
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 1, MPI_COMM_WORLD) ;
	
	delete [] send_buffer ;
      } 
    }
  }
}

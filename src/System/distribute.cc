#include <distribute.h>
#include <Tools/debug.h>
#include <stdlib.h>

//#define DEBUGGER
#ifdef DEBUGGER
#include "debugger.h"
#endif

#include <vector>
using std::vector ;

#include "metis.h"
#include <mpi.h>

#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
using std::ifstream ;

#include <algorithm>
using std::swap ;

//#define SCATTER_DIST
//#define UNITY_MAPPING

#ifdef SCATTER_DIST
#define UNITY_MAPPING
#endif

namespace Loci {
  int MPI_processes ;
  int MPI_rank ;
  int num_threads = 1 ;

  ofstream debugout[128] ;
  
  void Init(int* argc, char*** argv)  {
#ifdef DEBUGGER
    debug_execname = (*argv)[0] ;
#endif
    
    MPI_Init(argc, argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    for(int i=0;i<MPI_processes;++i) {
      ostringstream oss ;
      oss << "debug."<<MPI_rank ;
      string filename  = oss.str() ;
      debugout[i].open(filename.c_str(),ios::out) ;
    }
#ifdef DEBUGGER    
    chopsigs_() ;
#endif
  }

  void Finalize() {
    MPI_Finalize() ;
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
      if(vp->RepType() == STORE) {
	storeRepP sp = vp->getRep() ;
	map_entities += sp->domain() ;
      }
    }
    
    store<entitySet> dynamic_map ;
    dynamic_map.allocate(map_entities) ;
    for(vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;
      if(vp->RepType() == MAP) {
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

#ifdef SCATTER_DIST
    // Test code
    unsigned short int seed[3] = {0,0,1} ;

    for(int i=0;i<size_map;++i)
      part[i] = (nrand48(seed))%num_partitions ;
    // end test code
#endif
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
          std::vector<variableSet> vvs ;
	  for(int i = 0; i < vmsi->mapping.size(); ++i) {
              variableSet v ;
              for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
                vi != vmsi->mapping[i].end();
                ++vi) {
                v += variable(*vi,time_ident()) ;
              }
              vvs.push_back(v) ;
	  }
          maps.insert(vvs) ;
	}
      }
    
      for(vmsi = ri->get_info().desc.sources.begin();
          vmsi != ri->get_info().desc.sources.end();
          ++vmsi) {
        if(vmsi->mapping.size() != 0) {
          std::vector<variableSet> vvs ;
	  for(int i = 0; i < vmsi->mapping.size(); ++i) {
              variableSet v ;
              for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
                vi != vmsi->mapping[i].end();
                ++vi) {
                v += variable(*vi,time_ident()) ;
              }
              vvs.push_back(v) ;
	  }
          maps.insert(vvs) ;
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
    debugout[MPI_rank] << "expand_map maps = " << endl ;
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      entitySet image ;
      const std::vector<variableSet> &mv = *smi ;
      for(int i = 0; i < mv.size(); ++i) {
        variableSet v = mv[i] ;
        debugout[MPI_rank] << v << "->" ;
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
      debugout[MPI_rank] << endl ;
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
    entitySet::const_iterator ei, ti ;
    std::vector<entitySet> proc_entities ;
    categories(facts,iv) ;
    entitySet e ;
#ifdef DEBUG
    debugout[MPI_rank] << "categories size = " << iv.size()
                       << " {" << endl ;
    for(int i = 0; i < iv.size(); ++i) 
      debugout[MPI_rank] << iv[i] << endl ;

    debugout[MPI_rank] << "}" << endl ;
#endif
    for(int i = 0; i < iv.size(); ++i) {
      // Within each category:
      // 1) Number local processor entities first
      e = get_entities[myid][myid] & iv[i] ; 
      if(e != EMPTY){
	proc_entities.push_back(e) ;
	size += e.size() ;
      }
      // 2) Number clone region entities next
      for(int j = 0; j < num_procs; ++j) 
        if(myid != j) {
          e = get_entities[myid][j] & iv[i];
          if(e != EMPTY) {
            proc_entities.push_back(e) ;
            size += e.size() ;
          }
        }
    }


    entitySet g ;

#ifdef UNITY_MAPPING
    for(int i=0;i<proc_entities.size();++i)
      g+= proc_entities[i] ;
    l2g.allocate(g) ;
    for(entitySet::const_iterator ei=g.begin();ei!=g.end();++ei)
      l2g[*ei] = *ei ;
#else
    e = interval(0, size - 1) ;
    l2g.allocate(e) ;
    for(int i = 0; i < proc_entities.size(); ++i) {
      g += proc_entities[i] ;
      for(ei = proc_entities[i].begin(); ei != proc_entities[i].end(); ++ei ) {
	l2g[j] = *ei ;
	++j ;
      }
    }
#endif
    
    df->g2l.allocate(g) ;
    entitySet ldom = l2g.domain() ;
    for(entitySet::const_iterator ei=ldom.begin();ei!=ldom.end();++ei) {
      df->g2l[l2g[*ei]] = *ei ;
    }
      
    entitySet send_neighbour ;
    entitySet recv_neighbour ;
    store<entitySet> send_entities ;
    store<entitySet> recv_entities ;
    
    for(int i = 0 ; i < num_procs; ++i ) 
      if(myid != i )
	if(get_entities[myid][i] != EMPTY) 
	  recv_neighbour += i ; 
    
    for(int i = 0; i < num_procs; ++i)
      if(myid != i)
	if(get_entities[i][myid] != EMPTY)
	  send_neighbour += i ;
    
    send_entities.allocate(send_neighbour) ;
    recv_entities.allocate(recv_neighbour) ;
    
    for(ei = recv_neighbour.begin(); ei!= recv_neighbour.end(); ++ei) {
      for(ti =  get_entities[myid][*ei].begin(); ti != get_entities[myid][*ei].end(); ++ti)
	recv_entities[*ei] += df->g2l[*ti] ;
    }
    for(ei = send_neighbour.begin(); ei!= send_neighbour.end(); ++ei) {
      for(ti =  get_entities[*ei][myid].begin(); ti != get_entities[*ei][myid].end(); ++ti)
	send_entities[*ei] +=  df->g2l[*ti] ;
    }
    reorder_facts(facts, df->g2l) ;
    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    g = EMPTY ;
    for(ei = get_entities[myid][myid].begin(); ei != get_entities[myid][myid].end(); ++ei)
      g+= df->g2l[*ei] ;
    my_entities = g ;
    df->myid = myid ;
    df->my_entities = g ;
    for(ei=send_neighbour.begin(); ei != send_neighbour.end();++ei)
      df->xmit.push_back
        (fact_db::distribute_info::dist_data(*ei,send_entities[*ei])) ;
    for(ei=recv_neighbour.begin(); ei != recv_neighbour.end();++ei)
      df->copy.push_back
        (fact_db::distribute_info::dist_data(*ei,recv_entities[*ei])) ;
    int total = 0 ;
    for(int i=0;i<df->xmit.size();++i)
      total += df->xmit[i].size ;
    df->xmit_total_size = total ;
    total = 0 ;
    for(int i=0;i<df->xmit.size();++i)
      total += df->copy[i].size ;
    df->copy_total_size = total ;
    
      
    facts.put_distribute_info(df) ;
    facts.create_fact("l2g", l2g) ;
    facts.create_fact("my_entities", my_entities) ;
    debugout[MPI_rank] << "my_entities = " << my_entities << endl ;
    
  }
  
  

  entitySet fill_entitySet(const entitySet& e, fact_db &facts) {

    entitySet re ;

    if(facts.isDistributed()) {  
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      int **send_buffer, **recv_buffer ;
      int *recv_size ;

      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;

        recv_buffer[0] = new int[d->copy_total_size] ;
        recv_size[0] = d->copy[0].size ;
        for(int i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->copy[i-1].size ;
          recv_size[i] = d->copy[i].size ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size] ;
        for(int i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->xmit[i-1].size ;
      }
        
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;

      for(int i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(int i=0;i<d->xmit.size();++i) {
        entitySet temp = e & d->xmit[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;


          int send_size = temp.size() ;
        MPI_Send(send_buffer[i], send_size, MPI_INT, d->xmit[i].proc,
                 1, MPI_COMM_WORLD) ;
      }
        

      if(d->copy.size() > 0)
	MPI_Waitall(d->copy.size(), recv_request, status) ;
      
      for(int i = 0; i < d->copy.size(); ++i) {
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        for(int j = 0 ; j < recieved; ++j) 
          re += d->g2l[recv_buffer[i][j]] ;
      }
      

      if(d->copy.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->xmit.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] status ;
      delete [] recv_request ;
      
    } 
    return re ;
  }

  vector<entitySet> fill_entitySet(const vector<entitySet>& ev,
                                   fact_db &facts) {

    vector<entitySet> re(ev.size()) ; ;

    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      const int evsz = ev.size() ;
      int *evtmp = new int[evsz] ;
      int **send_buffer, **recv_buffer ;
      int *recv_size ;

      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;

        recv_buffer[0] = new int[d->copy_total_size*evsz+evsz*d->copy.size()] ;
        recv_size[0] = d->copy[0].size*evsz+evsz ;
        for(int i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->copy[i-1].size*evsz+evsz) ;
          recv_size[i] = d->copy[i].size*evsz+evsz ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        for(int i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
      }
        
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;

      for(int i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(int i=0;i<d->xmit.size();++i) {
        int j=evsz ;
        for(int k=0;k<evsz;++k) {
          entitySet temp = ev[k] & d->xmit[i].entities ;
          send_buffer[i][k] = temp.size() ;

          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
            send_buffer[i][j++] = l2g[*ei] ;
        }
        int send_size = j ;
        MPI_Send(send_buffer[i], send_size, MPI_INT, d->xmit[i].proc,
                 1, MPI_COMM_WORLD) ;
      }
        

      if(d->copy.size() > 0)
	MPI_Waitall(d->copy.size(), recv_request, status) ;

      for(int i = 0; i < d->copy.size(); ++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=evsz ;
        for(int k=0;k<evsz;++k) {
          for(int l=0;l<recv_buffer[i][k];++l)
            re[k] += d->g2l[recv_buffer[i][j++]] ;
        }
        WARN(j!=recieved) ;
      }
      

      if(d->copy.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->xmit.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] evtmp ;
      delete [] status ;
      delete [] recv_request ;
      
    } 
    return re ;
  }
  
  entitySet send_entitySet(const entitySet& e, fact_db &facts) {
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

      for(int i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      for(int i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;

        int send_size = temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }
      
      if(d->xmit.size() > 0)
	MPI_Waitall(d->xmit.size(), recv_request, status) ;
      
      
      for(int i=0;i<d->xmit.size();++i) {
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        for(int j=0;j<recieved;++j)
          re += d->g2l[recv_buffer[i][j]] ;
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

  vector<entitySet> send_entitySet(const vector<entitySet>& ev,
                                   fact_db &facts) {
    vector<entitySet> re(ev.size()) ;
    if(facts.isDistributed()) {  
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      const int evsz = ev.size() ;
      int **send_buffer, **recv_buffer ;
      int *recv_size ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        recv_size[0] = d->xmit[0].size*evsz + evsz ;

        for(int i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
          recv_size[i] = d->xmit[i].size*evsz+evsz ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size*evsz+evsz*d->copy.size()] ;
        for(int i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size*evsz+evsz ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(int i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      for(int i=0;i<d->copy.size();++i) {
        int j=evsz ;
        for(int k=0;k<evsz;++k) {
          entitySet temp = ev[k] & d->copy[i].entities ;
          send_buffer[i][k] = temp.size() ;

          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
            send_buffer[i][j++] = l2g[*ei] ;
        }
        int send_size = j ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }
      
      if(d->xmit.size() > 0)
	MPI_Waitall(d->xmit.size(), recv_request, status) ;
      
      
      for(int i=0;i<d->xmit.size();++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=evsz ;
        for(int k=0;k<evsz;++k)
          for(int l=0;l<recv_buffer[i][k];++l)
            re[k] += d->g2l[recv_buffer[i][j++]] ;
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
    return re ;
  }

  void print_global(entitySet e, fact_db &facts) {
    MPI_Status *status ;
    MPI_Request *recv_request ;
    int MAX = 100 ;
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

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

 
  dist_facts& distribute_info::get_dist_facts() {
    return distributed_facts ;
  }
  
  distribute_info::distribute_info(int myid) {
    distributed_facts.isDistributed.allocate(interval(myid, myid)) ;
    distributed_facts.send_neighbour.allocate(interval(myid, myid)) ;
    distributed_facts.recv_neighbour.allocate(interval(myid, myid)) ;
    distributed_facts.isDistributed[myid] = 0 ;
    distributed_facts.my_entities = EMPTY ;
  }
  
  void distribute_info::set_dist_facts(int myid, store<int> isDistributed, constraint my_entities, Map g2l, Map l2g, store<entitySet> send_neighbour, store<entitySet> recv_neighbour) {
    dist_facts d ;
    d = get_dist_facts() ;
    d.isDistributed[myid] = isDistributed[myid] ;
    d.my_entities = my_entities ;
    d.g2l = g2l ;
    d.l2g = l2g ;
    d.send_neighbour[myid] = send_neighbour[myid] ;
    d.recv_neighbour[myid] = recv_neighbour[myid] ;
  }
 
  
  void Init(int argc, char** argv, int& num_procs, int& myid)  {
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
  }
  
 void metis_facts(fact_db &facts,int num_partitions,
		  std::vector<entitySet> &ptn, store<int> &partition ) {
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


  
  void generate_distribution(fact_db &facts, rule_db &rdb, int num_procs,
			     std::vector<std::vector<entitySet> > &get_entities ) {
    std::vector<entitySet> ptn ;
    store<int> partition ;
    std::set<std::vector<Loci::variableSet> > maps ;
    std::vector<entitySet> copy(num_procs) ;
    std::vector<entitySet> image(num_procs) ;
    metis_facts(facts,num_procs,ptn,partition) ;
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
  
  void distribute_facts(std::vector<std::vector<entitySet> > &get_entities, fact_db &facts, int myid)  {
    int num_procs = get_entities.size() ;
    int j = 0 ;
    int size = 0 ;
    distribute_info dist(myid) ;
    Map g2l, l2g ;
    store<int> isDistributed;
    store<entitySet> send_neighbour, recv_neighbour ;
    isDistributed.allocate(interval(myid, myid)) ;
    send_neighbour.allocate(interval(myid, myid)) ;
    recv_neighbour.allocate(interval(myid, myid)) ;
    std::vector<Loci::interval> iv ;
    constraint my_entities ;
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
    
    g2l.allocate(g) ;
    j = 0 ;
    for(int i = 0; i < proc_entities.size(); ++i) {
      for(ei = proc_entities[i].begin(); ei != proc_entities[i].end(); ++ei ) {
	g2l[*ei] = j ;
	++j ;
      }
    }
    
    for(int i = 0 ; i < num_procs; ++i ) 
      if(myid != i )
	if(get_entities[myid][i] != EMPTY) 
	  recv_neighbour[myid] += i ; 
    
    for(int i = 0; i < num_procs; ++i)
      if(myid != i)
	if(get_entities[i][myid] != EMPTY)
	  send_neighbour[myid] += i ;
    
    reorder_facts(facts, g2l) ;
    isDistributed[myid] = 1 ;
    my_entities = get_entities[myid][myid] ;
    dist.set_dist_facts(myid, isDistributed, my_entities, g2l, l2g, send_neighbour, recv_neighbour) ;
    facts.create_fact("isDistributed", isDistributed) ;
    facts.create_fact("g2l", g2l) ;
    facts.create_fact("l2g", l2g) ;
    facts.create_fact("my_entities", my_entities) ;
    facts.create_fact("send_neighbour", send_neighbour) ;
    facts.create_fact("recv_neighbour", recv_neighbour) ;
    facts.write(cout) ;
  }
  
}

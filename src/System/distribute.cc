#include <distribute.h>
#include <Tools/debug.h>
#include <entitySet.h>

#include <stdlib.h>
#include <string.h>

#include "Tools/debugger.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <map>
using std::map ;

#include "metis.h"
#include <mpi.h>

#include <iostream>
using std::cout ;
using std::cerr ; 
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::ostream ;
using std::ostringstream ;

#include <algorithm>
using std::swap ;
using std::sort ;
//#define SCATTER_DIST
//#define UNITY_MAPPING
#ifdef SCATTER_DIST
#define UNITY_MAPPING
#endif

#include <Tools/stream.h>

namespace Loci {
  int MPI_processes = 1;
  int MPI_rank ;
  int num_threads = 1 ;
  ofstream debugout ;
  double barrier_time = 0 ;
  double total_memory_usage = 0 ;
  void Init(int* argc, char*** argv)  {
    
    char *execname = (*argv)[0] ;
    const char *hostname = "localhost" ;
    const char *debug = "gdb" ;
    MPI_Init(argc, argv) ;
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN) ;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;

    ostringstream oss ;
    oss << "debug."<<MPI_rank ;
    string filename  = oss.str() ;
    debugout.open(filename.c_str(),ios::out) ;
    
    bool debug_setup = false ;
    int i = 1 ;
    while(i<*argc) {
      if(!strcmp((*argv)[i],"--display")) {
        debug_setup = true ;
        hostname = (*argv)[i+1] ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--debug")) {
        debug = (*argv)[i+1] ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--threads")) {
        cerr << "warning --threads not yet implemented" << endl ;
        num_threads = atoi((*argv)[i+1]) ;
        i+=2 ;
      } else
        break ;
    }
    if(i!=1) {
      *argc -= (i-1) ;
      for(int k=1;k<*argc;++k)
        (*argv)[k] = (*argv)[k+i-1] ;
    }
    
    if(debug_setup) {
      setup_debugger(execname,debug,hostname) ;
      chopsigs_() ;
    }
    
  }
  
  void Finalize() {
    MPI_Finalize() ;
  }

  void Abort() {
    debugger_() ;
  }
  
  vector<entitySet> read_partition(const char *fname,int num_partitions) {
    vector<entitySet> ptn ;
    ifstream infile ;
    ostringstream oss ;
    int part ;
    oss << fname << "." << num_partitions ;
    string filename = oss.str() ; 
    infile.open(filename.c_str(), ios::in) ;
    if(infile.fail()) {
      cerr << "File " << filename <<  "   not found \n First create the file using -exit option \n " << endl ;
       Finalize() ;
      exit(0) ;
    }
    infile >> part ;
    FATAL(part != num_partitions) ;
    int *partition = new int[num_partitions] ;
    for(int i = 0; i < num_partitions; ++i) 
      infile >> partition[i] ;
    for(int i = 0; i < num_partitions; ++i) {
      entitySet parti ;
      for(int j = 0; j < partition[i]; ++j) {
        infile >> part ;
        parti += part ;
      }
      ptn.push_back(parti);
    }
    return ptn ;
  }    

  void write_partition(const char *fname, const vector<entitySet> &ptn) {
    if(MPI_rank == 0) {
      int num_partitions = ptn.size() ;
      ostringstream oss ;
      oss << fname << "." << num_partitions ;
      string filename = oss.str() ;
      ofstream ofile ;
      entitySet::const_iterator tei ;
      ofile.open(filename.c_str(), ios::out) ;
      ofile << num_partitions << endl ;
      for(int i = 0; i < ptn.size(); ++i) 
        ofile << ptn[i].size() << endl ;
      for(int i = 0; i < num_partitions; ++i) {
        entitySet temp = ptn[i];
        for(tei = temp.begin(); tei != temp.end(); ++tei)
          ofile << *tei << endl ;
      }
    }
  }
  
  void metis_facts(fact_db &facts, vector<entitySet> &ptn, int num_partitions) {
    if(num_partitions == 0)
      num_partitions = MPI_processes ;
    
    variableSet::const_iterator vi ;
    entitySet::const_iterator ei ;
    variableSet fact_vars ;
    fact_vars = facts.get_typed_variables() ;
    entitySet map_entities ;
    
    /*Initially a serial fact_database is set up on all the
      processors. We then split the fact_database into p parts if
      there are p processes. First step to partitioning is setting up
      the graph . To create the graph we need to extract the entities
      associated with the stores and the maps and their relationship
      with each other . This is probably not the most efficient way to 
      create the graph but it lets us do the partitioning without any
      prior knowledge about the problem( for example : whether there are
      faces, cells etc ...). A problem which might occur is poor load
      balancing. We are not assured that the partitioning is load
      balanced.  */
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
	map_entities += vp->domain() ;
      }
    }
    store<entitySet> dynamic_map ;
    dynamic_map.allocate(map_entities) ;
    int map_count = 0 ;
    for(vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;
      if(vp->RepType() == MAP) {
	map_count++ ;
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
    if(!map_count) {
      for(vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
	storeRepP vp = facts.get_variable(*vi) ;
	if(vp->RepType() == STORE) {
	  entitySet dom = vp->domain() ; 
	  for(ei=dom.begin();ei!=dom.end();++ei) {
	    dynamic_map[*ei] += *ei ;
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
    size_adj.allocate(dom_map) ;
    reverse.allocate(dom_map) ;
    count = 0 ;
    /* First the entities are grouped together by renumbering them. */
    for(ei = map_entities.begin(); ei!=map_entities.end(); ++ei) {
      entities[*ei] = count ;
      ++count ;
    }
    count = 0 ;
    for(ei = map_entities.begin(); ei != map_entities.end(); ++ei) {
      size_adj[count] = dynamic_map[*ei].size() ;  
      ++count ;
    }
    // Create a reverse mapping to revert to the original numbering 
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
    double t = MPI_Wtime() ;
    METIS_PartGraphKway(&size_map,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,&num_partitions,&options,&edgecut,part) ;
    double et = MPI_Wtime() ;
    debugout << "Time taken for METIS_PartGraphKway = " << et - t << "  seconds " << endl ;
    debugout << " Edge cut   " <<  edgecut << endl ;
    
    entitySet num_parts = interval(0, num_partitions-1) ;
    store<int> number ;
    store<int> dummy_number ;
    number.allocate(num_parts) ;
    dummy_number.allocate(num_parts) ;
    for(ei = num_parts.begin(); ei!= num_parts.end(); ++ei)
      number[*ei] = 0 ;
    
    for(int i = 0; i < size_map; i++) 
      number[part[i]] += 1 ;
    
    for(ei = num_parts.begin(); ei!=num_parts.end(); ++ei) {
      dummy_number[*ei] = 0 ;
    }
    multiMap epp ;
    epp.allocate(number) ;
    for(int i = 0; i < size_map; i++)
      epp[part[i]][dummy_number[part[i]]++] = reverse[i] ;
    for(ei=num_parts.begin();ei!=num_parts.end();++ei) {
      entitySet parti ;
      for(const int *ii=epp.begin(*ei);ii!= epp.end(*ei);++ii) {
        parti += *ii ; 
      }
      ptn.push_back(parti) ;
    }
    delete [] xadj ;
    delete [] part ;
    delete [] adjncy ;
  }



  
  /*This routine loops over all the rules in the database and extracts
  all the variables associated with the mappings in the head, body and
  the constraints of the rules. */

  void get_mappings(rule_db &rdb, fact_db &facts,
                    set<vector<variableSet> > &maps_ret) {
    ruleSet rules = rdb.all_rules() ;
    set<vector<variableSet> > maps ;
    
    for(ruleSet::const_iterator ri = rules.begin(); ri != rules.end(); ++ri) {
      set<vmap_info>::const_iterator vmsi ;
      for(vmsi = ri->get_info().desc.targets.begin();
          vmsi != ri->get_info().desc.targets.end();
          ++vmsi) {
        if(vmsi->mapping.size() != 0) {
          vector<variableSet> vvs ;
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
          vector<variableSet> vvs ;
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
              v += variable(*vi,time_ident()) ;
              vector<variableSet> vvs ;
              vvs.push_back(v) ;
              maps.insert(vvs) ;
            }
	  
          }
        }
      }
    }

    set<vector<variableSet> >::const_iterator mi ;

    variableSet vars = facts.get_typed_variables() ;
    maps_ret.clear() ;
    for(mi=maps.begin();mi!=maps.end();++mi) {
      const vector<variableSet> &vss = *mi ;
      // If the maps aren't in the fact database then exclude it
      variableSet notvars ;
      for(int i=0;i<vss.size();++i)
        notvars += vss[i]-vars ;
      if(notvars != EMPTY)
        continue ;
      // Now expand mapping of length up two 2
      if(vss.size() < 2) {
        for(variableSet::const_iterator vi=vss[0].begin();vi!=vss[0].end();++vi) {
          vector<variableSet> vs(1) ;
          vs[0] += *vi ;
          maps_ret.insert(vs) ;
        }
      } else {
        for(variableSet::const_iterator vi=vss[0].begin();
            vi!=vss[0].end();++vi) 
          for(variableSet::const_iterator vvi=vss[1].begin();
              vvi!=vss[1].end();++vvi) {
            vector<variableSet> vs(vss.size()) ;
            vs[0] += *vi ;
            vs[1] += *vvi ;
            for(int i=2;i<vss.size();++i)
              vs[i] = vss[i] ;
            maps_ret.insert(vs) ;
          }
      }
            
    }

  }
  
  /*The expand_map routine helps in identifying the clone regions. The
  entities which are related are found out by taking the image of the
  maps associated with the rules in the database. The entitySet which
  is usually passed on to the routine will contain the my_entities
  associated with a particular process. */
  
  entitySet expand_map(entitySet domain, fact_db &facts,
                       const set<vector<variableSet> > &maps) {
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    set<vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      const vector<variableSet> &mv = *smi ;
      for(int i = 0; i < mv.size(); ++i) {
        variableSet v = mv[i] ;
        v &= vars ;
        entitySet image ;
        for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
          storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() == MAP) {
            MapRepP mp = MapRepP(p->getRep()) ;
            image += mp->image(p->domain() & locdom) ;
          }
        }
	// The image of the map is added to the entitySet to be
	// returned 
        dom += image ;
        locdom = image ;
      }
    }
    return dom ;
  }
  
  entitySet dist_expand_map(entitySet domain, fact_db &facts,
			    const std::set<std::vector<variableSet> > &maps) {   
    std::vector<entitySet> ptn = facts.get_init_ptn() ;
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = ptn[i] ;
      ptn[i] = tmp & interval(0,  UNIVERSE_MAX) ;
    }
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      const vector<variableSet> &mv = *smi ;
      for(int i = 0; i < mv.size(); ++i) {
	variableSet v = mv[i] ;
	v &= vars ; 
	entitySet image ;
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    entitySet tmp_dom = p->domain() ;
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    entitySet out_of_dom = locdom - tmp_dom ; 
	    entitySet::const_iterator ei ;
	    entitySet tmp_out = out_of_dom ; 
	    storeRepP sp = mp->expand(tmp_out, ptn) ;
	    if(sp->domain() != tmp_dom) {
	      //debugout << " updated map " << *vi << endl ;
	      //debugout << " old_domain.Min() = " << tmp_dom.Min() << endl ;
	      //debugout << " old_domain.Max() = " << tmp_dom.Max() << endl ;
	      //debugout << " old_domain.size = " << tmp_dom.size() << endl ;
	      facts.update_fact(variable(*vi), sp) ; 
	      //debugout << " new_domain.Min() = " << sp->domain().Min() << endl ;
	      //debugout << " new_domain.Max() = " << sp->domain().Max() << endl ;
	      //debugout << " new_domain.size = " <<  sp->domain().size() << endl ;
	    }
	    image +=  MapRepP(sp)->image((sp->domain()) & locdom) ;
	  }
	}
	dom += image ;
	locdom = image ;
      }
    }
    return dom ;
  }
  
  void categories(fact_db &facts,std::vector<entitySet> &pvec) {
    entitySet active_set ;
    set<entitySet> set_of_sets ;
    set_of_sets.insert(~EMPTY) ;
    variableSet vars = facts.get_typed_variables(
) ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP p = facts.get_variable(*vi) ;
      if(p->RepType() == MAP) {
        MapRepP mp = MapRepP(p->getRep()) ;
        active_set += p->domain() ;
        active_set += mp->image(p->domain()) ;
        set_of_sets.insert(p->domain()) ;
      } else if(p->RepType() == STORE) {
        active_set += p->domain() ;
        set_of_sets.insert(p->domain()) ;
      } else {
        if(p->domain() != ~EMPTY) {
	  entitySet tmp = p->domain() ;
	  if(facts.is_distributed_start())
	    tmp = all_collect_entitySet(tmp) ;
	  active_set += tmp ;
	  set_of_sets.insert(p->domain()) ;
	}
      }
    }
    vector<int> vals,vals2 ;
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
    std::vector<interval> tmp_pvec ;
    FATAL(vals2.size() < 4) ;
    for(int i=1;i<vals2.size()-1;) {
      int i1 = vals2[i] ;
      int i2 = vals2[i+1] ;
      if(i1+1 == i2)
        i2 = i1 ;
      else 
        ++i ;
      ++i ;
      interval iv(i1,i2) ;
      if(facts.is_distributed_start())
	tmp_pvec.push_back(iv) ;
      else
	pvec.push_back(entitySet(iv)) ;
    }
    
    if(facts.is_distributed_start()) {
      std::map<variable, entitySet> vm ;
      std::map<variable, entitySet>::const_iterator mvi ;
      entitySet total_entities ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
	entitySet active_set ;
	storeRepP p = facts.get_variable(*vi) ;
	if((p->RepType() == MAP)) {
	  entitySet tmp = p->domain() ;
	  entitySet tmp_all_collect = Loci::all_collect_entitySet(tmp) ;
	  vm[*vi] = tmp_all_collect ; 
	  total_entities += tmp_all_collect ;
	}
	else if((p->RepType() == STORE)) {
	  entitySet tmp = p->domain() ;
	  entitySet all_collect = Loci::all_collect_entitySet(tmp) ;
	  vm[*vi] = all_collect ;
	  total_entities += all_collect ;
	}
	else {
	  if(p->domain() != ~EMPTY) {
	    entitySet tmp = p->domain() ;
	    entitySet all_collect = Loci::all_collect_entitySet(tmp) ;
	    vm[*vi] = all_collect ; 
	    total_entities += all_collect ;
	  }
	}
      }
      for(variableSet::const_iterator vsi = vars.begin(); vsi != vars.end(); ++vsi) {
	storeRepP p = facts.get_variable(*vsi) ;
	if((p->RepType() == MAP)) {
	  MapRepP mp = MapRepP(p->getRep()) ;
	  entitySet image_dom = mp->image(p->domain()) ;
	  entitySet out_of_set = image_dom - total_entities ;
	  entitySet left_out_categories = all_collect_entitySet(out_of_set) ;
	  if(left_out_categories != EMPTY) {
	    std::string name = "image_" ;
	    name.append(vsi->get_info().name) ;
	    variable v = variable(name) ;
	    Loci::debugout << " new variable  " << v << " is created "  << endl ;
	    vm[v] = left_out_categories ;
	  }
	}
      } 
      pvec = modified_categories(facts, vm, tmp_pvec) ;
      Loci::debugout << " pvec.size() = " << pvec.size() << endl ;
    }
  }

  vector<entitySet> generate_distribution(fact_db &facts, rule_db &rdb, int num_partitions) {
    if(num_partitions == 0)
      num_partitions = MPI_processes ;
    vector<entitySet> ptn ;
    if(num_partitions == 1)
      return ptn ;
    
    debugout << "Synchronising before metis_facts" << endl ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    double start = MPI_Wtime() ;

    metis_facts(facts,ptn,num_partitions) ;

    double end_time  = MPI_Wtime() ;
    debugout << "Time taken for metis_facts = " << end_time -start << endl ;
    return ptn ;
  } 
  
  void  distribute_facts(vector<entitySet> &ptn, fact_db &facts, rule_db &rdb) {
    if(ptn.size() == 0)
      return ;
    
    vector<vector<entitySet> > get_entities(MPI_processes) ;
    double start = MPI_Wtime() ;
    set<vector<variableSet> > maps ;
    get_mappings(rdb,facts,maps) ;
    double end_time  = MPI_Wtime() ;
    debugout << "Time taken for get_mappings =   = " << end_time -start << endl ; 
    set<vector< variableSet> >::const_iterator smi ;
    
    start = MPI_Wtime() ;
    int num_procs = MPI_processes ;
    vector<entitySet> copy(num_procs) ;
    vector<entitySet> image(num_procs) ;
    entitySet tmp ;
    //    for(int i = 0; i < num_procs; ++i)
    //      debugout << " partition[ " << i << " ] = " << ptn[i] << endl ;
    for(int pnum = 0; pnum < num_procs; pnum++) {
      image[pnum] = expand_map(ptn[pnum], facts, maps) ;
      // The clone region is obtained here 
      copy[pnum] = image[pnum] - ptn[pnum] ;
      if(pnum != MPI_rank)
	tmp += copy[pnum] ;
      for(int i = 0; i < num_procs; ++i) {
	// The information abt the clone region is found out here.  
	entitySet slice = copy[pnum] & ptn[i] ;
	get_entities[pnum].push_back(slice) ;
      }
      //The diagonal elements of the 2D vector (get_entities) contains
      //my_entities. 
      get_entities[pnum][pnum] = ptn[pnum] ;
    }
    //     debugout << " total clone information needed =  " << tmp << endl ;
    end_time  = MPI_Wtime() ;
    debugout << "Time taken for all the calls to expand_map is =   = " << end_time-start << endl ; 
    start = MPI_Wtime() ;
    int myid = MPI_rank ;
    int size = 0 ;
    fact_db::distribute_infoP df = new fact_db::distribute_info  ;
    Map l2g ;
    constraint my_entities ;
    int isDistributed ;
    vector<entitySet> iv ;
    entitySet::const_iterator ei, ti ;
    vector<entitySet> proc_entities ;
    categories(facts,iv) ;
    debugout << " initial categories.size() = " << iv.size() << endl ;
    entitySet e ;
#ifdef DEBUG
    //debugout << "categories size = " << iv.size()
    //                 << " {" << endl ;
    //for(int i = 0; i < iv.size(); ++i) 
    //debugout << iv[i] << endl ;
    
    //debugout << "}" << endl ;
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
    int j = 0 ;
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
    df->l2g = l2g ;
    df->g2l.allocate(g) ;
    entitySet ldom = l2g.domain() ;
    for(entitySet::const_iterator ei=ldom.begin();ei!=ldom.end();++ei) {
      df->g2l[l2g[*ei]] = *ei ;
    }
    entitySet send_neighbour ;
    entitySet recv_neighbour ;
    store<entitySet> send_entities ;
    store<entitySet> recv_entities ;
    
    for(int i = 0 ; i < num_procs; ++i) 
      if(myid != i )
	if(get_entities[myid][i] != EMPTY) 
	  recv_neighbour += i ; 
    
    for(int i = 0; i < num_procs; ++i)
      if(myid != i)
	if(get_entities[i][myid] != EMPTY)
	  send_neighbour += i ;
    

    send_entities.allocate(send_neighbour) ;
    recv_entities.allocate(recv_neighbour) ;

    entitySet recv, send ;
    for(ei = recv_neighbour.begin(); ei != recv_neighbour.end(); ++ei) {
      for(ti =  get_entities[myid][*ei].begin(); ti != get_entities[myid][*ei].end(); ++ti)
	recv_entities[*ei] += df->g2l[*ti] ;
      recv += recv_entities[*ei] ;
    }
    for(ei = send_neighbour.begin(); ei!= send_neighbour.end(); ++ei) {
      for(ti =  get_entities[*ei][myid].begin(); ti != get_entities[*ei][myid].end(); ++ti)
	send_entities[*ei] +=  df->g2l[*ti] ;
      send += send_entities[*ei] ;
    }
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for creating intitial info =  " << end_time - start << endl ;
    //debugout << "g2l = " << df->g2l << endl ;
    start = MPI_Wtime() ;
    reorder_facts(facts, df->g2l) ;
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for reordering =  " << end_time - start << endl ; 
    start = MPI_Wtime() ;
    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    g = EMPTY ;
    for(ei = get_entities[myid][myid].begin(); ei != get_entities[myid][myid].end(); ++ei)
      g += df->g2l[*ei] ;
    my_entities = g ;
    df->myid = myid ;
    df->my_entities = g ;
#ifdef DEBUG
    // debugout << "my_entities = " << g << endl ;
#endif
    /*xmit data structure contains the information as to what
      entities are to be send to what processor . The copy data
      structure contains the entities that are to be received from a
      particular processor(information regarding the clone region
      entities). All the entities are stored in their local
      numbering. A local to global numbering l2g  is provided to send the
      entities in their original global numbering.*/ 
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
    for(int i=0;i<df->copy.size();++i)
      total += df->copy[i].size ;
    df->copy_total_size = total ;
    
    
    facts.put_distribute_info(df) ;
    facts.create_fact("l2g", l2g) ;
    facts.put_l2g(l2g) ;
    facts.create_fact("my_entities", my_entities) ;
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for creating final info =  " << end_time - start << endl ;
    /*
    variableSet tmp_vars = facts.get_typed_variables();
    for(variableSet::const_iterator vi = tmp_vars.begin(); vi != tmp_vars.end(); ++vi) {
       storeRepP tmp_sp = facts.get_variable(*vi) ;
      if(tmp_sp->RepType() ==  MAP) {
	 debugout << " Map "  << *vi << " =  " << endl  ;
	entitySet tmp_dom = tmp_sp->domain() ;
	 debugout << " local domain = " << tmp_dom << endl ;
	entitySet global_dom =  all_collect_entitySet(tmp_dom) ;
	 debugout << " global_map domain = " << global_dom << endl ;
	entitySet tmp_image =  MapRepP(tmp_sp->getRep())->image(tmp_dom) ;
	global_dom =  all_collect_entitySet(tmp_image) ;
	 debugout << " global map image = " << global_dom << endl ;
      }
      }
    */
  }
  
  /*The fill_entitySet routine fills in the clone region entities
    . The send_buffer and the recv_buffer are allocated only once to
    contain the maximum clone region */ 
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
      /*The recv_size is the maximum possible, so that even if we
	receive a shorter message there won't be any problem */
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
        

      if(d->copy.size() > 0) {
	int err = MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }
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

  /*This is an optimization to the fill_entitySet routine to which we
    passed only an entitySet. In this case we pass in a vector of
    entitySet so that we can group the communication of entities. This 
    avoids the additional start up cost incurred when we send the
    entities corresponding to an entitySet
    each time . ie with one startup cost ts we can send all the
    entities required to a particular processor. */
  vector<entitySet> fill_entitySet(const vector<entitySet>& ev,
                                   fact_db &facts) {

    vector<entitySet> re(ev.size()) ;

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
        

      if(d->copy.size() > 0) {
	int err = MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }

      for(int i = 0; i < d->copy.size(); ++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=evsz ;
        WARN(recieved < evsz) ;
        for(int k=0;k<evsz;++k) {
          for(int l=0;l<recv_buffer[i][k];++l) {
       	    re[k] += d->g2l[recv_buffer[i][j++]] ;
	  }
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
  /*The send_entitySet routine is used to handle cases when there are
    mapping in the output. Sometimes we might compute entities in the
    clone region. Since these entities are not owned by the processor
    it needs to be send to the processor that actually owns them. */
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

      /*By intersecting the given entitySet with the clone region
	entities we can find out which entities are to be sent */
      for(int i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;

        int send_size = temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }
      
      if(d->xmit.size() > 0) {
	int err = MPI_Waitall(d->xmit.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }
      
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
      
      if(d->xmit.size() > 0) {
	int err = MPI_Waitall(d->xmit.size(), recv_request, status) ;
      	FATAL(err != MPI_SUCCESS) ;
      }
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
    if(facts.isDistributed()) {  
      MPI_Status *status ; 
      MPI_Request *recv_request ;
      int MAX = 100 ;
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = new fact_db::distribute_info ;
      d = facts.get_distribute_info() ;
      l2g = facts.get_variable("l2g") ;
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
	
	int err = MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
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
  
  /* This is a routine used for outputting the final result. Each
     processor send its entitySet (in the global numbering ) to
     processor 0.Processor 0 collects all the entitySet(including its
     local one )  and returns an entitySet in global
     numbering. */ 
  entitySet collect_entitySet(entitySet e, fact_db &facts) {
    if(!facts.isDistributed())
      return e ;
    entitySet re ;
    if(facts.isDistributed()) {  
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      d = facts.get_distribute_info() ;
      l2g = facts.get_variable("l2g") ;
      if(MPI_processes == 1) {
	entitySet temp = e ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  re += l2g[*ti] ;
        return re ;
      }
      if(d->myid == 0) {
	MPI_Status *status, *size_status ;
	MPI_Request *recv_request, *size_request ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	entitySet temp = e & d->my_entities ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  re += l2g[*ti] ;
        
	recv_size = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(k = 0; k < MPI_processes-1; k++) {
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,1, MPI_COMM_WORLD, &size_request[k]);
        }
	int err = MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	recv_buffer = new int*[MPI_processes-1] ;
        int total_size = 0 ;
        for(int i=0;i<MPI_processes-1;++i)
          total_size += recv_size[i] ;
        recv_buffer[0] = new int[total_size] ;
        
	for(int i = 1; i < MPI_processes-1; ++i)
	  recv_buffer[i] = recv_buffer[i-1]+recv_size[i-1] ;
        
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,2, MPI_COMM_WORLD, &recv_request[k] );  
	
	err = MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(k = 0; k < MPI_processes-1; ++k)       
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    FATAL(re.inSet(recv_buffer[k][i])) ;
	    re += recv_buffer[k][i] ;
	  }
        delete [] recv_buffer[0] ;
	delete [] recv_buffer ;
	delete [] recv_size ;
        delete[] size_request ;
        delete[] size_status ;
        
        delete[] recv_request ;
        delete[] status ;
      } else {
	int *send_buffer;
	int send_size ;
	entitySet temp = e & d->my_entities ;
	send_size = temp.size() ;
	send_buffer = new int[send_size] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[j] = l2g[*ti] ; 
	  re += send_buffer[j] ;
	  ++j ;
	}
	MPI_Send(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD) ;
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD) ;
	delete [] send_buffer ;
      }
    }
    return re ;
  }      

  /* This routine collects the store variables into processor
     0. Finally we have a single store allocated over the entities in
     global numbering in processor 0 */
  storeRepP collect_store(storeRepP &sp, fact_db &facts) {
    storeRepP nsp = sp ;
    if(facts.isDistributed())  {  
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      if(d->myid == 0) {
	std::vector<sequence> vseq(MPI_processes-1) ;
	MPI_Status *status, *size_status, *store_status ;
	MPI_Request *recv_request, *size_request, *store_request ;
	int **recv_buffer ;
	int *recv_size ;
	entitySet re ;
	sequence te ;
	entitySet temp = sp->domain() & d->my_entities ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  te += l2g[*ti] ;
	re += entitySet(te) ;
	recv_size = new int[MPI_processes-1] ;
	int *recv_size_bytes = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(int k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,11, MPI_COMM_WORLD, &size_request[k]);
	
	int err = MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;

	recv_buffer = new int*[MPI_processes-1] ;
        int recv_size_total = 0 ;
        for(int k=0;k<MPI_processes-1;++k) {
	  recv_size_total += recv_size[k] ;
	}
	recv_buffer[0] = new int[recv_size_total] ;
	for(int i = 1; i < MPI_processes-1; ++i)
	  recv_buffer[i] = recv_buffer[i-1] + recv_size[i-1] ;
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;
	
	for(int k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,12, MPI_COMM_WORLD, &recv_request[k] );  
	
	err = MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(int k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_size_bytes[k],1,MPI_INT, k+1,14,
		    MPI_COMM_WORLD, &size_request[k]) ;  
	
	err = MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	
	for(int k = 0; k < MPI_processes-1; ++k) {      
	  sequence tempseq ;
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    re += recv_buffer[k][i] ;
	    tempseq += recv_buffer[k][i] ;
	  }
	  vseq[k] = tempseq ;
	}
	int my_sz= sp->pack_size(temp) ;
	int my_unpack = 0 ;
	int loc_unpack = 0 ;
	int loc_pack = 0 ;
	int *r_size = new int[MPI_processes-1] ;
	store_request = new MPI_Request[MPI_processes-1] ;
	store_status = new MPI_Status[MPI_processes-1] ;
	
	int sz = 0 ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  r_size[i] = recv_size_bytes[i] ;
	  sz += r_size[i] ;
	}
	
	unsigned char **recv_ptr = new unsigned char*[MPI_processes-1] ;
	unsigned char* my_stuff = new unsigned char[my_sz] ;
	sp->pack(my_stuff, loc_pack,my_sz,temp) ;
	nsp = sp->new_store(re) ;
	//nsp = sp->new_store(EMPTY) ;
	//nsp->allocate(re) ;
	recv_ptr[0] = new unsigned char[sz] ;
	for(int i = 1; i < MPI_processes-1; i++)
	  recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
	
	for(int i = 0; i < MPI_processes-1; i++)
	  MPI_Irecv(recv_ptr[i],r_size[i] , MPI_PACKED, i+1, 13,
		    MPI_COMM_WORLD, &store_request[i]) ;
	
	err = MPI_Waitall(MPI_processes-1, store_request, store_status) ;
	FATAL(err != MPI_SUCCESS) ;
	nsp->unpack(my_stuff, my_unpack, my_sz, te) ; 
	for(int i = 0; i < MPI_processes-1; ++i) {
	  loc_unpack = 0 ;
	  nsp->unpack(recv_ptr[i], loc_unpack, r_size[i],
		      vseq[i]) ;
	}
	
	delete [] recv_size ;
        delete [] recv_size_bytes ; 
        delete [] recv_buffer[0] ;
	delete [] recv_buffer ;
	delete [] status ;
	delete [] recv_request ;
	delete [] size_request ;
	delete [] size_status ;
	delete [] store_request ;
	delete [] store_status ;
        delete [] r_size ;
        delete [] my_stuff ;
        delete [] recv_ptr[0] ;
	delete [] recv_ptr ;
      }
      else {
        
        entitySet dom = sp->domain() ;
	entitySet temp = dom & d->my_entities ;
	
	int send_size = temp.size() ;
	int *send_buffer = new int[send_size] ;
	int sz = sp->pack_size(temp) ;
	unsigned char *send_ptr = new unsigned char[sz] ;
	
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) 
	  send_buffer[j++] = l2g[*ti] ; 
	
	int loc_pack = 0;
	sp->pack(send_ptr, loc_pack, sz, temp) ;
	
	MPI_Send(&send_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 12, MPI_COMM_WORLD) ;
        MPI_Send(&sz,1,MPI_INT,0,14,MPI_COMM_WORLD) ;
	MPI_Send(send_ptr, sz, MPI_PACKED, 0, 13, MPI_COMM_WORLD) ;
	delete [] send_buffer ;
	delete [] send_ptr ;
      } 
    }
    MPI_Barrier(MPI_COMM_WORLD) ;
    return nsp ;
    
  }

  
  storeRepP collect_reorder_store(storeRepP &sp, Map& remap, fact_db &facts) {
    entitySet dom = sp->domain() ;
    dom =  collect_entitySet(dom) ;
    dom =  all_collect_entitySet(dom) ;
    std::vector<entitySet> chop_ptn = facts.get_chop_ptn() ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    entitySet remap_dom = remap.domain() ;
    Map l2g ;
    l2g = facts.get_variable("l2g") ;
    constraint my_entities ;
    my_entities = facts.get_variable("my_entities") ;
    dMap reverse ;
    remap_dom = init_ptn[ MPI_rank] & dom ;
    dmultiMap d_remap ;
    entitySet tot_remap_dom =  all_collect_entitySet(remap_dom) ;
    if(facts.is_distributed_start())
      distributed_inverseMap(d_remap, remap, tot_remap_dom, tot_remap_dom, init_ptn) ;
    else
      inverseMap(d_remap, remap, tot_remap_dom, tot_remap_dom) ;
    FORALL(remap_dom, ri) {
      if(d_remap[ri].size() == 1)
	reverse[ri] = d_remap[ri][0] ;
    } ENDFORALL ;
    Map tmp_remap ;
    entitySet tmp_remap_dom =  MapRepP(l2g.Rep())->preimage(dom&init_ptn[ MPI_rank]).first ;
    tmp_remap.allocate(tmp_remap_dom) ;
    FORALL(tmp_remap_dom, ri) {
      tmp_remap[ri] = l2g[ri] ;
    } ENDFORALL ;
    entitySet owned_entities = *my_entities & sp->domain() ;
    MapRepP(tmp_remap.Rep())->compose(reverse, tmp_remap_dom) ;
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(chop_ptn[ MPI_rank] & dom) ;
    entitySet global_owned =  MapRepP(tmp_remap.Rep())->image(owned_entities) ;
    remap_dom = chop_ptn[MPI_rank] & dom ; 
    owned_entities = tmp_remap.domain()  ;
    entitySet out_of_dom, filled_entities, local_entities ;
    
    FORALL(owned_entities, ii) {
      if(remap_dom.inSet(tmp_remap[ii])) { 
	filled_entities += tmp_remap[ii] ;
	local_entities += ii ;
      }
    } ENDFORALL ;
    storeRepP tmp_remap_sp = MapRepP(tmp_remap.Rep())->thaw() ;
    dMap d_tmp_remap(tmp_remap_sp) ; 
    qcol_rep->scatter(d_tmp_remap, sp, local_entities) ;
    out_of_dom = global_owned - filled_entities ;
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    int size_send = 0 ;
    std::vector<entitySet> copy( MPI_processes), send_clone( MPI_processes) ;
    std::vector<store<int> > scl( MPI_processes), rcl( MPI_processes) ;
    
    for(int i = 0; i <  MPI_processes; ++i) {
      send_clone[i] = out_of_dom & chop_ptn[i] ;
      entitySet tmp = MapRepP(tmp_remap.Rep())->preimage(send_clone[i]).first ;
      send_count[i] = 0 ;
      store<int> tmp_store ; 
      if(tmp.size())
	tmp_store.allocate(interval(0, tmp.size()-1)) ;
      else
	tmp_store.allocate(tmp) ;
      send_count[i] = tmp_store.Rep()->pack_size(tmp) ;
      ei = tmp.begin() ;
      for(int j = 0; j < tmp.size(); ++j) {
	tmp_store[j] = tmp_remap[*ei] ;
	++ei ;
      }
      scl[i] = tmp_store ;
      size_send += send_count[i];
    }
    unsigned char *send_buf = new unsigned char[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      scl[i].Rep()->pack(send_buf, loc_pack, size_send, scl[i].domain()) ; 
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    unsigned char *recv_buf = new unsigned char[size_send] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_PACKED,
		  recv_buf, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    std::vector< sequence> recv_dom( MPI_processes) ;
    int loc_unpack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      int tmp_size = recv_count[i] / sizeof(int) ;
      entitySet tmp ;
      if(tmp_size > 0)
	tmp = interval(0, tmp_size-1) ;
      else
	tmp = EMPTY ;
      recv_dom[i] =  sequence(tmp) ;
      store<int> tmp_store ;
      tmp_store.allocate(tmp) ;
      tmp_store.Rep()->unpack(recv_buf, loc_unpack, size_send, recv_dom[i]) ; 
      rcl[i] = tmp_store ;
    }
    size_send = 0 ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_dom[i] = MapRepP(tmp_remap.Rep())->preimage(send_clone[i]).first ;
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    }
    
    unsigned char *send_store = new unsigned char[size_send] ;
    std::vector<storeRepP> tmp_sp( MPI_processes) ;
    int size_recv = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] = 0 ;
      tmp_sp[i] = sp->new_store(entitySet(recv_dom[i])) ;
      recv_count[i] =  tmp_sp[i]->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    }
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_unpack = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tmp_sp[i]->unpack(recv_store, loc_unpack, size_recv, recv_dom[i]) ; 
    for(int i = 0; i < MPI_processes; ++i) {	
      dMap m;
      m.allocate(rcl[i].domain()) ;
      FORALL(m.domain(), mi) {
	m[mi] = rcl[i][mi] ;
      } ENDFORALL ;
      qcol_rep->scatter(m, tmp_sp[i], tmp_sp[i]->domain()) ;  
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return qcol_rep ;
    
    /*
      std::vector<entitySet> chop_ptn = facts.get_chop_ptn() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Map l2g ;
      l2g = facts.get_l2g() ;
      entitySet loc_dom = sp->domain() ;
      entitySet AR_gdom = MapRepP(l2g.Rep())->image(loc_dom) ;
      entitySet AR_tot_gdom = all_collect_entitySet(AR_gdom) ;
      entitySet final_gdom = AR_tot_gdom & chop_ptn[ MPI_rank] ;
      entitySet::const_iterator ei ;
      storeRepP g_sp = sp->remap(l2g) ;
      entitySet clone = g_sp->domain() - init_ptn[MPI_rank] ;
      entitySet remap_dom = remap.domain() ;
      dmultiMap d_remap ;
      entitySet tot_remap_dom = all_collect_entitySet(remap_dom) ;
      distributed_inverseMap(d_remap, remap, tot_remap_dom, tot_remap_dom, init_ptn) ;
      
      Map reverse ;
      remap_dom = d_remap.domain() ;
      reverse.allocate(remap_dom) ;
      
      FORALL(remap_dom, ri) {
      if(d_remap[ri].size() > 1)
      cout << " ERROR : " << endl ;
      reverse[ri] = d_remap[ri][0] ;
      } ENDFORALL ;
      
      storeRepP remap_gsp = g_sp->remap(reverse) ;
      entitySet tmp = remap_gsp->domain() ;
      clone = (chop_ptn[ MPI_rank]) - tmp ;
      entitySet tot_dom = tmp + clone ;
      storeRepP new_sp = remap_gsp->new_store(tot_dom) ;
      new_sp->copy(remap_gsp, tmp) ;
      fill_clone(new_sp, clone, init_ptn) ;
      Map new_map ;
      new_map.allocate(final_gdom) ;
      FORALL(final_gdom, fi) {
      new_map[fi] = fi ;
      } ENDFORALL ;
      storeRepP final_sp = new_sp->remap(new_map) ;
      return final_sp ;
    */
  }
  
  storeRepP distribute_reorder_store(storeRepP &sp, Map& remap, fact_db &facts) {
    return sp ;
  }
  
  storeRepP distribute_store(storeRepP &sp, fact_db &facts) {
    if(!facts.isDistributed()) {
      return sp ;
    }
    storeRepP nsp ;
    if(facts.isDistributed()) {  
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = new fact_db::distribute_info ;
      d = facts.get_distribute_info() ;
      l2g = facts.get_variable("l2g") ;
      if(d->myid == 0) {
	std::vector<entitySet> ent(MPI_processes-1) ;
	MPI_Status *status, *size_status, *store_status;
	MPI_Request *recv_request, *size_request, *store_request ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	sequence te ;
	entitySet temp, me ;
	entitySet my ;
	
	for(ti = d->my_entities.begin(); ti != d->my_entities.end(); ++ti) {
	  temp += l2g[*ti] ;
	}
	me = temp & sp->domain() ;
	for(ti = me.begin(); ti != me.end(); ++ti) {
	  te += d->g2l[*ti] ;
	  my += d->g2l[*ti] ;
	}
	
	recv_size = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(k = 0; k < MPI_processes-1; k++)
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,1, MPI_COMM_WORLD, &size_request[k]);  
	int err = MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	recv_buffer = new int*[MPI_processes-1] ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  recv_buffer[i] = new int[recv_size[i]] ;
	}
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,2, MPI_COMM_WORLD, &recv_request[k] );  
	
	err = MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(k = 0; k < MPI_processes-1; ++k) {      
	  entitySet re ;
	  for(int i = 0 ; i < recv_size[k]; ++i) 
	    re += recv_buffer[k][i] ;
	  
	  ent[k] = re ;
	}
	nsp = sp->new_store(my) ;
	//nsp = sp->new_store(EMPTY) ;
        //nsp->allocate(my) ;
	
	int sz = 0 ;
	int my_sz = sp->pack_size(me) ;
	int my_pack = 0 ;
	int my_unpack = 0 ;
	int *s_size = new int[MPI_processes-1] ;
	
	for(int i = 0; i < MPI_processes-1; ++i) { 
	  s_size[i] = sp->pack_size(ent[i]) ;
	  sz += sp->pack_size(ent[i]) ;
	}
	unsigned char **send_ptr = new unsigned char*[MPI_processes-1] ;
	unsigned char* my_stuff = new unsigned char[my_sz] ;
	store_request = new MPI_Request[MPI_processes-1] ;
	store_status = new MPI_Status[MPI_processes-1] ;
	sp->pack(my_stuff, my_pack,my_sz,me) ;
	send_ptr[0] = new unsigned char[sz] ;
	for(int i = 1; i < MPI_processes-1; i++)
	  send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
	for(int i = 0; i < MPI_processes-1; i++) {
	  int loc_pack = 0 ;
	  sp->pack(send_ptr[i], loc_pack, s_size[i], ent[i]) ;
	  MPI_Isend(send_ptr[i], s_size[i], MPI_PACKED, i+1, 3,
		    MPI_COMM_WORLD, &store_request[i]) ;
	}
	err = MPI_Waitall(MPI_processes-1, store_request, store_status) ;
	FATAL(err != MPI_SUCCESS) ;
	nsp->unpack(my_stuff, my_unpack, my_sz, te) ; 
	delete [] recv_size ;
	delete [] recv_buffer ;
	delete [] status ;
	delete [] recv_request ;
	delete [] size_request ;
	delete [] size_status ;
	delete [] store_request ;
	delete [] store_status ;
	delete [] send_ptr ;
	
      }
      else {
	int *send_buffer ;
	int send_size ;
	int loc_unpack = 0 ;
	unsigned char *recv_ptr ;
	MPI_Status stat ;
	entitySet re ;
	sequence tempseq ;
	entitySet my = sp->domain() & d->my_entities ;
	send_size = my.size() ;
	send_buffer = new int[send_size] ;
	int sz = sp->pack_size(my) ;
	
	recv_ptr = new unsigned char[sz] ;
	nsp = sp->new_store(EMPTY) ;
        nsp->allocate(my) ;
	int j = 0 ;
	
	for(ti = my.begin(); ti != my.end(); ++ti) {
	  re += l2g[*ti] ;
	}
	for(ti = re.begin(); ti != re.end(); ++ti) {
	  send_buffer[j] = *ti  ;
	  ++j ;
	}
	for(ti = re.begin(); ti != re.end(); ++ti)
	  tempseq += d->g2l[*ti] ;
	
	MPI_Send(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD) ;
	
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD) ;
	MPI_Recv(recv_ptr, sz, MPI_PACKED, 0, 3, MPI_COMM_WORLD, &stat) ;
	nsp->unpack(recv_ptr, loc_unpack, sz,tempseq) ;
	delete [] send_buffer ;
	delete [] recv_ptr ;
	
      } 
    }   
    return nsp ;
     
  }
  
  void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    
    entitySet preloop = input_preimage & input_map.domain() ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;
    std::vector<std::vector<int> > map_elems(MPI_processes) ;
    std::vector<int> tmp_vec ;
    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
      FORALL(local_input_image,i) {
	result[i] = tmp_vec ;
      }ENDFORALL ;
      
      FORALL(preloop,i) {
	int elem = input_map[i] ;
	if(input_image.inSet(elem)) 
	  for(int j = 0; j < MPI_processes; ++j)
	    if(init_ptn[j].inSet(elem)) {
	      map_elems[j].push_back(elem) ;
	      map_elems[j].push_back(i) ;
	    }
      }ENDFORALL ;
      for(int i = 0; i < MPI_processes; ++i) 
	send_count[i] = map_elems[i].size() ;
      
      int size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	size_send += send_count[i] ; 
      int *send_buf = new int[size_send] ;
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ; 
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i) {
	size_send += recv_count[i] ;
      }
      int *recv_buf = new int[size_send] ;
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	for(std::vector<int>::const_iterator vi = map_elems[i].begin(); vi != map_elems[i].end(); ++vi) {
	  send_buf[size_send] = *vi ;
	  ++size_send ;
	}
      send_displacement[0] = 0 ;
      recv_displacement[0] = 0 ;
      for(int i = 1; i < MPI_processes; ++i) {
	send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
	recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
      }
      MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		    recv_buf, recv_count, recv_displacement, MPI_INT,
		    MPI_COMM_WORLD) ;  
      HASH_MAP(int, std::set<int> ) hm ;
      for(int i = 0; i < MPI_processes; ++i) {
	for(int j = recv_displacement[i]; j <
	      recv_displacement[i]+recv_count[i]-1; ++j) {
	  hm[recv_buf[j]].insert(recv_buf[j+1]);
	  j++ ;
	}
      }
      for(HASH_MAP(int, std::set<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
	for(std::set<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	  result[hmi->first].push_back(*si) ;
      
      delete [] recv_count ;
      delete [] send_count ;
      delete [] send_displacement ;
      delete [] recv_displacement ;
      delete [] send_buf ;
      delete [] recv_buf ;
    }
  
  void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    
    entitySet preloop = input_preimage & input_map.domain() ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes];
    int *recv_displacement = new int[MPI_processes];
    std::vector<std::vector<int> > map_elems(MPI_processes) ;
    std::vector<int> tmp_vec ;
    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    FORALL(local_input_image,i) {
      result[i] = tmp_vec ;
    }ENDFORALL ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      std::vector<int> tmp_vec ;
      if(input_image.inSet(elem)) 
	for(int j = 0; j < MPI_processes; ++j)
	  if(init_ptn[j].inSet(elem)) {
	    map_elems[j].push_back(elem) ;
	    map_elems[j].push_back(i) ;
	  }
    }ENDFORALL ;
    
    for(int i = 0; i < MPI_processes; ++i) 
      send_count[i] = map_elems[i].size() ;
    
    int size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += map_elems[i].size() ;
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      size_send += recv_count[i] ;
    }
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      for(std::vector<int>::const_iterator vi = map_elems[i].begin(); vi != map_elems[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, std::set<int> ) hm ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	hm[recv_buf[j]].insert(recv_buf[j+1]);
	j++ ;
      }
    }
    for(HASH_MAP(int, std::set<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
      for(std::set<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	result[hmi->first].push_back(*si) ;
    
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
  }
  void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    entitySet preloop = input_preimage & input_map.domain() ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes];
    int *recv_displacement = new int[MPI_processes];
    std::vector<HASH_MAP(int, std::vector<int> > ) map_elems(MPI_processes); 
    std::vector<int> tmp_vec ;
    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    FORALL(local_input_image,i) {
      result[i] = tmp_vec ;
    }ENDFORALL ;
    
    FORALL(preloop,i) {
      for(std::vector<int>::const_iterator mi = input_map[i].begin(); mi != input_map[i].end(); ++mi) {
	int elem = *mi ;
	if(input_image.inSet(elem)) 
	  for(int j = 0; j < MPI_processes; ++j)
	    if(init_ptn[j].inSet(elem)) 
	      (map_elems[j])[elem].push_back(i) ;
      }
    } ENDFORALL ;
    
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] = 2*map_elems[i].size() ;
      for(HASH_MAP(int, std::vector<int> )::iterator hi = map_elems[i].begin(); hi != map_elems[i].end(); ++hi)
	send_count[i] += hi->second.size() ; 
    }
    int size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += send_count[i] ;
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      size_send += recv_count[i] ;
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(HASH_MAP(int, std::vector<int> )::const_iterator miv = map_elems[i].begin(); miv != map_elems[i].end(); ++miv) {
	send_buf[size_send] = miv->first ;
	++size_send ;
	send_buf[size_send] = miv->second.size() ;
	++size_send ;
	for(std::vector<int>::const_iterator vi = miv->second.begin(); vi != miv->second.end(); ++vi) { 
	  send_buf[size_send] = *vi ;
	  ++size_send ;
	}
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, std::set<int> ) hm ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	int count = recv_buf[j+1] ;
	for(int k = 0; k < count; ++k)
	  hm[recv_buf[j]].insert(recv_buf[j+k+2]);
	j += count + 1 ;
      }
    }
    for(HASH_MAP(int, std::set<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
      for(std::set<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	result[hmi->first].push_back(*si) ;
    
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
  }
  
  void distributed_inverseMap(dmultiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    entitySet preloop = input_preimage & input_map.domain() ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes];
    int *recv_displacement = new int[MPI_processes];
    std::vector<HASH_MAP(int, std::vector<int> > ) map_elems(MPI_processes); 
    std::vector<int> tmp_vec ;
    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    FORALL(local_input_image,i) {
      result[i] = tmp_vec ;
    }ENDFORALL ;
    
    FORALL(preloop,i) {
      for(const int *mi = input_map.begin(i); mi != input_map.end(i); ++mi) {
	int elem = *mi ;
	if(input_image.inSet(elem)) 
	  for(int j = 0; j < MPI_processes; ++j)
	    if(init_ptn[j].inSet(elem)) 
	      (map_elems[j])[i].push_back(elem) ;
      }
    }ENDFORALL ;
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] = 2*map_elems[i].size() ;
      for(HASH_MAP(int, std::vector<int> )::iterator hi = map_elems[i].begin(); hi != map_elems[i].end(); ++hi)
	send_count[i] += hi->second.size() ; 
    }
    int size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += send_count[i] ;
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      size_send += recv_count[i] ;
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(HASH_MAP(int, std::vector<int> )::const_iterator miv = map_elems[i].begin(); miv != map_elems[i].end(); ++miv) {
	send_buf[size_send] = miv->first ;
	++size_send ;
	send_buf[size_send] = miv->second.size() ;
	++size_send ;
	for(std::vector<int>::const_iterator vi = miv->second.begin(); vi != miv->second.end(); ++vi) { 
	  send_buf[size_send] = *vi ;
	  ++size_send ;
	}
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, std::set<int> ) hm ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	int count = recv_buf[j+1] ;
	for(int k = 0; k < count; ++k)
	  hm[recv_buf[j+k+2]].insert(recv_buf[j]);
	j += count + 1 ;
      }
    }
    for(HASH_MAP(int, std::set<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
      for(std::set<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	result[hmi->first].push_back(*si) ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
  }
  
 
  void fill_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei) {
	copy[i].push_back(*ei) ;
      }
      std::sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      std::sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_dom[i] += *vi ;
      }
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] =  sp->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
  }
  
  std::vector< storeRepP> fill_global_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec( MPI_processes) ;
    std::vector< storeRepP> tmp_sp( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      recv_count[i] =  sp->pack_size(entitySet(recv_dom[i])) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp[i]->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }
  
  std::pair< storeRepP, Map > send_clone( storeRepP& sp, Map &m,  entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 0 ;
      FORALL(send_dom[i], si) {
	entitySet tmp = interval(si, si) ;
	send_count[i] +=  sp->pack_size(tmp) ;
      } ENDFORALL ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    Map m_int ;
    int tot = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      tot += entitySet(recv_dom[i]).size() ; //+= entitySet(recv_dom[i]) ;
    
    tmp_dom = interval(0, tot-1) ;
    
     storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] = 0 ; 
      FORALL(entitySet(recv_dom[i]), ri) {
	entitySet tmp = interval(ri, ri) ; 
	recv_count[i] +=  tmp_sp->pack_size(tmp) ;
      } ENDFORALL ;
      size_recv += recv_count[i] ;
    }
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    dmultiMap inverse ;
    entitySet m_dom = m.domain() ;
    entitySet range =  MapRepP(m.Rep())->image(m_dom) ;
    inverseMap(inverse, m, range, m_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      FORALL(send_dom[i], si) {
	entitySet tmp = interval(inverse[si][0], inverse[si][0]) ;
	sp->pack(send_store, loc_pack, size_send, tmp) ;
      } ENDFORALL ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    m_int.allocate(tmp_dom) ;
    tot = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(sequence::const_iterator si = recv_dom[i].begin() ; si != recv_dom[i].end(); ++si) {
	m_int[tot] = *si ;
	++tot ;
      }
    }
    FORALL(tmp_dom, ti) {
      entitySet tmp = interval(ti, ti) ;
      sequence tmp_seq = sequence(tmp) ;
      tmp_sp->unpack(recv_store, loc_pack, size_recv, tmp_seq) ; 
    } ENDFORALL ;
    
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return std::make_pair(tmp_sp, m_int)  ;
  }
  
   storeRepP send_clone_non( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
  int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    for(int i = 0; i <  MPI_processes; ++i)
      tmp_dom += entitySet(recv_dom[i]) ;
     storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] =  tmp_sp->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }
   std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) {
  int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec( MPI_processes) ;
    std::vector< storeRepP> tmp_sp( MPI_processes) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      recv_count[i] =  tmp_sp[i]->pack_size(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    
     loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      tmp_sp[i]->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
 }
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	hm[recv_map[j]] = recv_map[j+1];
	j++ ;
      }
    }
    dMap tmp_dm ;
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) {
      tmp_dm[hmi->first] = hmi->second ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_dm ;
  }
 std::vector<dMap> send_global_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    std::vector<HASH_MAP(int, int) > hm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm ;
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	tmp_hm[recv_map[j]] = recv_map[j+1];
	j++ ;
      }
      hm[i] = tmp_hm ;
    }
    std::vector<dMap> v_dm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm = hm[i] ; 
      dMap tmp_dm ;
      for(HASH_MAP(int, int)::const_iterator hmi = tmp_hm.begin(); hmi != tmp_hm.end(); ++hmi) 
	tmp_dm[hmi->first] = hmi->second ;
      v_dm[i] = tmp_dm ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return v_dm ;
  }
  
 entitySet all_collect_entitySet(const entitySet &e) {
    entitySet collect ;
    entitySet::const_iterator ei ;
    if(MPI_processes > 1) {
      int *recv_count = new int[MPI_processes] ;
      int *send_count = new int[MPI_processes] ;
      int *send_displacement = new int[MPI_processes];
      int *recv_displacement = new int[MPI_processes];
      int size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i) {
	send_count[i] = 2 * e.num_intervals() ;
	size_send += send_count[i] ; 
      }
      int *send_buf = new int[size_send] ;
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ; 
      MPI_Barrier(MPI_COMM_WORLD) ;
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	size_send += recv_count[i] ;
      int *recv_buf = new int[size_send] ;
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	for(int j = 0; j < e.num_intervals(); ++j) {
	  send_buf[size_send] = e[j].first ;
	  ++size_send ;
	  send_buf[size_send] = e[j].second ;
	  ++size_send ;
	}
      send_displacement[0] = 0 ;
      recv_displacement[0] = 0 ;
      for(int i = 1; i < MPI_processes; ++i) {
	send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
	recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
      }
      MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		    recv_buf, recv_count, recv_displacement, MPI_INT,
		    MPI_COMM_WORLD) ;
      std::vector<entitySet> add(MPI_processes) ;
      for(int i = 0; i < MPI_processes; ++i) {
	for(int j = recv_displacement[i]; j <
	      recv_displacement[i]+recv_count[i]-1; ++j) { 
	  collect +=  interval(recv_buf[j], recv_buf[j+1]) ;
	  j++ ;
	}
      }
      delete [] send_count ;
      delete [] recv_count ;
      delete [] recv_displacement ;
      delete [] send_displacement ;
      delete [] send_buf ;
      delete [] recv_buf ;
    }
    else
      return e ;
    return collect ;
  }

  
  std::vector<int> all_collect_sizes(int size) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    for(int i = 0; i <  MPI_processes; ++i) 
      send_count[i] = size ;
    
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    std::vector<int> vset( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i)
      vset[i] = recv_count[i] ;
    
    delete [] send_count ;
    delete [] recv_count ;
    return vset ;
  }

 std::vector<entitySet> all_collect_vectors(entitySet &e) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes];
    int *recv_displacement = new int[ MPI_processes];
    int size_send = 0 ;
    entitySet::const_iterator ei ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = e.size() ;
      size_send += send_count[i] ; 
    }  
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      for(ei = e.begin(); ei != e.end(); ++ei) {
	send_buf[size_send] = *ei ;
	++size_send ;
      }
    } 
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    } 
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    std::vector<entitySet> vset( MPI_processes) ;
    
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	vset[i] += recv_buf[j] ;
    } 
    delete [] send_count ;
    delete [] recv_count ;
    delete [] recv_displacement ;
    delete [] send_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
    return vset ;
  }

  
class cat_tree : public CPTR_type {
public:
  typedef CPTR<cat_tree>     treeP ;
  typedef std::list<treeP> treeList ;
  variableSet                generalization ;
  treeList                    before, after ;
  variableSet get_below_vars() ;
  int is_swapped ;
  void Print(ostream &s) ;
  void fill_list(std::list<treeP> &list_tp) ;
  std::list<variableSet> get_categories(variableSet &vs) ;
  void recursive_reorder(std::list<variableSet> &lvs) ;
  void reorder(std::list<variableSet> &lvs) ;
  void down(int &i) ;
  bool operator==(const treeP &tp) const {
    if((generalization == tp->generalization) && (before.size() == tp->before.size()) && (after.size() == tp->after.size()))
      return 1 ;
    else
      return 0 ;
  }
  void before2after() ;
} ;

typedef cat_tree::treeP     treeP ;
typedef cat_tree::treeList  treeList ;

inline bool operator<(const treeP &tp1, const treeP &tp2) {
  
  if(tp1->generalization != tp2->generalization)
    return 0 ;
  else
    return  tp1->generalization < tp2->generalization ||
      (tp1->generalization==tp2->generalization &&
       (tp1->before.size()< tp2->before.size())) ; 
}

variableSet cat_tree::get_below_vars() {
  variableSet vs = generalization ;
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    vs += (*li)->get_below_vars() ;
  for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
    vs += (*li)->get_below_vars() ;
  return vs ;
}

void cat_tree::Print(ostream &s) {
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    (*li)->Print(s) ;
  s << "  generalization = " << generalization << "  before.size = " << before.size() << "  after.size()  =  " << after.size() << endl ;
  for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
    (*li)->Print(s) ;
}

void cat_tree::down(int &i ) {
  i++ ;
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    (*li)->down(i) ;
}

std::list<variableSet> cat_tree::get_categories(variableSet &vs) {
  vs += generalization ;
  variableSet tvs  = vs ;
  std::list<variableSet> lvs ;
  
  if(before.size() > 1) {
    for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
      variableSet tmp_vs = vs;
      tmp_vs += (*li)->generalization ;
      std::list<variableSet> tmp_lvs = (*li)->get_categories(tmp_vs) ;
      for(std::list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
	lvs.push_back(*lvi) ;
    }
  }
  else {
    for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
      variableSet tmp_vs = vs;
      std::list<variableSet> tmp_lvs ;
      if(after.size()) {
	tmp_vs += (*li)->generalization ;
	tmp_lvs = (*li)->get_categories(tmp_vs) ;
	for(std::list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
	  lvs.push_back(*lvi) ;
      }
      else {
	vs += (*li)->generalization ;
	tmp_lvs = (*li)->get_categories(vs) ;
	for(std::list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
	  lvs.push_back(*lvi) ;
      }  
    }
  }
  if(tvs != EMPTY)
    lvs.push_back(tvs) ;
  if(after.size() > 1) {
    for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li) {
      variableSet tmp_vs = tvs;
      tmp_vs += (*li)->generalization ;
      std::list<variableSet> tmp_lvs = (*li)->get_categories(tmp_vs) ;
      for(std::list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
	lvs.push_back(*lvi) ;
    }
  }
  else {
    for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li) {
      variableSet tmp_vs = vs ;
      std::list<variableSet> tmp_lvs ;
      if(before.size()) {
	tmp_vs += (*li)->generalization ;
	tmp_lvs = (*li)->get_categories(tmp_vs) ;
      }
      else {
	vs += (*li)->generalization ;
	tmp_lvs = (*li)->get_categories(vs) ;
      }
      for(std::list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
	lvs.push_back(*lvi) ;
    }  
  }
  
  return lvs ;
}

void cat_tree::fill_list(std::list<treeP> &list_tp) {
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    (*li)->fill_list(list_tp) ;
  list_tp.push_front(this) ;
  for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
    (*li)->fill_list(list_tp) ;
}

void cat_tree::before2after() {
  std::list<treeP> vt ;
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
    if((*li)->before.size())
      vt.push_back(*li) ;
    else
      after.push_back(*li) ;
  }
  before = vt ;
}

void cat_tree::reorder(std::list<variableSet> &lvs) {
  std::set<variableSet> tmp_vars ;
  std::list<treeP> lone_ltp, tmp_buf, tmp_before ;
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    if((*li)->before.size() || (*li)->after.size())
      tmp_buf.push_back(*li) ;
    else
      lone_ltp.push_back(*li) ;
  for(std::list<treeP>::iterator li = tmp_buf.begin(); li != tmp_buf.end(); ++li) {
    tmp_before.push_back(*li) ;
    tmp_vars.insert((*li)->generalization) ;
    variableSet vset = (*li)->get_below_vars() ;
    for(std::list<treeP>::iterator tmp_li = lone_ltp.begin(); tmp_li != lone_ltp.end(); ++tmp_li)
      if((((*tmp_li)->generalization) & vset) != EMPTY)
	if(tmp_vars.find((*tmp_li)->generalization) == tmp_vars.end()) {
	  tmp_before.push_back(*tmp_li) ;
	  tmp_vars.insert((*tmp_li)->generalization) ;
	}
  }
  for(std::list<treeP>::iterator tmp_li = lone_ltp.begin(); tmp_li != lone_ltp.end(); ++tmp_li)
    if(tmp_vars.find((*tmp_li)->generalization) == tmp_vars.end()) {
      tmp_vars.insert((*tmp_li)->generalization) ;
      tmp_before.push_back(*tmp_li) ;
    }
  before = tmp_before ;
  variableSet int_vset ;
  for(std::list<variableSet>::iterator li = lvs.begin(); li != lvs.end(); ++li) {
    for(variableSet::const_iterator vi = (*li).begin(); vi != (*li).end(); ++vi)
      int_vset += *vi ;
  }
  for(std::list<treeP>::iterator lii = before.begin(); lii != before.end(); ++lii) {
    entitySet vs = entitySet((*lii)->get_below_vars()) ;
    for(std::list<treeP>::iterator lij = lii; lij != before.end(); ++lij) {
      if(lii != lij) {
	entitySet next_vs =  entitySet((*lij)->get_below_vars()) ;
	entitySet tmp_vs = vs & next_vs ;
	tmp_vs -= (tmp_vs & entitySet(int_vset)) ;
	if(tmp_vs != EMPTY) {
	  lvs.push_front(variableSet(tmp_vs)) ;
	  int_vset += variableSet(tmp_vs) ;
	  std::list<treeP> tlp, tmp_tlp ;
	  for(std::list<treeP>::iterator tlpi = (*lij)->before.begin(); tlpi != (*lij)->before.end(); ++tlpi)
	    if(((*tlpi)->get_below_vars() & variableSet(tmp_vs)) == EMPTY)
	      (*lij)->after.push_back(*tlpi) ;
	    else
	      tlp.push_back(*tlpi) ;
	  (*lij)->before = tlp ; 
	  (*lij)->is_swapped = 1 ;
	}
      }
    }
    for(std::list<treeP>::iterator lij = before.begin(); lij != lii; ++lij) {
      if(!(*lij)->is_swapped) {
	entitySet vs = entitySet((*lij)->get_below_vars()) ;
	for(std::list<variableSet>::iterator lvi = lvs.begin(); lvi != lvs.end(); ++lvi) {
	  entitySet next_vs =  entitySet(*lvi) ;
	  entitySet tmp_vs = vs & next_vs ;
	  if(tmp_vs != EMPTY) {
	    std::list<treeP> tlp, tmp_tlp ;
	    
	    for(std::list<treeP>::iterator tlpi = (*lij)->before.begin(); tlpi != (*lij)->before.end(); ++tlpi)
	      if(((*tlpi)->get_below_vars() & variableSet(tmp_vs)) == EMPTY)
		tmp_tlp.push_back(*tlpi) ;
	      else
		(*lij)->after.push_front(*tlpi) ;
	    (*lij)->before = tmp_tlp ; 
	  
	  }
	}
      }
    }
  }
}


void cat_tree::recursive_reorder(std::list<variableSet> &tmp_lvs) {
  for(std::list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
    (*li)->recursive_reorder(tmp_lvs) ;
  reorder(tmp_lvs) ;
  for(std::list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
    (*li)->recursive_reorder(tmp_lvs) ;
} 
 
std::vector<std::vector<variableSet> > create_orig_matrix(std::vector<variableSet> &ves) {
  std::vector<variableSet>::const_iterator vesi, viter ;
  std::vector<std::vector<variableSet> > vvvs(ves.size()) ;
  int tmp = 0 ;
  viter = ves.begin() ;
  while(tmp < ves.size()) { 
    entitySet vs ;
    for(vesi = ves.begin(); vesi != ves.end(); ++vesi) { 
      vs = (entitySet(*viter) & entitySet(*vesi)) ;
      vvvs[tmp].push_back(variableSet(vs));
    }
    tmp++ ;
    viter++ ;
  }
  /*
  if(Loci::MPI_rank == 0) {
    cout << " original_matrix =  " << endl ;
    for(int i = 0; i < vvvs.size(); ++i) {
      for(int j = 0; j < vvvs[i].size(); ++j) 
	cout << vvvs[i][j] << "   "  ; 
      cout << endl ;
    }
  }
  */
  return vvvs ;
}

std::vector<std::vector<variableSet> > create_sub_matrix(std::vector<variableSet> &ves, variableSet v) {
  std::vector<variableSet>::const_iterator vesi, viter ;
  std::vector<std::vector<variableSet> > vvvs ;
  std::vector<variableSet>  tmp_ves(ves.size()) ;
  int tmp = 0 ;
  for(int i = 0; i < ves.size(); ++i) 
    tmp_ves[i] = variableSet(entitySet(ves[i])- entitySet(v)) ;
  viter = tmp_ves.begin() ;
  while(tmp < tmp_ves.size()) { 
    entitySet vs ;
    std::vector<variableSet> tmp_vvs ;
    for(vesi = tmp_ves.begin(); vesi != tmp_ves.end(); ++vesi) {
      vs = (entitySet(*viter) & entitySet(*vesi)) ;
      tmp_vvs.push_back(variableSet(vs));
    }
    vvvs.push_back(tmp_vvs) ;
    tmp++ ;
    viter++ ;
  }
  /*
  if(Loci::MPI_rank == 0) {
    cout << " Sub matrix =  " << endl ;
    for(int i = 0; i < vvvs.size(); ++i) {
      for(int j = 0; j < vvvs[i].size(); ++j) 
	cout << vvvs[i][j] << "   "  ; 
      cout << endl ;
    }
  }
  */
  return vvvs ;
}

treeP recursive_matrix_tree(std::vector<std::vector<variableSet> > &tmp_vvvs, variableSet vset, std::set<variableSet> total_varset) {
  std::vector<std::vector<variableSet> > vvvs ;
  treeP tp = new cat_tree ;
  variableSet vs ;
  tp->generalization = vset ;
  tp->is_swapped = 0 ;
  if(tmp_vvvs.size() == 0) {
    return tp ;
  }
  std::vector<variableSet> tmp_vvs ;
  for(int i = 0; i < tmp_vvvs.size(); ++i) {
    if(tmp_vvvs[i][i] != EMPTY) {
      tmp_vvs.push_back(tmp_vvvs[i][i]) ;
    }
  }
  if(tmp_vvs.size() < tmp_vvvs.size())
    vvvs = create_orig_matrix(tmp_vvs) ;
  else
    vvvs = tmp_vvvs ;

  for(int i = 0; i < vvvs.size(); ++i) {
    std::set<variableSet> set_vars ;
    std::vector<variableSet> vvs ;
    entitySet es ;
    es = entitySet(vvvs[i][i]) ;
    if(set_vars.find(vvvs[i][i]) == set_vars.end()) {
      vvs.push_back(vvvs[i][i]) ;
      set_vars.insert(vvvs[i][i]) ;
    }
    std::vector<std::vector<variableSet> > sub_vvvs ;
    for(int j = 0; j < vvvs[i].size(); ++j)
      if(j != i) {
	if(vvvs[i][j] != EMPTY) {
	  es &= entitySet(vvvs[i][j]) ;
	  vs = variableSet(es) ;
	  if(vs != EMPTY) {
	    if(set_vars.find(vvvs[j][j]) == set_vars.end()) {
	      vvs.push_back(vvvs[j][j]) ;
	      set_vars.insert(vvvs[j][j]) ;
	    }
	  }
	}
      }
    if(vvs.size()) {
      es = vvs[0] ;
      for(int j = 1; j < vvs.size(); ++j)
	es &= entitySet(vvs[j]) ;
      if(total_varset.find(variableSet(es)) == total_varset.end()) {
	total_varset.insert(variableSet(es)) ;
	sub_vvvs = create_sub_matrix(vvs, variableSet(es)) ;
	tp->before.push_back(recursive_matrix_tree(sub_vvvs,
						   variableSet(es), total_varset));
      }
    }
    if(set_vars.size() == vvvs.size() )
      break ;
  }
  return tp ;
}

std::vector<entitySet> modified_categories(fact_db &facts, std::map<variable, entitySet> &vm, std::vector<interval> &pvec) {
  entitySet::const_iterator ei, ci ;
  std::vector<entitySet> tmp_pvec ;
  std::vector<variableSet> vvs(pvec.size()) ;
  variableSet vars = facts.get_typed_variables() ;
  std::map<variableSet, entitySet> mve ;
  std::map<variableSet, entitySet>::const_iterator miter ;
  std::map<variable, entitySet>::const_iterator svi ;
  variableSet initial_varset ;
  Loci::debugout << " Size of the vector passed to modified categories = " << pvec.size() << endl ;
  for(int i = 0; i < pvec.size(); ++i) {
    for(svi = vm.begin(); svi != vm.end(); ++svi) {
      if((svi->second & entitySet(pvec[i])) != EMPTY) {
	vvs[i] += svi->first ;
	initial_varset += svi->first ;
      }
    }
    if(vvs[i] != EMPTY) {
      mve[vvs[i]] += pvec[i]; 
    }
  }
  
  std::vector<variableSet> tmp_vvs ;
  for(miter = mve.begin(); miter != mve.end(); ++miter) {
    tmp_vvs.push_back(miter->first) ;
    Loci::debugout << "******************************************************" << endl ;
    Loci::debugout << " grouping variables " << miter->first << endl ;
    //Loci::debugout << " Entities shared = " << miter->second << endl ;
    Loci::debugout << " Total Entities causing the grouping = " << miter->second.size() << endl ;  
    Loci::debugout << "******************************************************" << endl ;
  }
  
  Loci::debugout << " The number of variable sets grouped due to common categories = " << tmp_vvs.size() << endl ;
  std::vector<std::vector<variableSet> > vvvs = create_orig_matrix(tmp_vvs) ;
  vvs.clear() ;
  std::vector<variableSet> indep ;
  std::set<variableSet> total_vars, itotal_vars ;
  
  for(int i = 0; i < vvvs.size(); ++i) {
    entitySet vs, tmp_vs ;
    for(int j = 0; j < vvvs[i].size(); ++j) {
      if(j != i) {
	vs = entitySet(vvvs[i][i]) & entitySet(vvvs[i][j]) ;
	tmp_vs += vs ;
	if(vs != EMPTY) {
	  if(total_vars.find(vvvs[j][j]) == total_vars.end()) {
	    vvs.push_back(vvvs[j][j]) ;
	    total_vars.insert(vvvs[j][j]) ;
	  }
	}
      }
    }
    if(tmp_vs == EMPTY) 
      indep.push_back(vvvs[i][i]) ;
  }
  std::vector<std::vector<variableSet> > tmp_vvvs ;
  vvvs = create_orig_matrix(vvs) ;
  for(int i = 0; i < vvvs.size(); ++i) {
    vvs.clear() ;
    entitySet es = entitySet(vvvs[i][i]);
    for(int j = 0; j < vvvs[i].size(); ++j) {
      entitySet tmp_es ;
      if((vvvs[i][j] != EMPTY))  {
	if(itotal_vars.find(vvvs[j][j]) == itotal_vars.end()) {
	  tmp_es = es & entitySet(vvvs[i][j]) ;
	  if(tmp_es != EMPTY) {
	    es = tmp_es ;
	    itotal_vars.insert(vvvs[j][j]) ;
	    vvs.push_back(vvvs[j][j]) ;
	  }
	}
	else
	  j++ ;
      }
    }
    if(vvs.size())
      tmp_vvvs.push_back(vvs) ;
    if(itotal_vars.size() == vvvs.size())
      break ; 
  }
  treeP super_tp = new cat_tree ;
  variableSet vset ;
  super_tp->generalization = vset ;
  for(int i = 0; i < tmp_vvvs.size(); ++i) {
    treeP tp = new cat_tree ;
    variableSet vs ;
    vvvs = create_orig_matrix(tmp_vvvs[i]) ;
    std::set<variableSet> total_varset ;
    tp = recursive_matrix_tree(vvvs, vs, total_varset) ;
    super_tp->before.push_back(tp) ;
  }
  treeP tp = new cat_tree ;
  tp->generalization = vset ;
  for(int i = 0 ; i < indep.size(); ++i) {
    treeP tmp_tp = new cat_tree ;
    tmp_tp->generalization = indep[i] ;
    tp->before.push_back(tmp_tp) ;
  } 
  for(std::list<treeP>::iterator li = super_tp->before.begin(); li != super_tp->before.end(); ++li) 
    for(std::list<treeP>::iterator lj = (*li)->before.begin(); lj != (*li)->before.end(); ++lj)
      tp->before.push_back(*lj) ;
  
  std::list<variableSet> lvs = tp->get_categories(vset) ; 
  std::list<variableSet> tlist ;
  tp->recursive_reorder(tlist) ;
  
  lvs = tp->get_categories(vset) ;
  entitySet total_entities ;
  HASH_MAP(int, entitySet) map_vec ;
  std::vector<int> sort_vec ;
  for(std::list<variableSet>::iterator lvi = lvs.begin(); lvi != lvs.end(); ++lvi) {
    if(mve.find(*lvi) != mve.end()) {
      entitySet tmp = mve[*lvi] ; 
      if(tmp!= EMPTY) {
	//map_vec[*ei] = tmp ;
	//sort_vec.push_back(*ei) ;
	//Loci::debugout << " Correspoding to " << *lvi << " pushing back " << tmp << endl ;
	tmp -= total_entities ;
	if(tmp!= EMPTY) {
	  tmp_pvec.push_back(tmp) ;
	  //ei = tmp.begin() ;
	  //map_vec[*ei] = tmp ;
	  //sort_vec.push_back(*ei) ;
	}
	total_entities += tmp ;
      }
    }
  }
  
  //for(int i = 0; i < tmp_pvec.size(); ++i)
  //Loci::debugout << " tmp_pvec[" << i << " ] = " << tmp_pvec[i] << endl ;
  //std::vector<int>::const_iterator svci = sort_vec.begin();
  /*
    while(svci != sort_vec.end()) {
    entitySet tmp = map_vec[*svci] ; 
    ++svci ;
    while(svci != sort_vec.end()) 
    if(tmp.Max() > map_vec[*svci].Min()) {
    tmp += map_vec[*svci] ;
    ++svci ;
    }
    else 
    break ;
    tmp_pvec.push_back(tmp) ;
    }
  */
  //for(std::vector<int>::const_iterator svci = sort_vec.begin(); svci != sort_vec.end(); ++svci)
  //tmp_pvec.push_back(map_vec[*svci]) ;
  Loci::debugout << " total_entities = " << total_entities << endl ;
  
  for(variableSet::const_iterator vsi = vars.begin(); vsi != vars.end(); ++vsi) {
    storeRepP p = facts.get_variable(*vsi) ;
    if((p->RepType() == MAP)) {
      MapRepP mp = MapRepP(p->getRep()) ;
      entitySet image_dom = mp->image(p->domain()) ;
      entitySet out_of_set = image_dom - total_entities ;
      entitySet left_out_categories = all_collect_entitySet(out_of_set) ;
      if(left_out_categories != EMPTY) {
	Loci::debugout << " left out stuff  = " << left_out_categories  << endl ;
	//ei = left_out_categories.begin() ;
	//map_vec[*ei] = left_out_categories ;
	//sort_vec.push_back(*ei) ;
	tmp_pvec.push_back(left_out_categories) ;
	total_entities += left_out_categories ;
      }
    }
  }
  //std::sort(sort_vec.begin(), sort_vec.end()) ; 
  //for(std::vector<int>::const_iterator svci = sort_vec.begin(); svci != sort_vec.end(); ++svci)
  //tmp_pvec.push_back(map_vec[*svci]) ;
  
  return tmp_pvec ;
}
  
  
  int GLOBAL_OR(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LOR,MPI_COMM_WORLD) ;
    return result ;
  }
  
  int GLOBAL_AND(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LAND,MPI_COMM_WORLD) ;
    return result ;
  }
  
  int GLOBAL_MAX(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD) ;
    return result ;
  }
  int GLOBAL_MIN(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD) ;
    return result ;
   }
  
} 

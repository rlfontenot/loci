#include "dist_tools.h"
#include "loci_globs.h"
#include <Tools/debug.h>
#include <entitySet.h>

#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <unistd.h>

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
using std::cout ;
#include <algorithm>
using std::swap ;
using std::sort ;
#ifdef SCATTER_DIST
#define UNITY_MAPPING
#endif

#include <Tools/stream.h>

namespace Loci {
  int MPI_processes = 1;
  int MPI_rank ;
  int num_threads = 1 ;
  int method = 3 ;
  /////////////////////////////
  // flags to turn on/off the visualization feature
  bool show_graphs = false ;
  bool show_decoration = false ;
  // flag to enable/disable the dynamic memory management
  bool use_dynamic_memory = false ;
  // flag to enable/disable output of dynamic memory
  // and multilevel graph decoration information
  bool show_dmm_verbose = false ;
  // flag to enable/disable chomping
  bool use_chomp = false ;
  // flag to visualize the chomping graph
  bool show_chomp = false ;
  // flag to turn on the summary report on chomping
  bool chomp_verbose = false ;
  // flag to enable outputing schedule to file
  bool schedule_output = false ;
  // flag to enable dynamic scheduling feature
  bool use_dynamic_scheduling = false ;
  // chomping size(unit is KB), default is 128KB
  int chomping_size = 128 ;
  /////////////////////////////
  
  ofstream debugout ;

  double total_memory_usage = 0 ;

  extern int current_rule_id ;

  void debug_print_rule() {
    if(current_rule_id != 0) {
      rule r(current_rule_id) ;
      cerr << "crash occured in rule " << r << endl ;

      if(exec_current_fact_db != 0) {
        char buf[512] ;
        sprintf(buf,"crash_dump.%d",MPI_rank) ;
        ofstream cfile(buf,ios::out) ;
        cfile << "rule: " << r << endl ;
        
        variableSet v = r.sources() ;

        cfile << "facts = {" << endl ;
        for(variableSet::const_iterator vi=v.begin();vi!=v.end();++vi) {
          cfile << *vi << ':' ;
          storeRepP p = Loci::exec_current_fact_db->get_variable(*vi);
          p->Print(cfile) ;
        }
        cfile << "}" << endl ;
      }
    }
  }
  //This is the first call to be made for any Loci program be it
  //sequential or parallel. 
  void Init(int* argc, char*** argv)  {
    char *execname = (*argv)[0] ;
    const char *hostname = "localhost" ;
    const char *debug = "gdb" ;
    //Setting up of the global variables for processor ID and the
    //total number of processes.  
    MPI_Init(argc, argv) ;
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN) ;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
    //Create a debug file for each process
    ostringstream oss ;
    // if output directory doesn't exist, create one
    bool debug_is_directory = true ;
    struct stat statbuf ;
    if(GLOBAL_OR(stat("debug",&statbuf)!=0)) {
      if(MPI_rank == 0)
        mkdir("debug",0755) ;
      for(int i=0;i<1000;++i) {
        if(GLOBAL_AND(stat("debug",&statbuf)==0))
          break ;
      }
    } else {
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file 'debug' should be a directory!, rename 'output' and start again."
             << endl ;
        debug_is_directory = false ;

      }
    }

    if(debug_is_directory) {
      if(MPI_processes == 1)
        oss << "debug/debug" ;
      else
        oss << "debug/debug."<<MPI_rank ;
    } else
      oss << "debug."<< MPI_rank ;
    
    string filename  = oss.str() ;
    debugout.open(filename.c_str(),ios::out) ;
    // All the rules in an unnamed namespace are first copied into the 
    // global rule list. To add rules to the rule database we just
    // neeed to use the global_rule_list. Inititally when the rules
    // are registered using the register rule it gets pushed into the
    // register_rule_list which is a static rule list. 
    if(!register_rule_list.empty()) {
      global_rule_list.copy_rule_list(register_rule_list) ;
      register_rule_list.clear() ;
    }
    bool debug_setup = false ;
    int i = 1 ;
    while(i<*argc) {
      if(!strcmp((*argv)[i],"--display")) {
        debug_setup = true ;
        hostname = (*argv)[i+1] ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--debug")) {
        debug = (*argv)[i+1] ;
        schedule_output = true ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--threads")) {
        cerr << "warning --threads not yet implemented" << endl ;
        num_threads = atoi((*argv)[i+1]) ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--scheduleoutput")) {
        schedule_output = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--graphs")) {
        show_graphs = true ; // visualize the dependency graph &
                             // the decomposed graph & every supernode
        i++ ;
      } else if(!strcmp((*argv)[i],"--decoration")) {
        show_decoration = true ; // visualize the decorated multilevel graph
        i++ ;
      } else if(!strcmp((*argv)[i],"--dmm")) {
        use_dynamic_memory = true ; // use the dynamic memory management
        i++ ;
      } else if(!strcmp((*argv)[i],"--dmmverbose")) {
        // output some info about dmm
        show_dmm_verbose = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chomp")) {
        // use the chomping scheme
        use_chomp = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--showchomp")) {
        // visualize the chomp graph
        show_chomp = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chompverbose")) {
        // summary report
        chomp_verbose = true ;
        i++ ;
      } else if(!strcmp((*argv)[i],"--chompingsize")) {
        std::stringstream ss ;
        ss << (*argv)[i+1] ;        
        ss >> chomping_size ;
        i+=2 ;
      } else if(!strcmp((*argv)[i],"--method")) {
        method = atoi((*argv)[i+1]);
        i+=2;
      } else if(!strcmp((*argv)[i],"--balance")) {
        use_dynamic_scheduling = true ;
        i++ ;
      } else
        break ;
    }
    if(i!=1) {
      *argc -= (i-1) ;
      for(int k=1;k<*argc;++k)
        (*argv)[k] = (*argv)[k+i-1] ;
    }
    
    set_debug_callback(debug_print_rule) ;
    if(debug_setup) {
      setup_debugger(execname,debug,hostname) ;
    }
    chopsigs_() ;
  }
  //All Loci programs must end with this call. 
  void Finalize() {
    MPI_Finalize() ;
  }

  void Abort() {
    debugger_() ;
  }

  void get_clone(fact_db &facts, const rule_db &rdb) {
    fact_db::distribute_infoP df = facts.get_distribute_info()  ;
    std::vector<entitySet> &ptn = facts.get_init_ptn() ;
    entitySet bdom = ptn[Loci::MPI_rank] & interval(Loci::UNIVERSE_MIN, -1) ;
    entitySet global_bdom = Loci::all_collect_entitySet(bdom) ;
    int p = 0; 
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      entitySet tmp = ptn[i] ; 
      ptn[i] = tmp & interval(0, Loci::UNIVERSE_MAX) ;
    }
    FORALL(global_bdom, i) {
      int tmp = p % Loci::MPI_processes ;
      ptn[tmp] += i ;
      p++ ;
    } ENDFORALL ;
    variableSet tmp_vars = facts.get_typed_variables();
    for(variableSet::const_iterator vi = tmp_vars.begin(); vi != tmp_vars.end(); ++vi) {
      Loci::storeRepP tmp_sp = facts.get_variable(*vi) ;
      if(tmp_sp->RepType() == Loci::CONSTRAINT) {
        entitySet tmp_dom = tmp_sp->domain() ;
        if(tmp_dom != ~EMPTY) {
          entitySet global_tmp_dom = Loci::all_collect_entitySet(tmp_dom) ;
          constraint tmp ;
          *tmp = global_tmp_dom ;
          facts.update_fact(variable(*vi), tmp) ; 
        }
      }
      if(tmp_sp->RepType() == Loci::MAP) {
        storeRepP map_sp = Loci::MapRepP(tmp_sp->getRep())->thaw() ; 
        facts.replace_fact(*vi, map_sp) ;
      }
    }
    entitySet tmp_copy, image ;
    std::set<std::vector<variableSet> > maps ;

    entitySet::const_iterator ei ;
    Loci::get_mappings(rdb,facts,maps) ;
    std::vector<entitySet> copy(MPI_processes), send_clone(MPI_processes) ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes];
    int *recv_displacement = new int[MPI_processes];
    int size_send = 0 ;
    entitySet tmp_set = ptn[Loci::MPI_rank] ; 
    image = Loci::dist_expand_map(tmp_set, facts, maps) ;
    tmp_copy =  image - ptn[MPI_rank] ; 
    for(int i = 0; i < MPI_processes; ++i) {
      copy[i] = tmp_copy & ptn[i] ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
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
      for(ei = copy[i].begin(); ei != copy[i].end(); ++ei) {
        send_buf[size_send] = *ei ;
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
            recv_displacement[i]+recv_count[i]; ++j) 
        send_clone[i] += recv_buf[j] ;
    }
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
  
    variableSet vars = facts.get_typed_variables() ;
    double start = MPI_Wtime() ;
    int myid = MPI_rank ;
    int size = 0 ;
    Map l2g ;
    constraint my_entities ;
    int isDistributed ;
    std::vector<entitySet> iv ; 
    entitySet::const_iterator ti ;
    vector<entitySet> proc_entities ;
    Loci::categories(facts,iv) ;
    entitySet e ;
    if(Loci::MPI_rank == 0)
      cout << " initial_categories =  " << iv.size() << endl ;
  
    for(size_t i = 0; i < iv.size(); ++i) {
      //    Loci::debugout << " iv[ " << i << " ] = " << iv[i] << endl ;
      // Within each category:
      // 1) Number local processor entities first
      e = ptn[myid] & iv[i] ; 
      if(e != EMPTY) {
        proc_entities.push_back(e) ; 
        size += e.size() ;
      }
      // 2) Number clone region entities next
      for(int j = 0; j < MPI_processes; ++j) 
        if(myid != j) {
          e = copy[j] & iv[i];
          if(e != EMPTY) {
            proc_entities.push_back(e) ;
            size += e.size() ;
          }
        }
    }
    iv.clear() ;
    entitySet g ; 
    //#define UNITY_MAPPING
#ifdef UNITY_MAPPING
    cout << "Using Unity Mapping " << endl ;
    for(int i=0;i<proc_entities.size();++i)
      g+= proc_entities[i] ;
    l2g.allocate(g) ;
    for(entitySet::const_iterator ei=g.begin();ei!=g.end();++ei)
      l2g[*ei] = *ei ;
#else
    MPI_Barrier(MPI_COMM_WORLD) ;
    int j = 0 ;
    e = interval(0, size-1) ;
    l2g.allocate(e) ;
    for(size_t i = 0; i < proc_entities.size(); ++i) {
      g += proc_entities[i] ;
      for(ei = proc_entities[i].begin(); ei != proc_entities[i].end(); ++ei ) {
        l2g[j] = *ei ;
        ++j ;
      }
    }
    proc_entities.clear() ;
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
    for(int i = 0 ; i < MPI_processes; ++i) 
      if(myid != i )
        if(copy[i] != EMPTY) 
          recv_neighbour += i ; 
  
    for(int i = 0; i < MPI_processes; ++i)
      if(myid != i)
        if(send_clone[i] != EMPTY)
          send_neighbour += i ;
  
    send_entities.allocate(send_neighbour) ;
    recv_entities.allocate(recv_neighbour) ;
    for(ei = recv_neighbour.begin(); ei != recv_neighbour.end(); ++ei) {
      for(ti =  copy[*ei].begin(); ti != copy[*ei].end(); ++ti) {
        recv_entities[*ei] += df->g2l[*ti] ;
      }
    }
  
    for(ei = send_neighbour.begin(); ei!= send_neighbour.end(); ++ei) {
      for(ti =  send_clone[*ei].begin(); ti != send_clone[*ei].end(); ++ti)
        send_entities[*ei] +=  df->g2l[*ti] ;
    }
    double end_time =  MPI_Wtime() ;
    Loci::debugout << "  Time taken for creating intitial info =  " << end_time - start << endl ;
    start = MPI_Wtime() ;
    Loci::reorder_facts(facts, df->g2l) ;
    end_time =  MPI_Wtime() ;
    Loci::debugout << "  Time taken for reordering =  " << end_time - start << endl ; 
    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    g = EMPTY ;
    for(ei = ptn[myid].begin(); ei != ptn[myid].end(); ++ei)
      g += df->g2l[*ei] ;
    my_entities = g ;
    df->myid = myid ;
    df->my_entities = g ;
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
    for(size_t i=0;i<df->xmit.size();++i)
      total += df->xmit[i].size ;
    df->xmit_total_size = total ;
    total = 0 ;
    for(size_t i=0;i<df->copy.size();++i)
      total += df->copy[i].size ;
    df->copy_total_size = total ;
    facts.put_distribute_info(df) ;
    facts.create_fact("l2g", l2g) ;
    facts.put_l2g(l2g) ;
    facts.create_fact("my_entities",my_entities);
  }
  
  //This routine was first written to read in a partition of
  //entities. No longer needed if we are relying completely on the
  //scalable version.  
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
  //This routine writes out a generalized partition of entities. It
  //calls the generalized non scalable partitioning routine and writes 
  //out p partitions . 
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
      for(unsigned int i = 0; i < ptn.size(); ++i) 
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
    for(int i = 0; i < size_map; ++i) {
      xadj[i+1] = xadj[i] + size_adj[i] ;
    }
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

  void get_mappings(const rule_db &rdb, fact_db &facts,
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
	  for(unsigned int i = 0; i < vmsi->mapping.size(); ++i) {
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
	  for(unsigned int i = 0; i < vmsi->mapping.size(); ++i) {
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
          for(unsigned int i = 0; i < vmsi->mapping.size(); i++) {
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
      for(unsigned int i=0;i<vss.size();++i)
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
            for(unsigned int i=2;i<vss.size();++i)
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
  associated with a particular process. This routine doesn't need to
  perform any communication as the whole map is present on all the
  processors. */
  
  entitySet expand_map(entitySet domain, fact_db &facts,
                       const set<vector<variableSet> > &maps) {
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    set<vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      const vector<variableSet> &mv = *smi ;
      for(unsigned int i = 0; i < mv.size(); ++i) {
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
  //This routine  is similar to the expand map but it works for maps
  //which are distributed across processors. 
  entitySet dist_expand_map(entitySet domain, fact_db &facts,
			    const std::set<std::vector<variableSet> > &maps) {   
    std::vector<entitySet> ptn = facts.get_init_ptn() ;
    
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = ptn[i] ;
      ptn[i] = tmp & interval(0, UNIVERSE_MAX) ;
    }
    entitySet dom = domain ;
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet locdom = domain ;
      const vector<variableSet> &mv = *smi ;
      for(unsigned int i = 0; i < mv.size(); ++i) {
	variableSet v = mv[i] ;
	v &= vars ; 
	entitySet image ;
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    entitySet tmp_dom = p->domain() ;
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	    glob_dom &= dom ;
	    entitySet tmp_out = (glob_dom & locdom) - tmp_dom ; 
	    //Loci::debugout << " variable = " << *vi << " tmp_out = " << tmp_out << endl ;
	    storeRepP sp = mp->expand(tmp_out, ptn) ;
	    if(sp->domain() != tmp_dom) {
	      facts.update_fact(variable(*vi), sp) ; 
	      //Loci::debugout << "updated map = " << *vi <<  endl ;
	      //sp->Print(Loci::debugout) ;
	    }
	    image +=  MapRepP(sp)->image((sp->domain()) & locdom) ;
	    dom += tmp_out ;
	  }
	}
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    entitySet tmp_dom = p->domain() ;
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	    glob_dom &= dom ;
	    entitySet tmp_out = (glob_dom & locdom) - tmp_dom ; 
	    storeRepP sp = mp->expand(tmp_out, ptn) ;
	    if(sp->domain() != tmp_dom) {
	      facts.update_fact(variable(*vi), sp) ; 
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
  
  // This routine does a generalized partitioning of the entities in
  // the fact database and returns the partition. 
  vector<entitySet> generate_distribution(fact_db &facts, rule_db &rdb, int num_partitions) {
    if(num_partitions == 0)
      num_partitions = MPI_processes ;
    vector<entitySet> ptn ;
    if(num_partitions == 1)
      return ptn ;
    
    debugout << "Synchronising before metis_facts" << endl ;
    double start = MPI_Wtime() ;

    metis_facts(facts,ptn,num_partitions) ;

    double end_time  = MPI_Wtime() ;
    debugout << "Time taken for metis_facts = " << end_time -start << endl ;
    fact_db::distribute_infoP df = new fact_db::distribute_info  ;
    facts.put_init_ptn(ptn) ;
    facts.put_distribute_info(df) ;
    return ptn ;
  }
  //This  routine is almost similar to the above one but for the
  //introduction of the chopped partitioning and the remapping of
  //entities needed for scalable I/O
  vector<entitySet> generate_scalable_distribution(fact_db &facts, rule_db &rdb, int num_partitions) {
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
    std::vector<entitySet> partition(Loci::MPI_processes) ;
    std::vector<entitySet> chop_ptn(Loci::MPI_processes) ;
    dMap remap ;
    for(int i = 0 ; i < Loci::MPI_processes; ++i) {
      entitySet tmp_set = ptn[i] ;
      if(Loci::MPI_rank == i) {
	FORALL(tmp_set, ei) {
	  remap[ei] = ei ;
	} ENDFORALL ;
      }
    }
    int indx = 1 ;
    // This is needed for the new scalable I/O routines. 
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int ivl = ptn[i].size() ;
      chop_ptn[i] = Loci::interval(indx, indx + ivl-1) ;
      indx += ivl ;
    }
    fact_db::distribute_infoP df = new fact_db::distribute_info  ;
    //Remapping of entities keeps the partitioning contiguous . 
    df->remap = remap ;
    df->chop_ptn = chop_ptn ;
    facts.put_init_ptn(ptn) ;
    facts.put_distribute_info(df) ;
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
    start = MPI_Wtime() ;
    int num_procs = MPI_processes ;
    vector<entitySet> copy(num_procs) ;
    vector<entitySet> image(num_procs) ;
    entitySet tmp ;
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
    end_time  = MPI_Wtime() ;
    debugout << "Time taken for all the calls to expand_map is =   = " << end_time-start << endl ; 
    start = MPI_Wtime() ;
    int myid = MPI_rank ;
    int size = 0 ;
    fact_db::distribute_infoP df = facts.get_distribute_info()  ;
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
    for(unsigned int i = 0; i < iv.size(); ++i) {
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
    for(unsigned int i = 0; i < proc_entities.size(); ++i) {
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
    for(unsigned int i=0;i<df->xmit.size();++i)
      total += df->xmit[i].size ;
    df->xmit_total_size = total ;
    
    total = 0 ;
    for(unsigned int i=0;i<df->copy.size();++i)
      total += df->copy[i].size ;
    df->copy_total_size = total ;
    
    
    facts.put_distribute_info(df) ;
    facts.create_fact("l2g", l2g) ;
    facts.put_l2g(l2g) ;
    facts.create_fact("my_entities", my_entities) ;
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for creating final info =  " << end_time - start << endl ;
  }
  
  /*The fill_entitySet routine fills in the clone region entities
    . The send_buffer and the recv_buffer are allocated only once to
    contain the maximum clone region */ 
  entitySet fill_entitySet(const entitySet& e, fact_db &facts) {
    
    entitySet re ;
    
    if(facts.isDistributed()) {  
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;
      
      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;
        recv_buffer[0] = new int[d->copy_total_size] ;
        recv_size[0] = d->copy[0].size ;
        for(unsigned int i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->copy[i-1].size ;
          recv_size[i] = d->copy[i].size ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size] ;
        for(unsigned int i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->xmit[i-1].size ;
      }
      
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;
      /*The recv_size is the maximum possible, so that even if we
	receive a shorter message there won't be any problem */
      for(unsigned int i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(unsigned int i=0;i<d->xmit.size();++i) {
        entitySet temp = e & d->xmit[i].entities ;
	
        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;


	int send_size = temp.size() ;
        MPI_Send(send_buffer[i], send_size, MPI_INT, d->xmit[i].proc,
                 1, MPI_COMM_WORLD) ;
      }
        

      if(d->copy.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }
      for(unsigned int i = 0; i < d->copy.size(); ++i) {
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
      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;

        recv_buffer[0] = new int[d->copy_total_size*evsz+evsz*d->copy.size()] ;
        recv_size[0] = d->copy[0].size*evsz+evsz ;
        for(unsigned int i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->copy[i-1].size*evsz+evsz) ;
          recv_size[i] = d->copy[i].size*evsz+evsz ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        for(unsigned int i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
      }
        
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;

      for(unsigned int i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(unsigned int i=0;i<d->xmit.size();++i) {
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
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }

      for(unsigned int i = 0; i < d->copy.size(); ++i) {
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

      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[d->xmit_total_size] ;
        recv_size[0] = d->xmit[0].size ;

        for(unsigned int i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->xmit[i-1].size ;
          recv_size[i] = d->xmit[i].size ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size] ;
        for(unsigned int i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(unsigned int i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      /*By intersecting the given entitySet with the clone region
	entities we can find out which entities are to be sent */
      for(unsigned int i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei)
          send_buffer[i][j++] = l2g[*ei] ;

        int send_size = temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }
      
      if(d->xmit.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }
      
      for(unsigned int i=0;i<d->xmit.size();++i) {
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
      int **send_buffer = 0 ; 
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        recv_size[0] = d->xmit[0].size*evsz + evsz ;

        for(unsigned int i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
          recv_size[i] = d->xmit[i].size*evsz+evsz ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size*evsz+evsz*d->copy.size()] ;
        for(unsigned int i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size*evsz+evsz ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(unsigned int i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      for(unsigned int i=0;i<d->copy.size();++i) {
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
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
      	FATAL(err != MPI_SUCCESS) ;
      }
      for(unsigned int i=0;i<d->xmit.size();++i) {
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
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
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

#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
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
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
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

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
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

#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
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

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(int k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_size_bytes[k],1,MPI_INT, k+1,14,
		    MPI_COMM_WORLD, &size_request[k]) ;  

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
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
	recv_ptr[0] = new unsigned char[sz] ;
	for(int i = 1; i < MPI_processes-1; i++)
	  recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
	
	for(int i = 0; i < MPI_processes-1; i++)
	  MPI_Irecv(recv_ptr[i],r_size[i] , MPI_PACKED, i+1, 13,
		    MPI_COMM_WORLD, &store_request[i]) ;

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, store_request, store_status) ;
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
    return nsp ;
    
  }
  //This routine is used to collect the entire store on processor
  //0. The store is initially distributed across all the processors in 
  //their initial global numbering. This is not a scalable option for
  //large containers. Written for use in the couple_grids routine in
  //fvm_coupler.cc  
  storeRepP collect_global_store(storeRepP &sp) {
    storeRepP nsp = sp ;
    if(Loci::MPI_rank == 0) {
      std::vector<sequence> vseq(MPI_processes-1) ;
      MPI_Status *status, *size_status, *store_status ;
      MPI_Request *recv_request, *size_request, *store_request ;
      int **recv_buffer ;
      int *recv_size ;
      entitySet re ;
      entitySet temp = sp->domain() ;
      sequence te = sequence(temp) ;
      re += temp ;
      recv_size = new int[MPI_processes-1] ;
      int *recv_size_bytes = new int[MPI_processes-1] ;
      size_request = new MPI_Request[MPI_processes-1] ;
      size_status = new MPI_Status[MPI_processes-1] ;
      for(int k = 0; k < MPI_processes-1; k++) 
	MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,11, MPI_COMM_WORLD, &size_request[k]);
      
#ifdef DEBUG
      int err =
#endif 
	MPI_Waitall(MPI_processes-1, size_request, size_status) ;
      FATAL(err != MPI_SUCCESS) ;
      
      recv_buffer = new int*[MPI_processes-1] ;
      int recv_size_total = 0 ;
      for(int k=0;k<MPI_processes-1;++k) 
	recv_size_total += recv_size[k] ;
      
      recv_buffer[0] = new int[recv_size_total] ;
      for(int i = 1; i < MPI_processes-1; ++i)
	recv_buffer[i] = recv_buffer[i-1] + recv_size[i-1] ;
      recv_request = new MPI_Request[MPI_processes-1] ;
      status = new MPI_Status[MPI_processes-1] ;
      
      for(int k = 0; k < MPI_processes-1; k++)
	MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,12, MPI_COMM_WORLD, &recv_request[k] );  
      
#ifdef DEBUG
      err = 
#endif 
	MPI_Waitall(MPI_processes-1, recv_request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int k = 0; k < MPI_processes-1; k++) 
	MPI_Irecv(&recv_size_bytes[k],1,MPI_INT, k+1,14,
	 	  MPI_COMM_WORLD, &size_request[k]) ;  
      
#ifdef DEBUG
      err =
#endif 
 	MPI_Waitall(MPI_processes-1, size_request, size_status) ;
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
      recv_ptr[0] = new unsigned char[sz] ;
      for(int i = 1; i < MPI_processes-1; i++)
	recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
      
      for(int i = 0; i < MPI_processes-1; i++)
	MPI_Irecv(recv_ptr[i],r_size[i] , MPI_PACKED, i+1, 13,
		  MPI_COMM_WORLD, &store_request[i]) ;
      
#ifdef DEBUG
      err =
#endif 
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
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
      entitySet temp = dom ;
      int send_size = temp.size() ;
      int *send_buffer = new int[send_size] ;
      int sz = sp->pack_size(temp) ;
      unsigned char *send_ptr = new unsigned char[sz] ;
      
      int j = 0 ;
      for(entitySet::const_iterator ti = temp.begin(); ti != temp.end(); ++ti) 
	send_buffer[j++] = *ti ; 
      
      int loc_pack = 0;
      sp->pack(send_ptr, loc_pack, sz, temp) ;
      
      MPI_Send(&send_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
      MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 12, MPI_COMM_WORLD) ;
      MPI_Send(&sz,1,MPI_INT,0,14,MPI_COMM_WORLD) ;
      MPI_Send(send_ptr, sz, MPI_PACKED, 0, 13, MPI_COMM_WORLD) ;
      delete [] send_buffer ;
      delete [] send_ptr ;
    } 
    return nsp ;
  }
  //This is a generalized routine for writing out storeRepP's. Has
  //been tested for stores, storeVecs and multiStores. The initial
  //store in the local numbering is first redistributed such that it
  //ends up in a blocked partitioning format in the global numbering
  //across all the processors. The qrep passed to this routine is in
  //the chunked partitioning format in global numbering.   
void write_container(hid_t group_id, storeRepP qrep) {
    entitySet dom = qrep->domain() ;
    entitySet q_dom = all_collect_entitySet(dom) ;
    unsigned char* tmp_send_buf ;
    std::vector<int> sort_max ;
    int local_size = qrep->pack_size(dom) ;
    sort_max = all_collect_sizes(local_size) ;
    std::vector<int> interval_sizes = all_collect_sizes(dom.size()) ; 
    int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
    tmp_send_buf = new unsigned char[total_size] ;
    if(Loci::MPI_rank == 0)
      Loci::HDF5_WriteDomain(group_id, q_dom);
    frame_info fi = qrep->write_frame_info(group_id) ;
    int array_size = 0 ;
    if(fi.size) 
      if(fi.is_stat)   
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	array_size = fi.size * dom.size() ;
    else
      if(fi.is_stat) 
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	  array_size += *fvi ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size) ;
    int tot_arr_size = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_arr_size += arr_sizes[i] ;
    if(Loci::MPI_rank != 0) {
      MPI_Status status ;
      int send_size_buf ;
      send_size_buf = qrep->pack_size(dom) ;
      int tot_size = send_size_buf ;
      int loc_pack = 0 ;
      qrep->pack(tmp_send_buf, loc_pack, total_size, dom) ;
      int flag = 0 ;
      MPI_Recv(&flag,1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status) ;
      if(flag) {
	MPI_Send(&tot_size, 1, MPI_INT, 0, 11, MPI_COMM_WORLD) ;
	MPI_Send(tmp_send_buf, tot_size, MPI_PACKED, 0, 12, MPI_COMM_WORLD) ;
      }
    } else {
      int rank = 1 ;
      hsize_t dimension = 1 ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = arr_sizes[0] ;
      dimension =  tot_arr_size ;
      hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      dimension = count ;
      start += dimension ;
      hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace, H5P_DEFAULT) ;
      qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
      H5Dclose(dataset) ;
      
      int curr_indx = dom.size() ;
      for(int i = 1; i < Loci::MPI_processes; ++i) {
	MPI_Status status ;
	int recv_total_size ;
	entitySet tmpset = entitySet(interval(curr_indx, interval_sizes[i]+curr_indx-1)) ;
	curr_indx += interval_sizes[i] ;
	Loci::storeRepP t_qrep = qrep->new_store(tmpset) ;
	int loc_unpack = 0 ;
	int flag = 1 ;
	MPI_Send(&flag, 1, MPI_INT, i, 10, MPI_COMM_WORLD) ;
	MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, MPI_COMM_WORLD, &status) ;
	MPI_Recv(tmp_send_buf, recv_total_size, MPI_PACKED, i, 12, MPI_COMM_WORLD, &status) ;
	Loci::sequence tmp_seq = Loci::sequence(tmpset) ;
	t_qrep->unpack(tmp_send_buf, loc_unpack, total_size, tmp_seq) ;
	dimension = arr_sizes[i] ;
	count = dimension ; 
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
	start += count ;
	dataset = H5Dopen(group_id, "data") ;
	t_qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", tmpset) ;
        H5Dclose(dataset) ;
      }
      H5Sclose(dataspace) ;
    }
    delete [] tmp_send_buf ;
  }
  //This routine 
  void read_container(hid_t group_id, storeRepP qrep, entitySet &dom) {
    entitySet q_dom ;
    if(Loci::MPI_rank == 0)
      Loci::HDF5_ReadDomain(group_id, q_dom) ;
    unsigned char* tmp_buf = 0;
    std::vector<int> interval_sizes ;
    q_dom = all_collect_entitySet(q_dom) ;
    if(dom == EMPTY) {
      int sz = q_dom.size() / Loci::MPI_processes ;
      for(int i = 0; i < MPI_processes-1; ++i)
	interval_sizes.push_back(sz) ;
      interval_sizes.push_back(sz + (q_dom.size() % Loci::MPI_processes)) ;
      int tmp = q_dom.Min() ;
      int max = q_dom.Max() ;
      if(MPI_rank < (MPI_processes-1))
	dom = interval(tmp + MPI_rank*sz, tmp + (MPI_rank+1)*sz-1) ;
      else
	dom = interval(tmp + MPI_rank*sz, max) ;
    } else {
      interval_sizes = all_collect_sizes(dom.size()) ;
      entitySet tset = all_collect_entitySet(dom) ;
      if(tset != q_dom) {
	cerr << "The total domain of the container and the sum of domains across the processors does not match" << endl ;
	cerr << "q_dom = " << q_dom << endl ;
	cerr << "tset = " << tset << endl ;
      }
    }
    if(qrep->domain() == EMPTY)
      qrep->allocate(dom) ;
    FATAL(qrep->domain() != dom) ;
    frame_info fi = qrep->read_frame_info(group_id) ;
    int array_size = 0 ;
    int vec_size = 0 ;
    if(fi.size)
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
	vec_size = fi.second_level.size() ;
      }
      else {
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	array_size = fi.size * dom.size() ;
      }
    else { 
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi) 
	  array_size += *vi ;
	vec_size = fi.second_level.size() + dom.size() ;
      }
      else {
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi) 
	  array_size += *fvi ;
	vec_size = dom.size() ;
      }
    }
    std::vector<int> tmp_sizes = all_collect_sizes(vec_size) ;
    int max_tmp_size = *std::max_element(tmp_sizes.begin(), tmp_sizes.end()) ;
    int max_eset_size = *std::max_element(interval_sizes.begin(), interval_sizes.end()) ;
    int* tmp_int  ;
    tmp_int = new int[max_tmp_size] ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size) ;
    int tot_arr_size = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_arr_size += arr_sizes[i] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      int t = 0 ;
      if(fi.size) { 
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	if(fi.is_stat) 
	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	    tmp_int[t++] = *vi ;
      }
      else {
	if(fi.is_stat) {
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;
	  
	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi) 
	    tmp_int[t++] = *vi ;
	}
	else
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;
      }
      if(tmp_sizes[Loci::MPI_rank])
	MPI_Send(tmp_int, tmp_sizes[Loci::MPI_rank], MPI_INT, 0, 10, MPI_COMM_WORLD) ;
      int total_size = 0 ;
      MPI_Recv(&total_size, 1, MPI_INT, 0, 11,
	       MPI_COMM_WORLD, &status) ;  
      tmp_buf = new unsigned char[total_size] ;
      MPI_Recv(tmp_buf, total_size, MPI_PACKED, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      Loci::sequence tmp_seq = Loci::sequence(dom) ;
      int loc_unpack = 0 ;
      if(qrep->domain() == EMPTY)
	qrep->allocate(dom) ;
      FATAL(dom != qrep->domain()) ;
      qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
    } else {
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int curr_indx = 0 ;
      int total_size = 0 ;
      int tmp_total_size = 0 ;
      entitySet max_set = interval(0, max_eset_size-1) ;
      storeRepP tmp_sp ;
      if(fi.size) 
	tmp_sp = qrep->new_store(max_set) ; 
      
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	entitySet local_set = entitySet(interval(curr_indx, interval_sizes[p]+curr_indx-1)) ;
	curr_indx += interval_sizes[p] ;
	hsize_t dimension = arr_sizes[p] ;
	count = dimension ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	entitySet tmp_set = interval(0, local_set.size()-1) ;
	if(p && tmp_sizes[p]) {
	  MPI_Recv(tmp_int, tmp_sizes[p], MPI_INT, p, 10,
		   MPI_COMM_WORLD, &status) ;
	  std::vector<int> vint, fvint ;
	  int t = 0 ;
	  if(fi.size) {
	    if(fi.is_stat) {
	      for(int i = 0; i < tmp_sizes[p]; ++i)
		vint.push_back(tmp_int[t++]) ;
	      fi.second_level = vint ;
	    } 
	  }
	  else {
	    for(int i = 0; i < local_set.size(); ++i)
	      fvint.push_back(tmp_int[t++]) ;
	    for(int i = 0; i < tmp_sizes[p]-local_set.size(); ++i)
	      vint.push_back(tmp_int[t++]) ;
	    fi.first_level = fvint ;
	    fi.second_level = vint ;
	  }
	}
	storeRepP t_sp ;
	int t = 0 ;
	if(p == 0) 
	  if(!fi.size)
	    for(std::vector<int>::const_iterator vi = fi.first_level.begin(); vi != fi.first_level.end(); ++vi)
	      tmp_int[t++] = *vi ;
	if(fi.size) {
	  tmp_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ; 
	  tmp_total_size = tmp_sp->pack_size(tmp_set) ;
	}
	else {
	  t_sp = qrep->new_store(tmp_set, tmp_int) ;
	  t_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ; 
	  tmp_total_size = t_sp->pack_size(tmp_set) ;
	}
	if(tmp_total_size > total_size) {
	  total_size = tmp_total_size ;
	  if(p)
	    delete [] tmp_buf ;
	  tmp_buf = new unsigned char[total_size] ;
	}
	start += count ;
	int loc = 0 ;
	if(fi.size)
	  tmp_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
	else
	  t_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
	if(p == 0) {
	  int loc_unpack = 0 ;
	  Loci::sequence tmp_seq = Loci::sequence(dom) ;
	  if(fi.size) 
	    if(fi.size > 1) 
	      qrep->set_elem_size(fi.size) ;
	  if(qrep->domain() == EMPTY)
	    qrep->allocate(dom) ;
	  FATAL(dom != qrep->domain()) ;
	  qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
	} else { 
	  MPI_Send(&total_size, 1, MPI_INT, p, 11, MPI_COMM_WORLD) ;
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, p, 12, MPI_COMM_WORLD) ;
	}
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
    delete [] tmp_buf ;
    delete [] tmp_int ; 
  }
  void read_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size) {
    std::vector<int> vec_size = all_collect_sizes(dom_size) ;
    hsize_t dimension = 0 ;
    hid_t dataset = 0;
    hid_t dataspace = 0;
    if(Loci::MPI_rank == 0) {
      dataset = H5Dopen(group_id, name) ;
      dataspace = H5Dget_space(dataset) ;
      H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    }
    int dim = dimension ;
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    int rank = 1 ;
    hid_t datatype = H5T_NATIVE_INT ;
    int total_size = dim / MPI_processes ;
    std::vector<int> sizes = all_collect_sizes(dom_size) ;
    total_size = *std::max_element(sizes.begin(), sizes.end() );
    int *tmp_int = new int[total_size] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      MPI_Recv(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      for(int i = 0; i < sizes[MPI_rank]; ++i)
	vint.push_back(tmp_int[i]) ;
    } else { 
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	dimension = sizes[p] ;
	count = dimension ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_int) ;
        if(err < 0) {
          cerr << "H5Dread() failed" << endl ;
        }
	H5Sclose(memspace) ;
	start += count ;
	if(p == 0) {
	  for(int i = 0; i < sizes[p]; ++i) 
	    vint.push_back(tmp_int[i]) ;
	} else 
	  MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, MPI_COMM_WORLD) ;
      }
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
    }
    delete [] tmp_int ; 
  }
  void read_multi_vector_int(hid_t group_id, const char* name, int dim,  std::vector<int>& vint) {
    hsize_t dimension = 0 ;
    hid_t dataset = 0;
    hid_t dataspace = 0;
    if(Loci::MPI_rank == 0) {
      dataset = H5Dopen(group_id, name) ;
      dataspace = H5Dget_space(dataset) ;
      H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    }
    int rank = 1 ;
    hid_t datatype = H5T_NATIVE_INT ;
    int total_size = 0 ;
    std::vector<int> sizes = all_collect_sizes(dim) ;
    total_size = *std::max_element(sizes.begin(), sizes.end());
    int *tmp_int = new int[total_size] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      MPI_Recv(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      for(int i = 0; i < sizes[MPI_rank]; ++i)
	vint.push_back(tmp_int[i]) ;
    } else { 
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	dimension = sizes[p] ;
	count = dimension ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_int) ;
        if(err < 0) {
          cerr << "H5Dread() failed" << endl ;
        }
	H5Sclose(memspace) ;
	start += count ;
	if(p == 0) {
	  for(int i = 0; i < sizes[p]; ++i) 
	    vint.push_back(tmp_int[i]) ;
	} else 
	  MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, MPI_COMM_WORLD) ;
      }
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
    }
    delete [] tmp_int ; 
  }
  
  void write_vector_int(hid_t group_id, const char* name, std::vector<int>& vint) {
    std::vector<int> sort_max(MPI_processes) ;
    int tot_entities = vint.size() ;
    sort_max = all_collect_sizes(tot_entities) ;
    std::vector<int> sizes = sort_max ;
    std::sort(sort_max.begin(), sort_max.end()) ;
    tot_entities = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_entities += sort_max[i] ;
    int *tmp_int = new int[sort_max[MPI_processes-1]] ;
    int tmp = 0 ;
    for(std::vector<int>::iterator vi = vint.begin(); vi != vint.end(); ++vi)
      tmp_int[tmp++] = *vi ;
    if(Loci::MPI_rank != 0) {
      MPI_Status status ;
      int flag = 0 ;
      MPI_Recv(&flag,1, MPI_INT, 0, 11, MPI_COMM_WORLD, &status) ;
      if(flag)
	MPI_Send(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12, MPI_COMM_WORLD) ;
    } else {
      hid_t datatype = H5T_NATIVE_INT ;
      int rank = 1 ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      hsize_t dimension = tot_entities ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      count = sizes[0] ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
      dimension = sizes[0] ;
      start += dimension ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t dataset = H5Dcreate(group_id, name , datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
      H5Dclose(dataset) ;
      H5Sclose(memspace) ;
      for(int i = 1; i < Loci::MPI_processes; ++i) {
	MPI_Status status ;
	int flag = 1 ;
	MPI_Send(&flag, 1, MPI_INT, i, 11, MPI_COMM_WORLD) ;
	MPI_Recv(tmp_int, sizes[i], MPI_INT, i, 12, MPI_COMM_WORLD, &status) ;
	dimension = sizes[i] ;
	count = dimension ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
	start += count ;
	memspace = H5Screate_simple(rank, &dimension, NULL) ;
	dataset = H5Dopen(group_id, name) ;
	H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
	H5Dclose(dataset) ;
	H5Sclose(memspace) ;
      }
      H5Sclose(dataspace) ;  
    }
    delete [] tmp_int ;
  }
  
  storeRepP collect_reorder_store(storeRepP &sp, dMap& remap, fact_db &facts) {
    entitySet dom = sp->domain() ;
    dom =  collect_entitySet(dom) ;
    dom =  all_collect_entitySet(dom) ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    std::vector<entitySet> chop_ptn = d->chop_ptn ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    entitySet remap_dom ;
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
      else
        if(d_remap[ri].size() > 1)
          cerr << "d_remap has multiple entries!" << endl ;
    } ENDFORALL ;
    Map tmp_remap ;
    entitySet tmp_remap_dom =  MapRepP(l2g.Rep())->preimage(dom&init_ptn[ MPI_rank]).first ;
    tmp_remap.allocate(tmp_remap_dom) ;
    FORALL(tmp_remap_dom, ri) {
      tmp_remap[ri] = l2g[ri] ;
    } ENDFORALL ;
    entitySet owned_entities = *my_entities & sp->domain() ;

    entitySet tmp_remap_image = MapRepP(tmp_remap.Rep())->image(tmp_remap.domain());
    entitySet reverse_dom = reverse.domain() ;
    entitySet tmp_remap_preimage = MapRepP(tmp_remap.Rep())->preimage(reverse_dom).first ;

    if((tmp_remap_dom&tmp_remap_preimage) != tmp_remap_dom) {
      debugout << "tmp_remap.image() = " << tmp_remap_image << endl ;
      debugout << "reverse.domain() = " << reverse.domain() << endl ;
      debugout << "tmp_remap.preimage(reverse.domain()) = "
               << tmp_remap_preimage << endl ;
      cerr << "something fishy in collect_reorder_store" << endl ;
      cerr << "missing entities" << tmp_remap_dom-tmp_remap_preimage << endl ;
    }
    MapRepP(tmp_remap.Rep())->compose(reverse, tmp_remap_dom&tmp_remap_preimage) ;
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
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp[i] = sp->new_store(entitySet(recv_dom[i])) ;
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
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_unpack = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      tmp_sp[i]->unpack(recv_store, loc_unpack, size_recv, recv_dom[i]) ; 
    }
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
  } 
  
  void distribute_reorder_store(storeRepP &new_sp, storeRepP sp_init, dMap& remap, fact_db &facts) {
    entitySet our_dom = new_sp->domain() ;
    entitySet read_set = Loci::collect_entitySet(our_dom) ;
    read_set = all_collect_entitySet(read_set) ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    std::vector<entitySet> chop_ptn = d->chop_ptn ;
    Loci::constraint my_entities ; 
    Map l2g ; 
    l2g = facts.get_variable("l2g") ;
    entitySet remap_dom = remap.domain() & read_set ;
    entitySet local_own = d->g2l.domain() & remap.image(remap_dom) ;
    entitySet t_dom = remap.preimage(local_own).first ;
    dMap tmp_remap ; 
    FORALL(t_dom, ti) { 
      tmp_remap[ti] = remap[ti] ;
    } ENDFORALL ; 
    Loci::MapRepP(tmp_remap.Rep())->compose(d->g2l, t_dom) ;
    my_entities = d->my_entities ; 
    entitySet local_req = *my_entities & our_dom ;
    entitySet global_req = Loci::MapRepP(l2g.Rep())->image(local_req) ;
    global_req = Loci::MapRepP(remap.Rep())->preimage(global_req).first ;
    std::vector<int> sizes(Loci::MPI_processes) ;
    for(int i = 0; i < Loci::MPI_processes; ++i)
      sizes[i] = (read_set & chop_ptn[i]).size() ;
    std::sort(sizes.begin(), sizes.end()) ;
    int max_size = sizes[Loci::MPI_processes-1] ;
    entitySet max_set = interval(0, max_size-1) ;
    unsigned char *tmp_buf, *recv_buf ;
    int total_size = 0 ;
    total_size = sp_init->pack_size(max_set) ;
    tmp_buf = new unsigned char[total_size] ;
    recv_buf = new unsigned char[total_size] ;
    int loc = 0 ;
    entitySet tmp_dom = sp_init->domain() ;
    sp_init->pack(tmp_buf, loc, total_size, tmp_dom) ; 
    MPI_Status status, stat_1 ;
    if(Loci::MPI_rank != 0) {
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	MPI_Recv(recv_buf, total_size, MPI_PACKED, Loci::MPI_rank-1, 12,
		 MPI_COMM_WORLD, &status) ;  
	if(Loci::MPI_rank < Loci::MPI_processes-1) { 
	  MPI_Send(recv_buf, total_size, MPI_PACKED, Loci::MPI_rank+1, 12, MPI_COMM_WORLD) ;
	}
	if(p == 0)
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, 0, 12, MPI_COMM_WORLD) ;
	entitySet recv_set = read_set & chop_ptn[p] ;
	Loci::sequence tmp_seq = Loci::sequence(recv_set) ;
	int loc_unpack = 0 ;
	storeRepP sp = sp_init->new_store(recv_set) ;
	sp->unpack(recv_buf, loc_unpack, total_size, tmp_seq) ;
	entitySet local_set = recv_set & global_req ;
	new_sp->scatter(tmp_remap, sp, local_set) ; 
      }
    } else {
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	if(p > 0) 
	  MPI_Recv(recv_buf, total_size, MPI_PACKED, p, 12,
		   MPI_COMM_WORLD, &stat_1) ;
	if(p == 0)
	  MPI_Send(tmp_buf, total_size, MPI_PACKED, 1, 12, MPI_COMM_WORLD) ;
	else
	  MPI_Send(recv_buf, total_size, MPI_PACKED, 1, 12, MPI_COMM_WORLD) ;
	dMap tmp_map , new_map;
	entitySet local_set = chop_ptn[p] & read_set ;  
	int loc_unpack = 0 ;
	Loci::sequence tmp_seq = Loci::sequence(local_set) ;
	storeRepP sp = sp_init->new_store(local_set) ;
	if(p != 0)
	  sp->unpack(recv_buf, loc_unpack, total_size, tmp_seq) ;
	local_set &= global_req ;
	if(p == 0)
	  new_sp->scatter(tmp_remap, sp_init, local_set) ; 
	else 
	  new_sp->scatter(tmp_remap, sp, local_set) ; 
      }
    }
    delete [] recv_buf ;
    delete [] tmp_buf ;
  }
  
  
  storeRepP distribute_store(storeRepP &sp, fact_db &facts) {
    if(!facts.isDistributed()) {
      return sp ;
    }
    storeRepP nsp ;
    if(facts.isDistributed()) {  
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
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
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	recv_buffer = new int*[MPI_processes-1] ;
	for(int i = 0; i < MPI_processes-1; ++i) {
	  recv_buffer[i] = new int[recv_size[i]] ;
	}
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,2, MPI_COMM_WORLD, &recv_request[k] );  

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
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
#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, store_request, store_status) ;
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
  
  Map distribute_global_map(Map &m, const std::vector<entitySet> &init_ptn) {
    if(Loci::MPI_processes == 1)
      return m ;
    
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    entitySet image ;
    if(Loci::MPI_rank == 0)
      image = MapRepP(m.Rep())->image(dom) ;
    image = all_collect_entitySet(image) ;

    entitySet my ;
    if(Loci::MPI_rank == 0) {
      std::vector<entitySet> ent(MPI_processes) ;
      MPI_Status *store_status ;
      MPI_Request *store_request ;
      int sz = 0 ;
      int *s_size = new int[MPI_processes] ;
      
      for(int i = 0; i < MPI_processes; ++i) 
	ent[i] = init_ptn[i] & image ;
      for(int i = 0; i < MPI_processes-1; ++i) { 
	s_size[i] = 2 * ent[i+1].size() ;
	sz += s_size[i] ;
      }
      std::vector<entitySet> final_ent(MPI_processes) ;
      for(int i = 0; i < MPI_processes; ++i) {
	dom = MapRepP(m.Rep())->preimage(ent[i]).first ;
	final_ent[i] = dom ;
      }
      /*
      for(int i = 0; i < MPI_processes; ++i) {   
	for(int j = i+1 ; j < MPI_processes; ++j) {
	  entitySet  tmp = final_ent[i] & final_ent[j] ;
	  if(tmp != EMPTY)
	    cerr << " ERROR: Something is screwed up ???  " << endl ;
	  final_ent[j] -= tmp ;
	}
      }
      */
      int **send_ptr = new int*[MPI_processes] ;
      store_request = new MPI_Request[MPI_processes] ;
      store_status = new MPI_Status[MPI_processes] ;
      int tmp = 0 ;
      send_ptr[0] = new int[sz*10] ;
      for(int i = 1; i < MPI_processes-1; i++)
	send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
      for(int i = 0; i < MPI_processes-1; i++) {
	tmp = 0 ;
	for(entitySet::const_iterator ei = final_ent[i+1].begin(); ei != final_ent[i+1].end(); ++ei) {
	  send_ptr[i][tmp++] = *ei ;
	  send_ptr[i][tmp++] = m[*ei] ;
	  
	}
	MPI_Isend(send_ptr[i], s_size[i], MPI_INT, i+1, 3,
		  MPI_COMM_WORLD, &store_request[i]) ;
      }
#ifdef DEBUG
      int err =
#endif
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
      FATAL(err != MPI_SUCCESS) ;
      if(Loci::MPI_processes > 1)
	my = final_ent[0] ;
      else
	my = dom ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei)
	nm[*ei] = m[*ei] ;
      delete [] store_request ;
      delete [] store_status ;
      delete [] send_ptr ;
    }
    else {
      MPI_Status stat ;
      int sz = 2*(init_ptn[Loci::MPI_rank] & image).size();
      int *recv_ptr = new int[sz*10] ;
      MPI_Recv(recv_ptr, sz, MPI_INT, 0, 3, MPI_COMM_WORLD, &stat) ;
      int tmp = 1 ;
      dom = EMPTY ;
      for(int i = 0; i < sz; i += 2)
	dom += recv_ptr[i] ;
      nm.allocate(dom) ;
      for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) {
	nm[*ei] = recv_ptr[tmp] ;
	tmp += 2 ;
      }
      delete [] recv_ptr ;
    }
    return nm ;
  }
  
  Map distribute_gmap(Map &m, const std::vector<entitySet> &init_ptn) {
    if(Loci::MPI_processes == 1)
      return m ;
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    dom = Loci::all_collect_entitySet(dom) ;

    entitySet my ;
    if(Loci::MPI_rank == 0) {
      std::vector<entitySet> ent(MPI_processes-1) ;
      MPI_Status *store_status ;
      MPI_Request *store_request ;
      int sz = 0 ;
      int *s_size = new int[MPI_processes-1] ;
      for(int i = 0; i < MPI_processes-1; ++i) { 
	ent[i] = init_ptn[i+1] & dom ;
	s_size[i] = ent[i].size() ;
	sz += s_size[i] ;
      }
      int **send_ptr = new int*[MPI_processes-1] ;
      store_request = new MPI_Request[MPI_processes-1] ;
      store_status = new MPI_Status[MPI_processes-1] ;
      int tmp = 0 ;
      send_ptr[0] = new int[sz] ;
      for(int i = 1; i < MPI_processes-1; i++)
	send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
      for(int i = 0; i < MPI_processes-1; i++) {
	tmp = 0 ;
	for(entitySet::const_iterator ei = ent[i].begin(); ei != ent[i].end(); ++ei) 
	  send_ptr[i][tmp++] = m[*ei] ;
	MPI_Isend(send_ptr[i], s_size[i], MPI_INT, i+1, 3,
		  MPI_COMM_WORLD, &store_request[i]) ;
      }
#ifdef DEBUG
      int err =
#endif
	MPI_Waitall(MPI_processes-1, store_request, store_status) ;
      FATAL(err != MPI_SUCCESS) ;
      my = init_ptn[0] & dom ; ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei)
	nm[*ei] = m[*ei] ;
      delete [] store_request ;
      delete [] store_status ;
      delete [] send_ptr ;
    }
    else {
      MPI_Status stat ;
      my = init_ptn[Loci::MPI_rank] & dom ;
      int sz = my.size();
      int *recv_ptr = new int[sz] ;
      MPI_Recv(recv_ptr, sz, MPI_INT, 0, 3, MPI_COMM_WORLD, &stat) ;
      int tmp = 0 ;
      nm.allocate(my) ;
      for(entitySet::const_iterator ei = my.begin(); ei != my.end(); ++ei) 
	nm[*ei] = recv_ptr[tmp++] ;
      delete [] recv_ptr ;
    }
    return nm ;
   }
  
  Map distribute_whole_map(Map &m) {
    if(Loci::MPI_processes == 1)
      return m ;
    Map nm ;
    entitySet dom ;
    if(Loci::MPI_rank == 0)
      dom = m.domain() ;
    dom = Loci::all_collect_entitySet(dom) ;

    nm.allocate(dom) ;
    int sz = 0 ;
    sz = dom.size() ;
    int *recv_ptr = new int[sz] ;
    if(Loci::MPI_rank == 0) {
      int tmp = 0 ;
      for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) 
	recv_ptr[tmp++] = m[*ei] ;
    }
    int tmp = 0 ;
    MPI_Bcast(recv_ptr, sz, MPI_INT, 0, MPI_COMM_WORLD) ;
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei) 
      nm[*ei] = recv_ptr[tmp++] ;
    delete [] recv_ptr ;
    return nm ;
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
    std::vector<int> vset( MPI_processes) ;
    if(MPI_processes > 1) {
      int *recv_count = new int[ MPI_processes] ;
      int *send_count = new int[ MPI_processes] ;
      for(int i = 0; i <  MPI_processes; ++i) 
	send_count[i] = size ;
      
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ; 
      for(int i = 0; i <  MPI_processes; ++i)
	vset[i] = recv_count[i] ;
      
      delete [] send_count ;
      delete [] recv_count ;
    }
    else
      vset[0] = size ;
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


  entitySet broadcastEntitySet(entitySet s, int root) {
    int num_ivals = s.num_intervals() ;

    MPI_Bcast(&num_ivals,1,MPI_INT,root,MPI_COMM_WORLD) ;
    int *buf = new int[num_ivals*2] ;
    if(MPI_rank == root) {
      for(int i=0;i<num_ivals;++i) {
        buf[2*i] = s[i].first ;
        buf[2*i+1] = s[i].second ;
      }
    }
    MPI_Bcast(&buf,num_ivals*2,MPI_INT,root,MPI_COMM_WORLD) ;
    entitySet ret_val ;
    for(int i=0;i<num_ivals;++i)
      ret_val += interval(buf[2*i],buf[2*i+1]) ;
    delete[] buf ;
    return ret_val ;
  }


  // Pass in a set of entitySets one for each processor to broadcast
  // to other processors.  Return with a vector of entitySets as sent to
  // you by each other processor.
  vector<entitySet> Alltoall_entitySet(vector<entitySet> v) {
    WARN(int(v.size()) != MPI_processes) ;

    if(MPI_processes == 1)
      return v ;

    const int p = v.size() ;
    vector<int> ivals(p) ;
    for(int i=0;i<p;++i)
      ivals[i] = v[i].num_intervals() ;
    vector<int> ivalr(p) ;
    MPI_Alltoall(&(ivals[0]),1,MPI_INT,&(ivalr[0]),1,MPI_INT,MPI_COMM_WORLD) ;
    vector<int> sdispls(p),rdispls(p) ;
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+ivals[i-1]*2 ;
      rdispls[i] = rdispls[i-1]+ivalr[i-1]*2 ;
    }
    vector<int> sbuf(sdispls[p-1]+ivals[p-1]*2) ;
    vector<int> rbuf(rdispls[p-1]+ivalr[p-1]*2) ;
    vector<int> scounts(p),rcounts(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = 2*ivals[i] ;
      rcounts[i] = 2*ivalr[i] ;
      for(int j=0;j<ivals[i];++j) {
        sbuf[sdispls[i]+j*2] = v[i][j].first ;
        sbuf[sdispls[i]+j*2+1] = v[i][j].second ;
      }
    }
    MPI_Alltoallv(&(sbuf[0]),&(scounts[0]),&(sdispls[0]),MPI_INT,
                  &(rbuf[0]),&(rcounts[0]),&(rdispls[0]),MPI_INT,
                  MPI_COMM_WORLD) ;

    vector<entitySet> retv(p) ;
    for(int i=0;i<p;++i) {
      for(int j=0;j<ivalr[i];++j) {
        retv[i] += interval(rbuf[rdispls[i]+j*2],rbuf[rdispls[i]+j*2+1]) ;
      }
    }
    return retv ;
  }
  
  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) {
    if(MPI_processes == 1)
      return m ;

    const int p = MPI_processes ;
    entitySet dom = m.domain() ;
    std::vector<entitySet> send_slices(p) ;
    for(int i=0;i<p;++i) {
      send_slices[i] = dom & init_ptn[i] ;
      dom -= init_ptn[i] ;
    }
    WARN(dom != EMPTY) ;

    vector<entitySet> recv_slices = Alltoall_entitySet(send_slices) ;

    vector<int> scounts(p),rcounts(p), sdispls(p),rdispls(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = send_slices[i].size() ;
      rcounts[i] = recv_slices[i].size() ;
    }
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    vector<int> sbuf(sdispls[p-1]+scounts[p-1]) ;
    vector<int> rbuf(rdispls[p-1]+rcounts[p-1]) ;
    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = send_slices[i].begin();
          ei != send_slices[i].end();++ei,++k) {
        sbuf[sdispls[i]+k] = m[*ei] ;
      }
    }

    MPI_Alltoallv(&(sbuf[0]),&(scounts[0]),&(sdispls[0]),MPI_INT,
                  &(rbuf[0]),&(rcounts[0]),&(rdispls[0]),MPI_INT,
                  MPI_COMM_WORLD) ;

    dMap ret_map ;

    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = recv_slices[i].begin();
          ei != recv_slices[i].end();++ei,++k) {
        ret_map[*ei] = rbuf[rdispls[i]+k] ;
      }
    }

    return ret_map ;
    
  }
} 

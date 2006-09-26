#include <vector>
using std::vector;

#include <set>
using std::set;


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <algorithm>
using std::sort;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>
#include "loci_globs.h"

namespace Loci {
  void get_clone(fact_db &facts, const rule_db &rdb) {
    fact_db::distribute_infoP df = facts.get_distribute_info()  ;
    std::vector<entitySet> &ptn = facts.get_init_ptn() ;
    entitySet bdom = ptn[MPI_rank] & interval(UNIVERSE_MIN, -1) ;
    entitySet global_bdom = all_collect_entitySet(bdom) ;
    int p = 0; 
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = ptn[i] ; 
      ptn[i] = tmp & interval(0, UNIVERSE_MAX) ;
    }

    FORALL(global_bdom, i) {
      int tmp = p % MPI_processes ;
      ptn[tmp] += i ;
      if(duplicate_work)
	if(tmp == MPI_rank)
	  //Since we are manually adding entities to init_ptn, we also need to 
	  //add those entities in global_comp_entities
	  facts.global_comp_entities += i;
      p++ ;
    } ENDFORALL ;
    variableSet tmp_vars = facts.get_typed_variables();
    for(variableSet::const_iterator vi = tmp_vars.begin(); vi != tmp_vars.end(); ++vi) {
      storeRepP tmp_sp = facts.get_variable(*vi) ;
      if(tmp_sp->RepType() == CONSTRAINT) {
        entitySet tmp_dom = tmp_sp->domain() ;
        if(tmp_dom != ~EMPTY) {
          entitySet global_tmp_dom = all_collect_entitySet(tmp_dom) ;
          constraint tmp ;
          *tmp = global_tmp_dom ;
          //facts.update_fact(variable(*vi), tmp) ;
          facts.replace_fact(*vi,tmp) ;
        }
      }
      if(tmp_sp->RepType() == MAP) {
        storeRepP map_sp = MapRepP(tmp_sp->getRep())->thaw() ;
        facts.replace_fact(*vi, map_sp) ;
      }
    }

    entitySet tmp_set;
    if(duplicate_work) {
      std::set<std::vector<variableSet> > context_maps ;

      if(extended_duplication) 
	get_mappings(rdb, facts, context_maps);
      else 
	get_mappings(rdb, facts, context_maps, 1);

      context_maps = classify_moderate_maps(facts, context_maps);

      //Find out entities that can produce target variable entities owned by a processor
      facts.global_comp_entities += context_for_map_output(ptn[MPI_rank], facts, context_maps);
      
      //Add entities so that maximum depth of duplication can be achieved
      if(multilevel_duplication) {
	entitySet mySet = ptn[MPI_rank];
	bool continue_adding = true;
	int num_levels = 1;
	do{
	  mySet += dist_special_expand_map(facts.global_comp_entities,
						 facts, context_maps) ;
	  entitySet added_entities = context_for_map_output(mySet,  facts, context_maps);
	  added_entities -= facts.global_comp_entities;
	  if(all_collect_entitySet(added_entities) == EMPTY)
	    continue_adding = false;
	  else {
	    facts.global_comp_entities += added_entities;
	    num_levels++;
	  }
	}while(continue_adding);
	facts.global_comp_entities += mySet;
	debugout << "Number of Duplication Levels: " << num_levels << endl; 
      }
      tmp_set = facts.global_comp_entities;
    }
    else
      tmp_set = ptn[MPI_rank] ;
    std::set<std::vector<variableSet> > dist_maps ;
    get_mappings(rdb,facts,dist_maps) ;
    entitySet tmp_copy, image ;
    image = dist_expand_map(tmp_set, facts, dist_maps) ;
    tmp_copy =  image - ptn[MPI_rank] ; 
    std::vector<entitySet> copy(MPI_processes), send_clone(MPI_processes) ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes];
    int *recv_displacement = new int[MPI_processes];
    int size_send = 0 ;
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
    entitySet::const_iterator ei ;
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
    categories(facts,iv) ;
    for(size_t i=0;i < iv.size(); ++i) 
      //      iv[i] = all_collect_entitySet(iv[i]) ;
      iv[i] = dist_expand_entitySet(iv[i],tmp_copy,ptn) ;
    

    entitySet::const_iterator ti ;
    vector<entitySet> proc_entities ;
    entitySet e ;

    debugout << " initial_categories =  " << iv.size() << endl ;
  
    for(size_t i = 0; i < iv.size(); ++i) {
      //    debugout << " iv[ " << i << " ] = " << iv[i] << endl ;
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

    df->l2g = l2g.Rep() ;
    df->dl2g = MapRepP(l2g.Rep())->thaw() ;
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
    debugout << "  Time taken for creating initial info =  " << end_time - start << endl ;
    start = MPI_Wtime() ;
    reorder_facts(facts, df->g2l) ;
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for reordering =  " << end_time - start << endl ; 
    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    g = EMPTY ;
#ifdef DEBUG
    entitySet g2ldom = df->g2l.domain() ;
    if(ptn[myid]-g2ldom != EMPTY) {
      cerr << "problem with g2l " << ptn[myid]-g2ldom << endl ;
    }
#endif
    for(ei = ptn[myid].begin(); ei != ptn[myid].end(); ++ei)
      g += df->g2l[*ei] ;
    my_entities = g ;
    df->myid = myid ;
    df->my_entities = g ;

    //Add comp_entities
    if(duplicate_work) {
      g = EMPTY;
      for(ei = facts.global_comp_entities.begin();
	ei != facts.global_comp_entities.end(); ++ei)
	g += df->g2l[*ei] ;
      df->comp_entities = g;
    }
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
    // this needs to be an intensional fact
    facts.create_intensional_fact("l2g", l2g) ;
    facts.put_l2g(l2g) ;
    facts.create_intensional_fact("my_entities",my_entities);
  }
  
  /*This routine loops over all the rules in the database and extracts
    all the variables associated with the mappings in the head, body and
    the constraints of the rules. */
  //rule_part = 0: Output and Input of the rule mapping
  //          = 1: Only Outupt
  //          = 2: Only Input
  void get_mappings(const rule_db &rdb, fact_db &facts,
                    set<vector<variableSet> > &maps_ret, unsigned int rule_part) {
    ruleSet rules = rdb.all_rules() ;
    set<vector<variableSet> > maps ;
    
    for(ruleSet::const_iterator ri = rules.begin(); ri != rules.end(); ++ri) {
      set<vmap_info>::const_iterator vmsi ;
      if(rule_part != 2) {
	for(vmsi = ri->get_info().desc.targets.begin();
	    vmsi != ri->get_info().desc.targets.end();
	    ++vmsi) {
	  if(vmsi->mapping.size() != 0) {
	    vector<variableSet> vvs ;
	    for(size_t i = 0; i < vmsi->mapping.size(); ++i) {
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
      }
      if(rule_part != 1) {
	for(vmsi = ri->get_info().desc.sources.begin();
	    vmsi != ri->get_info().desc.sources.end();
	    ++vmsi) {
	  if(vmsi->mapping.size() != 0) {
	    vector<variableSet> vvs ;
	    for(size_t i = 0; i < vmsi->mapping.size(); ++i) {
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
	  vector<variableSet> vvs ;
	  if(vmsi->mapping.size() != 0) {
	    for(size_t i = 0; i < vmsi->mapping.size(); i++) {
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
      }
    }

    set<vector<variableSet> >::const_iterator mi ;

    variableSet vars = facts.get_typed_variables() ;
    maps_ret.clear() ;
    for(mi=maps.begin();mi!=maps.end();++mi) {
      const vector<variableSet> &vss = *mi ;
      // If the maps aren't in the fact database then exclude it
      variableSet notvars ;
      for(size_t i=0;i<vss.size();++i)
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
            for(size_t i=2;i<vss.size();++i)
              vs[i] = vss[i] ;
            maps_ret.insert(vs) ;
          }
      }
            
    }

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
      for(size_t i = 0; i < mv.size(); ++i) {
	variableSet v = mv[i] ;
	v &= vars ; 
	entitySet image ;
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    entitySet tmp_dom = p->domain() ;
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	    entitySet tmp_out = (glob_dom & locdom) - tmp_dom ; 
	    storeRepP sp = mp->expand(tmp_out, ptn) ;
	    if(sp->domain() != tmp_dom) {
	      //facts.update_fact(variable(*vi), sp) ;
              facts.replace_fact(*vi,sp) ;
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
        for(size_t i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->copy[i-1].size ;
          recv_size[i] = d->copy[i].size ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size] ;
        for(size_t i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->xmit[i-1].size ;
      }
      
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;
      /*The recv_size is the maximum possible, so that even if we
	receive a shorter message there won't be any problem */
      for(size_t i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(size_t i=0;i<d->xmit.size();++i) {
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
      for(size_t i = 0; i < d->copy.size(); ++i) {
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
        for(size_t i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->copy[i-1].size*evsz+evsz) ;
          recv_size[i] = d->copy[i].size*evsz+evsz ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        for(size_t i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
      }
        
      
      Map l2g ;
      l2g = facts.get_variable("l2g") ;
      
      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;

      for(size_t i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(size_t i=0;i<d->xmit.size();++i) {
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

      for(size_t i = 0; i < d->copy.size(); ++i) {
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

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->xmit[i-1].size ;
          recv_size[i] = d->xmit[i].size ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      /*By intersecting the given entitySet with the clone region
	entities we can find out which entities are to be sent */
      for(size_t i=0;i<d->copy.size();++i) {
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
      
      for(size_t i=0;i<d->xmit.size();++i) {
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

  //Similar to send_entitySet defined above.
  //Instead of returning set of entities, it returns
  //vector of entitySet describing exactly which entities 
  //are received from which processor.
  vector<entitySet> send_entitySetv(const entitySet& e, fact_db &facts) {
    vector<entitySet> re(MPI_processes);
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

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+d->xmit[i-1].size ;
          recv_size[i] = d->xmit[i].size ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      /*By intersecting the given entitySet with the clone region
	entities we can find out which entities are to be sent */
      for(size_t i=0;i<d->copy.size();++i) {
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
      
      for(size_t i=0;i<d->xmit.size();++i) {
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        for(int j=0;j<recieved;++j)
          re[d->xmit[i].proc] += d->g2l[recv_buffer[i][j]] ;
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

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(d->xmit[i-1].size*evsz+evsz) ;
          recv_size[i] = d->xmit[i].size*evsz+evsz ;
        }
      }
      
      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[d->copy_total_size*evsz+evsz*d->copy.size()] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+d->copy[i-1].size*evsz+evsz ;
      }
      Map l2g ;
      l2g = facts.get_variable("l2g") ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;  

      for(size_t i=0;i<d->copy.size();++i) {
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
      for(size_t i=0;i<d->xmit.size();++i) {
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

  std::set<std::vector<variableSet> > classify_moderate_maps(fact_db &facts, const std::set<std::vector<variableSet> > &maps) {
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    std::set<std::vector<variableSet> > return_maps;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      entitySet domain;
      entitySet image;
      entitySet tmp;
      const vector<variableSet> &mv = *smi ;
      for(unsigned int i = 0; i < mv.size(); i++) {
	variableSet v = mv[i] ;
	v &= vars ;
	if(i == 0) {
	  for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	    storeRepP p = facts.get_variable(*vi) ;
	    if(p->RepType() ==  MAP) {
	      MapRepP mp =  MapRepP(p->getRep()) ;
	      domain += mp->domain();
	    }
	  }
	  domain = all_collect_entitySet(domain);
	  tmp = domain;
	}
	  
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    image += mp->image(tmp);
	  }
	}

	tmp = all_collect_entitySet(image);
	image = EMPTY;
      }

      image = tmp;
      debugout << "map ratio: " << 1.0*domain.size()/image.size() << endl;
      for(unsigned int i = 0; i < mv.size(); i++) {
	debugout << "(" << mv[i] << ")  ->  "; 
      }
      debugout << endl;
      if(domain.size()*1.0 > image.size() *  10.0) {
	facts.intensive_output_maps.insert(*smi);
      }
      else {
	return_maps.insert(*smi);
      }
    }
    return return_maps;
  } 

  
  //Finds the context for maps whose final image will be in provided domain.
  entitySet context_for_map_output(entitySet domain, fact_db &facts,
	      const std::set<std::vector<variableSet> > &maps) {
    std::vector<entitySet> ptn = facts.get_init_ptn() ;
    entitySet context;
    variableSet vars = facts.get_typed_variables() ;
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = maps.begin(); smi != maps.end(); ++smi) {
      std::vector<entitySet>  preimage_vec = all_collect_vectors(domain);
      const vector<variableSet> &mv = *smi ;
      for(int i = mv.size() -1; i >= 0; --i) {
	variableSet v = mv[i] ;
	v &= vars ;
	std::vector<entitySet>  tmp_preimage_vec(MPI_processes);
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    for(int j = 0; j < MPI_processes; j++) {
	      tmp_preimage_vec[j] += mp->preimage(preimage_vec[j]).second;
	    }
	  }
	}
	for(int j = 0; j < MPI_processes; j++) {
	  preimage_vec[j] = all_collect_entitySet(tmp_preimage_vec[j]);
	}

	if(i == 0) {
	  context += preimage_vec[MPI_rank];
	}
      }
    }
    return context ;
  }

  //It works just like dist_expand_map except the return entities only contain the image
  //of the last map in the chain.  
  entitySet dist_special_expand_map(entitySet domain, fact_db &facts,
				    const std::set<std::vector<variableSet> > &maps) {
    entitySet special_return;
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
      for(size_t i = 0; i < mv.size(); ++i) {
	variableSet v = mv[i] ;
	v &= vars ; 
	entitySet image ;
	for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	  storeRepP p = facts.get_variable(*vi) ;
	  if(p->RepType() ==  MAP) {
	    entitySet tmp_dom = p->domain() ;
	    MapRepP mp =  MapRepP(p->getRep()) ;
	    entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	    entitySet tmp_out = (glob_dom & locdom) - tmp_dom ; 
	    storeRepP sp = mp->expand(tmp_out, ptn) ;
	    if(sp->domain() != tmp_dom) {
	      //facts.update_fact(variable(*vi), sp) ;
              facts.replace_fact(*vi,sp) ;
	    }
	    image +=  MapRepP(sp)->image((sp->domain()) & locdom) ;

	  }
	}
	dom += image ;
	locdom = image ;
	if(i == mv.size() -1)
	  special_return += image;
      }
    }
    return special_return ;
  }

  // function that restores the fact_db back to its global numbering
  void restore_global_facts(fact_db& facts) {
    if(MPI_rank == 0)
      cerr << "restore_global_facts" << endl ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    // not yet distributed, we don't need to do anything
    if(df == 0)
      return ;
    if(!df->isDistributed)
      return ;

    facts.remove_variable(variable("l2g")) ;
    facts.remove_variable(variable("my_entities")) ;
    df->xmit.clear() ;
    df->copy.clear() ;
    dMap l2g ;
    entitySet dom = df->l2g.domain() ;
    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei)
      l2g[*ei] = df->l2g[*ei] ;

    if(MPI_rank==0)
      cerr << "before reorder" << endl ;
    
    reorder_facts(facts, l2g) ;

    if(MPI_rank==0)
      cerr << "after reorder" << endl ;

    df->isDistributed = 0 ;
  }

}


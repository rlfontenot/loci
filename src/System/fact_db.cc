#include <fact_db.h>
#include <constraint.h>
#include <Tools/stream.h>
#include <Tools/debugger.h>
#include <DStore.h>

#include <typeinfo>

extern "C" {
#include <hdf5.h>
}

using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;
using std::make_pair ;

#include <Tools/parse.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;
  extern fact_db *exec_current_fact_db;

  fact_db::fact_db() {
    constraint EMPTY ;
    create_fact("EMPTY",EMPTY) ;
    constraint UNIVERSE ;
    UNIVERSE = ~EMPTY ;
    create_fact("UNIVERSE",UNIVERSE) ;
    distributed_info = 0 ;
    maximum_allocated = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      init_ptn.push_back(EMPTY) ;
    }
    exec_current_fact_db = this;
  }

  fact_db::~fact_db() {}

  void fact_db::synonym_variable(variable v, variable synonym) {

    // Find all variables that should be synonymous with v
    variableSet synonym_set ;
    std::map<variable,variable>::const_iterator mi ;
    while((mi=synonyms.find(v)) != synonyms.end()) {
      synonym_set += v ;
      v = mi->second ;
    }
    variable s = synonym;
    while((mi=synonyms.find(s)) != synonyms.end()) {
      synonym_set += s ;
      s = mi->second ;
    }
    synonym_set += s ;

    // If the two are already synonymous, we are done
    if(s == v)
      return ;

    // Sanity check, nake sure v exists
    std::map<variable,fact_info>::iterator vmi, vmj ;
    if((vmi = fmap.find(v)) == fmap.end()) {
      cerr << "variable type not known for target of synonym" << endl ;
      abort() ;
    }
    // If the synonym already points to a different variable instance,
    // remove it
    if((vmj = fmap.find(s)) != fmap.end()) {
      fact_info &finfo = vmj->second ;
      if((finfo.data_rep->domain() != EMPTY &&
          finfo.data_rep->RepType() != PARAMETER)) {
        cerr << "unable to define synonym variable " << synonym
             << " when varaiable already created in db. "  << endl ;
        cerr << "variable v = " << v << endl ;
        abort() ;
      }
      remove_variable(synonym) ;
    }

    // Add new synonyms so that they point to v
    for(variableSet::const_iterator vi = synonym_set.begin();
        vi!=synonym_set.end();
        ++vi) {
      synonyms[*vi] = v ;
    }
  }
  
  
  
  void fact_db::update_fact(variable v, storeRepP st) {
    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }
    warn(synonyms.find(v) != synonyms.end()) ;
    std::map<variable, fact_info>::iterator mi = fmap.find(v) ;
    
    if(mi != fmap.end()) {
      mi->second.data_rep->setRep(st->getRep()) ;
    } else
      cerr << "warning: update_fact: fact does not exist for variable " << v
           << endl ;
  }
  
  void fact_db::create_fact(variable v, storeRepP st) {

    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }

    if(synonyms.find(v) != synonyms.end()) {
      v = remove_synonym(v) ;
      std::map<variable, fact_info>::iterator mi = fmap.find(v) ;
      if(mi==fmap.end()) {
        fmap[v].data_rep = new store_ref ;
        fmap[v].data_rep->setRep(st->getRep()) ;
      } else {
        if(typeid(st->getRep()) != typeid(mi->second.data_rep->getRep())) {
          cerr << "set_variable_type() method of fact_db changing type for variable " << v << endl ;
        }
        mi->second.data_rep->setRep(st->getRep()) ;
      }
      return ;
    }
    
    std::map<variable, fact_info>::iterator mi = fmap.find(v) ;
    if(mi != fmap.end()) {
      cerr << "WARNING: fact_db::set_variable_type retyping variable "
	   << v << endl ;
      mi->second.data_rep->setRep(st->getRep()) ;
    } else {
      fmap[v].data_rep = new store_ref ;
      fmap[v].data_rep->setRep(st->getRep()) ;
    }
    
  } 

  void fact_db::remove_variable(variable v) {
    std::map<variable, variable>::iterator si ;
    std::map<variable, fact_info>::iterator mi ;
    if((si=synonyms.find(v)) != synonyms.end()) {
      variable real_var = remove_synonym(v) ;
      synonyms.erase(si) ;
      remove_variable(real_var) ;
    } else if((mi=fmap.find(v)) != fmap.end()) {
      // First remove any synonyms to this variable.
      variableSet syn_vars ;
      vector<map<variable,variable>::iterator > lrm ;
      for(si=synonyms.begin();si!=synonyms.end();++si)
        if(si->second == v)
          syn_vars += si->first ;
      for(variableSet::const_iterator vi=syn_vars.begin();
          vi!=syn_vars.end();++vi) {
        if((si=synonyms.find(*vi)) != synonyms.end()) {
          synonyms.erase(si) ;
        }
      }

      // Now erse the variable
      fmap.erase(mi) ;
    }
  }
  
  
  variableSet fact_db::get_typed_variables() const {
    std::map<variable, fact_info>::const_iterator mi ;
    std::map<variable, variable>::const_iterator si ;
    variableSet all_vars ;
    for(mi=fmap.begin();mi!=fmap.end();++mi)
      all_vars += mi->first ;
    
    for(si=synonyms.begin();si!=synonyms.end();++si)
      all_vars += si->first ;
    return all_vars ;
  }
    

  std::pair<interval, interval> fact_db::get_distributed_alloc(int size) {
    dist_from_start = 1 ;
    if(MPI_processes > 1) {
      int* send_buf = new int[MPI_processes] ;
      int* size_send = new int[MPI_processes] ;
      int* size_recv = new int[MPI_processes] ;
      int* recv_buf = new int[MPI_processes] ;
      for(int i = 0; i < MPI_processes; ++i) {
	send_buf[i] = maximum_allocated ;
	size_send[i] = size ;
      } 
      MPI_Alltoall(send_buf, 1, MPI_INT, recv_buf, 1, MPI_INT, MPI_COMM_WORLD) ;
      MPI_Alltoall(size_send, 1, MPI_INT, size_recv, 1, MPI_INT, MPI_COMM_WORLD) ;
      std::sort(recv_buf, recv_buf+MPI_processes) ;
      maximum_allocated = recv_buf[MPI_processes-1] ;
      int local_max = maximum_allocated ;
      int global_max = 0 ;
      for(int i = 0; i < MPI_rank; ++i)
	local_max += size_recv[i] ;
      for(int i = 0; i < MPI_processes; ++i) 
	global_max += size_recv[i] ;
      
      for(int i = 0 ; i < MPI_processes; ++i) {
	int local = maximum_allocated ;
	for(int j = 0; j < i; ++j)
	  local += size_recv[j] ;
	init_ptn[i] += interval(local, local+size_recv[i]-1) ;
      }
      
      interval local_ivl = interval(local_max, local_max + size - 1) ;
      interval global_ivl = interval(maximum_allocated, maximum_allocated+global_max-1) ; 
      maximum_allocated = local_max + size ;
      delete [] send_buf ;
      delete [] recv_buf ;
      delete [] size_send ;
      delete [] size_recv ;
      return(make_pair(local_ivl, global_ivl)) ;
    }
    interval alloc = interval(maximum_allocated,maximum_allocated+size-1) ;
    maximum_allocated += size ;
    return (make_pair(alloc, alloc)) ; ;
  }

  storeRepP fact_db::get_variable(variable v) {
    v = remove_synonym(v) ;
    std::map<variable, fact_info>::iterator mi =
      fmap.find(remove_synonym(v)) ;
    if(mi == fmap.end()) 
      return storeRepP(0) ;
    else
      return storeRepP(mi->second.data_rep) ;
  }
  
  fact_db::distribute_infoP fact_db::get_distribute_info() {
    return(distributed_info);
  }
  
  void fact_db::put_distribute_info(distribute_infoP dp) {
    distributed_info = dp ;
  }
 
  bool fact_db::isDistributed() {
    if(distributed_info == 0)
      return 0 ;
    else 
      return 1 ;
  }
  
  void fact_db::rotate_vars(const std::list<variable> &lvars) {
    list<variable>::const_iterator jj ;
    jj = lvars.begin() ;
    storeRepP cp = fmap[remove_synonym(*jj)].data_rep->getRep() ;
    ++jj ;
    if(jj != lvars.end()) {
      for(;jj!=lvars.end();++jj) {
        fact_info &fd = fmap[remove_synonym(*jj)] ;
        storeRepP tmp = fd.data_rep->getRep() ;
        fd.data_rep->setRep(cp) ;
        cp = tmp ;
      }
    }
    fmap[remove_synonym(lvars.front())].data_rep->setRep(cp) ;
  }

  ostream &fact_db::write(ostream &s) const {
    std::map<variable, fact_info>::const_iterator vmi ;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP storeRep = storeRepP(vmi->second.data_rep) ;
      entitySet en=storeRep->domain();
      std::string groupname = (v.get_info()).name;
      s << groupname << ":" ;
      storeRep->Print(s);
    }
    return s ;
  }

  istream &fact_db::read(istream &s) {
    string vname ;
    parse::kill_white_space(s) ;
    if(s.peek()!='{') {
      cerr << "format error in fact_db::read" << endl ;
      return s ;
    }
    s.get() ;
    
    for(;;) {
      parse::kill_white_space(s) ;
      if(s.peek() == '}') {
        s.get() ;
        break ;
      }
      if(s.peek() == char_traits<char>::eof()) {
        cerr << "unexpected EOF in fact_db::read" << endl ;
        exit(1) ;
      }
      parse::kill_white_space(s) ;
      if(parse::is_name(s)) 
        vname = parse::get_name(s) ;
      else {
        cerr << "syntax error in fact_db::read" << endl ;
        exit(1) ;
      }
      parse::kill_white_space(s) ;
      if(!parse::get_token(s,":")) {
        cerr << "syntax error in fact_db::read, no ':' separator"
             << endl ;
        exit(1) ;
      }


      variable var(vname) ;
      storeRepP vp = get_variable(var) ;
      if(vp == 0) {
        vp = get_variable_type(var) ;
        if(vp != 0) {
          create_fact(var,vp) ;
        }
        vp = get_variable(var) ;
      }
      if(vp == 0) {
        cerr << "variable named '" << vname
             << "' not found in database in fact_db::read." << endl
             << "Error not recoverable. " << endl ;
        exit(-1) ;
      }
      vp->Input(s) ;
    }
    return s ;
  }

  void fact_db::write_hdf5(const char *fname){
    hid_t  file_id, group_id;
    char   filename[500], str[100];
    std::map<variable, fact_info>::const_iterator vmi ;

    strcpy(filename, fname);
    strcat( filename, "_p");
    sprintf( str, "%d", Loci::MPI_rank);
    strcat( filename, str);
    strcat( filename, ".hdf5");

    file_id =  H5Fcreate( filename, H5F_ACC_TRUNC,
                          H5P_DEFAULT, H5P_DEFAULT);

    // Write Processor information ....
    hsize_t  dimension = 2;
    int      ibuf[]    = {0,1};
    int      rank      = 1;
    group_id         = H5Gcreate(file_id, "ProcessorID", 0);
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dcreate(group_id, "Processor", vDatatype, vDataspace,
                                  H5P_DEFAULT);
    ibuf[0] = Loci::MPI_rank;
    ibuf[1] = Loci::MPI_processes;
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuf);
    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

    entitySet eset;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v);
      eset += store_Rep->domain();
    }

    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v);
      entitySet en=store_Rep->domain();
      std::string groupname = (v.get_info()).name;
      group_id = H5Gcreate(file_id, groupname.c_str(), 0);
      (store_Rep->getRep())->writehdf5(group_id, en);
    }
    H5Fclose(file_id);
    H5Gclose(group_id);

  }

  void fact_db::read_hdf5(const char *fname){
    hid_t   group_id, group_id1, group_id2;
    char    filename[500], str[100];
    hid_t  *file_id;

    entitySet  fileSet;
    entitySet   eset, myLocalEntity, myGlobalEntity;
    entitySet::const_iterator ei, ci;

    file_id = new hid_t[Loci::MPI_processes];

    //-----------------------------------------------------------------
    // First know how many files, the previous execution created. This
    // is stored in the in the files. and since there was atleast one
    // processors, get this information from processor 0;
    //-----------------------------------------------------------------
    int     iproc, rank=1, ibuf[2], maxFiles;
    hsize_t dimension = 2;
    hid_t   vDataspace, vDatatype, vDataset;


    strcpy(filename, fname);
    if( Loci::MPI_rank == 0) {
      iproc = 0;
      strcat( filename, "_p");
      sprintf( str, "%d", iproc);
      strcat( filename, str);
      strcat( filename, ".hdf5");
      file_id[0] = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDatatype  = H5T_NATIVE_INT;
      group_id   = H5Gopen(file_id[0], "/ProcessorID");
      vDataset   = H5Dopen( group_id, "Processor");
      H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, ibuf);

      H5Sclose( vDataspace );
      H5Dclose( vDataset   );
      H5Fclose( file_id[0] );
      H5Gclose(group_id);

      maxFiles = ibuf[1];
    }

    //-----------------------------------------------------------------
    // Now decide which file(s) I can read from the pool of files, If
    // there are more processors than number of files, then everyone
    // can read only one file, otherwise some processors will read
    // more files than other, which can create redundant entities on
    // the processors, which need to be migrated to other processor.
    //-----------------------------------------------------------------

    if( Loci::MPI_processes > 1)
        MPI_Bcast( &maxFiles, 1, MPI_INT, 0, MPI_COMM_WORLD );

    std::map<variable, fact_info>::const_iterator vmi ;
    for( int iproc = Loci::MPI_rank; iproc < maxFiles; iproc+=Loci::MPI_processes) {
      strcpy(filename, fname);
      strcat( filename, "_p");
      sprintf( str, "%d", iproc);
      strcat( filename, str);
      strcat( filename, ".hdf5");
      file_id[iproc] = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    }

    Loci::fact_db::distribute_infoP d =
        Loci::exec_current_fact_db->get_distribute_info() ;


    entitySet usr_eset, currdom, localdom;
    Map       lg, currg2l, currl2l;

    if( isDistributed() ) currg2l = d->g2l;

    std::set<string>  donotRead;
    donotRead.insert("EMPTY");
    donotRead.insert("UNIVERSE");
    donotRead.insert("my_entities");
    donotRead.insert("l2g");

    map<string,entitySet>  gsetRead;

    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
        variable v=vmi->first;
        std::string vname = (v.get_info()).name;
        if( donotRead.find(vname) != donotRead.end()) continue;
     
        storeRepP store_Rep = get_variable(v)->getRep();
        currdom = EMPTY;
        std::string groupname = vname;
        for(int iproc=Loci::MPI_rank; iproc < maxFiles; iproc+=Loci::MPI_processes) {
            group_id1 = H5Gopen(file_id[iproc], "l2g");
            group_id2 = H5Gopen(file_id[iproc], groupname.c_str() );
            if( group_id1 > 0 && group_id2 > 0) {
                // Read the domain of the data. It is written with local numbering
                HDF5_ReadDomain(group_id2, eset);

                // Get the local->global numbering written in the file ..
                HDF5_Local2Global(group_id1, eset, lg);

                // We might have different global to local numbering in this run.
                // convert the global number to local number. and if the local number
                // is not assigned, use the global number for this run.
                for( ci = eset.begin(); ci != eset.end(); ++ci) {
                     gsetRead[vname] += lg[*ci];
                     if( currg2l.domain().inSet(lg[*ci]) )
                         currdom += currg2l[lg[*ci]];
                     else
                         currdom += lg[*ci];
                }
                H5Gclose(group_id1);
                H5Gclose(group_id2);
            }
        }
        store_Rep->allocate(currdom) ;
        for(int iproc=Loci::MPI_rank; iproc < maxFiles; iproc+=Loci::MPI_processes) {
          group_id1 = H5Gopen(file_id[iproc], "l2g");
          group_id2 = H5Gopen(file_id[iproc], groupname.c_str() );

          if( group_id1 > 0 && group_id2 > 0) {
            storeRepP tmpStore_Rep = store_Rep->new_store(usr_eset);
            localdom = ~EMPTY;
            tmpStore_Rep->readhdf5(group_id2, localdom);
            localdom = tmpStore_Rep->domain();
            HDF5_Local2Global(group_id1, localdom, lg);
            currl2l.allocate( localdom ); // local-local mapping from previous to current
            for( ci = localdom.begin(); ci != localdom.end(); ++ci) {
                if( currg2l.domain().inSet(lg[*ci]) )
                    currl2l[*ci] = currg2l[lg[*ci]];
                else
                    currl2l[*ci] = lg[*ci];
            }
            tmpStore_Rep = tmpStore_Rep->remap( currl2l );
            store_Rep->copy( tmpStore_Rep, tmpStore_Rep->domain());
            H5Gclose(group_id1);
            H5Gclose(group_id2);
          }
        }
        update_fact(v,store_Rep);
    }

    for( int iproc = Loci::MPI_rank; iproc < maxFiles; iproc+=Loci::MPI_processes) 
        H5Fclose(file_id[iproc]);

    delete [] file_id;

    if( Loci::MPI_processes == 1 || isDistributed() == 0) return;

    Map g2l, l2g;
    l2g = Loci::exec_current_fact_db->get_variable("l2g");
    g2l = d->g2l;
    myLocalEntity = d->my_entities;

    for( ei = myLocalEntity.begin(); ei != myLocalEntity.end(); ++ei)
        myGlobalEntity += l2g[*ei];

    //----------------------------------------------------------------
    // Now we are ready to allocate appropriate entities to each 
    // processors. Some of the entities might be read already and 
    // some will be extracted through redundant array.
    //----------------------------------------------------------------
    std::vector<int>  redundantSize(Loci::MPI_processes),
                      globalRedundant(Loci::MPI_processes),
                      displ(Loci::MPI_processes),
                      unpackSize(Loci::MPI_processes);
 
    std::vector<int>  isendbuf;
    std::vector<unsigned char> packbuf, unpackbuf;

    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
        variable v=vmi->first;
        std::string vname = (v.get_info()).name;
        if( donotRead.find(vname) != donotRead.end()) continue;


        storeRepP store_Rep = get_variable(v)->getRep();
        localdom = store_Rep->domain() & myLocalEntity;

        // get the global entity of the container

        // distribution will set the appropriate entities on this
        // processor. So find the redundant entities on this processor
        // which need to be migrated to other processor.

        entitySet myGlobalRedundant = gsetRead[vname] - myGlobalEntity;
        // Let everyone know that global entities each one sending
        int  numsend  =  myGlobalRedundant.size();
        MPI_Allgather( &numsend, 1, MPI_INT, &redundantSize[0], 1, 
                       MPI_INT, MPI_COMM_WORLD);

        int indx = 0;
        isendbuf.resize(numsend+1);
        for(ei=myGlobalRedundant.begin(); ei != myGlobalRedundant.end(); ++ei)
            isendbuf[indx++] = *ei;


        // All processor will be sending their "redundant entities". collect
        // all of them on each processor. Every processor will accept, if
        // the entity in the "redundant buffer" belongs to the processor,
        // otherwise, it will reject them..
        indx = 0;
        int  numTotal = 0;
        for( int i = 0; i < Loci::MPI_processes; i++) {
            displ[i]    = indx;
            indx       += redundantSize[i];
            numTotal   += redundantSize[i];
        }

       globalRedundant.resize(numTotal+1);
       MPI_Allgatherv( &isendbuf[0], numsend, MPI_INT, &globalRedundant[0], 
                       &redundantSize[0], &displ[0], MPI_INT, MPI_COMM_WORLD);

       entitySet localSet;
       sequence  localSeq;
       for( int i = 0; i < numTotal; i++)  {
            localSet  += g2l[globalRedundant[i]];
            localSeq  += g2l[globalRedundant[i]];
       }

       // Once the information about entities is known, we can pack the 
       // attributes of the container 
       int packSize = 0;
       packSize = store_Rep->pack_size(myGlobalRedundant);
       packbuf.resize(packSize+4);

       int location = 0;
       store_Rep->pack( &packbuf[0], location, packSize, localSet);
 
       MPI_Allgather( &packSize, 1, MPI_INT, &unpackSize[0], 1, MPI_INT, 
                      MPI_COMM_WORLD);
       indx     = 0;
       numTotal = 0;
       for( int i = 0; i < Loci::MPI_processes; i++) {
            displ[i]  = indx;
            indx     += unpackSize[i];
            numTotal += unpackSize[i];
       }

       unpackbuf.resize(numTotal+4);
       MPI_Allgatherv( &packbuf[0],  packSize, MPI_BYTE, &unpackbuf[0], &unpackSize[0], 
                       &displ[0], MPI_BYTE, MPI_COMM_WORLD);

       location = 0;
       store_Rep->allocate(myLocalEntity);
       store_Rep->unpack( &unpackbuf[0], location, numTotal, localSeq);
    }

  }



  
  void reorder_facts(fact_db &facts, Map &remap) {
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      if(facts.is_distributed_start())
	facts.replace_fact(*vi,p->remap(remap)) ;
      else
	facts.update_fact(*vi,p->remap(remap)) ;
    }
  }
  
  void serial_freeze(fact_db &facts) {
    variableSet vars = facts.get_typed_variables() ;
    entitySet map_entities ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      if(p->RepType() == MAP) {
        MapRepP mp = MapRepP(p->getRep()) ;
	entitySet dom = mp->domain() ;
        map_entities += dom ;
        map_entities += mp->image(dom) ;
      }
      if(p->RepType() == STORE) {
	map_entities += p->domain() ;
      }
    }
    Map m ;
    m.allocate(map_entities) ;
    for(entitySet::const_iterator ei = map_entities.begin(); ei != map_entities.end(); ++ei)
      m[*ei] = *ei ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      facts.replace_fact(*vi,p->remap(m)) ;
    }
  }

  void fact_db::set_variable_type(variable v, storeRepP st) {
    tmap[v] = storeRepP(st->new_store(EMPTY)) ;
  }

  storeRepP fact_db::get_variable_type(variable v) const {
    map<variable,storeRepP>::const_iterator mi ;
    if((mi=tmap.find(v)) != tmap.end())
      return storeRepP(mi->second->new_store(EMPTY)) ;
    else
      return storeRepP(0) ;
  }
}


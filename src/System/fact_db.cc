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

using std::istream ;
using std::ostream ;
using std::endl ;
using std::ios ;

#include <Tools/parse.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;
  extern fact_db *exec_current_fact_db;
  fact_db::fact_db() {
    distributed_info = 0 ;
    maximum_allocated = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      init_ptn.push_back(EMPTY) ;
    }
    exec_current_fact_db = this;
    dist_from_start = 0 ;
    constraint EMPTY_constraint ;
    EMPTY_constraint = EMPTY ;
    create_fact("EMPTY",EMPTY_constraint) ;
    constraint UNIVERSE_constraint ;
    UNIVERSE_constraint = ~EMPTY ;
    create_fact("UNIVERSE",UNIVERSE_constraint) ;
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

    // Sanity check, make sure v exists
    std::map<variable,fact_info>::iterator vmi, vmj ;
    if((vmi = fmap.find(v)) == fmap.end()) {
      cerr << "WARNING: synonym_variable("<<v<<","<<synonym<<")"<<endl ;
      cerr << "WARNING: type not known for target of synonym, ignoring operation" << endl ;
      //      abort() ;
      return ;
    }
    // If the synonym already points to a different variable instance,
    // remove it
    if((vmj = fmap.find(s)) != fmap.end()) {
      fact_info &finfo = vmj->second ;
      if((finfo.data_rep->domain() != EMPTY &&
          finfo.data_rep->RepType() != PARAMETER)) {
        cerr << "unable to define synonym variable " << synonym
             << " when variable already created in db. "  << endl ;
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
    variable tmp_v ;
    if(nspace_vec.size()) {  
      tmp_v = v ;
      for(unsigned int i = 0; i < nspace_vec.size(); ++i)
	tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    }
    else
      tmp_v = v ;
    
    warn(synonyms.find(tmp_v) != synonyms.end()) ;
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    
    if(mi != fmap.end()) {
      mi->second.data_rep->setRep(st->getRep()) ;
    } else
      cerr << "warning: update_fact: fact does not exist for variable " << tmp_v
	   << endl ;
  }
  
  void fact_db::create_fact(variable v, storeRepP st) {
    
    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }
    variable tmp_v ;
    if(nspace_vec.size()) {
      tmp_v = v ;
      for(unsigned int i = 0; i < nspace_vec.size(); ++i)
	tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    }
    else
      tmp_v = v ;
    if(synonyms.find(tmp_v) != synonyms.end()) {
      tmp_v = remove_synonym(tmp_v) ;
      std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
      if(mi==fmap.end()) {
        fmap[tmp_v].data_rep = new store_ref ;
        fmap[tmp_v].data_rep->setRep(st->getRep()) ;
      } else {
        if(typeid(st->getRep()) != typeid(mi->second.data_rep->getRep())) {
          cerr << "set_variable_type() method of fact_db changing type for variable " << tmp_v << endl ;
        }
        mi->second.data_rep->setRep(st->getRep()) ;
      }
      return ;
    }
    
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    if(mi != fmap.end()) {
      cerr << "WARNING: fact_db::set_variable_type retyping variable "
	   << tmp_v << endl ;
      mi->second.data_rep->setRep(st->getRep()) ;
    } else {
      fmap[tmp_v].data_rep = new store_ref ;
      fmap[tmp_v].data_rep->setRep(st->getRep()) ;
    }
    // cout << " tmp_v = " << tmp_v << endl ;
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
    

  std::pair<entitySet, entitySet> fact_db::get_distributed_alloc(int size) {
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
	if(size_recv[i] > 0 )
	  init_ptn[i] += interval(local, local+size_recv[i]-1) ;
      }
      entitySet local_ivl, global_ivl ;
      if(size > 0 )
	local_ivl = entitySet(interval(local_max, local_max + size - 1)) ;
      else {
	local_ivl = EMPTY ;
      }
      
      global_ivl = entitySet(interval(maximum_allocated, maximum_allocated+global_max-1)) ;
      if(size > 0)
	maximum_allocated = local_max + size ;
      delete [] send_buf ;
      delete [] recv_buf ;
      delete [] size_send ;
      delete [] size_recv ;
      return(make_pair(local_ivl, global_ivl)) ;
    }
    entitySet alloc = entitySet(interval(maximum_allocated,maximum_allocated+size-1)) ;
    maximum_allocated += size ;
    init_ptn[0] += alloc ;
    return (make_pair(alloc, alloc)) ;
  }

  storeRepP fact_db::get_variable(variable v) {
    variable tmp_v ;
    if(nspace_vec.size()) {  
      tmp_v = v ;
      for(unsigned int i = 0; i < nspace_vec.size(); ++i)
	tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    }
    else
      tmp_v = v ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, fact_info>::iterator mi =
      fmap.find(remove_synonym(tmp_v)) ;
    if(mi == fmap.end()) {
      if(Loci::MPI_rank == 0)
	cout << " returning null  storeRep for variable " << tmp_v<< endl ;
      return storeRepP(0) ;
      
    }
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

  /////////////////////////////////////////////////////////////////////////////
  
  void fact_db::write_hdf5(const char *fname, int num_partitions){
    hid_t  file_id, group_id;
    char   filename[500] ;
    std::map<variable, fact_info>::const_iterator vmi ;
    
    strcpy(filename, fname);
    //if(Loci::MPI_processes > 1) {
    //strcat( filename, "_");
    //sprintf( str, "%d", Loci::MPI_rank);
    //strcat( filename, str);
    //}
    //strcat( filename, ".hdf5");
    
    file_id =  H5Fcreate( filename, H5F_ACC_TRUNC,
                          H5P_DEFAULT, H5P_DEFAULT);
    num_partitions = max(num_partitions, Loci::MPI_processes);
    
    hsize_t  dimension = 1;
    int      rank      = 1;
    group_id         = H5Gcreate(file_id, "Partitions", 0);
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dcreate(group_id, "NumParts", vDatatype, vDataspace,H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &num_partitions);
    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

    //
    // The following facts are internal to Loci, which are not important from user's
    // point of view, so we will not write them into the facts database.
    //
    std::set<string>  donotWrite;
    donotWrite.insert("EMPTY");
    donotWrite.insert("UNIVERSE");
    donotWrite.insert("my_entities");

    // Domain of a contains contains both "MyEntitySet" and perhaphs clone
    // entities too. While writing we don't write clone entities, therefore
    // domain should be intersected with my "MyLocalEntitySet"

    entitySet myLocalEntity;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v);
      std::string groupname = (v.get_info()).name;
      if( donotWrite.find(groupname) != donotWrite.end()) continue;
      entitySet en=store_Rep->domain();
      if( isDistributed() )  {
        Loci::fact_db::distribute_infoP d;
        d   = Loci::exec_current_fact_db->get_distribute_info() ;
        myLocalEntity = d->my_entities;
        en = en & myLocalEntity;
      }
      group_id = H5Gcreate(file_id, groupname.c_str(), 0);
      (store_Rep->getRep())->writehdf5(group_id, en);
    }
    H5Fclose(file_id);
    H5Gclose(group_id);
  }

  /////////////////////////////////////////////////////////////////////////////

  void fact_db::read_hdf5(const char *fname){
    hid_t   group_id, group_id1, group_id2;
    char    filename[500], str[100];
    hid_t  *file_id;
    Map     entity_owner;
     
    entitySet::const_iterator ei;
    entitySet  fileSet, eset, myLocalEntity, globalEntitySet,
      myGlobalEntitySet, exportSet;
    sequence   importSeq;
   
    std::vector<int> displ(Loci::MPI_processes),
      unpacksize(Loci::MPI_processes), isendbuf, irecvbuf;
    
    MPI_Request  *request;
    MPI_Status   status;

    // Because of Asynchronous message passing, we have to store the buffer
    // till the communication is over...
    struct Message {
      int            isendbuf1, isendbuf2;
      int           *isendbuf;
      unsigned char *sendbuf;
    };
    Message  *message;

    std::vector<unsigned char> packbuf, unpackbuf;
    
    file_id = new hid_t[Loci::MPI_processes];
    message = new Message[Loci::MPI_processes];
    request = new MPI_Request[Loci::MPI_processes];

    H5Eset_auto (NULL, NULL);

    for( int i = 0; i < Loci::MPI_processes; i++) {
      message[i].isendbuf = NULL;
      message[i].sendbuf  = NULL;
    }
    //-----------------------------------------------------------------
    // first know how many files, the previous execution created. this
    // is stored in the in the files. and since there was atleast one
    // processors, get this information from processor 0;
    //-----------------------------------------------------------------
    int     ifile, rank=1, maxfiles=1, gid, numsend, numrecv, packsize, location, num_parts;
    hsize_t dimension = 2;
    hid_t   vdataspace, vdatatype, vdataset;

    strcpy(filename, fname);
    if( Loci::MPI_rank == 0) {
      ifile = 0;
      if(Loci::MPI_processes > 1) {
	strcat( filename, "_");
	sprintf( str, "%d", ifile);
	strcat( filename, str);
      }
      strcat( filename, ".hdf5");
      file_id[0] = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if( file_id[0] > 0) {
        vdataspace = H5Screate_simple(rank, &dimension, NULL);
        vdatatype  = H5T_NATIVE_INT;
        group_id   = H5Gopen(file_id[0], "/Partitions");
        if( group_id > 0) {
          vdataset   = H5Dopen( group_id, "NumParts");
          H5Dread( vdataset, vdatatype, H5S_ALL, vdataspace, H5P_DEFAULT, &num_parts);
          H5Sclose( vdataspace );
          H5Dclose( vdataset   );
          H5Fclose( file_id[0] );
          H5Gclose( group_id );

          maxfiles = num_parts;
 
        }
      } else {
        cout << "Warning: Unable to open file " << filename << endl;
        return;
      }
    }


    //-----------------------------------------------------------------
    // now decide which file(s) i can read from the pool of files, if
    // there are more processors than number of files, then everyone
    // can read only one file, otherwise some processors will read
    // more files than other, which can create redundant entities on
    // the processors, which need to be migrated to other processor.
    //-----------------------------------------------------------------

    if( Loci::MPI_processes > 1)
      MPI_Bcast( &maxfiles, 1, MPI_INT, 0, MPI_COMM_WORLD );

    // Which files, I am supposed to read from the pool of files.
    std::vector<int> files_assigned;
    for(int I=Loci::MPI_rank; I < maxfiles; I+=Loci::MPI_processes)
      files_assigned.push_back(I);

    std::map<variable, fact_info>::const_iterator vmi ;
    for(unsigned int I=0;I < files_assigned.size(); I++){
      strcpy(filename, fname);
      if(Loci::MPI_processes > 1) {
	strcat( filename, "_");
	sprintf( str, "%d", files_assigned[I]);
	strcat( filename, str);
      }
      strcat( filename, ".hdf5");
      cout << " reading file_name " << filename << endl ;
      file_id[I] = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if( file_id[I] < 0) {
        cout << "Warning: Couldn't open file " << filename << endl;
        return;
      }
    }
    entitySet  new_eset,currdom,localdom, gsetRead,retainSet;
    Map        l2g, lg, currl2l;
    dMap currg2l ;
    Loci::fact_db::distribute_infoP d;
    if( isDistributed() )  {
      d       = Loci::exec_current_fact_db->get_distribute_info() ;
      currg2l = d->g2l;
      l2g     = Loci::exec_current_fact_db->get_variable("l2g");
      myLocalEntity = d->my_entities;

      // currg2l's domain may contain "clone" entities too, therefore use
      // myLocalEntitySet to know what are the global entitySet assigned to
      // this process.
      eset = EMPTY;
      for( ei = myLocalEntity.begin(); ei != myLocalEntity.end(); ++ei)
        eset += l2g[*ei];

      numsend  =  2*eset.size();
      MPI_Allgather( &numsend, 1, MPI_INT, &unpacksize[0], 1, MPI_INT, 
                     MPI_COMM_WORLD);
      int indx = 0;
      isendbuf.resize(numsend+1);
      for( ei = eset.begin(); ei != eset.end(); ++ei){
        isendbuf[indx++] = *ei;
        isendbuf[indx++] = Loci::MPI_rank;
      }

      int  numtotal = 0;
      for( int i = 0; i < Loci::MPI_processes; i++) {
        displ[i]    = numtotal;
        numtotal   += unpacksize[i];
      }

      irecvbuf.resize(numtotal+1);
      MPI_Allgatherv( &isendbuf[0], numsend, MPI_INT, &irecvbuf[0], 
                      &unpacksize[0], &displ[0], MPI_INT, MPI_COMM_WORLD);

      eset = EMPTY;
      for( int i = 0; i < numtotal/2; i++) eset += irecvbuf[2*i];

      entity_owner.allocate(eset);
      int iproc;
      for( int i = 0; i < numtotal/2; i++){
        gid     = irecvbuf[2*i];
        iproc   = irecvbuf[2*i+1];
        entity_owner[gid] = iproc; 
      }
    }

    std::set<string>  donotRead;
    donotRead.insert("EMPTY");
    donotRead.insert("UNIVERSE");
    donotRead.insert("my_entities");
    donotRead.insert("l2g");

    Map  localmap, scatter_map;
    storeRepP tmpstore_Rep;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      std::string vname = (v.get_info()).name;
      if( donotRead.find(vname) != donotRead.end()) continue;
     
      storeRepP store_Rep = get_variable(v)->getRep();

      // parameter and constraints are specical case and handle separately

      if( store_Rep->RepType() == PARAMETER || 
          store_Rep->RepType() == CONSTRAINT ) continue;
      gsetRead   = EMPTY;
      std::string groupname = vname;

      // first get information about entitset and assign local number to them
   
      for(unsigned int I = 0; I < files_assigned.size(); I++){
        group_id2 = H5Gopen(file_id[I], groupname.c_str() );  
        if( group_id2 < 0) continue;

        HDF5_ReadDomain(group_id2, eset);
        H5Gclose(group_id2);
        // get the local->global numbering written in the file ..
	if( maxfiles > 1) {
	group_id1 = H5Gopen(file_id[I], "l2g");           
	if( group_id1 > 0) {
	  HDF5_Local2Global(group_id1, eset, lg);
	  H5Gclose(group_id1);         
	  for( ei = eset.begin(); ei != eset.end(); ++ei)
	    gsetRead += lg[*ei];
	} else
	  gsetRead += eset;
	} else 
	  gsetRead += eset;
      }
      // First of all, if the facts are distributed, then not all "globalEntitySet"
      // need to be allocated for each container. Let everyone know what are the 
      // global entities for this container and the required entitySet for this
      // container will be difference of "globalEntitySet" and "gsetRead".
      if( isDistributed() ){       
        numsend  =  gsetRead.size();
        MPI_Allgather( &numsend, 1, MPI_INT, &unpacksize[0], 1, MPI_INT, 
                       MPI_COMM_WORLD);

        int indx = 0;
        isendbuf.resize(numsend+1);
        for(ei=gsetRead.begin(); ei != gsetRead.end(); ++ei)
          isendbuf[indx++] = *ei;

        int  numtotal = 0;
        for( int i = 0; i < Loci::MPI_processes; i++) {
          displ[i]    = numtotal;
          numtotal   += unpacksize[i];
        }

        irecvbuf.resize(numtotal+1);
        MPI_Allgatherv( &isendbuf[0], numsend, MPI_INT, &irecvbuf[0], 
                        &unpacksize[0], &displ[0], MPI_INT, MPI_COMM_WORLD);

        globalEntitySet =  EMPTY;
        for( int i = 0; i < numtotal; i++) 
          globalEntitySet +=  irecvbuf[i];

        myGlobalEntitySet = EMPTY;
        for( ei = globalEntitySet.begin(); ei != globalEntitySet.end(); ++ei){
          if(getEntityIdentity( entity_owner[*ei] ) == Loci::MPI_rank)
            myGlobalEntitySet += *ei;
        }
        
      } else {
        globalEntitySet   = gsetRead;
        myGlobalEntitySet = gsetRead;
      }
 
      // For all entitities of the container find out the global to local numbering.
      // If the fact is already distributed, then some of the entities might already
      // get the local numbering on a given processor. In case, a global entitySet is
      // not assigned a local numbering, we are at liberty to choose of our own for
      // this module. ( This might change later by calling global to local numbering
      // module.
      currdom  = EMPTY;
      new_eset = globalEntitySet;
      localmap.allocate(new_eset);
      if( isDistributed() ){       
        int indx  = 0;
        for(ei = new_eset.begin(); ei != new_eset.end(); ++ei) {
          if(currg2l.domain().inSet(*ei) ) {
            currdom       += currg2l[*ei];
            localmap[*ei]  = currg2l[*ei];
          }else
            localmap[*ei]  = indx++;
        }
      } else {
        myLocalEntity     = gsetRead;
        currdom           = gsetRead;
        for(ei = new_eset.begin(); ei != new_eset.end(); ++ei) 
          localmap[*ei] = *ei;
      }

      // divide the entity set into two group, retainSet and forwardSet
      retainSet  = gsetRead & myGlobalEntitySet;

      store_Rep->allocate(currdom);
      new_eset = globalEntitySet;
      int has_file, file_counter;

      std::vector<int> :: const_iterator file_iter,  file_end;
      file_iter  = files_assigned.begin();
      file_end   = files_assigned.end();
      
      ifile = -1;
      while( 1 ) {
        
        has_file = 0;
        if( file_iter != file_end ) {
          ifile++;
          has_file = 1;
          ++file_iter;
        }

        if( Loci::MPI_processes > 1)
          MPI_Allreduce( &has_file, &file_counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        else
          file_counter = has_file;
	
        if( file_counter == 0) break;

        localdom = EMPTY; 
        if( has_file ) {
          group_id2 = H5Gopen(file_id[ifile], groupname.c_str() );
          // read the attribute data in temporary buffer. Remember that entities are in
          // local numbering in the file.     
          if( group_id2 > 0) {
            tmpstore_Rep = store_Rep->new_store(new_eset);
            localdom = ~EMPTY;
            tmpstore_Rep->readhdf5(group_id2, localdom);
            localdom = tmpstore_Rep->domain();            
            H5Gclose(group_id2);
	    lg.allocate(localdom);
            // if more than one file, then local-global mapping has to be read
	    group_id1 = H5Gopen(file_id[ifile], "l2g");           
	    if( group_id1 > 0) {
	      HDF5_Local2Global(group_id1, localdom, lg);
	      H5Gclose(group_id1);         
	    } else {
              for( ei = localdom.begin(); ei != localdom.end(); ++ei)
                lg[*ei] = *ei;
            }

            // A Composition is needed which can map entities from the local numbering in
            // the previous run to the local numbering in the present run..
          
            scatter_map.allocate( localdom );
          
            currdom = EMPTY;
            for( ei = localdom.begin(); ei != localdom.end(); ++ei) {
              if( retainSet.inSet(lg[*ei])) {
                currdom          +=  *ei;
                scatter_map[*ei]  =  localmap[lg[*ei]];
              }
            }
	    storeRepP scatter_sp = MapRepP(scatter_map.Rep())->thaw() ;
	    dMap d_scatter_map(scatter_sp) ;
	    store_Rep->scatter( d_scatter_map, tmpstore_Rep, currdom );
	  }
        }
	
        if( !isDistributed() ) continue;

        // Now send the entities which are not owned by this process. There
        // four messages with tags(1...4)
        for(int jproc=0; jproc < Loci::MPI_processes; jproc++) {
          exportSet = EMPTY;
          for( ei = localdom.begin(); ei != localdom.end(); ++ei) 
            if(getEntityIdentity(entity_owner[lg[*ei]]) == jproc) 
	      exportSet += *ei;
    
          // First message with Tag=1 contains number of entitySet send to
          // remote processor
          numsend = exportSet.size();
          message[jproc].isendbuf1 = numsend;
          MPI_Isend( &message[jproc].isendbuf1, 1, MPI_INT, jproc, 1, 
                     MPI_COMM_WORLD, &request[jproc] );

          if( numsend ) {
            // Second message contains the entitySet (with global number)
            // send to remote processor.
            message[jproc].isendbuf = new int[numsend];
            int indx = 0;
            for(ei=exportSet.begin(); ei != exportSet.end(); ++ei)
              message[jproc].isendbuf[indx++] = lg[*ei];
            MPI_Isend( message[jproc].isendbuf, numsend, MPI_INT, 
                       jproc, 2, MPI_COMM_WORLD, &request[jproc] );

            // Third message with TAG=3 contains the size of packet that will
            // be send to the processor. This cann't be calculated by the remote
            // processor, because it may not know the contains of the objects.
            // so this message will be used to memory allocation on the remote
            // processor. In fact, MPI_PROBE could be used to decide the message
            // size, but possibly, Probing is inefficient, so I am avoiding this
            // function.
            packsize = tmpstore_Rep->pack_size(exportSet);
          
            message[jproc].isendbuf2 = packsize;
            MPI_Isend( &message[jproc].isendbuf2, 1, MPI_INT, 
                       jproc, 3, MPI_COMM_WORLD, &request[jproc]);

            // Fourth Message with TAG=4 value actual contains the attribute data of
            // the container which need to be passed to the remote processor. At
            // present we are sending data as BYTES, which may change in future.
            // 
            if( packsize ) {
              location = 0;
              message[jproc].sendbuf = new unsigned char[packsize];
              tmpstore_Rep->pack( message[jproc].sendbuf, location, 
                                  packsize, exportSet);
              MPI_Isend(message[jproc].sendbuf, packsize, MPI_BYTE, 
                        jproc, 4, MPI_COMM_WORLD, &request[jproc]);
            }
          }
        }
        //
        // Now receive the data from the remote process.
        //
        for(int jproc=0; jproc < Loci::MPI_processes; jproc++) {
          MPI_Recv( &numrecv, 1, MPI_INT, jproc, 1, MPI_COMM_WORLD, 
                    &status);
          if( numrecv ) {
            irecvbuf.resize(numrecv);
            MPI_Recv( &irecvbuf[0], numrecv, MPI_INT, jproc, 2, MPI_COMM_WORLD, 
                      &status);
            importSeq = EMPTY;
            for( int i = 0; i < numrecv; i++)
              importSeq += currg2l[irecvbuf[i]];
   
            MPI_Recv( &unpacksize[jproc], 1, MPI_INT, jproc, 3, MPI_COMM_WORLD, 
                      &status);
            if( unpacksize[jproc] ) {
              unpackbuf.resize(unpacksize[jproc]);

              MPI_Recv( &unpackbuf[0], unpacksize[jproc], MPI_BYTE, jproc, 4, MPI_COMM_WORLD, 
                        &status);
              location = 0;
              store_Rep->unpack( &unpackbuf[0], location, unpacksize[jproc], importSeq);
            }
          }
        }
        for(int jproc=0; jproc < Loci::MPI_processes; jproc++) {
          if( message[jproc].isendbuf ) {
            delete [] message[jproc].isendbuf;
            message[jproc].isendbuf = NULL;
          }
          if( message[jproc].sendbuf ) {
            delete [] message[jproc].sendbuf;
            message[jproc].sendbuf = NULL;
          }
        }
       
      } // Complete all the containers
      update_fact(v,store_Rep);
    }
    
     
    // parameter and constraints are specical case and handle separately
    // All the entitySet are with respect to global numbering evern for distributed
    // memory
  
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      std::string vname = (v.get_info()).name;
      if( donotRead.find(vname) != donotRead.end()) continue;

      storeRepP store_Rep = get_variable(v)->getRep();

      if( store_Rep->RepType() == PARAMETER || 
          store_Rep->RepType() == CONSTRAINT ) {

        std::string groupname = vname;
      
        gsetRead = EMPTY;
        for(unsigned int I = 0; I < files_assigned.size(); I++){
          group_id2 = H5Gopen(file_id[I], groupname.c_str() );  
          if( group_id2 < 0) continue;
          HDF5_ReadDomain(group_id2, eset);
          H5Gclose(group_id2);
          gsetRead += eset;
        }
	
        new_eset = EMPTY;
        store_Rep = store_Rep->new_store(new_eset);
        store_Rep->allocate(gsetRead);
      
        
	for(int proc = 0; proc < Loci::MPI_processes; ++proc) {
          group_id2 = H5Gopen(file_id[proc], groupname.c_str() );
          if( group_id2 > 0)
            store_Rep->readhdf5(group_id2, gsetRead);
        }
	// Only one processor is sufficient to read the parameter or constraints
        // and we can broadcast these balues.
	
	/*
        if(Loci::MPI_rank == 0 || Loci::MPI_processes == 1 ) {
          group_id2 = H5Gopen(file_id[0], groupname.c_str() );
          if( group_id2 > 0)
            store_Rep->readhdf5(group_id2, gsetRead);
        }
	
        if( Loci::MPI_processes > 1 ) {
          packsize = store_Rep->pack_size(gsetRead);
          MPI_Bcast(&packsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
          
          location = 0;
          packbuf.resize(packsize);
          store_Rep->pack( &packbuf[0], location, packsize, ~EMPTY);
          MPI_Bcast(&packbuf[0], packsize, MPI_BYTE,  0, MPI_COMM_WORLD );
          
          location = 0;
          store_Rep->unpack( &packbuf[0], location, packsize, ~EMPTY);
	}
	*/
        update_fact(v,store_Rep);        
      }
    }

  

    for(unsigned int I = 0; I < files_assigned.size(); I++)
      H5Fclose(file_id[I]);

    delete [] file_id;
    delete [] message;
    delete [] request;
  }


  void fact_db::Print_diagnostics() {
    std::map<variable, fact_info>::iterator mi ;
    ostringstream oss ;
    oss << "memory." ;
    oss << MPI_rank ;
    string file_name = oss.str() ;
    std::ofstream ofile(file_name.c_str(), ios::out) ;
    double total_size = 0 ;
    double total_wasted = 0 ;
    entitySet dom, total, unused ;
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(finfo.data_rep->RepType() == STORE) {
	dom = finfo.data_rep->domain() ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += finfo.data_rep->pack_size(dom) ; 
	total_wasted += finfo.data_rep->pack_size(unused) ;
      }
    }
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(finfo.data_rep->RepType() == STORE) {
	dom = finfo.data_rep->domain() ;
	double size = finfo.data_rep->pack_size(dom) ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	double wasted_space = finfo.data_rep->pack_size(unused) ;
	ofile << " ****************************************************" << endl ;
	ofile << " Total_size = " << total_size << endl ;
	ofile << "Variable = "  << mi->first << endl ;
	ofile << "Domain = " << dom << endl ;
	ofile << "Size allocated = " << size << endl ; 
	if( isDistributed() )  {
	  Loci::fact_db::distribute_infoP d ;
	  d   = Loci::exec_current_fact_db->get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = finfo.data_rep->pack_size(clone) ;
	  ofile << "----------------------------------------------------" << endl;
	  ofile << " My_entities = " << my_entities << endl ;
	  ofile << " Clone entities = " << clone << endl ;
	  ofile << "Memory required for the  clone region  = " << clone_size << endl ;
	  ofile << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  ofile << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  ofile << "----------------------------------------------------" << endl;
	}
	ofile << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	ofile << "----------------------------------------------------" << endl;
	ofile << "Total wasted size = " << total_wasted << endl ;
	ofile << "Unused entities = " << unused << endl ;
	ofile << "Wasted space = " << wasted_space << endl ;
	ofile << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	ofile << " ***************************************************" << endl << endl << endl ;
      }
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  
  void reorder_facts(fact_db &facts, dMap &remap) {
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      if(facts.is_distributed_start())
        facts.replace_fact(*vi,p->remap(remap)) ;
      else
        facts.update_fact(*vi,p->remap(remap)) ;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  
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
    dMap m ;
    // No need to allocate, dMap is self allocating
    //    m.allocate(map_entities) ;
    for(entitySet::const_iterator ei = map_entities.begin();
        ei != map_entities.end();
        ++ei) {
      m[*ei] = *ei ;
    }

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

  

#include <fact_db.h>
#include <constraint.h>
#include <Tools/stream.h>
#include <Tools/debugger.h>

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
  }

  fact_db::~fact_db() {}

  void fact_db::synonym_variable(variable v, variable synonym) {
    v = remove_synonym(v) ;
    std::map<variable, fact_info>::iterator vmi ;
    if((vmi = fmap.find(synonym)) != fmap.end()) {
      fact_info &finfo = vmi->second ;
      if((finfo.data_rep->domain() != EMPTY &&
          finfo.data_rep->RepType() != PARAMETER)) {
        cerr << "unable to define synonym variable " << synonym
             << " when varaiable already created in db. "  << endl ;
        cerr << "variable v = " << v << endl ;
        abort() ;
      }
      // This is a hack
      remove_variable(synonym) ;
    }
    synonyms[synonym] = v ;
  }
  
  
  void fact_db::create_fact(variable v, storeRepP st) {
    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }
    set_variable_type(v,st) ;
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
  
  void fact_db::set_variable_type(variable v, storeRepP st) {
    //    warn(synonyms.find(v) != synonyms.end()) ;
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
      synonyms.erase(si) ;
    } else if((mi=fmap.find(v)) != fmap.end()) {
      mi->second.data_rep = 0 ;
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
	// debugout << "Partition for processor  " << i <<"   = " << init_ptn[i] << endl ;
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

      storeRepP vp = get_variable(vname) ;
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

  void fact_db::write_hdf5(const char *filename){
    hid_t  file_id, group_id;
    file_id =  H5Fcreate( filename, H5F_ACC_TRUNC,
                          H5P_DEFAULT, H5P_DEFAULT);
    std::map<variable, fact_info>::const_iterator vmi ;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v);
      entitySet en=store_Rep->domain();
      std::string groupname = (v.get_info()).name;
      cout<<"Write "<<groupname<<" to HDF5 file "<<endl;
      group_id = H5Gcreate(file_id, groupname.c_str(), 0);
      (store_Rep->getRep())->writehdf5(group_id, en);
    }

  }

  void fact_db::read_hdf5(const char *filename){
    hid_t  file_id, group_id;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    std::map<variable, fact_info>::const_iterator vmi ;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v)->getRep();
      std::string groupname = (v.get_info()).name;
      cout<<"Read "<<groupname<<" from HDF5 file "<<endl;
      if(H5Gopen(file_id,groupname.c_str())>0){
        group_id = H5Gopen(file_id, groupname.c_str());
        entitySet dom = store_Rep->domain() ;
        store_Rep->readhdf5(group_id, dom);
        update_fact(v,store_Rep);
      }
      else
        cerr<<("Warning: variable \""+groupname+"\" is not found in file \""+filename+"\"")<<endl;
      //store_Rep->Print(cout);
    }
  }

  
  void reorder_facts(fact_db &facts, Map &remap) {
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      //facts.create_fact(*vi,p->remap(remap)) ;
      facts.update_fact(*vi,p->remap(remap)) ;
    }
  }
}


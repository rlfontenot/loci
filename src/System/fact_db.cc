#include <fact_db.h>
#include <constraint.h>
#include <Tools/stream.h>
#include <Tools/debugger.h>
#include <DStore.h>
#include "dist_tools.h"
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
    minimum_allocated = -1 ;
    for(int i = 0; i < MPI_processes; ++i) {
      init_ptn.push_back(EMPTY) ;
    }
    //    exec_current_fact_db = this;
    dist_from_start = 0 ;
    constraint EMPTY_constraint ;
    EMPTY_constraint = EMPTY ;
    create_fact("EMPTY",EMPTY_constraint) ;
    constraint UNIVERSE_constraint ;
    UNIVERSE_constraint = ~EMPTY ;
    create_fact("UNIVERSE",UNIVERSE_constraint) ;
  }

  fact_db::~fact_db() {}
  void fact_db::set_maximum_allocated(int i) {
    maximum_allocated = i ;
  }
  void fact_db::set_minimum_allocated(int i) {
    minimum_allocated = i ;
  }
  
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
      for(size_t i = 0; i < nspace_vec.size(); ++i)
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
      for(size_t i = 0; i < nspace_vec.size(); ++i)
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
      for(size_t i = 0; i < nspace_vec.size(); ++i)
	tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    }
    else
      tmp_v = v ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, fact_info>::iterator mi =
      fmap.find(remove_synonym(tmp_v)) ;
    if(mi == fmap.end()) {
      //      if(Loci::MPI_rank == 0)
      //	cout << " returning null  storeRep for variable " << tmp_v<< endl ;
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
  void fact_db::write_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    write_hdf5(filename, vars) ;
  }
  void fact_db::read_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    read_hdf5(filename, vars) ; 
  }
  void fact_db::write_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0, group_id=0;
    if(Loci::MPI_rank == 0) 
      file_id =  H5Fcreate(filename, H5F_ACC_TRUNC,
			   H5P_DEFAULT, H5P_DEFAULT) ;
    
    
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(p->RepType() == STORE) {
	if(MPI_rank == 0)
	  group_id = H5Gcreate(file_id, (variable(*vi).get_info().name).c_str(), 0) ;
	
	if(isDistributed()) {
	  fact_db::distribute_infoP df = get_distribute_info() ;
	  dMap remap = df->remap ;
	  storeRepP reorder_sp = collect_reorder_store(p, remap, *this) ;
	  write_container(group_id, reorder_sp) ;
	} else 
	  write_container(group_id, p->getRep()) ;
	if(MPI_rank == 0)
	  H5Gclose(group_id) ;
      }
    }
    if(Loci::MPI_rank == 0) 
      H5Fclose(file_id) ;
  }
  
  
  void fact_db::read_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0, group_id=0;
    if(Loci::MPI_rank == 0) 
      file_id =  H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT) ;
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(p->RepType() == STORE) {
	if(Loci::MPI_rank == 0) 
	  group_id = H5Gopen(file_id, (variable(*vi).get_info().name).c_str()) ;
	entitySet dom = p->domain() ;
	read_container(group_id, p, dom) ;
	if(Loci::MPI_rank == 0) 
	  H5Gclose(group_id) ;
      }
    }
    if(Loci::MPI_rank == 0) 
      H5Fclose(file_id) ;
    
  }
}

  

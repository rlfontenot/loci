#include <fact_db.h>
#include <constraint.h>
#include <Tools/stream.h>
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
  
  fact_db::fact_db() {
    constraint EMPTY ;
    create_fact("EMPTY",EMPTY) ;
    constraint UNIVERSE ;
    UNIVERSE = ~EMPTY ;
    create_fact("UNIVERSE",UNIVERSE) ;
    distributed_info == 0 ;
  }

  fact_db::~fact_db() {}

  namespace {
    inline bool offset_sort(const variable &v1, const variable &v2)
    { return v1.get_info().offset > v2.get_info().offset ; }
  }
  
  void fact_db::register_variable(variable v) {
    time_ident vtime = v.time() ;
    if(!all_vars.inSet(v) && vtime != time_ident() &&
       !v.get_info().assign && (v.get_info().priority.size() == 0)) {
      const string &name = v.get_info().name ;
      list<variable> s ;
      s.push_back(v) ;
      time_map[vtime][name].merge(s,offset_sort) ;
    }
    all_vars += v ;
  }
  
  void fact_db::install_fact_info(variable v, fact_info info) {
    vmap_type::iterator mi = vmap.find(v) ;
    if(mi != vmap.end()) {
      cerr << "error:  reinstalling fact in database for variable "
           << v << endl ;
      abort() ;
    }
    info.synonyms += v ;
    vmap[v] = info ;
    fact_infov[info.fact_info_ref].aliases += v ;
    register_variable(v) ;
  }
  
  void fact_db::install_fact_data(variable v, fact_data data) {
    if(free_set == EMPTY) {
      fact_infov.push_back(data) ;
      install_fact_info(v,fact_info(fact_infov.size()-1)) ;
    } else {
      int val = *(free_set.begin()) ;
      free_set -= val ;
      fact_infov[val] = data ;
      install_fact_info(v,fact_info(val)) ;
    }
  }

  
  fact_db::fact_info &fact_db::get_fact_info(variable v) {
    vmap_type::iterator mi = vmap.find(remove_synonym(v)) ;
    if(mi == vmap.end()) {
      //    cerr << "error: variable " << v << " does not exist in fact database"
      //         << endl ;
      //      exit(-1) ;
      //    return vmap.begin()->second ;
      return get_fact_info(variable("EMPTY")) ;
    }
    return mi->second ;
  }

  void fact_db::create_fact(variable v, storeRepP st) {
    set_variable_type(v,st) ;
    variable_is_fact_at(v,st->domain()) ;
  }

  void fact_db::update_fact(variable v, storeRepP st) {
    if(all_vars.inSet(v)) {
      fact_data &fd = fact_infov[get_fact_info(v).fact_info_ref] ;
      (*fd.data_rep) = st ;
      fd.imageMap.clear() ;
      fd.preimageMap.clear() ;
      fd.ismap = (st->RepType() == Loci::MAP) ;
      if(fd.ismap)
        fd.minfo = MapRepP(fd.data_rep->getRep()) ; 
      //    (get_fact_data(v)) = fact_data(v,st) ;
      fact_info &fi = get_fact_info(v) ;
      fi.fact_installed = st->domain() ;
      fi.existence = st->domain() ;
    } else
      cerr << "warning: update_fact: fact does not exist for variable " << v
           << endl ;
  }

  void fact_db::set_variable_type(variable v, storeRepP st) {
    if(all_vars.inSet(v)) {
      cerr << "WARNING: fact_db::set_variable_type retyping variable "
           << v << endl ;
      get_fact_data(v).data_rep->setRep(st) ;
    } else
      install_fact_data(v,fact_data(v,st)) ;
  }

  void fact_db::allocate_variable(variable v, entitySet s) {
    cerr << "allocate_variable not implemented " << endl ;
    exit(-1) ;
  }
  
  void fact_db::variable_is_fact_at(variable v,entitySet s) {
    fact_info &fi = get_fact_info(v) ;
    fi.fact_installed += s ;
    fi.existence += s ;
    variableSet aliases = fact_infov[fi.fact_info_ref].aliases ;
    aliases -= v ;
    for(variableSet::const_iterator vi=aliases.begin();vi!=aliases.end();++vi) 
      get_fact_info(v).fact_installed -= s ;
  }

  void fact_db::alias_variable(variable v, variable alias) {
    if(all_vars.inSet(v)) {
      if(all_vars.inSet(alias)) {
        fact_info &vinfo = get_fact_info(v) ;
        fact_info &ainfo = get_fact_info(alias) ;
        if(vinfo.fact_info_ref != ainfo.fact_info_ref) {
          // merge alias and v
          warn(fact_infov[ainfo.fact_info_ref].data_rep->domain() != EMPTY) ;
          int del = ainfo.fact_info_ref ;
          // move aliases from record to be deleted
          variableSet move_aliases = fact_infov[del].aliases ;
          // change reference number to that of current record
          variableSet::const_iterator vi ;
          for(vi=move_aliases.begin();vi!=move_aliases.end();++vi) 
            get_fact_info(*vi).fact_info_ref = vinfo.fact_info_ref ;
          // update new record alias set
          fact_infov[vinfo.fact_info_ref].aliases += move_aliases ;
          // delete record and save for recovery
          fact_infov[del] = fact_data() ;
          free_set += del ;
        }
        return ;
      }

      int ref = get_fact_info(v).fact_info_ref ;
      fact_infov[ref].aliases += alias ;
      install_fact_info(alias,fact_info(ref)) ;
      
    } else if(all_vars.inSet(alias)) {
      alias_variable(alias,v) ;
    } else {
      cerr << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
      abort() ;
    }
  }

  void fact_db::synonym_variable(variable v, variable synonym) {
    v = remove_synonym(v) ;
    vmap_type::iterator vmi ;
    if((vmi = vmap.find(synonym)) != vmap.end()) {
      const fact_info &finfo = vmi->second ;
      fact_data &fdata = fact_infov[finfo.fact_info_ref] ;
      if(finfo.fact_installed != EMPTY ||
         (fdata.data_rep->domain() != EMPTY &&
          fdata.data_rep->RepType() != PARAMETER)) {
        cerr << "unable to define synonym variable " << synonym
             << " when varaiable already created in db. "  << endl ;
        cerr << "variable v = " << v << endl ;
        cerr << "finfo.fact_installed == " << finfo.fact_installed << endl ;
        cerr << "fdata.aliases = " ;
        for(variableSet::const_iterator vi=fdata.aliases.begin();
            vi!=fdata.aliases.end();++vi) {
          cerr << *vi << " " << endl ;
        }
        
        //      abort() ;
      }
      fact_info &vfinfo = get_fact_info(v) ;
      if(vfinfo.fact_info_ref != finfo.fact_info_ref) {
        // if variables exist and they are not synonymous, alias them first
        alias_variable(v,synonym) ;
      }

      // Join synonym mappings into one. = v ;

      // Remove variable from vmap structure
      fact_infov[vfinfo.fact_info_ref].aliases -= synonym ;
      vmap.erase(vmi) ;
    }
    fact_info &vfinfo = get_fact_info(v) ;

    vfinfo.synonyms += synonym ;
    synonyms[synonym] = v ;
    register_variable(synonym) ;
  }
    
  
  storeRepP fact_db::get_variable(variable v) {
    if(all_vars.inSet(v))
      return storeRepP(get_fact_data(v).data_rep) ;
    else
      return storeRepP(0) ;
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
  
  fact_db::time_infoP fact_db::get_time_info(time_ident tl) {
    time_infoP ti = new time_info ;
    if(tl == time_ident())
      return ti ;
    variable time_var(tl) ;
    
    if(get_typed_variables().inSet(time_var))
      ti->time_var = get_variable(time_var) ;
    else
      set_variable_type(time_var,ti->time_var.Rep()) ;
    typedef map<std::string,std::list<variable> >  maptype ;
    maptype &tinfo = time_map[tl] ;
    maptype::const_iterator ii ;
    for(ii=tinfo.begin();ii!=tinfo.end();++ii) {
      if(ii->second.size() < 2)
        continue ;
      std::list<variable>::const_iterator jj ;
      bool overlap = false ;
      variableSet vtouch ;
      for(jj=ii->second.begin();jj!=ii->second.end();++jj) {
        variableSet aliases = get_aliases(*jj) ;
        variableSet as = aliases ;
        variableSet::const_iterator vi ;
        
        for(vi=aliases.begin();vi!=aliases.end();++vi) 
          as += get_synonyms(*vi) ;
        if((as & vtouch) != EMPTY)
          overlap = true ;
        vtouch += as ;
      }
      if(!overlap) {
        //        cout << "memory variable " << ii->first << "{" << tl << "}"
        //             << " has a history of " << ii->second.size() << endl ;
        ti->rotate_lists.push_back(ii->second) ;
      } else {
        if(ii->second.size() !=2) {
          cerr << "unable to have history on variables aliased in time"
               << endl
               << "error occured on variable " << ii->first << "{" << tl << "}"
               << endl ;
          exit(-1) ;
        }
      }
    }
    return ti ;
  }
  
  void fact_db::initialize_time(fact_db::time_infoP ti) {
    *(ti->time_var) = 0 ;
  }
  
  void fact_db::advance_time(fact_db::time_infoP ti) {
    (*(ti->time_var))++ ;
    list<list<variable> >::const_iterator ii ;
    list<variable>::const_iterator jj ;

    for(ii=ti->rotate_lists.begin();ii!=ti->rotate_lists.end();++ii) {
      jj=ii->begin() ;
      storeRepP cp = get_fact_data(*jj).data_rep->getRep() ;
      ++jj ;
      if(jj != ii->end()) {
        for(;jj!=ii->end();++jj) {
          fact_data &fd = get_fact_data(*jj) ;
          storeRepP tmp = fd.data_rep->getRep() ;
          *(fd.data_rep) = cp ;
          cp = tmp ;
        }
      }
      *(get_fact_data(ii->front()).data_rep) = cp ;
    }
  }

  void fact_db::close_time(fact_db::time_infoP ti) {
    *(ti->time_var) = 0 ;
  }
  
  void fact_db::set_existential_info(variable v, rule f, entitySet x) {
    fact_info &finfo = get_fact_info(v) ;
    if((finfo.existence & x) != EMPTY) {
      ruleSet conflicts ;
      map<rule,existential_info>::iterator mi ;
      for(mi=finfo.exist_map.begin();mi!=finfo.exist_map.end();++mi) {
        if(mi->first == f) {
          continue ;
        }
        if((mi->second.exists & x) != EMPTY) {
          const vector<string> &p1 = v.get_info().priority ;
          const vector<string> &p2 = mi->second.v.get_info().priority ;
          for(int i=0;i<min(p1.size(),p2.size());++i)
            if(p1[i] != p2[i])
              conflicts += mi->first ;
          if(p1.size() == p2.size())
            conflicts += mi->first ;
          else if(p1.size() > p2.size()) {
            mi->second.exists -= x ;
            //cerr << f << " has priority over " << mi->first << endl ;
          } else {
            x -= mi->second.exists ;
            //cerr << mi->first << " has priority over " << f << endl ;
          }
        }
      }
      if(conflicts != EMPTY && v.get_info().name != string("OUTPUT")) {
        cerr << "rule " << f << " conflicts with " << conflicts << endl ;
        exit(-1) ;
      }
    }
    existential_info &einfo = finfo.exist_map[f] ;
    einfo.v = v ;
    einfo.exists += x ;
    finfo.existence += x ;
  }
  
  entitySet fact_db::variable_existence(variable v) {
    if(v.get_info().tvar)  // if time variable, variable exists everywhere
      return ~EMPTY ;
    return get_fact_info(v).existence ;
  }

  ruleSet fact_db::get_existential_rules(variable v) {
    std::map<rule,existential_info>::const_iterator mi ;
    fact_info &finfo = get_fact_info(v) ;
    ruleSet rules ;
    for(mi=finfo.exist_map.begin();mi!=finfo.exist_map.end();++mi)
      rules += mi->first ;
    return rules ;
  }
  
  void fact_db::variable_request(variable v, entitySet e) {
    if(v.get_info().tvar)
      return ;
    get_fact_info(v).requested += e ;
  }

  entitySet fact_db::get_variable_request(rule f, variable v) {
    fact_info &finfo = get_fact_info(v) ;
    map<rule,existential_info>::iterator mi = finfo.exist_map.find(f) ;
    if(mi == finfo.exist_map.end())
      return EMPTY ;
    return mi->second.exists & finfo.requested ;
  }
  
  bool fact_db::is_a_Map(variable v) {
    return get_fact_data(v).ismap ;
  }

  entitySet fact_db::image(variable v, entitySet e) {
    fact_data &fdata = get_fact_data(v) ;
    if(!fdata.ismap)
      return EMPTY ;
    map<entitySet,entitySet>::const_iterator ii = fdata.imageMap.find(e) ;
    if(ii != fdata.imageMap.end())
      return ii->second ;
    else
      return fdata.imageMap[e] = get_fact_data(v).minfo->image(e) ;
  }

  pair<entitySet,entitySet> fact_db::preimage(variable v, entitySet e) {
    fact_data &fdata = get_fact_data(v) ;
    if(!fdata.ismap)
      return make_pair(EMPTY,EMPTY) ;
    map<entitySet,pair<entitySet,entitySet> >::const_iterator ii ;
    ii = fdata.preimageMap.find(e) ;
    if(ii != fdata.preimageMap.end())
      return ii->second ;
    else
      return fdata.preimageMap[e] = get_fact_data(v).minfo->preimage(e) ;
  }
  
  void fact_db::printSummary(ostream &s) const {
    vmap_type::const_iterator vmi ;
    for(vmi=vmap.begin();vmi!=vmap.end();++vmi) {
      s << "--------------------------------------------------------------"
        << endl ;
      s << "variable: " << vmi->first << ", installed = " <<
        vmi->second.fact_installed << endl ;
      s << "synonyms: " << vmi->second.synonyms << endl ;
      s << "existence: " << vmi->second.existence << endl ;
      s << "requested: " << vmi->second.requested << endl ;
      const variableSet &aliases = fact_infov[vmi->second.fact_info_ref].aliases ;
      s << "aliases: " << aliases << endl ;

    }

  }

  ostream &fact_db::write(ostream &s) const {
    vmap_type::const_iterator vmi ;
    for(vmi=vmap.begin();vmi!=vmap.end();++vmi) {
      variable v=vmi->first;
      const fact_data &fd = fact_infov[vmi->second.fact_info_ref] ;
      storeRepP store_Rep = storeRepP(fd.data_rep) ;
      entitySet en=store_Rep->domain();
      std::string groupname = (v.get_info()).name;
      s << groupname << ":" ;
      store_Rep->Print(s);
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
    H5::H5File file(filename, H5F_ACC_TRUNC);

    vmap_type::const_iterator vmi ;
    for(vmi=vmap.begin();vmi!=vmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v);
      entitySet en=store_Rep->domain();
      std::string groupname = (v.get_info()).name;
      cout<<"Write "<<groupname<<" to HDF5 file "<<endl;
      //store_Rep->Print(cout);
      H5::Group group = file.createGroup("/"+groupname);
      (store_Rep->getRep())->writehdf5(group,en);
    }
  }

  void fact_db::read_hdf5(const char *filename){
    H5::H5File file(filename, H5F_ACC_RDONLY);

    vmap_type::const_iterator vmi ;
    for(vmi=vmap.begin();vmi!=vmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP store_Rep = get_variable(v)->getRep();
      std::string groupname = (v.get_info()).name;
      cout<<"Read "<<groupname<<" from HDF5 file "<<endl;
      if(H5Gopen(file.getId(),groupname.c_str())>0){
        H5::Group group = file.openGroup("/"+groupname);
        entitySet dom = store_Rep->domain() ;
        store_Rep->readhdf5(group, dom);
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
      storeRepP p = facts.get_variable(*vi) ;
      facts.update_fact(*vi,p->remap(remap)) ;
    }
  }
}


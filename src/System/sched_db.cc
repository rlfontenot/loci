#include <sched_db.h>
using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;
using std::make_pair ;

#include <Tools/parse.h>

//extern "C" {
//#include <hdf5.h>
//}

#include <Tools/stream.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;

  extern ofstream debugout ;
  
  sched_db::sched_db(fact_db &facts) {
    variableSet tmp_all_vars = facts.get_typed_variables() ;
    for(int i = 0; i < tmp_all_vars.size(); i++)
      sched_infov.push_back(sched_data()) ;
    int sched_alloc = 0 ;
    for(variableSet::const_iterator vi = tmp_all_vars.begin(); vi != tmp_all_vars.end(); ++vi) {
      sched_info si ; 
      si.sched_info_ref = sched_alloc++ ;
      storeRepP rp = facts.get_variable(*vi) ;
      si.existence = rp->domain() ;
      si.fact_installed = si.existence ;
      si.synonyms  += *vi ;
      vmap[*vi] = si ;
      variable v = variable(*vi) ;
      storeRepP st = rp->getRep() ;
      sched_data sd(v, st) ;
      sd.aliases += *vi ;
      sched_infov[si.sched_info_ref] = sd ;
    }

    //    synonyms = facts.synonyms ;
    //    free_set = facts.free_set ;

    all_vars = tmp_all_vars ;
  }
  
  sched_db::sched_db() {}
  sched_db::~sched_db() {}
  
  void sched_db::install_sched_info(variable v, sched_info info) {
    vmap_type::iterator mi = vmap.find(v) ;
    if(mi != vmap.end()) {
      cerr << "error:  reinstalling fact in database for variable "
	   << v << endl ;
      abort() ;
    }
    info.synonyms += v ;
    vmap[v] = info ;
    sched_infov[info.sched_info_ref].aliases += v ;
    all_vars += v ;
  }
  
  void sched_db::install_sched_data(variable v, sched_data data) {
    if(free_set == EMPTY) {
      sched_infov.push_back(data) ;
      install_sched_info(v,sched_info(sched_infov.size()-1)) ;
    } else {
      int val = *(free_set.begin()) ;
      free_set -= val ;
      sched_infov[val] = data ;
      install_sched_info(v,sched_info(val)) ;
    }
  }
  sched_db::sched_info &sched_db::get_sched_info(variable v) {
    vmap_type::iterator mi = vmap.find(remove_synonym(v)) ;
    if(mi == vmap.end()) 
      return get_sched_info(variable("EMPTY")) ;
    return mi->second ;
  }
  
  void sched_db::variable_is_fact_at(variable v,entitySet s, fact_db &facts) {
    sched_info &fi = get_sched_info(v) ;
    fi.fact_installed += s ;
    fi.existence += s ;
    variableSet aliases = sched_infov[fi.sched_info_ref].aliases ;
    aliases -= v ;
    for(variableSet::const_iterator vi=aliases.begin();vi!=aliases.end();++vi) 
      get_sched_info(v).fact_installed -= s ;
  }
  
  void sched_db::alias_variable(variable v, variable alias, fact_db &facts) {

    facts.synonym_variable(v,alias) ;

    if(all_vars.inSet(v)) {
      if(all_vars.inSet(alias)) {
        cerr << "alias already in fact_db!" << endl ;
        cerr << "error found in alias_variable("<<v<<","<< alias<<")" << endl ;
        abort() ;
      }
      
      int ref = get_sched_info(v).sched_info_ref ;
      sched_infov[ref].aliases += alias ;
      install_sched_info(alias,sched_info(ref)) ;
      
    } else if(all_vars.inSet(alias)) {
      alias_variable(alias,v, facts) ;
    } else {
      cerr << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
      abort() ;
    }
    
  }

  void sched_db::synonym_variable(variable v, variable synonym, fact_db &facts) {
    facts.synonym_variable(v, synonym) ;
    v = remove_synonym(v) ;
    vmap_type::iterator vmi ;
    if((vmi = vmap.find(synonym)) != vmap.end()) {
      cerr << "synonym already in fact_db!" << endl ;
      cerr << "error found in synonym_variable("<<v<<","<<synonym<<")"<<endl;
      abort() ;
    }

    sched_info &vfinfo = get_sched_info(v) ;
    
    vfinfo.synonyms += synonym ;
    synonyms[synonym] = v ;
    all_vars += synonym ;
  }
  
  void sched_db::set_existential_info(variable v, rule f, entitySet x) {
    sched_info &finfo = get_sched_info(v) ;
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
          for(int i=0,j=p1.size()-1,k=p2.size()-1;
              i<min(int(p1.size()),int(p2.size()));
              ++i,--j,--k)
            if(p1[j] != p2[k]) {
              conflicts += mi->first ;
              cerr << "adding to conflicts because " << p1[i]
                   << "!=" << p2[i]<<endl ;
            }
          if(p1.size() == p2.size()) {
            conflicts += mi->first ;
          } else if(p1.size() > p2.size()) {
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
        cerr << "conflicting entities are " << (finfo.existence & x) << endl ;
        //        debugger_() ;
        //        exit(-1) ;
      }
    }
    existential_info &einfo = finfo.exist_map[f] ;
    einfo.v = v ;
    einfo.exists += x ;
    finfo.existence += x ;
  }
  
  entitySet sched_db::variable_existence(variable v) {
    if(v.get_info().tvar)  // if time variable, variable exists everywhere
      return ~EMPTY ;
    return get_sched_info(v).existence ;
  }
  
  ruleSet sched_db::get_existential_rules(variable v) {
    std::map<rule,existential_info>::const_iterator mi ;
    sched_info &finfo = get_sched_info(v) ;
    ruleSet rules ;
    for(mi=finfo.exist_map.begin();mi!=finfo.exist_map.end();++mi)
      rules += mi->first ;
    return rules ;
  }
  
  void sched_db::variable_request(variable v, entitySet e) {
    if(v.get_info().tvar)
      return ;
    get_sched_info(v).requested += e ;
  }
  
  entitySet sched_db::get_variable_request(rule f, variable v) {
    sched_info &finfo = get_sched_info(v) ;
    map<rule,existential_info>::iterator mi = finfo.exist_map.find(f) ;
    if(mi == finfo.exist_map.end())
      return EMPTY ;
    return mi->second.exists & finfo.requested ;
  }

  void sched_db::set_variable_type(variable v, storeRepP st, fact_db &facts) {
    facts.create_fact(v, st) ;
    if(!all_vars.inSet(v)) {
      install_sched_data(v, sched_data(v, st)) ;
    }
  }

  entitySet sched_db::image(variable v, entitySet e) {
    sched_data &fdata = get_sched_data(v) ;
    if(!fdata.ismap)
      return EMPTY ;
    map<entitySet,entitySet>::const_iterator ii = fdata.imageMap.find(e) ;
    if(ii != fdata.imageMap.end())
      return ii->second ;
    else
      return fdata.imageMap[e] = get_sched_data(v).minfo->image(e) ;
  }
  
  pair<entitySet,entitySet> sched_db::preimage(variable v, entitySet e) {
    sched_data &fdata = get_sched_data(v) ;
    if(!fdata.ismap)
      return make_pair(EMPTY,EMPTY) ;
    map<entitySet,pair<entitySet,entitySet> >::const_iterator ii ;
    ii = fdata.preimageMap.find(e) ;
    if(ii != fdata.preimageMap.end())
      return ii->second ;
    else
      return fdata.preimageMap[e] = get_sched_data(v).minfo->preimage(e) ;
  }

  std::ostream &sched_db::print_summary(fact_db &facts, std::ostream &s) {
    s << "Summary of Existential deduction:" << endl ;
    std::map<variable,sched_info>::const_iterator mi ;
    for(mi=vmap.begin();mi!=vmap.end();++mi) {
      storeRepP sp = facts.get_variable(mi->first) ;
      if(sp->RepType() == MAP)
	s << " Container = MAP " << endl ;
      else if(sp->RepType() == STORE)
	s << " Container = STORE " << endl ;
      else if(sp->RepType() == PARAMETER)
	s << " Container = PARAMETER " << endl ;
      else if(sp->RepType() == CONSTRAINT)
	s << " Container = CONSTRAINT " << endl ;
      
      double size = 0 ;
      size = sp->pack_size(mi->second.requested) / 1000000 ;
      s << mi->first << " " <<mi->second.synonyms << " "<< mi->second.existence
        << " request= " << mi->second.requested << "  size requested = " << size << " MB " << endl ;
    }
    return s ;
  }
  
}

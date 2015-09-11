//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <typeinfo>
#include <Tools/except.h>

#include <gfact_db.h>
//#include <rule.h>
//#include <constraint.h>
#include <Tools/debugger.h>
#include <gstore.h>
#include "dist_tools.h"

#include <Loci_Datatypes.h>

using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;


using std::istream ;
using std::ostream ;
using std::ostringstream ;
using std::endl ;
using std::ios ;

#include <Tools/parse.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;
    
  gfact_db::gfact_db() {
    //set default key manager as distKeyManager 
    key_manager = CPTR<gKeyManager>(new distKeyManager());
  }

  gfact_db::~gfact_db() {}

  //this function is not tested
  void
  gfact_db::copy_all_from(const gfact_db& f) {
    synonyms = f.synonyms ;
    fmap = f.fmap ;
    nspace_vec = f.nspace_vec ;

  


    // deep copy the key manager
    if(f.key_manager != 0) {
      key_manager = f.key_manager->clone() ;
    }
    
  }

  variable gfact_db::remove_synonym(variable v) const {
    std::map<variable,variable>::const_iterator mi ;
    while((mi=synonyms.find(v)) != synonyms.end())
      v = mi->second ;
    return v ;
  }
  
  variable gfact_db::get_real_var(variable v) const {
    //  variable v = add_namespace(variable(vname)) ;
    variable real_v = v;
    if(synonyms.find(v) != synonyms.end()) {
      real_v = remove_synonym(v) ;
    }
    return real_v;
  }
    

  
  void gfact_db::synonym_variable(variable v, variable synonym) {
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
    std::map<variable,gStoreRepP>::iterator vmi, vmj ;
    if((vmi = fmap.find(v)) == fmap.end()) {
      cerr << "WARNING: synonym_variable("<<v<<","<<synonym<<")"<<endl ;
      cerr << "WARNING: type not known for target of synonym, ignoring operation" << endl ;
      //      abort() ;
      return ;
    }
    // If the synonym already points to a different variable instance,
    // remove it
    if((vmj = fmap.find(s)) != fmap.end()) {
      gStoreRepP &finfo = vmj->second ;
      if((finfo->domain() != genIntervalSet<gEntity>::EMPTY &&
          finfo->RepType() != GPARAMETER)) {
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
  
  
 

  variable gfact_db::add_namespace(variable v) const {
    variable tmp_v ;
    if((nspace_vec.size() > 0) && (v.get_info().namespac.size() == 0)) {
      tmp_v = v ;
      for(size_t i = 0; i < nspace_vec.size(); ++i)
        tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    } else
      tmp_v = v ;
    return tmp_v ;
  }
  // this is the basic method that creats a fact
  // in the gfact_db. It is served as the basis
  // for create_fact methods
  //this method adds st to fmap with duplication check
  //this method will not process keyspace info
  void gfact_db::create_pure_fact(const variable& v, gStoreRepP st) {
    //get real variable
    variable real_v = get_real_var(v);
    std::map<variable, gStoreRepP>::iterator mi = fmap.find(real_v) ;
    if(mi==fmap.end()) {
      fmap[real_v]= st ;
    } else {
      if(typeid(st->getRep()) != typeid(mi->second->getRep())) {
        cerr << "set_variable_type() method of gfact_db changing type for variable " << v << endl ;
      }
      mi->second = st ;
    }
    return ;
  }

  void gfact_db::set_variable_domain_space(const variable& v, gStoreRepP st, gKeySpaceP space){
    if(space != gKeySpaceP(0)){
      //get real variable   
      variable real_v = get_real_var(v) ;
      st->set_domain_space(&(*space));
      space->add_out_var(real_v);
    }
    return;
  }
  void gfact_db::set_variable_image_space(const variable& v, gMapRepP st, gKeySpaceP space){
    if(space != 0){
      //get real variable   
      variable real_v = get_real_var(v) ;
      st->set_image_space(&(*space));
      space->add_in_var(real_v); 
      return;
    }
  }
  void gfact_db::remove_variable_domain_space(const variable& v){
     
    gStoreRepP st = get_variable(v);
    gKeySpaceP space = st->get_domain_space();
    if(space != 0 ){
      variable real_v = get_real_var(v);
      space->remove_out_var(real_v);
      st->set_domain_space(0);
    }
  }
  
  void gfact_db::remove_variable_image_space(const variable& v){
    //get real varaible 
    gStoreRepP st = get_variable(v);
    gKeySpaceP space = gMapRepP(st)->get_image_space();
    if(space != 0){
      variable real_v = get_real_var(v);
      space->remove_in_var(real_v);
      gMapRepP(st)->set_image_space(0);
    }
  }
    

  
  void gfact_db::create_fact(const variable& v, gStoreRepP st,
                             gKeySpaceP domain_space,
                             gKeySpaceP image_space) {
    variable tmp_v =add_namespace(v) ;
    create_pure_fact(tmp_v,st) ;
    if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
  }
  void gfact_db::create_fact(const std::string& vname, gStoreRepP st,
                             gKeySpaceP domain_space ,
                             gKeySpaceP image_space ) {
    variable v = variable(vname) ;
    create_fact(v,st, domain_space, image_space) ;
  }
  
  void gfact_db::create_fact(const variable& v, gstore_instance &si,
                             gKeySpaceP domain_space,
                             gKeySpaceP image_space) {
    variable tmp_v = add_namespace(v) ;
    gStoreRepP st = si.Rep(); 
    create_pure_fact(tmp_v,st) ;
    if(domain_space != 0)set_variable_domain_space(v, st, domain_space);
    if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    si.setRep(get_variable(v)) ;
  }

  void gfact_db::create_fact(const std::string& vname, gstore_instance &si,
                             gKeySpaceP domain_space ,
                             gKeySpaceP image_space ) {
    variable v = variable(vname) ;
    create_fact(v,si, domain_space, image_space) ;
  }


  //update fmap,
  //if new keyspaces are provided, remove v from old spaces, and connect v with  new ones
  //otherwise, assume st already has the domain space and image space set up
  void gfact_db::update_fact(variable v, gStoreRepP st,
                             gKeySpaceP domain_space ,
                             gKeySpaceP image_space ) {
    //get real variable
    variable tmp_v = add_namespace(v) ;
      
    //tmp_v should not have synonyms
    warn(synonyms.find(tmp_v) != synonyms.end()) ;

    //update fmap and keyspaces
    std::map<variable, gStoreRepP>::iterator mi = fmap.find(tmp_v) ;
     
    //tmp_v should be in fmap
    if(mi != fmap.end()) {
      if(domain_space != 0){
        remove_variable_domain_space(v);
      }
      if(image_space != 0){
        remove_variable_image_space(v);
      }
         
      mi->second = st ;

      if(domain_space != 0)set_variable_domain_space(v, st, domain_space);
      if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    } else{
      cerr << "warning: update_fact: fact does not exist for variable " << tmp_v
           << endl ;
    }
  }

  
  /*! remove from synonym and fmap */
  //also remove from keyspaces
  void gfact_db::remove_variable(variable v) {
    std::map<variable, variable>::iterator si ;
    std::map<variable, gStoreRepP>::iterator mi ;
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
      //before erase the variable, remove it from keyspaces
      remove_variable_domain_space(v);
      remove_variable_image_space(v);
      // Now erase the variable
      fmap.erase(mi) ;
     
    }
  }
  
  variableSet gfact_db::get_typed_variables() const {
    std::map<variable, gStoreRepP>::const_iterator mi ;
    std::map<variable, variable>::const_iterator si ;
    variableSet all_vars ;
    for(mi=fmap.begin();mi!=fmap.end();++mi)
      all_vars += mi->first ;
    
    for(si=synonyms.begin();si!=synonyms.end();++si)
      all_vars += si->first ;
    return all_vars ;
  }


 
  
  gStoreRepP gfact_db::get_variable(variable v) {
    variable tmp_v = add_namespace(v) ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, gStoreRepP>::iterator mi =
      fmap.find(remove_synonym(tmp_v)) ;
    if(mi == fmap.end()) {
      if(Loci::MPI_rank == 0)
      	cerr << " returning null  gStoreRep for variable " << tmp_v<< endl ;
      return gStoreRepP(0) ;
    }
    else
      return mi->second ;
  }

  
 
}

  

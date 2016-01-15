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

#ifndef GFACT_DB_H
#define GFACT_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/except.h>
#include <rule.h>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <gstore_rep.h>
#include <variable.h>
#include <gmap_rep.h>
#include <gmap.h>
#include <gkey_manager.h>
#include <gkeyspace.h>
#include <distribute_long.h>
#include <gparameter.h>
#include <Tools/simple_partition_long.h>
using std::string;


namespace Loci {
  class rule_db ;

  /*
    Developer's notes: To repalce fact_db by gfact_db, first rename everything in gfact_db related to gcontainers so that
    the member data and member method related to tradititional containers and other info will remains the same as fact_db. 
  */
  class gfact_db {
  private:
    // key manager
    gKeyManagerP gkey_manager ;
    
    //map each variable to its gcontainer storage
    std::map<variable,gStoreRepP> gfmap ;
    
    //map each variable to its traditional container storage
    //this map is not constructed by create_gfact method
    //instead, when variables freeze, they are removed from gfmap
    // and added to fmap
    std::map<variable,storeRepP> fmap ; 
    
    //map the variables to their types
    //tmap and gtmap are used for optional rules and default rules
    //the rules will specify the type of variables
    //and .vars file will provide the data
    //The target variables of these rules are gParams in UniverseSpace,
    //not over subset of a specific keyspace
    std::map<variable,storeRepP> tmap ;
    


    // support for multiple queries and experimental
    // extensions to the fact_db to distinguish
    // extensional facts and intensional facts
    variableSet extensional_facts ;

    
    std::vector<std::string> nspace_vec ;//allow namespaces in front of variable name
    /*! all variables that point to the same gStoreRepP, the second of the pair at the end of chain
      is the variable suppose to appear in gfmap*/
    std::map<variable,variable> synonyms ; 
    
    
    // //don't know what its for yet
    // std::set<std::vector<variableSet> > intensive_output_maps;
  private:
    // a copy function
    void copy_all_from(const gfact_db& f) ;
    variable remove_synonym(variable v) const;

    void set_variable_domain_space(const variable& v, gStoreRepP st, gKeySpaceP space);
    void set_variable_image_space(const variable& v, gMapRepP st, gKeySpaceP space);
    void remove_variable_domain_space(const variable& v);
    void remove_variable_image_space(const variable& v);
    
    // this is the basic method that creats a fact
    // in the gfact_db. It is served as the basis
    // for create_gfact methods
    //this method adds st to gfmap with duplication check
    //this method will not process keyspace info
    void create_pure_gfact(const variable& v, gStoreRepP st) ;
  public:
    //constructor
    gfact_db() ;
    // copy constructor
    gfact_db(const gfact_db &f) {
      copy_all_from(f) ;
    }
    // the assignment operator
    gfact_db &operator=(const gfact_db &f) {
      if(&f != this) {
        copy_all_from(f) ;
      }
      return *this ;
    }
    //destructor
    ~gfact_db() ;

    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname, storeRepP st)
    { set_variable_type(variable(vname),st) ;}
    void set_variable_type(variable v, store_instance &si)
    { set_variable_type(v,si.Rep()) ; }
    void set_variable_type(std::string vname, store_instance &si)
    { set_variable_type(variable(vname),si) ; }
    
    storeRepP get_variable_type(variable v) const ;
    storeRepP get_variable_type(std::string vname) const
    { return get_variable_type(variable(vname)) ;}
    
    //it is safe to be public?
    void set_key_manager( gKeyManagerP km){gkey_manager=km;}
    
    //real_var is the var used in paces such as gfmap, in_vars and out_vars of keyspaces, etc.
    //this method remove remove_synonym
    variable get_real_var(variable  v) const; 
    
    variable add_namespace(variable v) const ;
    //copy all variables in gfmap to traditional containers
    void copy_facts(fact_db& facts) const ;
    
  public:

    // create_gfact now defaults to create an extensional fact
    // as this is the primary interface for users of Loci
    void create_gfact(const variable& v, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0) ;
    
    void create_gfact(const std::string& vname, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0); 
    
    void create_gfact(const variable& v, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0) ;
    
    void create_gfact(const std::string& vname, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0);

    //update gfmap,
    //if new keyspaces are provided, remove v from old spaces, and connect v with  new ones
    //otherwise, connect the old spaces to st
    void update_gfact(variable v, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0);
    
    void update_gfact(std::string vname, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_gfact(variable(vname),st, domain_space, image_space) ;}
    void update_gfact(variable v, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_gfact(v,si.Rep(), domain_space, image_space) ; si.setRep(get_gvariable(v)) ; }
    void update_gfact(std::string vname, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_gfact(variable(vname),si, domain_space, image_space) ; }

    void replace_gfact(variable v, gStoreRepP st,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0) 
    { remove_gvariable(v) ;
      create_gfact(v,st, domain_space, image_space) ; }
    void replace_gfact(std::string vname, gStoreRepP st,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0)
    { replace_gfact(variable(vname),st, domain_space, image_space) ;}
    void replace_gfact(variable v, gstore_instance &si,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0)
    { replace_gfact(v,si.Rep(), domain_space, image_space) ;
      si.setRep(get_gvariable(v)) ; }
    void replace_gfact(std::string vname, gstore_instance &si,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0 )
    { replace_gfact(variable(vname),si, domain_space, image_space) ; }

    gStoreRepP get_gfact(variable &v) { return get_gvariable(v); }
    gStoreRepP get_gfact(std::string vname)
    { return get_gvariable(variable(vname)) ; }

    gStoreRepP get_gvariable(variable v) ;
    gStoreRepP get_gvariable(std::string vname)
    { return get_gvariable(variable(vname)) ; }

    /////////////////////////////////////////////////////////
    // support methods for extensional & intensional facts //
    /////////////////////////////////////////////////////////
    variableSet get_extensional_facts() const {
      return extensional_facts ;
    }
    variableSet get_intensional_facts() const {
      return variableSet(get_typed_variables()-extensional_facts) ;
    }

    // // we no longer provide create_extensional_fact methods with explicit name to
    // // create extentional facts, they are just as the same as
    // // the default create_gfact methods
    // void create_extensional_fact(const variable& v, gStoreRepP st) {
    //   create_gfact(v,st) ;
    // }
    // void create_extensional_fact(const std::string& vname, gStoreRepP st) {
    //   create_gfact(vname,st) ;
    // }
    // void create_extensional_fact(const variable& v, gstore_instance &si) {
    //   create_gfact(v,si) ;
    // }
    // void create_extensional_fact(const std::string& vname,
    //                              gstore_instance &si) {
    //   create_gfact(vname,si) ;
    // }
    // this method will convert all intensional facts (if any) in
    // the fact database to extensional facts
    void make_all_extensional() {
      variableSet intensional_facts = get_intensional_facts() ;
      extensional_facts += intensional_facts ;
    }

    // // we have the corresponding intensional facts creation, which create fact in gfmap
    // //Don't know if these methods are needed or not
    // void create_intensional_gfact(const variable& v, gStoreRepP st,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   variable v_tmp = add_namespace(v) ;
    //   create_pure_gfact(v_tmp,st) ;
    //   if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    //   if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    // }
    // void create_intensional_gfact(const std::string& vname, gStoreRepP st,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   create_intensional_fact(variable(vname), st, domain_space, image_space);
      
    // }
    // void create_intensional_gfact(const variable& v, gstore_instance &si,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   variable v_tmp = add_namespace(v) ;
    //   create_pure_gfact(v_tmp,si.Rep()) ;
    //   gStoreRepP st = si.Rep(); 
    //   if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    //   if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    //   si.setRep(get_gvariable(v_tmp)) ;
    // }
    // void create_intensional_gfact(const std::string& vname,
    //                              gstore_instance &si,
    //                               gKeySpaceP domain_space = 0,
    //                               gKeySpaceP image_space = 0)
    // {
    //   create_intensional_fact(variable(vname), si, domain_space, image_space);
    // }





    
    // // we have the corresponding intensional facts creation, which create fact in fmap
    // void create_intensional_fact(const variable& v, gStoreRepP st,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   variable v_tmp = add_namespace(v) ;
    //   create_pure_fact(v_tmp,st) ;
    //   if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    //   if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    // }
    // void create_intensional_fact(const std::string& vname, gStoreRepP st,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   create_intensional_fact(variable(vname), st, domain_space, image_space);
      
    // }
    // void create_intensional_fact(const variable& v, gstore_instance &si,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   variable v_tmp = add_namespace(v) ;
    //   create_pure_fact(v_tmp,si.Rep()) ;
    //   gStoreRepP st = si.Rep(); 
    //   if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    //   if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    //   si.setRep(get_gvariable(v_tmp)) ;
    // }
    // void create_intensional_fact(const std::string& vname,
    //                              gstore_instance &si,
    //                              gKeySpaceP domain_space = 0,
    //                              gKeySpaceP image_space = 0)
    // {
    //   create_intensional_fact(variable(vname), si, domain_space, image_space);
    // }

    
    // // this method erases all intensional facts, assume it is traditional facts for now
    // void erase_intensional_facts() {
    //   variableSet intensional_facts = get_intensional_facts() ;
    //   for(variableSet::const_iterator vi=intensional_facts.begin();
    //       vi!=intensional_facts.end();++vi)
    //     remove_variable(*vi) ;
    // }
    // // this method will convert an intensional fact to
    // // a extensional fact
    // void make_extensional_fact(const variable& v) ;
    // void make_extensional_fact(const std::string& vname) {
    //   make_extensional_fact(variable(vname)) ;
    // }
    // // this method will convert a extensional fact
    // // to an intensional one
    // void make_intensional_fact(const variable& v) ;
    // void make_intensional_fact(const std::string& vname) {
    //   make_intensional_fact(variable(vname)) ;
    // }
    
    //the following two methods return the storeRepP in fmap
    storeRepP get_variable(variable v) ;
    storeRepP get_variable(std::string vname)
    { return get_variable(variable(vname)) ; }

    void remove_gvariable(variable v) ;
    
    void synonym_variable(variable v, variable synonym) ;

   
    
    void set_namespace(std::string name_space){
      nspace_vec.insert(nspace_vec.begin(), 1, name_space) ; 
    }
    void remove_namespace() {
      nspace_vec.pop_back() ;
    }  
    void unset_namespace() {
      nspace_vec.clear() ;
    }
    
    //this method returns all variables in gfmap and their synonyms
    //has nothing to do with tmap
    variableSet get_typed_variables() const ;

    std::istream& read_vars(std::istream& s, const rule_db& rdb);
    void read_vars(std::string filename, const rule_db &rdb) {
      std::ifstream ifile(filename.c_str(),std::ios::in) ;
      if(ifile.fail()) {
        std::string error = std::string("Can't open '") +filename + "'" ;
        throw(StringError(error)) ;
      }
      read_vars(ifile,rdb) ;
      ifile.close() ;
    }
    gKeyManagerP get_key_manager() const {return gkey_manager ;}
  
   
  };
}
    

#endif

    

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
  class gfact_db {
  private:
    // key manager
    gKeyManagerP key_manager ;
   
   
    std::map<variable,gStoreRepP> fmap ; //map each variable to its storage
    std::vector<std::string> nspace_vec ;//allow namespaces in front of variable name
    /*! all variables that point to the same gStoreRepP, the second of the pair at the end of chain
      is the variable suppose to appear in fmap*/
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
    // for create_fact methods
    //this method adds st to fmap with duplication check
    //this method will not process keyspace info
    void create_pure_fact(const variable& v, gStoreRepP st) ;
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
    
    //it is safe to be public?
    void set_key_manager( gKeyManagerP km){key_manager=km;}
    
    //real_var is the var used in paces such as fmap, in_vars and out_vars of keyspaces, etc.
    //this method remove remove_synonym
    variable get_real_var(variable  v) const; 
    
    variable add_namespace(variable v) const ;
    //copy all variables in fmap to traditional containers
    void copy_facts(fact_db& facts) const ;
    
  public:
       
    void create_fact(const variable& v, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0) ;
    
    void create_fact(const std::string& vname, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0); 
    
    void create_fact(const variable& v, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0) ;
    
    void create_fact(const std::string& vname, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0);

    //update fmap,
    //if new keyspaces are provided, remove v from old spaces, and connect v with  new ones
    //otherwise, connect the old spaces to st
    void update_fact(variable v, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0);
    
    void update_fact(std::string vname, gStoreRepP st,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_fact(variable(vname),st, domain_space, image_space) ;}
    void update_fact(variable v, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_fact(v,si.Rep(), domain_space, image_space) ; si.setRep(get_variable(v)) ; }
    void update_fact(std::string vname, gstore_instance &si,
                     gKeySpaceP domain_space = 0,
                     gKeySpaceP image_space = 0)
    { update_fact(variable(vname),si, domain_space, image_space) ; }

    void replace_fact(variable v, gStoreRepP st,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0) 
    { remove_variable(v) ;
      create_fact(v,st, domain_space, image_space) ; }
    void replace_fact(std::string vname, gStoreRepP st,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0)
    { replace_fact(variable(vname),st, domain_space, image_space) ;}
    void replace_fact(variable v, gstore_instance &si,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0)
    { replace_fact(v,si.Rep(), domain_space, image_space) ;
      si.setRep(get_variable(v)) ; }
    void replace_fact(std::string vname, gstore_instance &si,
                      gKeySpaceP domain_space = 0,
                      gKeySpaceP image_space = 0 )
    { replace_fact(variable(vname),si, domain_space, image_space) ; }

    gStoreRepP get_fact(variable &v) { return get_variable(v); }
    gStoreRepP get_fact(std::string vname)
    { return get_variable(variable(vname)) ; }

    gStoreRepP get_variable(variable v) ;
    gStoreRepP get_variable(std::string vname)
    { return get_variable(variable(vname)) ; }
    
    void remove_variable(variable v) ;
    
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
   
    variableSet get_typed_variables() const ;
    
    gKeyManagerP get_key_manager() const {return key_manager ;}
  
   
  };
}
    

#endif

    

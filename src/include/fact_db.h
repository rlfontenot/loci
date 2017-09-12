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

#ifndef FACT_DB_H
#define FACT_DB_H

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
////////////////////////////////////
//the following is temporary, 
////////////////////////////////////
#include <distribute.h>
#include <store_rep.h>
#include <Map_rep.h>
#include <Map.h>
#include <DMap.h>
#include <key_manager.h>
#include <keyspace.h>

/////////////////////////////////////////////
using std::string;


namespace Loci {
  class rule_db ;

  /*
    Developer's notes: To repalce fact_db by fact_db, first rename everything in fact_db related to gcontainers so that
    the member data and member method related to tradititional containers and other info will remains the same as fact_db.
    At present time, fact_db includes everything in fact_db, clean-up need to be done later. 
  */
  class fact_db {
  private:
    // key manager
    gKeyManagerP gkey_manager ;
    
    //map each variable to its gcontainer storage
    //when variables freeze, they are removed from gfmap
    // and added to fmap
    std::map<variable,gStoreRepP> gfmap ;


    struct fact_info {
      store_refP data_rep ;
    } ;
    //map each variable to its traditional container storage
    std::map<variable,fact_info> fmap ; 

    //map the variables to their types
    //tmap are used for optional rules and default rules
    //these rules will specify the type of variables
    //while .vars file will provide the data
    //The target variables of these rules are gParams in UniverseSpace,
    //not over a subset of a specific keyspace
    std::map<variable,storeRepP> tmap ;
    


    // support for multiple queries and experimental
    // extensions to the fact_db to distinguish
    // extensional facts and intensional facts
    variableSet extensional_facts ;

    
    std::vector<std::string> nspace_vec ;//allow namespaces in front of variable name
    /*! all variables that point to the same gStoreRepP, the second of the pair at the end of chain
      is the variable suppose to appear in gfmap*/
    std::map<variable,variable> synonyms ; 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //the following data is temporary, will clean up later
    /////////////////////////////////////////////////////////////////////////////////////////////////////
  public:
    struct distribute_info : public CPTR_type {
      struct dist_data {
        int proc ;
        entitySet entities ;
        int size ;
        dist_data(int p, entitySet &e) :proc(p),entities(e)
        { size = e.size() ; }
      } ;
      
      int myid ;
      int isDistributed ;
      Map l2g ; // local numbering to global numbering
      dMap g2l ; // global numbering to local numbering
      
      entitySet my_entities ;
      entitySet comp_entities;//local numbering of global_comp_entities
      
      std::vector<dist_data> copy ; // Data describing which processors own
      // the entities in the clone region
      std::vector<dist_data> xmit ; // Data describing which processors to
      // send entities that this processor owns
      int copy_total_size ;
      int xmit_total_size ;
      //      dMap remap ;
      dMap g2f ; // Global to file numbering
      distribute_info() {} ;
    }  ;
    std::vector<entitySet> init_ptn ;

    //Related to Duplication: global_comp_entities set defines
    //entities that can be computed on a processor
    entitySet global_comp_entities; 

    typedef CPTR<distribute_info> distribute_infoP ;
    distribute_infoP distributed_info ;
    
    std::set<std::vector<variableSet> > intensive_output_maps;

    // experimental inclusion of keyspace partitions
    std::map<std::string,KeySpaceP> keyspace ;
    // experimental inclusion of a key manager
    KeyManagerP key_manager ;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
   
    
    
  private:
    // a copy function
    void copy_all_from(const fact_db& f) ;
    

    void set_variable_domain_space(const variable& v, gStoreRepP st, gKeySpaceP space);
    void set_variable_image_space(const variable& v, gMapRepP st, gKeySpaceP space);
    void remove_variable_domain_space(const variable& v);
    void remove_variable_image_space(const variable& v);
    
    // this is the basic method that creats a fact
    // in the fact_db. It is served as the basis
    // for create_gfact methods
    //this method adds st to gfmap with duplication check
    //this method will not process keyspace info
    void create_pure_gfact(const variable& v, gStoreRepP st) ;
    variable add_namespace(variable v) const ;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //the following private methods are temporary
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // this is the basic method that creats a fact
    // in the fact_db. It is served as the basis
    // for create_fact & the create_intensional_fact methods
    void create_pure_fact(const variable& v, storeRepP st) ;
   
    std::pair<entitySet, entitySet> get_dist_alloc(int size) ;
    int maximum_allocated ;
    int minimum_allocated ;
    int dist_from_start ;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
  public:
    //constructor
    fact_db() ;
    // copy constructor
    fact_db(const fact_db &f) {
      copy_all_from(f) ;
    }
    // the assignment operator
    fact_db &operator=(const fact_db &f) {
      if(&f != this) {
        copy_all_from(f) ;
      }
      return *this ;
    }
    //destructor
    ~fact_db() ;

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
    
    //real_var is the var used in places such as gfmap, in_vars and out_vars of keyspaces, etc.
    variable get_real_var(variable  v) const; 
    
  
    //copy all variables in gfmap to traditional containers in fmap
    void copy_facts();
    
    variable remove_synonym(variable v) const;
    void remove_variable(variable v) ;
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
    
    //this method returns all variables and their synonyms
    //has nothing to do with tmap
    //should get from gfmap or fmap? assume it is fmap now
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
    gKeyManagerP get_gkey_manager() const {return gkey_manager ;}


    ////////////////////////////////////////////////////////////////////////////////// 
    //the following public methods are temporary, need reconsider later
    /////////////////////////////////////////////////////////////////////////////////

    void set_maximum_allocated(int i) ;
    void set_minimum_allocated(int i) ;
    int get_max_alloc() {
      return maximum_allocated ;
    }
    int get_min_alloc() {
      return minimum_allocated ;
    }

    // create_fact now defaults to create an extensional fact
    // as this is the primary interface for users of Loci
    //? Be careful, don't add namespace twice
    void create_fact(const variable& v, storeRepP st) {
      variable tmp_v =add_namespace(v) ;
      create_pure_fact(tmp_v,st) ;
      extensional_facts += tmp_v ;
    }
    void create_fact(const std::string& vname, storeRepP st) {
      variable v = variable(vname) ;
      create_fact(v,st) ;
    }
    void create_fact(const variable& v, store_instance &si) {
      variable tmp_v = add_namespace(v) ;
      create_pure_fact(tmp_v,si.Rep()) ;
      si.setRep(get_variable(v)) ;
      extensional_facts += tmp_v ;
    }
    void create_fact(const std::string& vname, store_instance &si) {
      variable v = add_namespace(variable(vname)) ;
      create_pure_fact(v,si.Rep()) ;
      si.setRep(get_variable(vname)) ;
      extensional_facts += v ;
    }
    


    
    void update_fact(variable v, storeRepP st);
    void update_fact(std::string vname, storeRepP st)
    { update_fact(variable(vname),st) ;}
    void update_fact(variable v, store_instance &si)
    { update_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void update_fact(std::string vname, store_instance &si)
    { update_fact(variable(vname),si) ; }

    void replace_fact(variable v, storeRepP st)
    { remove_variable(v) ; create_fact(v,st) ; }
    void replace_fact(std::string vname, storeRepP st)
    { replace_fact(variable(vname),st) ;}
    void replace_fact(variable v, store_instance &si)
    { replace_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void replace_fact(std::string vname, store_instance &si)
    { replace_fact(variable(vname),si) ; }
    storeRepP get_fact(variable &v) { return get_variable(v); }
    storeRepP get_fact(std::string vname)
    { return get_variable(variable(vname)) ; }
    
    storeRepP get_variable(variable v) ;
    storeRepP get_variable(std::string vname)
    { return get_variable(variable(vname)) ; }


    /////////////////////////////////////////////////////////
    // support methods for extensional & intensional facts //
    /////////////////////////////////////////////////////////
    variableSet get_extensional_facts() const {
      return extensional_facts ;
    }
    variableSet get_intensional_facts() const {
      return variableSet(get_typed_variables()-extensional_facts) ;
    }

    // we still provide these methods with explicit name to
    // create extentional facts, they are just as the same as
    // the default create_fact methods
    void create_extensional_fact(const variable& v, storeRepP st) {
      create_fact(v,st) ;
    }
    void create_extensional_fact(const std::string& vname, storeRepP st) {
      create_fact(vname,st) ;
    }
    void create_extensional_fact(const variable& v, store_instance &si) {
      create_fact(v,si) ;
    }
    void create_extensional_fact(const std::string& vname,
                                 store_instance &si) {
      create_fact(vname,si) ;
    }
    // this method will convert all intensional facts (if any) in
    // the fact database to extensional facts
    void make_all_extensional() {
      variableSet intensional_facts = get_intensional_facts() ;
      extensional_facts += intensional_facts ;
    }
    // and then we have the corresponding intensional facts creation
    void create_intensional_fact(const variable& v, storeRepP st) {
      variable v_tmp = add_namespace(v) ;
      create_pure_fact(v_tmp,st) ;
    }
    void create_intensional_fact(const std::string& vname, storeRepP st) {
      create_pure_fact(add_namespace(variable(vname)),st) ;
    }
    void create_intensional_fact(const variable& v, store_instance &si) {
      variable v_tmp = add_namespace(v) ;
      create_pure_fact(v_tmp,si.Rep()) ;
      si.setRep(get_variable(v_tmp)) ;
    }
    void create_intensional_fact(const std::string& vname,
                                 store_instance &si) {
      variable v = add_namespace(variable(vname)) ;
      create_pure_fact(v,si.Rep()) ;
      si.setRep(get_variable(v)) ;
    }
    // this method erases all intensional facts
    void erase_intensional_facts() {
      variableSet intensional_facts = get_intensional_facts() ;
      for(variableSet::const_iterator vi=intensional_facts.begin();
          vi!=intensional_facts.end();++vi)
        remove_variable(*vi) ;
    }
    // this method will convert an intensional fact to
    // a extensional fact
    void make_extensional_fact(const variable& v) ;
    void make_extensional_fact(const std::string& vname) {
      make_extensional_fact(variable(vname)) ;
    }
    // this method will convert a extensional fact
    // to an intensional one
    void make_intensional_fact(const variable& v) ;
    void make_intensional_fact(const std::string& vname) {
      make_intensional_fact(variable(vname)) ;
    }

    
    void rotate_vars(const std::list<variable> &lvars) ;
    // this method will adjust the rotation vars
    // so that the all the history part will be
    // having the same domain as the current part
    // for example, we have A{n+1}, A{n}, A{n-1}
    // that is rotating. this method will make
    // A{n} and A{n-1} have the same domain as
    // A{n+1}, for those out of A{n+1}'s domain,
    // they are erased, for those missing, they
    // are copied from A{n+1}. this method is
    // mainly for those dynamic rule's rotation
    // list because during the execution, there
    // might have been erase and insertion occurred
    // and we therefore must adjust the history
    // variables accordingly.
    void adjust_rotation_vars(const std::list<variable>& lvars) ;
    void update_remap(const std::vector<std::pair<int, int> > &remap_update);
    entitySet get_allocation(int size) {
      entitySet alloc = interval(maximum_allocated,maximum_allocated+size-1) ;
      maximum_allocated += size ;
      return alloc ;
    }
    entitySet negative_allocation(int size) {
      entitySet alloc = interval(minimum_allocated-size+1, minimum_allocated) ;
      minimum_allocated -= size ;
      return alloc ;
    }

    std::pair<entitySet, entitySet> get_distributed_alloc(int size) ;
    std::pair<entitySet, entitySet> get_distributed_alloc(const std::vector<int> &remap_entities);
    std::pair<entitySet, entitySet> get_distributed_alloc(int size, int offset);
    int is_distributed_start() {return dist_from_start ;}
    std::vector<entitySet>& get_init_ptn() {return init_ptn ;}
    void  put_init_ptn(std::vector<entitySet> &t_init ) {init_ptn = t_init ;}
    
    fact_db::distribute_infoP get_distribute_info() ;
    void put_distribute_info(fact_db::distribute_infoP dp) ;
    bool isDistributed() ;

    void setupDefaults(const rule_db &rdb) ;

    void write_all_hdf5(const char *filename) ;
    void read_all_hdf5(const char *filename) ;
    void write_hdf5(const char *filename, variableSet &vars) ;
    void read_hdf5(const char *filename, variableSet &vars) ;
    void Print_diagnostics() ;
    //the following functions have nothing to do with gKeySpace, gKeyManager
    // experimental code to create keyspace from the
    // global registered keyspace list
    // returns "true" to indicate the methods succeeded,
    // "false" to indicate an error.
    bool create_keyspace(KeySpaceList& global_list) ;
    KeySpaceP get_keyspace(const std::string& kname) ;
    KeyManagerP get_key_manager() const {return key_manager ;}
    void init_key_manager() ;

   
    
  };
  ///////////////////////////////////////////////////////////////////////
  // temporary function
  /////////////////////////////////////////////////////////////////////////
  void reorder_facts(fact_db &facts, dMap &remap) ;
  void serial_freeze(fact_db &facts) ; 
}
    

#endif

    

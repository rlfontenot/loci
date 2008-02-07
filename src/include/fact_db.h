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

#include <string>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <store_rep.h>
#include <variable.h>
#include <Map_rep.h>
#include <Map.h>
#include <DMap.h>

namespace Loci {
  class rule_db ;
  
  class fact_db {
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
    } ;
    std::vector<entitySet> init_ptn ;

    //Related to Duplication: global_comp_entities set defines
    //entities that can be computed on a processor
    entitySet global_comp_entities; 

    typedef CPTR<distribute_info> distribute_infoP ;
    distribute_infoP distributed_info ;

    std::set<std::vector<variableSet> > intensive_output_maps;
  private:
    struct fact_info {
      store_refP data_rep ;
    } ;
    
    std::map<variable,variable> synonyms ;
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      while((mi=synonyms.find(v)) != synonyms.end())
        v = mi->second ;
      return v ;
    }

    // a copy function
    void copy_all_from(const fact_db& f) {
      init_ptn = f.init_ptn ;
      global_comp_entities = f.global_comp_entities ;
      synonyms = f.synonyms ;
      maximum_allocated = f.maximum_allocated ;
      minimum_allocated = f.minimum_allocated ;
      dist_from_start = f.dist_from_start ;
      fmap = f.fmap ;
      tmap = f.tmap ;
      nspace_vec = f.nspace_vec ;
      extensional_facts = f.extensional_facts ;
      /* we cannot use the following direct assignment
         to copy the distributed_info from f since
         distributed_info is a NPTR pointer and is
         reference counted this would not be a true copy
      */
      // distributed_info = f.distributed_info ;
      distribute_infoP& df = distributed_info ;
      const distribute_infoP& fdf = f.distributed_info ;
      if(fdf == 0) {
        df = 0 ;
        return ;
      }
      df = new distribute_info ;
      df->myid = fdf->myid ;
      df->isDistributed = fdf->isDistributed ;
      // we make a deep copy of the maps
      entitySet l2g_alloc = fdf->l2g.domain() ;
      entitySet g2l_alloc = fdf->g2l.domain() ;
      df->l2g.allocate(l2g_alloc) ;
      df->g2l.allocate(g2l_alloc) ;
      for(entitySet::const_iterator ei=l2g_alloc.begin();
          ei!=l2g_alloc.end();++ei) {
        df->l2g[*ei] = fdf->l2g[*ei] ;
      }
      for(entitySet::const_iterator ei=g2l_alloc.begin();
          ei!=g2l_alloc.end();++ei)
        df->g2l[*ei] = fdf->g2l[*ei] ;
      df->my_entities = fdf->my_entities ;
      df->comp_entities = fdf->comp_entities ;
      df->copy = fdf->copy ;
      df->xmit = fdf->xmit ;
      df->copy_total_size = fdf->copy_total_size ;
      df->xmit_total_size = fdf->xmit_total_size ;
      //      df->remap = fdf->remap ;
      df->g2f = fdf->g2f ;
    }
    
    std::pair<entitySet, entitySet> get_dist_alloc(int size) ;
    
    int maximum_allocated ;
    int minimum_allocated ;
    int dist_from_start ;
    
    std::map<variable,fact_info> fmap ;
    std::map<variable,storeRepP> tmap ;
    std::vector<std::string> nspace_vec ;
    // support for multiple queries and experimental
    // extensions to the fact_db to distinguish
    // extensional facts and intensional facts
    variableSet extensional_facts ;

    variable add_namespace(variable v) const ;
    // this is the basic method that creats a fact
    // in the fact_db. It is served as the basis
    // for create_fact & the create_intensional_fact methods
    void create_pure_fact(const variable& v, storeRepP st) ;
  public:
    fact_db() ;
    ~fact_db() ;
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
    void set_maximum_allocated(int i) ;
    void set_minimum_allocated(int i) ;
    int get_max_alloc() {
      return maximum_allocated ;
    }
    int get_min_alloc() {
      return minimum_allocated ;
    }
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

      
    
    // create_fact now defaults to create an extensional fact
    // as this is the primary interface for users of Loci
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
    
    void update_fact(variable v, storeRepP st) ;
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
    
    void remove_variable(variable v) ;
    
    void synonym_variable(variable v, variable synonym) ;

    void rotate_vars(const std::list<variable> &lvars) ;
    
    void set_namespace(std::string name_space){
      nspace_vec.insert(nspace_vec.begin(), 1, name_space) ; 
    }
    void remove_namespace() {
      nspace_vec.pop_back() ;
    }  
    void unset_namespace() {
      nspace_vec.clear() ;
    }
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

    void update_remap(const std::vector<std::pair<int, int> > &remap_update);

    std::pair<entitySet, entitySet> get_distributed_alloc(int size) ;
    std::pair<entitySet, entitySet> get_distributed_alloc(const std::vector<int> &remap_entities);
    std::pair<entitySet, entitySet> get_distributed_alloc(int size, int offset);
    int is_distributed_start() {return dist_from_start ;}
    std::vector<entitySet>& get_init_ptn() {return init_ptn ;}
    void  put_init_ptn(std::vector<entitySet> &t_init ) {init_ptn = t_init ;}
    
    fact_db::distribute_infoP get_distribute_info() ;
    void put_distribute_info(fact_db::distribute_infoP dp) ;
    bool isDistributed() ;
    
    
    variableSet get_typed_variables() const ;
    std::ostream &write(std::ostream &s) const ;
    std::istream &read(std::istream &s) ;
    void setupDefaults(const rule_db &rdb) ;
    std::istream& read_vars(std::istream& s, const rule_db& rdb) ;

    void read_vars(std::string filename, const rule_db &rdb) {
      std::ifstream ifile(filename.c_str(),std::ios::in) ;
      if(ifile.fail()) {
        std::string error = std::string("Can't open '") +filename + "'" ;
        throw(StringError(error)) ;
      }
      read_vars(ifile,rdb) ;
      ifile.close() ;
    }

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

    void write_all_hdf5(const char *filename) ;
    void read_all_hdf5(const char *filename) ;
    void write_hdf5(const char *filename, variableSet &vars) ;
    void read_hdf5(const char *filename, variableSet &vars) ;
    void Print_diagnostics() ;
  } ;
  
  void reorder_facts(fact_db &facts, dMap &remap) ;
  void serial_freeze(fact_db &facts) ;
  
}

#endif

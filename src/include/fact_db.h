#ifndef FACT_DB_H
#define FACT_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <string>
#include <map>
#include <vector>
#include <list>
#include <ostream>
#include <store_rep.h>
#include <variable.h>
//#include <constraint.h>
//#include <rule.h>
//#include <parameter.h>
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
      Map l2g ;
      dMap g2l ; 
      
      entitySet my_entities ;
      entitySet comp_entities;
      
      std::vector<dist_data> copy ; // Data describing which processors own
      // the entities in the clone region
      std::vector<dist_data> xmit ; // Data describing which processors to
      // send entities that this processor owns
      int copy_total_size ;
      int xmit_total_size ;
      dMap remap ;
      distribute_info() {} ;
    } ;
    std::vector<entitySet> init_ptn ;
    entitySet global_comp_entities;
    typedef CPTR<distribute_info> distribute_infoP ;
    distribute_infoP distributed_info ;

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

    fact_db(const fact_db &f) ;
    fact_db &operator=(const fact_db &f) ;
    std::pair<entitySet, entitySet> get_dist_alloc(int size) ;
    
    Map l2g ;
    int maximum_allocated ;
    int minimum_allocated ;
    int dist_from_start ;
    
    std::map<variable,fact_info> fmap ;
    std::map<variable,storeRepP> tmap ;
    std::vector<std::string> nspace_vec ;
  public:
    fact_db() ;
    ~fact_db() ;
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
    
    void create_fact(variable v, storeRepP st) ;
    void create_fact(std::string vname, storeRepP st)
      { create_fact(variable(vname),st) ;}
    void create_fact(variable v, store_instance &si)
      { create_fact(v,si.Rep()) ;si.setRep(get_variable(v)) ; }
    void create_fact(std::string vname, store_instance &si)
      { create_fact(variable(vname),si) ; }
    
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
    Map& get_l2g() { return l2g; } 
    void put_l2g(Map& lg) { l2g = lg ; }
    
    fact_db::distribute_infoP get_distribute_info() ;
    void put_distribute_info(fact_db::distribute_infoP dp) ;
    bool isDistributed() ;
    
    
    variableSet get_typed_variables() const ;
    std::ostream &write(std::ostream &s) const ;
    std::istream &read(std::istream &s) ;
    std::istream& read_vars(std::istream& s, const rule_db& rdb) ;
    
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

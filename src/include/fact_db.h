#ifndef FACT_DB_H
#define FACT_DB_H

#include <string>
#include <map>
#include <vector>
#include <list>
#include <ostream>
#include <store_rep.h>
#include <variable.h>
#include <constraint.h>
#include <rule.h>
#include <parameter.h>
#include <Map_rep.h>
#include <Map.h>

namespace Loci {
  class fact_db {
  public:
    class time_info : public CPTR_type {
      friend class fact_db ;
    protected:
      time_info() {}
    public:
      param<int> time_var ;
      std::list<std::list<variable> > rotate_lists ;
    } ;
    
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
      Map g2l ;

      entitySet my_entities ;

      std::vector<dist_data> copy ; // Data describing which processors own
      // the entities in the clone region
      std::vector<dist_data> xmit ; // Data describing which processors to
      // send entities that this processor owns
      int copy_total_size ;
      int xmit_total_size ;

      distribute_info() {} ;
    } ;
    
    typedef CPTR<time_info> time_infoP ;
    typedef CPTR<distribute_info> distribute_infoP ;
    distribute_infoP distributed_info ;
  private:
    fact_db(const fact_db &f) ;
    fact_db &operator=(const fact_db &f) ;
    struct fact_data {
      variableSet aliases ;
      store_refP data_rep ;
      bool ismap ;
      MapRepP minfo ;
      std::map<entitySet,entitySet> imageMap ;
      std::map<entitySet,std::pair<entitySet,entitySet> > preimageMap
 ;
      fact_data() {} 
      fact_data(variable v, storeRepP &st)
      { aliases += v, data_rep = new store_ref ; (*data_rep) = st ; 
      ismap = (st->RepType() == Loci::MAP);
      if(ismap) minfo = MapRepP(data_rep->getRep()) ; }
    } ;
    
    struct existential_info {
      variable v ;
      entitySet exists ;
      existential_info() {} 
      existential_info(const variable &vin,const entitySet &e) :
        v(vin), exists(e) {}
    } ;
    
    struct fact_info ;
    friend struct fact_info ;
    struct fact_info {
      int fact_info_ref ;
      entitySet fact_installed ;
      variableSet synonyms ;
      std::map<rule,existential_info> exist_map ;
      entitySet existence ;
      entitySet requested ;
      entitySet shadow ; // Used by distributed memory apply rules
      
      fact_info() {fact_info_ref = -1 ; }
      fact_info(int ref)
      { fact_info_ref = ref ; }
    } ;
    
    
    std::map<variable,variable> synonyms ;
    
    typedef std::map<variable, fact_info> vmap_type ;
    vmap_type vmap ;
    std::vector<fact_data> fact_infov ;
    intervalSet free_set ;
    
    variableSet all_vars ;
    
    std::map<time_ident,std::map<variable,std::list<variable> > > time_map;
    
    
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      if((mi=synonyms.find(v)) != synonyms.end())
        return mi->second ;
      return v ;
    }
    
    void install_fact_info(variable v, fact_info info) ;
    void install_fact_data(variable v, fact_data data) ;
    
    const fact_info & get_fact_info(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end()) {
        cerr << "get_fact_info: variable " << v << " does not exist" << endl ;
        abort() ;
      }
      return mi->second ;
    }
    
    const fact_data & get_fact_data(variable v) const 
      { return fact_infov[get_fact_info(v).fact_info_ref] ; }
      
    fact_info & get_fact_info(variable v) ;
    fact_data & get_fact_data(variable v) 
      { return fact_infov[get_fact_info(v).fact_info_ref] ; }
    void register_variable(variable v) ;
    std::vector<entitySet> init_ptn ;
    int maximum_allocated ;
    int dist_from_start ;
  public:

    fact_db() ;
    ~fact_db() ;

    void create_fact(variable v, storeRepP st) ;
    void create_fact(std::string vname, storeRepP st)
      { create_fact(variable(vname),st) ;}
    void create_fact(variable v, store_instance &si)
      { create_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void create_fact(std::string vname, store_instance &si)
      { create_fact(variable(vname),si) ; }
    
    void update_fact(variable v, storeRepP st) ;
    void update_fact(std::string vname, storeRepP st)
      { update_fact(variable(vname),st) ;}
    
    void update_fact(variable v, store_instance &si)
      { update_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void update_fact(std::string vname, store_instance &si)
      { update_fact(variable(vname),si) ; }

    storeRepP get_fact(variable &v) { return get_variable(v); }
    storeRepP get_fact(std::string vname)
      { return get_variable(variable(vname)) ; }
    
    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname,const storeRepP st)
      { set_variable_type(variable(vname),st) ; }
    
    void set_variable_type(variable v, store_instance &si)
      { set_variable_type(v,si.Rep()) ; }
    void set_variable_type(std::string vname, store_instance &si)
      { set_variable_type(variable(vname),si) ; }
    
    void allocate_variable(variable v, entitySet s) ;
    void allocate_variable(std::string vname, const entitySet s)
      { allocate_variable(variable(vname),s) ; }
    
    void variable_is_fact_at(variable v,entitySet s) ;
    void variable_is_fact_at(std::string vname, const entitySet s)
      { variable_is_fact_at(variable(vname),s) ; }
    
    void alias_variable(variable v, variable alias) ;
    void alias_variable(std::string vname, std::string alias)
      { alias_variable(variable(vname),variable(alias)) ; }

    void synonym_variable(variable v, variable synonym) ;
    void synonym_variable(std::string vname, std::string synonym)
      { synonym_variable(variable(vname),variable(synonym)) ; }

    
    entitySet get_allocation(int size) {
      entitySet alloc = interval(maximum_allocated,maximum_allocated+size-1) ;
      maximum_allocated += size ;
      return alloc ;
    }
    std::pair<interval, interval> get_distributed_alloc(int size) ;
    int is_distributed_start() {
      return dist_from_start ;
    }
    std::vector<entitySet>& get_init_ptn() {
      return init_ptn ;}
    storeRepP get_variable(variable v) ;
    storeRepP get_variable(std::string vname)
      { return get_variable(variable(vname)) ; }
    
    fact_db::time_infoP get_time_info(time_ident tl) ;
    fact_db::distribute_infoP get_distribute_info() ;
    void put_distribute_info(fact_db::distribute_infoP dp) ;
    bool isDistributed() ;
    void initialize_time(time_infoP ti) ;
    void advance_time(time_infoP ti) ;
    void close_time(time_infoP ti) ;
    
    variableSet get_typed_variables() const { return all_vars ; }
    
    variableSet get_synonyms(variable v) const 
      { return get_fact_info(v).synonyms ; }
    
    variableSet get_aliases(variable v) const
      { return get_fact_data(v).aliases ; }
    
    void set_existential_info(variable v,rule f,entitySet x) ;
    ruleSet get_existential_rules(variable v) ;

    entitySet get_existential_info(variable v, rule f) {
      fact_info &finfo = get_fact_info(v) ;
      std::map<rule,existential_info>::const_iterator mi ;
      mi = finfo.exist_map.find(f) ;
      if(mi!=finfo.exist_map.end()) {
        return mi->second.exists ;
      } else
        return EMPTY ;
    }

    entitySet variable_existence(variable v) ;
    void variable_request(variable v, entitySet e) ;
    void set_variable_shadow(variable v, entitySet e)
    { get_fact_info(v).shadow = e ; }
    void variable_shadow(variable v, entitySet e)
    { get_fact_info(v).shadow += e ; }
    entitySet get_variable_shadow(variable v) const
    { return get_fact_info(v).shadow ; }
    entitySet get_variable_request(rule f, variable v) ;
    entitySet get_variable_requests(variable v) const 
    { return  get_fact_info(v).requested ; }
    
    bool is_a_Map(variable v) ;
    entitySet image(variable v, entitySet e) ;
    std::pair<entitySet,entitySet> preimage(variable v, entitySet e) ;

    ostream &write(ostream &s) const ;
    istream &read(istream &s) ;
    
    void write_hdf5(const char *filename);
    void read_hdf5(const char *filename);
    void printSummary(std::ostream &s) const ;
  } ;
  
  void reorder_facts(fact_db &facts, Map &remap) ;
  
  
}

#endif

#ifndef SCHED_DB_H
#define SCHED_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <fact_db.h>
namespace Loci {
  class sched_db {
    
    struct sched_data {
      variableSet aliases ;
      bool ismap ;
      MapRepP minfo ;
      std::map<entitySet,entitySet> imageMap ;
      std::map<entitySet,std::pair<entitySet,entitySet> > preimageMap
      ;
      sched_data() {} 
      sched_data(variable v, storeRepP &st)
	{ aliases += v ; 
	ismap = (st->RepType() == Loci::MAP);
	if(ismap) minfo = MapRepP(st->getRep()) ; }
    } ;
    
    struct existential_info {
      variable v ;
      entitySet exists ;
      existential_info() {} 
      existential_info(const variable &vin,const entitySet &e) :
        v(vin), exists(e) {}
    } ;
    struct sched_info ;
    friend struct sched_info ;
    struct sched_info {
      int sched_info_ref ;
      entitySet fact_installed ;
      variableSet synonyms ;
      std::map<rule,existential_info> exist_map ;
      entitySet existence ;
      entitySet requested ;
      entitySet shadow ; // Used by distributed memory apply rules
      
      sched_info() {sched_info_ref = -1 ; }
      sched_info(int ref)  
	{ sched_info_ref = ref ; }
    } ;
    void register_variable(variable v) ;
  
    variableSet all_vars ;
    std::map<variable,variable> synonyms ;
    typedef std::map<variable, sched_info> vmap_type ;
    vmap_type vmap ;
    std::vector<sched_data> sched_infov ;
    intervalSet free_set ;
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      if((mi=synonyms.find(v)) != synonyms.end()) 
	return mi->second ;
      return v ;
    }
  public:
    sched_db() ;
    ~sched_db() ;
    sched_db(fact_db &facts) ; 
    void install_sched_data(variable v, sched_data data) ;
    void install_sched_info(variable v, sched_info info) ;
    
    variableSet get_synonyms(variable v) const 
      { return get_sched_info(v).synonyms ; }
    
    variableSet get_aliases(variable v) const
      { return get_sched_data(v).aliases ; }
    
    void set_variable_type(variable v, storeRepP st, fact_db &facts) ;
    void set_variable_type(std::string vname,const storeRepP st, fact_db &facts)
      { set_variable_type(variable(vname),st, facts) ; }
    
    void set_variable_type(variable v, store_instance &si, fact_db &facts)
      { set_variable_type(v,si.Rep(), facts) ; }
    void set_variable_type(std::string vname, store_instance &si, fact_db &facts)
      { set_variable_type(variable(vname),si, facts) ; }
    
    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname,const storeRepP st)
      { set_variable_type(variable(vname),st) ; }
    
    void set_variable_type(variable v, store_instance &si)
      { set_variable_type(v,si.Rep()) ; }
    void set_variable_type(std::string vname, store_instance &si)
      { set_variable_type(variable(vname),si) ; } 
    
    void variable_is_fact_at(variable v,entitySet s, fact_db &facts) ;
    void variable_is_fact_at(std::string vname, const entitySet s, fact_db &facts)
      { variable_is_fact_at(variable(vname),s, facts) ; }
    
    void variable_is_fact_at(variable v,entitySet s) ;
    void variable_is_fact_at(std::string vname, const entitySet s)
      { variable_is_fact_at(variable(vname),s) ; }
    
    void alias_variable(variable v, variable alias, fact_db &facts) ;
    void alias_variable(std::string vname, std::string alias, fact_db &facts)
      { alias_variable(variable(vname),variable(alias), facts) ; }
   
    void alias_variable(variable v, variable alias) ;
    void alias_variable(std::string vname, std::string alias)
      { alias_variable(variable(vname),variable(alias)) ; }
    
    void synonym_variable(variable v, variable synonym, fact_db &facts) ;
    void synonym_variable(std::string vname, std::string synonym, fact_db &facts)
      { synonym_variable(variable(vname),variable(synonym), facts) ; }
    
    void synonym_variable(variable v, variable synonym) ;
    void synonym_variable(std::string vname, std::string synonym)
      { synonym_variable(variable(vname),variable(synonym)) ; }
    
    void set_existential_info(variable v,rule f,entitySet x) ;
    ruleSet get_existential_rules(variable v) ;
    const sched_info & get_sched_info(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end()) {
        cerr << "const get_sched_info: variable " << v << " does not exist" << endl ;
        abort() ;
      }
      return mi->second ;
    }
    const sched_data & get_sched_data(variable v) const 
  { return sched_infov[get_sched_info(v).sched_info_ref] ; }
    
    sched_info & get_sched_info(variable v) ;
    sched_data & get_sched_data(variable v) 
      { return sched_infov[get_sched_info(v).sched_info_ref] ; }

    bool is_a_Map(variable v) {
      return get_sched_data(v).ismap ;
    }
    entitySet get_existential_info(variable v, rule f) {
      sched_info &finfo = get_sched_info(v) ;
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
    { get_sched_info(v).shadow = e ; }
    void variable_shadow(variable v, entitySet e)
    { get_sched_info(v).shadow += e ; }
    entitySet get_variable_shadow(variable v) const
    { return get_sched_info(v).shadow ; }
    entitySet get_variable_request(rule f, variable v) ;
    entitySet get_variable_requests(variable v) const 
    { return  get_sched_info(v).requested ; }
    variableSet get_typed_variables() const { return all_vars ; }
    entitySet image(variable v, entitySet e) ;
    std::pair<entitySet,entitySet> preimage(variable v, entitySet e) ;

    std::ostream &print_summary(fact_db &facts, std::ostream &s) ;
  } ;
}

#endif

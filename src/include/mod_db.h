#ifndef MOD_DB_H
#define MOD_DB_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#ifdef HAVE_DLOPEN
#include <dlfcn.h>
#endif
#include <rule.h>
#include <fact_db.h>
#include <string>
namespace Loci {
  std::string remove_space(const std::string &str) ; 
  void parse_str(const std::string& str, std::vector<std::string> &str_vec) ;
  class mod {
  public: 
    struct mod_info {
     rule_impl_list loaded_rule_list ;
      std::string mod_name ;
      void *m_library ;
      void (*m_init_model)(fact_db &facts, const char *problem_name) ;
      std::string name() { return mod_name ; } 
      mod_info(rule_impl_list& rl, std::string str) {
	loaded_rule_list.copy_rule_list(rl) ;
	mod_name = str ;
      }
      mod_info(rule_impl_list& rl) { 
	loaded_rule_list.copy_rule_list(rl) ;
	mod_name = std::string("") ;
      }
      mod_info() {
	mod_name = "" ;
	m_library = 0 ;
	m_init_model = 0 ;
      }
      mod_info(const mod_info& mi) {
	loaded_rule_list.copy_rule_list(mi.loaded_rule_list) ;
	mod_name = mi.mod_name ;
        m_library = mi.m_library ;
        m_init_model = mi.m_init_model ;
      }
      ~mod_info() { } 
    } ;
  private: 
    friend class mod_info ;
    struct mod_db {
      std::map<std::string,mod_info> mod_map ;
      typedef std::map<std::string,mod_info>::iterator MI ;
      void put_info(mod_info& md) {
	MI msi ;
	if((msi = mod_map.find(md.name())) == mod_map.end())
	  mod_map.insert(std::make_pair(md.name(),md ));
      }
      
      mod_info &get_info(const std::string &str) ;
      mod_info &get_info(const std::string &str, const std::string &to_str, const char* problem_name, fact_db &facts) ;
    } ;
  
    static mod_db *mdb ;
    void create_mod_db() {if(0==mdb) mdb = new mod::mod_db ; }
  public:
    mod(const std::string& name) {
      create_mod_db() ;
      get_info(name) ;
    }
    mod(const std::string& name, const std::string& to_str, const char* problem_name, fact_db &facts) {
      create_mod_db() ;
      get_info(name, to_str, problem_name, facts) ;
    }
    mod() {
     create_mod_db() ;
     mod_info md = mod_info(init_rule_list) ;
    }
    mod_info &get_info(const std::string name) { 
      return mdb->get_info(name) ; }
    mod_info &get_info(const std::string name, const std::string &to_str, const char* problem_name, fact_db &facts) {
      return mdb->get_info(name, to_str, problem_name, facts); }
    void put_info(mod_info &md) {mdb->put_info(md) ; }
  } ;
  
  void load_module(const std::string from_str, const std::string to_str, rule_db& rdb, std::set<std::string> &str_set) ;
  
  void load_module(const std::string from_str, const std::string to_str, const char* problem_name, fact_db &facts, rule_db& rdb, std::set<std::string> &str_set) ;	 
}

#endif

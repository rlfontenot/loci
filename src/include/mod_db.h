#ifndef MOD_DB_H
#define MOD_DB_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
namespace Loci {
  using Loci::rule_impl_list ;
  using Loci::variable ;
  std::string remove_space(const std::string &str) {
    std::string tmp_str ;
    for(std::string::const_iterator si = str.begin(); si != str.end(); ++si) {
      if(*si != ' ')
	tmp_str += *si ;
    }
    return tmp_str ;
  }
  void parse_str(const std::string& str, std::vector<std::string> &str_vec) {
    if(!str.empty()) {
      int tmp = 0 ;
      std::string sub_str = remove_space(str) ;
      while(tmp != std::string::npos) {
	tmp = sub_str.find(",") ;
	str_vec.push_back(sub_str.substr(0,tmp)) ;
      sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
      }
    }
  }
  class mod {
  public: 
   struct mod_info {
     rule_impl_list loaded_rule_list ;
     std::string mod_name ;
     std::vector<std::string> default_ns_vec, using_ns_vec;
     std::vector<variable> export_vars, import_vars ;
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
     mod_info() {mod_name = "" ;
     void *m_library = 0 ;
     void (*m_init_model)(fact_db &facts, const char *problem_name) = 0 ;
     }
     ~mod_info() { } ;
     mod_info(const mod_info& mi) {
       loaded_rule_list.copy_rule_list(mi.loaded_rule_list) ;
       mod_name = mi.mod_name ;
     }
     
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
     mod_info &get_info(const std::string &str) { 
       MI msi ;
       mod_info md ;
       std::string tmp_str ;
       tmp_str.append(str) ;
       tmp_str.append("_m.so") ;
       if((msi = mod_map.find(tmp_str)) == mod_map.end()) {
	 md.default_ns_vec.push_back(str) ;
	 if(Loci::MPI_rank == 0)
	   cout << "Loading  in rules from  " << tmp_str << endl ;
	 md.m_library = dlopen(tmp_str.c_str(),RTLD_GLOBAL|RTLD_NOW) ;
	 if(md.m_library == 0) {
	   if(Loci::MPI_rank == 0)
	     cerr << "unable to open " << tmp_str.c_str() << endl ;
	   const char *error = dlerror() ;
	   if(error)
	     if(Loci::MPI_rank == 0)
	       cerr << "reason for failure is " << error << endl ;
	   exit(-1) ;
	 }
	 //for(rule_impl_list::iterator gi = global_rule_list.begin(); gi !=global_rule_list.end(); ++gi) {
	 //if(!(gi.get_p())->rr->is_module_rule()) 
	 //  cout << str << " global rule added = " << rule(*gi) << endl ; 
	 //}
	 md.loaded_rule_list.copy_rule_list(global_rule_list) ;
	 md.mod_name = tmp_str ;
	 global_rule_list.clear() ;
	 put_info(md) ;
	 return mod_map[tmp_str] ;
       }
       return msi->second ;
     } 
     mod_info &get_info(const std::string &str, const std::string &to_str, const char* problem_name, fact_db &facts) { 
       MI msi ;
       mod_info md ;
       std::string tmp_str ;
       tmp_str.append(str) ;
       tmp_str.append("_m.so") ;
       if((msi = mod_map.find(tmp_str)) == mod_map.end()) {
	 md.default_ns_vec.push_back(str) ;
	 if(Loci::MPI_rank == 0)
	   cout << "Loading  in rules from  " << tmp_str << endl ;
	 md.m_library = dlopen(tmp_str.c_str(),RTLD_GLOBAL|RTLD_NOW) ;
	 if(md.m_library == 0) {
	   if(Loci::MPI_rank == 0)
	     cerr << "unable to open " << tmp_str.c_str() << endl ;
	   const char *error = dlerror() ;
	   if(error)
	     if(Loci::MPI_rank == 0)
	       cerr << "reason for failure is " << error << endl ;
	   exit(-1) ;
	 }
	 md.m_init_model = (void (*)(fact_db &facts, const char *problem_name))
	   dlsym(md.m_library,"init_model") ;
	 if(md.m_init_model != 0) {
	   if(!to_str.empty()) {
	     int tmp = 0 ;
	     std::string sub_str = to_str ;
	     while(tmp != std::string::npos) {
	       tmp = sub_str.find("_") ;
	       facts.set_namespace(sub_str.substr(0,tmp)) ;
	     sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
	     }
	   }
	   md.m_init_model(facts,problem_name) ;
	   facts.unset_namespace() ;
	 }
	 md.loaded_rule_list.copy_rule_list(global_rule_list) ;
	 md.mod_name = tmp_str ;
	 global_rule_list.clear() ;
	 put_info(md) ;
	 return mod_map[tmp_str] ; 
       }
       return msi->second ;
     } 
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
  mod::mod_db *mod::mdb = 0 ;
  
  void load_module(const std::string from_str, const std::string to_str, rule_db& rdb, std::set<std::string> &str_set) {
    str_set.insert(from_str) ;
    mod md(from_str) ;
    mod::mod_info m = md.get_info(from_str) ;
    for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
      if(!(gi.get_p())->rr->is_module_rule()) {
	if(!to_str.empty()) {
	  int tmp = 0 ;
	  std::string sub_str = to_str ;
	  while(tmp != std::string::npos) {
	    tmp = sub_str.find("_") ;
	    m.using_ns_vec.push_back(sub_str.substr(0,tmp)) ;
	    sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
	  }
	  rule_implP rp = *gi;
	  variableSet vars = rp->get_var_list() ;
	  std::map<variable,variable> new_vars;
	  for(variableSet::variableSetIterator i=vars.begin();i!=vars.end();++i) {
	    variable tmp_var = *i ;
	    for(int j = m.using_ns_vec.size()-1; j >= 0 ; --j)
	      tmp_var =  tmp_var.add_namespace(m.using_ns_vec[j]); 
	    new_vars[*i] = tmp_var ;
	  }
	  rp->rename_vars(new_vars) ;
	  rdb.add_rule(Loci::rule(rp)) ; 
	} 
	else
	  rdb.add_rule(Loci::rule(*gi)) ; 
      }
      else {
	if(Loci::MPI_rank == 0)
	  cerr << "Module rule found in " << from_str << endl ;
	std::string load  =  ((Loci::register_module*)(gi.get_p()->rr))->using_nspace() ;
	std::vector<std::string> str_vec ;
	parse_str(load, str_vec) ;
	for(int i = 0; i < str_vec.size(); ++i) 
	  if(str_set.find(str_vec[i]) == str_set.end()) {
	    if(Loci::MPI_rank == 0)
	      cout << "loading in rules from " << str_vec[i] <<"  for module " << to_str << endl ; 
	    load_module(str_vec[i], to_str, rdb, str_set) ;
	  }
      }
    }	 
  }
  
  void load_module(const std::string from_str, const std::string to_str, const char* problem_name, fact_db &facts, rule_db& rdb, std::set<std::string> &str_set) {
    str_set.insert(from_str) ;
    mod md(from_str, to_str, problem_name, facts) ;
    mod::mod_info m = md.get_info(from_str, to_str, problem_name, facts) ;
    variableSet input_vars, output_vars ;
    if(!to_str.empty()) {
      int tmp = 0 ;
      std::string sub_str = to_str ;
      while(tmp != std::string::npos) {
	tmp = sub_str.find("_") ;
	m.using_ns_vec.push_back(sub_str.substr(0,tmp)) ;
	sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
      }
      for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
	if((gi.get_p())->rr->is_module_rule()) {
	  if(Loci::MPI_rank == 0)
	    cerr << "Module rule found in " << from_str << endl ;
	  std::string load  =  ((Loci::register_module*)(gi.get_p()->rr))->using_nspace() ;
	  std::vector<std::string> str_vec ;
	  parse_str(load, str_vec) ;
	  for(int i = 0; i < str_vec.size(); ++i) 
	    if(str_set.find(str_vec[i]) == str_set.end()) {
	      if(Loci::MPI_rank == 0)
		cout << "loading in rules from " << str_vec[i] <<"  for module " << to_str << endl ; 
	      load_module(str_vec[i], to_str, problem_name, facts, rdb, str_set) ;
	    }
	  load = ((Loci::register_module*)(gi.get_p()->rr))->input_vars() ;
	  input_vars = variableSet(expression::create(load)) ;
	  load = ((Loci::register_module*)(gi.get_p()->rr))->output_vars() ;
	  output_vars = variableSet(expression::create(load)) ;
	  break ;
	}
      }
      for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
	if(!(gi.get_p())->rr->is_module_rule()) {
	  rule_implP rp = *gi ;
	  variableSet vars = rp->get_var_list() ;
	  std::map<variable,variable> new_vars;
	  for(variableSet::variableSetIterator i=vars.begin();i!=vars.end();++i)
	    if(input_vars.inSet(*i) || output_vars.inSet(*i) || i->is_time_variable())
	      new_vars[*i] = *i ;
	    else {
	      variable tmp_var = *i ;
	      for(int j = m.using_ns_vec.size()-1; j >= 0; --j) 
		tmp_var = tmp_var.add_namespace(m.using_ns_vec[j]); 
	      new_vars[*i] = tmp_var ;
	    }
	  rp->rename_vars(new_vars) ;
	  rdb.add_rule(Loci::rule(rp)) ;
	}
      }
    }
    
  }	 
}

#endif

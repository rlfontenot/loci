#ifndef EXECUTE_H
#define EXECUTE_H

#include <ostream>
#include <vector>

#include <Config/conf.h>
#include <Tools/cptr.h>

namespace Loci {
  class fact_db ;
  class execute_modules : public CPTR_type {
  protected:
    bool control_thread ;
  public:
    execute_modules() {control_thread = false ; }
    virtual void execute(fact_db &facts) = 0 ;
    virtual void Print(std::ostream &s) const = 0 ;
    bool is_control_thread() { return control_thread; }
  } ;

  typedef CPTR<execute_modules> executeP ;

  class execute_list : public execute_modules {
    std::vector<executeP> elist ;
  public:
    execute_list() {control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };
    
  class execute_sequence : public execute_modules {
    std::vector<executeP> elist ;
  public:
    execute_sequence() {control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };

  class execute_par : public execute_modules {
    std::vector<executeP> elist ;
  public:
    execute_par() {control_thread = true ;}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };

  class execute_create_threads : public execute_modules {
    int num_threads ;
  public:
    execute_create_threads() {control_thread = true ; }
    execute_create_threads(int nth) : num_threads(nth) {} 
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  class execute_destroy_threads : public execute_modules {
  public:
    execute_destroy_threads() {control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  class execute_thread_sync : public execute_modules {
  public:
    execute_thread_sync() {control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  
  const int max_threads = 32 ;

}

#endif

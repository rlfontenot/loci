#ifndef EXECUTE_H
#define EXECUTE_H

#include <ostream>
#include <list>

#include <Tools/cptr.h>

namespace Loci {
  class fact_db ;
  class execute_modules : public CPTR_type {
  public:
    virtual void execute(fact_db &facts) = 0 ;
    virtual void Print(std::ostream &s) const = 0 ;
  } ;

  typedef CPTR<execute_modules> executeP ;

  class execute_list : public execute_modules {
    std::list<executeP> elist ;
  public:
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };

}

#endif

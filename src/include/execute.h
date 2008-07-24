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
#ifndef EXECUTE_H
#define EXECUTE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <ostream>
#include <vector>
#include <string>
using std::string;

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
	virtual string getName() = 0;
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
	virtual string getName() { return "execute_list";};
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
	virtual string getName() { return "execute_sequence";};
  };

}
#endif

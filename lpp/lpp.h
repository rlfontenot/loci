//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#ifndef LPP_H
#define LPP_H

#include "variable.h"
#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>

extern bool prettyOutput ;

struct parseError {
  std::string error_type ;
  parseError(std::string errs) : error_type(errs) {}
} ;
  
class parseFile {
  int cnt ;
  std::string filename ;
  int line_no ;
  std::ifstream is ;
  std::map<Loci::variable,std::pair<std::string,std::string> > type_map ;
  int killsp() ;
  int killspout(std::ostream &outputFile) ;

  void syncFile(std::ostream &outputFile) {
    if(!prettyOutput)
      outputFile << "#line " << line_no << " \"" << filename << "\"" << std::endl ;
  }

  void validate_VariableAccess(Loci::variable v, 
			       const std::list<Loci::variable> &vlist,
			       bool first_name,
			       const std::map<Loci::variable,std::string> &vnames,
			       const std::set<std::list<Loci::variable> > &validate_set) ;

  std::string process_String(std::string instring,
			     const std::map<Loci::variable,std::string> &vnames,
			     const std::set<std::list<Loci::variable> > &validate_set) ;
  
  void process_SpecialCommand(std::ostream &outputFile,
                              const std::map<Loci::variable,std::string> &vnames,
                              int &openbrace) ;
  void process_Prelude(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Compute(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Calculate(std::ostream &outputFile,
                         const std::map<Loci::variable,std::string> &vnames,
                         const std::set<std::list<Loci::variable> > & validate_set) ;
  void setup_Type(std::ostream &outputFile) ;
  void setup_Rule(std::ostream &outputFile) ;
public:
  parseFile() {
    line_no = 0 ;
    cnt = 0 ;
    Loci::variable OUTPUT("OUTPUT") ;
    type_map[OUTPUT] = std::pair<std::string,std::string>("param","<bool>") ;
  }
  void processFile(std::string file, std::ostream &outputFile) ;
} ;

extern std::list<std::string> include_dirs ;
#endif

#ifndef LPP_H
#define LPP_H

#include <Loci.h>
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

  std::string process_String(std::string instring,
                       const std::map<Loci::variable,std::string> &vnames) ;
                        
  void process_Prelude(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Compute(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Calculate(std::ostream &outputFile,
                         const std::map<Loci::variable,std::string> &vnames) ;
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

#ifndef LPP_H
#define LPP_H

#include <Loci.h>
#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>

struct parseError {
  std::string error_type ;
  parseError(std::string errs) : error_type(errs) {}
} ;
  
class parseFile {
  std::string filename ;
  int line_no ;
  std::ifstream is ;
  std::map<Loci::variable,std::pair<std::string,std::string> > type_map ;
  void killsp() {
    while(is.peek() == ' ' || is.peek() == '\t' || is.peek() == '\n'
          || is.peek() == '\r') {
      if(is.peek() == '\n') line_no++ ;
      is.get();
    }
  }
  void killspout(std::ostream &outputFile) {
    while(is.peek() == ' ' || is.peek() == '\t' || is.peek() == '\n'
          || is.peek() == '\r') {
      if(is.peek() == '\n') line_no++ ;
      char c = is.get();
      outputFile << c ;
    }
  }
  void syncFile(std::ostream &outputFile) {
    outputFile << "#line " << line_no << " \"" << filename << "\"" << std::endl ;
  }

  void process_Compute(std::ostream &outputFile,
                         const std::map<Loci::variable,std::string> &vnames) ;
  void process_Calculate(std::ostream &outputFile,
                         const std::map<Loci::variable,std::string> &vnames) ;
  void setup_Type(std::ostream &outputFile) ;
  void setup_Rule(std::ostream &outputFile) ;
public:
  parseFile() {
    line_no = 0 ;
  }
  void processFile(std::string file, std::ostream &outputFile) ;
} ;

#endif

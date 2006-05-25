#include "lpp.h"

using namespace std ;

int main(int argc, char *argv[]) {

  string filename = "test.loci" ;
  
  if(argc >= 2) {
    filename = string(argv[1]) ;
  }

  parseFile parser ;
  try {
    parser.processFile(filename,cout) ;
  } catch(parseError pe) {
    if(pe.error_type != "syntax error")
      cerr << pe.error_type << endl ;
    exit(-1) ;
  }

  return 0 ;
}

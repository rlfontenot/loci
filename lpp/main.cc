#include "lpp.h"

using namespace std ;

int main(int argc, char *argv[]) {
  Loci::Init(&argc,&argv) ;
  string filename = "test.loci" ;
  
  if(argc >= 2) {
    filename = string(argv[1]) ;
  }
  parseFile parser ;

  parser.processFile(filename,cout) ;

  Loci::Finalize() ;
  return 0 ;
}

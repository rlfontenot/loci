#include "lpp.h"
#include <unistd.h>

using namespace std ;

list<string> include_dirs ;

bool prettyOutput = false ;

void Usage(int argc, char *argv[]) {
  cerr << "Loci Pre-Processor Usage:" << endl ;
  cerr << argv[0] <<" -I<dir> filename.loci -o filename.cc" << endl ;
  exit(-1) ;
}

int main(int argc, char *argv[]) {


  bool file_given = false; 
  string filename ;
  bool out_given = false ;
  string outfile ;
  for(int i=1;i<argc;++i) {
    if(argv[i][0] == '-') {
      if(argv[i][1] == 'I') {
        string dir = &argv[i][2] ;
        include_dirs.push_back(dir) ;
      }
      if(argv[i][1] == 'p') {
        prettyOutput = true ;
      }
      if(argv[i][1] == 'o') {
        if(i+1>argc || out_given)
          Usage(argc,argv) ;
        outfile = argv[i+1] ;
        i++ ;
        out_given = true ;
      }
      if(argv[i][1] == 'v') {
        cout << "Loci version: " << Loci::version() << endl ;
      }
      if(argv[i][1] == 'V') {
        cout << "Loci version: " << Loci::version() << endl ;
      }
    } else {
      if(file_given == true) {
        cerr << "multiple filenames given" << endl ;
        Usage(argc,argv) ;
      }
      filename = argv[i] ;
      file_given = true ;
    }
  }
  if(!file_given) {
    cerr << "no filename" << endl ;
    Usage(argc,argv) ;
  }
  parseFile parser ;
  try {
    if(out_given) {
      ofstream file(outfile.c_str(),ios::out) ;
      if(file.fail()) {
        cerr << "unable to open file " << outfile << " for writing!" << endl ;
        exit(-1) ;
      }
      parser.processFile(filename,file) ;
    } else {
      parser.processFile(filename,cout) ;
    }
  } catch(parseError pe) {
    if(pe.error_type != "syntax error")
      cerr << pe.error_type << endl ;
    if(out_given)
      unlink(outfile.c_str()) ;
    exit(-1) ;
  } catch(...) {
    cerr << "Unknown exception caught!" << endl ;
    if(out_given)
      unlink(outfile.c_str()) ;
    exit(-1) ;
  }

  return 0 ;
}

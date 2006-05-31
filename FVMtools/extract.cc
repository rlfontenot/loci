#include <Loci.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
using std::istream ;
using std::ostream ;
using std::endl ;
using std::cin ;
using std::cout ;
using std::cerr ;
using std::ios ;
using std::ofstream ;
using std::ifstream ;

using std::istringstream ;
using std::ostringstream ;
using std::string ;
using std::char_traits ;

#include <vector>
using std::vector ;


void usage(char *s) {
  cout << "Usage: " << endl
       << string(s) << " <package> [package options] <problem_name> <time_step> <key(s)> " << endl << endl ;
  cout << "where <package> may be:" << endl
       << "-2d :  extract for the 2dgv plotting package" << endl
       << "-fv :  extract for the FieldView post-processing package" << endl
       << "-en :  extract for the Ensight post-processing package" << endl
       << "-tec:  extract for the TecPlot post-procesing package" << endl
       << "-ascii: extract to an ascii file" << endl
       << endl ;
  cout << "where <key> may be: " << endl
       << "r - extract density " << endl
       << "p - extract log10 pressure" << endl
       << "P - extract actual pressure" << endl 
       << "u - extract velocity magnitude" << endl
       << "m - extract mach number" << endl
       << "t - extract temperature" << endl
       << "a - extract soundspeed" << endl
       << "G - extract Gamma" << endl
       << "R - extract Gas Constant R~" << endl
       << "f<species> - extract mass fractions for <species>" << endl
       << "n<X> - extract nodal data dump file: output/X_hdf5.0"<<endl 
       << endl ;

  cout << "extra options for -2d postprocessing package" << endl 
       << "  -bc num  : specify which boundary to extract for (multiple -bc's"
       << endl 
       << "             are merged)" << endl 
       << "  -xy : project boundary onto z=0 plane (default)" << endl
       << "  -yz : project boundary onto x=0 plane" << endl
       << "  -xz : project boundary onto y=0 plane" << endl
       << "  -xr : project boundary onto x,radius plane (for axisymmetric grids)"
       << endl << endl ;
  cout << "example:  to extract OH species from time step 0 of ssme simulation for visualization with 2dgv use:" << endl
       << string(s) << " -2d -bc 1 -xr ssme 0 fOH" << endl ;
}


inline double square(double x) { return x*x; }

int FieldViewExtract(int ac, char *av[]) ;

int TecplotExtract(int ac, char *av[]) ;
int twoDExtract(int ac, char *av[]) ;
int asciiExtract(int ac, char *av[]) ;
int EnsightExtract(int ac, char *av[]) ;

int main(int ac, char *av[])
{
  Loci::Init(&ac, &av) ;
  
  if(ac>1 && !strcmp(av[1],"-fv")) {
    av[1] = av[0] ;
    av++ ;
    ac-- ;
    bool success = !FieldViewExtract(ac,av) ;
    Loci::Finalize() ;
    return success ;
   
  } else if(ac>1 && !strcmp(av[1],"-tec")) {
    av[1] = av[0] ;
    av++ ;
    ac-- ;
    bool success = !TecplotExtract(ac,av) ;
    Loci::Finalize() ;
    return success ;
  } else if(ac>1 && !strcmp(av[1],"-2d")) {
    av[1] = av[0] ;
    av++ ;
    ac-- ;
    bool success = !twoDExtract(ac,av) ;
    Loci::Finalize() ;
    return success ;
  } else if(ac>1 && !strcmp(av[1],"-en")) {
    av[1] = av[0] ;
    av++ ;
    ac-- ;
    bool success = !EnsightExtract(ac,av) ;
    Loci::Finalize() ;
    return success ;
  } else if(ac >1 && !strcmp(av[1],"-ascii")) {
    av[1] = av[0] ;
    av++ ;
    ac-- ;
    bool success = !asciiExtract(ac,av) ;
    Loci::Finalize() ;
    return success ;
    
  } 
      
  
  usage(av[0]) ;
  
  Loci::Finalize() ;
  return -1 ;
}

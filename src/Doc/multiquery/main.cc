// This is an example file (main.cc) to illustrate the new Loci's feature
// of making multiple queries for stationary time facts. Basically everything
// is unchanged. See the comments before the queries at the end of the file 
// for detail.

#include <Loci.h>
#include <Tools/fpe.h>
#include "gridReader.h"

#include <sys/stat.h>
#include <unistd.h>
#include <string>

using std::string ;
using std::cout ;
using std::endl ;
using std::cerr ;

int main(int ac, char *av[]) {
  Loci::Init(&ac, &av) ;
  
  // default query 
  string query = "solution" ;

  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 2 && !strcmp(av[1],"-q")) {
      query = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
        Loci::Finalize() ;
        exit(0) ;
      }
      ac-- ;
      av++ ;
    } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }
  
  // if output directory doesn't exist, create one
  struct stat statbuf ;
  if(stat("output",&statbuf))
    mkdir("output",0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file 'output' should be a directory!, rename 'output' and start again."
           << endl ;
      Loci::Abort() ;
    }
  
  // set up floating point exception handlers
  set_fpe_abort() ;
  

  if(ac <= 1) {
    cout << "Loci version: " << Loci::version() << endl ;
    Loci::Finalize() ;
    exit(0) ;
  }


  int av1sz = strlen(av[1]) ;
  if(av1sz>0 && (av[1][av1sz-1] == '.'))
    av[1][av1sz-1] = '\0' ;

  cout << "creating rule database" <<endl ;
  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;

  cout << "reading vars file" << endl ;
  // read grid, connectivity information, and user supplied information into
  // the fact database
  fact_db facts ;
  heat::read_grid(facts,rdb,av[1]) ;

#define MULTIQUERY
  // Making multiple queries are essentially the same as making
  // a single query, except for that we can issue many of them
  // right now. If the user is quering an extensional fact, then
  // Loci won't generate an execution plan, it will simply report
  // that you are quering an extensional fact. For any query, the
  // queried results are maintained in the fact_db until the next
  // query. Before every query, all intensional facts in the 
  // fact_db will be "erased."
#ifdef MULTIQUERY
  if(!Loci::makeQuery(rdb,facts,"area")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"vol")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"ic_stime")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"boundary_node")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"nodalw_sum")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"Kcoeff")) {
    cerr << "query failed!" << endl ;
    //Loci::Abort() ;
    Loci::Finalize() ;
    return -1 ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"area")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,"area")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
#endif
  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  cout << endl ;
  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }

  
  Loci::Finalize() ;
}


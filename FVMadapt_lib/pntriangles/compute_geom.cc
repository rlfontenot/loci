#include <Loci.h>
#include <Tools/fpe.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

using std::istringstream ;
using std::map ;

using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::ifstream ;
using std::ios ;
using std::cout;

void Usage() {
  cerr << "Usage:" << endl ;
  cerr << " computegeom <options> <casename>" << endl
       << " where <options> are:" << endl
       << " -sym <bc>: specify symmetry boundary conditions" << endl
       << " -ignore <bc>: specify boundary to ignore" << endl
       << " -o <filename>: specify output filename" << endl ;
  exit(-1) ;
}

int main(int argc, char *argv[]) {
  Loci::Init(&argc, &argv) ;
  string query;
  vector<string> none_bcs ;
  vector<string> sym_bcs ;
  string outfile ;
  while(argc>=2 && argv[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(argc >= 3 && !strcmp(argv[1],"-q")) {
      query = argv[2] ;
      argc -= 2 ;
      argv += 2 ;
    } else if(!strcmp(argv[1],"-sym")) {
      string s = argv[2] ;
      argc -= 2 ;
      argv += 2 ;
      sym_bcs.push_back(s) ;
    } else if(!strcmp(argv[1],"-ignore")) {
      string s = argv[2] ;
      argc -= 2 ;
      argv += 2 ;
      none_bcs.push_back(s) ;
    } else if(!strcmp(argv[1],"-o")) {
      string s = argv[2] ;
      argc -= 2 ;
      argv += 2 ;
      outfile = s ;
    } else {
      cerr << "argument " << argv[1] << " is not understood." << endl ;
      argc-- ;
      argv++ ;
      Usage() ;
    }
  }
  if(argc != 2) {
    cerr << "must pass in argument for grid to process"<< endl ;
    Usage() ;
    exit(-1) ;
  }
  string base = argv[1] ;
  string meshfile = base + ".vog" ;
  if(outfile == "")
    outfile = base +".geom";
  string varfile = base +".vars";
 
  fact_db facts ;

  
  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;
  Loci::load_module("fvm",rdb) ;

  //read in grid file
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    cerr << "unable to read file " << meshfile << endl ;
    Loci::Abort() ;
  }

  store<string> boundary_names ;
  boundary_names = facts.get_variable("boundary_names") ;

  map<string,string> binfo ;
  entitySet dom = boundary_names.domain() ;
  FORALL(dom,cc) {
    binfo[boundary_names[cc]] = string("geometry") ;
  } ENDFORALL ;
  for(size_t i = 0; i<none_bcs.size();++i)
    binfo[none_bcs[i]] = string("none") ;
  for(size_t i = 0; i<sym_bcs.size();++i)
    binfo[sym_bcs[i]] = string("symmetry") ;
  
  
  try {
    string var_str ;
    var_str = "{ boundary_conditions: <" ;
    map<string,string>::const_iterator si ;
    for(si = binfo.begin();si!=binfo.end();) {
      var_str += si->first + "=" + si->second ;
      ++si ;
      if(si != binfo.end())
        var_str += "," ;
    }
    var_str += "> } " ;
    cout << "var_str="<<var_str <<endl ;
    istringstream ifile(var_str) ;

    facts.read_vars(ifile,rdb) ;
  } catch(const Loci::BasicException &err) {
    err.Print(cerr) ;
    cerr << "aborted in vars setup" << endl ;
    Loci::Abort() ;
  }


  param<std::string> outfile_par ;

  *outfile_par = outfile;
  facts.create_fact("outfile_par",outfile_par) ;

  //set up boundary
  setupBoundaryConditions(facts) ;
  
  //get boundary nodes(boundary_nodes) and boundary faces(bc1)
  constraint bc1 ;
  bc1 = facts.get_fact("geometry_BC") ;
  multiMap face2node ;
  face2node = facts.get_fact("face2node") ;
  Loci::MapRepP mp = Loci::MapRepP(face2node.Rep()) ;
  entitySet bcnodes = mp->image(*bc1) ;
  constraint boundary_nodes ;
  *boundary_nodes = bcnodes ;
  facts.create_fact("boundary_nodes",boundary_nodes) ;
 
 
  //compute and output boundary nodes    
  if(!Loci::makeQuery(rdb,facts,"bnode_output")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }

  //compute and output face2node and geometry of boundary faces
  if(!Loci::makeQuery(rdb,facts,"bface_output")) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }

  cout << "geometry written to " << outfile << endl ;
  Loci::Finalize() ;
}

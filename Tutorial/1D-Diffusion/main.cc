//*******************************************************************
//  Solution of 1D Diffusion equation using LOCI Framework:
// This was written as a part of "LOCI Beginner's Tutorial".
//
// Modified By:
// Chaman Singh Verma
// Mississippi State University.
// 22 Feb. 2001
//
//*******************************************************************
#include <Loci.h>
#include <iostream>

using namespace std ;

int main(int argc, char *argv[])
{
  Loci::Init(&argc,&argv) ;


  set_fpe_abort() ;
  
  // query for this variable by default
  std::string query = "solution" ;
  // do not write out facts by default
  bool write_facts = false;
  // do not write out rules by default
  bool write_rules = false;

  // Parse out the LOCI specific arguments.
  //	-q queryvar	make alternative queries
  //	-fact		print out the fact database during program execution
  //	-rule		print out the rule database during program execution
  for (int i=1; i<argc; i++) {
    if (!strcmp(argv[i], "-q") && (i+1) < argc) {
      query = argv[i+1];
      for (int j=0; j<argc-i-2; j++) {
	argv[i+j] = argv[i+j+2];
      }
      argc -= 2;
      i--;
    }
    else if (!strcmp(argv[i], "-fact")) {
      write_facts = !write_facts;
      for (int j=0; j<argc-i-1; j++) {
	argv[i+j] = argv[i+j+1];
      }
      argc -= 1;
      i--;
    }
    else if (!strcmp(argv[i], "-rule")) {
      write_rules = !write_rules;
      for (int j=0; j<argc-i-1; j++) {
	argv[i+j] = argv[i+j+1];
      }
      argc -= 1;
      i--;
    } else {
      std::cerr << "argument " << argv[i] << " is not understood."
		<< std::endl;
    }
  }      

  //-----------------------------------------------------------------
  // Create Rule Database
  //-----------------------------------------------------------------
  
  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;

  //-----------------------------------------------------------------
  // Create Fact Database
  //-----------------------------------------------------------------
  fact_db facts ;


  ifstream file("heat.vars",ios::in) ;
  if(ifile.fail()) {
    cerr << "can't open 'heat.vars'" << endl ;
    exit(-1) ;
  }

  param<int> Nin ;
  *Nin = 50 ;
  
  facts.read_vars(ifile,rdb) ;
  
  Nin = facts.get_variable("N") ;
  
  const int N = *Nin ; // Number of points in grid.

  //-----------------------------------------------------------------
  // Create a 1-d unstructured grid ; Allocate space for nodes and
  // cells.  Distribute to processors using a simple block partition
  //-----------------------------------------------------------------
  pair<entitySet,entitySet> node_alloc =
    facts.get_distributed_alloc(block_partition(N)) ;
  pair<entitySet,entitySet> cell_alloc = 
    facts.get_distributed_alloc(block_partition(N-1)) ;

  //-----------------------------------------------------------------
  // Generate 1D grid positions at the nodes.
  //-----------------------------------------------------------------
  store<float> x; 
  x.allocate(nodes);

  entitySet::const_iterator ei ; // Generic iterator

  for(ei=nodes.begin();ei!=nodes.end();++ei)
    x[*ei] = float(*ei)/float(N);
  //-----------------------------------------------------------------
  // Create mapping from interface to cells
  // cl = cell left , cr = cell right
  //-----------------------------------------------------------------
  
  Map cl,cr ;
  cl.allocate(nodes-interval(0,0)) ; // do not allocate for leftmost interface
  cr.allocate(nodes-interval(N,N)) ; // do not allocate for rightmost interface

  // Assign maps from nodes to cells
  // cl = {(i,l) | i \in [1,N], l = i+N}
  // cr = {(i,r) | i \in [0,N-1], r = i+N+1}
  for(ei=cl.domain().begin();ei!=cl.domain().end();++ei) 
    cl[*ei] = *ei + N;
  for(ei=cr.domain().begin();ei!=cr.domain().end();++ei) 
    cr[*ei] = *ei + N + 1;

  //-----------------------------------------------------------------
  // Create mapping from cells to interface
  // il = interface left, ir = interface right
  // il = {(c,l) | c \in cells, l = c-N-1},   
  // ir = {(c,r) | c \in cells, l = c-N}
  //-----------------------------------------------------------------
  Map il,ir ;
  il.allocate(cells) ;
  ir.allocate(cells) ;

  for(ei=cells.begin();ei!=cells.end();++ei) {
    il[*ei] = *ei - N - 1 ;
    ir[*ei] = *ei - N ;
  }


  facts.create_fact("il",il) ;
  facts.create_fact("ir",ir) ;
  facts.create_fact("x", x) ;
  facts.create_fact("cl",cl) ;
  facts.create_fact("cr",cr);

  // Diffusion constant
  param<float> nu ;
  *nu = 1.0 ;
  facts.create_fact("nu",nu) ;

  // Number of iterations to run simulation
  param<int> max_iteration ;
  *max_iteration = 100 ;
  facts.create_fact("max_iteration",max_iteration) ;

  // Minimum L1 norm for convergence test
  param<double> error_tolerance;
  *error_tolerance = 1.0E-03;
  facts.create_fact( "error_tolerance", error_tolerance);

  // Identify boundary conditions
  constraint left_boundary ;
  constraint right_boundary ;
  *right_boundary = cl.domain() - cr.domain() ;
  *left_boundary = cr.domain() - cl.domain() ;

  facts.create_fact("left_boundary",left_boundary) ;
  facts.create_fact("right_boundary",right_boundary) ;

  if (Loci::MPI_processes == 1) {
    // Write out the initial fact database if -fact
    if (write_facts) {
      char fn_facts[11];
      strcpy(fn_facts, "facts_init");
      std::ofstream of_facts(fn_facts);

      facts.write(of_facts);

      of_facts.close();
    }

    // Write out the rule database if -rule
    if (write_rules) {
      char fn_rules[6];
      strcpy(fn_rules, "rules");
      std::ofstream of_rules(fn_rules);

      Loci::ruleSet rset = rdb.all_rules();
      Loci::ruleSet::const_iterator itr;
      for (itr = rset.begin(); itr!=rset.end(); itr++) {
	of_rules << *itr << std::endl;
      }

      of_rules.close();
    }
  }

  // Query Loci for fact derived fact 'solution'
  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  
  //-----------------------------------------------------------------
  // Final Step: Query the database for solution:
  //-----------------------------------------------------------------

  if(query == "solution") {
    store<float> usol ;
    usol = facts.get_variable("solution") ;

    cout << "The solution is : " <<endl;
    for(ei=cells.begin();ei!=cells.end();++ei) 
      cout << ""<< *ei<<" "<<usol[*ei]<<endl ;
  }
  //-----------------------------------------------------------------
  // End of computations::
  //-----------------------------------------------------------------

  Loci::Finalize() ;
  return 0;

}
//*******************************************************************

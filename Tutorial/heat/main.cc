#include <Loci.h>
#include <iostream>

using namespace std ;


int main(int argc, char *argv[])
{
  // Initialize Loci
  Loci::Init(&argc,&argv) ;

  // Setup exceptions so program aborts on floating point exceptions
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

  // Load module of finite-volume helper rules
  Loci::load_module("fvm",rdb) ;

  //-----------------------------------------------------------------
  // Create Fact Database
  //-----------------------------------------------------------------
  fact_db facts ;


  char *filename = "test" ;
  // Read in the vars file
  char buf[512] ;
  sprintf(buf,"%s.vars",filename) ;
  ifstream ifile(buf,ios::in) ;
  if(ifile.fail()) {
    cerr<<"can't open " << buf << endl ;
    exit(-1) ;
  }
  facts.read_vars(ifile,rdb) ;
  ifile.close() ;

  // Read in the XDR file
  string file = string(filename) + string(".xdr") ;
  if(!Loci::setupFVMGrid(facts,file)) {
    cerr << "unable to read grid file '" << file << "'" << endl ;
    Loci::Abort() ;
  }

  setupBoundaryConditions(facts) ;

  createLowerUpper(facts) ;

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
    entitySet::const_iterator ei ;
    entitySet dom = usol.domain() ;
    for(ei=dom.begin();ei!=dom.end();++ei) 
      cout << ""<< *ei<<" "<<usol[*ei]<<endl ;
  } else {
    // If we queried for something else, print out the results
    using Loci::variableSet ;
    variableSet queries = variableSet(Loci::expression::create(query)) ;
    for(variableSet::const_iterator vi=queries.begin();vi!=queries.end();++vi){
      Loci::storeRepP rep = facts.get_variable(*vi) ;
      cout << "variable " << *vi << "= " << endl ;
      rep->Print(cout) ;
    }
  }
   
  //-----------------------------------------------------------------
  // End of computations::
  //-----------------------------------------------------------------

  Loci::Finalize() ;
  return 0;
}
  

//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
//#include <fstream>
#include <iostream>
#include <string>
//#include <utility>
//#include <vector>
//#include <list>
//#include <algorithm>
#include <Loci.h>
//#include "hexcell.h"
#include "defines.h"
//#include <rpc/xdr.h>
//#include <rpc/rpc.h>
using std::string;
using std::cout;
using std::endl;
using std::cerr;
//using std::ofstream;
//using std::ifstream;
using std::ios;

namespace Loci{
  void parallelClassifyCell(fact_db &facts) ;
   void createEdgesParallel(fact_db &facts);
}


int main(int argc, char ** argv) {
  // Let Loci initialize itself.
  // This also gives Loci first dibs on the command line arguments.
  Loci::Init(&argc, &argv);
  

  // This is the name of the mesh file that you want to read in.
  // This may be overridden by the command line argument "-g file.xdr"
  //  string meshfile = "testGrid.xdr";
  string meshfile = "testGrid.vog";
  //This is the name of the refinement plan file
  string planFile = "out.plan";
  //This is the name of the output gridfile
  string outFile  = "out.xdr";
  
  // Here's where we parse out the command line arguments that are
  // relevant to this program.
  int j=1;
  //string pathname = "../gridgen/restore/mixedcell/";
  // string pathname = "/var/tmp/qxue/grid/";
   string pathname = "";
   //  bool restart = false;
   // bool xdr = true;


  
    //print out help info
  if( (argc == 1)||(argc==2) ){
    if(Loci::MPI_rank == 0){
      cout<<"command line:" << endl;
      cout <<"refmesh <options> <filename> <options> <filename> ... "<< endl;
      cout << endl;
      cout << "options:" << endl;
      cout <<"-g <file> -- original grid file, refinement plans are based on this grid" << endl;
      cout <<"-r <file> -- input refinement plan file" <<endl;
      cout <<"-o <file> -- output grid file" << endl;
     
    }
    return 0;
  }
  
  for (int i=1; i<argc; i++) {
    // Let's look at the i'th argument
    string arg(argv[i]);
    
    if (arg == "-g" && (i+1) < argc) {
      // Replace the mesh filename with the next argument
      meshfile =  argv[++i];      
    }
    else if(arg == "-o" && (i+1) < argc){
       //replace the output filename with the next argument
      outFile =  argv[++i];
    }
    
    else if(arg == "-r" && (i+1) < argc){
      //replace the input refinement plan filename with the next argument
      planFile =  argv[++i];
      
    }
     
    
    else{
      // Anything that we don't recognize gets recycled back into argv
      argv[j++] = argv[i];
    }
  }
  // Set the new size of the argument list to however many arguments
  // we put back in argv.
  argc = j;
  meshfile = pathname + meshfile;
  outFile = pathname + outFile;
  planFile = pathname + planFile;
  



 // Setup the fact database.
  fact_db facts;
  

  
  
  param<std::string> planfile_par ;
  *planfile_par = planFile;
  facts.create_fact("planfile_par",planfile_par) ;

  param<std::string> meshfile_par ;
  *meshfile_par = meshfile;
  facts.create_fact("meshfile_par",meshfile_par) ;

  param<std::string> outfile_par ;
  *outfile_par = outFile;
  facts.create_fact("outfile_par",outfile_par) ;

  
  // Setup the rule database.
  // Add all registered rules.
  rule_db rules;
  rules.add_rules(global_rule_list);
  
 
 if(Loci::MPI_rank == 0) std::cout <<"reading in meshfile" << std::endl;
  // Read in the mesh file.  Setup Loci datastructures
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    std::cerr << "unable to read grid file '" << meshfile << "'" << std::endl ;
    Loci::Abort() ;
  }
  
  Loci::createLowerUpper(facts) ;
  Loci::createEdgesParallel(facts) ;
 
  Loci:: parallelClassifyCell(facts);
  

  
   if(!Loci::makeQuery(rules, facts, "node_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
   }
   
   if(!Loci::makeQuery(rules, facts, "face_output")) {
     std::cerr << "query failed!" << std::endl;
     Loci::Abort();
   }
 

 
  
  
  // Tell Loci to cleanup after itself and exit gracefully.
  Loci::Finalize();
  
}

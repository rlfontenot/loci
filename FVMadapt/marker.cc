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
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                        test8.cc
//
//   this file read in a coarse grid file, a refinement plan file(based on the coarse grid),
//and a tag file(if input refinement file is provided, the tag file is for new grid generated using the input
//refinement plan, otherwise, the tag file is for the coarse grid file),
//and then write out a new  refinement plan file(the plan also based on the coarse grid).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <string>

#include <Loci.h>
#include <stdlib.h>

#include "globals.h"
#include <Loci>

using std::cout ;


using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using Loci::MPI_processes;
namespace Loci{
  void parallelClassifyCell(fact_db &facts) ;
   void createEdgesParallel(fact_db &facts) ;
}

  
int main(int argc, char ** argv) {
  // Let Loci initialize itself.
  // This also gives Loci first dibs on the command line arguments.
  Loci::Init(&argc, &argv);
  
  
  // This is the name of the mesh file that you want to read in.
  // This may be overridden by the command line argument "-g file.xdr"
  //  string meshfile = "testGrid.xdr";
  string meshfile = "testGrid.vog";
  //This is the name of the input refinement plan file
  string planFile = "out1.plan"; //default file can be empty file
  //This is the name of the output refinement file
  string outFile  = "out.plan";
  //this is the name of the tag file
  string tagFile = "r1.tag";
  //this is the name of xml file
  string xmlFile = "region.xml";
  // Here's where we parse out the command line arguments that are
  // relevant to this program.
  int j=1;
  //  string pathname = "/scratch/qxue/output/";
  // string pathname = "/var/tmp/qxue/grid/";
  string pathname = "";
  //if the input refinement filename is provided in command line, restart is true,
  //otherwise, restart is false
  bool restart = false;
  //if xml_input is true, cell plan is defines by a region and tol value. otherwise,
  //it's decided by tagfile
  bool xml_input = false;
  
  int split_mode = 0;//default mode, hexcells and prisms split according to edge length
  //print out help info
  if( (argc == 1)||(argc==2) ){
  
    if(Loci::MPI_rank == 0){
      cout<<"Usuage:" << endl;
      cout <<"marker <options> <filename/input value> <options> <filename/input value> ... "<< endl;
      cout << endl;
      cout << "options:" << endl;
      cout <<"-g <file> -- original grid file(required)," << endl;
      cout <<"             refinement plans are based on this grid" << endl;
      cout <<"-r <file> -- input refinement plan file(optional)" <<endl;
      cout << "-xml <file> -- input xml file(optional), the file define a geometric region" << endl;
      cout << "            --all cells inside the region will be refined to a tolerance value" << endl;
      cout << "            --xml option and -tag option can not be selected at the same time" <<endl;
      
      cout <<"-tag <file> -- input tag file(optional), " << endl;
      cout <<"               if there is an input refinement plan file,"<<endl;
      cout <<"               the tag file is for the refined grid" << endl;
      cout <<"               otherwise, the tag file is for the original grid" << endl;
      cout <<"-o <file> -- output refinement plan file(default)" << endl;
      cout <<"-tol <double> -- tolerance, minimum grid spacing allowed(default)" << endl;
      cout <<"-fold <double> -- twist value, maximum face folding allowed(default)" << endl;
      cout <<"-levels <int> --  levels of refinement(default)" << endl;
      cout << "-mode <int: 0, 1, 2, 3> -- split mode 0: anisotropic split according to edge length" << endl;
      cout << "                        -- split mode 1: don't refine in z direction" << endl;
      cout << "                        -- split mode 2: fully isotropic split" << endl;
    
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
      //replace the input filename with the next argument
      planFile =  argv[++i];
      restart = true;
    }
    
    else if(arg == "-tag" && (i+1) < argc){
      //replace the tag filename with the next argument
      tagFile =  argv[++i];
    }
     else if(arg == "-xml" && (i+1) < argc){
       //replace the xml filename with the next argument
       xmlFile =  argv[++i];
       xml_input = true;
     }

     
     else if(arg == "-tol" && (i+1) < argc){
       //replace the tolarence(minimum grid spacing) with the next argument
      char** endptr = 0;
      Globals::tolerance = strtod(argv[++i], endptr);
      
    }
    
    else if(arg == "-fold" && (i+1) < argc){
      //replace the fold(maximum face fold value)with the next argument
      char** endptr = 0;
      Globals::fold = strtod(argv[++i], endptr);
      
    }


     else if(arg == "-mode" && (i+1) < argc){
       //replace the fold(maximum face fold value)with the next argument
       split_mode = atoi(argv[++i]);
      
     }
       
     else if(arg == "-levels" && (i+1) < argc){
       //replace the fold(maximum face fold value)with the next argument
       Globals::levels = atoi(argv[++i]);
       
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
  tagFile = pathname + tagFile;
  xmlFile = pathname + xmlFile;
  
  if(Loci::MPI_rank == 0){
    cout <<"Marker running" <<endl;
    cout << "fold: " << Globals::fold << endl;
    cout << "tolerance: "<<Globals::tolerance <<endl;
    cout << "levels: " <<Globals::levels << endl;
    cout << "restart: " << restart << endl;
    cout << "split_mode: " << split_mode << endl;
        
  }

  if(split_mode == 0 && Globals::levels > 1){
    cerr << "WARNING: multi-level refinement is not allowed in anisotropic refinement"<< endl;
    Loci::Abort();
  }

  if(xml_input && Globals::tolerance == 0.0){
    cerr <<"WARNING: zero tolerance will cause infinite loop at -xml option"<< endl;
    Loci::Abort();
  }
  

  


  // Setup the rule database.
  // Add all registered rules.  
  rule_db rules;
  rules.add_rules(global_rule_list);
  
  // Setup the fact database.
  fact_db facts;
  
  if(Loci::MPI_rank == 0) cout <<"reading in meshfile" << std::endl;

  // Read in the mesh file.  
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    std::cerr << "unable to read grid file '" << meshfile << "'" << std::endl ;
    Loci::Abort() ;
  }
  
  //Setup Loci datastructures
  Loci::createLowerUpper(facts) ;
  Loci::createEdgesParallel(facts) ;
  Loci:: parallelClassifyCell(facts);
 
  
 
  
  param<std::string> tagfile_par ;
  
  *tagfile_par = tagFile;
  facts.create_fact("tagfile_par",tagfile_par) ;
  
  param<std::string> outfile_par ;
  *outfile_par = outFile;
  facts.create_fact("outfile_par",outfile_par) ;
  
  
 

  if(restart && xml_input){
    param<std::string> planfile_par ;
    *planfile_par = planFile;
    facts.create_fact("planfile_par",planfile_par) ;
    
    param<int> restart_xml_par;
    *restart_xml_par = 1;
    facts.create_fact("restart_xml_par",restart_xml_par);
  }
  else if(restart && (!xml_input)){
    param<std::string> planfile_par ;
    *planfile_par = planFile;
    facts.create_fact("planfile_par",planfile_par) ;
    
    param<int> restart_no_xml_par;
    *restart_no_xml_par = 1;
    facts.create_fact("restart_no_xml_par",restart_no_xml_par);
  }

  else if((!restart) && xml_input){
    param<std::string> xmlfile_par ;
    *xmlfile_par = xmlFile;
    facts.create_fact("xmlfile_par",xmlfile_par) ;
    
    param<int> no_restart_xml_par;
    *no_restart_xml_par = 1;
    facts.create_fact("no_restart_xml_par",no_restart_xml_par);
  }
  else if((!restart) && (!xml_input)){
    param<std::string> xmlfile_par ;
    *xmlfile_par = xmlFile;
    facts.create_fact("xmlfile_par",xmlfile_par) ;
    
    param<int> no_restart_no_xml_par;
    *no_restart_no_xml_par = 1;
    facts.create_fact("no_restart_no_xml_par",no_restart_no_xml_par);
  }

  
  param<int> split_mode_par;
  *split_mode_par = split_mode;
  facts.create_fact("split_mode_par", split_mode_par);

  
  
      if(!Loci::makeQuery(rules, facts, "cellplan_output")) {
  //  if(!Loci::makeQuery(rules, facts, "balancedCellPlan")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }
  
  // Tell Loci to cleanup after itself and exit gracefully.
  Loci::Finalize();
}
  

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

int main()
{
  const int N = 50 ; // Number of points in grid.

  //-----------------------------------------------------------------
  // Create a 1-d unstructured grid ; Node and Cells are identified
  // separately.
  //-----------------------------------------------------------------
  entitySet nodes  = interval(0,N) ;
  entitySet cells  = interval(N+1,2*N);

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

  //-----------------------------------------------------------------
  // Create Facts Database
  //-----------------------------------------------------------------
  fact_db facts ;

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

  //-----------------------------------------------------------------
  // Create Rule database ...
  //-----------------------------------------------------------------
  
  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;

  //-----------------------------------------------------------------
  // Create and execute the schedule to  obtains the solution
  //-----------------------------------------------------------------
  executeP schedule = create_execution_schedule(rdb,facts,"solution") ;

  //schedule->Print(cout) ;                   // Display schedule

  schedule->execute(facts) ;           // Execute the schedule

  //-----------------------------------------------------------------
  // Final Step: Query the database for solution:
  //-----------------------------------------------------------------

  store<float> usol ;
  usol = facts.get_variable("solution") ;

  cout << "The solution is : " <<endl;
  for(ei=cells.begin();ei!=cells.end();++ei) 
    cout << ""<< *ei<<" "<<usol[*ei]<<endl ;

  //-----------------------------------------------------------------
  // End of computations::
  //-----------------------------------------------------------------

  return(0);

}
//*******************************************************************

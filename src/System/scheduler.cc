#include "sched_tools.h"
#include "distribute.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

#define PROFILE_CODE

namespace Loci {

  double get_timer() {
#ifdef PROFILE_CODE
    clock_t tc ;
    static double to = 0;
    double tn,t ;
    tc = clock() ;
    tn = tc/1000000.0 ;
    t = tn - to ;
    to = tn ;
    return t ;
#else
    return -1.0 ;
#endif
  }


  
  executeP create_execution_schedule(rule_db &rdb,
                                     fact_db &facts,
                                     std::string target_string,
                                     int nth) {
    num_threads = min(nth,max_threads) ;

    //double timer = get_timer() ;
    
    variableSet given = facts.get_typed_variables() ;
    variableSet target(expression::create(target_string)) ;

    cout << "generating dependency graph..." << endl ;
    double start_time = MPI_Wtime() ;
    digraph gr = dependency_graph(rdb,given,target).get_graph() ;
    // If graph is empty, return a null schedule 
    if(gr.get_target_vertices() == EMPTY)
      return executeP(0) ;

    cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr) ;
    
    cout << "decomposing graph..." << endl ;

    
    decomposed_graph decomp(gr,given,target) ;
    graph_compiler compile_graph(decomp) ;
    double end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for graph processing  = " << end_time  - start_time << "  seconds " << endl ;
#ifdef PROFILE_CODE
    //timer = get_timer() ;
    //cout << "Graph Processing Time: "<<timer << " seconds" << endl ;
#endif
    
    cout << "existential analysis..." << endl ;
    start_time = MPI_Wtime() ;
    compile_graph.existential_analysis(facts) ;
    end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for existential_analysis  = " << end_time  - start_time << "  seconds " << endl ;
    cout << "creating execution schedule..." << endl;
    start_time = MPI_Wtime() ;
    executeP sched =  compile_graph.execution_schedule(facts,num_threads) ;
    end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for schedule generation  = " << end_time  - start_time << "  seconds " << endl ;
#ifdef PROFILE_CODE    
    //timer = get_timer() ;
    //cout << "Schedule Generation Time: " << timer << " seconds" << endl ;
#endif
    return sched ;
  }
}

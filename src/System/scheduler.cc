#include "sched_tools.h"
#include "distribute.h"
#include "param_rule.h"
#include <Tools/stream.h>
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ofstream ;

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
    sched_db scheds(facts) ;
    variableSet given = facts.get_typed_variables() ;
    variableSet target(expression::create(target_string)) ;
    if(Loci::MPI_rank==0)
      cout << "generating dependency graph..." << endl ;
    double start_time = MPI_Wtime() ;
    rule_db par_rdb ;
    par_rdb = parametric_rdb(rdb,target) ;
    digraph gr = dependency_graph(par_rdb,given,target).get_graph() ;
    // If graph is empty, return a null schedule 
    if(gr.get_target_vertices() == EMPTY)
      return executeP(0) ;
    
    if(Loci::MPI_rank==0)
      cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr, scheds) ;
    if(Loci::MPI_rank==0)
      cout << "decomposing graph..." << endl ;
    decomposed_graph decomp(gr,given,target) ;
    variableSet fact_vars, initial_vars ;
    fact_vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      if(variable(*vi).time().level_name() == "*" ) {
	storeRepP vp = facts.get_variable(*vi) ;
	if(vp->RepType() == STORE) {
	  ostringstream oss ;
	  oss << "source(" <<"EMPTY"<<')' ;
	  oss << ",target(" << *vi << ')' ;
	  string sig = oss.str() ;
	  rule r(sig) ;
	  entitySet t ;
	  if(par_rdb.rules_by_target(*vi) == EMPTY) {
	    if(facts.isDistributed()) {
	      fact_db::distribute_infoP d = facts.get_distribute_info() ;
	      for(int i = 0; i < d->copy.size(); ++i)
		t += d->copy[i].entities ;
	      initial_vars += *vi ;
	      scheds.set_existential_info(*vi, r, t) ;
	    }
	  }
	}
      }
    }
    if(Loci::MPI_rank==0)
      cout << " initial_vars = " << initial_vars << endl ;
    graph_compiler compile_graph(decomp, initial_vars) ;
    double end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for graph processing  = "
                   << end_time  - start_time << "  seconds " << endl ;
#ifdef PROFILE_CODE
    //timer = get_timer() ;
    //cout << "Graph Processing Time: "<<timer << " seconds" << endl ;
#endif
    if(Loci::MPI_rank==0)
      cout << "existential analysis..." << endl ;
    start_time = MPI_Wtime() ;
    compile_graph.existential_analysis(facts, scheds) ;
    end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for existential_analysis  = "
                   << end_time  - start_time << "  seconds " << endl ;
    if(Loci::MPI_rank==0)
      cout << "creating execution schedule..." << endl;
    start_time = MPI_Wtime() ;
    executeP sched =  compile_graph.execution_schedule(facts,scheds, num_threads) ;
    end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for schedule generation  = " << end_time  - start_time << "  seconds " << endl ;

    //scheds.print_summary(Loci::debugout) ;
#ifdef PROFILE_CODE    
    //timer = get_timer() ;
    //cout << "Schedule Generation Time: " << timer << " seconds" << endl ;
#endif
    return sched ;
  }
}

#include "comp_tools.h"

// flag to use the new memory scheme
// disable it and turn on the USE_ALLOCATE_ALL_VARS (below, sched_comp.cc)
// to use the original memory allocation
//#define USE_MEMORY_SCHEDULE
#define USE_ALLOCATE_ALL_VARS

namespace Loci {
  void promote_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
#ifdef USE_MEMORY_SCHEDULE
    /////////////////////////////////////////////////////////////////////
    // set_var_existence is from given to the goal, this is
    // executed first than the process_var_requests.
    // in this phase, disable the de-allocation of all the rule's sources
    // the de-allocation compilers will remove any variables from being
    // de-allocated set here

    scheds.add_disabled(r.sources()) ;
    
    /////////////////////////////////////////////////////////////////////
#endif
  }

  void promote_compiler::process_var_requests(fact_db &facts, sched_db &scheds)
  {
#ifdef USE_MEMORY_SCHEDULE
    /////////////////////////////////////////////////////////////////////
    // process_var_requests is from the goal to the given, this is
    // executed after than the set_var_requests.
    // in this phase, disable the allocation of all the rule's targets
    // the allocation compilers will remove any variables from being
    // allocated set here

    scheds.add_disabled(r.targets()) ;

    ////////////////////////////////////////////////////////////////////
#endif
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef USE_ALLOCATE_ALL_VARS
    return executeP(0) ;
#endif
#ifdef USE_MEMORY_SCHEDULE
    std::ostringstream oss ;
    oss << r ;
    return executeP(new execute_msg(oss.str())) ;
#endif
  }

  void generalize_compiler::set_var_existence(fact_db &facts, sched_db &scheds)
  {
#ifdef USE_MEMORY_SCHEDULE
    ////////////////////////////////////////////////////////////////////
    // the same as promote_compiler

    scheds.add_disabled(r.sources()) ;

    ////////////////////////////////////////////////////////////////////
#endif
  }

  void generalize_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
#ifdef USE_MEMORY_SCHEDULE
    ////////////////////////////////////////////////////////////////////
    // the same as promote_compiler

    scheds.add_disabled(r.targets()) ;

    ////////////////////////////////////////////////////////////////////
#endif
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef USE_ALLOCATE_ALL_VARS
    return executeP(0) ;
#endif
#ifdef USE_MEMORY_SCHEDULE
    std::ostringstream oss ;
    oss << r ;
    return executeP(new execute_msg(oss.str())) ;
#endif
  }

  void priority_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void priority_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP priority_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    return executeP(0) ;
  }
}

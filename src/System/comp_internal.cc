#include "comp_tools.h"

// flags to enable(disable) writing out the
// these internals to the schedule file
#define WRITE_OUT
//#define DONT_WRITE_OUT

namespace Loci {
  void promote_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void promote_compiler::process_var_requests(fact_db &facts, sched_db &scheds)
  {
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
    std::ostringstream oss ;
    oss << r ;
    return executeP(new execute_msg(oss.str())) ;
#endif
  }

  void generalize_compiler::set_var_existence(fact_db &facts, sched_db &scheds)
  {
  }

  void generalize_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
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
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
    std::ostringstream oss ;
    oss << r ;
    return executeP(new execute_msg(oss.str())) ;
#endif
  }
}

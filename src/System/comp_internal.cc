#include "comp_tools.h"

namespace Loci {
  void promote_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void promote_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    return executeP(0) ;
  }

  void generalize_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void generalize_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    return executeP(0) ;
  }

  void priority_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void priority_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP priority_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    return executeP(0) ;
  }
}

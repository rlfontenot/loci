#include "comp_tools.h"

namespace Loci {
  void promote_compiler::set_var_existence(fact_db &facts) {
  }

  void promote_compiler::process_var_requests(fact_db &facts) {
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }

  void generalize_compiler::set_var_existence(fact_db &facts) {
  }

  void generalize_compiler::process_var_requests(fact_db &facts) {
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }

  void priority_compiler::set_var_existence(fact_db &facts) {
  }

  void priority_compiler::process_var_requests(fact_db &facts) {
  }

  executeP priority_compiler::create_execution_schedule(fact_db &facts) {
    return executeP(0) ;
  }
}

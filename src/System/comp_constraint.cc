#include "comp_tools.h"
#include "dist_tools.h"
using std::ostream ;
using std::endl ;

namespace Loci {

  extern int current_rule_id ;

  class execute_constraint_rule: public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    bool do_run ;
  public:
    execute_constraint_rule(rule fi, sequence seq,
                            fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  execute_constraint_rule::
  execute_constraint_rule(rule fi, sequence seq,
                          fact_db& facts, sched_db& scheds) {
    do_run = true ;
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  void execute_constraint_rule::execute(fact_db &facts) {
    current_rule_id = rule_tag.ident() ;
    if(do_run) {
      rp->compute(exec_seq) ;
    }
  }
  
  void execute_constraint_rule::Print(ostream &s) const {
    s << rule_tag << "  over sequence " << exec_seq << endl ;
  }
  
  void constraint_compiler::set_var_existence(fact_db& facts,
                                              sched_db& scheds)
  {}
  
  void constraint_compiler::process_var_requests(fact_db& facts,
                                                 sched_db& scheds) {
    variableSet var_requests = constraint_rule.targets() ;
    variableSet::const_iterator vi ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      scheds.variable_request(*vi,~EMPTY) ;
    }
    
    process_rule_requests(constraint_rule,facts, scheds) ;
  }
  
  executeP constraint_compiler::create_execution_schedule(fact_db& facts,
                                                          sched_db& scheds) {
    return new execute_constraint_rule(constraint_rule,
                                       ~EMPTY, facts,scheds) ;
  }

}

#include "sched_tools.h"

namespace Loci {
  dynamic_schedule_rule::dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_set = eset ;
  }


  void dynamic_schedule_rule::execute(fact_db &facts) {
    rp->compute(sequence(exec_set)) ;
  }

  void dynamic_schedule_rule::Print(std::ostream &s) const {
    s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl ;
  }
}

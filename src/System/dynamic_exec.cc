#include "sched_tools.h"
using std::cout ;
using std::cerr ;
using std::endl ;

namespace Loci {
  dynamic_schedule_rule::dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    local_compute = rp->new_rule_impl() ;
    entitySet in = rule_tag.sources() ;
    outputs = rule_tag.targets() ;
    //    cout << "inputs = " << inputs << endl ;
    for(variableSet::const_iterator vi=in.begin();vi!=in.end();++vi) {
      storeRepP store_ptr = rp->get_store(*vi) ;
      if((store_ptr != 0) && store_ptr->RepType() == Loci::STORE) {
        inputs += *vi ;
        local_facts.create_fact(*vi,store_ptr->new_store(EMPTY)) ;
      } else {
        local_facts.create_fact(*vi,facts.get_variable(*vi)) ;
      }
    }
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP store_ptr = rp->get_store(*vi)->new_store(EMPTY) ;
      local_facts.create_fact(*vi,store_ptr) ;
    }
    
    local_compute->initialize(local_facts) ;
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

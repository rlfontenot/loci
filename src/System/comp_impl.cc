#include "comp_tools.h"
#include "dist_tools.h"
#include "loci_globs.h"
using std::ostream ;
using std::endl ;

namespace Loci {

  int current_rule_id = 0 ;
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
   
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p, const sched_db &scheds) {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  void execute_rule::execute(fact_db &facts) {
    current_rule_id = rule_tag.ident() ;
    rp->compute(exec_seq) ;
  }
  
  void execute_rule::Print(ostream &s) const {
    printIndent(s) ;
    s << rule_tag << " over sequence " ;
    if(verbose || exec_seq.num_intervals() < 4) {
      s << exec_seq << endl ;
    } else {
      s << "[ ... ]" << endl ;
    }
  }
  
  void impl_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(impl,facts, scheds) ;
  }
  
  void impl_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    exec_seq = process_rule_requests(impl,facts, scheds) ;
  }
  
  executeP impl_compiler::create_execution_schedule(fact_db &facts,sched_db &scheds ) {
    //    if(GLOBAL_AND(exec_seq.size()==0)) {
    //      return executeP(0) ;
    //    }
    variableSet targets = impl.targets() ;
    WARN(targets.size() == 0) ;
    if(num_threads > 1 &&
       exec_seq.size() > num_threads*30 &&
       !impl.get_info().output_is_parameter &&
       (impl.get_info().rule_impl->get_rule_class()==rule_impl::POINTWISE||
        impl.get_info().rule_impl->get_rule_class()==rule_impl::UNIT)&&
       impl.get_info().rule_impl->thread_rule() &&
       (targets.begin()->get_info()).name != "OUTPUT") {
      execute_par *ep = new execute_par ;
      parallel_schedule(ep,exec_seq,impl,facts, scheds) ;
      return ep ;
    }
    if((targets.begin()->get_info()).name == "OUTPUT") {
      CPTR<execute_list> el = new execute_list ;
      execution_factory ef(impl,sequence(exec_seq), facts, scheds);
      el->append_list(ef.create_product());
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
      return executeP(el) ;
    }
    if(impl.get_info().rule_impl->dynamic_schedule_rule() && use_dynamic_scheduling) 
      return new dynamic_schedule_rule(impl,exec_seq,facts, scheds) ;
    else {
      execution_factory ef(impl,sequence(exec_seq),facts, scheds);
      return ef.create_product();
    }
  }

  // blackbox_compiler code
  void
  blackbox_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    // set UNIVERSE existence for all targets
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.set_existential_info(*vi, impl, ~EMPTY) ;
  }

  void
  blackbox_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    // everyone will need to request for their existence
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;

    variableSet sources = impl.sources() ;
    for(variableSet::const_iterator vi=sources.begin();
        vi!=sources.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
  }

  executeP
  blackbox_compiler::create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) {
    return new execute_rule(impl, ~EMPTY, facts, scheds);
  }

}
 
 

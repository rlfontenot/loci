#include "comp_tools.h"
#include "dist_tools.h"
using std::ostream ;
using std::endl ;

namespace Loci {

  int current_rule_id = 0 ;
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts, sched_db &scheds)  {
    do_run = true ;
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p, sched_db &scheds)
  {
    do_run = true ;
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  execute_rule::execute_rule(bool output_empty, rule fi, sequence seq, fact_db &facts, sched_db &scheds)
  {
    do_run = false ;
    if(output_empty) {
      do_run = true ;
      rp = fi.get_rule_implP() ;
      rule_tag = fi ;
      rp->initialize(facts) ;
      exec_seq = seq ;
      control_thread = false ; 
    }
  }
  
  void execute_rule::execute(fact_db &facts) {
    current_rule_id = rule_tag.ident() ;
    if(do_run) {
      rp->compute(exec_seq) ;
    }
  }
  
  void execute_rule::Print(ostream &s) const {
    s << rule_tag << "  over sequence " << exec_seq << endl ;
  }
  
  void impl_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(impl,facts, scheds) ;
  }
  
  void impl_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    exec_seq = process_rule_requests(impl,facts, scheds) ;
  }
  
  executeP impl_compiler::create_execution_schedule(fact_db &facts,sched_db &scheds ) {
#ifndef DEBUG
    //if(exec_seq.size() == 0)
    //return executeP(0) ;
#endif
    if(GLOBAL_AND(exec_seq.size()==0))
      return executeP(0) ;
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
      if(exec_seq == EMPTY) {
	bool output_empty = true ;
	el->append_list(new execute_rule(output_empty,impl,sequence(exec_seq),facts, scheds)) ;
      }
      else
	el->append_list(new execute_rule(impl,sequence(exec_seq),facts, scheds)) ;
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
      return executeP(el) ;
    }
    if(impl.get_info().rule_impl->dynamic_schedule_rule() && use_dynamic_scheduling) 
      return new dynamic_schedule_rule(impl,exec_seq,facts, scheds) ;
    else 
      return new execute_rule(impl,sequence(exec_seq),facts, scheds) ;
  }

}
 
 

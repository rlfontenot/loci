#include "comp_tools.h"
#include <distribute.h>
namespace Loci {
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts)  {
    do_run = true ;
    if(seq.num_intervals() == 0)
      do_run = false ;
    
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p)
  {
    do_run = true ;
    if(seq.num_intervals() == 0)
      do_run = false ;
    
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  void execute_rule::execute(fact_db &facts) {
    if(do_run)
      rp->compute(exec_seq) ;
  }
  
  void execute_rule::Print(ostream &s) const {
    if(MPI_processes > 1)
      s <<  "processor : " << MPI_rank  << "\t" ;
    s << rule_tag << " over sequence " << exec_seq << endl ;
  }
  
  void impl_compiler::set_var_existence(fact_db &facts) {
    existential_rule_analysis(impl,facts) ;
  }
  
  void impl_compiler::process_var_requests(fact_db &facts) {
    exec_seq = process_rule_requests(impl,facts) ;
  }
  
  executeP impl_compiler::create_execution_schedule(fact_db &facts) {
    
#ifndef DEBUG
    if(exec_seq.size() == 0)
      return executeP(0) ;
#endif
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
      parallel_schedule(ep,exec_seq,impl,facts) ;
      return ep ;
    }
    if((targets.begin()->get_info()).name == "OUTPUT") {
      CPTR<execute_list> el = new execute_list ;
      el->append_list(new execute_rule(impl,sequence(exec_seq),facts)) ;
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
      return executeP(el) ;
    }
    return new execute_rule(impl,sequence(exec_seq),facts) ;
  }

}
 
 

#include "sched_tools.h"
#include "loci_globs.h"
#include <distribute.h>

namespace Loci {
  execution_factory::execution_factory(rule fi, sequence seq, fact_db &ft, 
				       const sched_db &sd)
    :scheds(sd), facts(ft){
    rule_tag = fi ;
    exec_seq = seq ; 
  } 
  
  execute_modules* execution_factory::create_product() {
    if((rule_tag.get_info().output_is_parameter ||
        rule_tag.get_info().rule_impl->thread_rule())) {
      if(GLOBAL_AND(exec_seq.size() == 0)) {
        return new execute_rule_null(rule_tag) ;
      }
    } else if(exec_seq.size() == 0) {
      return new execute_rule_null(rule_tag) ;
    }
      
    if(Loci::collect_timings)
      return new timed_execute_rule(rule_tag, exec_seq, facts, scheds,
				    Loci::time_duration_to_collect_data);
    else
      return new execute_rule(rule_tag, exec_seq, facts, scheds);
  }

  execute_modules* execution_factory::create_product(variable v, const storeRepP &p) {
    //    if(GLOBAL_AND(exec_seq.size() ==0)) {
    //      return new execute_rule_null(rule_tag) ;
    //    }
    if(Loci::collect_timings)
      return new timed_execute_rule(rule_tag, exec_seq, facts, v, p, scheds,
				    Loci::time_duration_to_collect_data);
    else
      return new execute_rule(rule_tag, exec_seq, facts, v, p, scheds);
  }
 
  timed_execute_rule::timed_execute_rule(rule fi, sequence seq, fact_db &facts,
		     const sched_db &scheds, double td)
    :execute_rule(fi, seq, facts, scheds), min_time_duration(td) {
  }
  
  timed_execute_rule::timed_execute_rule(rule fi, sequence seq, fact_db &facts, variable v,
		     const storeRepP &p, const sched_db &scheds, double td)
    :execute_rule(fi, seq, facts, v, p, scheds), min_time_duration(td){
  }
  
  void timed_execute_rule::execute(fact_db &facts) {
    int count = 0;
    double et;
    double st = MPI_Wtime();
    sequence empty_seq;
    do{
      rp->compute(empty_seq);
      et = MPI_Wtime();
      count++;
    }while(et - st < min_time_duration); 
    
    timeout << rule_tag <<"\t" << count << "\t" << et-st << "\t";
    
    count = 0;
    st = MPI_Wtime();
    do{
      rp->compute(exec_seq);
      et = MPI_Wtime();
      count++;
    }while(et - st < min_time_duration);
      
    timeout << exec_seq.size() <<"\t" << count << "\t" << et-st << endl;
  }
}

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
        !rule_tag.get_info().rule_impl->thread_rule())) {
      if(GLOBAL_AND(exec_seq.size() == 0)) {
        return new execute_rule_null(rule_tag) ;
      }
    }
    
    if(Loci::collect_timings && rule_tag.targets().begin()->get_info().name != "OUTPUT")
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

    //rule_tag.PrintName(timeout);
    timeout << rule_tag.get_info().name();
    timeout << "\t" << count << "\t" << et-st << "\t";
    count = 0;
    st = MPI_Wtime();
    do{
      rp->compute(exec_seq);
      et = MPI_Wtime();
      count++;
    }while(et - st < min_time_duration);
      
    timeout << exec_seq.size() <<"\t" << count << "\t" << et-st << "\t";
    
    std::map<variable, std::pair<int, unsigned char *> > storeMap;
    variableSet targets = rule_tag.targets();
    for(variableSet::const_iterator vi = targets.begin(); 
	vi != targets.end(); vi++) {
      storeRepP vRep = facts.get_variable(*vi);
      int size = vRep->pack_size(vRep->domain());
      unsigned char * buf = new unsigned char[size];
      int position = 0;
      vRep->pack(buf, position, size, vRep->domain());
      std::pair<int, unsigned char *> myPair(size, buf);
      storeMap[*vi] = myPair;
    }
    
    int size_factor[4];
    if(exec_seq.size() <= 100) {
      size_factor[0] = 2;
      size_factor[1] = 4;
      size_factor[2] = 8;
      size_factor[3] = 16;
    }
    else if(exec_seq.size() <= 1000) {
      size_factor[0] = -2;
      size_factor[1] = 2;
      size_factor[2] = 4;
      size_factor[3] = 8;
    }
    else if(exec_seq.size() <= 10000) {
      size_factor[0] = -4;
      size_factor[1] = -2;
      size_factor[2] = 2;
      size_factor[3] = 4;
    }
    else if(exec_seq.size() <= 50000) {
      size_factor[0] = -8;
      size_factor[1] = -4;
      size_factor[2] = -2;
      size_factor[3] = 2;
    }
    else {
      size_factor[0] = -16;
      size_factor[1] = -8;
      size_factor[2] = -4;
      size_factor[3] = -2;
    }

    for(int i = 0; i < 4; i++) {
      sequence my_seq;
      if(size_factor[i] < 0) {
	int size = exec_seq.size()/size_factor[i]*(-1);
	sequence::const_iterator si = exec_seq.begin();
	for(int j = 0; j < size; j++) {
	  my_seq += *si;
	  si++;
	}
      }
      else {
	for(int j = 0; j < size_factor[i]; j++) {
	  my_seq += exec_seq;
	}
      }
      count = 0;
      st = MPI_Wtime();
      do{
	rp->compute(my_seq);
	et = MPI_Wtime();
	count++;
      }while(et - st < min_time_duration);
      timeout << my_seq.size() <<"\t" << count << "\t" << et-st << "\t";
    }
    timeout << endl;
    for(std::map<variable, std::pair<int, unsigned char*> >::iterator mi = storeMap.begin();
	mi != storeMap.end(); mi++) {
      storeRepP vRep = facts.get_variable(mi->first);
      int position = 0;
      vRep->unpack(mi->second.second, position, mi->second.first, vRep->domain());
      delete [] mi->second.second;
    }
  }
}

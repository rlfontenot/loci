#include "dist_tools.h" // for use of Loci::debugout
#include "comp_tools.h"
#include "visitorabs.h"
#include <vector>
using std::vector ;
#include <deque>
using std::deque ;
#include <map>
using std::map ;

#include <malloc.h>

namespace Loci {
  extern int chomping_size ;
  extern bool profile_memory_usage ;

  extern double LociAppPeakMemory ;
  extern double LociAppAllocRequestBeanCounting ;
  extern double LociAppFreeRequestBeanCounting ;
  extern double LociAppPeakMemoryBeanCounting ;
  extern double LociAppPMTemp ;
  namespace {
    // memory profile function
    int currentMem(void) {
#ifdef LINUX
      struct mallinfo info = mallinfo() ;
      return info.arena+info.hblkhd ;
#else
      cerr << "currentMem not supported" << endl ;
      return 0 ;
#endif
    }
  }

  class execute_chomp: public execute_modules {
    entitySet total_domain ;
    vector<pair<rule,rule_compilerP> > chomp_comp ;
    vector<rule_implP> chomp_compP ;
    deque<entitySet> rule_seq ;
    variableSet chomp_vars ;
    vector<vector<entitySet> > seq_table ;
    int_type chomp_size ;
    int_type chomp_iter ;
    vector<int_type> chomp_offset ;
    vector<storeRepP> chomp_vars_rep ;
    int_type D_cache_size ;
  public:
    execute_chomp(const entitySet& td,
                  const vector<pair<rule,rule_compilerP> >& comp,
                  const deque<entitySet>& seq,
                  const variableSet& cv,
                  fact_db& facts):
      total_domain(td),chomp_comp(comp),rule_seq(seq),
      chomp_vars(cv),chomp_size(0),chomp_iter(0),D_cache_size(0) {

      for(vector<pair<rule,rule_compilerP> >::const_iterator vi=comp.begin();
          vi!=comp.end();++vi)
        chomp_compP.push_back((vi->first).get_rule_implP()) ;

      for(vector<rule_implP>::iterator vi=chomp_compP.begin();
          vi!=chomp_compP.end();++vi)
        (*vi)->initialize(facts) ;

      for(variableSet::const_iterator vi=chomp_vars.begin();
          vi!=chomp_vars.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        chomp_vars_rep.push_back(srp) ;
      }

      // we'll need to set up the chomp_size
      // and the chomping sequence table
      entitySet test_alloc = interval(1,1) ;
      int_type total_obj_size = 0 ;
      int_type min_store_size = UNIVERSE_MAX ;
      D_cache_size = chomping_size * 1024 ; // assume 128KB now
      for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
          vi!=chomp_vars_rep.end();++vi) {
        int_type my_size = (*vi)->pack_size(test_alloc) ;
        total_obj_size += my_size ;
        min_store_size = min(min_store_size,my_size) ;
      }

      double scale = (double)D_cache_size / (double)total_obj_size ;

      if(scale <= 1.0)
        chomp_size = 1 ;
      else
        chomp_size = (int)scale ;

      entitySet copy_total_domain = total_domain ;
      chomp_offset.clear() ;
      int_type low_pos ;
      int_type high_pos ;
      while(copy_total_domain != EMPTY) {
        low_pos = copy_total_domain.Min() ;
        high_pos = low_pos + chomp_size - 1 ;
        entitySet seg = interval(low_pos, high_pos) ;
        vector<entitySet> seq_vec ;
        for(deque<entitySet>::const_iterator vi=rule_seq.begin();
            vi!=rule_seq.end();++vi)
          seq_vec.push_back(*vi & seg) ;
        
        seq_table.push_back(seq_vec) ;

        copy_total_domain -= seg ;

        if(copy_total_domain != EMPTY)
          chomp_offset.push_back(copy_total_domain.Min() - low_pos) ;
        else
          chomp_offset.push_back(0) ;
      }

      chomp_iter = seq_table.size() ;
    }
    
    virtual void execute(fact_db& facts) ;
    virtual void Print(std::ostream& s) const ;
  } ;

  void execute_chomp::execute(fact_db& facts) {
    if(total_domain == EMPTY) {
      // call the compute() method at least once
      vector<rule_implP>::iterator vri ;
      for(vri=chomp_compP.begin();vri!=chomp_compP.end();++vri) {
        (*vri)->compute(sequence(EMPTY)) ;
      }
      return ;
    }

    {
      entitySet first_alloc =
        entitySet(interval(total_domain.Min(),
                           total_domain.Min()+chomp_size-1)
                  ) ;
      // first we need to allocate the chunks of variables
      for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
          vi!=chomp_vars_rep.end();++vi) {
        (*vi)->allocate(first_alloc) ;
      }
    }
    // do memory profiling
    if(profile_memory_usage) {
      double currmen = currentMem() ;
      if(currmen > LociAppPeakMemory)
        LociAppPeakMemory = currmen ;
      LociAppAllocRequestBeanCounting += D_cache_size ;
      LociAppPMTemp += D_cache_size ;
      if(LociAppPMTemp > LociAppPeakMemoryBeanCounting)
        LociAppPeakMemoryBeanCounting = LociAppPMTemp ;
    }

    // begin execution, the loop number would be seq_table.size()
    vector<vector<entitySet> >::const_iterator vvi ;
    int count = 0 ;
    for(vvi=seq_table.begin();vvi!=seq_table.end();++vvi,++count) {
      vector<entitySet>::const_iterator vei ;
      vector<rule_implP>::iterator vri ;
      for(vri=chomp_compP.begin(),vei=vvi->begin();
          vri!=chomp_compP.end();++vri,++vei) {
        (*vri)->compute(sequence(*vei)) ;
      }
      // we shift the alloc domain for each chomp_vars_repS
      // first get the offset
      int_type offset = chomp_offset[count] ;
      if(offset != 0) {
        for(vector<storeRepP>::iterator vsi=chomp_vars_rep.begin();
            vsi!=chomp_vars_rep.end();++vsi)
          (*vsi)->shift(offset) ;
      }
    }
    // at last, we deallocate all the chomp_vars
    for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
        vi!=chomp_vars_rep.end();++vi) {
      (*vi)->allocate(EMPTY) ;
    }
    
    // do memory profiling
    if(profile_memory_usage) {
      LociAppFreeRequestBeanCounting += D_cache_size ;
      LociAppPMTemp -= D_cache_size ;
      if(LociAppPMTemp < 0)
        std::cout << "MEMORY PROFILING WARNING: negative memory size"
                  << endl ;
    }
  }

  void execute_chomp::Print(std::ostream& s) const {
    s << "--Start chomping: (chomping interval size: "
      << chomp_size << ", iter number: " << chomp_iter
      << ", total domain: " << total_domain
      << ")" << endl ;
    s << "--Perform chomping for the following rule sequence: " << endl ;
    for(vector<pair<rule,rule_compilerP> >::const_iterator
          vi=chomp_comp.begin();vi!=chomp_comp.end();++vi)
      s << "-- " << vi->first << endl ;
    s << "--End chomping" << endl ;
  }
  
  chomp_compiler::chomp_compiler(const digraph& cgraph,
                                 const variableSet& cvars,
                                 const map<rule,rule>& a2u)
    :chomp_graph(cgraph),chomp_vars(cvars),apply2unit(a2u) {}

  void chomp_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    vector<pair<rule,rule_compilerP> >::iterator i ;
    for(i=chomp_comp.begin();i!=chomp_comp.end();++i) {
      rule r = i->first ;
      rule_compilerP bc = i->second ;

      if(r.get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
        existential_rule_analysis(r,facts,scheds) ;
      else
        existential_applyrule_analysis(r,facts,scheds) ;
    }
  }
  
  void chomp_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    vector<pair<rule,rule_compilerP> >::reverse_iterator ri ;
    for(ri=chomp_comp.rbegin();ri!=chomp_comp.rend();++ri) {
      rule r = ri->first ;
      rule_compilerP bc = ri->second ;
      
      entitySet exec_seq ;
      if(r.get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
        exec_seq = process_rule_requests(r,facts,scheds) ;
      else {
        map<rule,rule>::const_iterator mi ;
        mi = apply2unit.find(r) ;
        FATAL(mi == apply2unit.end()) ;
        exec_seq = process_applyrule_requests(r,mi->second,facts,scheds) ;
      }
      rule_seq.push_front(exec_seq) ;
    }
  }

  executeP chomp_compiler::create_execution_schedule(fact_db& facts,
                                                     sched_db& scheds) {
    // we first union all the rule sequence
    entitySet total ;
    
    deque<entitySet>::size_type i ;
    
    for(i=0;i<rule_seq.size();++i)
      total += rule_seq[i] ;
    
    return new execute_chomp(total,chomp_comp,rule_seq,chomp_vars,facts) ;
  }
  
} // end of namespace Loci

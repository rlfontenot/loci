#include "comp_tools.h"
#include "dist_tools.h"
using std::ostream ;
using std::endl ;

namespace Loci {

  extern int current_rule_id ;
  
  class execute_map_rule: public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    variableSet sources ; // source vars of the map rule
    variableSet targets ; // target vars of the map rule
  public:
    // tsv is the target and source variables of this map rule
    execute_map_rule(rule fi, sequence seq,
                     const variableSet& svars,
                     const variableSet& tvars,
                     fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  execute_map_rule::
  execute_map_rule(rule fi, sequence seq,
                   const variableSet& svars,
                   const variableSet& tvars,
                   fact_db& facts, sched_db& scheds) {
    rp = fi.get_rule_implP() ;
    sources = svars ;
    targets = tvars ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    control_thread = false ;
  }
  
  void execute_map_rule::execute(fact_db &facts) {
    current_rule_id = rule_tag.ident() ;
    // before the execute begins, we need to restore
    // all the facts associated with this map rule
    // to their global numbering scheme
    if(facts.is_distributed_start()) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      for(variableSet::const_iterator vi=sources.begin();
          vi!=sources.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        // Note, we can use facts.replace_fact here
        // But doing so will require rebind the rule
        // to the fact_db
        facts.update_fact(*vi,srp->remap(df->dl2g)) ;
        //facts.replace_fact(*vi,srp->remap(l2g)) ;
      }
    }
    rp->compute(exec_seq) ;
    // and then after computing, we'll convert all
    // all the facts associated with this map rule
    // to their local numbering scheme
    // But since we have generated new maps
    // we could possibly end up having new entities
    // included in the clone region on local process.
    // Therefore we should only transform all the source
    // facts since transforming the generated maps
    // at last will likely to lose data.
    if(facts.is_distributed_start()) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;

      for(variableSet::const_iterator vi=sources.begin();
          vi!=sources.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        facts.update_fact(*vi,srp->remap(df->g2l)) ;
        //facts.replace_fact(*vi,srp->remap(df->g2l)) ;
      }
    }
  }
  
  void execute_map_rule::Print(ostream &s) const {
    s << rule_tag << "  over sequence " << exec_seq << endl ;
  }
  
  void map_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    existential_rule_analysis(map_impl, facts, scheds) ;
  }
  
  void map_compiler::process_var_requests(fact_db& facts,
                                          sched_db& scheds) {
    //variableSet sources = map_impl.sources() ;
    //variableSet::const_iterator vi ;
    //for(vi=sources.begin();vi!=sources.end();++vi) {
    //  scheds.variable_request(*vi,scheds.variable_existence(*vi)) ;
    //}
    exec_seq = process_rule_requests(map_impl, facts, scheds) ;
  }
  
  executeP map_compiler::create_execution_schedule(fact_db& facts,
                                                   sched_db& scheds) {
    //return new execute_map_rule(map_impl, ~EMPTY, facts,scheds) ;
    return new execute_map_rule(map_impl, exec_seq,
                                map_impl.sources(),
                                map_impl.targets(),
                                facts,scheds) ;
  }
  
}

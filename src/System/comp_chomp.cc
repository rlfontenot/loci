#include "comp_tools.h"
#include "visitorabs.h"
#include <vector>
using std::vector ;
#include <deque>
using std::deque ;


namespace Loci {

  class execute_chomp: public execute_modules {
    entitySet total_domain ;
    vector<rule> chomp_comp ;
    vector<rule_implP> chomp_compP ;
    deque<entitySet> rule_seq ;
    variableSet chomp_vars ;
    vector<vector<entitySet> > seq_table ;
    bool is_table_set ;
    int_type chomp_size ;
    vector<int_type> chomp_offset ;
    vector<storeRepP> chomp_vars_rep ;
  public:
    execute_chomp(const entitySet& td,
                  const vector<rule>& comp,
                  const deque<entitySet>& seq,
                  const variableSet& cv,
                  fact_db& facts):
      total_domain(td),chomp_comp(comp),rule_seq(seq),
      chomp_vars(cv),is_table_set(false) {

      for(vector<rule>::const_iterator vi=comp.begin();
          vi!=comp.end();++vi)
        chomp_compP.push_back(vi->get_rule_implP()) ;

      for(vector<rule_implP>::iterator vi=chomp_compP.begin();
          vi!=chomp_compP.end();++vi)
        (*vi)->initialize(facts) ;

      for(variableSet::const_iterator vi=chomp_vars.begin();
          vi!=chomp_vars.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        chomp_vars_rep.push_back(srp) ;
      }
      
      int_type dom_min = total_domain.Min() ;
      int_type dom_max = total_domain.Max() ;
      int_type tmp = (dom_max - dom_min)/chomp_vars.size() ;
      chomp_size = (tmp<=50)?tmp:50 ;
    }
    
    virtual void execute(fact_db& facts) ;
    virtual void Print(std::ostream& s) const ;
  } ;

  void execute_chomp::execute(fact_db& facts) {
    if(total_domain == EMPTY)
      return ;

    // first of all, we build the table, if necessary
    if(!is_table_set) {
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
        
        chomp_offset.push_back(copy_total_domain.Min() - low_pos) ;
      }

      is_table_set = true ;
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
        vi!=chomp_vars_rep.end();++vi)
      (*vi)->allocate(EMPTY) ;
  }

  void execute_chomp::Print(std::ostream& s) const {
    s << "--Start chomping: (chomping interval size: "
      << chomp_size << ")" << endl ;
    s << "--Perform chomping for the following rule sequence: " << endl ;
    for(vector<rule>::const_iterator vi=chomp_comp.begin();
        vi!=chomp_comp.end();++vi)
      s << "-- " << *vi << endl ;
    s << "--End chomping" << endl ;
  }
  
  chomp_compiler::chomp_compiler(const digraph& cgraph,
                                 const variableSet& cvars)
    :chomp_graph(cgraph),chomp_vars(cvars) {}

  void chomp_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    vector<rule>::iterator i ;
    for(i=chomp_comp.begin();i!=chomp_comp.end();++i)
      existential_rule_analysis(*i,facts,scheds) ;
  }
  
  void chomp_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    vector<rule>::reverse_iterator ri ;
    for(ri=chomp_comp.rbegin();ri!=chomp_comp.rend();++ri)
      rule_seq.push_front(process_rule_requests(*ri,facts,scheds)) ;
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

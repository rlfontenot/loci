#include "sched_tools.h"
using std::cout ;
using std::cerr ;
using std::endl ;
namespace Loci {
dynamic_schedule_rule::dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds)  {

  rp = fi.get_rule_implP() ; //get rule from rule database (rhs)
  rule_tag = fi ;            //store rule tag in rule_tag
  local_compute = rp->new_rule_impl() ; //another instance of rule rhs     
  entitySet in = rule_tag.sources() ; //inputs from rule rhs
  outputs = rule_tag.targets() ;      //outputs as in rhs  
  exec_set = eset ;    

  // Setup local facts input variables (types only no allocation)
  for(variableSet::const_iterator vi=in.begin();vi!=in.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    if((store_ptr != 0) && store_ptr->RepType() == Loci::STORE) {	  
      inputs += *vi ; 
      local_facts.create_fact(*vi,store_ptr->new_store(EMPTY)) ;      
    } else {
      local_facts.create_fact(*vi,facts.get_variable(*vi)) ;
    }
  }
  // Setup local facts output variables 
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    local_facts.create_fact(*vi,store_ptr->new_store(EMPTY)) ;
  }

  // Initialize both functions for remote and local execution.
  local_compute->initialize(local_facts) ;
  rp->initialize(facts) ;
}

void dynamic_schedule_rule::execute(fact_db &facts) { 
  for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
    storeRepP sp = local_facts.get_variable(*vi) ;
    sp->allocate(exec_set) ;
  }
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP sp = local_facts.get_variable(*vi) ;
    sp->allocate(exec_set) ;
  }
  
  // Pack inputs from facts
  int position = 0 ;
  int size = 0 ;
  for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
    storeRepP s_ptr = facts.get_variable(*vi) ;
    size += s_ptr->pack_size(exec_set) ;
  }
  unsigned char *buf = new unsigned char[size] ;
  for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
    storeRepP s_ptr = facts.get_variable(*vi) ;
    s_ptr->pack(buf,position,size,exec_set) ;
  }

  // unpack inputs into local facts
  position = 0 ;
  for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
    storeRepP s_ptr = local_facts.get_variable(*vi) ;
    s_ptr->unpack(buf,position,size,sequence(exec_set)) ;
  }

  //    local_facts.Print_diagnostics();
  
  //  rp->compute(sequence(exec_set)) ;

  // Now we will compute in seperate memory and results will be in output
  local_compute->compute(sequence(exec_set)) ;

  // Pack output
  position = 0 ;
  size = 0 ;
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP s_ptr = local_facts.get_variable(*vi) ;
    size += s_ptr->pack_size(exec_set) ;
  }
  delete[] buf ;
  buf = new unsigned char[size] ;
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP s_ptr = local_facts.get_variable(*vi) ;
    s_ptr->pack(buf,position,size,exec_set) ;
  }

  // unpack outputs into facts
  position = 0 ;
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP s_ptr = facts.get_variable(*vi) ;
    s_ptr->unpack(buf,position,size,sequence(exec_set)) ;
  }

  // Delete allocated temporaries
  delete[] buf ;

  for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
    storeRepP sp = local_facts.get_variable(*vi) ;
    sp->allocate(EMPTY) ;
  }
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP sp = local_facts.get_variable(*vi) ;
    sp->allocate(EMPTY) ;
  }
}

void dynamic_schedule_rule::Print(std::ostream &s) const {

 s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl ; 
  
}

}

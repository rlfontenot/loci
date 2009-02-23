//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include "comp_tools.h"
#include "dist_tools.h"
#include "loci_globs.h"
using std::ostream ;
using std::endl ;

namespace Loci {

  int current_rule_id = 0 ;
  int rule_count = 0;

  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
  }

  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p, const sched_db &scheds) {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
  }

  void execute_rule::execute(fact_db &facts) {
    stopWatch s ;
    s.start() ;
    current_rule_id = rule_tag.ident() ;
    rp->compute(exec_seq);
    current_rule_id = 0 ;
    timer.addTime(s.stop(),exec_size) ;
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

  
  void execute_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "rule: "<<rule_tag ;

    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }
  
  execute_modules_decorator_factory* impl_compiler::decoratorFactory = NULL;

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

    if (impl.get_info().rule_impl->dynamic_schedule_rule() && use_dynamic_scheduling) {
      executeP execute_dynamic = new dynamic_schedule_rule(impl,exec_seq,facts, scheds) ;
      if(decoratorFactory != NULL)
        execute_dynamic = decoratorFactory->decorate(execute_dynamic);
      return execute_dynamic;
    }
    executeP exec_rule = new execute_rule(impl,sequence(exec_seq),facts, scheds);
    if(decoratorFactory != NULL)
      exec_rule = decoratorFactory->decorate(exec_rule);
    return exec_rule;
  }

  // blackbox_compiler code
  void
  blackbox_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    // set UNIVERSE existence for all targets
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();++vi)
      scheds.set_existential_info(*vi, impl, ~EMPTY) ;
  }

  void blackbox_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    // everyone will need to request for their existence
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi = targets.begin(); vi != targets.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;

    variableSet sources = impl.sources() ;
    for(variableSet::const_iterator vi=sources.begin(); vi != sources.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
  }

  execute_modules_decorator_factory* blackbox_compiler::decoratorFactory = NULL;

  executeP blackbox_compiler::create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_rule(impl, ~EMPTY, facts, scheds);
    if (decoratorFactory != NULL)
      execute = decoratorFactory->decorate(execute);
    return execute;
  }

  // superRule_compiler code
  void
  superRule_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_existential(impl,facts,scheds) ;
  }

  void superRule_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_requests(impl,facts,scheds) ;
  }

  execute_modules_decorator_factory* superRule_compiler::decoratorFactory = NULL;

  executeP superRule_compiler::create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_rule(impl, ~EMPTY, facts, scheds);
    if (decoratorFactory != NULL)
      execute = decoratorFactory->decorate(execute);
    return execute;
  }


  // Lets set up some common super rule functions

  class NOT_rule : public super_rule {
    param<bool> NOT ;
  public:
    NOT_rule() {
      name_store("NOT(X)",NOT) ;
      constraint("X") ;
      output("NOT(X)") ;
    }
    void compute(const sequence &seq) {
      *NOT = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = ~EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      debugout << "constraints = " << constraints << endl ;
      // Now complement (entities we own)
      constraints = (~constraints) & my_entities ;
      debugout << "constraints & my_entities =" << constraints << endl ;

      variableSet output = r.targets() ;
      debugout << "setting " << output << endl ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<NOT_rule> register_NOT_rule ;

  class OR_rule : public super_rule {
    param<bool> OR ;
  public:
    OR_rule() {
      name_store("OR(X,Y)",OR) ;
      constraint("X,Y") ;
      output("OR(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR_rule> register_OR_rule ;

  class OR3_rule : public super_rule {
    param<bool> OR ;
  public:
    OR3_rule() {
      name_store("OR(X,Y,Z)",OR) ;
      constraint("X,Y,Z") ;
      output("OR(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR3_rule> register_OR3_rule ;

  class OR4_rule : public super_rule {
    param<bool> OR ;
  public:
    OR4_rule() {
      name_store("OR(W,X,Y,Z)",OR) ;
      constraint("W,X,Y,Z") ;
      output("OR(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR4_rule> register_OR4_rule ;

  class AND2_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND2_rule() {
      name_store("AND(X,Y)",AND) ;
      constraint("X,Y") ;
      output("AND(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND2_rule> register_AND2_rule ;

  class AND3_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND3_rule() {
      name_store("AND(X,Y,Z)",AND) ;
      constraint("X,Y,Z") ;
      output("AND(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND3_rule> register_AND3_rule ;

  class AND4_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND4_rule() {
      name_store("AND(W,X,Y,Z)",AND) ;
      constraint("W,X,Y,Z") ;
      output("AND(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND4_rule> register_AND4_rule ;
    
  
}

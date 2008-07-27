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
    rp->compute(exec_seq);
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

}

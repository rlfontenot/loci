//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include "loci_globs.h"
#include "sched_tools.h"
#include "dist_tools.h"
#include "param_rule.h"
#include "thread.h"
#include <Tools/except.h>
#include <constraint.h>
#include <new>
#include "comp_tools.h"

using std::bad_alloc ;
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ios ;
using std::ofstream ;
using std::istream ;
using std::ostream ;

///////////////////////////////////
#include <sstream>
#include <algorithm>
#include <sys/time.h> // for gettimeofday function
///////////////////////////////////

namespace Loci {
  namespace {
    rule create_rule(variable sv, variable tv, string qualifier) {
      ostringstream oss ;
      oss << "source(" << sv << ')' ;
      oss << ",target(" << tv << ')' ;
      oss << ",qualifier(" << qualifier << ')' ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
  }

  inline variable makeGPUVAR(variable v) {
    variable::info vinfo = v.get_info() ;
    string vname = string("__GPU__") + vinfo.name ;
    vinfo.name = vname ;
    return variable(vinfo) ;
  }
  rule_db rename_gpu_containers(fact_db  &facts,
				const rule_db &rdb)  {

    rule_db gpu_rdb ;

    variableSet gpuInputs ;
    variableSet gpuOutputs ;
    variableSet inputs ;
    variableSet outputs ;
    ruleSet rset = rdb.all_rules() ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end();++rsi) {
      rule_implP rp = rsi->get_rule_implP() ;
      if(rp == 0) {
	// Not a rule with an implementation so no GPU container
	// types
	gpu_rdb.add_rule(*rsi) ;
        continue ;
      }
      variableSet targets = rsi->targets() ;
      variableSet sources = rsi->sources() ;
      bool hasGPUVar = false ;
      bool hasCPUVar = false ;
      map<variable,variable> vm ;
      for(variableSet::const_iterator vsi = targets.begin(); vsi !=
	    targets.end(); ++vsi) {
	storeRepP sp = rp->get_store(*vsi) ;
	if(isGPU(sp)) {
	  gpuOutputs += *vsi ;
	  hasGPUVar = true ;
	  vm[*vsi] = makeGPUVAR(*vsi) ;
	} else {
	  outputs += *vsi ;
	  hasCPUVar = true ;
	  vm[*vsi] = *vsi ;
	}
      }
      for(variableSet::const_iterator vsi = sources.begin(); vsi !=
	    sources.end(); ++vsi) {
	storeRepP sp = rp->get_store(*vsi) ;
	if(isGPU(sp)) {
	  gpuInputs += *vsi ;
	  hasGPUVar = true ;
	  vm[*vsi] = makeGPUVAR(*vsi) ;
	} else {
	  inputs += *vsi ;
	  hasCPUVar = true ;
	  vm[*vsi] = *vsi ;
	}
      }
      if(hasGPUVar) {
	rp->rename_vars(vm) ;
	gpu_rdb.add_rule(rule(rp)) ;
      } else {
	gpu_rdb.add_rule(*rsi) ;
      }
      
    }

    cout << "inputs = " << inputs << endl ;
    cout << "outputs = " << outputs << endl ;
    cout << "gpuInputs = " << gpuInputs << endl ;
    cout << "gpuOutputs = " << gpuOutputs << endl ;

    variableSet cpu2gpu = gpuInputs ;
    cpu2gpu -= gpuOutputs ;
    for(variableSet::const_iterator vsi = cpu2gpu.begin(); vsi !=
	  cpu2gpu.end(); ++vsi) {
      gpu_rdb.add_rule(create_rule(*vsi,makeGPUVAR(*vsi),"cpu2gpu")) ;
    }
    variableSet gpu2cpu = gpuOutputs ;
    gpuOutputs -= gpuInputs ;
    for(variableSet::const_iterator vsi = gpu2cpu.begin(); vsi !=
	  gpu2cpu.end(); ++vsi) {
      gpu_rdb.add_rule(create_rule(makeGPUVAR(*vsi),*vsi,"gpu2cpu")) ;
    }    

    return gpu_rdb ;
  }

  void gpu2cpu_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(r,facts, scheds) ;
  }

  void gpu2cpu_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = process_rule_requests(r,facts, scheds) ;
    scheds.update_exec_seq(r, exec_seq);
  }

  executeP gpu2cpu_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(r) ;
    variable vgpu = *r.sources().begin() ;
    storeRepP p = facts.get_variable(vgpu)->getRep() ;
    gpuRepP gp = gpuRepP(p) ; 
    variable vcpu = *r.targets().begin() ;
    storeRepP cp = facts.get_variable(vcpu) ;
    executeP execute = executeP(new execute_gpu2cpu_copy(r,gp,cp,exec_seq)) ;
    return execute;
  }

  void cpu2gpu_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(r,facts, scheds) ;
  }

  void cpu2gpu_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = process_rule_requests(r,facts, scheds) ;
    scheds.update_exec_seq(r, exec_seq);
  }

  executeP cpu2gpu_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(r) ;
    variable vgpu = *r.targets().begin() ;
    storeRepP p = facts.get_variable(vgpu)->getRep() ;
    gpuRepP gp = gpuRepP(p) ; 
    variable vcpu = *r.sources().begin() ;
    storeRepP cp = facts.get_variable(vcpu);

    executeP execute = executeP(new execute_cpu2gpu_copy(r,gp,cp,exec_seq)) ;
    return execute;
  }

  void execute_gpu2cpu_copy::execute(fact_db &facts, sched_db &scheds) {
    gpuvar->copyTo(cpuvar,copyset) ;
  }

  void execute_gpu2cpu_copy::Print(ostream &s) const {
    printIndent(s) ;
    s << r << " over sequence " ;
    if(verbose || copyset.num_intervals() < 4) {
      s << copyset << endl ;
    } else {
      s << "[ ... ], l=" << copyset.size() << endl ;
    }
  }

  void execute_gpu2cpu_copy::dataCollate(collectData &data_collector) const {
    //    ostringstream oss ;
    //    oss << "rule: "<<rule_tag ;
    //
    //    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  void execute_cpu2gpu_copy::execute(fact_db &facts, sched_db &scheds) {
    gpuvar->copyFrom(cpuvar,copyset) ;
  }

  void execute_cpu2gpu_copy::Print(ostream &s) const {
    printIndent(s) ;
    s << r << " over sequence " ;
    if(verbose || copyset.num_intervals() < 4) {
      s << copyset << endl ;
    } else {
      s << "[ ... ], l=" << copyset.size() << endl ;
    }
  }

  void execute_cpu2gpu_copy::dataCollate(collectData &data_collector) const {
    //    ostringstream oss ;
    //    oss << "rule: "<<rule_tag ;
    //
    //    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }
  

}

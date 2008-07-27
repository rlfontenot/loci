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

#ifndef PERFORMANCE_ANALYSIS_H
#define PERFORMANCE_ANALYSIS_H

#include <iostream>
using std::iostream;
using std::endl;
#include <stdio.h>
#include <map>
using std::map;
using std::multimap;
#include <string>
using std::string;
#include <iomanip>
#include <mpi.h>
using namespace std;
#include <fact_db.h>
#include <sched_db.h>
#include <execute.h>

#ifdef USE_PAPI
#include <papi.h>
#include <papiStdEventDefs.h>
#endif

namespace Loci {
#ifdef USE_PAPI
  typedef long_long timeType;
#else
  typedef double timeType;
#endif

  typedef void * timer_token;
	
  class performance_data;
	
  class performance_analysis {
    performance_data *data; // internal data objects
  public:
    performance_analysis(void);
    void add2RuleTimingsTable(string id, string rule_name, int size, timeType time);
    void * start_timer(string key);
    void stop_timer(void * token);
    void create_report(void);
  private:
    void processPerformanceData();
    void printRuleTimingsDetails();
    void printRuleTimingsSummary();
    void printExecuteModuleTable();
    void printCompilerTimingTable();
    void printBasicExecuteModuleTable();
    void printBasicCompilerTimingTable();
    void printPhaseTimingTable();
  };

  class execute_modules_timer : public execute_modules {
    executeP wrapped;
    performance_analysis* perf_analysis;
    string compiler_key, module_key, timer_key;
  public:
    execute_modules_timer(executeP module, performance_analysis* analysis, string compiler_key);
    virtual ~execute_modules_timer() {};
    virtual void execute(fact_db &facts);
    virtual void Print(std::ostream &s) const;
    virtual string getName() { return "execute_modules_timer";};
    virtual void dataCollate(collectData &data_collector) const {}

  };
	
  class execute_modules_decorator_factory {
  public:
    virtual executeP decorate(executeP module) = 0;
    virtual ~execute_modules_decorator_factory() {}
  };
	
  class execute_modules_timer_factory : public execute_modules_decorator_factory {
    performance_analysis* perf_analysis;
    string key;
  public: 
    execute_modules_timer_factory(performance_analysis* analysis, string module_key);
    virtual ~execute_modules_timer_factory() {};
    virtual executeP decorate(executeP module);
  };
} //end of namespace Loci
#endif

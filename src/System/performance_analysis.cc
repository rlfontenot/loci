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

#include "performance_analysis.h"
#include "loci_globs.h"

namespace Loci {
	
	struct rule_timing_info {
		int context_size;
		timeType total_time;
	};
	
	struct timer_data {
		string key;
		timeType start_time, stop_time, time;
		int count;
	};
	
	class performance_data {
	public :
		map<string, rule_timing_info*> rule_timings_table;
		multimap<string, timeType> execute_timings;
		map<string, timeType> phases_timing;
		multimap<string, timer_data*> execute_module_timing;
		multimap<string, timer_data*> compiler_timing;
	};

	performance_analysis::performance_analysis(void) {
		data = new performance_data();
	}

	timer_token performance_analysis::start_timer(string key) {
		timer_data *d = new timer_data;
		d->key = key;
#ifdef USE_PAPI
	    d->start_time = PAPI_get_real_usec();
#else
	    d->start_time = MPI_Wtime();
#endif
		return d;
	}
  
	void performance_analysis::stop_timer(timer_token token) {
		timer_data *d = (timer_data *) token;
		timeType time_elapsed;
#ifdef USE_PAPI
	    d->stop_time = PAPI_get_real_usec();
	    time_elapsed = (d->stop_time - d->start_time) / 10000;
#else
	    d->stop_time = MPI_Wtime();
	    time_elapsed = d->stop_time - d->start_time;
#endif  
	    data->execute_timings.insert(pair<string, timeType> (d->key, abs(time_elapsed)));
      }
	
	void performance_analysis::add2RuleTimingsTable(string rule_name, int size, timeType time) {
	      rule_timing_info *rule_info = new rule_timing_info;
	      rule_info->context_size = size;
	      rule_info->total_time = time;
	      data->rule_timings_table[rule_name] = rule_info;
	}

	void performance_analysis::create_report() {
		perfReportOut << setiosflags(ios::left) << endl;
#ifdef USE_PAPI
		perfReportOut << "The following performance timings were collected using PAPI" <<endl;
#endif
		processPerformanceData();
		if(minimum_performance_analysis) {
			printPhaseTimingTable();
			printBasicExecuteModuleTable();
		} else if(full_performance_analysis) {
			printPhaseTimingTable();
			printBasicExecuteModuleTable();
			printBasicCompilerTimingTable();
			printRuleTimingsSummary();
		} else if(detailed_performance_analysis) {
			printPhaseTimingTable();
			printExecuteModuleTable();
			printCompilerTimingTable();
			printRuleTimingsSummary();
			printRuleTimingsDetails();
		} else if(measure_rule_timings) {
			printRuleTimingsSummary();
			printRuleTimingsDetails();
		}
	}

	void performance_analysis::processPerformanceData() {
		multimap<string, timeType>::iterator itr;
		timer_data *execute_data, *compiler_data;
		timeType total_exec_time = 0, comm_time = 0, comp_time = 0, io_time = 0;
		string exec_module_name, compiler_name;
		
		for(itr = data->execute_timings.begin(); itr != data->execute_timings.end(); itr++) {
			int count = (int) data->execute_timings.count((*itr).first);
			string key = (*itr).first;
			size_t pos = key.find("-");
			if (pos > 0 && pos < key.length()) {
				exec_module_name = key.substr(0, pos);
				compiler_name = key.substr(pos+1);
			} else {
				exec_module_name = key;
				compiler_name = "";
			}
			total_exec_time = 0;

			for(int i = 0; i < count; i++) {
				if(exec_module_name == "execute_comm" || exec_module_name == "execute_comm_reduce") {
					comm_time += (*itr).second;
				} else if(exec_module_name == "execute_rule" || exec_module_name == "execute_rule_null") {
					comp_time += (*itr).second;
				} else if(exec_module_name == "execute_msg") {
					io_time +=  (*itr).second;
				}
				total_exec_time += (*itr).second;
				itr++;
			}
			// Put data in separate tables
			if (compiler_name == "") {
				data->phases_timing[exec_module_name] = total_exec_time;
			} else {
				execute_data = new timer_data;
				execute_data->key = compiler_name;
				execute_data->time = total_exec_time;
				execute_data->count = count;
				compiler_data = new timer_data;
				compiler_data->key = exec_module_name;
				compiler_data->time = total_exec_time;
				compiler_data->count = count;
				data->execute_module_timing.insert(pair<string, timer_data*>(exec_module_name, execute_data));
				data->compiler_timing.insert(pair<string, timer_data*>(compiler_name, compiler_data));
			}
			itr--;
		}
		data->phases_timing["Communication Time"] = comm_time;
		data->phases_timing["Computation Time"] = comp_time;
		data->phases_timing["I/O Time"] = io_time;		
	}
	
	void performance_analysis::printPhaseTimingTable() {
		multimap<string, timeType>::iterator itr;

		perfReportOut << endl;
		perfReportOut << "Execution Time Summary" << endl;
		perfReportOut << "================================" << endl;
		perfReportOut << endl;

		perfReportOut << setw(15) << "" << setw(20) << "Phase" << setw(20) << "Execution Time" << endl;
		perfReportOut << setw(5)  << "" << setw(30) << "---------------------------" << setw(25) << "--------------" << endl;
		
		timeType phase_time = 0;
		for(itr = data->phases_timing.begin(); itr != data->phases_timing.end(); itr++) {
			string phase_name = (*itr).first;
			phase_time = (*itr).second;
			perfReportOut << setw(5) << "" << setw(33) << phase_name << "" << setw(20) << phase_time << endl;
		}
		perfReportOut << setw(5) << "" << "============================================" << endl;		
		perfReportOut << endl;
	}

	void performance_analysis::printBasicExecuteModuleTable() {
		multimap<string, timer_data*>::iterator itr;
		
		perfReportOut << endl;
		perfReportOut << "Execution Time of execution_modules::execute()" << endl;
		perfReportOut << "========================================================" << endl;
		perfReportOut << endl;
		
		perfReportOut << setw(5) << "" << setw(25) << "Execute Module" << setw(12) << "Call Count" << setw(20) << "Execution Time" << endl;
		perfReportOut << setw(5) << "" << setw(85) << "======================== =========== ===============" << endl;
		for(itr = data->execute_module_timing.begin(); itr != data->execute_module_timing.end(); itr++) {
			int count = (int) data->execute_module_timing.count((*itr).first);
			perfReportOut << setw(5) << "" << setw(25) << (*itr).first;
			timeType total_time = 0;
			int total_count = 0;
			for (int i = 0; i < count; i++) {
				timer_data *info = (*itr).second;
				total_time += info->time;
				total_count += info->count;
				itr++;
			}
			perfReportOut << setw(12) << total_count << setw(20) << total_time << endl;
			if (itr == data->execute_module_timing.end())
				perfReportOut << setw(5) << "" << setw(85) << "====================================================" << endl;
			itr--;
		}
		perfReportOut << endl;
	}
	
	void performance_analysis::printExecuteModuleTable() {
		multimap<string, timer_data*>::iterator itr;
		
		perfReportOut << endl;
		perfReportOut << "Execution Time of execution_modules::execute()" << endl;
		perfReportOut << "========================================================" << endl;
		perfReportOut << endl;
		
		perfReportOut << setw(5) << "" << setw(25) << "Execute Module" << setw(20) << "Compiler" 
						<< setw(12) << "Call Count" << setw(20) << "Execution Time" << endl;
		perfReportOut << setw(5) << "" << setw(85) << "======================== =================== =========== ===============" << endl;
		for(itr = data->execute_module_timing.begin(); itr != data->execute_module_timing.end(); itr++) {
			int count = (int) data->execute_module_timing.count((*itr).first);
			perfReportOut << setw(5) << "" << setw(25) << (*itr).first;
			timeType total_time = 0;
			for (int i = 0; i < count; i++) {
				timer_data *info = (*itr).second;
				if (i == 0)
					perfReportOut << setw(20) << info->key << setw(12) << info->count << setw(20) << info->time << endl;
				else 
					perfReportOut << setw(30) << "" << setw(20) << info->key << setw(12) << info->count << setw(20) << info->time << endl;
				total_time += info->time;	
				itr++;
			}
			perfReportOut << setw(40) << "" << setw(22) << "Total execution time:" << setw(20) << total_time << endl;
			if (itr != data->execute_module_timing.end())
				perfReportOut << setw(5) << "" << setw(85) << "------------------------ ------------------- ----------- ---------------" << endl;
			else 
				perfReportOut << setw(5) << "" << setw(85) << "========================================================================" << endl;
			itr--;
		}
		perfReportOut << endl;
	}

	void performance_analysis::printBasicCompilerTimingTable() {
		multimap<string, timer_data*>::iterator itr;
		
		perfReportOut << endl;
		perfReportOut << "Execution Time per Compiler" << endl;
		perfReportOut << "=====================================" << endl;
		perfReportOut << endl;

		perfReportOut << setw(5) << "" << setw(20) << "Compiler" << setw(20) << "Execution Time" << endl;
		perfReportOut << setw(5) << "" << setw(85) << "=================== ===============" << endl;
		for(itr = data->compiler_timing.begin(); itr != data->compiler_timing.end(); itr++) {
			int count = (int) data->compiler_timing.count((*itr).first);
			perfReportOut << setw(5) << "" << setw(20) << (*itr).first;
			timeType total_time = 0;
			int total_count = 0;
			for (int i = 0; i < count; i++) {
				timer_data *info = (*itr).second;
				total_time += info->time;
				total_count += info->count;
				itr++;
			}
			perfReportOut << setw(20) << total_time << endl;
			if (itr == data->compiler_timing.end())
				perfReportOut << setw(5) << "" << setw(85) << "===================================" << endl;
			itr--;		
		}
		perfReportOut << endl;
	}
		
	void performance_analysis::printCompilerTimingTable() {
		multimap<string, timer_data*>::iterator itr;
		
		perfReportOut << endl;
		perfReportOut << "Execution Time per Compiler for execution_modules::execute()" << endl;
		perfReportOut << "=====================================================================" << endl;
		perfReportOut << endl;

		perfReportOut << setw(5) << "" << setw(20) << "Compiler" << setw(25) << "Execution Module" 
						<< setw(12) << "Call Count" << setw(20) << "Execution Time" << endl;
		perfReportOut << setw(5) << "" << setw(85) << "=================== ======================== =========== ===============" << endl;
		for(itr = data->compiler_timing.begin(); itr != data->compiler_timing.end(); itr++) {
			int count = (int) data->compiler_timing.count((*itr).first);
			perfReportOut << setw(5) << "" << setw(20) << (*itr).first;
			timeType total_time = 0;
			for (int i = 0; i < count; i++) {
				timer_data *info = (*itr).second;
				if (i == 0)
					perfReportOut << setw(25) << info->key << setw(12) << info->count << setw(20) << info->time << endl;
				else 
					perfReportOut << setw(25) << "" << setw(25) << info->key << setw(12) << info->count << setw(20) << info->time << endl;
				total_time += info->time;	
				itr++;
			}
			perfReportOut << setw(40) << "" << setw(22) << "Total execution time:" << setw(20) << total_time << endl;
			if (itr != data->compiler_timing.end())
				perfReportOut << setw(5) << "" << setw(85) << "------------------- ------------------------ ----------- ---------------" << endl;
			else 
				perfReportOut << setw(5) << "" << setw(85) << "========================================================================" << endl;
			itr--;		
		}
		perfReportOut << endl;
	}
	
	void performance_analysis::printRuleTimingsSummary() {
		map<string, rule_timing_info*>::iterator it;
		string rule_name;
		timeType max_rule_time = -1, min_rule_time = 100000, total_rule_time = 0;
		string max_rule_name, min_rule_name;
		int num_rules = 0, min_rule_size = 0, max_rule_size = 0;
		rule_timing_info rule_info;

		perfReportOut << endl;
		perfReportOut << "Rule Timings Summary:" << endl;
		perfReportOut << "===============================" << endl;
		perfReportOut << endl;
		
		perfReportOut << "Rule Name" << endl;
		perfReportOut << "Context Size\tTime Taken\tAverage Time" << endl;
		perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		perfReportOut << endl;
		
		for(it = data->rule_timings_table.begin(); it != data->rule_timings_table.end(); it++) { 
			rule_name = (*it).first;
		    rule_info.context_size = (*it).second->context_size;
		    rule_info.total_time = (*it).second->total_time;
			if (rule_info.total_time >= max_rule_time){
				max_rule_name = rule_name;
				max_rule_time = rule_info.total_time;
				max_rule_size = rule_info.context_size;
			}
			if (rule_info.total_time <= min_rule_time) {
				min_rule_name = rule_name;
				min_rule_time = rule_info.total_time;
				min_rule_size = rule_info.context_size;
			}
			total_rule_time += rule_info.total_time;
			num_rules++;
		}

		perfReportOut << "Average time spent on a rule       " << total_rule_time / num_rules << endl;
		perfReportOut << "Total number of rules executed     " << num_rules << endl;
		perfReportOut << "Total execution time of all rules  " << total_rule_time << endl;
		perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		perfReportOut << endl;
		 
		perfReportOut << "Maximum time spent on a rule       " << max_rule_time << endl;
		perfReportOut << "--------------------------------------" << endl;
		perfReportOut << "Rule: " << max_rule_name << endl;
		perfReportOut << max_rule_size << "\t" <<  max_rule_time << "\t" <<  max_rule_time / max_rule_size << endl;
		perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		perfReportOut << endl;
					
		perfReportOut << "Minimum time spent on a rule       " << min_rule_time << endl;
		perfReportOut << "--------------------------------------" << endl;
		perfReportOut << "Rule: " << min_rule_name << endl;
		perfReportOut << min_rule_size << "\t" << min_rule_time << "\t" <<  min_rule_time / min_rule_size << endl;
		perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		perfReportOut << endl;

		perfReportOut << endl;
	}

	void performance_analysis::printRuleTimingsDetails() {
		map<string, rule_timing_info*>::iterator it;
		string rule_name;
		rule_timing_info rule_info;

		perfReportOut << endl;
		perfReportOut << "Complete Rule Timings Analysis:" << endl;
		perfReportOut << "===============================================" << endl;
		perfReportOut << endl;
		
		perfReportOut << "Rule Name " << endl;
		perfReportOut << "Context Size\tTime Taken\tAverage Time" << endl;
		perfReportOut << endl;
		perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		perfReportOut << endl;
		
	    for(it = data->rule_timings_table.begin(); it != data->rule_timings_table.end(); it++) { 
			rule_name = (*it).first;
		    rule_info.context_size = (*it).second->context_size;
		    rule_info.total_time = (*it).second->total_time;

		    perfReportOut << rule_name << "\n";
			if( rule_info.context_size != 0) 
			      perfReportOut << rule_info.context_size  << "\t" << rule_info.total_time << "\t" << rule_info.total_time / rule_info.context_size << endl;
			else 
			      perfReportOut << rule_info.context_size  << "\t" << rule_info.total_time << "\t" << endl;
			perfReportOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	    }
		perfReportOut << endl;
	}
				
	execute_modules_timer::execute_modules_timer(executeP module, performance_analysis* analysis, string key) {
		wrapped = module;
		compiler_key = key;
		module_key = module->getName();
		// Tracking the time that each compiler spends in a specific execute_modules execute()
		timer_key = module_key + "-" + compiler_key;
		perf_analysis = analysis;
	}
	
	void execute_modules_timer::execute(fact_db &facts) {
		timer_token module_timer = perf_analysis->start_timer(timer_key);
		wrapped->execute(facts);
		perf_analysis->stop_timer(module_timer);
	}

	void execute_modules_timer::Print(std::ostream &s) const {
		wrapped->Print(s);
	}
	
	execute_modules_timer_factory::execute_modules_timer_factory(performance_analysis* analysis, string module_key){
		perf_analysis = analysis;
		key = module_key;
	}

	executeP execute_modules_timer_factory::decorate(executeP module) {
		return new execute_modules_timer(module, perf_analysis, key);
	}
} // end namespace Loci

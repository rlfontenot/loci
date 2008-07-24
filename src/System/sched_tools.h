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
#ifndef SCHED_TOOLS_H
#define SCHED_TOOLS_H

#include <map>
#include <vector>
#include <set>
#include <list>

#include <Tools/cptr.h>
#include <scheduler.h>
#include <Tools/digraph.h>
#include <fact_db.h>
#include <sched_db.h>
#include <execute.h>
#include <depend_graph.h>
#include <Map.h>

#ifdef PROFILE_CODE
#include <time.h>
#endif
#include "performance_analysis.h"
#include "sched_mlg.h"
using std::vector;

namespace Loci {
  void extract_rule_sequence(std::vector<rule> &rule_seq,
                             const std::vector<digraph::vertexSet> &v) ;
  void set_var_types(fact_db &facts, const digraph &dg, sched_db &scheds) ;
  rule make_rename_rule(variable new_name, variable old_name) ;

  class execute_rule : public execute_modules {
  protected:
    rule_implP rp ;
    rule rule_tag ; 
    sequence exec_seq ;
  public:
    execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds);
    execute_rule(rule fi, sequence seq, fact_db &facts, variable v, const storeRepP &p, const sched_db &scheds);
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
	virtual string getName() {return "execute_rule";};
  } ;

  class execute_rule_null : public execute_modules {
  protected:
    rule rule_tag ; 
  public:
    execute_rule_null(rule fi) : rule_tag(fi) {}
    virtual void execute(fact_db &facts) {}
    virtual void Print(std::ostream &s) const {s << rule_tag << " over empty sequence."<< endl ;}
	virtual string getName() {return "execute_rule_null";};
  } ;

  class dynamic_schedule_rule: public execute_modules {
    rule_implP rp ;
    variableSet inputs, outputs ;
    fact_db local_facts;
    rule_implP local_compute1;
    rule rule_tag ;
    entitySet pre_exec_set ;
    entitySet exec_set ;
    
  public:
    dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds) ;
    virtual ~dynamic_schedule_rule() ;
  
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
	virtual string getName() {return "dynamic_schedule_rule";};
  } ;
  
  class visitor ;
  
  class rule_compiler : public CPTR_type {
  public:
    ////////////////////
    virtual void accept(visitor& v) = 0 ;//method to accept a visitor
    ////////////////////
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) = 0 ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) = 0 ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) = 0;
  } ;

  typedef CPTR<rule_compiler> rule_compilerP ;
  typedef std::map<rule, rule_compilerP > rulecomp_map ;

  struct decomposed_graph {
    multiLevelGraph mlg ;
    digraph::vertexSet loops, recursive, conditional ;
    decomposed_graph(digraph dg, digraph::vertexSet sources,
                     digraph::vertexSet targets) ;
  } ;

  struct graph_compiler {
    rule_compilerP fact_db_comm ;
    rulecomp_map rule_process ;
    rule baserule ;
    /////////////////////
    std::vector<rule> super_rules ;
    /////////////////////
    graph_compiler(decomposed_graph &deco, variableSet initial_vars) ;
    ///////////////////
    // visit order
    void top_down_visit(visitor& v) ;
    void bottom_up_visit(visitor& v) ;
    void compile(fact_db& facts, sched_db& scheds,
                 const variableSet& given, const variableSet& target) ;
    ///////////////////
    void existential_analysis(fact_db &facts, sched_db &scheds) ;
    executeP execution_schedule(fact_db &facts, sched_db &scheds,
                                const variableSet& alloc) ;
  } ;
  
  struct comm_info {
    variable v ;
    int processor ;
    entitySet send_set ;
    sequence recv_set ;
  } ;

   class execute_param_red : public execute_modules {
     vector<variable> reduce_vars ;
     vector<rule> unit_rules ;
     MPI_Op create_join_op ;
     vector<CPTR<joiner> >join_ops ; 
   public:
     execute_param_red(vector<variable> reduce_vars, vector<rule> unit_rules,
                       vector<CPTR<joiner> > join_ops) ; 
     ~execute_param_red() ;
     virtual void execute(fact_db &facts) ;
     virtual void Print(std::ostream &s) const ;
 	virtual string getName() {return "execute_param_red";};
  } ;

  // experimental dynamic scheduling function
  void dynamic_scheduling(digraph& gr, fact_db& facts,
                          variableSet& given,
                          const variableSet& target) ;
  // this version will construct a graph internally and then
  // will throw it away before exiting the function. this is
  // safer than the above version.
  void
  dynamic_scheduling2(rule_db&, fact_db&, const variableSet&) ;
  // experimental dynamic mapping generation
  // in the stationary time level
  void stationary_relation_gen(rule_db&, fact_db&, const variableSet&) ;
  
  // experimental code to process static & dynamic constraints in
  // a unified way, this is the stage 1 --- mainly to compute the
  // static constraints and also to do some pre-process to those
  // dynamic ones
  variableSet
  constraint_process_stage1(rule_db&, fact_db&, const variableSet&) ;
  // stage2 --- generate new rule_db and setting up things for
  // dynamic constraints
  rule_db
  constraint_process_stage2(const rule_db&, fact_db&, const variableSet&) ;
}
#endif


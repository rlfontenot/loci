#ifndef SCHED_TOOLS_H
#define SCHED_TOOLS_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#include "sched_mlg.h"



namespace Loci {
  void extract_rule_sequence(std::vector<rule> &rule_seq,
                             const std::vector<digraph::vertexSet> &v) ;
  void set_var_types(fact_db &facts, const digraph &dg, sched_db &scheds) ;
  rule make_super_rule(variableSet sources, variableSet targets,
                       variable cond = variable()) ;
  rule make_rename_rule(variable new_name, variable old_name) ;

  class execute_rule : public execute_modules {
    rule_implP rp ;
    rule rule_tag ; 
    sequence exec_seq ;
    bool do_run ;
  public:
    execute_rule(rule fi, sequence seq, fact_db &facts, sched_db &scheds) ;
     execute_rule(rule fi, sequence seq, fact_db &facts, variable v, const storeRepP &p, sched_db &scheds) ;
    execute_rule(bool output_empty, rule fi, sequence seq, fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  
  class rule_compiler : public CPTR_type {
  public:
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
    graph_compiler(decomposed_graph &deco, variableSet initial_vars ) ;
    void existential_analysis(fact_db &facts, sched_db &scheds) ;
    executeP execution_schedule(fact_db &facts, sched_db &scheds, int nth) ;
  } ;
  
  struct comm_info {
    variable v ;
    int processor ;
    entitySet send_set ;
    sequence recv_set ;
  } ;

   class execute_param_red : public execute_modules {
     variable reduce_var ;
     rule unit_rule ;
     CPTR<joiner> join_op ; 
   public:
     execute_param_red(variable reduce_var, rule unit_rule,
                       CPTR<joiner> join_op) ; 
     virtual void execute(fact_db &facts) ;
     virtual void Print(std::ostream &s) const ;
   } ;
}
#endif


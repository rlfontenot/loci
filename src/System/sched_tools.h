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

#include "sched_mlg.h"
using std::vector;
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

class dynamic_schedule_rule: public execute_modules {
    rule_implP rp ;
    variableSet inputs, outputs ;
    fact_db local_facts;
    rule_implP local_compute1;
    // rule_implP local_compute2;
    rule rule_tag ;
    entitySet exec_set ;
   
     
  public:
    dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
   
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
                                const variableSet& alloc, int nth) ;
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
   } ;

   // experimental dynamic scheduling function
   void dynamic_scheduling(digraph& gr, fact_db& facts,
                           variableSet& given,
                           const variableSet& target) ;
}
#endif


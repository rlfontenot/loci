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
  void set_var_types(fact_db &facts, const digraph &dg) ;
  rule make_super_rule(variableSet sources, variableSet targets,
                       variable cond = variable()) ;
  rule make_rename_rule(variable new_name, variable old_name) ;

  class execute_rule : public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    bool do_run ;
  public:
    execute_rule(rule fi, sequence seq, fact_db &facts) ;
    execute_rule(rule fi, sequence seq, fact_db &facts, variable v, const storeRepP &p) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  
  class rule_compiler : public CPTR_type {
  public:
    virtual void set_var_existence(fact_db &facts) = 0 ;
    virtual void process_var_requests(fact_db &facts) = 0 ;
    virtual executeP create_execution_schedule(fact_db &facts) = 0;
  } ;

  typedef std::map<rule, CPTR<rule_compiler> > rulecomp_map ;

  struct decomposed_graph {
    multiLevelGraph mlg ;
    digraph::vertexSet loops, recursive, conditional ;
    decomposed_graph(digraph dg, digraph::vertexSet sources,
                     digraph::vertexSet targets) ;
  } ;

  struct graph_compiler {
    rulecomp_map rule_process ;
    rule baserule ;
    graph_compiler(decomposed_graph &deco) ;
    void existential_analysis(fact_db &facts) ;
    executeP execution_schedule(fact_db &facts, int nth) ;
  } ;
  struct proc_details {
    int processor ;
    entitySet send_set ;
    sequence recv_set ;
  } ;
  
  struct comm_info {
    variable v ;
    std::vector<proc_details> send_info ;
    std::vector<proc_details> recv_info ;
  } ;
  class execute_precomm : public execute_modules {
    std::list<comm_info> comm_list ;
    sequence exec_sequence ;
  public:
    execute_precomm(std::list<comm_info> plist , sequence seq, fact_db &facts) ; 
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  class execute_postcomm : public execute_modules {
    std::list<comm_info> comm_list ;
    sequence exec_sequence ;
  public:
    execute_postcomm(std::list<comm_info> clist , sequence seq, fact_db &facts) ; 
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
}
#endif

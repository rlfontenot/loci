#ifndef COMP_TOOLS_H
#define COMP_TOOLS_H

#include <fact_db.h>
#include <execute.h>
#include <Tools/digraph.h>

#include "sched_tools.h"
#include <vector>
#include <map>

namespace Loci {
  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts) ;
  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute) ;
  void existential_rule_analysis(rule f, fact_db &facts) ;
  entitySet process_rule_requests(rule f, fact_db &facts) ;

  void parallel_schedule(execute_par *ep,const entitySet &exec_set,
                         const rule &impl, fact_db &facts) ;
  std::vector<entitySet> partition_set(const entitySet &s,int nthreads) ;

  typedef std::map<variable,entitySet> vdefmap ;
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts) ;
  entitySet vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                                 entitySet compute) ;

  std::vector<digraph::vertexSet> schedule_dag(const digraph &g,
                                          digraph::vertexSet start_vertices = EMPTY,
                                          digraph::vertexSet only_vertices =
                                          interval(UNIVERSE_MIN,UNIVERSE_MAX)) ;



  class loop_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    CPTR<rule_compiler> calc(const rule &r) 
    {return rule_process[r] ;}

    digraph dag ;
  
    std::vector<rule> rule_schedule ;
    std::vector<rule> collapse ;
    std::vector<digraph::vertexSet> collapse_sched ;
    variable cond_var ;
    std::vector<rule> advance ;
    std::vector<digraph::vertexSet> advance_sched ;
    variableSet advance_vars ;
    time_ident tlevel ;
    variable output ;
    bool output_present ;
  public:
    loop_compiler(rulecomp_map &rp, digraph gin) ;
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  // rule compiler for rule with concrete implementation
  class impl_compiler : public rule_compiler {
    rule impl ;  // rule to implement

    // existential analysis info
    entitySet exec_seq ;
  public:
    impl_compiler(rule r)  { impl=r;}
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  // rule compiler for single rule recursion
  class impl_recurse_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    struct fcontrol {
      std::list<entitySet> control_list ;
      std::map<variable,entitySet> generated ;
      bool use_constraints ;
      entitySet nr_sources, constraints ;
      struct mapping_info {
        entitySet total ;
        std::vector<MapRepP> mapvec ;
        std::vector<variable> mapvar ;
        variable v ;
      } ;
      std::vector<mapping_info> recursion_maps, target_maps ;
    } ;
    fcontrol  control_set ;
    sequence fastseq ;
    std::vector<entitySet > par_schedule ;
  public:
    impl_recurse_compiler(rule r)
    { impl = r ;}
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class recurse_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    CPTR<rule_compiler> calc(const rule &r) 
    {return rule_process[r] ;}

    ruleSet recurse_rules ;
    struct fcontrol {
      std::list<entitySet> control_list ;
      std::map<variable,entitySet> generated ;
      bool use_constraints ;
      entitySet nr_sources, constraints ;
      struct mapping_info {
        entitySet total ;
        std::vector<MapRepP> mapvec ;
        variable v ;
      } ;
      std::vector<mapping_info> recursion_maps, target_maps ;
    } ;
    std::map<rule,fcontrol > control_set ;
  public:
    recurse_compiler(rulecomp_map &rp, ruleSet rs) : rule_process(rp)
    { recurse_rules = rs ; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class dag_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    CPTR<rule_compiler> calc(const rule &r) 
    {return rule_process[r] ;}
  
    digraph dag ;
    std::vector<rule> rule_schedule ;
    std::vector<digraph::vertexSet> dag_sched ;
  public:
    dag_compiler(rulecomp_map &rp, digraph gin) ;
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class conditional_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    CPTR<rule_compiler> calc(const rule &r) 
    {return rule_process[r] ;}

    digraph dag ;
    std::vector<rule> rule_schedule ;
    std::vector<digraph::vertexSet> dag_sched ;
    variable cond_var ;
  public:
    conditional_compiler(rulecomp_map &rp, digraph gin,
                           variable conditional) ;
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  // apply rule compiler 
  class apply_compiler : public rule_compiler {
    rule apply,unit_tag ;  // rule to applyement
  
    // existential analysis info
    entitySet exec_seq ;
    bool output_mapping ;
  public:
    apply_compiler(rule r, rule ut)
    { apply=r; unit_tag = ut ; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;





  class promote_compiler : public rule_compiler {
    rule r ;
  public:
    promote_compiler(rule rin)
    { r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  class generalize_compiler : public rule_compiler {
    rule r ;
  public:
    generalize_compiler(rule rin)
    { r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  
  class priority_compiler : public rule_compiler {
    rule r ;
  public:
    priority_compiler(rule rin)
    { r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  
}

#endif

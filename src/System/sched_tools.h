#ifndef SCHED_TOOLS_H
#define SCHED_TOOLS_H

#include <map>
#include <vector>
#include <set>
#include <list>

#include <scheduler.h>
#include <Tools/digraph.h>
#include <fact_db.h>
#include <execute.h>
#include <depend_graph.h>
#include <Map.h>

#ifdef PROFILE_CODE
#include <time.h>
#endif

#define DEVELOP

namespace Loci {
  void extract_rule_sequence(std::vector<rule> &rule_seq,
                             const std::vector<digraph::vertexSet> &v) ;
  void set_var_types(fact_db &facts, const digraph &dg) ;
  rule make_super_rule(variableSet sources, variableSet targets,
                       variable cond = variable()) ;
  rule make_rename_rule(variable new_name, variable old_name) ;
  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts) ;
  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute) ;

  digraph::vertexSet visit_vertices(digraph dg,digraph::vertexSet begin) ;
  digraph::vertexSet visit_vertices_exclusive(digraph dg, digraph::vertexSet begin) ;
  class decompose_graph ;
  class rule_compiler ;
  class loop_compiler ;

  class decompose_graph {
    friend class rule_compiler ;
    friend class loop_compiler ;
    variable start,finish ;
    struct reduce_info {
      rule unit_rule ;
      ruleSet apply_rules ;
    } ;
    std::map<variable,reduce_info> reduce_map ;
    variableSet reduce_vars ;
    std::map<rule,rule> apply_to_unit ;

    ruleSet looping_rules ;
    ruleSet conditional_rules ;
    variableSet conditional_vars ;
    std::map<variable,ruleSet> conditional_map ;
  
    rule top_level_rule ;
    struct super_node_info {
      digraph graph ;
      variableSet sources,targets ;
    } ;
    std::map<rule, super_node_info> supermap ;
    ruleSet supernodes ;
    rule create_supernode(digraph &g, digraph::vertexSet ns,
                          variable cond_var = variable()) ;
    ruleSet recursive_supernodes ;
    ruleSet looping_supernodes ;
  
    std::map<rule, rule_compiler *> rule_process ;
  
  public:
    decompose_graph(digraph dg,digraph::vertexSet sources, digraph::vertexSet targets) ;
    void existential_analysis(fact_db &facts) ;
    executeP execution_schedule(fact_db &facts, int nth) ;
  } ;

  class rule_compiler {
  protected:
    decompose_graph *graph_ref ;
  
    rule_compiler *calc(const rule &r) 
    {return graph_ref->rule_process[r] ;}
  public:
    virtual void set_var_existence(fact_db &facts) = 0 ;
    virtual void process_var_requests(fact_db &facts) = 0 ;
    virtual executeP create_execution_schedule(fact_db &facts) = 0;
  } ;

  class loop_compiler : public rule_compiler {
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
    loop_compiler(decompose_graph *gr, digraph gin) ;
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
    impl_compiler(decompose_graph *gr, rule r)  { graph_ref = gr; impl=r;}
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
    impl_recurse_compiler(decompose_graph *gr, rule r)
    { graph_ref = gr; impl = r ;}
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class recurse_compiler : public rule_compiler {
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
    recurse_compiler(decompose_graph *gr, ruleSet rs)
    { graph_ref = gr, recurse_rules = rs ; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class dag_compiler : public rule_compiler {
    digraph dag ;
    std::vector<rule> rule_schedule ;
    std::vector<digraph::vertexSet> dag_sched ;
  public:
    dag_compiler(decompose_graph *gr, digraph gin) ;
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;

  class conditional_compiler : public rule_compiler {
    digraph dag ;
    std::vector<rule> rule_schedule ;
    std::vector<digraph::vertexSet> dag_sched ;
    variable cond_var ;
  public:
    conditional_compiler(decompose_graph *gr, digraph gin,
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
    apply_compiler(decompose_graph *gr, rule r, rule ut)
    { graph_ref = gr; apply=r; unit_tag = ut ; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;





  class promote_compiler : public rule_compiler {
    rule r ;
  public:
    promote_compiler(decompose_graph *gr, rule rin)
    { graph_ref = gr, r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  class generalize_compiler : public rule_compiler {
    rule r ;
  public:
    generalize_compiler(decompose_graph *gr, rule rin)
    { graph_ref = gr, r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  
  class priority_compiler : public rule_compiler {
    rule r ;
  public:
    priority_compiler(decompose_graph *gr, rule rin)
    { graph_ref = gr, r = rin; }
    virtual void set_var_existence(fact_db &facts) ;
    virtual void process_var_requests(fact_db &facts) ;
    virtual executeP create_execution_schedule(fact_db &facts) ;
  } ;
  
  class error_compiler : public rule_compiler {
  public:
    error_compiler() {}
    virtual void set_var_existence(fact_db &facts)
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual void process_var_requests(fact_db &facts) 
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual executeP create_execution_schedule(fact_db &facts)
    { cerr << "Internal consistency error" << endl ; exit(-1);
    return executeP(0);}
  } ;


  class allocate_all_vars : public execute_modules {
  public:
    allocate_all_vars() { control_thread = true ; }
    virtual void execute(fact_db &facts) ;
    virtual void Print(ostream &s) const ;
  } ;


  void test_decompose(digraph dg, digraph::vertexSet sources,
                      digraph::vertexSet targets) ;

}

//#define VERBOSE
#endif

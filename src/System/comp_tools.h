#ifndef COMP_TOOLS_H
#define COMP_TOOLS_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <fact_db.h>
#include <execute.h>
#include <Tools/digraph.h>

#include "sched_tools.h"
#include <vector>
#include <map>
#include <deque>

#include <mpi.h>
using std::vector;

// for the mallinfo function
#include <malloc.h>

#ifdef LINUX
#define HAS_MALLINFO
#endif

#ifdef SPARC
#define HAS_MALLINFO
#endif

#ifdef SGI
#define HAS_MALLINFO
#endif

namespace Loci {

  bool rule_has_mapping_in_output(rule r);
  variableSet input_variables_with_mapping(rule r);
  variableSet input_variables(rule r);
  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts, sched_db &scheds) ;
  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute, sched_db &scheds) ;

  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var, sched_db &scheds) ;
 
  
  void existential_rule_analysis(rule f, fact_db &facts, sched_db &scheds) ;
  entitySet process_rule_requests(rule f, fact_db &facts, sched_db &scheds) ;
  
  void existential_applyrule_analysis(rule apply, fact_db &facts, sched_db &scheds) ;
  entitySet process_applyrule_requests(rule apply, rule unit_tag, fact_db &facts, sched_db &scheds) ;
  
  std::vector<std::pair<variable,entitySet> >
    barrier_existential_rule_analysis(variableSet vlst, fact_db &facts, sched_db &scheds) ;
  std::vector<std::pair<variable,entitySet> >
    send_ent_for_plist(variableSet vlst, fact_db &facts, sched_db &scheds);
  std::list<comm_info>
  barrier_process_rule_requests(variableSet vars, fact_db &facts, sched_db &scheds) ;

  entitySet send_requests(const entitySet& e, variable v, fact_db &facts,
                          std::list<comm_info> &clist) ;
  
  std::list<comm_info>
  put_precomm_info(std::vector<std::pair<variable,entitySet> > send_entities,
                   fact_db &facts) ;
  
  std::list<comm_info> sort_comm(std::list<comm_info> slist, fact_db &facts) ;
  
  void parallel_schedule(execute_par *ep,const entitySet &exec_set,
                         const rule &impl, fact_db &facts, sched_db &scheds) ;
  std::vector<entitySet> partition_set(const entitySet &s,int nthreads) ;
  
  void create_user_function(unsigned char* , unsigned char* , int*,
                            MPI_Datatype* ) ;
  
  typedef std::map<variable,entitySet> vdefmap ;
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts, sched_db &scheds) ;
  entitySet vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                                 entitySet compute, sched_db &scheds) ;

  std::vector<digraph::vertexSet> schedule_dag(const digraph &g,
                           digraph::vertexSet start_vertices = EMPTY,
                           digraph::vertexSet only_vertices =
                           interval(UNIVERSE_MIN,UNIVERSE_MAX)) ;

  class visitor ;
  
  class loop_compiler : public rule_compiler {
  public:
    int cid ; // id number of this compiler
    std::vector<rule_compilerP> collapse_comp ; 
    variable cond_var ;
    std::vector<rule_compilerP> advance_comp ;
    variableSet advance_vars ;
    time_ident tlevel ;
    variable output ;
    bool output_present ;
    variableSet all_loop_vars ;
    ////////////
    variableSet collapse_vars ; // the collapse variables
    std::vector<digraph::vertexSet> collapse_sched ;
    std::vector<digraph::vertexSet> advance_sched ;
    // the rotate list
    std::list<std::list<variable> > rotate_lists ;
    ////////////////
    // the internal loop graph and the rulecomp map
    digraph loop_gr ;
    digraph collapse_gr ;
    digraph advance_gr ;
    rulecomp_map rule_compiler_map ;
    ////////////////
  public:
    loop_compiler(rulecomp_map &rp, digraph gin, int id) ;
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  // rule compiler for rule with concrete implementation
  class impl_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    // existential analysis info
    entitySet exec_seq ;
  public:
    impl_compiler(rule r)  { impl=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
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

    std::vector<std::pair<variable,entitySet> > pre_send_entities ;
    std::list<comm_info> clist ;
  public:
    int cid ; // id number of this compiler
    impl_recurse_compiler(rule r, int id):cid(id)
      { impl = r ;}
    ruleSet get_rules() const
      {ruleSet rs; rs += impl; return rs ;}
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class recurse_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    rule_compilerP calc(const rule &r) 
    {return rule_process[r] ;}

    ruleSet recurse_rules ;
    variableSet recurse_vars ;
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
    std::list<std::vector<std::pair<variable,entitySet> > > recurse_send_entities ;
    std::map<variable,std::vector<std::list<comm_info> > > send_req_var ;
    std::list<std::list<comm_info> > recurse_clist ;
    std::list<std::list<comm_info> > recurse_plist ;
      

    std::vector<std::pair<variable,entitySet> > pre_send_entities ;
    std::list<comm_info> pre_clist ;
    std::list<comm_info> post_clist ;
    std::list<comm_info> pre_plist ;

  public:
    int cid ; //id number of this compiler
    recurse_compiler(rulecomp_map &rp, ruleSet rs, int id) : rule_process(rp),cid(id)
    {
      recurse_rules = rs ;
      for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) 
        recurse_vars += ri->targets() ;
    }
    ruleSet get_rules() const {return recurse_rules ;}
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class dag_compiler : public rule_compiler {
  public:
    int cid ; // id number of this compiler
    std::vector<rule_compilerP> dag_comp ;
    ////////////
    std::vector<digraph::vertexSet> dag_sched ;
    /////////////////
    digraph dag_gr ;
    rulecomp_map rule_compiler_map ;
    /////////////////
  public:
    dag_compiler(rulecomp_map &rp, digraph dag, int id) ;
    /////////////////
    virtual void accept(visitor& v) ;
    /////////////////
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class barrier_compiler : public rule_compiler {
    variableSet barrier_vars ;
    std::vector<std::pair<variable,entitySet> > send_entities ;
    std::list<comm_info> clist ;
    std::list<comm_info> plist ;
  public:
    barrier_compiler(variableSet &vars)
      : barrier_vars(vars) {}
    /////////////
    virtual void accept(visitor& v) {}
    /////////////
    virtual void set_var_existence(fact_db &facts,  sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class singleton_var_compiler : public rule_compiler {
    variableSet barrier_vars ;
  public:
    singleton_var_compiler(variableSet &vars)
      :  barrier_vars(vars) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class reduce_param_compiler : public rule_compiler {
    vector<variable> reduce_vars ;
    vector<rule> unit_rules ;
    vector<CPTR<joiner> > join_ops ;
  public:
    reduce_param_compiler(const vector<variable> &v, const vector<rule> &ur,
                          vector<CPTR<joiner> >&jop) :
      reduce_vars(v), unit_rules(ur), join_ops(jop) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class reduce_store_compiler : public rule_compiler {
    variable reduce_var ;
    rule unit_rule ;
    CPTR<joiner> join_op ;

    std::list<comm_info> rlist ;  // reduction communication
    std::list<comm_info> clist ;  // comm from owner to requester
  public:
    reduce_store_compiler(const variable &v, const rule &ur,
                          CPTR<joiner> &jop) :
      reduce_var(v), unit_rule(ur), join_op(jop) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class conditional_compiler : public rule_compiler {
  public:
    int cid ; //id number of this compiler
    std::vector<rule_compilerP> dag_comp ;
    variable cond_var ;
    /////////
    std::vector<digraph::vertexSet> dag_sched ;
    ///////////
    digraph cond_gr ;
    rulecomp_map rule_compiler_map ;
    ///////////
  public:
    conditional_compiler(rulecomp_map &rp, digraph gin,
                         variable conditional, int id) ;
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
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
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class promote_compiler : public rule_compiler {
    rule r ;
  public:
    promote_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  class generalize_compiler : public rule_compiler {
    rule r ;
  public:
    generalize_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class priority_compiler : public rule_compiler {
    rule r ;
  public:
    priority_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class execute_msg : public execute_modules {
    std::string msg ;
  public:
    execute_msg(std::string m) : msg(m) {}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;

  struct send_var_info {
    variable v ;
    entitySet set ;
    send_var_info(variable iv, const entitySet &iset) : v(iv),set(iset) {}
  } ;
  struct recv_var_info {
    variable v ;
    sequence seq ;
    recv_var_info(variable iv, const sequence &iseq) : v(iv),seq(iseq) {}
  } ;

  class execute_comm : public execute_modules {
    std::vector<std::pair<int,std::vector<send_var_info> > > send_info ;
    std::vector<std::pair<int,std::vector<recv_var_info> > > recv_info ;
    std::vector<std::vector<storeRepP> > send_vars ;
    std::vector<std::vector<storeRepP> > recv_vars ;
    int *maxr_size, *maxs_size, *r_size, *s_size, *recv_sizes ;
    unsigned char **recv_ptr , **send_ptr ;
    MPI_Request *request;
    MPI_Status *status ;
  public:
    execute_comm(std::list<comm_info> &plist, fact_db &facts) ;
    ~execute_comm() ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ; 
  
  class allocate_var_compiler : public rule_compiler {
    variableSet allocate_vars ;
  public:
    allocate_var_compiler(const variableSet& vars)
      : allocate_vars(vars) {}
    allocate_var_compiler(const variable& var)
    {allocate_vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class free_var_compiler : public rule_compiler {
    variableSet free_vars ;
  public:
    free_var_compiler(const variableSet& vars)
      : free_vars(vars) {}
    free_var_compiler(const variable& var)
    {free_vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class execute_memProfileAlloc : public execute_modules {
    variableSet vars ;
  public:
    execute_memProfileAlloc(const variableSet& vars): vars(vars) {}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    // memory profile function
    int currentMem(void) {
#ifdef HAS_MALLINFO
      struct mallinfo info = mallinfo() ;
      return info.arena+info.hblkhd ;
#else
      cerr << "memProfile not implemented" << endl;
      return 0 ;
#endif
    }    
  } ;

  class execute_memProfileFree : public execute_modules {
    variableSet vars ;
  public:
    execute_memProfileFree(const variableSet& vars) : vars(vars) {}
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
    // memory profile function
    int currentMem(void) {
#ifdef HAS_MALLINFO
      struct mallinfo info = mallinfo() ;
      return info.arena+info.hblkhd ;
#else
      cerr << "memProfile not implemented" << endl ;
      return 0 ;
#endif
    }    
  } ;

  class memProfileAlloc_compiler : public rule_compiler {
    variableSet vars ;
  public:
    memProfileAlloc_compiler(const variableSet& vars): vars(vars) {}
    memProfileAlloc_compiler(const variable& var)
    {vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class memProfileFree_compiler : public rule_compiler {
    variableSet vars ;
  public:
    memProfileFree_compiler(const variableSet& vars): vars(vars) {}
    memProfileFree_compiler(const variable& var)
    {vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class chomp_compiler: public rule_compiler {
  public:
    digraph chomp_graph ;
    variableSet chomp_vars ;
    std::vector<digraph::vertexSet> chomp_sched ;
    std::vector<std::pair<rule,rule_compilerP> > chomp_comp ;
    std::deque<entitySet> rule_seq ;
    std::map<rule,rule> apply2unit ;
  public:
    chomp_compiler(const digraph& cgraph,const variableSet& cvars,
                   const std::map<rule,rule>& a2u) ;
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts,sched_db& scheds) ;
    virtual void process_var_requests(fact_db& facts,sched_db& scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,sched_db& scheds) ;
  } ;

  // compiler to handle the dynamic scheduling stuff
  class dynamic_compiler: public rule_compiler {
  public:
    int cid ;
    std::vector<rule_compilerP> comp ;
    std::vector<digraph::vertexSet> sched ;
    digraph gr ;
    rulecomp_map rule_compiler_map ;
    variableSet all_reduce_vars ;
    std::map<variable,std::pair<rule,CPTR<joiner> > > reduce_info ;
  public:
    dynamic_compiler(rulecomp_map& rp, const digraph& g, int id) ;
    virtual void accept(visitor& v) {}
    void collect_reduce_info() ;
    void schedule() ;
    void compile() ;
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) ;
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;
  
  // rule compiler for constraint rules
  class constraint_compiler : public rule_compiler {
    rule constraint_rule ;  // the constraint rule
  public:
    constraint_compiler(rule r)  { constraint_rule=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;


}

#endif

#include "sched_tools.h"
#include "comp_tools.h"
#include "dist_tools.h"
#include "visitor.h"

#include <Tools/stream.h>

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <set>
using std::set ;
#include <string>
using std::string ;
#include <sstream>
using std::stringstream ;

//#define HACK ; 

namespace Loci {
  // Loci option flags
  extern double total_memory_usage ;
  extern bool show_decoration ;
  extern bool use_dynamic_memory ;
  extern bool show_dmm_verbose ;
  extern bool use_chomp ;
  extern bool show_chomp ;
  extern bool chomp_verbose ;
  extern int chomping_size ;
  extern bool profile_memory_usage ;
  extern bool memory_greedy_schedule ;

  // memory profiling counters
  extern double LociAppPeakMemory ;
  extern double LociAppAllocRequestBeanCounting ;
  extern double LociAppFreeRequestBeanCounting ;
  extern double LociAppPeakMemoryBeanCounting ;
  extern double LociAppLargestAlloc ;
  extern variable LociAppLargestAllocVar ;
  extern double LociAppLargestFree ;
  extern variable LociAppLargestFreeVar ;
  extern double LociAppPMTemp ;
  namespace {
    // used to pre-process preallocation memory profiling
    variableSet LociRecurrenceVarsRealloc ;
  }
  
  class error_compiler : public rule_compiler {
  public:
    error_compiler() {}
    virtual void accept(visitor& v)
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds)
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) 
    { cerr << "Internal consistency error" << endl ; exit(-1);}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds)
    { cerr << "Internal consistency error" << endl ; exit(-1);
    return executeP(0);}
  } ;

  int get_supernode_num(int rid) {
    rule r(rid) ;
    string rqualifier = r.get_info().qualifier() ;
    string head = rqualifier.substr(0,2) ;
    if(head != "SN") return -1 ;
    
    string number = rqualifier.substr(2,rqualifier.size()-2) ;
    stringstream ss ;
    ss << number ;
    int ret ;
    ss >> ret ;

    return ret ;
  }
 
  graph_compiler::graph_compiler(decomposed_graph &deco,variableSet
				 initial_vars) {
    fact_db_comm = new barrier_compiler(initial_vars) ;
    multiLevelGraph &mlg = deco.mlg ;
    vector<int> levels ;
    baserule = rule(mlg.toplevel) ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }
    for(vector<int>::reverse_iterator ii=levels.rbegin();
        ii!=levels.rend();
        ++ii) {
      multiLevelGraph::subGraph *p = mlg.find(*ii) ;
      ruleSet level_rules = extract_rules(p->graph_v-mlg.subgraphs) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;
      gr.remove_dangling_vertices() ;

      digraph grt = gr.transpose() ;

      // Collect information on unit rules on this level
      map<rule,rule> apply_to_unit ;
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
          variable unit_var = *(ri->targets().begin()) ;
          ruleSet apply_rules = extract_rules(grt[unit_var.ident()]) ;
          for(ruleSet::const_iterator rii=apply_rules.begin();
              rii!=apply_rules.end();
              ++rii) {
            apply_to_unit[*rii] = *ri ;
          }
        }
      }

      // Assign compilers to rules on this level
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->type() == rule::INTERNAL) {
          if(ri->get_info().qualifier() == "promote")
            rule_process[*ri] = new promote_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "generalize")
            rule_process[*ri] = new generalize_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "priority")
            rule_process[*ri] = new priority_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "ALLOCATE")
            rule_process[*ri] = new allocate_var_compiler(ri->targets()) ;
          else
            rule_process[*ri] = new error_compiler ;
        } else {
          if(ri->get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
            rule_process[*ri] = new impl_compiler(*ri) ;
          else
            rule_process[*ri] = new apply_compiler(*ri,apply_to_unit[*ri]) ;
        }
      }

      // We have created compilers for every rule on this level, now
      // we create a compiler for this supernode

      rule snrule = rule(*ii) ;
      // push each snrule into the vector super_rules for later use
      super_rules.push_back(snrule) ;
      int id = get_supernode_num(*ii) ;

      if(deco.loops.inSet(*ii)) {
        // Looping supernode
        rule_process[snrule] = new loop_compiler(rule_process,gr,id) ;
      } else if(deco.recursive.inSet(*ii)) {
        // recursive supernode
        ruleSet recurse_rules = extract_rules(p->graph_v) ;
        if(recurse_rules.size() == 1 &&
           recurse_rules.begin()->get_info().desc.constraints.size() == 0
           && (MPI_processes ==1 || (*recurse_rules.begin()).get_rule_implP()->is_relaxed()))
          // Single rule recursion
          rule_process[snrule] =
            new impl_recurse_compiler(*(recurse_rules.begin()),id) ;
        else
          // Multi-rule recursion
          rule_process[snrule] =
            new recurse_compiler(rule_process, recurse_rules, id) ;
      } else if(deco.conditional.inSet(*ii)) {
        // conditional supernode
        variable cond_var = *(snrule.get_info().desc.conditionals.begin()) ;
        rule_process[snrule] =
          new conditional_compiler(rule_process,gr,cond_var,id) ;
      } else {
        // DAG supernode
        rule_process[snrule] = new dag_compiler(rule_process,gr,id) ;
      }
    } 
  }

  void graph_compiler::top_down_visit(visitor& v) {
    for(vector<rule>::reverse_iterator ii=super_rules.rbegin();
        ii!=super_rules.rend();++ii) {
      rulecomp_map::iterator ri ;
      ri = rule_process.find(*ii) ;
      FATAL(ri == rule_process.end()) ;
      (ri->second)->accept(v) ;
    }
  }
  
  void graph_compiler::bottom_up_visit(visitor& v) {
    for(vector<rule>::const_iterator ii=super_rules.begin();
        ii!=super_rules.end();++ii) {
      rulecomp_map::iterator ri ;
      ri = rule_process.find(*ii) ;
      FATAL(ri == rule_process.end()) ;
      (ri->second)->accept(v) ;
    }
  }
  
  fact_db *exec_current_fact_db = 0 ;

  void graph_compiler::compile(fact_db& facts,sched_db& scheds,
                               const variableSet& given,
                               const variableSet& target) {
    /***********************************
    1. unordered visitor phase (1..n)
         top -> down
         bottom -> up
       graph editing phase (1..n)

    2. ordering phase (1..n)

    3. assembly phase.
    ************************************/
    // timing variables
    double dst1=0,det1=0,dst2=0,det2=0,cst=0,cet=0 ;
    dst1 = MPI_Wtime() ;
    // must do this visitation
    // the loop rotate_lists is computed here
    rotateListVisitor rotlv(scheds) ;
    top_down_visit(rotlv) ;

    // get all the recurrence information
    // i.e. the generalize, promote, priority and rename info
    recurInfoVisitor recv ;
    top_down_visit(recv) ;

    // set the reallocated recurrence vars for
    // memory profiling preallocation
    if(!use_dynamic_memory && profile_memory_usage) {
      variableSet sources = variableSet(recv.get_recur_source_vars()
                                        - recv.get_recur_target_vars()) ;
      map<variable,variableSet> s2t_table = recv.get_recur_vars_s2t() ;
      LociRecurrenceVarsRealloc +=
        get_recur_target_for_vars(sources,s2t_table) ;
    }
    
    // get the input variables
    // get all the variables in the multilevel graph
    getAllVarVisitor allvarV ;
    top_down_visit(allvarV) ;
    variableSet input = variableSet(given & allvarV.get_all_vars_ingraph()) ;

    // check preconditions of recurrence
    check_recur_precondition(recv,input) ;

    // get the unit apply info (reduce info)
    unitApplyMapVisitor reduceV ;
    top_down_visit(reduceV) ;

    det1 = MPI_Wtime() ;
    if(use_dynamic_memory) {
      if(Loci::MPI_rank == 0)
        cout << "USING DYNAMIC MEMORY MANAGEMENT" << endl ;

      if(use_chomp) {
        if(Loci::MPI_rank == 0) 
          cout << "USING CHOMPING"
               << " (chomping size: "
               << chomping_size << "KB)"
               << endl ;
        
        cst = MPI_Wtime() ;
        chompPPVisitor cppv(facts,
                            rotlv.get_rotate_vars_table(),
                            rotlv.get_loop_shared_table(),
                            variableSet(recv.get_rename_source_vars() +
                                        recv.get_rename_target_vars())
                            ) ;
        top_down_visit(cppv) ;
        
        chompRuleVisitor crv(cppv.get_good_vars(),
                             cppv.get_bad_vars(),
                             reduceV.get_apply2unit()) ;
        top_down_visit(crv) ;
        cet = MPI_Wtime() ;

        if(show_chomp)
          crv.visualize(cout) ;
        if(chomp_verbose)
          crv.summary(cout) ;

        dagCheckVisitor dagcV1(true) ;
        top_down_visit(dagcV1) ;

      }
      
      dst2 = MPI_Wtime() ;
      // get inter/intra supernode information
      snInfoVisitor snv ;
      top_down_visit(snv) ;
      
      // get all untyped variables
      unTypedVarVisitor untypevarV(recv.get_recur_vars_s2t(),
                                   recv.get_recur_vars_t2s(),
                                   input) ;
      top_down_visit(untypevarV) ;

      // get the representitive variable for each rename target cluster
      variableSet only_targets = recv.get_rename_target_vars() ;
      only_targets -= recv.get_rename_source_vars() ;
      
      variableSet cluster_remaining ;
      cluster_remaining = pick_rename_target(recv.get_rename_s2t(),
                                             recv.get_rename_t2s(),
                                             only_targets
                                             ) ;

      // get representitive for each promoted variable cluster
      promotePPVisitor pppV(recv.get_promote_t2s(),
                            recv.get_promote_s2t(),
                            recv.get_promote_source_vars(),
                            recv.get_promote_target_vars(),
                            snv.get_graph_sn(),
                            snv.get_cond_sn(),input
                            ) ;
      top_down_visit(pppV) ;

      // reserved variableSet for the deleteInfoVisitor
      variableSet delInfoV_reserved ;

      variableSet promoted_rep =
        variableSet(pppV.get_rep() - recv.get_rename_source_vars()) ;

      delInfoV_reserved += recv.get_recur_source_vars() ;
      delInfoV_reserved += pppV.get_remaining() ;
      delInfoV_reserved -= promoted_rep ;
      delInfoV_reserved += input ;
      delInfoV_reserved += get_recur_target_for_vars(input,
                                                     recv.get_recur_vars_s2t()
                                                     ) ;
      delInfoV_reserved += target ;
      delInfoV_reserved += untypevarV.get_untyped_vars() ;
      delInfoV_reserved += cluster_remaining ;
      
      // compute how to do allocation
      allocInfoVisitor aiv(snv.get_graph_sn(),
                           recv.get_recur_target_vars(),
                           snv.get_loop_sn(),
                           rotlv.get_rotate_vars_table(),
                           rotlv.get_loop_shared_table(),
                           untypevarV.get_untyped_vars()
                           ) ;
      top_down_visit(aiv) ;
      
      // compute how to do deletion
      deleteInfoVisitor div(get_loop_alloc_table(aiv.get_alloc_table(),
                                                 snv.get_subnode_table(),
                                                 snv.get_loop_sn(),
                                                 snv.get_graph_sn(),
                                                 recv.get_recur_vars_s2t()),
                            recv.get_recur_vars_t2s(),
                            recv.get_recur_vars_s2t(),
                            snv.get_graph_sn(),
                            get_parentnode_table(snv.get_subnode_table()),
                            snv.get_loop_col_table(),
                            snv.get_loop_sn(),
                            snv.get_cond_sn(),
                            rotlv.get_rotate_vars_table(),
                            rotlv.get_loop_shared_table(),
                            delInfoV_reserved
                            ) ;
      top_down_visit(div) ;
      
      // decorate the graph to insert allocation rules
      allocGraphVisitor agv(aiv.get_alloc_table(),
                            snv.get_loop_sn(),
                            rotlv.get_rotate_vars_table(),
                            recv.get_priority_s2t(),
                            recv.get_priority_source_vars()
                            ) ;
      top_down_visit(agv) ;

      if(profile_memory_usage) {
        // decorate the graph to include
        // the memory profiling alloc compiler
        memProfileAllocDecoVisitor mpadv(aiv.get_alloc_table(),
                                         rotlv.get_rotate_vars_table()
                                         ) ;
        top_down_visit(mpadv) ;
      }
      
      // decorate the graph to insert deletion rules
      deleteGraphVisitor dgv(div.get_delete_table(),
                             div.get_recur_source_other_rules()) ;
      top_down_visit(dgv) ;
      det2 = MPI_Wtime() ;
      
      // check if the decorated graphs are acyclic
      dagCheckVisitor dagcV(true) ;
      top_down_visit(dagcV) ;

      if(use_chomp) {
        // compile all the chomp_compilers
        compChompVisitor compchompv ;
        top_down_visit(compchompv) ;
      }

      if(show_dmm_verbose) {
        // the output stream
        ostream& os = cout ;
        //print out the loop rotate list variables
        os << endl ;
        os << "---------loop rotate list variables-----------" << endl ;
        os << rotlv.get_rotate_vars_table() ;
        //print out the loop shared variables
        os << endl ;
        os << "---------loop shared variables-----------" << endl ;
        os << rotlv.get_loop_shared_table() ;
        os << endl ;
        //print out the untyped variables
        os << "----------untyped variables--------------" << endl ;
        os << untypevarV.get_untyped_vars() ;
        os << endl ;

        os << endl ;
        // get all typed variables in the fact database
        variableSet vars = facts.get_typed_variables() ;
        variableSet allgraphv = allvarV.get_all_vars_ingraph() ;
        // the difference
        variableSet typed_not_in_graph = variableSet(vars - allgraphv) ;
        if(typed_not_in_graph != EMPTY) {
          os << "These are typed variables in fact database, "
             << "but are not in the multilevel graph: "
             << typed_not_in_graph << endl << endl ;
        }
        
        variableSet allocated_vars ;
        map<int,variableSet> alloc_table = aiv.get_alloc_table() ;
        map<int,variableSet>::const_iterator miter ;
        for(miter=alloc_table.begin();miter!=alloc_table.end();++miter)
          allocated_vars += miter->second ;
        vars -= allocated_vars ;
        vars -= recv.get_recur_target_vars() ;
        // input variables
        rulecomp_map::iterator ri ;
        ri = rule_process.find(*super_rules.rbegin()) ;
        FATAL(ri == rule_process.end()) ;
        dag_compiler* dc = dynamic_cast<dag_compiler*>(&(*(ri->second))) ;
        digraph gr = dc->dag_gr ;
        digraph grt = gr.transpose() ;
        
        digraph::vertexSet input_vars ;
        digraph::vertexSet all_vertices = grt.get_all_vertices() ;  
        for(digraph::vertexSet::const_iterator vi=all_vertices.begin();
            vi!=all_vertices.end();++vi) {
          if(grt[*vi] == EMPTY)
            input_vars += *vi ;
        }
        vars -= extract_vars(input_vars) ;
        vars &= allgraphv ;
        
        if(vars != EMPTY) {
          os << "These are typed variables in the multilevel "
             << "graph, but are not scheduled to be allocated: "
             << vars << endl << endl ;
        }
        
        // report some allocation and deletion info
        allocDeleteStat adstat(aiv.get_alloc_table(),
                               div.get_delete_table(),
                               recv.get_recur_vars_t2s(),
                               recv.get_recur_vars_s2t(),
                               recv.get_recur_source_vars(),
                               recv.get_recur_target_vars()
                               ) ;
        adstat.report(os) ;
        os << endl ;

        allocDelNumReportVisitor adnrv(cout) ;
        top_down_visit(adnrv) ;
        os << endl ;
      } // end of if(show_dmm_verbose)
      
    } // end of if(use_dynamic_memory)
    
    // visualize each decorated graph
    if(show_decoration) {
      graphVisualizeVisitor visV ;
      top_down_visit(visV) ;
    }
    
    //orderVisitor ov ;    
    //bottom_up_visit(ov) ;
    if(!memory_greedy_schedule) {
      simLazyAllocSchedVisitor lazyschedv ;
      top_down_visit(lazyschedv) ;
    } else {
      if(Loci::MPI_rank == 0)
        cout << "memory greedy scheduling" << endl ;
      memGreedySchedVisitor mgsv(facts) ;
      top_down_visit(mgsv) ;
    }
    
    assembleVisitor av(reduceV.get_all_reduce_vars(),
                       reduceV.get_reduceInfo());
    bottom_up_visit(av) ;

    exec_current_fact_db = &facts ;

    if(use_dynamic_memory)
      Loci::debugout << "Time taken for dmm graph decoration = "
                     << (det1-dst1) + (det2-dst2) << " seconds " << endl ;
    if(use_chomp)
      Loci::debugout << "Time taken for chomping subgraph searching = "
                     << cet-cst << " seconds " << endl ;
  }
  
  void graph_compiler::existential_analysis(fact_db &facts, sched_db &scheds) {
    fact_db_comm->set_var_existence(facts, scheds) ;
    (rule_process[baserule])->set_var_existence(facts, scheds) ;
    variableSet var_requests = baserule.targets() ;
    variableSet::const_iterator vi ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi) ;
      scheds.variable_request(*vi,vexist) ;
    }
    (rule_process[baserule])->process_var_requests(facts, scheds) ;
    fact_db_comm->process_var_requests(facts, scheds) ;
  }
  
  class allocate_all_vars : public execute_modules {
    variableSet vars ;
    bool is_alloc_all ;
    std::map<variable,entitySet> v_requests, v_existence ;
  public:
    allocate_all_vars() { control_thread = true ; }
    allocate_all_vars(fact_db &facts, sched_db &scheds,
                      const variableSet& alloc,
                      bool is_alloc_all) ;
    void fill_in_requests(fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts) ;
    virtual void Print(std::ostream &s) const ;
  } ;
  
  
  //fact_db *exec_current_fact_db = 0 ;
  
  allocate_all_vars::allocate_all_vars(fact_db &facts, sched_db &scheds,
                                       const variableSet& alloc,
                                       bool is_alloc_all)
    :vars(alloc),is_alloc_all(is_alloc_all) {
    control_thread = true ;
    fill_in_requests(facts, scheds) ;
  }
  void allocate_all_vars::fill_in_requests(fact_db &facts, sched_db &scheds) {
    //    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi,vii ;
    
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      variableSet aliases = scheds.get_aliases(*vi) ;
      entitySet requests, existence ;
      for(vii=aliases.begin();vii!=aliases.end();++vii) {
	existence += scheds.variable_existence(*vii) ;
	requests += scheds.get_variable_requests(*vii) ;
      }
      v_requests[*vi] = requests ;
      v_existence[*vi] = existence ;
    }
  }
  
  void allocate_all_vars::execute(fact_db &facts) {
    //exec_current_fact_db = &facts ;
    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi ;  
    double total_size = 0 ;
    entitySet dom, total, unused ;
    double total_wasted = 0 ;
    bool facts_distributed = facts.is_distributed_start() ;
    
    //#define DIAGNOSTICS
#ifdef DIAGNOSTICS
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->RepType() == Loci::STORE) {
	dom = v_requests[*vi] ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += srp->pack_size(dom) ; 
	total_wasted += srp->pack_size(unused) ;
      }
    }
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->RepType() == Loci::STORE) {
	dom = v_requests[*vi] ;
	double size, wasted_space ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	size = srp->pack_size(dom) ;
	wasted_space = srp->pack_size(unused) ;
	Loci::debugout << " ****************************************************" << endl ;
	Loci::debugout << " Total_size = " << total_size << endl ;
	Loci::debugout << "Variable = "  << *vi << endl ;
	Loci::debugout << "Domain = " << dom << endl ;
	Loci::debugout << "Size allocated = " << size << endl ; 
	if(facts.isDistributed() )  {
	  Loci::fact_db::distribute_infoP d ;
	  d   = facts.get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = srp->pack_size(clone) ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	  Loci::debugout << " My_entities = " << my_entities << endl ;
	  Loci::debugout << " Clone entities = " << clone << endl ;
	  Loci::debugout << "Memory required for the  clone region  = " << clone_size << endl ;
	  Loci::debugout << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  Loci::debugout << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	}
	Loci::debugout << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	Loci::debugout << "----------------------------------------------------" << endl;
	Loci::debugout << "Total wasted size = " << total_wasted << endl ;
	Loci::debugout << "Unused entities = " << unused << endl ;
	Loci::debugout << "Wasted space = " << wasted_space << endl ;
	Loci::debugout << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	Loci::debugout << " ***************************************************" << endl << endl << endl ;
      }
    }
#endif

    //#define HACK
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
#ifdef HACK
      entitySet alloc_dom = v_existence[*vi]  + srp->domain() ;
#else
      entitySet alloc_dom = v_requests[*vi] + srp->domain() ;
#endif
      if(srp->domain() == EMPTY) {
	srp->allocate(alloc_dom) ;
      }
      else {
	if(srp->RepType() == Loci::STORE) {
	  entitySet tmp = interval(alloc_dom.Min(), alloc_dom.Max()) ;
	  if(tmp.size() >= 2*srp->domain().size())
	    Loci::debugout << "Variable = " << *vi << "  more than twice the space allocated :  allocated over " << alloc_dom << " size = " << tmp.size()  << "  while domain is only  " << srp->domain() << " size = " << srp->domain().size() << endl ;
	  if(alloc_dom != srp->domain()) {
	    Loci::debugout << "reallocating " << *vi << "  over  " << alloc_dom << " initially it was over  " << srp->domain() << endl ;
	    srp->allocate(alloc_dom) ;
	  }
	} 
      }
      // do memory profiling
      if(profile_memory_usage) {
        // we only do profiling for parallel start
        // since only in parallel, input vars are resized
        if(facts_distributed) {
          int packsize = srp->pack_size(alloc_dom) ;
          LociAppAllocRequestBeanCounting += packsize ;
          LociAppPMTemp += packsize ;
          if(LociAppPMTemp > LociAppPeakMemoryBeanCounting)
            LociAppPeakMemoryBeanCounting = LociAppPMTemp ;
          if(packsize > LociAppLargestAlloc) {
            LociAppLargestAlloc = packsize ;
            LociAppLargestAllocVar = *vi ;
          }
        }
      }
      
    }
    total_memory_usage = total_size + total_wasted ;
  }
  
  void allocate_all_vars::Print(std::ostream &s) const {
    if(is_alloc_all)
      s << "allocate all variables" << endl ;
    else
      s << "reallocate all given variables" << endl ;
  }


  executeP graph_compiler::
  execution_schedule(fact_db &facts, sched_db &scheds,
                     const variableSet& alloc, int nth) {

    CPTR<execute_list> schedule = new execute_list ;

    if(!use_dynamic_memory)
      schedule->append_list(new allocate_all_vars(facts,scheds,facts.get_typed_variables(),true)) ;
    else
      if(facts.is_distributed_start())
        if((MPI_processes > 1))
          schedule->append_list(new allocate_all_vars(facts,scheds,alloc,false)) ;

    schedule->append_list(new execute_create_threads(nth)) ;
    schedule->append_list(fact_db_comm->create_execution_schedule(facts, scheds));
    executeP top_level_schedule = (rule_process[baserule])->
      create_execution_schedule(facts, scheds) ;
    if(top_level_schedule == 0) 
      return executeP(0) ;
    schedule->append_list(top_level_schedule) ;
    schedule->append_list(new execute_destroy_threads) ;
    if(!use_dynamic_memory)
      if(profile_memory_usage) {
        variableSet profile_vars = facts.get_typed_variables() ;
        // given variables don't need to be profiled
        profile_vars -= alloc ;
        // we also need to exclude those recurrence variables
        profile_vars -= LociRecurrenceVarsRealloc ;
        
        schedule->append_list(new execute_memProfileAlloc(profile_vars)) ;
      }
    
    return executeP(schedule) ;
  }
}

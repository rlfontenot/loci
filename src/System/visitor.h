#ifndef VISITOR_H
#define VISITOR_H

#include <vector>
#include <map>
#include <set>
#include <Tools/intervalSet.h>
#include <Tools/digraph.h>
#include "sched_tools.h"

namespace Loci {

  class loop_compiler ;
  class dag_compiler ;
  class conditional_compiler ;
  
  class visitor {
  public:
    virtual ~visitor() {}
    virtual void visit(loop_compiler&) = 0 ;
    virtual void visit(dag_compiler&) = 0 ;
    virtual void visit(conditional_compiler&) = 0 ;
  } ;

  class orderVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    std::vector<digraph::vertexSet>
    order_dag(const digraph &g,
              digraph::vertexSet start_vertices = EMPTY,
              digraph::vertexSet only_vertices =
              interval(UNIVERSE_MIN,UNIVERSE_MAX)
              ) ;
  } ;

  class assembleVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    void compile_dag_sched(std::vector<rule_compilerP> &dag_comp,
                           const std::vector<digraph::vertexSet> &dag_sched,
                           const rulecomp_map& rcm,
                           const digraph& dag) ;
  } ;

  // obsolete visitor, need to be removed later
  class topDownInfoVisitor: public visitor {
  public:
    topDownInfoVisitor():inter_node_vars(EMPTY) {}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    void print_info(std::ostream& os) ;
    
    std::map<variable,int> get_ownership_table()
    {return ownership_table ;}
    std::map<variable,variable> get_recur_table()
    {return recurrence_vars_table ;}
    variableSet get_keep()
    {return keep_vars ;}
  protected:
    std::map<variable,int> ownership_table ;
    std::map<variable,variable> recurrence_vars_table ;
    variableSet keep_vars ;
  private:
    void visit_gr(const digraph& gr,int id) ;
    variableSet inter_node_vars ;
  } ;

  // obsolete visitor, need to be removed later
  class graphEditVisitor: public visitor {
  public:
    graphEditVisitor(const std::map<variable,int>& ot,
                     const std::map<variable,variable>& rt,
                     const variableSet& kv) ;
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  protected:
    std::map<variable,int> ownership_table ;
    std::map<variable,variable> recurrence_vars_table ;//area{n} -> area
                                                       //area -> area{n}
    variableSet keep_vars ;
  private:
    variableSet ot_vars ; // keys in the ownership table
    variableSet rvt_vars ;// keys in the recurrence_vars_table ;
    void edit_gr(digraph& gr,rulecomp_map& rcm,int id) ;
  } ;

  class graphVisualizeVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  } ;

  // generate allocation information table
  // used in top - down order
  class allocInfoVisitor: public visitor {
  public:
    // need graph_sn information
    allocInfoVisitor(const std::set<int>& gsn,
                     const variableSet& rtv,
                     const std::set<int>& lsn,
                     const std::map<int,variableSet>& rot_vt,
                     const std::map<int,variableSet>& lcommont,
                     const variableSet& untyped_vars)
      :graph_sn(gsn),recur_target_vars(rtv),
       loop_sn(lsn),rotate_vtable(rot_vt),
       loop_common_table(lcommont),
       allocated_vars(untyped_vars){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<int,variableSet> get_alloc_table()
    {return alloc_table ;}
    std::map<int,variableSet> get_loop_alloc_table()
    {return loop_alloc_table ;}
  protected:
    // return the all the variables that allocated
    // in the graph
    variableSet gather_info(const digraph& gr,int id) ;
    std::map<int,variableSet> alloc_table ;
    // holds info about the allocation of every loop super node
    std::map<int,variableSet> loop_alloc_table ;
    // variables that have been allocated up to now
    // or reserved not to be allocated
    variableSet allocated_vars ; 
    // super nodes that have a graph inside it
    std::set<int> graph_sn ;
    // the set of all the recurrence target variable
    variableSet recur_target_vars ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // used to decorate the graph to include allocation rules
  class allocGraphVisitor: public visitor {
  public:
    allocGraphVisitor(const std::map<int,variableSet>& t,
                      const std::set<int>& lsn,
                      const std::map<int,variableSet>& rot_vt,
                      const std::map<int,variableSet>& lcommont)
      :alloc_table(t),loop_sn(lsn),rotate_vtable(rot_vt),
       loop_common_table(lcommont){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  protected:
    void edit_gr(digraph& gr,rulecomp_map& rcm,int id) ;
    std::map<int,variableSet> alloc_table ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // generate delete information table
  // used in top - down order
  class deleteInfoVisitor: public visitor {
  public:
    // need loop allocation info and recurrence variable info
    deleteInfoVisitor(const std::map<int,variableSet>& lat,
                      const std::map<variable,variableSet>& rvt2s,
                      const std::map<variable,variableSet>& rvs2t,
                      const variableSet& rsv,
                      const variableSet& rtv,
                      const std::set<int>& gsn,
                      const std::map<int,int>& pnt,
                      const std::map<int,int>& lct,
                      const std::set<int>& lsn,
                      const std::map<int,variableSet>& rot_vt,
                      const std::map<int,variableSet>& lcommont,
                      const variableSet& target,
                      const variableSet& untyped_vars)
      :loop_alloc_table(lat),recur_vars_t2s(rvt2s),
       recur_vars_s2t(rvs2t),recur_source_vars(rsv),
       recur_target_vars(rtv),graph_sn(gsn),
       pnode_table(pnt),loop_ctable(lct),loop_sn(lsn),
       rotate_vtable(rot_vt),loop_common_table(lcommont),
       queried_target(target),deleted_vars(untyped_vars){
      
      for(std::map<int,int>::const_iterator mi=loop_ctable.begin();
          mi!=loop_ctable.end();++mi)
        col_sn.insert(mi->second) ;
      for(std::map<int,variableSet>::const_iterator mi=rotate_vtable.begin();
          mi!=rotate_vtable.end();++mi)
        all_rot_vars += mi->second ;
    }
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<int,variableSet> get_delete_table()
    {return delete_table ;}
    std::map<variable,ruleSet> get_recur_source_other_rules()
    {return recur_source_other_rules ;}
    std::map<variable,int> get_redirect_table()
    {return redirect_table ;}
  protected:
    // determine whether variables in working_vars
    // should be deleted in this graph

    // return all the variables that deleted
    // in the graph
    variableSet gather_info(const digraph& gr,
                            const variableSet& working_vars) ;

    // gather information for gather_info function use
    // get working_vars in the graph
    variableSet get_start_info(const digraph& gr,int id) ;
    // post processing of working variables
    void post_proc_working_vars(variableSet& working_vars) ;
    // constrain deletions of graphs that are inside n number of loops
    int only_loop_alloc(variableSet& working_vars,int id,int loop_num) ;
    // the looping algorithm for schedule deletion of variables
    // in the multilevel graph
    void looping_algr(const variableSet& working_vars,const digraph& dg,
                      int id,int start_loop_num) ;
    // function that checks the deletion of recurrence target vars
    void check_recur_target(const digraph& gr,int id) ;
    // given a rule set and a graph
    // return true if there is only one super node and
    // all other rules have path to the super node, or if
    // it is only one super node
    // return false otherwise
    bool let_it_go(const digraph& gr, ruleSet rules, const variable& v) ;
    // let_it_go for check_recur_target use (a different version)
    bool let_it_go(const digraph& gr, ruleSet rules) ;
    // the delete information table
    std::map<int,variableSet> delete_table ;
    // variables that have been deleted up to now or reserved
    // not to be deleted
    variableSet deleted_vars ;
    // info about recurrence variables (target to source map)
    std::map<variable,variableSet> recur_vars_t2s ;
    // source to target map
    std::map<variable,variableSet> recur_vars_s2t ;
    // all recurrence source variables
    variableSet recur_source_vars ;
    // all recurrence target variables
    variableSet recur_target_vars ;
    // if a variable is a recurrence target variable,
    // and there are other rules (other than the recurrence rule)
    // that its recurrence source variables reach in the same graph,
    // we then record such rules in this table
    std::map<variable,ruleSet> recur_source_other_rules ;
    // info about every loop's allocation
    std::map<int,variableSet> loop_alloc_table ;
    // super nodes that have a graph inside it
    std::set<int> graph_sn ;
    // table that records the parent node 
    std::map<int,int> pnode_table ;
    // table that holds the collapse node of each loop
    std::map<int,int> loop_ctable ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
    // these variables are for additional checking of the
    // possibility of deletion of recurrence target variables
    variableSet processed_targets ;
    std::map<variable,int> redirect_table ;
    // the queried variable set by the user
    variableSet queried_target ;
    // all the collapse node id
    std::set<int> col_sn ;
    // all the loop rotate variables
    variableSet all_rot_vars ;
  } ;

  // visitor that get all the recurrence variables in the
  // multilevel graph
  class recurInfoVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<variable,variableSet> get_recur_vars_t2s()
    {return recur_vars_t2s ;}
    std::map<variable,variableSet> get_recur_vars_s2t()
    {return recur_vars_s2t ;}
    variableSet get_recur_source_vars()
    {return recur_source_vars ;}
    variableSet get_recur_target_vars()
    {return recur_target_vars ;}
  protected:
    void gather_info(const digraph& gr) ;
    // from x{n} -> x, i.e. from target -> source
    std::map<variable,variableSet> recur_vars_t2s ; 
    // from source -> target, e.g. x -> x{n}
    std::map<variable,variableSet> recur_vars_s2t ;
    // set of recurrence source and target variables
    variableSet recur_source_vars ;
    variableSet recur_target_vars ;
  } ;

  // used to decorate the graph to include deletion rules
  class deleteGraphVisitor: public visitor {
  public:
    deleteGraphVisitor(const std::map<int,variableSet>& t,
                       const std::map<variable,ruleSet>& rsor)
      :delete_table(t),recur_source_other_rules(rsor){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  protected:
    void edit_gr(digraph& gr,rulecomp_map& rcm,int id) ;
    std::map<int,variableSet> delete_table ;
    // if a variable is a recurrence target variable,
    // and there are other rules (other than the recurrence rule)
    // that its recurrence source variables reach in the same graph,
    // this table holds such rules
    std::map<variable,ruleSet> recur_source_other_rules ;
  } ;

  // visitor to get some inter- super node information
  class snInfoVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::set<int> get_graph_sn()
    {return graph_sn ;}
    std::set<int> get_loop_sn()
    {return loop_sn ;}
    std::map<int,std::set<int> > get_subnode_table()
    {return subnode_table ;}
    std::map<int,int> get_loop_col_table()
    {return loop_col_table ;}
  protected:
    void fill_subnode_table(const digraph& gr, int id) ;
    // super node that has a graph in it
    // i.e. the loop, dag, and conditional super node
    std::set<int> graph_sn ;
    // all the loop super node id
    std::set<int> loop_sn ;
    // table that holds all the super nodes' id inside a super node
    std::map<int,std::set<int> > subnode_table ;
    // table that holds the collapse node of each loop
    std::map<int,int> loop_col_table ;
  } ;

  // check the if the allocation and deletion table are
  // consistent
  class allocDeleteStat {
    variableSet allocated_vars ;
    variableSet deleted_vars ;
    std::map<variable,variableSet> recur_vars_t2s ;
    std::map<variable,variableSet> recur_vars_s2t ;
    variableSet recur_source_vars ;
    variableSet recur_target_vars ;
  public:
    allocDeleteStat(const std::map<int,variableSet>& alloc_table,
                    const std::map<int,variableSet>& delete_table,
                    const std::map<variable,variableSet>& t2s,
                    const std::map<variable,variableSet>& s2t,
                    const variableSet& rsv,const variableSet& rtv) ;
    std::ostream& report(std::ostream& s) ;
  } ;

  // function to get the multilevel graph parent hierarchy table
  std::map<int,int>
  get_parentnode_table(const std::map<int,std::set<int> >& subnode_table) ;
  // function to get the allocation of each loop
  std::map<int,variableSet>
  get_loop_alloc_table(const std::map<int,variableSet>& alloc_table,
                       const std::map<int,std::set<int> >& subnode_table,
                       const std::set<int>& loop_sn,
                       const std::set<int>& graph_sn,
                       const std::map<variable,variableSet>& rvs2t) ;

  // visitor to compute the loop_rotate lists
  class rotateListVisitor: public visitor {
  public:
    rotateListVisitor(const sched_db& sd):scheds(sd) {}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) {}
    virtual void visit(conditional_compiler& cc) {}
    std::map<int,variableSet> get_rotate_vars_table()
    {return rotate_vars_table ;}
    std::map<int,variableSet> get_loop_common_table()
    {return loop_common_table ;}
  private:
    // reference to the schedule database
    const sched_db& scheds ;
    // table that holds variables in each loop's
    // rotate list
    std::map<int,variableSet> rotate_vars_table ;
    // table holds the common varibles between adv & col part of loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // visitor that checks if a graph has cycle
  class dagCheckVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    bool check_dag(digraph gr) ;
  } ;

  // visitor that discover all the un-typed variables
  // in the multilevel graph
  class unTypedVarVisitor: public visitor {
  public:
    unTypedVarVisitor(const std::map<variable,variableSet>& s2t,
                      const std::map<variable,variableSet>& t2s)
      :recur_vars_s2t(s2t),recur_vars_t2s(t2s) {}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    variableSet get_untyped_vars()
    {return untyped_vars ;}
  private:
    void discover(const digraph& gr) ;
    variableSet untyped_vars ;
    variableSet typed_vars ;
    // recurrence variable mapping table
    // source to target
    std::map<variable,variableSet> recur_vars_s2t ;
    // target to source
    std::map<variable,variableSet> recur_vars_t2s ;
  } ;

  // visitor that collects all the variables
  // in the multilevel graph
  class getAllVarVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    variableSet get_all_vars_ingraph()
    {return all_vars ;}
  private:
    void collect_vars(const digraph& gr) ;
    variableSet all_vars ;
  } ;

  // overload "<<" to print out an std::set
  template<typename T>
  inline std::ostream& operator<<(std::ostream& s, const std::set<T>& ss) {
    typename std::set<T>::const_iterator si ;
    si=ss.begin() ;
    s << *si ;
    ++si ;
    for(;si!=ss.end();++si) {
      s << " " ;
      s << *si ;
    }

    return s ;
  }
  
  // overload "<<" to print out an std::map
  template<typename T1, typename T2>
  inline std::ostream& operator<<(std::ostream& s, const std::map<T1,T2>& m) {
    typename std::map<T1,T2>::const_iterator mi ;
    for(mi=m.begin();mi!=m.end();++mi) {
      s << mi->first
        << ": " << mi->second
        << std::endl ;
    }

    return s ;
  }

}

#endif

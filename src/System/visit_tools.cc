#include "visitor.h"
#include "visit_tools.h"
#include "comp_tools.h"
#include <Tools/stream.h>

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include<map>
using std::map ;
#include <list>
using std::list ;
#include <algorithm>
using std::sort ;
using std::find_if ;
using std::copy ;
#include <functional>
using std::bind2nd ;
using std::ptr_fun ;
#include <iterator>
using std::ostream_iterator ;
#include <utility>
using std::pair ;
using std::make_pair ;
using std::min ;

namespace Loci {

  //////////////////////////////////////////////////////////////
  // function to get all the recur targets of a variableSet
  // from a given recur mapping table
  //////////////////////////////////////////////////////////////
  variableSet
  get_recur_target_for_vars(const variableSet& vars,
                            const map<variable,variableSet>& t) {
    variableSet ret ;
    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      ret += get_all_recur_vars(t,*vi) ;
    }

    return ret ;
  }


  //////////////////////////////////////////////////////////////
  // recurInfoVisitor
  //////////////////////////////////////////////////////////////
  void recurInfoVisitor::gather_info(const digraph& gr) {
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // looping over all the rules and gather info
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      // recurrence variables are the target and source of
      // these internal rules
      if(ruleIter->type() == rule::INTERNAL) {
        if(ruleIter->get_info().qualifier() == "generalize") {
          variable target,source ;
          target = *(ruleIter->targets().begin()) ;
          source = *(ruleIter->sources().begin()) ;

          recur_vars_t2s[target] += source ;
          recur_vars_s2t[source] += target ;
          recur_source_vars += source ;
          recur_target_vars += target ;

          generalize_t2s[target] += source ;
          generalize_s2t[source] += target ;
          generalize_source_vars += source ;
          generalize_target_vars += target ;
        }
        else if(ruleIter->get_info().qualifier() == "promote") {
          variable target,source ;
          target = *(ruleIter->targets().begin()) ;
          source = *(ruleIter->sources().begin()) ;

          recur_vars_t2s[target] += source ;
          recur_vars_s2t[source] += target ;
          recur_source_vars += source ;
          recur_target_vars += target ;

          promote_t2s[target] += source ;
          promote_s2t[source] += target ;
          promote_source_vars += source ;
          promote_target_vars += target ;
        }
        else if(ruleIter->get_info().qualifier() == "priority") {
          variable target,source ;
          target = *(ruleIter->targets().begin()) ;
          source = *(ruleIter->sources().begin()) ;

          recur_vars_t2s[target] += source ;
          recur_vars_s2t[source] += target ;
          recur_source_vars += source ;
          recur_target_vars += target ;

          priority_t2s[target] += source ;
          priority_s2t[source] += target ;
          priority_source_vars += source ;
          priority_target_vars += target ;
        }
      } else {
        // we check for rename (inplace update rules)
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ruleIter->get_info().desc.targets.begin();
            vmsi!=ruleIter->get_info().desc.targets.end(); ++vmsi) {
          if(vmsi->assign.size() != 0)
            for(unsigned int i=0;i<vmsi->assign.size();++i) {
              variable new_name = vmsi->assign[i].first ;
              variable old_name = vmsi->assign[i].second ;

              recur_vars_t2s[new_name] += old_name ;
              recur_vars_s2t[old_name] += new_name ;
              recur_source_vars += old_name ;
              recur_target_vars += new_name ;

              rename_t2s[new_name] += old_name ;
              rename_s2t[old_name] += new_name ;
              rename_source_vars += old_name ;
              rename_target_vars += new_name ;
            }
        }
      }
    }
    
  }

  void recurInfoVisitor::visit(loop_compiler& lc) {
    gather_info(lc.collapse_gr) ;
    gather_info(lc.advance_gr) ;
  }

  void recurInfoVisitor::visit(dag_compiler& dc) {
    gather_info(dc.dag_gr) ;
  }

  void recurInfoVisitor::visit(conditional_compiler& cc) {
    gather_info(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////
  // snInfoVisitor
  /////////////////////////////////////////////////////////////
  void snInfoVisitor::fill_subnode_table(const digraph& gr, int id) {
    // first get all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    set<int> subnode ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter))
        subnode.insert(get_supernode_num(*ruleIter)) ;
    }
    subnode_table[id] = subnode ;
  }

  void snInfoVisitor::visit(loop_compiler& lc) {
    graph_sn.insert(lc.cid) ;
    loop_sn.insert(lc.cid) ;

    fill_subnode_table(lc.loop_gr,lc.cid) ;

    rule collapse_node ;
    ruleSet all_rules = extract_rules(lc.collapse_gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri) {
      if(ri->target_time().before(ri->source_time())) {
        collapse_node = *ri ;
        break ;
      }
    }
    if(!is_super_node(collapse_node)) {
      cerr << "Internal error, the collapse part of loop compiler: "
           << lc.cid << " does not have a collapse rule" << endl ;
      exit(-1) ;
    }
    int collapse_id = get_supernode_num(collapse_node) ;
    if(collapse_id == -1) {
      cerr << "Error: conditional node has wrong id number" << endl ;
      exit(-1) ;
    }

    loop_col_table[lc.cid] = collapse_id ;

  }
  
  void snInfoVisitor::visit(dag_compiler& dc) {
    graph_sn.insert(dc.cid) ;
    fill_subnode_table(dc.dag_gr,dc.cid) ;
  }
  
  void snInfoVisitor::visit(conditional_compiler& cc) {
    graph_sn.insert(cc.cid) ;
    fill_subnode_table(cc.cond_gr,cc.cid) ;
  }

  ///////////////////////////////////////////////////////////
  // allocDeleteStat
  ///////////////////////////////////////////////////////////
  allocDeleteStat::allocDeleteStat(const map<int,variableSet>& alloc_table,
                                   const map<int,variableSet>& delete_table,
                                   const map<variable,variableSet>& t2s,
                                   const map<variable,variableSet>& s2t,
                                   const variableSet& rsv,
                                   const variableSet& rtv)
    :recur_vars_t2s(t2s),recur_vars_s2t(s2t),
     recur_source_vars(rsv),recur_target_vars(rtv){
    map<int,variableSet>::const_iterator mi ;
    for(mi=alloc_table.begin();mi!=alloc_table.end();++mi)
      allocated_vars += mi->second ;
    for(mi=delete_table.begin();mi!=delete_table.end();++mi)
      deleted_vars += mi->second ;
  }

  ostream& allocDeleteStat::report(ostream& s) const {
    variableSet d_a, a_d ;
    variableSet::const_iterator vi,vii ;

    s << "***** Begin allocation and deletion info report *****" << endl ;
    variableSet miss_alloc ;
    for(vi=allocated_vars.begin();vi!=allocated_vars.end();++vi) {
      if(recur_target_vars.inSet(*vi))
        miss_alloc += *vi ;
    }
    if(miss_alloc != EMPTY) {
      s << endl ;
      s << "These variables should not been allocated, because they are "
        << "recurrence target variables: " << endl << miss_alloc << endl ;
    }


    d_a = deleted_vars - allocated_vars ;

    // some recurrence variable processing
    // NEEDS SOME COMMENTS LATER HERE!!!!!
    variableSet remove ;
    for(vi=d_a.begin();vi!=d_a.end();++vi) {
      variableSet init_sources = get_leaf_recur_vars(recur_vars_t2s,*vi);
      
      if(init_sources != EMPTY) {
        bool all_allocated = true ;
        for(vii=init_sources.begin();vii!=init_sources.end();++vii)
          if(!allocated_vars.inSet(*vii)) {
            all_allocated = false ;
            break ;
          }

        if(all_allocated)
          remove += *vi ;
      }
    }
    d_a -= remove ;

    remove = variableSet(EMPTY) ;
    a_d = allocated_vars - deleted_vars ;
    for(vi=a_d.begin();vi!=a_d.end();++vi) {
      variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
      
      if(all_targets != EMPTY) {
        for(vii=all_targets.begin();vii!=all_targets.end();++vii)
          if(deleted_vars.inSet(*vii)) {
            remove += *vi ;
            break ;
          }
      }
    }
    a_d -= remove ;

    if(a_d == EMPTY) {
      s << endl ;
      s << "All allocated variables are scheduled to be deleted" << endl ;
    }else {
      s << endl ;
      s << "These variables are scheduled to be allocated, but not to be deleted:" << endl
        << a_d << endl ;
    }

    if(d_a != EMPTY) {
      s << endl ;
      s << "These variables are scheduled to be deleted, but they are input variables: " << endl << d_a << endl ;
    }
    s << endl ;
    s << "***** Finish allocation and deletion info report *****" << endl ;

    return s ;
  }

  /////////////////////////////////////////////////////////////////////
  // function to get the parent table
  /////////////////////////////////////////////////////////////////////  
  map<int,int>
  get_parentnode_table(const map<int,set<int> >& subnode_table) {
    map<int,int> ret ;

    // first get the parent node table
    map<int,set<int> >::const_iterator mi ;
    for(mi=subnode_table.begin();mi!=subnode_table.end();++mi) {
      // get all the subnodes of a node
      set<int> subnode = mi->second ;
      for(set<int>::const_iterator si=subnode.begin();
          si!=subnode.end();++si) {
        map<int,int>::const_iterator found ;
        found = ret.find(*si) ;
        if(found != ret.end()) {
          cerr << "multilevel graph error!" << endl ;
          exit(-1) ;
        }
        ret[*si] = mi->first ;
      }
    }
    // super node 0 has no parent node
    ret[0] = -1 ;

    return ret ;
  }

  ///////////////////////////////////////////////////////////////////////
  // function to get the allocation of each loop
  ///////////////////////////////////////////////////////////////////////
  map<int,variableSet>
  get_loop_alloc_table(const map<int,variableSet>& alloc_table,
                       const map<int,set<int> >& subnode_table,
                       const set<int>& loop_sn,
                       const set<int>& graph_sn,
                       const map<variable,variableSet>& rvs2t) {
    map<int,variableSet> ret ;

    for(set<int>::const_iterator si=loop_sn.begin();
        si!=loop_sn.end();++si) {
      // get all the subnode of a loop
      set<int> all_loop_subnode ;
      set<int> working ;

      all_loop_subnode.insert(*si) ;
      working.insert(*si) ;
      
      while(!working.empty()) {
        set<int> tmp ;
        for(set<int>::const_iterator si2=working.begin();
            si2!=working.end();++si2) {
          // skip recursive super node
          if(!inSet(graph_sn,*si2))
            continue ;
          map<int,set<int> >::const_iterator found ;
          found = subnode_table.find(*si2) ;
          FATAL(found == subnode_table.end()) ;
          tmp.insert( (found->second).begin(), (found->second).end()) ;
        }
        working = tmp ;
        all_loop_subnode.insert(tmp.begin(),tmp.end()) ;
      }

      // compute all the allocations of this loop
      variableSet loop_alloc ;
      for(set<int>::const_iterator si3=all_loop_subnode.begin();
          si3!=all_loop_subnode.end();++si3) {
        // skip recursive super node
        if(!inSet(graph_sn,*si3))
          continue ;
        map<int,variableSet>::const_iterator found ;
        if(inSet(loop_sn,*si3)) {
          // collapse part
          found = alloc_table.find(-(*si3)) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
          // advance part
          found = alloc_table.find(*si3) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
        } else {
          found = alloc_table.find(*si3) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
        }
      }
      // we need to include any recurrence variables
      variableSet loop_alloc_recur ;
      for(variableSet::const_iterator vi=loop_alloc.begin();
          vi!=loop_alloc.end();++vi) {
        // get all the recurrence target variables
        variableSet targets = get_all_recur_vars(rvs2t,*vi) ;
        if(targets != EMPTY)
          loop_alloc_recur += targets ;
      }
      loop_alloc += loop_alloc_recur ;
      // fill the ret table
      ret[*si] = loop_alloc ;
    }

    return ret ;
  }

  /////////////////////////////////////////////////////////////////////
  // rotateListVisitor
  /////////////////////////////////////////////////////////////////////
  namespace {
    inline bool offset_sort(const variable& v1, const variable& v2)
    { return v1.get_info().offset > v2.get_info().offset ; }
  }
  void rotateListVisitor::visit(loop_compiler& lc) {
    // compute the rotate list first
    
    map<variable,list<variable> > vlist ;
    for(variableSet::const_iterator vi=lc.all_loop_vars.begin();
        vi!=lc.all_loop_vars.end();++vi) {
      if(vi->time() == lc.tlevel && !vi->assign) {
        variable var_stationary(*vi,time_ident()) ;
        
        list<variable> s ;
        s.push_back(*vi) ;
        vlist[var_stationary].merge(s,offset_sort) ;
      }
    }
    
    map<variable,list<variable> >::const_iterator ii ;
    for(ii=vlist.begin();ii!=vlist.end();++ii) {
      if(ii->second.size() < 2) 
        continue ;
      std::list<variable>::const_iterator jj ;
      bool overlap = false ;
      variableSet vtouch ;
      for(jj=ii->second.begin();jj!=ii->second.end();++jj) {
        variableSet aliases = scheds.get_aliases(*jj) ;
        variableSet as = aliases ;
        variableSet::const_iterator vi ;
        
        for(vi=aliases.begin();vi!=aliases.end();++vi) 
          as += scheds.get_synonyms(*vi) ;
        if((as & vtouch) != EMPTY)
          overlap = true ;
        vtouch += as ;
      }
      if(!overlap) {
        lc.rotate_lists.push_back(ii->second) ;
      } else {
        if(ii->second.size() !=2) {
          cerr << "unable to have history on variables aliased in time"
               << endl
               << "error occured on variable " << ii->first
               << "{" << lc.tlevel << "}"
               << endl ;
          exit(-1) ;
        }
      }
    }

    // then get all the variables in the list and fill the table
    variableSet rotate_vars ;
    for(list<list<variable> >::const_iterator li=lc.rotate_lists.begin();
        li!=lc.rotate_lists.end();++li) {
      for(list<variable>::const_iterator lii=li->begin();
          lii!=li->end();++lii) {
        rotate_vars += *lii ;
      }
    }

    rotate_vars_table[lc.cid] = rotate_vars ;

    // and then compute the shared variables between
    // the advance part and the collapse part
    
    // we then get the common variable shared between the
    // advance and collapse part of the loop
    digraph::vertexSet collapse_gr_ver = lc.collapse_gr.get_all_vertices() ;
    digraph::vertexSet advance_gr_ver = lc.advance_gr.get_all_vertices() ;
    variableSet common = extract_vars(collapse_gr_ver & advance_gr_ver) ;
    // exclude time level variable
    variable tvar = variable(lc.tlevel) ;
    common -= tvar ;

    loop_common_table[lc.cid] = common ;
      
  }

  ////////////////////////////////////////////////////////////////
  // dagCheckVisitor
  ////////////////////////////////////////////////////////////////
  // return true if gr has NO cycle
  // return false if gr has cycle
  bool dagCheckVisitor::check_dag(digraph gr) {
    // if graph is empty, it is dag
    if(is_dg_empty(gr))
      return true ;
    // get all the leaf nodes of this graph
    digraph::vertexSet leaf_nodes ;
    do {
      digraph::vertexSet allv ;

      leaf_nodes = EMPTY ;
      allv = gr.get_all_vertices() ;
      for(digraph::vertexSet::const_iterator vi=allv.begin();
          vi!=allv.end();++vi) {
        if(gr[*vi] == EMPTY)
          leaf_nodes += *vi ;
      }
      gr.remove_vertices(leaf_nodes) ;
      
    }while(leaf_nodes != EMPTY) ;

    return is_dg_empty(gr) ;
  }

  void dagCheckVisitor::visit(loop_compiler& lc) {
    if(!check_dag(lc.collapse_gr)) {
      cerr << "ERROR: the collapse graph of loop super node("
           << lc.cid << ") has cycle(s)" << endl ;
      exit(-1) ;
    }
    if(!check_dag(lc.advance_gr)) {
      cerr << "ERROR: the advance graph of loop super node("
           << lc.cid << ") has cycle(s)" << endl ;
      exit(-1) ;
    }
  }
  
  void dagCheckVisitor::visit(dag_compiler& dc) {
    if(!check_dag(dc.dag_gr)) {
      cerr << "ERROR: the graph of dag super node("
           << dc.cid << ") has cycle(s)" << endl ;
      exit(-1) ;
    }
  }
  
  void dagCheckVisitor::visit(conditional_compiler& cc) {
    if(!check_dag(cc.cond_gr)) {
      cerr << "ERROR: the graph of conditional super node("
           << cc.cid << ") has cycle(s)" << endl ;
      exit(-1) ;
    }
  }

  ////////////////////////////////////////////////////////////////
  // unTypedVarVisitor
  ////////////////////////////////////////////////////////////////
  unTypedVarVisitor::
  unTypedVarVisitor(const std::map<variable,variableSet>& s2t,
                    const std::map<variable,variableSet>& t2s,
                    const variableSet& input)
    :recur_vars_s2t(s2t),recur_vars_t2s(t2s) {

    typed_vars += input ;
    variableSet add ;
    for(variableSet::const_iterator vi=typed_vars.begin();
        vi!=typed_vars.end();++vi) {
      variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
      add += all_targets ;
    }
    typed_vars += add ;
  }
  
  void unTypedVarVisitor::discover(const digraph& gr) {
    digraph grt = gr.transpose() ;
    // typed variables in this graph
    variableSet typed_here ;
    // get all target variables
    variableSet allvars ;
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      allvars += ri->targets() ;

    // typed vars don't need to be processed
    allvars -= typed_vars ;
    for(variableSet::const_iterator vi=allvars.begin();
        vi!=allvars.end();++vi) {
      // process one variable each time
      ruleSet source_rules ;
      
      ruleSet working ;
      working = extract_rules(grt[vi->ident()]) ;
      source_rules += working ;
      while(working != EMPTY) {
        ruleSet next ;
        for(ruleSet::const_iterator ri=working.begin();
            ri!=working.end();++ri) {
          variableSet source_vars ;
          if(!is_super_node(ri)) {
            source_vars = extract_vars(grt[ri->ident()]) ;
            for(variableSet::const_iterator vi2=source_vars.begin();
                vi2!=source_vars.end();++vi2)
              next += extract_rules(grt[vi2->ident()]) ;
          }
        }
        source_rules += next ;
        working = next ;
      } // end of while

      bool is_untyped = true ;
      for(ruleSet::const_iterator ri=source_rules.begin();
          ri!=source_rules.end();++ri) {
        if(!is_internal_rule(ri)) {
          is_untyped = false ;
          break ;
        }
      }
      if(is_untyped)
        untyped_vars += *vi ;
      else
        typed_here += *vi ;
    } // end of for
    // we make all typed variables' recurrence target variable(if any)
    // and all the recurrence source variable(if any) typed
    typed_vars += typed_here ;
    for(variableSet::const_iterator vi=typed_here.begin();
        vi!=typed_here.end();++vi) {
      variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
      typed_vars += all_targets ;
    }
    for(variableSet::const_iterator vi=typed_here.begin();
        vi!=typed_here.end();++vi) {
      variableSet all_sources = get_all_recur_vars(recur_vars_t2s,*vi) ;
      typed_vars += all_sources ;
    }
    for(variableSet::const_iterator vi=typed_vars.begin();
        vi!=typed_vars.end();++vi)
      if(untyped_vars.inSet(*vi))
        untyped_vars -= *vi ;
  }

  void unTypedVarVisitor::visit(loop_compiler& lc) {
    discover(lc.collapse_gr) ;
    discover(lc.advance_gr) ;
  }
  
  void unTypedVarVisitor::visit(dag_compiler& dc) {
    discover(dc.dag_gr) ;
  }
  
  void unTypedVarVisitor::visit(conditional_compiler& cc) {
    discover(cc.cond_gr) ;
  }

  ////////////////////////////////////////////////////////////////
  // getAllVarVisitor
  ////////////////////////////////////////////////////////////////

  void getAllVarVisitor::collect_vars(const digraph& gr) {
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    all_vars += extract_vars(allvertices) ;
  }

  void getAllVarVisitor::visit(loop_compiler& lc) {
    collect_vars(lc.loop_gr) ;
  }
  
  void getAllVarVisitor::visit(dag_compiler& dc) {
    collect_vars(dc.dag_gr) ;
  }
  
  void getAllVarVisitor::visit(conditional_compiler& cc) {
    collect_vars(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // promotePPVisitor
  /////////////////////////////////////////////////////////////////
  promotePPVisitor::
  promotePPVisitor(const map<variable,variableSet>& pt2s,
                   const map<variable,variableSet>& ps2t,
                   const variableSet& psource,
                   const variableSet& ptarget,
                   const set<int>& gsn,
                   variableSet& input):
    promote_source_vars(psource), promote_target_vars(ptarget),
    promote_t2s(pt2s), promote_s2t(ps2t), graph_sn(gsn) {
    
    reserved_vars += input ;
    reserved_vars += get_recur_target_for_vars(input,promote_s2t) ;
  }

  // given a rule set a variable can reach, if the variable
  // must be deleted in this graph, then it is a rep
  bool promotePPVisitor::is_rep(const digraph& gr, ruleSet rules) {
    if(rules == EMPTY) {
      return true ;
    }
    
    int supernode_num = 0 ;
    int supernode_id = 0 ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        // we eliminate the recursive super node
        // by looking up the graph_sn set
        if(inSet(graph_sn,get_supernode_num(*ruleIter))) {
          supernode_id = ruleIter->ident() ;
          ++supernode_num ;
        }
      }
    }

    if( (rules.size() == 1) && is_internal_rule(rules.begin()))
      return false ;
    if(supernode_num != 1)
      return true ;
    // the remaining case must be there are only one super node
    // and at least one non super node, we check if there is path
    // from each non super node to the super node in the graph
    rules -= rule(supernode_id) ;
    bool ret = false ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(!has_path(gr,ruleIter->ident(),supernode_id)) {
        ret = true ;
        break ;
      }
    }
    
    return ret ;
  }

  void promotePPVisitor::pick_rep(const digraph& gr) {
    // get variables to evaluate
    digraph::vertexSet allv = gr.get_all_vertices() ;
    
    ruleSet rules = extract_rules(allv) ;
    variableSet vars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      vars += ri->sources() ;

    vars &= variableSet(promote_source_vars + promote_target_vars) ;
    vars -= reserved_vars ;

    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      if(processed.inSet(*vi))
        continue ;
      // get the final target
      variableSet targets = get_leaf_recur_vars(promote_s2t,*vi) ;
      
      if(targets.size() == 1) {
        if(allv.inSet(targets.begin()->ident())) {
          processed += *vi ;
          continue ;
        }
      }

      // get the rules all_source can reach
      ruleSet others ;
      ruleSet tmp = extract_rules(gr[vi->ident()]) ;
      others += tmp ;
      for(ruleSet::const_iterator ruleIter=tmp.begin();
          ruleIter!=tmp.end();++ruleIter)
        if(!is_super_node(ruleIter))
          others += extract_rules(gr[ruleIter->ident()]) ;
      
      FATAL(others == EMPTY) ;

      if(is_rep(gr,others)) {
        rep += *vi ;
        variableSet all_source = get_all_recur_vars(promote_t2s,*vi) ;
        variableSet all_target = get_all_recur_vars(promote_s2t,*vi) ;
        processed += *vi ;
        processed += all_source ;
        processed += all_target ;
        remaining += all_source ;
        remaining += all_target ;
      }
    }
  }
  
  void promotePPVisitor::visit(loop_compiler& lc) {
    pick_rep(lc.loop_gr) ;
  }

  void promotePPVisitor::visit(dag_compiler& dc) {
    pick_rep(dc.dag_gr) ;
  }

  void promotePPVisitor::visit(conditional_compiler& cc) {
    pick_rep(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // function to analyze the renamed target variable cluster
  // and pick out a representitive variable for each cluster,
  // only this variable is allowed to be deleted
  /////////////////////////////////////////////////////////////////
  variableSet pick_rename_target(const map<variable,variableSet>& s2t,
                                 const map<variable,variableSet>& t2s,
                                 const variableSet& allvars) {
    variableSet filtered_vars ;
    
    variableSet processed ;
    for(variableSet::const_iterator vi=allvars.begin();
        vi!=allvars.end();++vi) {
      if(processed.inSet(*vi))
        continue ;

      variableSet all_other_target = get_all_leaf_target(t2s,s2t,*vi) ;

      if(all_other_target == EMPTY)
        continue ;
      
      processed += all_other_target ;
      
      vector<variable> all_targets ;
      all_targets.push_back(*vi) ;
      for(variableSet::const_iterator vi2=all_other_target.begin();
          vi2!=all_other_target.end();++vi2)
        all_targets.push_back(*vi2) ;

      sort(all_targets.begin(),all_targets.end(),time_before) ;

      /*
      cerr << "found target cluster: " ;
      copy(all_targets.begin(),all_targets.end(),
           ostream_iterator<variable>(cerr,",")) ;
      cerr << endl ;
      */
      
      vector<variable>::const_iterator pos ;
      pos = find_if(all_targets.begin()+1,all_targets.end(),
                    bind2nd(ptr_fun(time_after),all_targets[0])
                    ) ;

      if(pos - all_targets.begin() > 1) {
        variableSet time_ident_vars ;
        for(vector<variable>::const_iterator vi2=all_targets.begin();
            vi2!=pos;++vi2)
          time_ident_vars += *vi2 ;

        cerr << "WARNING: These renamed variables coexist in the same time level, and they refer to the same memory location, this is dangerous!: " << time_ident_vars << endl ;
      }

      for(vector<variable>::const_iterator vi2=all_targets.begin()+1;
          vi2!=all_targets.end();++vi2)
        filtered_vars += *vi2 ;
    }

    return filtered_vars ;
  }

  ///////////////////////////////////////////////////////////////////
  // function that checks some preconditions that all the recurrence
  // variables should meet.
  // generalize and priority rules should not form branches
  // promoted variables should not been renamed, and input variables
  // also should not been renamed
  ///////////////////////////////////////////////////////////////////
  void check_recur_precondition(const recurInfoVisitor& v,
                                const variableSet& input) {
    variableSet::const_iterator vi ;
    map<variable,variableSet>::const_iterator found ;

    variableSet generalize_source = v.get_generalize_source_vars() ;
    variableSet promote_source = v.get_promote_source_vars() ;
    variableSet priority_source = v.get_priority_source_vars() ;
    variableSet rename_source = v.get_rename_source_vars() ;
    
    map<variable,variableSet> generalize_s2t = v.get_generalize_s2t() ;
    map<variable,variableSet> promote_s2t = v.get_promote_s2t() ;
    map<variable,variableSet> priority_s2t = v.get_priority_s2t() ;
    map<variable,variableSet> rename_s2t = v.get_rename_s2t() ;
    
    // generalize rules should not form branches
    for(vi=generalize_source.begin();vi!=generalize_source.end();++vi) {
      found = generalize_s2t.find(*vi) ;
      FATAL(found == generalize_s2t.end()) ;
      if(found->second.size() > 1) {
        cerr << "WARNING: " << *vi << " is in the chain of "
             << "generalize rules, but it forms multiple targets: "
             << found->second << endl ;
      }
      if(promote_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of promote rules " ;
        found = promote_s2t.find(*vi) ;
        FATAL(found == promote_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
      if(priority_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of priority rules " ;
        found = priority_s2t.find(*vi) ;
        FATAL(found == priority_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
      if(rename_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of rename rules " ;
        found = rename_s2t.find(*vi) ;
        FATAL(found == rename_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
    }

    // priority rules should not form branches
    for(vi=priority_source.begin();vi!=priority_source.end();++vi) {
      found = priority_s2t.find(*vi) ;
      FATAL(found == priority_s2t.end()) ;
      if(found->second.size() > 1) {
        cerr << "WARNING: " << *vi << " is in the chain of "
             << "priority rules, but it forms multiple targets: "
             << found->second << endl ;
      }
      if(generalize_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of generalize rules " ;
        found = generalize_s2t.find(*vi) ;
        FATAL(found == generalize_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
      if(promote_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of promote rules " ;
        found = promote_s2t.find(*vi) ;
        FATAL(found == promote_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
      if(rename_source.inSet(*vi)) {
        cerr << "\tit is also in the chain of rename rules " ;
        found = rename_s2t.find(*vi) ;
        FATAL(found == rename_s2t.end()) ;
        cerr << "forms targets: " << found->second << endl ;
      }
    }

    // promoted variables should not be renamed
    variableSet promote_vars = v.get_promote_target_vars() ;
    for(vi=promote_vars.begin();vi!=promote_vars.end();++vi) {
      found = rename_s2t.find(*vi) ;
      if(found != rename_s2t.end()) {
        variableSet rt = found->second ;
        cerr << "WARNING: promoted variable: " << *vi
             << " is renamed to: " << rt << endl ;
      }
    }

    // input variables should not be renamed
    for(vi=input.begin();vi!=input.end();++vi) {
      found = rename_s2t.find(*vi) ;
      if(found != rename_s2t.end()) {
        variableSet rt = found->second ;
        cerr << "WARNING: input variable: " << *vi
             << " is renamed to: " << rt << endl ;
      }
      variableSet::const_iterator vi2 ;
      variableSet target ;
      
      target = get_all_recur_vars(promote_s2t,*vi) ;
      for(vi2=target.begin();vi2!=target.end();++vi2) {
        found = rename_s2t.find(*vi) ;
        if(found != rename_s2t.end()) {
          cerr << "WARNING: variable " << *vi2
               << " is promoted from input variable " << *vi
               << " but is renamed to: " << found->second
               << endl ;
        }
      }
    }

  }

  ///////////////////////////////////////////////////////////////
  // chopRuleVisitor
  ///////////////////////////////////////////////////////////////
  namespace {
    // given a contrete rule, return all the sources
    // that have maps
    inline variableSet map_sources(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet sources ;
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        if(!si->mapping.empty())
          sources += si->var ;
      }
      return sources ;
    }
    
    // given a contrete rule, return all the targetss
    // that have maps
    inline variableSet map_targets(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(!si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }
    
    // given a contrete rule, return all the targets
    // that do not have maps
    inline variableSet non_map_targets(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }
    
    inline variableSet non_map_targets(const ruleSet::const_iterator& ri) {
      const rule_impl::info& rinfo = ri->get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }

    // return the vertexSet of the rule's sources, targets
    // and the rule itself
    inline digraph::vertexSet get_rule_vertexSet(const rule& r) {
      digraph::vertexSet ret ;
      variableSet vars ;
      vars += r.sources() ;
      vars += r.targets() ;
      ret = get_vertexSet(vars) ;
      ret += r.ident() ;
      return ret ;
    }

    digraph::vertexSet get_ruleSet_vertexSet(const ruleSet& rules) {
      digraph::vertexSet ret ;
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
        ret += get_rule_vertexSet(*ri) ;
      }
      return ret ;
    }
    
    // determine if the two chop chains can be merged
    bool
    if_merge_chop_chain(const chop_chain& c1, const chop_chain& c2) {
      // get the two chop_vars from c1 and c2
      variableSet c1vars = c1.second ;
      variableSet c2vars = c2.second ;

      return ( (c1vars & c2vars) != EMPTY) ;
    }

  } // end of unamed namespace

  void chopRuleVisitor::find_chain(const digraph& gr, int id) {
    // get all vertices in the graph
    digraph::vertexSet allv = gr.get_all_vertices() ;
    // get the transpose of the graph
    digraph grt = gr.transpose() ;
    // find out the init source vertices
    digraph::vertexSet init ;
    for(digraph::vertexSet::const_iterator vi=allv.begin();
        vi!=allv.end();++vi)
      if(grt[*vi] == EMPTY)
        init += *vi ;
    // find out the init contrete rule set in the graph
    digraph::vertexSet tmpv ;
    for(digraph::vertexSet::const_iterator vi=init.begin();
        vi!=init.end();++vi)
      tmpv += gr[*vi] ;
    tmpv += init ;
    ruleSet init_crules ;
    ruleSet tmpr = extract_rules(tmpv) ;
    for(ruleSet::const_iterator ri=tmpr.begin();
        ri!=tmpr.end();++ri)
      if(!is_internal_rule(ri))
        init_crules += *ri ;

    //search the graph to find out the chain, if any
    list<chop_chain> chop_chain_list ;
    for(ruleSet::const_iterator ri=init_crules.begin();
        ri!=init_crules.end();++ri) {
      ruleSet allrules ;
      ruleSet working ;
      variableSet processed_vars ;
      
      working += *ri ;

      digraph::vertexSet vertices ;
      variableSet chop_vars ;
      while(working != EMPTY) {
        allrules += working ;
        vertices += get_ruleSet_vertexSet(working) ;
        ruleSet next ;
        for(ruleSet::const_iterator ri2=working.begin();
            ri2!=working.end();++ri2) {
          // get non map targets of each rule
          variableSet nmt = non_map_targets(ri2) ;
          nmt -= processed_vars ;
          processed_vars += nmt ;

          // for each variable, get the possible rules that can be reached
          for(variableSet::const_iterator vi=nmt.begin();
              vi!=nmt.end();++vi) {
            // check if it is a store
            storeRepP srp = facts.get_variable(*vi) ;
            if(srp->RepType() != Loci::STORE)
              continue ;

            digraph::vertexSet next_vertices ;
            next_vertices = gr[vi->ident()] ;
            if(next_vertices == EMPTY)
              continue ;

            ruleSet tmp = extract_rules(next_vertices) ;
            ruleSet tmp2 ;
            variableSet all_sources_have_map ;
            bool has_supernode = false ;
            for(ruleSet::const_iterator ri3=tmp.begin();
                ri3!=tmp.end();++ri3) {
              // a possible rule should be a contrete rule
              if(!is_internal_rule(ri3)) {
                // and the sources that have map should
                // not be in chop_vars
                
                // get the rule's source that involve maps
                all_sources_have_map += map_sources(*ri3) ;
                tmp2 += *ri3 ;
              } else {
                has_supernode = true ;
                break ;
              }
            }

            if(has_supernode)
              continue ;
            if(all_sources_have_map.inSet(*vi))
              continue ;

            // we also need to check if all the rules that
            // generate this variable have maps for
            // this variable
            digraph::vertexSet parent_vertices ;
            parent_vertices = grt[vi->ident()] ;
            tmp = extract_rules(parent_vertices) ;
            variableSet all_targets_have_map ;
            for(ruleSet::const_iterator ri3=tmp.begin();
                ri3!=tmp.end();++ri3)
              if(!is_internal_rule(ri3))
                all_targets_have_map += map_targets(*ri3) ;
            if(all_targets_have_map.inSet(*vi))
              continue ;

            if(tmp2 != EMPTY) {
              chop_vars += *vi ;
              next += tmp2 ;
            }
          }
        }
        working = next ;
      }
      if(allrules.size() > 1) {
        digraph sub_gr = gr.subgraph(vertices) ;
        chop_chain_list.push_back(make_pair(sub_gr,chop_vars)) ;
      }
    }
    
    // we got the chop_chain_list now, and we do some
    // post processing to merge the list, if it can be
    list<chop_chain>::iterator fix,move ;
    fix = chop_chain_list.begin() ;
    while(fix != chop_chain_list.end()) {
      move = fix ;
      ++move ;
      bool no_merge = true ;
      while(move != chop_chain_list.end()) {
        if(if_merge_chop_chain(*fix,*move)) {
          digraph::vertexSet new_gr_vertices =
            fix->first.get_all_vertices() +
            move->first.get_all_vertices() ;
          digraph new_gr = gr.subgraph(new_gr_vertices) ;
          
          variableSet new_chop_vars = fix->second ;
          new_chop_vars += move->second ;
          
          *fix = make_pair(new_gr,new_chop_vars) ;

          list<chop_chain>::iterator old_move = move ;
          ++move ;
          chop_chain_list.erase(old_move) ;

          no_merge = false ;
        }
        else
          ++move ;
      }
      if(no_merge)
        ++fix ;
    }

    if(!chop_chain_list.empty())
      all_chains[id] = chop_chain_list ;
  }
  
  void chopRuleVisitor::visit(loop_compiler& lc) {
    find_chain(lc.collapse_gr,-lc.cid) ;
    find_chain(lc.advance_gr,lc.cid) ;
  }

  void chopRuleVisitor::visit(dag_compiler& dc) {
    find_chain(dc.dag_gr,dc.cid) ;
  }

  void chopRuleVisitor::visit(conditional_compiler& cc) {
    find_chain(cc.cond_gr,cc.cid) ;
  }
  
} // end of namespace Loci


#include "visitor.h"
#include "visit_tools.h"
#include "comp_tools.h"
#include <Tools/stream.h>

using std::vector ;
using std::map ;
using std::set ;


namespace Loci {

  ////////////////////////////////////////////////////////////////
  // allocInfoVisitor
  ////////////////////////////////////////////////////////////////
  variableSet allocInfoVisitor::gather_info(const digraph& gr,int id) {
    // obtain the transpose of the graph
    digraph grt = gr.transpose() ;
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // variables processed in this graph
    variableSet working_vars ;
    
    // looping over all the rules and gather variables
    // need to be processed in this graph for allocation
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      // get the targets of the rule
      variableSet targets ;
      targets = ruleIter->targets() ;

      working_vars += targets ;
    }

    //////////////////////////////////////////////////
    // NOTE: THIS IS TEMPORARY
    // but we add the adjusts
    variableSet remove = variableSet(recur_target_vars - adjust_vars) ;
    working_vars -= remove ;
    // we remove recurrence target variables from the working_vars.
    //working_vars -= recur_target_vars ;

    // variables that allocated in this level (graph)
    variableSet alloc_here ;
    
    // looping over working_vars and gather allocation info
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // if this variable has been allocated before
      // we don't need to process it here
      if(allocated_vars.inSet(*varIter))
        continue ;      
      // get the rules that produce this variable
      ruleSet source_rules = extract_rules(grt[varIter->ident()]) ;
      // see if these rules are all super node
      bool all_super_nodes = true ;
      for(ruleSet::const_iterator ruleIter=source_rules.begin();
          ruleIter!=source_rules.end();++ruleIter) {
        if(!is_super_node(ruleIter)) {
          all_super_nodes = false ;
          break ;
        }
        else {
          // we eliminate the recursive super node
          // by look up the graph_sn set
          set<int>::const_iterator found ;
          found = graph_sn.find(get_supernode_num(*ruleIter)) ;
          if(found == graph_sn.end()) {
            all_super_nodes = false ;
            break ;
          }
        }
      }
      // if all the source rules are super node
      // there are two possible cases:
      
      // the first is there is only ONE source super node
      // then we don't allocate it in this level (graph)
      // the second is there are more than one super node
      // that produce this variable, then we have to allocate
      // this variable in this level (graph)
      if( (all_super_nodes) &&
          (source_rules.size() == 1)
          )
        continue ;
      // else we allocate it in this level
      // add the variable into the alloc_here set
      alloc_here += *varIter ;
    }
    // we check if there is any loop super node in the graph,
    // if yes, we also allocate the rotate list variables
    // and the shared variables for the loop in this graph.
    variableSet all_rot_vars ;
    variableSet all_shared_vars ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        int id = get_supernode_num(*ruleIter) ;
        if(inSet(loop_sn,id)) {
          // get the rotate list variables
          map<int,variableSet>::const_iterator found ;
          found = rotate_vtable.find(id) ;
          FATAL(found == rotate_vtable.end()) ;
          all_rot_vars += found->second ;

          found = loop_shared_table.find(id) ;
          FATAL(found == loop_shared_table.end()) ;
          all_shared_vars += found->second ;
        }
      }
    }
    all_rot_vars -= recur_target_vars ;
    all_rot_vars -= allocated_vars ;
    all_shared_vars -= recur_target_vars ;
    all_shared_vars -= allocated_vars ;
    
    alloc_here += all_rot_vars ;
    alloc_here += all_shared_vars ;

    
    // add the variables into the allocated_vars set
    allocated_vars += alloc_here ;
    // edit alloc_table
    alloc_table[id] = alloc_here ;
    // return the alloc_here
    return alloc_here ;
  }

  void allocInfoVisitor::visit(loop_compiler& lc) {
    variableSet tmp,loop_alloc_vars ;

    // suppose the id of loop_compiler is x,
    // then -x refers to the collapse part
    // x refers to the advance part
    tmp=gather_info(lc.collapse_gr,-lc.cid) ;
    loop_alloc_vars += tmp ;
    
    tmp=gather_info(lc.advance_gr,lc.cid) ;
    loop_alloc_vars += tmp ;

    loop_alloc_table[lc.cid]=loop_alloc_vars ;
  }

  void allocInfoVisitor::visit(dag_compiler& dc) {
    gather_info(dc.dag_gr,dc.cid) ;
  }

  void allocInfoVisitor::visit(conditional_compiler& cc) {
    gather_info(cc.cond_gr,cc.cid) ;
  }

  /////////////////////////////////////////////////////////
  // allocGraphVisitor
  /////////////////////////////////////////////////////////

  // we decorate the graph to include the allocation rules
  void allocGraphVisitor::edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // first get the transposed graph
    digraph grt = gr.transpose() ;
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // variables processed in this graph
    variableSet working_vars ;

    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      working_vars += ruleIter->targets() ;
    }

    // first we get the info in alloc_table, find out which variables
    // need to be allocated in the graph
    map<int,variableSet>::const_iterator found ;
    found = alloc_table.find(id) ;
    FATAL(found == alloc_table.end()) ;
    variableSet alloc_vars = found->second ;

    if(alloc_vars == EMPTY) return ;
    
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // we check if the variable is in the alloc_vars
      // if it is not in the alloc_vars, we don't allocate it here
      if(!alloc_vars.inSet(*varIter))
        continue ;
      // now we decorate the graph to include the allocation
      // we first obtain the ident of rules that produce the variable
      digraph::vertexSet source_rules = grt[varIter->ident()] ;
      // if source_rules is empty, then the variable is not
      // present in this graph, we need a search to find out
      // the source_rules
      if(source_rules == EMPTY) {
        for(ruleSet::const_iterator ruleIter=rules.begin();
            ruleIter!=rules.end();++ruleIter) {
          variableSet targets = ruleIter->targets() ;
          if(targets.inSet(*varIter))
            source_rules += ruleIter->ident() ;
        }
      }
      
      // we create a rule for allocation
      variable sv("CREATE") ;
      rule alloc_rule = create_rule(sv,*varIter,"ALLOCATE") ;
      // edit graph
      gr.add_edges(alloc_rule.ident(),source_rules) ;
      // modify the rulecomp_map
      rcm[alloc_rule] = new allocate_var_compiler(*varIter) ;
      // remove variable in alloc_vars
      alloc_vars -= *varIter ;
    }
    // if alloc_vars is not empty, it must contain loop rotate list
    // variables and loop shared variables. we allocate them here
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      variableSet rotate_vars ;
      variableSet shared_vars ;
      variableSet allocated ;
      
      if(is_super_node(ruleIter)) {
        int id = get_supernode_num(*ruleIter) ;
        if(inSet(loop_sn,id)) {
          // get the rotate list variables
          map<int,variableSet>::const_iterator found ;
          found = rotate_vtable.find(id) ;
          FATAL(found == rotate_vtable.end()) ;
          rotate_vars += found->second ;

          // get the shared variables
          found = loop_shared_table.find(id) ;
          FATAL(found == loop_shared_table.end()) ;
          shared_vars += found->second ;

          for(variableSet::const_iterator vi=alloc_vars.begin();
              vi!=alloc_vars.end();++vi) {
            // rotate list vars & shared vars
            if(rotate_vars.inSet(*vi) || shared_vars.inSet(*vi)) {
              allocated += *vi ;
              
              // we create a rule for allocation
              variable sv("CREATE") ;
              rule alloc_rule = create_rule(sv,*vi,"ALLOCATE") ;
              // edit graph
              gr.add_edge(alloc_rule.ident(),ruleIter->ident()) ;
              // modify the rulecomp_map
              rcm[alloc_rule] = new allocate_var_compiler(*vi) ;
            }
          }
          alloc_vars -= allocated ;
        }
      }
    }
    
    // if alloc_vars is not empty at this stage, it means
    // there still have some variables need to be allocated in this
    // graph but these variables do not belong to this graph

    // generally it is not possible to have such case at this stage, we
    // leave it for future development. We make a warning here
    if(alloc_vars != EMPTY) {
      cerr << "Warning: allocation not complete in graph: "
           << id << endl
           << "Remaining variables are: " << alloc_vars << endl ;
    }
  }

  void allocGraphVisitor::visit(loop_compiler& lc) {
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,-lc.cid) ;
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }

  void allocGraphVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }

  void allocGraphVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }

  //////////////////////////////////////////////////////////
  // deleteInfoVisitor
  //////////////////////////////////////////////////////////
  deleteInfoVisitor::
  deleteInfoVisitor(const map<int,variableSet>& lat,
                    const map<variable,variableSet>& rvt2s,
                    const map<variable,variableSet>& rvs2t,
                    const set<int>& gsn,
                    const map<int,int>& pnt,
                    const map<int,int>& lct,
                    const set<int>& lsn,
                    const map<int,variableSet>& rot_vt,
                    const map<int,variableSet>& lsharedt,
                    const variableSet& reserved_vars)
    :recur_vars_t2s(rvt2s), recur_vars_s2t(rvs2t),
     loop_alloc_table(lat), graph_sn(gsn),
     pnode_table(pnt), loop_ctable(lct),
     loop_sn(lsn), rotate_vtable(rot_vt),
     loop_shared_table(lsharedt){
      
    for(std::map<int,int>::const_iterator mi=loop_ctable.begin();
        mi!=loop_ctable.end();++mi)
      col_sn.insert(mi->second) ;
    for(std::map<int,variableSet>::const_iterator mi=rotate_vtable.begin();
        mi!=rotate_vtable.end();++mi)
      all_rot_vars += mi->second ;

    deleted_vars += reserved_vars ;
  }

  variableSet deleteInfoVisitor::get_start_info(const digraph& gr,int id) {

    variableSet working_vars ;

    // obtain all the possible variables in the graph
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      working_vars += ri->sources() ;
    
    working_vars += extract_vars(allvertices) ;

    // remove any deleted or reserved variables
    working_vars -= deleted_vars ;
    
    return working_vars ;    
  }

  variableSet deleteInfoVisitor::gather_info(const digraph& gr,
                                             const variableSet& working_vars) {
    if(working_vars == EMPTY)
      return variableSet(EMPTY) ;
    // get all vertices in the graph
    digraph::vertexSet allv = gr.get_all_vertices() ;
    // variables that been deleted in this graph
    variableSet delete_here ;
    // looping over all the working variables to gather info
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {

      // first some recurrence variable processing

      // all the recurrence source variables of *varIter
      variableSet all_recur_sources ;

      all_recur_sources = get_all_recur_vars(recur_vars_t2s,*varIter) ;

      // get the rules that this variable reaches
      ruleSet target_rules = extract_rules(gr[varIter->ident()]) ;
      // if this variable is in the recurrence variable table,
      // we'll also need to get the rules that the recurrence variable
      // reach
      if(all_recur_sources != EMPTY) {
        // get the recurrence source variables that are in this graph
        digraph::vertexSet in_graph_recur_vars ;
        for(variableSet::const_iterator vi=all_recur_sources.begin();
            vi!=all_recur_sources.end();++vi) {
          if(allv.inSet(vi->ident()))
            in_graph_recur_vars += vi->ident() ;
        }
        if(in_graph_recur_vars != EMPTY) {
          ruleSet others ;
          for(digraph::vertexSet::const_iterator
                vi=in_graph_recur_vars.begin();
              vi!=in_graph_recur_vars.end();
              ++vi)
              others = extract_rules(gr[*vi]) ;
          // but we must exclude those "generalize", "promote" and
          // "priority" internal rules from others
          variableSet only_vars = extract_vars(in_graph_recur_vars) ;
          only_vars += *varIter ;
          ruleSet subtract ;
          for(ruleSet::const_iterator ruleIter=others.begin();
              ruleIter!=others.end();++ruleIter) {
            variable target = *(ruleIter->targets().begin()) ;
            if(only_vars.inSet(target))
              subtract += *ruleIter ;
          }
          others -= subtract ;
          if(others != EMPTY) {
            target_rules += others ;
            recur_source_other_rules[*varIter] = others ;
          }
        }
      }

      // see if we can delete the variable in this graph
      if(let_it_go(gr,target_rules,*varIter))
        continue ;
      // otherwise we can delete the variable here
      delete_here += *varIter ;
    }

    // return the delete_here
    return delete_here ;    
  }

  void deleteInfoVisitor::visit(loop_compiler& lc) {
    // first in the loop, we don't delete the time_level variable
    variable tvar = variable(lc.tlevel) ;
    deleted_vars += tvar ;

    // we then get the variables shared between the
    // advance and collapse part of the loop
    map<int,variableSet>::const_iterator found ;
    variableSet shared ;
    found = loop_shared_table.find(lc.cid) ;
    FATAL(found == loop_shared_table.end()) ;
    shared = found->second ;

    // gather variables in the rotate list
    variableSet rotate_vars ;

    found = rotate_vtable.find(lc.cid) ;
    FATAL(found == rotate_vtable.end()) ;
    rotate_vars = found->second ;

    // we then gather deletion information for the collapse part
    variableSet working_vars = get_start_info(lc.collapse_gr,-lc.cid) ;
    // we delete the shared variables until the last iteration
    // finishes and inside the conditional super node
    working_vars -= shared ;
    // we don't delete variabls in rotate lists, we defer that until
    // the last iteration finishes
    working_vars -= rotate_vars ;
    
    looping_algr(working_vars,lc.collapse_gr,-lc.cid,0) ;
    
    // we gather deletion information for the advance part
    working_vars = get_start_info(lc.advance_gr,lc.cid) ;
    working_vars -= shared ;
    working_vars -= rotate_vars ;

    looping_algr(working_vars,lc.advance_gr,lc.cid,0) ;

    // we now schedule the deletion of shared and rotate_vars

    // first get the collapse node of this loop
    map<int,int>::const_iterator found2 ;
    found2 = loop_ctable.find(lc.cid) ;
    FATAL(found2 == loop_ctable.end()) ;
    int collapse_id = found2->second ;

    working_vars = variableSet(shared + rotate_vars) ;
    working_vars -= deleted_vars ;

    looping_algr(working_vars,lc.loop_gr,collapse_id,2) ;

  }

  void deleteInfoVisitor::visit(dag_compiler& dc) {
    variableSet working_vars = get_start_info(dc.dag_gr,dc.cid) ;
    looping_algr(working_vars,dc.dag_gr,dc.cid,1) ;    
  }

  
  void deleteInfoVisitor::visit(conditional_compiler& cc) {
    variableSet working_vars = get_start_info(cc.cond_gr,cc.cid) ;
    looping_algr(working_vars,cc.cond_gr,cc.cid,1) ;
  }

  // return the found loop node id, -1 if not found
  // loop_num specifies how many loops to find and include
  // along in the parent node path
  int deleteInfoVisitor::only_loop_alloc(variableSet& working_vars,
                                         int id, int loop_num) {
    // the top level 
    if(id == 0)
      return -1 ;

    if(loop_num < 0) {
      cerr << "ERROR calling only_loop_alloc() with loop_num < 0" << endl ;
      exit(-1) ;
    }
    // we check if we need to flip the id, if we are processing
    // the collapse graph of the loop
    if(id < 0) {
      if(inSet(loop_sn,-id))
        id = -id ;
      else {
        cerr << "ERROR calling only_loop_alloc() function,"
             << "with negative id number, only loop compiler "
             << "is allowed to have negative id number for its"
             << "collapse graph, the caller(id): " << id
             << " is not a loop compiler" << endl ;
        exit(-1) ;
      }
    }

    // if loop_num == 0, the only possible pass of such argument
    // to this function is inside the loop compiler, we do a check
    if(loop_num == 0) {
      if(!inSet(loop_sn,id)) {
        cerr << "ERROR calling only_loop_alloc() function,"
             << "caller(id): " << id << " is not a loop compiler"
             << " and therefore cannot use loop_num with 0!" << endl ;
        exit(-1) ;
      }
    }

    map<int,int>::const_iterator found ;
    map<int,variableSet>::const_iterator found_alloc ;
    int ret_id = -1 ;

    int pnode = id ;
    int found_loop_num = 0 ;

    // we search for loop_num loops
    while( (pnode != 0) && (found_loop_num != loop_num)) {
      found = pnode_table.find(pnode) ;
      FATAL(found == pnode_table.end()) ;
      pnode = found->second ;
      
      if(inSet(loop_sn,pnode))
        ++found_loop_num ;
    }
    
    if(found_loop_num == loop_num) {
      // get the loop's allocation
      found_alloc = loop_alloc_table.find(pnode) ;
      FATAL(found_alloc == loop_alloc_table.end()) ;
      
      working_vars &= found_alloc->second ;
      
      ret_id = pnode ;
    }
    
    return ret_id ;
  }

  // the looping algorithm for schedule deletion of variables
  // in the multilevel graph

  // working_vars is the initial unconstrained variableSet to
  // be evaluated, dg is the graph of a particular compiler
  // id is the compiler's id, start_loop_num is the starting
  // number for loop allocation constrain (0 for loop compiler,
  // 1 for dag & conditional compiler).
  void deleteInfoVisitor::looping_algr(const variableSet& working_vars,
                                       const digraph& dg,int id,
                                       int start_loop_num) {
    // it is important to use "+=" here, because of the loop,
    // we may have something in the delete_table[id] already
    // use "=" only will erase the previous contents.

    // *************************************************** //
    // IN ALL CODE HERE, WHENEVER WE WANT TO USE delete_table
    // AS LEFT VALUE, WE SHOULD ALWAYS USE "+=".
    // *************************************************** //
    
    // we can also skip this checking, it is just slightly
    // faster to do it here
    if(working_vars == EMPTY) {
      delete_table[id] += variableSet(EMPTY) ;
      return ;
    }
    
    // the looping algorithm
    int loop_num = start_loop_num ;
    int found_loop_id, last_loop_id ;
    variableSet processed_vars = working_vars ;
    variableSet ret ;
    
    found_loop_id = only_loop_alloc(processed_vars,id,loop_num) ;
    ret = gather_info(dg,processed_vars) ;
    last_loop_id = found_loop_id ;

    // edit table
    deleted_vars += ret ;
    delete_table[id] += ret ;

    variableSet remaining_vars = variableSet(working_vars - processed_vars) ;

    while(remaining_vars != EMPTY) {
      last_loop_id = found_loop_id ;
      ++loop_num ;
      found_loop_id = only_loop_alloc(remaining_vars,id,loop_num) ;
      ret = gather_info(dg,remaining_vars) ;

      if(last_loop_id != -1) {
        // put the ret into the collapse node of last_loop_id loop
        if(ret != EMPTY) {
          // we then find out the collapse super node inside this loop node
          map<int,int>::const_iterator found ;
          found = loop_ctable.find(last_loop_id) ;
          FATAL(found == loop_ctable.end()) ;
          int collapse_id = found->second ;
          
          deleted_vars += ret ;
          delete_table[collapse_id] += ret ;
        }
      }
      processed_vars += remaining_vars ;
      remaining_vars = variableSet(working_vars - processed_vars) ;
    }
  }
                                              
  // given a rule set and a graph
  // return true if there is only one super node and
  // all other rules have path to the super node, or if
  // it is only one super node
  // return false otherwise
  bool deleteInfoVisitor::let_it_go(const digraph& gr, ruleSet rules,
                                    const variable& v) {
    if(rules == EMPTY) {
      // if the target rules are empty, we look back to the rules
      // that generate the variable, if it is only one super node
      // we let it go, meaning we don't delete the variable here
      // this only applies to non loop rotate variables
      digraph grt = gr.transpose() ;
      ruleSet source_rules = extract_rules(grt[v.ident()]) ;
      if(source_rules.size() == 1) {
        if(is_super_node(source_rules.begin()))
          if(inSet(graph_sn,get_supernode_num(*source_rules.begin())))
            if(!all_rot_vars.inSet(v))
              return true ;
      }
      return false ;
    }
    
    int supernode_num = 0 ;
    int supernode_id = -1 ;
    rule supernode ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        // we eliminate the recursive super node
        // by looking up the graph_sn set
        if(inSet(graph_sn,get_supernode_num(*ruleIter))) {
          supernode = *ruleIter ;
          supernode_id = ruleIter->ident() ;
          ++supernode_num ;
        }
      }
    }
    if( (rules.size() == 1) && (supernode_num == 1)) {
      // check to see if it is a collapse node
      // if it is, we delete it in the graph
      if(inSet(col_sn,get_supernode_num(supernode)))
        return false ;
      return true ;
    }
    if(supernode_num != 1)
      return false ;
    // the remaining case must be there are only one super node
    // and at least one non super node, we check if there is path
    // from each non super node to the super node in the graph
    FATAL(supernode_id < 0) ;
    rules -= rule(supernode_id) ;
    bool ret = true ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(!has_path(gr,ruleIter->ident(),supernode_id)) {
        ret = false ;
        break ;
      }
    }
    
    return ret ;
  }

  //////////////////////////////////////////////////////////
  // deleteGraphVisitor
  //////////////////////////////////////////////////////////
  void deleteGraphVisitor::edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // first we get the info in delete_table, find out which variables
    // need to be deleted in the graph
    map<int,variableSet>::const_iterator found ;
    found = delete_table.find(id) ;
    FATAL(found == delete_table.end()) ;
    variableSet delete_vars = found->second ;

    if(delete_vars == EMPTY) return ;
    
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;
    
    // get all possible variables in the graph
    variableSet working_vars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      working_vars += ri->sources() ;
    
    working_vars += extract_vars(allvertices) ;

    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // we check if the variable is in the delete_vars
      // if it is not in the delete_vars, we don't delete it here
      if(!delete_vars.inSet(*varIter))
        continue ;
      // now we decorate the graph to include the deletion
      // we first obtain the ident of rules that are reachable
      // from the variable
      digraph::vertexSet target_rules = gr[varIter->ident()] ;

      // if the variable is a recurrence target variable,
      // we should also get possible rules that the recurrence
      // source variable leads to 
      map<variable,ruleSet>::const_iterator found ;
      found = recur_source_other_rules.find(*varIter) ;
      
      if(found != recur_source_other_rules.end()) {
        ruleSet others = found->second ;
        for(ruleSet::const_iterator ruleIter=others.begin();
            ruleIter!=others.end();++ruleIter)
          target_rules += ruleIter->ident() ;
      }

      // if target_rules is empty, then the variable is not
      // present in this graph, we need a search to find out
      // the target_rules
      if(target_rules == EMPTY) {
        for(ruleSet::const_iterator ruleIter=rules.begin();
            ruleIter!=rules.end();++ruleIter) {
          variableSet sources = ruleIter->sources() ;
          if(sources.inSet(*varIter))
            target_rules += ruleIter->ident() ;
        }
      }
      // if target_ruls is still empty, then the variable must
      // be in the graph, but it is a leaf node, we append the
      // deletion rule after the variable itself.
      if(target_rules == EMPTY)
        target_rules += varIter->ident() ;
       
      // we create a rule for deletion
      variable sv("DESTROY") ;
      rule delete_rule = create_rule(sv,*varIter,"DELETE") ;
      // edit graph
      gr.add_edges(target_rules,delete_rule.ident()) ;
      // modify the rulecomp_map
      rcm[delete_rule] = new free_var_compiler(*varIter) ;
      // remove variable in alloc_vars
      delete_vars -= *varIter ;
    }
    // if delete_vars is not empty at this stage, it means
    // there must have some variables need to be deleted in this
    // graph but these variables do not belong to this graph

    if(delete_vars != EMPTY) {
      // first find out all the leaf node in the graph
      digraph::vertexSet leaf_node ;
      for(digraph::vertexSet::const_iterator vi=allvertices.begin();
          vi!=allvertices.end();++vi) {
        if(gr[*vi] == EMPTY)
          leaf_node += *vi ;
      }
      
      for(variableSet::const_iterator vi=delete_vars.begin();
          vi!=delete_vars.end();++vi) {
        // we create a rule for deletion
        variable sv("DESTROY") ;
        rule delete_rule = create_rule(sv,*vi,"DELETE") ;
        // we make the rule happen at last in the graph
        // we connect all the leaf node in the graph
        // to this deletion rule
        gr.add_edges(leaf_node,delete_rule.ident()) ;
        // modify the rulecomp_map
        rcm[delete_rule] = new free_var_compiler(*vi) ; 
      }
    }
  }

  void deleteGraphVisitor::visit(loop_compiler& lc) {
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,-lc.cid) ;
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }
  
  void deleteGraphVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }
  
  void deleteGraphVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }

} // end of namespace Loci

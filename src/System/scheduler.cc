#include "scheduler.h"
#include <digraph.h>
#include <fact_db.h>
#include <execute.h>
#include <depend_graph.h>
#include <Map.h>

using namespace Loci ;

#include <map>
using std::map ;
#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <list>
using std::list ;

using std::pair ;
using std::make_pair ;

//#define VERBOSE

namespace {
class decompose_graph ;

class decompose_graph {
  friend class rule_calculator ;
  friend class loop_calculator ;
  variable start,finish ;
  struct reduce_info {
    rule unit_rule ;
    ruleSet apply_rules ;
  } ;
  map<variable,reduce_info> reduce_map ;
  variableSet reduce_vars ;
  map<rule,rule> apply_to_unit ;

  ruleSet looping_rules ;
  ruleSet conditional_rules ;
  variableSet conditional_vars ;
  map<variable,ruleSet> conditional_map ;
  
  rule top_level_rule ;
  struct super_node_info {
    digraph graph ;
    variableSet sources,targets ;
  } ;
  map<rule, super_node_info> supermap ;
  ruleSet supernodes ;
  rule create_supernode(digraph &g, digraph::nodeSet ns,
                        variable cond_var = variable()) ;
  ruleSet recursive_supernodes ;
  ruleSet looping_supernodes ;

  map<rule, rule_calculator *> rule_process ;
  
public:
  decompose_graph(digraph dg,digraph::nodeSet sources, digraph::nodeSet targets) ;
  void existential_analysis(fact_db &facts) ;
  executeP execution_schedule(fact_db &facts) ;
} ;

class rule_calculator {
protected:
  decompose_graph *graph_ref ;
  rule_calculator *calc(const rule &r) { return graph_ref->rule_process[r] ; }
public:
  virtual void set_var_existence(fact_db &facts) = 0 ;
  virtual void process_var_requests(fact_db &facts) = 0 ;
  virtual executeP create_execution_schedule(fact_db &facts) = 0 ;
} ;

class loop_calculator : public rule_calculator {
  digraph dag ;
  vector<rule> rule_schedule ;
  vector<rule> collapse ;
  variable cond_var ;
  vector<rule> advance ;
  variableSet advance_vars ;
  time_ident tlevel ;
  variable output ;
  bool output_present ;
public:
  loop_calculator(decompose_graph *gr, digraph gin) ;
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;


digraph::nodeSet visit_nodes(digraph dg,digraph::nodeSet begin) {

  digraph::nodeSet visit = begin ;
  digraph::nodeSet visited ;
  digraph::nodeSet::const_iterator ni ;

  // keep visiting nodes until no new nodes are found
  while(visit != EMPTY) {
    digraph::nodeSet newnodes ;
    // visit all the nodes that this graph leads to
    for(ni=visit.begin();ni!=visit.end();++ni)
      newnodes += dg[*ni] ;
    // update visit, but don't re-visit a node
    visit = newnodes - visited ;
    visited = visited + newnodes ;
  }
  return visited ;
}

digraph::nodeSet visit_nodes_exclusive(digraph dg, digraph::nodeSet begin) {
  digraph dgt = dg.transpose() ;
  digraph::nodeSet visited, visit ;

  visit = begin ;
  visited = visit ;

  // keep visiting nodes until no new nodes are found
  while(visit != EMPTY) {
    digraph::nodeSet::const_iterator ni,nii,no ;
    digraph::nodeSet newnodes, not_visited = ~visited ;
    // visit all new nodes and check to see if they contain paths to
    // new candidate nodes
    for(ni=visit.begin();ni!=visit.end();++ni) 
      for(nii=dg[*ni].begin();nii!=dg[*ni].end();++nii) { // loop over edges
        // if this edge leads to a node that can only be reached through
        // the visited set, then it is a new node to be visited
        digraph::nodeSet out_bound_nodes = dgt[*nii] & not_visited ;

        // Check to see if out of bound nodes loop exclusively back to this
        // node, if so then add it to the list of new nodes.
        // This is a little bit of a hack since the graphs analyzed may
        // contain loops no larger than 2 edges and 2 vertexes.
        bool flg = true ;
        digraph::nodeSet local_not_visit = ~(visited + interval(*nii,*nii));
        for(no=out_bound_nodes.begin();no!=out_bound_nodes.end();++no) {
          flg = flg && (dgt[*no] & local_not_visit) == EMPTY ;
        }
        if(flg) 
          newnodes += *nii ;
      }
    // next time we will visit the nodes found in this iteration
    // update our list of visited nodes.
    visit = newnodes - visited ;
    visited += newnodes ;
  }
  return visited ;
}


// Create a schedule for traversing a directed acyclic graph.  This schedule
// may be concurrent, or many nodes of the graph may be visited at each
// step of the schedule  If the graph contains cycles, the schedule may
// not include all of the nodes in the graph.
vector<digraph::nodeSet> schedule_dag(const digraph &g,
                                      digraph::nodeSet start_nodes = EMPTY,
                                      digraph::nodeSet only_nodes =
                                      interval(UNIVERSE_MIN,UNIVERSE_MAX)) {
  vector<digraph::nodeSet> schedule ;
  digraph gt = g.transpose() ;
  // First schedule any nodes that have no edges leading into them and have
  // not been scheduled previously (in start nodes)
  digraph::nodeSet working = g.get_source_nodes() -
    (g.get_target_nodes()+start_nodes) ;
  if(working != EMPTY)
    schedule.push_back(working) ;
  // visited nodes are all nodes that have already been scheduled
  digraph::nodeSet visited_nodes = start_nodes + working ;
  // In the beginning our working set are all scheduled nodes
  working = visited_nodes ;
  while(working != EMPTY) {
    // While we have nodes to work on, compute additional nodes that
    // can be scheduled
    digraph::nodeSet new_nodes ;
    digraph::nodeSet::const_iterator ni ;
    // loop over working set and create a list of candidate nodes
    for(ni=working.begin();ni != working.end(); ++ni)
      new_nodes += g[*ni] ;
    // If a node has already been scheduled it can't be scheduled again,
    // so remove visited nodes
    new_nodes = new_nodes - visited_nodes    ;
    // We only schedule nodes that are also in the only_nodes set
    working = new_nodes & only_nodes ;
    new_nodes = EMPTY ;
    // Find any node from this working set that has had all nodes leading
    // to it scheduled
    for(ni=working.begin();ni != working.end(); ++ni) 
      if((gt[*ni] & visited_nodes) == gt[*ni])
        new_nodes += *ni ;
    working = new_nodes ;
    // and these new nodes to the schedule
    if(new_nodes != EMPTY)
      schedule.push_back(new_nodes) ;
    // update visited nodes set to include scheduled nodes
    visited_nodes += new_nodes ;
  }
  return schedule ;
}


// Create a sequential ordering of rules based on the concurrent dag schedule
void extract_rule_sequence(vector<rule> &rule_seq,
                           const vector<digraph::nodeSet> &v) {
  vector<digraph::nodeSet>::const_iterator i ;
  for(i=v.begin();i!=v.end();++i) {
    ruleSet rules = extract_rules(*i) ;
    ruleSet::const_iterator ri ;
    for(ri=rules.begin();ri!=rules.end();++ri)
      rule_seq.push_back(*ri) ;
  }
}
        
rule make_super_rule(variableSet sources, variableSet targets,
                     variable cond = variable()) {
  FATAL(targets == EMPTY) ;
  static int super_node_number = 0 ;
  ostringstream oss ;
  oss << "source("<<sources << "),target(" << targets << ")," ;
  if(cond != variable()) 
    oss<< "conditional(" << cond << ")," ;
  oss << "qualifier(SN" << super_node_number++ << ")" ;
   
  return rule(oss.str()) ;
}

rule make_rename_rule(variable new_name, variable old_name) {
  ostringstream oss ;
  oss << "source(" << old_name << "),target(" << new_name
      << "),qualifier(rename)"  ;
  return rule(oss.str()) ;
}

// Create variable types in fact database
// This is a necessary precursor to binding rules to the database
void set_var_types(fact_db &facts, const digraph &dg) {

  // Get all the variables and rules represented in the graph
  variableSet all_vars = extract_vars(dg.get_all_nodes()) ;
  ruleSet     all_rules = extract_rules(dg.get_all_nodes()) ;

  // extract the qualified rules, these are rules that are
  // automatically generated by the system.  Since these rules
  // can not provide type information directly, they are
  // singled out so that they can be handled as a separate case.
  ruleSet qualified_rules ;
  for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
    if(ri->type() == rule::INTERNAL)
      qualified_rules += *ri ;
  }

  // We need the transpose of the graph in order to find the rules that
  // generate a particular variable
  digraph dgt = dg.transpose() ;
  
  vector<pair<variable,variable> > rename_vars ;

  // Loop through the variables and type them when there is an appropriate
  // rule specification.  Deal with typing rules that rename their targets
  // separately
  for(variableSet::const_iterator
        vi=all_vars.begin();vi!=all_vars.end();++vi) {

    // only deal with rules given to the system by the user
    ruleSet rs = extract_rules(dgt[(*vi).ident()] - qualified_rules) ;
    if(rs == EMPTY)
      continue ;

    // storage for rename specification if one is found
    pair<variable,variable> rename_pair ;
    bool rename = false ;

    // Check all rules generating this variable, if the rules specify
    // renaming a variable for the target, then by default the target
    // will be an identical type to that of the previous variable.
    for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) {
      set<vmap_info>::const_iterator vmsi ;
      for(vmsi=ri->get_info().desc.targets.begin();
          vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
        if(vmsi->assign.size() != 0) 
          for(int i=0;i<vmsi->assign.size();++i) {
            variable new_name = vmsi->assign[i].first ;
            
            if(new_name == *vi) {
              variable old_name = vmsi->assign[i].second ;
              if(rename && old_name != rename_pair.second) {
                cerr << "rule " << *ri << endl 
                     << " is not consistent with other rules w/respect to "
                     << "renaming." << endl ;
                exit(-1) ;
              }

              rename_pair = make_pair(new_name,old_name) ;
              rename = true ;
            }
          }
    }

    // If a rename variable was found, save this information for later
    // use and continue ;
    if(rename) {
      rename_vars.push_back(rename_pair) ;
      continue ;
    }

    // Otherwise pick the first rule in the set of rules that generates it
    // and query the rule for the type of variable *vi.  Set the variable
    // type appropriately in the fact database.

    rule pick = *rs.begin() ;
    const rule_impl::info &info = pick.get_info().desc ;
    storeRepP st = pick.get_info().rule_impl->get_store(*vi) ;
    facts.set_variable_type(*vi,st) ;
    
  }

  // We now have to deal with internally generated rules that do operations
  // such as promoting variables from one time level to another,
  // generalizing a variable from a specific time to general time
  // (e.g. q{n}<-q{n=0}) and rules that create a priority relationship
  // from priority variables (e.g. qf<-degenerate::qf).
  // In addition we will deal with variables that may be renamed.
  // To do this, we must tell the fact database that variables are
  // synonyms (meaning that the two variable names are completely identical)
  // and variables that are aliases (meaning that the variable is exclusive
  // by name).  Since these variables have to be defined from
  // variables that are already typed and that these rules may
  // become chained, we need to make sure that we proceed in the
  // proper order.  We do this by creating a subgraph of these
  // internal rules, then we perform a topological sort on this
  // sub-graph to determine the order of typing

  // This is the internal rule graph
  digraph irg ;
  // place all promote, generalize, and priority rules into irg
  for(ruleSet::const_iterator
        ri = qualified_rules.begin() ; ri != qualified_rules.end();++ri) {
    if(ri->get_info().qualifier() == "promote" ||
       ri->get_info().qualifier() == "generalize" ||
       ri->get_info().qualifier() == "priority") {
      irg.add_edges((*ri).ident(),dg[(*ri).ident()]) ;
      irg.add_edges(dgt[(*ri).ident()],(*ri).ident()) ;
    } else if(ri->get_info().qualifier() == "looping") {
      // create the iteration variable type
      time_ident tl = ri->source_time() ;
      variable v(tl) ;
      param<int> timevar ;
      storeRepP st = timevar ;
      facts.set_variable_type(v,st) ;
    }
  }

  // Also add rules to represent renaming relationships
  for(int i=0;i<rename_vars.size();++i) {
    variable new_var = rename_vars[i].first ;
    variable old_var = rename_vars[i].second ;
    
    rule rr = make_rename_rule(new_var,old_var) ;
    irg.add_edge(old_var.ident(),rr.ident()) ;
    irg.add_edge(rr.ident(),new_var.ident()) ;
  }  

  // Perform a topological sort on the interal rule graph
  vector<digraph::nodeSet> components =
    component_sort(irg).get_components() ;

  for(int i=0;i<components.size();++i) {
    // Note, recursion should not occur ammong the internal rules
    FATAL(components[i].size()!=1) ;
    ruleSet rs = extract_rules(components[i]) ;
    // Note: components may be rules or variables, thus this check
    if(rs != EMPTY) {
      rule r = *rs.begin() ;
      FATAL(r.sources().size() != 1) ;
      FATAL(r.targets().size() != 1) ;
      // get source and target variables
      variable s = (*(r.sources().begin())) ;
      variable t = (*(r.targets().begin())) ;
      // A rename rule uses alias, while all others use synonym relationships
      if(r.get_info().qualifier() == "rename")
        facts.alias_variable(s,t) ;
      else
        facts.synonym_variable(s,t) ;
    }
  }

  // Check to make sure there are no type conflicts, loop over all
  // rules that are not internally generated.
  bool type_error = false ;
  ruleSet::const_iterator ri ;
  ruleSet rs = all_rules ;
  rs -= qualified_rules ;
  for(ri=rs.begin();ri!=rs.end();++ri) {
    variableSet varcheck ;
    const rule_impl::info &finfo = ri->get_info().desc ;

    // Collect all variables for which are actually read or written in the class
    set<vmap_info>::const_iterator i ;
    for(i=finfo.sources.begin();i!=finfo.sources.end();++i) {
      for(int j=0;j<i->mapping.size();++j)
        varcheck += i->mapping[j] ;
      varcheck += i->var ;
    }
    for(i=finfo.targets.begin();i!=finfo.targets.end();++i) {
      for(int j=0;j<i->mapping.size();++j)
        varcheck += i->mapping[j] ;
      varcheck += i->var ;
      for(int k=0;k<i->assign.size();++k) {
        varcheck -= i->assign[k].first ;
        varcheck += i->assign[k].second ;
      }
    }

    variableSet::const_iterator vi ;
    for(vi = varcheck.begin();vi!=varcheck.end();++vi) {
      storeRepP rule_type = ri->get_rule_implP()->get_store(*vi)->getRep() ;
      storeRepP fact_type = facts.get_variable(*vi)->getRep() ;
      if(typeid(*rule_type) != typeid(*fact_type)) {
        cerr << "variable type mismatch for variable " << *vi << " in rule "
             << *ri << endl ;
        cerr << "fact database has type " << typeid(*fact_type).name() << endl ;
        cerr << "rule has type " << typeid(*rule_type).name() << endl ;
        type_error = true ;
      }
    }
  }
  if(type_error)
    exit(-1) ;
}


entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts) {
  variableSet::const_iterator vi ;
  entitySet sources = ~EMPTY ;
  for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
    sources &= facts.variable_existence(*vi) ;
  vector<variableSet>::const_reverse_iterator mi ;
  for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
    entitySet working = ~EMPTY ;
    for(vi=mi->begin();vi!=mi->end();++vi) {
      FATAL(!facts.is_a_Map(*vi)) ;
      working &= facts.preimage(*vi,sources).first ;
    }
    sources = working ;
  }
  return sources ;
}


entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                            entitySet compute) {
  vector<variableSet>::const_iterator mi ;
  for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
    WARN(mi->size() != 1) ;
    variable v = *(mi->begin()) ;
    FATAL(!facts.is_a_Map(v)) ;
    compute = facts.image(v,compute) ;
  }
  return compute ;
}


void existential_rule_analysis(rule f, fact_db &facts) {

  FATAL(f.type() == rule::INTERNAL) ;

  entitySet sources = ~EMPTY ;
  entitySet constraints = ~EMPTY ;
  entitySet targets = EMPTY ;
  const rule_impl::info &finfo = f.get_info().desc ;
  set<vmap_info>::const_iterator si ;
  for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
    sources &= vmap_source_exist(*si,facts) ;
  for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
    constraints &= vmap_source_exist(*si,facts) ;
  if(finfo.constraints.begin() != finfo.constraints.end())
    if((sources & constraints) != constraints) {
      cerr << "Warning, rule " << f <<
        "can not supply all entities of constraint" << endl ;
      exit(-1) ;
    }
  sources &= constraints ;
  for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
    targets = vmap_target_exist(*si,facts,sources) ;
    const variableSet &tvars = si->var ;
    variableSet::const_iterator vi ;
    for(vi=tvars.begin();vi!=tvars.end();++vi)
      facts.set_existential_info(*vi,f,targets) ;
  }
#ifdef VERBOSE
  cout << "rule " << f << " generating " << targets << endl ;
#endif
}
typedef map<variable,entitySet> vdefmap ;

entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                               fact_db &facts) {
  variableSet::const_iterator vi ;
  entitySet targets ;
  for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
    FATAL(tvarmap.find(*vi) == tvarmap.end()) ;
    targets |= tvarmap.find(*vi)->second ;
  }
  for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
    facts.variable_request(*vi,targets) ;
  
  vector<variableSet>::const_reverse_iterator mi ;
  for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
    entitySet working = EMPTY ;
    for(vi=mi->begin();vi!=mi->end();++vi) {
      FATAL(!facts.is_a_Map(*vi)) ;
      working |= facts.preimage(*vi,targets).second ;
    }
    targets = working ;
  }
  return targets ;
}

void vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                          entitySet compute) {
  vector<variableSet>::const_iterator mi ;
  variableSet::const_iterator vi ;
  for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
    entitySet working ;
    for(vi=mi->begin();vi!=mi->end();++vi) {
      FATAL(!facts.is_a_Map(*vi)) ;
      working |= facts.image(*vi,compute) ;
    }
    compute = working ;
  }
  for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
    facts.variable_request(*vi,compute) ;
}

entitySet process_rule_requests(rule f, fact_db &facts) {

  FATAL(f.type() == rule::INTERNAL) ;

  variableSet targets = f.targets() ;
  
  variableSet::const_iterator vi ;
  vdefmap tvarmap ;
  for(vi=targets.begin();vi!=targets.end();++vi) {
    if(vi->get_info().name == "OUTPUT") 
      facts.variable_request(*vi,facts.variable_existence(*vi)) ;
    tvarmap[*vi] = facts.get_variable_request(f,*vi) ;
  }

  const rule_impl::info &finfo = f.get_info().desc ;
  set<vmap_info>::const_iterator si ;
  entitySet compute,isect = ~EMPTY ;
  for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
    entitySet tmp = vmap_target_requests(*si,tvarmap,facts) ;
    compute |= tmp ;
    isect &= tmp ;
  }

  if(isect != compute) {
    entitySet working = compute ;
    vector<variableSet>::const_reverse_iterator mi ;
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      for(mi=si->mapping.rbegin();mi!=si->mapping.rend();++mi) {
        entitySet tmp ;
        for(vi=mi->begin();vi!=mi->end();++vi)
          tmp |= facts.image(*vi,working) ;
        working = tmp ;
      }
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        facts.variable_request(*vi,working) ;
        WARN(facts.get_variable_request(f,*vi) != working) ;
      }
    }
  }

  for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
    vmap_source_requests(*si,facts,compute) ;

  for(vi=finfo.conditionals.begin();vi!=finfo.conditionals.end();++vi)
    facts.variable_request(*vi,compute) ;

#ifdef VERBOSE
  cout << "rule " << f << " computes over " << compute << endl ;
#endif
  return compute ;
}

class execute_rule : public execute_modules {
  rule_implP rp ;
  rule rule_tag ;
  sequence exec_seq ;
public:
  execute_rule(rule fi, sequence seq, fact_db &facts) ;
  virtual void execute(fact_db &facts) ;
  virtual void Print(std::ostream &s) const ;
} ;

execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts)
{
  rp = fi.get_rule_implP() ;
  rule_tag = fi ;
  rp->initialize(facts) ;
  exec_seq = seq ;
}

void execute_rule::execute(fact_db &facts) {
  rp->compute(exec_seq) ;
}

void execute_rule::Print(ostream &s) const {
  s << rule_tag << " over sequence " << exec_seq << endl ;
}

class execute_loop : public execute_modules {
  executeP collapse, advance ;
  variable cvar ;
  time_ident tlevel ;
public:
  execute_loop(const variable &cv,
               const executeP &col, const executeP &adv,
               const time_ident &tl) :
    cvar(cv),collapse(col),advance(adv),tlevel(tl)
  { warn(col==0 || advance==0) ;}
  virtual void execute(fact_db &facts) ;
  virtual void Print(std::ostream &s) const ;
} ;

void execute_loop::execute(fact_db &facts) {
  param<bool> test ;
  test = facts.get_variable(cvar) ;
  fact_db::time_infoP tinfo = facts.get_time_info(tlevel) ;
  facts.initialize_time(tinfo) ;
  for(;;) {
    collapse->execute(facts) ;
    if(*test) {
      facts.close_time(tinfo) ;
      return ;
    }
    advance->execute(facts) ;
    facts.advance_time(tinfo) ;
  }
}

void execute_loop::Print(ostream &s) const {
  s << "Perform loop for time level "<< tlevel << endl ;
  s << "--compute collapse rule, conditional on " << cvar << endl ;
  s << "--collapse iteration {" << tlevel << "}"<< endl ;
  collapse->Print(s) ;
  s << "--advance iteration {" << tlevel << "}" << endl ;
  advance->Print(s) ;
  s << "end of loop for time level " << tlevel << endl ;
}


  
class execute_conditional : public execute_modules {
  executeP conditional ;
  variable cvar ;
public:
  execute_conditional(const executeP &cond, const variable &cv) :
    conditional(cond),cvar(cv)
  { warn(cond==0) ;}
  virtual void execute(fact_db &facts) ;
  virtual void Print(std::ostream &s) const ;
} ;

void execute_conditional::execute(fact_db &facts) {
  param<bool> test ;
  test = facts.get_variable(cvar) ;

  if(*test) {
    conditional->execute(facts) ;
  }
}

void execute_conditional::Print(ostream &s) const {
  s << "--compute rule if conditional " << cvar << " true." << endl ;
  conditional->Print(s) ;
  s << "--end conditional" << endl ;
}


  
class allocate_all_vars : public execute_modules {
public:
  virtual void execute(fact_db &facts) ;
  virtual void Print(ostream &s) const ;
} ;

void allocate_all_vars::execute(fact_db &facts) {
  variableSet vars = facts.get_typed_variables() ;
  variableSet::const_iterator vi,vii ;
  set<time_ident> time_set ;
  for(vi=vars.begin();vi!=vars.end();++vi)
    time_set.insert(vi->time()) ;
  set<time_ident>::const_iterator si ;
  for(si=time_set.begin();si!=time_set.end();++si) {
    fact_db::time_infoP tip = facts.get_time_info(*si) ;
    if(!tip->rotate_lists.empty()) {
      for(list<list<variable> >::const_iterator llv = tip->rotate_lists.begin();
          llv != tip->rotate_lists.end();++llv) {
        entitySet time_space ;
        for(list<variable>::const_iterator lv=llv->begin();lv!=llv->end();++lv) {
          variableSet aliases = facts.get_aliases(*lv) ;
          for(vii=aliases.begin();vii!=aliases.end();++vii)
            time_space += facts.get_variable_requests(*vii) ;
        }
        for(list<variable>::const_iterator lv=llv->begin();lv!=llv->end();++lv) {
          facts.variable_request(*lv,time_space) ;
        }
      }
    }
  }
  
  for(vi=vars.begin();vi!=vars.end();++vi) {
    storeRepP srp = facts.get_variable(*vi) ;
    if(srp->domain() == EMPTY) {
      variableSet aliases = facts.get_aliases(*vi) ;
      entitySet all_requests ;
      for(vii=aliases.begin();vii!=aliases.end();++vii)
        all_requests += facts.get_variable_requests(*vii) ;
      srp->allocate(all_requests) ;
    }
  }
}

void allocate_all_vars::Print(ostream &s) const {
  s << "allocate all variables" << endl ;
}


// rule calculator for rule with concrete implementation
class impl_calculator : public rule_calculator {
  rule impl ;  // rule to implement

  // existential analysis info
  sequence exec_seq ;
public:
  impl_calculator(decompose_graph *gr, rule r)  { graph_ref = gr; impl=r;}
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

void impl_calculator::set_var_existence(fact_db &facts) {
  existential_rule_analysis(impl,facts) ;
}

void impl_calculator::process_var_requests(fact_db &facts) {
  exec_seq = sequence(process_rule_requests(impl,facts)) ;
}

executeP impl_calculator::create_execution_schedule(fact_db &facts) {
#ifndef DEBUG
  if(exec_seq.size() == 0)
    return executeP(0) ;
  else
#endif
    return new execute_rule(impl,exec_seq,facts) ;
}

// apply rule calculator 
class apply_calculator : public rule_calculator {
  rule apply,unit_tag ;  // rule to applyement

  // existential analysis info
  sequence exec_seq ;
public:
  apply_calculator(decompose_graph *gr, rule r, rule ut)
  { graph_ref = gr; apply=r; unit_tag = ut ; }
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

void apply_calculator::set_var_existence(fact_db &facts) {
}

void apply_calculator::process_var_requests(fact_db &facts) {

  vdefmap tvarmap ;
  variableSet targets = apply.targets() ;
  variableSet sources = apply.sources() ;

  fatal(targets.size() != 1) ;
  variable tvar = *(targets.begin()) ;

  if(facts.get_variable(tvar)->RepType() == Loci::PARAMETER) 
    tvarmap[tvar] = facts.variable_existence(tvar) ;
  else
    tvarmap[tvar] = facts.get_variable_request(unit_tag,tvar) ;
    

  const rule_impl::info &finfo = apply.get_info().desc ;
  set<vmap_info>::const_iterator si ;
  entitySet compute ;
  for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
    compute |= vmap_target_requests(*si,tvarmap,facts) ;
  }

  for(si=finfo.targets.begin();si!=finfo.targets.end(); ++si) {
    variableSet::const_iterator vi ;
    entitySet comp = compute ;
    vector<variableSet>::const_iterator mi ;
    for(mi=si->mapping.begin();mi!=si->mapping.end();++mi) {
      entitySet working ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!facts.is_a_Map(*vi)) ;
        working |= facts.image(*vi,comp) ;
      }
      comp = working ;
    }
    for(vi=si->var.begin();vi!=si->var.end();++vi) {
      if((comp - facts.variable_existence(*vi)) != EMPTY) {
        cerr << "ERROR: Apply rule " << apply <<  endl
             << " output mapping forces application to entities where unit does not exist." << endl ;
        cerr << "error occurs for entities " <<
          entitySet(comp-facts.variable_existence(*vi)) << endl ;
        cerr << "error occurs when applying to variable " << *vi << endl;
        cerr << "error is not recoverable, terminating scheduling process"
             << endl ;
        exit(-1) ;
      }
      facts.variable_request(*vi,comp) ;
    }
  }

  
  entitySet srcs = ~EMPTY ;
  entitySet cnstrnts = ~EMPTY ;
  for(si=finfo.sources.begin();si!=finfo.sources.end();++si)
    srcs &= vmap_source_exist(*si,facts) ;
  for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si)
    cnstrnts &= vmap_source_exist(*si,facts) ;
  if(finfo.constraints.begin() != finfo.constraints.end())
    if((srcs & cnstrnts) != cnstrnts) {
      cerr << "Warning, reduction rule:" << apply
           << "cannot supply all entities of constraint" << endl ;
      cerr << "constraints = " <<cnstrnts << endl ;
      entitySet sac = srcs & cnstrnts ;
      cerr << "srcs & constraints = " << sac << endl ;
      exit(-1) ;
    }
  srcs &= cnstrnts ;

  // now trim compute to what can be computed.
  compute &= srcs ;
  
  exec_seq = sequence(compute) ;
    
  for(si=finfo.sources.begin();si!=finfo.sources.end();++si) 
    vmap_source_requests(*si,facts,compute) ;
    
#ifdef VERBOSE
  cout << "rule " << apply << " computes over " << compute << endl ;
#endif
}

executeP apply_calculator::create_execution_schedule(fact_db &facts) {
#ifndef DEBUG
  if(exec_seq.size() == 0)
    return executeP(0) ;
  else
#endif
    return new execute_rule(apply,exec_seq,facts) ;
}

// rule calculator for single rule recursion
class impl_recurse_calculator : public rule_calculator {
  rule impl ;  // rule to implement
  struct fcontrol {
    list<entitySet> control_list ;
    map<variable,entitySet> generated ;
    bool use_constraints ;
    entitySet nr_sources, constraints ;
    struct mapping_info {
      entitySet total ;
      vector<MapRepP> mapvec ;
      vector<variable> mapvar ;
      variable v ;
    } ;
    vector<mapping_info> recursion_maps, target_maps ;
  } ;
  fcontrol  control_set ;
  sequence fastseq ;
public:
  impl_recurse_calculator(decompose_graph *gr, rule r)
  { graph_ref = gr; impl = r ;}
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

void impl_recurse_calculator::set_var_existence(fact_db &facts) {
  variableSet::const_iterator vi ;
  variableSet tvars ;
  tvars = impl.targets() ;
  fcontrol &fctrl = control_set ;
  const rule_impl::info &finfo = impl.get_info().desc ;
  warn(impl.type() == rule::INTERNAL) ;
  entitySet sources = ~EMPTY ;
  entitySet constraints = ~EMPTY ;
  set<vmap_info>::const_iterator si ;
  for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
    if((si->var & tvars) == EMPTY)
      sources &= vmap_source_exist(*si,facts) ;
    else {
      tvars += si->var ;
      int num_maps = si->mapping.size() ;
      vector<variableSet::const_iterator> miv(num_maps) ;
      for(int i=0;i<num_maps;++i)
        miv[i] = si->mapping[i].begin() ;
        
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        int i = 0 ;
        do {
          fcontrol::mapping_info minfo ;
          for(int j=0;j!=num_maps;++j) {
            FATAL(!facts.is_a_Map(*(miv[j]))) ;
            MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;  
            minfo.mapvec.push_back(m) ;
            minfo.mapvar.push_back(*(miv[j])) ;
          }
          minfo.v = *vi ;
            
          for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
            miv[i] = si->mapping[i].begin() ;
          }
          fctrl.recursion_maps.push_back(minfo) ;
        } while(i!=num_maps) ;
      }
    }
  } 
  for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
    int num_maps = si->mapping.size() ;
    vector<variableSet::const_iterator> miv(num_maps) ;
    for(int i=0;i<num_maps;++i)
      miv[i] = si->mapping[i].begin() ;
      
    for(vi=si->var.begin();vi!=si->var.end();++vi) {
      int i = 0 ;
      do {
        fcontrol::mapping_info minfo ;
        for(int j=0;j!=num_maps;++j) {
          FATAL(!facts.is_a_Map(*(miv[j]))) ;
          MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
          minfo.mapvec.push_back(m) ;
          minfo.mapvar.push_back(*(miv[j])) ;
        }
        minfo.v = *vi ;
          
        for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
          miv[i] = si->mapping[i].begin() ;
        }
        fctrl.target_maps.push_back(minfo) ;
      } while(i!=num_maps) ;
    }
  }
  fctrl.use_constraints = false ;
  for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
    warn((si->var & tvars) != EMPTY) ;
    constraints &= vmap_source_exist(*si,facts) ;
    fctrl.use_constraints = true ;
  }
  fctrl.nr_sources = sources ;
  fctrl.constraints = constraints ;
  warn(fctrl.recursion_maps.size() == 0) ;

  fatal(tvars.size() != 1 ||
        control_set.recursion_maps.size() != 1 ||
        control_set.target_maps.size() != 1 ||
        fctrl.use_constraints) ;
    
  fcontrol::mapping_info &rmap = fctrl.recursion_maps[0] ;
  fcontrol::mapping_info &tmap = fctrl.target_maps[0] ;
  vector<multiMap> read_maps,read_map_inv ;

  for(int j=0;j<rmap.mapvec.size();++j) {
    read_maps.push_back(rmap.mapvec[j]->get_map()) ;
  }
  for(int j=0;j<read_maps.size();++j) {
    read_map_inv.push_back(multiMap()) ;
  }
  variable rvar = *(tvars.begin()) ;
  entitySet sdelta = facts.variable_existence(rvar) ;
  entitySet domain = fctrl.nr_sources + sdelta ;

  for(int j=read_maps.size()-1;j>=0;--j) {
    entitySet newdomain = facts.preimage(rmap.mapvar[j],domain).first ;
#ifdef VERBOSE
    cout << "j = " << j << ", domain = " << domain << ", newdomain = "
         << newdomain << endl ;
#endif
    inverseMap(read_map_inv[j],read_maps[j],domain,newdomain) ;
    domain = newdomain ;
  }
  
  for(int j=rmap.mapvec.size()-1;j>=0;--j)
    sdelta = rmap.mapvec[j]->preimage(sdelta).first ;
  sdelta &= fctrl.nr_sources ;
  entitySet tdelta = sdelta ;
  for(int j=0;j<tmap.mapvec.size();++j)
    tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
  cout << "sdelta_init = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif
  fastseq += sequence(sdelta) ;
  entitySet generated = tdelta ;
  const entitySet nr_sources = fctrl.nr_sources ;
  store<bool> exists ;
  exists.allocate(nr_sources) ;
  for(entitySet::const_iterator
        ei=nr_sources.begin();ei!=nr_sources.end();++ei)
    exists[*ei] = false ;
  for(entitySet::const_iterator
        ei=tdelta.begin();ei!=tdelta.end();++ei)
    exists[*ei] = true ;
  vector<const int *> mi(read_maps.size()) ;
  vector<const int *> me(read_maps.size()) ;

  do {
    for(int j=read_map_inv.size()-1;j>=0;--j) {
      entitySet candidates ;
      const int *mi ;
      for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
#ifdef DEBUG
        if(!read_map_inv[j].domain().inSet(*di)) {
          cerr << " di = " << *di << ", j = " << j << endl ;
        } else
#endif
          //        fatal(!read_map_inv[j].domain().inSet(*di)) ;
        for(mi=read_map_inv[j].begin(*di); mi!=read_map_inv[j].end(*di);++mi) 
          candidates += *mi ;
      }
      sdelta = candidates ;
    }
    // we now have candidate sdeltas, check them to see if they
    // are satisfied.
    entitySet satisfied ;
    for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
      int c = *di ;
      mi[0] = read_maps[0].begin(c) ;
      me[0] = read_maps[0].end(c) ;
      int j=0 ;
      const int je=read_maps.size();
      bool chk = true ;
      // Check that all paths from candidate go to generated cells
      while(chk && j>=0) {
        c = *mi[j] ;
        j++ ;
        if(j==je) {
          chk = chk && exists[c] && nr_sources.inSet(c) ;
          do {
            j-- ;
          } while(j>=0 && (++mi[j]) == me[j]) ;
        } else {
          mi[j] = read_maps[j].begin(c) ;
          me[j] = read_maps[j].end(c) ;
        }
        
      }
      if(chk)
        satisfied += *di ;
    }
    sdelta = satisfied ;
    entitySet tdelta = sdelta ;
    for(int j=0;j<tmap.mapvec.size();++j)
      tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
    cout << "sdelta = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif
    fastseq += sequence(sdelta) ;
    for(entitySet::const_iterator
          ei=tdelta.begin();ei!=tdelta.end();++ei) 
      exists[*ei] = true ;
    
  } while(sdelta != EMPTY) ;

  for(entitySet::const_iterator
        ei=nr_sources.begin();ei!=nr_sources.end();++ei)
    if(exists[*ei]) {
      const int start = *ei ;
      const int end = nr_sources.Max() ;
      int finish = start ;
      for(;*ei!=end && exists[*ei];++ei)
        finish = *ei ;
      if(*ei == end && exists[end])
        finish = end ;
      generated += interval(start,finish) ;
    }
  
  
  fctrl.generated[rvar] = generated ;
#ifdef VERBOSE
  cout << "recursive rule " << impl << " generating " << generated << endl ;
#endif
    
  for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
      mi!=fctrl.generated.end();++mi)
    facts.set_existential_info(mi->first,impl,mi->second) ;
    
}
  
void impl_recurse_calculator::process_var_requests(fact_db &facts) {
  process_rule_requests(impl,facts) ;
}

executeP impl_recurse_calculator::create_execution_schedule(fact_db &facts) {
  return new execute_rule(impl,fastseq,facts) ;
}

class recurse_calculator : public rule_calculator {
  ruleSet recurse_rules ;
  struct fcontrol {
    list<entitySet> control_list ;
    map<variable,entitySet> generated ;
    bool use_constraints ;
    entitySet nr_sources, constraints ;
    struct mapping_info {
      entitySet total ;
      vector<MapRepP> mapvec ;
      variable v ;
    } ;
    vector<mapping_info> recursion_maps, target_maps ;
  } ;
  map<rule,fcontrol > control_set ;
public:
  recurse_calculator(decompose_graph *gr, ruleSet rs)
  { graph_ref = gr, recurse_rules = rs ; }
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

void recurse_calculator::set_var_existence(fact_db &facts) {
  control_set.clear() ;
  
  ruleSet::const_iterator fi ;
  variableSet::const_iterator vi ;
  variableSet tvars ;
  const ruleSet &fset = recurse_rules ;
  for(fi=fset.begin();fi!=fset.end();++fi) 
    tvars += fi->targets() ;
  for(fi=fset.begin();fi!=fset.end();++fi) {
    fcontrol &fctrl = control_set[*fi] ;
    const rule_impl::info &finfo = fi->get_info().desc ;
    warn(fi->type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
      if((si->var & tvars) == EMPTY)
        sources &= vmap_source_exist(*si,facts) ;
      else {
        tvars += si->var ;
        int num_maps = si->mapping.size() ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;
        
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!facts.is_a_Map(*(miv[j]))) ;
              minfo.mapvec.push_back(
                                     MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
            }
            minfo.v = *vi ;
              
            for(i=0;(i!=num_maps)&&(++(miv[i]) == si->mapping[i].end());++i)
              miv[i] = si->mapping[i].begin() ;

            fctrl.recursion_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
    } 
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      int num_maps = si->mapping.size() ;
      vector<variableSet::const_iterator> miv(num_maps) ;
      for(int i=0;i<num_maps;++i)
        miv[i] = si->mapping[i].begin() ;

      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        int i = 0 ;
        do {
          fcontrol::mapping_info minfo ;
          for(int j=0;j!=num_maps;++j) {
            FATAL(!facts.is_a_Map(*(miv[j]))) ;
            minfo.mapvec.push_back(
                                   MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
          }
          minfo.v = *vi ;
          
          for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
            miv[i] = si->mapping[i].begin() ;
          }
          fctrl.target_maps.push_back(minfo) ;
        } while(i!=num_maps) ;
      }
    }
    fctrl.use_constraints = false ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
      warn((si->var & tvars) != EMPTY) ;
      constraints &= vmap_source_exist(*si,facts) ;
      fctrl.use_constraints = true ;
    }
    fctrl.nr_sources = sources ;
    fctrl.constraints = constraints ;
    warn(fctrl.recursion_maps.size() == 0) ;
  }
  

  map<variable,entitySet> tvar_computed,tvar_update ;
  for(vi=tvars.begin();vi!=tvars.end();++vi) {
    tvar_computed[*vi] = facts.variable_existence(*vi) ;
    tvar_update[*vi] = tvar_computed[*vi] ;
  }
  map<variable,entitySet>::iterator mi ;
  for(mi=tvar_update.begin();mi!=tvar_update.end();++mi)
    mi->second = EMPTY ;

  bool finished = false;

  while(!finished) {
    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      // partial vmap_souce_exist
      entitySet srcs = fctrl.nr_sources ;
      for(int i=0;i<fctrl.recursion_maps.size();++i) {
        entitySet sdelta = tvar_computed[fctrl.recursion_maps[i].v] ;
        for(int j=fctrl.recursion_maps[i].mapvec.size()-1;j>=0;--j)
          sdelta = fctrl.recursion_maps[i].mapvec[j]->preimage(sdelta).first ;
        fctrl.recursion_maps[i].total += sdelta ;
        srcs &= fctrl.recursion_maps[i].total ;
      }
      if(fctrl.use_constraints) {
        if((srcs & fctrl.constraints) != fctrl.constraints) {
          cerr << "recursive rule: " << *fi
               << " cannot supply all entitites in constraint" << endl ;
          cerr << "constraints = " << fctrl.constraints << endl ;
          entitySet sac = srcs & fctrl.constraints ;
          cerr << "srcs & constraints = " << sac << endl ;
          exit(-1) ;
        }
        srcs &= fctrl.constraints ;
      }

      fctrl.control_list.push_back(srcs) ;

      for(int i=0;i<fctrl.target_maps.size();++i) {
        entitySet trgts = srcs ;
        for(int j=0;j<fctrl.target_maps[i].mapvec.size();++j)
          trgts = fctrl.target_maps[i].mapvec[j]->image(trgts) ;
        tvar_update[fctrl.target_maps[i].v] += trgts ;
        fctrl.generated[fctrl.target_maps[i].v] += trgts ;
      }
#ifdef VERBOSE
      cout << "recursive rule " << *fi << " generating " << srcs << endl ;
#endif
    }
    finished = true ;
    for(vi=tvars.begin();vi!=tvars.end();++vi) {
      entitySet tvar_delta = tvar_update[*vi] - tvar_computed[*vi] ;
      if(tvar_delta!=EMPTY)
        finished = false ;
      tvar_computed[*vi] = tvar_update[*vi] ;
    }
  }

  for(fi=fset.begin();fi!=fset.end();++fi) {
    fcontrol &fctrl = control_set[*fi] ;
    for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
        mi!=fctrl.generated.end();++mi)
      facts.set_existential_info(mi->first,*fi,mi->second) ;
  }
}

void recurse_calculator::process_var_requests(fact_db &facts) {
  ruleSet::const_iterator fi ;
  for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
    fcontrol &fctrl = control_set[*fi] ;
    entitySet control = process_rule_requests(*fi,facts) ;
    list<entitySet>::iterator ci ;
    entitySet total ;
    for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
      *ci -= total ;
      total += *ci ;
    }
    do {
      if(fctrl.control_list.back() == EMPTY)
        fctrl.control_list.pop_back() ;
      if(!fctrl.control_list.empty())
        fctrl.control_list.back() &= control ;
    } while(!fctrl.control_list.empty() && fctrl.control_list.back() == EMPTY) ;
  }
}

executeP recurse_calculator::create_execution_schedule(fact_db &facts) {
  CPTR<execute_list> el = new execute_list ;
    
  map<rule, list<entitySet>::const_iterator> rpos ;
  ruleSet::const_iterator ri ;
    
  for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
    rpos[*ri] = control_set[*ri].control_list.begin() ;

  bool finished ;
  do {
    finished = true ;
      
    for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
      const fcontrol &fctrl = control_set[*ri] ;
      list<entitySet>::const_iterator &li = rpos[*ri] ;
      if(li != fctrl.control_list.end()) {
        finished = false ;
        if(li->size() != 0)
          el->append_list(new execute_rule(*ri,sequence(*li),facts)) ;
        li++ ;
      }
    }
  } while(!finished) ;
    
  if(el->size() == 0)
    return 0 ;
  else
    return executeP(el) ;
}

class dag_calculator : public rule_calculator {
  digraph dag ;
  vector<rule> rule_schedule ;
public:
  dag_calculator(decompose_graph *gr, digraph gin) ;
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

dag_calculator::dag_calculator(decompose_graph *gr, digraph gin) {
  graph_ref = gr ;
  dag = gin ;
  vector<digraph::nodeSet> dag_sched = schedule_dag(dag) ;
  extract_rule_sequence(rule_schedule,dag_sched) ;
#ifdef DEBUG
  // sanity check, all nodes should be scheduled
  digraph::nodeSet allnodes ;
  for(int i=0;i< dag_sched.size();++i) 
    allnodes += dag_sched[i] ;
  warn(allnodes != dag.get_all_nodes()) ;
#endif
}

  
void dag_calculator::set_var_existence(fact_db &facts) {
  for(int i=0;i<rule_schedule.size();++i) 
    calc(rule_schedule[i])->set_var_existence(facts) ;
}

void dag_calculator::process_var_requests(fact_db &facts) {
  vector<rule>::reverse_iterator ri ;
  for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
    calc(*ri)->process_var_requests(facts) ;
}

executeP dag_calculator::create_execution_schedule(fact_db &facts) {
  CPTR<execute_list> elp = new execute_list ;
  for(int i=0;i<rule_schedule.size();++i) 
    elp->append_list(calc(rule_schedule[i])->create_execution_schedule(facts)) ;
  return executeP(elp) ;
}

loop_calculator::loop_calculator(decompose_graph *gr, digraph gin) {
  graph_ref = gr ;
  dag = gin ;
  ruleSet loopset =
    extract_rules(gin.get_all_nodes() & gr->looping_rules) ;
  if(loopset.size() != 1) {
    cerr << "internal consistency error, loopset size != 1" << endl ;
    exit(-1) ;
  }

  tlevel = loopset.begin()->source_time() ;
  dag.remove_node((*loopset.begin()).ident()) ;
  digraph dagt = dag.transpose() ;
  ruleSet collapse_rules, all_rules ;
  variableSet collapse_vars ;
  
  all_rules = extract_rules(dag.get_all_nodes()) ;
  for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) 
    if(ri->target_time().before(ri->source_time()))
      collapse_rules += *ri ;
  
  if(collapse_rules.size() != 1 ||
     collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
    cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
         << endl << "error is not recoverable" << endl ;
    exit(-1) ;
  }

  variableSet all_vars = extract_vars(dag.get_all_nodes()) ;
  for(variableSet::const_iterator vi=all_vars.begin();vi!=all_vars.end();++vi)
    if(vi->get_info().offset == 1)
      advance_vars += *vi ;
  
  collapse_vars = collapse_rules.begin()->targets() ;
  cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

  // create Output Variable
  output = variable(variable("OUTPUT"),tlevel) ;
  variableSet outputSet ;
  outputSet = interval(output.ident(),output.ident()) ;

  // Schedule part of graph that leads to collapse
  vector<digraph::nodeSet> dag_sched
    = schedule_dag(dag, EMPTY,visit_nodes(dagt,collapse_vars)) ;
  extract_rule_sequence(rule_schedule,dag_sched) ;
  extract_rule_sequence(collapse,dag_sched) ;

  // Schedule advance part of loop.  First try to schedule any output, then
  // schedule the advance
  digraph::nodeSet visited ;
  for(int i = 0;i<dag_sched.size();++i)
    visited += dag_sched[i] ;
  dag_sched = schedule_dag(dag,visited, visit_nodes(dagt,outputSet)) ;
  if(dag_sched.size() == 0)
    output_present = false ;
  else 
    output_present = true ;
  
  extract_rule_sequence(rule_schedule,dag_sched) ;
  extract_rule_sequence(advance,dag_sched) ;
  for(int i = 0;i<dag_sched.size();++i)
    visited += dag_sched[i] ;
  // now schedule everything that hasn't been scheduled
  dag_sched = schedule_dag(dag,visited) ;
  extract_rule_sequence(rule_schedule,dag_sched) ;
  extract_rule_sequence(advance,dag_sched) ;
  
#ifdef DEBUG
  // sanity check, all nodes should be scheduled
  digraph::nodeSet allnodes = visited ;
  for(int i=0;i< dag_sched.size();++i) 
    allnodes += dag_sched[i] ;
  warn(allnodes != dag.get_all_nodes()) ;
  if(allnodes != dag.get_all_nodes()) {
    cerr << " rules NOT scheduled = " << endl
         << extract_rules(dag.get_all_nodes() - allnodes) << endl ;
    cerr << " variables NOT scheduled = "
         << extract_vars(dag.get_all_nodes() - allnodes) << endl ;
    vector<digraph::nodeSet> components =
      component_sort(dag).get_components() ;
    for(int i=0;i<components.size();++i) {
      if(components[i].size() > 1) {
        cerr << "reason: graph not a dag, strongly connected component found"
             << endl ;
        cerr << "component rules = " << endl ;
        cerr << extract_rules(components[i]) << endl ;
      }
    }
    exit(-1) ;
  }
#endif
}

void loop_calculator::set_var_existence(fact_db &facts) {
  for(int i=0;i<rule_schedule.size();++i) 
    calc(rule_schedule[i])->set_var_existence(facts) ;
}

void loop_calculator::process_var_requests(fact_db &facts) {
  variableSet var_requests = advance_vars ;
  variableSet::const_iterator vi ;
  if(output_present)
    var_requests += output ;
  for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
    entitySet vexist = facts.variable_existence(*vi) ;
    facts.variable_request(*vi,vexist) ;
  }
  vector<rule>::reverse_iterator ri ;
  for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
    calc(*ri)->process_var_requests(facts) ;

}

executeP loop_calculator::create_execution_schedule(fact_db &facts) {
  CPTR<execute_list> col = new execute_list ;
  for(int i=0;i<collapse.size();++i) 
    col->append_list(calc(collapse[i])->create_execution_schedule(facts)) ;
  CPTR<execute_list> adv = new execute_list ;
  for(int i=0;i<advance.size();++i) 
    adv->append_list(calc(advance[i])->create_execution_schedule(facts)) ;
  return new execute_loop(cond_var,executeP(col),executeP(adv),tlevel) ;
}

class conditional_calculator : public rule_calculator {
  digraph dag ;
  vector<rule> rule_schedule ;
  variable cond_var ;
public:
  conditional_calculator(decompose_graph *gr, digraph gin,
                         variable conditional) ;
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;

conditional_calculator::conditional_calculator(decompose_graph *gr,
                                               digraph gin,
                                               variable conditional) {
  graph_ref = gr ;
  dag = gin ;
  cond_var = conditional ;
  vector<digraph::nodeSet> dag_sched = schedule_dag(dag) ;
  extract_rule_sequence(rule_schedule,dag_sched) ;
#ifdef DEBUGG
  // sanity check, all nodes should be scheduled
  digraph::nodeSet allnodes ;
  for(int i=0;i< dag_sched.size();++i) 
    allnodes += dag_sched[i] ;
  warn(allnodes != dag.get_all_nodes()) ;
#endif
}

    
void conditional_calculator::set_var_existence(fact_db &facts) {
  for(int i=0;i<rule_schedule.size();++i) 
    calc(rule_schedule[i])->set_var_existence(facts) ;
}

void conditional_calculator::process_var_requests(fact_db &facts) {
  vector<rule>::reverse_iterator ri ;
  for(ri=rule_schedule.rbegin();ri!=rule_schedule.rend();++ri)
    calc(*ri)->process_var_requests(facts) ;

}

executeP conditional_calculator::create_execution_schedule(fact_db &facts) {
  CPTR<execute_list> elp = new execute_list ;
  for(int i=0;i<rule_schedule.size();++i) 
    elp->append_list(calc(rule_schedule[i])->create_execution_schedule(facts)) ;
  return new execute_conditional(executeP(elp),cond_var) ;
}


class promote_calculator : public rule_calculator {
  rule r ;
public:
  promote_calculator(decompose_graph *gr, rule rin)
  { graph_ref = gr, r = rin; }
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;
  
void promote_calculator::set_var_existence(fact_db &facts) {
}

void promote_calculator::process_var_requests(fact_db &facts) {
}

executeP promote_calculator::create_execution_schedule(fact_db &facts) {
  return executeP(0) ;
}

class generalize_calculator : public rule_calculator {
  rule r ;
public:
  generalize_calculator(decompose_graph *gr, rule rin)
  { graph_ref = gr, r = rin; }
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;
  
void generalize_calculator::set_var_existence(fact_db &facts) {
}

void generalize_calculator::process_var_requests(fact_db &facts) {
}

executeP generalize_calculator::create_execution_schedule(fact_db &facts) {
  return executeP(0) ;
}

class priority_calculator : public rule_calculator {
  rule r ;
public:
  priority_calculator(decompose_graph *gr, rule rin)
  { graph_ref = gr, r = rin; }
  virtual void set_var_existence(fact_db &facts) ;
  virtual void process_var_requests(fact_db &facts) ;
  virtual executeP create_execution_schedule(fact_db &facts) ;
} ;
  
void priority_calculator::set_var_existence(fact_db &facts) {
}

void priority_calculator::process_var_requests(fact_db &facts) {
}

executeP priority_calculator::create_execution_schedule(fact_db &facts) {
  return executeP(0) ;
}

class error_calculator : public rule_calculator {
public:
  error_calculator() {}
  virtual void set_var_existence(fact_db &facts)
  { cerr << "Internal consistency error" << endl ; exit(-1);}
  virtual void process_var_requests(fact_db &facts) 
  { cerr << "Internal consistency error" << endl ; exit(-1);}
  virtual executeP create_execution_schedule(fact_db &facts)
  { cerr << "Internal consistency error" << endl ; exit(-1);
    return executeP(0);}
} ;

decompose_graph::decompose_graph(digraph dg,
                                 digraph::nodeSet sources,
                                 digraph::nodeSet targets) {


  // Add nodes to the graph so that the source variables and target variables
  // aren't left dangling.
  variableSet vars = extract_vars(dg.get_all_nodes()) ;
  start = variable("__START__") ;
  finish = variable("__FINISH__") ;
  dg.add_edges(start.ident(),sources) ;
  dg.add_edges(targets,finish.ident()) ;

  // In the first step we decompose the graph by sorting the graph
  // into iteraton levels.
  map<time_ident,digraph::nodeSet> time_sort_nodes ;

  for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
    time_sort_nodes[vi->time()] += (*vi).ident() ;

  ruleSet rules = extract_rules(dg.get_all_nodes()) ;

  for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
    // sort rules based on time level
    const time_ident rule_time =
      ri->type()==rule::COLLAPSE?ri->source_time():ri->target_time() ;
    time_sort_nodes[rule_time] += (*ri).ident() ;
    // collect information about unit rules
    if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
      if(ri->targets().size() != 1) {
        cerr << "unit rule must have only one variable as target: "
             << endl << *ri << endl ;
        exit(-1) ;
      }
      variable unit_var = *(ri->targets().begin()) ;
      if(reduce_vars.inSet(unit_var)) {
        cerr << "only one unit rule may be defined for reduction variable "
             << unit_var << endl ;
        exit(-1) ;
      }
      reduce_vars += unit_var ;
      reduce_info &rinfo = reduce_map[unit_var] ;
      rinfo.unit_rule = *ri ;
      dg.add_edge(unit_var.ident(),(*ri).ident()) ;
      rinfo.apply_rules = extract_rules(dg.transpose()[unit_var.ident()]) ;
      rinfo.apply_rules -= rinfo.unit_rule ;
      for(ruleSet::const_iterator
            rii=rinfo.apply_rules.begin();rii!=rinfo.apply_rules.end();++rii) {
        if(rii->get_info().rule_impl->get_rule_class() != rule_impl::APPLY) {
          cerr << "pointwise rule is producing reduction variable: "
               << unit_var << endl ;
          cerr << "This rule should be an apply rule. (offending rule is):"
               << endl << *rii << endl ;
          exit(-1) ;
        }
        if(rii->targets().size() != 1) {
          cerr << "apply rule should apply to only one reduction variable"
               << endl << "offending rule is:" << endl << *rii << endl ;
          exit(-1) ;
        }
        dg.add_edge(unit_var.ident(),(*rii).ident()) ;
      }
    }
    // collect information about looping rules
    if(ri->get_info().qualifier() == "looping")
      looping_rules += *ri ;
    // collect information about rules that are conditionally executed.
    if(ri->get_info().desc.conditionals != EMPTY) {
      variableSet conds = ri->get_info().desc.conditionals ;
      if(conds.size() != 1) {
        cerr << "improper rule: " << *ri << endl
             << "Rule Improperly specifies more than one conditional" << endl ;
        cerr << "Error not recoverable." << endl ;
        exit(-1) ;
      }
      variable cond = *(conds.begin()) ;

      if(ri->type() != rule::COLLAPSE &&
         (ri->targets().size() != 1 ||
          (ri->targets().begin())->get_info().name != "OUTPUT")) {
        cerr << "conditional rules must either be collapse rules or " << endl
             << "have a OUTPUT as the sole argument" << endl ;
        cerr << "offending rule: " <<*ri << endl ; 
        exit(-1) ;
      }

      conditional_rules += *ri ;
      conditional_vars += cond ;
      conditional_map[cond] += *ri ;
    }
  }

  // Create apply rule to unit rule mapping
  for(variableSet::const_iterator
        vi = reduce_vars.begin(); vi != reduce_vars.end(); ++vi) {

    reduce_info &rinfo = reduce_map[*vi] ;

    for(ruleSet::const_iterator
          it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) {
      apply_to_unit[*it] = rinfo.unit_rule ;
    }
  }

#ifdef VERBOSE
  cout << "reduction variables = " << reduce_vars << endl;
  cout << "looping rules = " << endl << looping_rules << endl ;
  cout << "conditional rules = " << endl << conditional_rules << endl ;
#endif
  
  // compute a depth first ordering of time levels
  vector<time_ident> time_stack,time_sequence ;
  time_stack.push_back(time_ident()) ;
  for(;!time_stack.empty();) {
    time_ident ctime = time_stack.back() ;
    time_stack.pop_back() ;
    time_sequence.push_back(ctime) ;
    const vector<time_ident> &children = ctime.children()  ;
    for(vector<time_ident>::const_iterator
          i=children.begin();i!=children.end();++i)
      time_stack.push_back(*i) ;
  }
  // proceed from the bottom of the time hierarchy tree, extracting each
  // time level and constructing a supernode based on that time level.
  digraph main_graph = dg ;
  for(vector<time_ident>::reverse_iterator
        ti = time_sequence.rbegin();ti!=time_sequence.rend();++ti) 
    if(time_sort_nodes[*ti] != EMPTY) {
      rule snode = create_supernode(main_graph,time_sort_nodes[*ti]) ;
      if(*ti != time_ident()) {
        const time_ident rule_time =snode.target_time() ;
        time_sort_nodes[rule_time] += snode.ident() ;
      } else
        top_level_rule = snode ; 

      if((snode.sources() & snode.targets()) != EMPTY) {
        cerr << "warning, time level " << snode << " is recursive." << endl ;
        cerr << "build and collapse rules not allowed to form a "
             << "recursive block" << endl ;
        exit(-1) ;
      }
   }

  // Go over each time level and decompose it into strongly connected components
  // At this level, loops remain as a single component.
  ruleSet time_levels = supernodes ;
  for(ruleSet::const_iterator
        ri=time_levels.begin();ri!=time_levels.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    if(*ri != top_level_rule) {
      ruleSet loop_rules = ruleSet(sni.graph.get_all_nodes() & looping_rules) ;
      
      if(loop_rules.size() != 1) {
        time_ident ti ;
        variableSet vars = extract_vars(sni.graph.get_all_nodes()) ;
        for(variableSet::const_iterator vi= vars.begin();vi!=vars.end();++vi) 
          if(ti.before(vi->time()))
            ti = vi->time() ;
        
        cerr << "incorrect specification of iteration on time level " << ti
             << endl << "no loop specification was produced." << endl
             << "Perhaps no advance rules could be applied at this level?"
             << endl ;
        exit(-1) ;
      }
      // Insure that collapse rules are included in loop supernode
      sni.graph.add_edges(sni.targets,(*(loop_rules.begin())).ident()) ;
    }
    vector<digraph::nodeSet> components =
      component_sort(sni.graph).get_components() ;
    vector<digraph::nodeSet>::const_iterator ci ;
    for(ci=components.begin();ci!=components.end();++ci)
      if(ci->size() > 1) {
        variableSet vars = variableSet(*ci & reduce_vars) ;
        ruleSet loop_rules = ruleSet(sni.graph.get_all_nodes() &
                                     looping_rules) ;
        if(vars== EMPTY || loop_rules != EMPTY)
          create_supernode(sni.graph,*ci) ;
        else {
          if(vars.size() > 1) {
            cerr << "reductions can't mix with recursions." << endl ;
            cerr << "error occured when parsing variables " << vars << endl ;
            exit(-1) ;
          }
          // remove the self-referential loops and make all pointwise rules
          // depend on the unit rule.
          variable rvar = *(vars.begin()) ;
          reduce_info &rinfo = reduce_map[rvar] ;

          for(ruleSet::const_iterator
                it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
            sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
          sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
          sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
        }
      }
  }

  ruleSet new_rules = ruleSet(supernodes - time_levels) ;
  // Now go over remaining supernodes to decompose time loops
  for(ruleSet::const_iterator
        ri=new_rules.begin();ri!=new_rules.end();++ri) {
    super_node_info &sni = supermap[*ri] ;

    ruleSet loop_rules = ruleSet(sni.graph.get_all_nodes() & looping_rules) ;

    digraph gtemp = sni.graph ;

    warn(loop_rules.size() > 1) ;

    if(loop_rules.size() == 1)
      gtemp.remove_node((*loop_rules.begin()).ident()) ;
    else continue ;
    
    vector<digraph::nodeSet> components =
      component_sort(gtemp).get_components() ;
    vector<digraph::nodeSet>::const_iterator ci ;
    for(ci=components.begin();ci!=components.end();++ci)
      if(ci->size() > 1) {
        variableSet vars = variableSet(*ci & reduce_vars) ;
        if(vars== EMPTY)
          create_supernode(sni.graph,*ci) ;
        else {
          if(vars.size() > 1) {
            cerr << "reductions can't mix with recursions." << endl ;
            cerr << "error occured when parsing variables " << vars << endl ;
            exit(-1) ;
          }
          // remove the self-referential loops and make all pointwise rules
          // depend on the unit rule.
          variable rvar = *(vars.begin()) ;
          reduce_info &rinfo = reduce_map[rvar] ;

          for(ruleSet::const_iterator
                it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
            sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
          sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
          sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
        }
      }
  }
  // Create conditional execution super nodes.
  // A conditional rule may exclusively depend on other rules for execution
  // In this case, we only need to execute these other rules when the condition
  // is satisified.  This is accomplished by creating a conditional supernode
  // that encapsulates this grouped behavior.
  new_rules = supernodes ;
  ruleSet cond_supernodes  ;
  for(ruleSet::const_iterator
        ri = new_rules.begin();ri!=new_rules.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    digraph::nodeSet graph_nodes = sni.graph.get_all_nodes() &
      ~sni.sources & ~sni.targets ;
    
    variableSet cond_vars = variableSet(graph_nodes & conditional_vars) ;
    // If there are conditional variables in this graph, then decompose them
    if(cond_vars != EMPTY) {
      // loop over conditional variables and extract the functions that
      // are contingent on them
      variableSet::const_iterator vi ;
      for(vi = cond_vars.begin();vi!=cond_vars.end();++vi) {

        ruleSet cond_rules = conditional_map[*vi] ;
        // Note, all of the rules that are conditional of this variable
        // should be in this graph, if not report this error and exit.
        if((graph_nodes & cond_rules) != cond_rules) {
          cerr << "problem with conditional functions of variable " << *vi
               << endl ;
          cerr << "rule(s) not in same level of graph:" << endl ;
          ruleSet except = ruleSet(cond_rules - (graph_nodes & cond_rules)) ;
          cerr << except << endl ;
          cerr << "error not recoverable." << endl ;
          exit(-1) ;
        }
        // Find the set of nodes that are exclusive in generating the
        // rules conditional on this variable (*vi)
        digraph::nodeSet cond_nodes =
          visit_nodes_exclusive(sni.graph.transpose(), cond_rules) ;
        // remove the rules that are exclusive to the conditional evaluation
        // itself
        cond_nodes = cond_nodes -
          visit_nodes_exclusive(sni.graph.transpose(),
                                interval((*vi).ident(),(*vi).ident())) ;

        // create a super node for this conditional set of rules
        cond_supernodes += create_supernode(sni.graph,cond_nodes,*vi) ;
      }
    }
  }
  // collect info on supernodes
  // remove recursive dependencies from rules
  ruleSet all_rules = supernodes ;
  for(ruleSet::const_iterator
        ri = supernodes.begin();ri!=supernodes.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    digraph::nodeSet graph_nodes = sni.graph.get_all_nodes() &
      ~sni.sources & ~sni.targets ;
    all_rules += extract_rules(graph_nodes) ;
    ruleSet loop_rules = ruleSet(sni.graph.get_all_nodes() & looping_rules) ;
    if(loop_rules != EMPTY)
      looping_supernodes += *ri ;
    if(loop_rules == EMPTY && (sni.sources & sni.targets) != EMPTY)
      recursive_supernodes += *ri ;
    else {
      // fix graph so that it is a directed acyclic graph.
      // At this stage the graph should have all cycles except
      // cycles of length two removed from the graph.  We now
      // remove these cycles and add the appropriate dependencies
      // to ensure that these cycles are scheduled last.
      digraph g = sni.graph ;
      digraph gt = sni.graph.transpose() ;
      ruleSet grules = extract_rules(graph_nodes - looping_rules) ;
      ruleSet cycle_rule ;
      ruleSet non_cycle_rule ;
      for(ruleSet::const_iterator
            rii = grules.begin(); rii != grules.end(); ++rii) {
        int id = (*rii).ident() ;
        if((g[id] & gt[id]) != EMPTY) {  // found a cycle
          if(!supernodes.inSet(*rii)) {
            cerr << "internal consistency error, offending rule = "
                 << *rii << endl ;
            exit(-1) ;
          }
          variableSet cvars = extract_vars(g[id] & gt[id]) ;
          cycle_rule += *rii ;
          // get other rules that generate cycle variables
          ruleSet other_rules ;
          for(variableSet::const_iterator
                vi = cvars.begin();vi!=cvars.end();++vi)
            other_rules += extract_rules(gt[(*vi).ident()]) ;
          other_rules -= *rii ;
          // remove edges causing loop
          for(variableSet::const_iterator
                vi = cvars.begin();vi!=cvars.end();++vi)
            sni.graph.remove_edge((*vi).ident(),id) ;
          // replace these edges with rule dependencies
          sni.graph.add_edges(other_rules,id) ;
        }
      }
      if((cycle_rule & non_cycle_rule) != EMPTY) {
        cerr << "internal consistency error, rule cycles collide for rules:"
             << endl ;
        cerr << extract_rules(cycle_rule & non_cycle_rule) ;
        exit(-1) ;
      }

    }
  }

  // Create schedule calculators for all graph nodes
  ruleSet working_set = ruleSet(all_rules - supernodes) ;
  for(ruleSet::const_iterator
        ri = working_set.begin();ri != working_set.end(); ++ri) {
    if(ri->type() != rule::INTERNAL) {
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
        rule_process[*ri] = new apply_calculator(this,*ri,apply_to_unit[*ri]) ;
      else
        rule_process[*ri] = new impl_calculator(this,*ri) ;
    } else if(ri->get_info().qualifier() == "promote")
      rule_process[*ri] = new promote_calculator(this,*ri) ;
    else if(ri->get_info().qualifier() == "generalize")
      rule_process[*ri] = new generalize_calculator(this,*ri) ;
    else if(ri->get_info().qualifier() == "priority")
      rule_process[*ri] = new priority_calculator(this,*ri) ;
    else
      rule_process[*ri] = new error_calculator ;
  }
  ruleSet reduction_corrections ;
  // now do same for supernodes
  for(ruleSet::const_iterator
        ri = recursive_supernodes.begin();ri!=recursive_supernodes.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    if((sni.targets & reduce_vars) != EMPTY) {
      // this is really a reduction block, so remove recursion dependencies
      // and save to schedule as a dag.
      warn((sni.targets & reduce_vars) != sni.targets) ;
      // remove the self-referential loops and make all pointwise rules
      // depend on the unit rule.
      warn(sni.targets.size() != 1) ;
      variable rvar = *(sni.targets.begin()) ;
      reduce_info &rinfo = reduce_map[rvar] ;
      
      for(ruleSet::const_iterator
            it=rinfo.apply_rules.begin();it!=rinfo.apply_rules.end();++it) 
        sni.graph.remove_edge(rvar.ident(),(*it).ident()) ;
      sni.graph.add_edges(rinfo.unit_rule.ident(),rinfo.apply_rules) ;
      sni.graph.remove_edge(rvar.ident(),rinfo.unit_rule.ident()) ;
      reduction_corrections += *ri ;
      
    } else {
      // recursive rule block
      ruleSet recurse_rules = extract_rules(sni.graph.get_all_nodes()) ;
      if(recurse_rules.size() == 1 &&
         recurse_rules.begin()->get_info().desc.constraints.size() == 0) {
        rule_process[*ri] =
          new impl_recurse_calculator(this,*(recurse_rules.begin())) ;
      } else {
        rule_process[*ri] = new recurse_calculator(this,recurse_rules) ;
      }
    }
  }
  recursive_supernodes -= reduction_corrections ;
  
  for(ruleSet::const_iterator
        ri = looping_supernodes.begin();ri!=looping_supernodes.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    rule_process[*ri] = new loop_calculator(this,sni.graph) ;
  }
  for(ruleSet::const_iterator
        ri = cond_supernodes.begin();ri!=cond_supernodes.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    variable cond_var = *(ri->get_info().desc.conditionals.begin()) ;
    rule_process[*ri] = new conditional_calculator(this,sni.graph,cond_var) ;
  }

  ruleSet dag_rules = ruleSet(supernodes -
                              (recursive_supernodes +
                               looping_supernodes +
                               cond_supernodes)) ;

  for(ruleSet ::const_iterator
        ri = dag_rules.begin();ri!=dag_rules.end();++ri) {
    super_node_info &sni = supermap[*ri] ;
    rule_process[*ri] = new dag_calculator(this,sni.graph) ;
  }

}

rule decompose_graph::create_supernode(digraph &g,
                                       const digraph::nodeSet nodes,
                                       variable cond_var) {
  // extract rules in this supernode 
  ruleSet srules = extract_rules(nodes) ;
  warn(srules == EMPTY) ;

  // include all variables that are sources or sinks to any rule in this
  // supernode in the supernode graph.
  digraph::nodeSet ns = nodes ;
  variableSet sources,targets,locals,local_sources,local_targets ;

  // Calculate all sources and targets for rules in this supernode
  digraph gt = g.transpose() ;
  for(ruleSet::const_iterator ii =srules.begin();ii!=srules.end();++ii) {
    int id = (*ii).ident() ;
    sources += extract_vars(gt[id]) ;
    targets += extract_vars(g[id]) ;
  }

  // Add these new nodes to the supernode graph
  ns += sources + targets ;

  variableSet all_vars = extract_vars(ns) ;
  // find all variables that are referenced exclusively from within
  // the supernode, that is find all local variables
  digraph::nodeSet not_ns = ~ns ;
  for(variableSet::const_iterator
        vi=all_vars.begin();vi!=all_vars.end();++vi) {
    int id = (*vi).ident() ;
    if((gt[id] & not_ns) == EMPTY)
      local_sources += *vi ;
    if((g[id] & not_ns) == EMPTY)
      local_targets += *vi ;
  }

  locals = local_sources & local_targets ;
  // remove local variables from the super node rule source and target lists
  sources -= local_sources ;
  targets -= local_targets ;

  super_node_info sninfo ;
  
  // create supernode info
  sninfo.graph = g.subgraph(ns) ;
  sninfo.graph.add_edges(start.ident(),sources) ;
  sninfo.graph.add_edges(targets,finish.ident()) ;
  
  sninfo.sources = sources ;
  sninfo.targets = targets ;
  // remove components of graph that comprise supernode from parent graph
  ns -= (sources + targets) ;
  g.remove_nodes( ns ) ;
  // create a rule for the supernode
  rule node_rule = make_super_rule(sources,targets,cond_var) ;
  // Add supernode rule to parent graph
  g.add_edges(sources,node_rule.ident()) ;
  g.add_edges(node_rule.ident(),targets) ;
  
  // keep a map of supernodes
  supermap[node_rule] = sninfo ;
  supernodes += node_rule ;

#ifdef VERBOSE
  cout << "node_rule = " << node_rule << endl ;
  cout << "locals  = " << locals << endl ;
  cout << "rules in block = " << endl << srules << endl ;
  if((sources & targets) != EMPTY)
    cout << "sources & targets = " << variableSet(sources & targets) << endl ;
#endif

  return node_rule ;
}

void decompose_graph::existential_analysis(fact_db &facts) {

  (rule_process[top_level_rule])->set_var_existence(facts) ;
  variableSet var_requests = top_level_rule.targets() ;
  variableSet::const_iterator vi ;
  for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
    entitySet vexist = facts.variable_existence(*vi) ;
    facts.variable_request(*vi,vexist) ;
  }
  (rule_process[top_level_rule])->process_var_requests(facts) ;
}

executeP decompose_graph::execution_schedule(fact_db &facts) {

  CPTR<execute_list> schedule = new execute_list ;
  schedule->append_list(new allocate_all_vars) ;
  executeP top_level_schedule = (rule_process[top_level_rule])->
    create_execution_schedule(facts) ;
  if(top_level_schedule == 0) 
    return executeP(0) ;
  schedule->append_list(top_level_schedule) ;
  return executeP(schedule) ;
}
  
}

namespace Loci {
  executeP create_execution_schedule(rule_db &rdb,
                                     fact_db &facts,
                                     std::string target_string) {
  
    variableSet given = facts.get_typed_variables() ;
    variableSet target(expression::create(target_string)) ;
    
    cout << "generating dependency graph..." << endl ;
    digraph gr = dependency_graph(rdb,given,target).get_graph() ;

    // If graph is empty, return a null schedule 
    if(gr.get_target_nodes() == EMPTY)
      return executeP(0) ;

    cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr) ;
    
    cout << "decomposing graph..." << endl ;
    decompose_graph decomp(gr,given,target) ;
    
    cout << "existential analysis..." << endl ;
    decomp.existential_analysis(facts) ;
    
    cout << "creating execution schedule..." << endl;
    return decomp.execution_schedule(facts) ;
  }
}

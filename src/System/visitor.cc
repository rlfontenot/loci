#include <sys/stat.h>
#include "visitor.h"
#include "comp_tools.h"

#include <Tools/stream.h>
#include <distribute.h>
#include <parameter.h>

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <string>
using std::string ;
#include <sstream>
using std::ostringstream ;
#include <iostream>
using std::ostream ;
#include <fstream>
using std::ofstream ;
#include <sstream>
using std::stringstream ;
#include <list>
using std::list ;
#include <set>
using std::set ;

namespace Loci {

  ///////////////////////////////////////////////////////////////////
  // orderVisitor
  ///////////////////////////////////////////////////////////////////
  // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet>
  orderVisitor::order_dag(const digraph &g,
                          digraph::vertexSet start_vertices,
                          digraph::vertexSet only_vertices) {
    digraph gt = g.transpose() ;
    
    vector<digraph::vertexSet> schedule ; 
    // First schedule any vertices that have no edges leading into them
    //and have not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    if(working != EMPTY) 
      schedule.push_back(working) ;
    
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = start_vertices + working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni) 
        new_vertices += g[*ni] ;
      
      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      // We only schedule vertices that are also in the only_vertices set
      working = new_vertices & only_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all vertices leading
      // to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni) 
        if((gt[*ni] & visited_vertices) == gt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // and these new vertices to the schedule
      if(new_vertices != EMPTY)
        schedule.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
    return schedule ;
  }  

  void orderVisitor::visit(loop_compiler& lc) {
    
    lc.collapse_sched = order_dag(lc.collapse_gr) ;
    lc.advance_sched = order_dag(lc.advance_gr) ;
    /*
    lc.collapse_sched = order_dag(lc.loop_gr,EMPTY,lc.collapse_vertices) ;
    digraph::vertexSet visited ;
    for(int i=0;i!=lc.collapse_sched.size();++i)
      visited += lc.collapse_sched[i] ;
    lc.advance_sched = order_dag(lc.loop_gr,visited) ;
    */
  }

  void orderVisitor::visit(dag_compiler& dc) {
    dc.dag_sched = order_dag(dc.dag_gr) ;
  }

  void orderVisitor::visit(conditional_compiler& cc) {
    cc.dag_sched = order_dag(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // assembleVisitor
  /////////////////////////////////////////////////////////////////
  namespace {
    class check_dump_on_startup {
      bool checked_dump ;
      bool do_dump ;
    public :
      check_dump_on_startup() {
        do_dump = false ;
        checked_dump = false ;
      }
      bool ok() {
        if(!checked_dump) {
          do_dump=true ;
          struct stat statbuf ;
          if(stat("dump_vars",&statbuf))
            do_dump = false ;
          else if(!S_ISDIR(statbuf.st_mode)) 
            do_dump = false ;
          checked_dump = true ;
        } 
        return do_dump ;
      }
      
    } ;
    
    check_dump_on_startup check_dump_vars ;
    
    class execute_dump_var : public execute_modules {
      variableSet dump_vars ;
    public:
      execute_dump_var(variableSet vars) : dump_vars(vars) {}
      virtual void execute(fact_db &facts) ;
      virtual void Print(std::ostream &s) const ;
    } ;
    
    map<variable, int> dump_var_lookup ;
    
    void execute_dump_var::execute(fact_db &facts) {
      
      for(variableSet::const_iterator vi=dump_vars.begin();
          vi!=dump_vars.end();++vi) {
        variable v = *vi ;
        debugout << "dumping variable " << v << endl ;
        storeRepP st = facts.get_variable(v) ;
        storeRepP sc = st ;
        if(facts.isDistributed() && st->RepType() == STORE) {
          sc = collect_store(st,facts) ;
        }
        if(MPI_rank == 0) {
          ostringstream oss ;
          oss << "dump_vars/"<<v ;
          if(dump_var_lookup.find(v) ==dump_var_lookup.end())
            dump_var_lookup[v] = 0 ;
          if(dump_var_lookup[v] != 0)
            oss << "."<<dump_var_lookup[v] ;
          dump_var_lookup[v]++ ;
          
          string filename = oss.str() ;
          ofstream ofile(filename.c_str(),ios::out) ;
          ofile.precision(6) ;
          sc->Print(ofile) ;
        }
      }
    }
    
    void execute_dump_var::Print(std::ostream &s) const {
      s << "dumping variables " << dump_vars << endl ;
    }
    
    
    class dump_vars_compiler : public rule_compiler {
      variableSet dump_vars ;
    public:
      dump_vars_compiler(variableSet &vars) : dump_vars(vars) {}
      virtual void compile() {}
      virtual void accept(visitor& v) {}
      virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
      virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
      virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) {
        return executeP(new execute_dump_var(dump_vars)) ;
      }
    } ;
    
  } // end of namespace
  
  void assembleVisitor::compile_dag_sched
  (std::vector<rule_compilerP> &dag_comp,
   const std::vector<digraph::vertexSet> &dag_sched,
   const rulecomp_map &rcm,
   const digraph &dag) {

    digraph dagt = dag.transpose() ;
    if(dag_sched.size() == 0)
      return ;

    for(unsigned int i=0;i<dag_sched.size();++i) {
      //Loci::debugout << " in comp_dag.cc dag_sched[i] = " << dag_sched[i] << endl ;
      variableSet vars = extract_vars(dag_sched[i]) ;
      //Loci::debugout << " in comp_dag.cc vars = " << vars  << endl ;
      ruleSet rules = extract_rules(dag_sched[i]) ;
      if(rules == EMPTY && i+1<dag_sched.size()) {
        ++i ;
        vars += extract_vars(dag_sched[i]) ;
        rules = extract_rules(dag_sched[i]) ;
      }

      variableSet barrier_vars, reduce_vars,singleton_vars,all_vars ;
      variableSet::const_iterator vi ;
      
      std::map<variable,ruleSet> reduce_info ;
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(dagt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        bool recursive = false ;
        bool unit_rule_exists = false ;
        for(ri=var_rules.begin();ri!=var_rules.end();++ri)
          if(ri->get_info().rule_class != rule::INTERNAL) {
            use_rules += *ri ;
            rule_implP rimp = ri->get_rule_implP() ;
            if(rimp->get_rule_class() == rule_impl::POINTWISE)
              pointwise = true ;
            if(rimp->get_rule_class() == rule_impl::UNIT ||
               rimp->get_rule_class() == rule_impl::APPLY)
              reduction = true ;
            if(rimp->get_rule_class() == rule_impl::UNIT)
              unit_rule_exists = true ;
            if(rimp->get_rule_class() == rule_impl::SINGLETON)
              singleton = true ;
          } else {
            if((ri->sources() & ri->targets()) != EMPTY)
              recursive = true ;
          }
        
        WARN(reduction && pointwise || pointwise && singleton ||
             reduction && singleton) ;

        if((use_rules != EMPTY)) {
          if(pointwise && !recursive && (vi->get_info().name != "OUTPUT")) {
            barrier_vars += *vi ;
          }
          if(pointwise && recursive && (vi->get_info().name != "OUTPUT")) {
            all_vars += *vi ;
          }
          if(reduction) {
            reduce_info[*vi] = use_rules ;
            if(reduction && unit_rule_exists)
              reduce_vars += *vi ;
          } if(singleton) {
            singleton_vars += *vi ;
          }
        }

      }

      all_vars += barrier_vars ;
      dag_comp.push_back(new barrier_compiler(barrier_vars)) ;
      
      all_vars += singleton_vars ;
      
      if(singleton_vars != EMPTY)
        dag_comp.push_back(new singleton_var_compiler(singleton_vars)) ;

      all_vars += reduce_vars;
      
      if(reduce_vars != EMPTY) {
        std::map<variable,ruleSet>::const_iterator xi ;
        variableSet vars ;
        for(xi=reduce_info.begin();xi!=reduce_info.end();++xi) {
          vars += xi->first ;
          rule unit_rule ;
          CPTR<joiner> join_op = CPTR<joiner>(0) ;
          ruleSet::const_iterator ri ;
          for(ri=xi->second.begin();ri!=xi->second.end();++ri) {
            if(ri->get_rule_implP()->get_rule_class() == rule_impl::UNIT)
              unit_rule = *ri ;
            else if(ri->get_rule_implP()->get_rule_class() == rule_impl::APPLY){
              //              if(join_op == 0)
                join_op = ri->get_rule_implP()->get_joiner() ;
              //#define TYPEINFO_CHECK
#ifdef TYPEINFO_CHECK
              else
                if(typeid(*join_op) !=
                   typeid(*(ri->get_rule_implP()->get_joiner()))) {
                  cerr << "Warning:  Not all apply rules for variable " << xi->first << " have identical join operations!" << endl ;
                }
#endif
            } else {
              cerr << "Warning: reduction variable " << xi->first
                   << " has a non-reduction rule contributing to its computation,"
                   << endl << "offending rule is " << *ri << endl ;
            }
          }
          if(join_op == 0) {
            cerr << "unable to find any apply rules to complete the reduction defined by rule"
                 << endl
                 << unit_rule << endl ;
          }
          FATAL(join_op == 0) ;
          storeRepP sp = join_op->getTargetRep() ;
          if(sp->RepType() == PARAMETER) {
            dag_comp.push_back(new reduce_param_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
          else if (sp->RepType() == BLACKBOX) {
	    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
          }else {
            dag_comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
        }
      }

      
      if(check_dump_vars.ok()) {
        if(all_vars != EMPTY) 
          dag_comp.push_back(new dump_vars_compiler(all_vars)) ;
      }
      
      if(rules != EMPTY) {
        ruleSet::const_iterator ri ;
        for(ri=rules.begin();ri!=rules.end();++ri) {
          rulecomp_map::const_iterator rmi ;
          rmi = rcm.find(*ri) ;
          FATAL(rmi == rcm.end()) ;
          dag_comp.push_back(rmi->second) ;
        }
      }
    }
  }

  void assembleVisitor::visit(loop_compiler& lc) {
    compile_dag_sched(lc.collapse_comp,lc.collapse_sched,
                      lc.rule_compiler_map,lc.loop_gr) ;
    compile_dag_sched(lc.advance_comp,lc.advance_sched,
                      lc.rule_compiler_map,lc.loop_gr) ;
  }
  
  void assembleVisitor::visit(dag_compiler& dc) {
    compile_dag_sched(dc.dag_comp,dc.dag_sched,
                      dc.rule_compiler_map,dc.dag_gr) ;
  }
  
  void assembleVisitor::visit(conditional_compiler& cc) {
    compile_dag_sched(cc.dag_comp,cc.dag_sched,
                      cc.rule_compiler_map,cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////
  // topDownInfoVisitor
  void topDownInfoVisitor::visit_gr(const digraph& gr,int id) {
    // first we get the targets of the graph
    digraph::vertexSet targets ;
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    for(digraph::vertexSet::const_iterator vi=allvertices.begin();
        vi!=allvertices.end();++vi) {
      digraph::vertexSet outcome = gr[*vi] ;
      if(outcome == EMPTY)
        targets += *vi ;
    }
    // extract the variables
    variableSet target_vars = extract_vars(targets) ;

    // fill the ownership table if possible
    if(inter_node_vars == EMPTY)
      inter_node_vars = target_vars ;
    else {
      for(variableSet::const_iterator vi=target_vars.begin();
          vi!=target_vars.end();++vi) {
        if(inter_node_vars.inSet(*vi))
          ownership_table[*vi] = id ;
        else
          inter_node_vars += *vi ;
      }
    }
    // fill the recurrence table
    ruleSet rules = extract_rules(allvertices) ;
    for(ruleSet::const_iterator ri=rules.begin();
        ri!=rules.end();++ri) {
      if(ri->type() == rule::INTERNAL) {
        if( (ri->get_info().qualifier() == "promote") ||
            (ri->get_info().qualifier() == "generalize")) {
          recurrence_vars_table[*(ri->targets().begin())] =
            *(ri->sources().begin()) ;
          recurrence_vars_table[*(ri->sources().begin())] =
            *(ri->targets().begin()) ;
        }
      }
    }
    // fill the keep_vars
    if(id == 0) {
      digraph grt = gr.transpose() ;
      digraph::vertexSet sources ;
      for(digraph::vertexSet::const_iterator vi=allvertices.begin();
          vi!=allvertices.end();++vi) {
        digraph::vertexSet income = grt[*vi] ;
        if(income == EMPTY)
          sources += *vi ;
      }
      keep_vars = extract_vars(sources) ;
      keep_vars += target_vars ;
    }
  }

  void topDownInfoVisitor::visit(loop_compiler& lc) {
    visit_gr(lc.loop_gr,lc.cid) ;
  }

  void topDownInfoVisitor::visit(dag_compiler& dc) {
    visit_gr(dc.dag_gr,dc.cid) ;
  }

  void topDownInfoVisitor::visit(conditional_compiler& cc) {
    visit_gr(cc.cond_gr,cc.cid) ;
  }

  void topDownInfoVisitor::print_info(ostream& os) {
    os << "variable:" << "\towner:" << endl ;
    for(map<variable,int>::const_iterator mi=ownership_table.begin();
        mi!=ownership_table.end();++mi) {
      os << mi->first << ":\t" << mi->second << endl ;
    }
    os << endl ;
    os << "variable:" << "\trecurrence variable:" << endl ;
    for(map<variable,variable>::const_iterator mi=recurrence_vars_table.begin();mi!=recurrence_vars_table.end();++mi) {
      os << mi->first << ":\t" << mi->second << endl ;
    }
    os << endl ;
    os << "do not delete the following variables:" << endl ;
    os << keep_vars << endl ;
  }

  ///////////////////////////////////////////////////////////////////
  // graphEditVisitor
  graphEditVisitor::graphEditVisitor(const map<variable,int>& ot,
                                     const map<variable,variable>& rt,
                                     const variableSet& kv):
    ownership_table(ot),recurrence_vars_table(rt),keep_vars(kv)
  {
    for(map<variable,int>::const_iterator mi=ownership_table.begin();
        mi!=ownership_table.end();++mi) {
      ot_vars += mi->first ;
    }
    for(map<variable,variable>::const_iterator mi=recurrence_vars_table.begin();mi!=recurrence_vars_table.end();++mi) {
      rvt_vars += mi->first ;
    }
    /*
    keep_vars += variable("simulation_finished{n}") ;
    keep_vars += variable("$n{n}") ;
    keep_vars += variable("T{n}") ;

    keep_vars += variable("gradDotN(T){n}") ;
    keep_vars += variable("stable_sum") ;
    keep_vars += variable("Tface{n}") ;
    keep_vars += variable("nodal_temp_sum{n}") ;
    keep_vars += variable("Tnode{n}") ;
    keep_vars += variable("heat_flux{n}") ;
    keep_vars += variable("flux_sum{n}") ;
    keep_vars += variable("rhs{n}") ;
    */
    //keep_vars += variable("heat_flux{n}") ;
  }

  namespace {
    rule create_rule(variable sv, variable tv, string qualifier) {
      ostringstream oss ;
      oss << "source(" << sv << ')' ;
      oss << ",target(" << tv << ')' ;
      oss << ",qualifier(" << qualifier << ')' ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
    
    rule create_rule(variableSet source, variableSet target, string qualifier) {
      ostringstream oss ;
      oss << "source(" << source << ')' ;
      oss << ",target(" << target << ')' ;
      oss << ",qualifier(" << qualifier << ")" ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }

    void create_digraph_dot_file(const digraph &dg, const char* fname)
    {
      digraph dgt = dg.transpose() ;
      digraph::vertexSet allvertices = dg.get_all_vertices() ;
      digraph::vertexSet::const_iterator ri ;
      ofstream outf(fname) ;
      
      outf<<"digraph G {\n" ;
      //outf<<"size = \"8.5,11\";\n" ;
      
      for(ri=allvertices.begin();ri!=allvertices.end();++ri) {

        if(*ri == 0)
          continue ;
        
        digraph::vertexSet outvertices = dg[*ri] ;
        digraph::vertexSet incomevertices = dgt[*ri] ;

        if(*ri < 0) {
          rule r(*ri) ;
          if(r.type() == rule::INTERNAL) {
            outf<<"\""<<r<<"\""
                <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
          }
          else {
            outf<<"\""<<r<<"\""<<"[shape=box,style=filled,color=gold];\n" ;
          }
        }
        else {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
        }
        if(incomevertices == EMPTY) {
          if(*ri > 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<r<<"\""<<"[style=filled,color=red];\n" ;
          }
        }
        if(outvertices == EMPTY) {
          if(*ri > 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<r<<"\""<<"[style=filled,color=blue];\n" ;
          }
          continue ;
        }
        
        if(*ri < 0) {
          rule r(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r2(*ii) ;
              outf<<"\""<<r<<"\""<<" -> "<<"\""<<r2<<"\""<<";\n" ;
            }else {
              variable v(*ii) ;
              outf<<"\""<<r<<"\""<<" -> "<<"\""<<v<<"\""<<";\n" ;
            }
          }
        }else {
          variable v(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<r<<"\""<<";\n" ;
            }else {
              variable v2(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""<<";\n" ;
            }
          }        
        }
      }
      outf<<"}\n" ;
      outf.close() ;
    }
  } // end of un-named namespace
  
  void graphEditVisitor::edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;
    
    // this map is to keep track which variables have been allocated
    // in case of multiple rule generates the same target
    map<variable,rule> allocated_vars ;
    
    // this map is to keep track which variables have been freed
    map<variable,rule> freed_vars ;

    // looping over all the rules and edit the graph accordingly
    for(ruleSet::const_iterator ruleIter=rules.begin();ruleIter!=rules.end();++ruleIter) {
      // we do not process these internal rules here
      if(ruleIter->type() == rule::INTERNAL) {
        if( (ruleIter->get_info().qualifier() == "promote") ||
            (ruleIter->get_info().qualifier() == "generalize") ||
            (ruleIter->get_info().qualifier() == "priority")
            )
          continue ;
      }

      // get the sources and targets of the rule
      variableSet sources, targets ;
      sources = ruleIter->sources() ;
      targets = ruleIter->targets() ;
      // first we allocate the targets, (loop over it)
      for(variableSet::const_iterator tvarIter=targets.begin();
          tvarIter!=targets.end();++tvarIter) {
        
        bool found_in_table ;
        int id_in_table ;

        map<variable,int>::const_iterator ownerIter ;
        ownerIter = ownership_table.find(*tvarIter) ;

        if(ownerIter == ownership_table.end()) {
          found_in_table = false ;
          id_in_table = -1 ;
        } else {
          found_in_table = true ;
          id_in_table = ownerIter->second ;
        }
        
        // we allocate the variable if it is not in ownership_table
        // or if it is in the table, we compare the id with that in
        // the table, if they match, we allocate it
        if( (!found_in_table) ||
            ((found_in_table)&&(id_in_table==id))
            ) {
          map<variable,rule>::const_iterator smi ;
          smi = allocated_vars.find(*tvarIter) ;
          // if it hasn't been allocated
          if(smi == allocated_vars.end()) {
            // create a rule that points to the rule
            // and edit the rulecomp_map to insert a new allocate compiler
            // for the rule
            variable sv("SOURCE") ;
            variable tv("SINK") ;
            
            rule alloc_rule = create_rule(sv,*tvarIter,"ALLOCATE") ;
            //            rule alloc_rule = create_rule(sv,tv,"ALLOCATE") ;
            gr.add_edge(alloc_rule.ident(),ruleIter->ident()) ;
            // modify the rulecomp_map
            rcm[alloc_rule] = new allocate_var_compiler(*tvarIter) ;
            
            allocated_vars[*tvarIter] = alloc_rule ;
          }
          else {
            // this means the variable has also been allocated by another
            // rule we therefore connect the already existed allocation
            // rule to this rule
            gr.add_edge( (smi->second).ident(),ruleIter->ident()) ;
          }
        }
      } // end of for (targets)

      // now we decorate the graph to include deallocation of variables
      for(variableSet::const_iterator svarIter=sources.begin();
          svarIter!=sources.end();++svarIter) {
        // if a variable belongs to the keep_vars set, then don't delete it
        if(keep_vars.inSet(*svarIter)) continue ;
        // if a variable is a recurrence variable, we don't delet it
        if(rvt_vars.inSet(*svarIter)) continue ;
        
        // get the ident of the targets of the rule
        digraph::vertexSet targets_ident ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi)
          targets_ident += vi->ident() ;

        map<variable,rule>::const_iterator smi ;
        smi = freed_vars.find(*svarIter) ;

        // if it hasn't been freed
        if(smi == freed_vars.end()) {
          // create a rule that from the rule's targets
          // and edit the rulecomp_map to insert a new free compiler
          // for the rule
          variableSet tv ;
          tv += *svarIter ;
          variableSet sv(EMPTY) ;
          
          rule free_rule = create_rule(sv,tv,"FREE") ;
          
          gr.add_edges(targets_ident,free_rule.ident()) ;
          // modify the rulecomp_map
          rcm[free_rule] = new free_var_compiler(*svarIter) ;

          // save the rule
          freed_vars[*svarIter] = free_rule ;
        }
        else {
          // this means another rule also deletes this variable
          // we make an edge to maintain the order
          gr.add_edges(targets_ident,(smi->second).ident()) ;
        }
        
      } // end of for (sources)

    } // end of for (ruleSet::const_iterator) loop
  }

  void graphEditVisitor::visit(loop_compiler& lc) {
    // first find out if there are any variables
    // in both collapse and advance graph
    digraph::vertexSet collapse_gr_ver = lc.collapse_gr.get_all_vertices() ;
    digraph::vertexSet advance_gr_ver = lc.advance_gr.get_all_vertices() ;
    variableSet common = extract_vars(collapse_gr_ver & advance_gr_ver) ;
    // we do not delete common in collapse and let the advance part
    // delete it, since in every loop, we first execute collapse part
    // and then the advance part
    keep_vars += common ;
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,lc.cid) ;
    //keep_vars -= common ; 
    
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }
  
  void graphEditVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }
  
  void graphEditVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }

  namespace {
    bool is_super_node(int rid)
    {
      if(rid >= 0) // a variable
        return false ;
      rule r(rid) ;
      string rqualifier = r.get_info().qualifier() ;
      
      return (rqualifier.substr(0,2) == "SN") ;
    }

    bool is_super_node(const ruleSet::const_iterator& ruleIter) {
      string rqualifier = ruleIter->get_info().qualifier() ;

      return (rqualifier.substr(0,2) == "SN") ;
    }

    bool is_super_node(const rule& r) {
      string rqualifier = r.get_info().qualifier() ;

      return (rqualifier.substr(0,2) == "SN") ;
    }

    int get_supernode_num(const rule& r) {
      string rqualifier = r.get_info().qualifier() ;
      string head = rqualifier.substr(0,2) ;
      if(head != "SN") {
        cerr << "get_supernode_num error! pass in rule is not a super node!"
             << endl ;
        exit(-1) ;
      }
      
      string number = rqualifier.substr(2,rqualifier.size()-2) ;
      stringstream ss ;
      ss << number ;
      int ret ;
      ss >> ret ;
      
      return ret ;
    }

    template<typename T>
    inline bool inSet(const set<T>& s, const T& elem) {
      typename set<T>::const_iterator si ;
      si = s.find(elem) ;
      if(si == s.end())
        return false ;
      return true ;
    }

    // is there a path from source vertex to target vertex in gr?
    bool has_path(const digraph& gr, int source, int target) {
      bool path = false ;
      digraph::vertexSet working = gr[source] ;
      while(working != EMPTY) {
        if(working.inSet(target)) {
          path = true ;
          break ;
        }
        digraph::vertexSet children ;
        for(digraph::vertexSet::const_iterator vi=working.begin();
            vi!=working.end();++vi) {
          children += gr[*vi] ;
        }
        working = children ;
      }

      return path ;
    }

    // given a variableSet, convert them to vertexSet
    inline digraph::vertexSet get_vertexSet(const variableSet& vars) {
      digraph::vertexSet ret ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
        ret += vi->ident() ;

      return ret ;
    }

    // given a ruleSet, convert them to vertexSet
    inline digraph::vertexSet get_vertexSet(const ruleSet& rules) {
      digraph::vertexSet ret ;
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
        ret += ri->ident() ;

      return ret ;
    }

    // get the init source variableSet from a target variable and
    // the recurrence target to source table, or
    // get the final target variableSet from a source variable and
    // the recurrence source to target table
    inline
    variableSet get_leaf_recur_vars(const map<variable,variableSet>& t,
                                    const variable& v) {
      map<variable,variableSet>::const_iterator found ;
      variableSet ret ;
      list<variable> working ;

      // first check if v is in the table
      found = t.find(v) ;
      // if v is not recurrence variable, return empty set
      if(found == t.end())
        return variableSet(EMPTY) ;

      // do a depth first search
      working.push_front(v) ;
      while(!working.empty()) {
        variable cur = working.front() ;
        working.pop_front() ;
        
        variableSet tmp ;
        found = t.find(cur) ;
        
        if(found != t.end())
          tmp = found->second ;
        else
          ret += cur ;
        
        for(variableSet::const_iterator vi=tmp.begin();
            vi!=tmp.end();++vi)
          working.push_front(*vi) ;
      }

      return ret ;
    }

    // get all the source variables from a target variable and
    // the recurrence target to source table, or
    // get all the target variables from a source variable and
    // the recurrence source to target table
    inline
    variableSet get_all_recur_vars(const map<variable,variableSet>& t,
                                   const variable& v) {
      map<variable,variableSet>::const_iterator found ;
      variableSet ret, working ;
      
      // first check if v is in the table
      found = t.find(v) ;
      // if v is not recurrence variable, return empty set
      if(found == t.end())
        return variableSet(EMPTY) ;

      // do a bread first search, and add up all the results
      working += v ;
      while(working != EMPTY) {
        variableSet tmp ;
        for(variableSet::const_iterator vi=working.begin();
            vi!=working.end();++vi) {
          found = t.find(*vi) ;
          if(found != t.end())
            tmp += found->second ;
        }
        ret += tmp ;
        working = tmp ;
      }

      return ret ;
    }

    // given a variable, if it is a recurrence variable, then get
    // all other leaf target variables that all its source can reach
    inline
    variableSet get_all_leaf_target(const map<variable,variableSet>& t2s,
                                    const map<variable,variableSet>& s2t,
                                    const variable& v) {
      map<variable,variableSet>::const_iterator found ;
      variableSet ret, working ;

      working += v ;
      while(working != EMPTY) {
        variableSet tmp ;
        for(variableSet::const_iterator vi=working.begin();
            vi!=working.end();++vi) {
          found = t2s.find(*vi) ;
          if(found != t2s.end()) {
            variableSet cur = found->second ;
            tmp += cur ;
            for(variableSet::const_iterator vi2=cur.begin();
                vi2!=cur.end();++vi2)
              ret += get_leaf_recur_vars(s2t,*vi2) ;
          }
        }
        working = tmp ;
      }

      ret -= v ;
      return ret ;
    }

  } // end-of-unnamed namespace
  
  /////////////////////////////////////////////////////////////
  // graphVisualizeVisitor
  /////////////////////////////////////////////////////////////
  void graphVisualizeVisitor::visit(loop_compiler& lc) {
    if(Loci::MPI_rank==0) {
      cerr << "the looping part graph (collapse and advance), id: "
           << lc.cid << endl ;
      create_digraph_dot_file(lc.loop_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;

      cerr << "\tthe collapse part graph" << endl ;
      create_digraph_dot_file(lc.collapse_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
      
      cerr << "\tthe advance part graph" << endl ;
      create_digraph_dot_file(lc.advance_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

  void graphVisualizeVisitor::visit(dag_compiler& dc) {
    if(Loci::MPI_rank==0) {
      cerr << "the dag graph, id: " << dc.cid << endl ;
      create_digraph_dot_file(dc.dag_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

  void graphVisualizeVisitor::visit(conditional_compiler& cc) {
    if(Loci::MPI_rank==0) {
      cerr << "the conditional graph, id: " << cc.cid << endl ;
      create_digraph_dot_file(cc.cond_gr,"visual_graph.dot") ;
      system("dotty visual_graph.dot") ;
      system("rm -fr visual_graph.dot") ;
    }
  }

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
    // we remove recurrence target variables from the working_vars.
    working_vars -= recur_target_vars ;

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
    // and the common variables for the loop in this graph.
    variableSet all_rot_vars ;
    variableSet all_common_vars ;
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

          found = loop_common_table.find(id) ;
          FATAL(found == loop_common_table.end()) ;
          all_common_vars += found->second ;
        }
      }
    }
    all_rot_vars -= recur_target_vars ;
    all_rot_vars -= allocated_vars ;
    all_common_vars -= recur_target_vars ;
    all_common_vars -= allocated_vars ;
    
    alloc_here += all_rot_vars ;
    alloc_here += all_common_vars ;

    
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
    // variables and loop common variables. we allocate them here
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      variableSet rotate_vars ;
      variableSet common_vars ;
      variableSet allocated ;
      
      if(is_super_node(ruleIter)) {
        int id = get_supernode_num(*ruleIter) ;
        if(inSet(loop_sn,id)) {
          // get the rotate list variables
          map<int,variableSet>::const_iterator found ;
          found = rotate_vtable.find(id) ;
          FATAL(found == rotate_vtable.end()) ;
          rotate_vars += found->second ;

          // get the common variables
          found = loop_common_table.find(id) ;
          FATAL(found == loop_common_table.end()) ;
          common_vars += found->second ;

          for(variableSet::const_iterator vi=alloc_vars.begin();
              vi!=alloc_vars.end();++vi) {
            // rotate list vars & common vars
            if(rotate_vars.inSet(*vi) || common_vars.inSet(*vi)) {
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
  variableSet deleteInfoVisitor::get_start_info(const digraph& gr,int id) {

    // first check the deletion of recurrence target variables
    check_recur_target(gr,id) ;
    
    variableSet working_vars ;

    // obtain all the possible variables in the graph
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      working_vars += ri->sources() ;
    
    working_vars += extract_vars(allvertices) ;
    
    // post processing (filter) of working_vars
    post_proc_working_vars(working_vars) ;

    return working_vars ;    
  }

  // post processing of working variables
  void deleteInfoVisitor::post_proc_working_vars(variableSet& working_vars) {
    // we reserve the user queried target variables, we don't delete them
    working_vars -= queried_target ;
    // we don't process recurrence source variables
    // that are possible in the working_vars
    // e.g. x -> x{n}, x is the recurrence source variables
    // i.e. we always don't delete recurrence source variables
    working_vars -= recur_source_vars ;
    
    // we exclude any recurrence target variables that is
    // in the working_vars and has been deleted already
    variableSet target_cand ;
    target_cand = variableSet(working_vars & recur_target_vars) ;
    target_cand -= deleted_vars ;
    for(variableSet::const_iterator vi=target_cand.begin();
        vi!=target_cand.end();++vi) {
      // search the redirect table
      map<variable,int>::const_iterator found ;
      found = redirect_table.find(*vi) ;
      if(found != redirect_table.end()) {
        if(!deleted_vars.inSet(*vi)) {
          deleted_vars += *vi ;
          delete_table[found->second] += *vi ;
        }
      }
    }

    // and we don't process variables that have been deleted
    // or been reserved
    working_vars -= deleted_vars ;
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
      // all the other same target variables
      variableSet other_targets ;

      all_recur_sources = get_all_recur_vars(recur_vars_t2s,*varIter) ;
      other_targets = get_all_leaf_target(recur_vars_t2s,
                                          recur_vars_s2t,
                                          *varIter) ;

      // some verification to the other_targets
      if(other_targets != EMPTY) {
        // if any of the other_targets is deleted, we don't
        // need to process this varible
        bool if_continue = false ;
        for(variableSet::const_iterator vi=other_targets.begin();
            vi!=other_targets.end();++vi)
          if(deleted_vars.inSet(*vi)) {
            deleted_vars += *varIter ;
            if_continue = true ;
            break ;
          }
        if(if_continue)
          continue ;
      }
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

    // we then get the common variable shared between the
    // advance and collapse part of the loop
    map<int,variableSet>::const_iterator found ;
    variableSet common ;
    found = loop_common_table.find(lc.cid) ;
    FATAL(found == loop_common_table.end()) ;
    common = found->second ;

    // gather variables in the rotate list
    variableSet rotate_vars ;

    found = rotate_vtable.find(lc.cid) ;
    FATAL(found == rotate_vtable.end()) ;
    rotate_vars = found->second ;

    // we then gather deletion information for the collapse part
    variableSet working_vars = get_start_info(lc.collapse_gr,-lc.cid) ;
    // we delete the common shared variables until the last iteration
    // finishes and inside the conditional super node
    working_vars -= common ;
    // we don't delete variabls in rotate lists, we defer that until
    // the last iteration finishes
    working_vars -= rotate_vars ;
    
    looping_algr(working_vars,lc.collapse_gr,-lc.cid,0) ;
    
    // we gather deletion information for the advance part
    working_vars = get_start_info(lc.advance_gr,lc.cid) ;
    working_vars -= common ;
    working_vars -= rotate_vars ;

    looping_algr(working_vars,lc.advance_gr,lc.cid,0) ;

    // we now schedule the deletion of common and rotate_vars

    // first get the collapse node of this loop
    map<int,int>::const_iterator found2 ;
    found2 = loop_ctable.find(lc.cid) ;
    FATAL(found2 == loop_ctable.end()) ;
    int collapse_id = found2->second ;

    working_vars = variableSet(common + rotate_vars) ;
    post_proc_working_vars(working_vars) ;

    looping_algr(working_vars,lc.loop_gr,collapse_id,2) ;

  }

  void deleteInfoVisitor::visit(dag_compiler& dc) {
    // if the dag has id 0, then this is the top level
    // graph, i.e. the computation begins from this dag.
    // we gather the input variables, we don't delete them.
    if(dc.cid == 0) {
      // get the transposed graph
      digraph grt = dc.dag_gr.transpose() ;
      digraph::vertexSet all_vertices = grt.get_all_vertices() ;
      
      digraph::vertexSet input ;
      for(digraph::vertexSet::const_iterator vi=all_vertices.begin();
          vi!=all_vertices.end();++vi) {
        if(grt[*vi] == EMPTY)
          input += *vi ;
      }

      variableSet input_vars = extract_vars(input) ;

      deleted_vars += input_vars ;

      // we also reserve all the recurrence targets of input variables
      variableSet input_vars_target ;
      for(variableSet::const_iterator vi=input_vars.begin();
          vi!=input_vars.end();++vi) {
        input_vars_target += get_leaf_recur_vars(recur_vars_s2t,*vi) ;
      }

      deleted_vars += input_vars_target ;
    }

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
                                              
  // function that checks the deletion of recurrence target vars
  void deleteInfoVisitor::check_recur_target(const digraph& gr,int id) {
    // we will also need to do an additional search
    // for one special case for the recurrence variables.
    // if any variable in the graph is a recurrence source
    // variable and it final recurrence target variable
    // is not in the graph. we need to check if there is
    // any rules that are reachable from these recurrence
    // source variables. if yes, we need to delete the
    // recurrence target variable here in this graph

    digraph::vertexSet allv = gr.get_all_vertices() ;
    // get all the leaf node
    digraph::vertexSet leaf_node ;
    for(digraph::vertexSet::const_iterator vi=allv.begin();
        vi!=allv.end();++vi) {
      if(gr[*vi] == EMPTY)
        leaf_node += *vi ;
    }
    
    ruleSet rules = extract_rules(allv) ;
    variableSet allvars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      allvars += ri->sources() ;

    allvars &= recur_source_vars ;

    // we don't need to process leaf variables
    allvars -= extract_vars(leaf_node) ;
    
    for(variableSet::const_iterator vi=allvars.begin();
        vi!=allvars.end();++vi) {
      // find out the final recurrence target variable
      variableSet targets = get_leaf_recur_vars(recur_vars_s2t,*vi) ;

      for(variableSet::const_iterator vi2=targets.begin();
          vi2!=targets.end();++vi2) {
        // see if it is in the graph
        // if in the graph, we don't do anything
        if(allv.inSet(vi2->ident()))
          continue ;
        // see if it is been seen before
        if(processed_targets.inSet(*vi2))
          continue ;
        // if it has already been deleted, we don't delete it
        if(deleted_vars.inSet(*vi2))
          continue ;
        // not in the graph, have not been processed yet
        variableSet all_sources = get_all_recur_vars(recur_vars_t2s,*vi2) ;
        
        // get all the source in the graph
        digraph::vertexSet all_source_vertices = get_vertexSet(all_sources) ;
        all_source_vertices &= allv ;
        
        // get the rules all_source can reach
        ruleSet others ;
        for(digraph::vertexSet::const_iterator
              veri=all_source_vertices.begin();
            veri!=all_source_vertices.end();++veri) {
          ruleSet tmp = extract_rules(gr[*veri]) ;
          others += tmp ;
          for(ruleSet::const_iterator ruleIter=tmp.begin();
              ruleIter!=tmp.end();++ruleIter)
            if(!is_super_node(ruleIter))
              others += extract_rules(gr[ruleIter->ident()]) ;
        }
        FATAL(others == EMPTY) ;
        // see if we must delete the target variable here
        if(!let_it_go(gr,others)) {
          // find out where to put the deletion

          // get all the parent loop super node (include itself if apply)
          vector<int> all_loop_pnode ;
          int pnode ;
          if(id < 0)
            pnode = -id ;
          else
            pnode = id ;
          if(inSet(loop_sn,pnode))
            all_loop_pnode.push_back(pnode) ;
          
          // we search for all loops
          map<int,int>::const_iterator found ;
          while(pnode != 0) {
            found = pnode_table.find(pnode) ;
            FATAL(found == pnode_table.end()) ;
            pnode = found->second ;
            
            if(inSet(loop_sn,pnode))
              all_loop_pnode.push_back(pnode) ;
          }

          int delete_graph_id ;
          if(all_loop_pnode.empty())
            delete_graph_id = id ;
          else {
            map<int,variableSet>::const_iterator found_alloc ;
            int index, found_loop_id ;
            for(index=0;index!=all_loop_pnode.size();++index) {
              // get the loop's allocation
              found_alloc = loop_alloc_table.find(all_loop_pnode[index]) ;
              FATAL(found_alloc == loop_alloc_table.end()) ;
              if(found_alloc->second.inSet(*vi2))
                break ;
            }
            if(index == 0)
              found_loop_id = id ;
            else
              found_loop_id = all_loop_pnode[index-1] ;
            
            // find out the loop collapse node id
            map<int,int>::const_iterator found ;
            found = loop_ctable.find(found_loop_id) ;
            FATAL(found == loop_ctable.end()) ;
            delete_graph_id = found->second ;
          }
                    
          processed_targets += *vi2 ;
          redirect_table[*vi2] = delete_graph_id ;
        }
      }
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
    int supernode_id ;
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

  // let_it_go for check_recur_target
  bool deleteInfoVisitor::let_it_go(const digraph& gr, ruleSet rules) {
    if(rules == EMPTY) {
      return false ;
    }
    
    int supernode_num = 0 ;
    int supernode_id ;
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
    if( (rules.size() == 1) && (supernode_num == 1))
      return true ;
    if(supernode_num != 1)
      return false ;
    // the remaining case must be there are only one super node
    // and at least one non super node, we check if there is path
    // from each non super node to the super node in the graph
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
        if( (ruleIter->get_info().qualifier() == "promote") ||
            (ruleIter->get_info().qualifier() == "generalize") ||
            (ruleIter->get_info().qualifier() == "priority")
            ) {
          variable target,source ;
          target = *(ruleIter->targets().begin()) ;
          source = *(ruleIter->sources().begin()) ;

          recur_vars_t2s[target] += source ;
          recur_vars_s2t[source] += target ;
          recur_source_vars += source ;
          recur_target_vars += target ;
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

  ostream& allocDeleteStat::report(ostream& s) {
    variableSet d_a, a_d ;
    variableSet::const_iterator vi,vii ;
    map<variable,variable>::const_iterator mi ;

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

    variableSet miss_delete ;
    for(vi=deleted_vars.begin();vi!=deleted_vars.end();++vi) {
      if(recur_source_vars.inSet(*vi))
        miss_delete += *vi ;
    }
    if(miss_delete != EMPTY) {
      s << endl ;
      s << "These variables should not been deleted, because they are "
        << "recurrence source variables: " << endl << miss_delete << endl ;
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
      variableSet final_targets = get_leaf_recur_vars(recur_vars_s2t,*vi) ;
      
      if(final_targets != EMPTY) {
        for(vii=final_targets.begin();vii!=final_targets.end();++vii)
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

    // and then compute the common variables between
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
  namespace {
    inline bool is_dg_empty(const digraph& gr) {
      digraph::vertexSet allv = gr.get_all_vertices() ;
      return (allv == EMPTY) ;
    }
  }

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
  namespace {
    // return if a rule is a generalize, promote or priority rule
    inline bool is_recur_rule(const ruleSet::const_iterator& ruleIter) {
      if(ruleIter->type() == rule::INTERNAL) {
        if( (ruleIter->get_info().qualifier() == "promote") ||
            (ruleIter->get_info().qualifier() == "generalize") ||
            (ruleIter->get_info().qualifier() == "priority")
            )
          return true ;
      }
      return false ;
    }

    inline bool is_internal_rule(const ruleSet::const_iterator& ruleIter) {
      return (ruleIter->type() == rule::INTERNAL) ;
    }

  } // end-of-namespace

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
    if(dc.cid == 0) {
      // get the transposed graph
      digraph grt = dc.dag_gr.transpose() ;
      digraph::vertexSet all_vertices = grt.get_all_vertices() ;
      
      digraph::vertexSet input ;
      for(digraph::vertexSet::const_iterator vi=all_vertices.begin();
          vi!=all_vertices.end();++vi) {
        if(grt[*vi] == EMPTY)
          input += *vi ;
      }

      typed_vars += extract_vars(input) ;

      variableSet add ;
      for(variableSet::const_iterator vi=typed_vars.begin();
          vi!=typed_vars.end();++vi) {
        variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
        add += all_targets ;
      }
      typed_vars += add ;
    }
    
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

  //////////////////////////////////////////////////////////
  // function to decorate the dependency graph
  //////////////////////////////////////////////////////////
  void deco_depend_gr(digraph& gr,const variableSet& given) {
    // first get all the vertices
    digraph::vertexSet allv = gr.get_all_vertices() ;
    // then we obtain the transpose of the graph
    digraph grt = gr.transpose() ;
    
    // then we collect the recurrence variables information
    variableSet recur_source_vars, recur_target_vars ;
    map<variable,variableSet> recur_vars_s2t ;
    map<variable,variableSet> recur_vars_t2s ;

    ruleSet rules = extract_rules(allv) ;
    // looping over all the rules and gather info
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      // recurrence variables are the target and source of
      // these internal rules
      if(ruleIter->type() == rule::INTERNAL) {
        if( (ruleIter->get_info().qualifier() == "promote") ||
            (ruleIter->get_info().qualifier() == "generalize") ||
            (ruleIter->get_info().qualifier() == "priority")
            ) {
          variable target,source ;
          target = *(ruleIter->targets().begin()) ;
          source = *(ruleIter->sources().begin()) ;

          recur_vars_t2s[target] += source ;
          recur_vars_s2t[source] += target ;
          recur_source_vars += source ;
          recur_target_vars += target ;
        }
      }
    }

    // variables to be evaluated during allocation and deletion
    variableSet alloc_working_vars ;
    variableSet delete_working_vars ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      alloc_working_vars += ruleIter->targets() ;
      delete_working_vars += ruleIter->sources() ;
    }
    alloc_working_vars -= recur_target_vars ;
    delete_working_vars -= recur_source_vars ;
    delete_working_vars -= given ;

    // we first do allocation
    for(variableSet::const_iterator vi=alloc_working_vars.begin();
        vi!=alloc_working_vars.end();++vi) {
      
      digraph::vertexSet source_rules = grt[vi->ident()] ;
      
      // we create a rule for allocation
      variable sv("CREATE") ;
      rule alloc_rule = create_rule(sv,*vi,"ALLOCATE") ;
      // edit graph
      gr.add_edges(alloc_rule.ident(),source_rules) ;
    }

  }

} // end of namespace Loci

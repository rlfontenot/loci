#include "comp_tools.h"
#include <vector>
using std::vector ;

#include <sys/stat.h>

#include <distribute.h>
#include <parameter.h>

#include <map>
using std::map ;


namespace Loci {

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
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) {
      return executeP(new execute_dump_var(dump_vars)) ;
    }
  } ;
  
  void compile_dag_sched(std::vector<rule_compilerP> &dag_comp,
                         const std::vector<digraph::vertexSet> &dag_sched,
                         const rulecomp_map &rcm,
                         const digraph &dag) {
    digraph dagt = dag.transpose() ;
    if(dag_sched.size() == 0)
      return ;
#ifdef TEST
    variableSet acvars = extract_vars(dag_sched[0]) ;
    vector<variableSet> vars_alloc(dag_sched.size()) ;
    vector<variableSet> vars_free(dag_sched.size()) ;
    acvars = EMPTY ;
    for(int i=0;i<dag_sched.size();++i) {
      ruleSet rules = extract_rules(dag_sched[i]) ;
      ruleSet::const_iterator ri ;
      variableSet vs ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        vs += ri->targets() ;
      vars_alloc[i] = vs - acvars ;
      //      cerr << "vars_alloc["<<i<<"] = " << vars_alloc[i] << endl ;
      acvars += vs ;
    }
#endif
    for(int i=0;i<dag_sched.size();++i) {
      variableSet vars = extract_vars(dag_sched[i]) ;
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

  
  dag_compiler::dag_compiler(rulecomp_map &rule_process, digraph dag) {

    std::vector<digraph::vertexSet> dag_sched = schedule_dag(dag) ;
    compile_dag_sched(dag_comp,dag_sched,rule_process,dag) ;
    
#ifdef DEBUG
    // sanity check, all vertices should be scheduled
    digraph::vertexSet allvertices ;
    for(int i=0;i< dag_sched.size();++i) 
      allvertices += dag_sched[i] ;
    warn(allvertices != dag.get_all_vertices()) ;
    if(allvertices != dag.get_all_vertices()) {
      digraph::vertexSet leftout = allvertices ^ dag.get_all_vertices() ;
      cerr << "leftout rules= " << extract_rules(leftout) << endl ;
      cerr << "leftout vars = " << extract_vars(leftout) << endl ;
    }
#endif
  }

  void dag_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void dag_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP dag_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    return executeP(elp) ;
  }
}


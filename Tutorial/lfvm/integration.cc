#include<Loci.h>

// time integration

class ic_rule : public pointwise_rule {
  const_param<double> T_initial ;
  const_param<double> heat_capacity ;
  store<double> e0 ;
public:
  ic_rule() {
    name_store("T_initial",T_initial) ;
    name_store("heat_capacity",heat_capacity) ;
    name_store("e0{n=0}",e0) ;
    input("T_initial,heat_capacity") ;
    output("e0{n=0}") ;
    constraint("area") ;
  }
  void calculate(Entity e) {
    e0[e] = T_initial[e]*heat_capacity[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<ic_rule> register_ic_rule ;


class advance_heat : public pointwise_rule {
  const_store<double> e0n ;
  const_param<double> time_step ;
  const_store<double> rhs ;
  store<double> e0np1 ;
public:
  advance_heat() {
    name_store("e0{n}",e0n) ;
    name_store("time_step{n}",time_step) ;
    name_store("rhs{n}",rhs) ;
    name_store("e0{n+1}",e0np1) ;
    input("time_step{n},rhs{n},e0{n}") ;
    output("e0{n+1}") ;
    constraint("e0{n}") ;
  }
  void calculate(Entity e) {
    e0np1[e] = e0n[e] + time_step[e]*rhs[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<advance_heat> register_advance_heat ;

class collapse_time : public pointwise_rule {
  const_store<double> Tnode ;
  store<double> solution ;
public:
  collapse_time() {
    name_store("Tnode{n}",Tnode) ;
    name_store("solution",solution) ;
    input("Tnode{n}") ;
    output("solution") ;
    conditional("simulation_finished{n}") ;
  }
  void calculate(Entity e) {
    solution[e] = Tnode[e] ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<collapse_time> register_collapse_time ;

class collapse_info : public singleton_rule {
  const_param<int> n ;
  const_param<double> time_step ;
  param<double> solution_time ;
public:
  collapse_info() {
    name_store("$n{n}",n) ;
    name_store("time_step{n}",time_step) ;
    name_store("solution_time",solution_time) ;
    input("$n{n},time_step{n}") ;
    output("solution_time") ;
    conditional("simulation_finished{n}") ;
  }

  virtual void compute(const sequence &seq) {
    *solution_time = double(*n)*(*time_step) ;
  }
} ;


register_rule<collapse_info> register_collapse_info ;

// Condition that determines when iteration is complete
class collapse_condition : public singleton_rule {
  const_param<int> n, max_iteration ;
  const_param<double> max_time ;
  const_param<double> time_step ;
  param<bool> simulation_finished ;
public:
  collapse_condition() {
    name_store("$n",n) ;
    name_store("max_iteration",max_iteration) ;
    name_store("time_step",time_step) ;
    name_store("max_time",max_time) ;
    name_store("simulation_finished",simulation_finished) ;
    input("$n,max_iteration,max_time,time_step") ; 
    output("simulation_finished") ;
  }
  virtual void compute(const sequence &seq) {
    *simulation_finished = (*n >= *max_iteration) 
      ||((double(*n) * *time_step) >= *max_time);
  }
} ;

register_rule<collapse_condition> register_collapse_condition ;

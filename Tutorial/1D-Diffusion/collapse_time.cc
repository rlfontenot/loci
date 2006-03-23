#include <Loci.h>

//*******************************************************************
// When simulation is finished, copy current iteration results to
// solution
//*******************************************************************

class collapse_time : public pointwise_rule {
  const_store<float> u ;
  store<float> solution ;
public:
  collapse_time() {
    name_store("u{n}",u) ;
    name_store("solution",solution) ;
    input("u{n}") ;
    output("solution") ;
    conditional("simulation_finished{n}") ;
  }
  void calculate(Entity e) {
    solution[e] = u[e] ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;
  
register_rule<collapse_time> register_collapse_time ;

//*******************************************************************
// Condition that determines when iteration is complete
//*******************************************************************

class collapse_condition : public singleton_rule {
  const_param<int> n, max_iteration ;
  param<bool> simulation_finished ;
public:
  collapse_condition() {
    name_store("$n",n) ;
    name_store("max_iteration",max_iteration) ;
    name_store("simulation_finished",simulation_finished) ;
    input("$n,max_iteration") ;
    output("simulation_finished") ;
  }
  virtual void compute(const sequence &seq) {
    *simulation_finished = (*n >= *max_iteration) ;
  }
} ;

register_rule<collapse_condition> register_collapse_condition ;

//*******************************************************************

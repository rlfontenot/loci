#include <Loci.h>

class life_step : public pointwise_rule {
  const_Map N, S, E, W ;
  const_store<int> life_n ;
  store<int> life_nP1 ;
public:
  life_step() {
    name_store("N",N) ;
    name_store("S",S) ;
    name_store("E",E) ;
    name_store("W",W) ;
    name_store("life{N}",life_n) ;
    name_store("life{N+1}",life_nP1) ;
    input("N->life{N},E->N->life{N},W->N->life{N}") ;
    input("S->life{N},E->S->life{N},W->S->life{N}") ;
    input("E->life{N}") ;
    input("W->life{N}") ;
    input("life{N}") ;
    output("life{N+1}") ;
  }
  inline void calculate(Entity e) {
    int sum = 
      life_n[N[e]]+life_n[N[E[e]]]+life_n[N[W[e]]]+
      life_n[S[e]]+life_n[S[E[e]]]+life_n[S[W[e]]]+
      life_n[E[e]]+
      life_n[W[e]] ;
    if(life_n[e]==0) {
      life_nP1[e] = (sum == 3) ;
    } else {
      life_nP1[e] = (sum == 2 || sum == 3) ;
    }
  }
  void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<life_step> register_life_step ;

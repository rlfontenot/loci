#include <Loci.h>

//-------------------------------------------------------------------
// Objective :  Compute u{n+1} based on explicit euler time integration method
//-------------------------------------------------------------------

class advance_time : public pointwise_rule {
  const_store<float> x,ux,un ;
  const_Map il,ir ;
  const_param<float> dt,nu ;
  store<float> unp1 ;
public:
  advance_time() {
    name_store("x{n}",x) ;
    name_store("ux{n}",ux) ;
    name_store("u{n}",un) ;
    name_store("il{n}",il) ;
    name_store("ir{n}",ir) ;
    name_store("nu{n}",nu) ;
    name_store("dt{n}",dt) ;
    name_store("u{n+1}",unp1) ;
    input("u{n},nu{n},dt{n},(ir{n},il{n})->(ux{n},x{n})") ;
    output("u{n+1}") ;

    constraint("(il{n},ir{n})->x{n}") ; // Assert update for all valid cells
  }
    
  void calculate(Entity e) {
    unp1[e] = un[e] + nu[e]*dt[e]*(ux[ir[e]]-ux[il[e]])/(x[ir[e]]-x[il[e]]) ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&advance_time::calculate) ;
  }
} ;

register_rule<advance_time> register_advance_time ;

//*******************************************************************

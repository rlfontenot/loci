#include <Loci.h>

#include <Loci.h>

//*******************************************************************
//Initializing values at time= 0 using u(x,0) = f(x). Since "xc" is
//used, it will calculate the values are cell-centers.
//*******************************************************************

float f(float x) {
  return 0 ;
}

class initial_condition : public pointwise_rule {
  const_store<float> xc ;
  store<float> u ;
public:
  initial_condition() {
    name_store("xc",xc) ;
    name_store("u{n=0}",u) ;
    input("xc") ;
    output("u{n=0}") ; // u{n=0}<-xc
  }
  void calculate(Entity e) {
    u[e] = f(xc[e]) ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<initial_condition> register_initial_condition ;

//*******************************************************************


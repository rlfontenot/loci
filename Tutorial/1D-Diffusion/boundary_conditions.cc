#include <Loci.h>

//*******************************************************************
// Compute boundary condition at leftmost interface. At left we are
// imposing Neumann Boundary Condition.
//*******************************************************************

class left_bc : public pointwise_rule {
  store<float> ux ;
public:
  left_bc() {
    name_store("ux",ux) ;
    constraint("left_boundary") ;
    output("ux") ; 
  }
  void calculate(Entity e) {
    ux[e] = -1 ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&left_bc::calculate) ;
  }
} ;

register_rule<left_bc> register_left_bc ;

//*******************************************************************
// Compute boundary condition at rightmost interface. At right we are
// imposing Neumann Boundary Condition.
//*******************************************************************

class right_bc : public pointwise_rule {
  const_store<float> u,xc,x ;
  const_Map cl ;
  store<float> ux ;
public:
  right_bc() {
    name_store("u",u) ;
    name_store("cl",cl) ;
    name_store("x",x) ;
    name_store("xc",xc) ;
    name_store("ux",ux) ;
    input("x,cl->(u,xc)") ;
    output("ux") ;
    constraint("right_boundary") ;
  }
  void calculate(Entity e) {
    ux[e] = (u[cl[e]])/(xc[cl[e]]-x[e]) ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&right_bc::calculate) ;
  }
} ;

register_rule<right_bc> register_right_bc ;

//*******************************************************************

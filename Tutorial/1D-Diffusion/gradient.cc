#include <Loci.h>

//*******************************************************************
// Compute gradient of u at all internal interfaces. An internal face
// is a face having left and right cell mapping defined. Since a
// boundary face has only one cell adjacent to it, it willn't calculate
// gradient at those points.
//*******************************************************************

class interface_gradient : public pointwise_rule {
  const_store<float> u ;
  const_Map cl,cr ;
  const_store<float> xc ;
  store<float> ux ;
public:
  interface_gradient() {
    name_store("u",u) ;
    name_store("cl",cl) ;
    name_store("cr",cr) ;
    name_store("xc",xc) ;
    name_store("ux",ux) ;
    input("(cl,cr)->(u,xc)") ;
    output("ux") ;      
  }
  void calculate(Entity e) {
    ux[e] = (u[cl[e]]-u[cr[e]])/(xc[cl[e]]-xc[cr[e]]) ;
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&interface_gradient::calculate) ;
  }
} ;

register_rule<interface_gradient> register_interface_gradient ;

//*******************************************************************

#include <Loci.h>

class calc_gradXdotn : public pointwise_rule {
  const_Map cl,cr ;
  const_store<vector2d<double> > centroid ;
  const_store<double> X ;
  store<double> gradDotN ;
public:
  calc_gradXdotn() {
    name_store("cl",cl) ;
    name_store("cr",cr) ;
    name_store("centroid",centroid) ;
    name_store("X",X) ;
    name_store("gradDotN(X)",gradDotN) ;
    input("(cl,cr)->(centroid,X)") ;
    output("gradDotN(X)") ;
  }
  void calculate(Entity e) {
    gradDotN[e] = (X[cr[e]]-X[cl[e]])/norm(centroid[cr[e]]-centroid[cl[e]]) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<calc_gradXdotn> register_calc_gradXdotn ;

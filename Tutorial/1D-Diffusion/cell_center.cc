#include <Loci.h>

//*******************************************************************
// Compute center of cells: Also make sure to compute the centroid of
// only those cells, where left and right interface entities are
// defined.
//*******************************************************************

class cell_center : public pointwise_rule {
  const_store<float> x ;
  const_Map il,ir ;
  store<float> xc ;
public:
  cell_center() {
    name_store("x",x) ;
    name_store("il",il) ;
    name_store("ir",ir) ;
    name_store("xc",xc) ;
    input("(il,ir)->x") ;
    output("xc") ;         // xc <- (il,ir)->x
  }
  void calculate(Entity e) {
    xc[e] = 0.5*(x[il[e]]+x[ir[e]]);
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&cell_center::calculate) ;
  }
} ;

register_rule<cell_center> register_cell_center ;

//*******************************************************************


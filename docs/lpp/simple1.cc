//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci.h>
#include "TSM_param.h"
using std::cerr ;
using std::endl ;

namespace Loci {

  
  typedef vector3d<real_t> vect3d ;
  typedef tensor3d<real_t> tens3d ;
  typedef real_t real ;
  
  // $type firstOrderCells store<char> 
  // $type geom_cells Constraint
  // $type vol store<Loci::real_t> 
  // $type cl Map
  // $type cr Map
  //comments 1
  class file_simple000_1608568559m193 : public Loci::unit_rule {
    store<char>  firstOrderCells ; 
public:
    file_simple000_1608568559m193() {
       name_store("firstOrderCells",firstOrderCells) ;
       output("firstOrderCells") ;
       constraint("geom_cells") ;
    }
    void calculate(Loci::Entity e) { 
    firstOrderCells[e]= 0 ;
  }    void compute(const Loci::sequence &seq) { 
      do_loop(seq,this) ;
    }
} ;
register_rule<file_simple000_1608568559m193> register_file_simple000_1608568559m193 ;


  
  class file_simple001_1608568559m193 : public Loci::apply_rule< store<char> ,Loci::Maximum<char> >  {
    const_store<Loci::real_t>  vol ; 
    const_Map cl ; 
    const_Map cr ; 
    store<char>  firstOrderCells ; 
public:
    file_simple001_1608568559m193() {
       name_store("firstOrderCells",firstOrderCells) ;
       name_store("vol",vol) ;
       name_store("cl",cl) ;
       name_store("cr",cr) ;
       input("(cl,cr)->vol") ;
       output("(cl,cr)->firstOrderCells") ;
       constraint("(cl,cr)->geom_cells") ;
    }
    void calculate(Loci::Entity e) { 
//comments 2
    if (max (vol[cl[e]],vol[cr[e]]) > 50.*min (vol[cl[e]],vol[cr[e]])) {
      char tmp = 1 ;
      join (firstOrderCells[cl[e]],tmp ) ;
      join (firstOrderCells[cr[e]],tmp ) ;
    }
  }    void compute(const Loci::sequence &seq) { 
      do_loop(seq,this) ;
    }
} ;
register_rule<file_simple001_1608568559m193> register_file_simple001_1608568559m193 ;


  inline real vlimit(real Xcc, real qmin, real qmax, real qdif, real eps2,
		     real ref) {
#ifdef REFLIM
    // old way of preventing div zero errors
    if(fabs(qdif)<=ref*1e-4) return 1.0 ;
#endif
    // delta +
    const real delp = (qdif>0.0)?qmin-Xcc:qmax-Xcc;
    // delta -
    const real delm = -qdif ;
    // numerator of limiter
    const real num = ((delp*delp+eps2)*delm+ 2.0*delm*delm*delp)  ;
    // denominator of limiter
    const real den = (delm*(delp*delp+2.0*delm*delm+delm*delp+eps2)) ;
    // make the limiting case of 0/0 work as expected
    const real e = (den >= 0.0?1.0e-30:-1.0e-30) ;
    return (num+e)/(den+e) ;
  }

  // $type X store<real> 
  // $type limiters(X0) store<Loci::real_t> 
  class file_simple002_1608568559m194 : public Loci::pointwise_rule {
    const_store<real>  X ; 
    store<Loci::real_t>  limitersX ; 
public:
    file_simple002_1608568559m194() {
       name_store("X",X) ;
       name_store("limiters(X)",limitersX) ;
       input("X") ;
       output("limiters(X)") ;
    }
    void calculate(Loci::Entity e) { 
    limitersX[e]= 1.0 ;
  }    void compute(const Loci::sequence &seq) { 
      do_loop(seq,this) ;
    }
} ;
register_rule<file_simple002_1608568559m194> register_file_simple002_1608568559m194 ;
class file_simple002_1608568559m194 : public PointwiseRule {
public:
    file_simple002_1608568559m194() {}
    file_simple002_1608568559m194(int ctxId)
    :PointwiseRule(ctxId)  {}


    class Compute { 
      public: 
      Loci::sequence domain ;
      const real* X ; 
      Loci::real_t* limitersX ; 

      void bind(StoreDB<GI, T> & db) {
       X =  db.X;
       limitersX =  db.limiters(X);
      }

      __host__ __device__
      Loci::sequence getDomain() const {
        return domain ;
      }

      __host__ __device__
      void setDomain(Loci::sequence& dom)  {
       domain = dom ;
      }

      __host__ __device__ 
      void operator()(Loci::Entity e) { 
    limitersX[e]= 1.0 ;
  }} ;
} ;

  
  // $type face2node multiMap
  // $type NGTNodalv3dMax(X0) store<vector3d<real_t> >  
 
  class file_simple003_1608568559m194 : public Loci::apply_rule< store<vector3d<real_t> > ,Loci::Maximum<vector3d<real_t> > >  {
    const_Map cl ; 
    const_store<real>  X ; 
    const_multiMap face2node ; 
    store<vector3d<real_t> >  NGTNodalv3dMaxX ; 
public:
    file_simple003_1608568559m194() {
       name_store("cl",cl) ;
       name_store("X",X) ;
       name_store("face2node",face2node) ;
       name_store("NGTNodalv3dMax(X)",NGTNodalv3dMaxX) ;
       input("cl->X") ;
       output("face2node->NGTNodalv3dMax(X)") ;
    }
    void calculate(Loci::Entity e) { 
    int nsz = face2node[e].size () ;
    for (int i =0;i <nsz ;++i )
      join (NGTNodalv3dMaxX[face2node[e][i ]],X[cl[e]]) ;
  }    void compute(const Loci::sequence &seq) { 
      do_loop(seq,this) ;
    }
} ;
register_rule<file_simple003_1608568559m194> register_file_simple003_1608568559m194 ;
class file_simple003_1608568559m194 : public ApplyRule< store<vector3d<real_t> > ,Loci::Maximum<vector3d<real_t> > >  {
public:
    file_simple003_1608568559m194() {}
    file_simple003_1608568559m194(int ctxId)
    :ApplyRule(ctxId)  {}


    class Compute { 
      public: 
      Loci::sequence domain ;
      const Loci::int_type* cl ; 
      const real* X ; 
      const Loci::int_type* face2node ; 
      const Loci::int_type* face2nodeOffset  ; 
      vector3d<real_t> * NGTNodalv3dMaxX ; 

      void bind(StoreDB<GI, T> & db) {
       cl =  db.cl;
       X =  db.X;
       face2node =  db.face2node;
       NGTNodalv3dMaxX =  db.NGTNodalv3dMax(X);
      }

      __host__ __device__
      Loci::sequence getDomain() const {
        return domain ;
      }

      __host__ __device__
      void setDomain(Loci::sequence& dom)  {
       domain = dom ;
      }

      __host__ __device__ 
      void operator()(Loci::Entity e) { 
    int nsz = face2node[e].size () ;
    for (int i =0;i <nsz ;++i )
      join (NGTNodalv3dMaxX[face2node[e][i ]],X[cl[e]]) ;
  }} ;
} ;


 
}//end of namespace Loci


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
  class file_simple002_1608568559m174 : public PointwiseRule {
public:
    file_simple002_1608568559m174() {}
    file_simple002_1608568559m174(int ctxId)
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
 
  class file_simple003_1608568559m174 : public ApplyRule< store<vector3d<real_t> > ,Loci::Maximum<vector3d<real_t> > >  {
public:
    file_simple003_1608568559m174() {}
    file_simple003_1608568559m174(int ctxId)
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


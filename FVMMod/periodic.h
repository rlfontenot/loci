#ifndef PERIODIC_H
#define PERIODIC_H
#ifdef OLD
#include <Loci.h>
#include <string>
#include "FVMtypes.h"

namespace Loci {

  struct rigid_transform {
    vect3d t1,t2 ;
    tens3d R,Rinv ;
    rigid_transform() {}
    rigid_transform(vect3d center, vect3d v, real angle, vect3d translate) {
      t1 = -1.*center ;
      t2 = center + translate ;
      R.x = vect3d(1,0,0) ;
      R.y = vect3d(0,1,0) ;
      R.z = vect3d(0,0,1) ;
      Rinv = R ;
      if(angle == 0)
        return ;
      real s = sin(angle) ;
      real c = cos(angle) ;
      real C = 1.-c ;
      real x = v.x ;
      real y = v.y ;
      real z = v.z ;
      R.x += vect3d(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      R.y += vect3d(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      R.z += vect3d(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
      s = -s ;
      Rinv.x += vect3d(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      Rinv.y += vect3d(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      Rinv.z += vect3d(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
    }
    vect3d rotate_vec(vect3d v) const {
      return dot(R,v) ;
    }
    tens3d rotate_tensor(const tens3d &t) const {
      return product(R,product(t,Rinv)) ;
    }
    vect3d transform(vect3d v) const {
      return rotate_vec(v+t1)+t2 ;
    }
  } ;
  
  struct periodic_info {
    std::string name ;
    bool master,processed ;
    vect3d center, v, translate ;
    real angle ;
    entitySet bset ;
    Entity bc_num ;
    periodic_info() {
      name = "PERIODIC" ;
      master = false ;
      processed = false ;
      center = vect3d(0.,0.,0.) ;
      v = vect3d(1.,0.,0.) ;
      translate = center ;
      angle = 0 ;
      bset = EMPTY ;
      bc_num = 0 ;
    }
  } ;

  void setup_periodic_bc(std::list<pair<periodic_info,periodic_info> >
                    &periodic_list,fact_db &facts) ;
  template<> struct data_schema_traits<Loci::rigid_transform> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::rigid_transform()) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,t1) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,t2) ;
      LOCI_INSERT_TYPE(ct,Loci::rigid_transform,R) ;
      return DatatypeP(ct) ;
      
    }
  } ;
}
#endif

#endif

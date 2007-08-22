#ifndef LOCI_GRID_READERS_H
#define LOCI_GRID_READERS_H
#include <fact_db.h>
#include <string>

namespace Loci {

  typedef double real_t ;

  //Define struct Area which data members of normal vector and area of the area
  struct Area {
    vector3d<real_t> n ;  //normal vector of the face
    real_t sada ; //area of the face
  } ;

  //Overload ostream and istream (Input/Output) operators for struct Area
  inline std::ostream & operator<<(std::ostream &s, const Area &v)
  {
    s << v.n << ' ' << v.sada << ' ' ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, Area &v)
  {
    s >> v.n >> v.sada  ;
    return s ;
  }

  template<> struct data_schema_traits<Loci::Area> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::Area()) ;
      LOCI_INSERT_TYPE(ct,Loci::Area,n) ;
      LOCI_INSERT_TYPE(ct,Loci::Area,sada) ;
      return DatatypeP(ct) ;
    }
  } ;
  
  
  struct rigid_transform {
    vector3d<real_t> t1,t2 ;
    tensor3d<real_t> R,Rinv ;
    rigid_transform() {}
    rigid_transform(vector3d<real_t> center, vector3d<real_t> v, real_t angle, vector3d<real_t> translate) {
      t1 = -1.*center ;
      t2 = center + translate ;
      R.x = vector3d<real_t>(1,0,0) ;
      R.y = vector3d<real_t>(0,1,0) ;
      R.z = vector3d<real_t>(0,0,1) ;
      Rinv = R ;
      if(angle == 0)
        return ;
      real_t s = sin(angle) ;
      real_t c = cos(angle) ;
      real_t C = 1.-c ;
      real_t x = v.x ;
      real_t y = v.y ;
      real_t z = v.z ;
      R.x += vector3d<real_t>(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      R.y += vector3d<real_t>(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      R.z += vector3d<real_t>(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
      s = -s ;
      Rinv.x += vector3d<real_t>(C*(x*x-1.),C*x*y-z*s,C*x*z+y*s) ;
      Rinv.y += vector3d<real_t>(C*x*y+z*s,C*(y*y-1.),C*y*z-x*s) ;
      Rinv.z += vector3d<real_t>(C*x*z-y*s,C*y*z+x*s,C*(z*z-1.)) ;
    }
    vector3d<real_t> rotate_vec(vector3d<real_t> v) const {
      return dot(R,v) ;
    }
    tensor3d<real_t> rotate_tensor(const tensor3d<real_t> &t) const {
      return product(R,product(t,Rinv)) ;
    }
    vector3d<real_t> transform(vector3d<real_t> v) const {
      return rotate_vec(v+t1)+t2 ;
    }
  } ;
  
  struct periodic_info {
    std::string name ;
    bool master,processed ;
    vector3d<real_t> center, v, translate ;
    real_t angle ;
    entitySet bset ;
    Entity bc_num ;
    periodic_info() {
      name = "PERIODIC" ;
      master = false ;
      processed = false ;
      center = vector3d<real_t>(0.,0.,0.) ;
      v = vector3d<real_t>(1.,0.,0.) ;
      translate = center ;
      angle = 0 ;
      bset = EMPTY ;
      bc_num = 0 ;
    }
  } ;

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

  bool readFVMGrid(fact_db &facts, std::string filename) ;
  bool setupFVMGrid(fact_db &facts, std::string filename) ;
  bool readBCfromVOG(std::string filename,
                     std::vector<std::pair<int,std::string> > &boundary_ids) ;
  void setupBoundaryConditions(fact_db &facts) ;
  void createLowerUpper(fact_db &facts) ;
  void createEdgesPar(fact_db& facts) ;


  extern void writeVOG(std::string filename,store<vector3d<double> > &pos,
                       Map &cl, Map &cr, multiMap &face2node,
                       std::vector<std::pair<int,std::string> > surfaceids) ;
}

#endif

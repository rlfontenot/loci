#ifndef METRICS_H
#define METRICS_H
#include "FVMtypes.h"

namespace Loci {

  
  //Define struct Area which data members of normal vector and area of the area
  struct Area {
    vect3d n ;  //normal vector of the face
    real sada ; //area of the face
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

}

#endif

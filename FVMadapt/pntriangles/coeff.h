#include <Loci.h>
#include <iostream>

typedef Loci::vector3d<double> vect3d;

//Define struct GeomCoeff
namespace Loci{
  struct GeomCoeff {
   
    vect3d b300;
    vect3d b030;
    vect3d b003;
    vect3d b210;
    vect3d b120;
    vect3d b021;
    vect3d b012;
    vect3d b102;
    vect3d b201;
    vect3d b111;
  } ;

  //Overload ostream and istream (Input/Output) operators for struct GeomCoeff
  inline std::ostream & operator<<(std::ostream &s, const GeomCoeff &g)
  {
    s << g.b300 << ' ' << g.b030 << ' ' <<g.b003 << ' ' << g.b210 << ' '
      << g.b120 << ' ' << g.b021 << ' ' <<g.b012 << ' ' << g.b102 << ' '
      << g.b201 << ' ' << g.b111 <<endl;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, GeomCoeff &g)
  {
    s >> g.b300  >> g.b030  >>g.b003  >> g.b210 
      >> g.b120  >> g.b021  >>g.b012  >> g.b102 
      >> g.b201  >> g.b111 ;
    return s ;
  }

  template<> struct data_schema_traits<Loci::GeomCoeff> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::GeomCoeff()) ;
    
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b300) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b030) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b003) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b210) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b120) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b021) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b012) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b102) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b201) ;
      LOCI_INSERT_TYPE(ct,Loci::GeomCoeff,b111) ;
      
      return DatatypeP(ct) ;
    }
  } ;
}

#ifndef LOCI_TYPES_H
#define LOCI_TYPES_H
#include <data_traits.h>
#include <Tools/stream.h>
namespace Loci {

  //-----------STD pair-------------------------------//
  template<class T1,class T2> std::ostream &
    operator<<(std::ostream &s, const std::pair<T1,T2> &v) {
    s<<"["<<v.first<<","<<v.second<<"]";
    return s;
  }

  template<class T1,class T2> std::istream &
    operator>>(std::istream &s, std::pair<T1,T2> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a '[' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a ',' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a ']' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    return s;
  }

  //---------------------vector3d------------------//
  template <class T> 
    struct vector3d {
      T x,y,z ;
      vector3d() {} 
      vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
      vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
    {
      s >> v.x >> v.y >> v.z ;
      return s ;
    }

  template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
  }

  template <class T> inline T dot(const vector3d<T> &v1, const T ra2[]) {
    return v1.x*ra2[0] + v1.y*ra2[1] + v1.z*ra2[2] ;
  }

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
  }

  template <class T> inline T dot(const T ra1[], const vector3d<T> &v2) {
    return ra1[0]*v2.x + ra1[1]*v2.y + ra1[2]*v2.z ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
    return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                       v1.z*ra2[0]-v1.x*ra2[2],
                       v1.x*ra2[1]-v1.y*ra2[0]) ;
  }
  template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
    return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                       ra1[2]*v2.x-ra1[0]*v2.z,
                       ra1[0]*v2.y-ra1[1]*v2.x) ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    target.z += val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    target.z -= val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(float r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(long double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template <class T>
    struct data_schema_traits< vector3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector3d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;
        ct->insert("x",offsetof(vector3d<T>,x),getLociType(t.x)) ;
        ct->insert("y",offsetof(vector3d<T>,y),getLociType(t.y)) ;
        ct->insert("z",offsetof(vector3d<T>,z),getLociType(t.z)) ;
        return DatatypeP(ct) ;
      }
    };
  
  //---------------------vector2d------------------//
  template <class T> 
    struct vector2d {
      T x,y ;
      vector2d() {} 
      vector2d(T xx,T yy) : x(xx),y(yy) {}
      vector2d(const vector2d &v) {x=v.x;y=v.y;}
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector2d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector2d<T> &v)
    {
      s >> v.x >> v.y ;
      return s ;
    }

  template <class T> inline T dot(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y ;
  }

  template <class T> inline T dot(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[0] + v1.y*ra2[1] ;
  }

  template <class T> inline T norm(const vector2d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y) ;
  }

  template <class T> inline T dot(const T ra1[], const vector2d<T> &v2) {
    return ra1[0]*v2.x + ra1[1]*v2.y ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.y-v1.y*v2.x ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[1]-v1.y*ra2[0] ;
  }

  template<class T> inline T cross(const T ra1[], const vector2d<T> &v2) {
    return ra1[0]*v2.y-ra1[1]*v2.x ;
  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+=(vector2d<T> &target, const vector2d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator-=(vector2d<T> &target, const vector2d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x+v2.x,v1.y+v2.y) ;
  }

  template<class T> inline vector2d<T> operator-(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x-v2.x,v1.y-v2.y) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, float r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(float r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, float r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, double r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(double r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, double r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, long double r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(long double r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, long double r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  
  template <class T>
    struct  data_schema_traits< vector2d<T> > {
    public:
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector2d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;
        ct->insert("x",offsetof(vector2d<T>,x),getLociType(t.x)) ;
        ct->insert("y",offsetof(vector2d<T>,y),getLociType(t.y)) ;
        return DatatypeP(ct) ;
      }
    };
  
  
  
  //---------------------Array----------------------//
  template <class T,unsigned int n> class Array {
    T x[n] ;
  public:
    typedef T * iterator ;
    typedef const T * const_iterator ;
    
    Array() {} ;
    Array(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] = v.x[i] ; } 
    Array<T,n> &operator=(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 

    Array<T,n> &operator +=(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
    Array<T,n> &operator -=(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
    Array<T,n> &operator *=(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
    Array<T,n> &operator /=(const Array<T,n> &v)
    { for(unsigned int i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }

    T &operator[](unsigned int indx) { return x[indx]; }
    const T &operator[](unsigned int indx) const { return x[indx] ; }

    iterator begin() { return &x[0] ; }
    iterator end() { return begin()+n ; }
    const_iterator begin() const { return &x[0] ; }
    const_iterator end() const { return begin()+n ; }

    unsigned int size() const  { return n ; }
  } ;

  template <class T,unsigned int n> inline std::ostream &
    operator<<(std::ostream &s, const Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,unsigned int n> inline std::istream &
    operator>>(std::istream &s, Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }

  
  template <class T,unsigned int n> 
  class data_schema_traits< Array<T,n> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = n ;
      return new ArrayType(getLociType(T()),sizeof(Array<T,n>),1,&dim) ;
    }
  };
}
#endif

#ifndef TYPES_H
#define TYPES_H

#include <Tools/stream.h>
namespace Loci {

  //------------STL vector----------------------------//
  template<class T> ostream & operator<<(ostream &s, const std::vector<T> &v)
    {
      std::vector<T>::const_iterator ii=v.begin();
      for(ii=v.begin();ii!=v.end();ii++){
	s<<*ii<<endl;
      }
      return s;
    }
  template<class T> istream & operator>>(istream &s, std::vector<T> &v)
    {
      std::vector<T>::iterator ii;
      for(ii=v.begin();ii!=v.end();ii++){
	s>>*ii;
      }
      return s;
    }

  //-----------STD pair-------------------------------//
 template<class T1,class T2> ostream & operator<<(ostream &s, const std::pair<T1,T2> &v)
    {
     s<<"["<<v.first<<","<<v.second<<"]";
      return s;
    }
  template<class T1,class T2> istream & operator>>(istream &s, std::pair<T1,T2> &i)
    {
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
  
   template <class T> inline ostream & operator<<(ostream &s, const vector3d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
      return s ;
    }

 template <class T> inline istream &operator>>(istream &s, vector3d<T> &v)
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

//---------------------vec----------------------//
template <class T,int n> class vec {
    T x[n] ;
  public:
    T &operator[](int indx) { return x[indx]; }
    const T &operator[](int indx) const { return x[indx] ; }
    vec<T,n> operator +=(const vec<T,n> &v)
      { for(int i=0;i<n;++i) x[i] += v.x[i] ; return v ; }
    vec<T,n> operator -=(const vec<T,n> &v)
      { for(int i=0;i<n;++i) x[i] -= v.x[i] ; return v ; }
    vec<T,n> operator *=(const vec<T,n> &v)
      { for(int i=0;i<n;++i) x[i] *= v.x[i] ; return v ; }
    vec<T,n> operator /=(const vec<T,n> &v)
      { for(int i=0;i<n;++i) x[i] /= v.x[i] ; return v ; }
} ;

template <class T,int n> inline ostream & operator<<(ostream &s, const vec<T,n> &v)
{
    for(int i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
}

template <class T,int n> inline istream &operator>>(istream &s, vec<T,n> &v)
{
    for(int i=0;i<n;++i)
      s >> v[i] ;
    return s ;
}

}
#endif

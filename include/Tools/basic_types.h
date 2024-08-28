//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
#ifndef LOCI_BASIC_TYPES_H
#define LOCI_BASIC_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/autodiff.h>


namespace Loci {

#ifdef USE_AUTODIFF
#ifdef AUTODIFF2ND
  typedef Loci::FAD2d real_t ;
#else
#ifdef MULTIFAD
  typedef Loci::MFADd real_t ;
#else
  typedef Loci::FADd real_t ;
#endif
#endif
#else
  // No autodiff
  typedef double real_t ;
#endif
  //---------------------Array----------------------//
  template <class T,size_t n> class Array {
    T x[n] ;
  public:
    typedef T * iterator ;
    typedef const T * const_iterator ;
    
    //    Array() {} ;
    //    Array(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; } 
    //    Array<T,n> &operator=(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 

    Array<T,n> &operator +=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
    Array<T,n> &operator -=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
    Array<T,n> &operator *=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
    Array<T,n> &operator /=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }

    T &operator[](size_t indx) { return x[indx]; }
    const T &operator[](size_t indx) const { return x[indx] ; }
    T &operator[](int indx) { return x[indx]; }
    const T &operator[](int indx) const { return x[indx] ; }
    T &operator[](unsigned char indx) { return x[indx]; }
    const T &operator[](unsigned char indx) const { return x[indx] ; }
    T &operator[](unsigned int indx) { return x[indx]; }
    const T &operator[](unsigned int indx) const { return x[indx] ; }

    iterator begin() { return &x[0] ; }
    iterator end() { return begin()+n ; }
    const_iterator begin() const { return &x[0] ; }
    const_iterator end() const { return begin()+n ; }

    size_t size() const  { return n ; }
    
  } ;

  //---------------------vector3d------------------//
  template <class T> 
    struct vector3d {
      T x,y,z ;
      vector3d() {} 
      vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
      vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
      template <class S> vector3d(const vector3d<S> &v) {x=v.x;y=v.y;z=v.z;}
      template <class S> vector3d(const Array<S,3> &a) {x=a[0];y=a[1];z=a[2];}
      template <class S> operator Array<S,3>() {
	Array<S,3> a ;
	a[0] = x ;
	a[1] = y ;
	a[2] = z ;
	return a;
      }
    T &operator[](int i) {
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      case 2:
        return z ;
      default:
        return z ;
      }
    }
    const T &operator[](int i) const {
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      case 2:
        return z ;
      default:
        return z ;
      }
    }
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

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
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

  template<class T,class S> inline vector3d<T> &operator*=(vector3d<T> &target, S val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T,class S> inline vector3d<T> &operator/=(vector3d<T> &target, S val) {
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

  template<class T,class S> inline vector3d<T> operator*(const vector3d<T> &v1, S r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }
  template<class T,class S> inline vector3d<T> operator*(S r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T,class S> inline vector3d<T> operator/(const vector3d<T> &v1, S r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T>  struct tensor3d : public vector3d<vector3d< T > > {
    tensor3d() {}
    tensor3d(vector3d<T> xx,vector3d<T> yy, vector3d<T> zz)
      : vector3d<vector3d< T> > (xx,yy,zz) {}
    tensor3d(const tensor3d &v) : vector3d<vector3d< T> >(v) {}
  } ;

  template<class T> inline vector3d<T> dot(const tensor3d<T> &t,
                                           const vector3d<T> &v) {
    return vector3d<T>(dot(t.x,v),dot(t.y,v),dot(t.z,v)) ;
  }

  template<class T> inline tensor3d<T> product(const tensor3d<T> &t1,
                                               const tensor3d<T> &t2) {
    tensor3d<T> temp ;
    temp.x.x = t1.x.x*t2.x.x+t1.x.y*t2.y.x+t1.x.z*t2.z.x ;
    temp.y.x = t1.y.x*t2.x.x+t1.y.y*t2.y.x+t1.y.z*t2.z.x ;
    temp.z.x = t1.z.x*t2.x.x+t1.z.y*t2.y.x+t1.z.z*t2.z.x ;

    temp.x.y = t1.x.x*t2.x.y+t1.x.y*t2.y.y+t1.x.z*t2.z.y ;
    temp.y.y = t1.y.x*t2.x.y+t1.y.y*t2.y.y+t1.y.z*t2.z.y ;
    temp.z.y = t1.z.x*t2.x.y+t1.z.y*t2.y.y+t1.z.z*t2.z.y ;

    temp.x.z = t1.x.x*t2.x.z+t1.x.y*t2.y.z+t1.x.z*t2.z.z ;
    temp.y.z = t1.y.x*t2.x.z+t1.y.y*t2.y.z+t1.y.z*t2.z.z ;
    temp.z.z = t1.z.x*t2.x.z+t1.z.y*t2.y.z+t1.z.z*t2.z.z ;

    return temp ;
  }

  inline vector3d<float> realToFloat(vector3d<double> v) { return vector3d<float>(float(v.x),float(v.y),float(v.z)); }
  inline vector3d<double> realToDouble(vector3d<double> v) { return v ; }

  inline vector3d<float> realToFloat(vector3d<MFADd> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<MFADd> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }

  inline vector3d<float> realToFloat(vector3d<FADd> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<FADd> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }

  inline vector3d<float> realToFloat(vector3d<FAD2d> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<FAD2d> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }

  //---------------------vector2d------------------//
  template <class T> 
    struct vector2d {
      T x,y ;
      vector2d() {} 
      vector2d(T xx,T yy) : x(xx),y(yy) {}
      vector2d(const vector2d &v) {x=v.x;y=v.y;}
    T &operator[](int i) {
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      default:
        return y ;
      }
    }
    const T &operator[](int i) const {
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      default:
        return y ;
      }
    }
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

  template<class T> inline T cross(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.y-v1.y*v2.x ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[1]-v1.y*ra2[0] ;
  }


  template<class T,class S> inline vector2d<T> &operator*=(vector2d<T> &target, S val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T,class S> inline vector2d<T> &operator/=(vector2d<T> &target, S val) {
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

  template<class T,class S> inline vector2d<T> operator*(const vector2d<T> &v1, S r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T,class S> inline vector2d<T> operator*(S r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T,class S> inline vector2d<T> operator/(const vector2d<T> &v1, S r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  // For allocating temporary arrays of small size
  const int tmp_array_internal_SIZE=25 ;
  template <class T> class tmp_array {
    int sz ;
    T data[tmp_array_internal_SIZE] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > tmp_array_internal_SIZE)
        p = new T[sz] ;
    }
    void free() {
      if(sz > tmp_array_internal_SIZE)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { return p[i] ; }
    T & operator[](int i) const { return p[i] ; }
    T & operator[](size_t i) { return p[i] ; }
    T & operator[](size_t i) const { return p[i] ; }
    T & operator[](unsigned int i) { return p[i] ; }
    T & operator[](unsigned int i) const { return p[i] ; }
    T & operator[](unsigned char i) { return p[i] ; }
    T & operator[](unsigned char i) const { return p[i] ; }
    T & operator[](unsigned short i) { return p[i] ; }
    T & operator[](unsigned short i) const { return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;
}

#endif

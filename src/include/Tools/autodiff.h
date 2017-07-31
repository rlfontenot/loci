//#############################################################################
//#
//# Copyright 2017, Mississippi State University
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
#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <iostream> 
#include <Tools/tools.h>
#include <algorithm>

namespace Loci {
  
  class FADd {
  public:

    double value;
    double grad;
#ifdef LIVING_DANGEROUS
    /* MPGCOMMENT [05-12-2017 13:50] ---> type-cast */
    explicit operator char() {return (char)value;}
    explicit operator bool() {return (bool)value;}
    explicit operator int() {return (int)value;}
    explicit operator float() {return (float)value;}
    operator double() {return (double)value;}
    explicit operator FADd() {return FADd(value,grad);}
    explicit operator long double() const {return (long double)value;}

    //    operator double*()  {return &value;}/* MPGCOMMENT [05-15-2017 10:51] ---> only one used, others set for completeness */
    //    operator const double*()  {return (const double*)&value;}/* MPGCOMMENT [05-15-2017 10:51] ---> only one used, others set for completeness */
#endif
    
    /* MPGCOMMENT [05-12-2017 13:49] ---> constructor */
#define gradcheck() //if(grad!=0) debugger_() ; //fatal(grad!=0.0) ;
    template< class T1>
      FADd(T1 a0): value(a0), grad(0.0){ gradcheck();}/* MPGCOMMENT [05-15-2017 14:46] ---> used in chem, chem::Temperature */
  FADd(int a0         ,int b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,float b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,double b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,long double b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}

    //    explicit FADd(bool &u) : value((u)?1.0:0.0), grad(0.0) { gradcheck();}
    FADd(const FADd &u) : value(u.value), grad(u.grad) { gradcheck();}
    //  explicit FADd(int a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(float a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(double a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(long double a0)        : value(a0), grad(0.0) { gradcheck();}
  FADd() : value(0.0), grad(0.0) { gradcheck();}


    //    friend void swap(double& a, FADd& b) {
    //      using std::swap;
    //      swap(a, b.value);
    //    }

    /**************************************************************/
    /* MPGCOMMENT [05-12-2017 10:38] ---> FADd --- FADd OPERATORS */

    /* MPGCOMMENT [05-12-2017 10:35] ---> ASSIGNMENT */ 
    /* template<class T> Loci::const_store<Loci::vector3d<FADd> >& operator =(Loci::const_store<Loci::vector3d<T> > &u) { */
    /*   value = u; */
    /*   grad = 0.0; */
    /*   return *this; */
    /* } */

    //    template<class T> FADd& operator =(const T &u) {
    //      value = u;
    //      grad = 0.0;
    //      return *this;
    //    }
    FADd& operator =(const FADd &u) {
      value = u.value;
      grad = u.grad;
      return *this;
    }
    FADd& operator =(const int &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const float &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const double &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const long double &u) {
      value = u;
      grad = 0.0;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:35] ---> ADDITION */
    //    template<class T> FADd operator +(const T &v) const{
    //      return FADd(value+(double)v, grad);
    //    }
    FADd operator +(const FADd &v) const{
      return FADd(value+v.value, grad+v.grad);
    }
    FADd operator +() const{
      return FADd(value, grad);
    }
    FADd operator +(const int &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const float &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const double &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const long double &v) const{
      return FADd(value+(long double)v, grad);
    }
    //    template<class T> FADd& operator +=(const T &u) {
    //      value+=(double)u;
    //      return *this;
    //    }
    FADd& operator +=(const FADd &u) {
      value+=u.value;
      grad+=u.grad;
      return *this;
    }
    FADd& operator +=(const int &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const float &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const double &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const long double &u) {
      value+=u;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:35] ---> SUBTRACTION */
    //    template<class T> FADd operator -(const T &v) const{
    //      return FADd(value-(double)v, grad);
    //    }
    FADd operator -(const FADd &v) const{
      return FADd(value-v.value, grad-v.grad);
    }
    FADd operator -() const {
      return FADd(-value, -grad);
    }
    FADd operator -(const int &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const float &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const double &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const long double &v) const{
      return FADd(value-(long double)v, grad);
    }
    //    template<class T> FADd& operator -=(const T &u) {
    //      value-=(double)u;
    //      return *this;
    //    }
    FADd& operator -=(const FADd &u) {
      value-=u.value;
      grad-=u.grad;
      return *this;
    }
    FADd& operator -=(const int &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const float &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const double &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const long double &u) {
      value-=u;
      return *this;
    }

    /* MPGCOMMENT [05-12-2017 10:36] ---> MULTIPLICATION */
    //    template<class T> FADd operator *(const T &v) const{
    //      return FADd(value*(double)v, grad);
    //    }
    FADd operator *(const FADd &v) const {
      return FADd(value*v.value,  grad*v.value + v.grad*value);
    }
    FADd operator *(const int &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const float &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const double &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const long double &v) const{
      return FADd(value*(long double)v,  grad*(long double)v);
    }
    //    template<class T> FADd& operator *=(const T &u) {
    //      value*=(double)u;
    //      return *this;
    //    }
    FADd& operator *=(const FADd &u) {
      grad = grad*u.value + u.grad * value;
      value*=u.value;
      return *this;
    }
    FADd& operator *=(const int &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const float &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const double &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const long double &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }

    /* MPGCOMMENT [05-12-2017 10:36] ---> DIVISION */
    //    template<class T> FADd operator /(const T &v) const{
    //      return FADd(value/(double)v, grad);
    //    }
    FADd operator /(const FADd &v) const{
      return FADd(value/v.value, (grad*v.value - value*v.grad)/v.value/v.value);
    }
    FADd operator /(const int &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const float &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const double &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const long double &v) const{
      return FADd(value/v, (grad/(long double)v));
    }
    //    template<class T> FADd& operator /=(const T &u) {
    //      value/=(double)u;
    //      return *this;
    //    }
    FADd& operator /=(const FADd &u) {
      grad = (grad*u.value - value * u.grad) / u.value / u.value;
      value/=u.value;
      return *this;
    }
    FADd& operator /=(const int &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const float &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const double &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const long double &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:37] ---> COMPARISON OPERATORS */
    //    template<class T> bool operator ==(const T &u) const {
    //      return ((value==(double)u)?true:false);
    //    }
    //    template<class T> bool operator !=(const T &u) const {
    //      return ((value!=(double)u)?true:false);
    //    }
    //    template<class T> bool operator >(const T &u) const {
      //      return ((value>(double)u)?true:false);
    //    }
    //    template<class T> bool operator <(const T &u) const {
    //      return ((value<(double)u)?true:false);
    //    }
    //    template<class T> bool operator >=(const T &u) const {
    //      return ((value>=(double)u)?true:false);
    //    }
    //    template<class T> bool operator <=(const T &u) const {
    //      return ((value<=(double)u)?true:false);
    //    }
    bool operator ==(const FADd &u)const  {
      return ((value==u.value)?true:false);
    }
    bool operator !=(const FADd &u) const {
      return ((value!=u.value)?true:false);
    }
    bool operator >(const FADd &u) const {
      return ((value>u.value)?true:false);
    }
    bool operator <(const FADd &u) const {
      return ((value<u.value)?true:false);
    }
    bool operator >=(const FADd &u) const {
      return ((value>=u.value)?true:false);
    }
    bool operator <=(const FADd &u) const {
      return ((value<=u.value)?true:false);
    }
    bool operator ==(const int &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const int &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const int &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const int &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const int &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const int &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const float &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const float &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const float &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const float &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const float &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const float &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const double &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const long double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const long double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const long double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const long double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const long double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const long double &u) const {
      return ((value<=u)?true:false);
    }



    /* Array<T,n> &operator +=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; } */
    /* Array<T,n> &operator -=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; } */
    /* Array<T,n> &operator *=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; } */
    /* Array<T,n> &operator /=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ;} */


  };

  inline std::ostream &operator <<(std::ostream& stream, const FADd &u) {
    stream << u.value << '^' << u.grad ;
    return stream;
  }
  inline std::istream &operator >> (std::istream& stream, FADd &u) {
    u.grad = 0 ;
    stream >> u.value;
    if(stream.peek() == '^') {
      stream.get() ;
      stream >> u.grad ;
    }
    return stream;
  }
#ifdef NO_CMATH
  inline FADd ceil(const FADd u) {
    return FADd(::ceil(u.value),0.0) ;
  }
  inline FADd floor(const FADd u) {
    return FADd(::floor(u.value),0.0) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(::sin(u.value), u.grad*::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(::cos(u.value), -u.grad*::sin(u.value));
  }
  inline FADd exp(const FADd u) {
    return FADd(::exp(u.value), u.grad*::exp(u.value));
  }
  inline FADd log(const FADd u) {
    return FADd(::log(u.value), u.grad/u.value);
  }
  inline FADd log10(const FADd u) {
    return FADd(::log(u.value), u.grad/u.value)/::log(10.0);
  }

  inline FADd fabs(FADd u) {
    return FADd(::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  inline FADd abs(FADd u) {
    return FADd(::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FADd pow(const FADd u, const int k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const float k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const double k) {
    return FADd(::pow(u.value, k), 
		k * ::pow(u.value, k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const long double k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value,  (double)k-1.0)*u.grad );
  }
  inline FADd sqrt(const FADd u) {
    double su = ::sqrt(u.value) ;
    return FADd(su,0.5*u.grad/max(su,1e-30)) ;
  }
  inline FADd pow(const FADd k, const FADd u) {
    double kpu = ::pow(k.value,u.value) ;
    return FADd(kpu,::pow(k.value,u.value-1.)*k.grad*u.value +
		kpu*::log(k.value)*u.grad) ;
  }
  inline FADd pow(const int k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(double(k))*u.grad) ;
  }
  inline FADd pow(const float k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(double(k))*u.grad) ;
  }
  inline FADd pow(const double k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(k)*u.grad) ;
  }
  inline FADd pow(const long double k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(k)*u.grad) ;
  }


  inline FADd sinh(const FADd u) {
    return FADd(sinh(u.value),
		0.5*u.grad*(exp(u.value) + exp(-u.value))) ;
  }
  inline FADd cosh(const FADd u) {
    return FADd(cosh(u.value),
		0.5*u.grad*(exp(u.value) - exp(-u.value))) ;
  }
  inline FADd tanh(const FADd u) {
    double ex = exp(min(u.value,350.0)) ;
    double exm = exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    return FADd(tanh(u.value),u.grad*(1.-dex*dex/(sex*sex))) ;
  }
  inline FADd asin(const FADd u) {
    return FADd(::asin(u.value), u.grad/::sqrt(1.0-u.value*u.value) );
  }
  inline FADd acos(const FADd u) {
    return FADd(::acos(u.value), -u.grad/::sqrt(1.0-u.value*u.value) );
  }
  inline FADd atan(const FADd u) {
    return FADd(::atan(u.value), u.grad/(1.0+u.value*u.value) );
  }
  // This will not work in general
  inline FADd atan2(const FADd u, const FADd v) {  
    return atan(u/v);
  }

  inline FADd asinh(const FADd u) {
    return FADd(::asinh(u.value), u.grad/::sqrt(1.0+u.value*u.value) );
  }
  inline FADd acosh(const FADd u) {
    return FADd(::acosh(u.value), u.grad/::sqrt(-1.0+u.value*u.value) );
  }
  inline FADd atanh(const FADd u) {
    return FADd(::atanh(u.value), u.grad/(1.0-u.value*u.value) );
  }


#else
  inline FADd ceil(const FADd u) {
    return FADd(std::ceil(u.value),0.0) ;
  }
  inline FADd floor(const FADd u) {
    return FADd(std::floor(u.value),0.0) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(std::sin(u.value), u.grad*std::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(std::cos(u.value), -u.grad*std::sin(u.value));
  }
  inline FADd exp(const FADd u) {
    return FADd(std::exp(u.value), u.grad*std::exp(u.value));
  }
  inline FADd log(const FADd u) {
    return FADd(std::log(u.value), u.grad/u.value);
  }
  inline FADd log10(const FADd u) {
    return FADd(std::log(u.value), u.grad/u.value)/std::log(10.0);
  }
  inline FADd abs(const FADd u) {
    return FADd(std::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  inline FADd fabs(const FADd u) {
    return FADd(std::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FADd pow(const FADd u, const int k) {
    return FADd(std::pow(u.value, k), 
		(double)k * std::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const float k) {
    return FADd(std::pow(u.value, k), 
		(double)k * std::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const double k) {
    return FADd(std::pow(u.value, k), 
		k * std::pow(u.value, k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const long double k) {
    return FADd(std::pow(u.value,  k), 
		(double)k * std::pow(u.value,  (double)k-1.0)*u.grad );
  }
  inline FADd sqrt(const FADd u) {
    double su = std::sqrt(u.value) ;
    return FADd(su,0.5*u.grad/max(su,1e-30)) ;
  }
  inline FADd pow(const FADd k, const FADd u) {
    double kpu = std::pow(k.value,u.value) ;
    return FADd(kpu,std::pow(k.value,u.value-1.)*k.grad*u.value +
		kpu*std::log(k.value)*u.grad) ;
  }
  inline FADd pow(const int k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(double(k))*u.grad) ;
  }
  inline FADd pow(const float k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(double(k))*u.grad) ;
  }
  inline FADd pow(const double k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(k)*u.grad) ;
  }
  inline FADd pow(const long double k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(k)*u.grad) ;
  }


  inline FADd sinh(const FADd u) {
    return FADd(std::sinh(u.value),
		0.5*u.grad*(std::exp(u.value) + std::exp(-u.value))) ;
  }
  inline FADd cosh(const FADd u) {
    return FADd(std::cosh(u.value),
		0.5*u.grad*(std::exp(u.value) - std::exp(-u.value))) ;
  }
  inline FADd tanh(const FADd u) {
    double ex = std::exp(min(u.value,350.0)) ;
    double exm = std::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    return FADd(std::tanh(u.value),
		u.grad*(1.-dex*dex/(sex*sex))) ;
  }

  inline FADd asin(const FADd u) {
    return FADd(std::asin(u.value), u.grad/std::sqrt(1.0-u.value*u.value) );
  }
  inline FADd acos(const FADd u) {
    return FADd(std::acos(u.value), -u.grad/std::sqrt(1.0-u.value*u.value) );
  }
  inline FADd atan(const FADd u) {
    return FADd(std::atan(u.value), u.grad/(1.0+u.value*u.value) );
  }
  // This will not work in general
  inline FADd atan2(const FADd u, const FADd v) {  
    return atan(u/v);
  }

  inline FADd asinh(const FADd u) {
    return FADd(::asinh(u.value), u.grad/std::sqrt(1.0+u.value*u.value) );
  }
  inline FADd acosh(const FADd u) {
    return FADd(::acosh(u.value), u.grad/std::sqrt(-1.0+u.value*u.value) );
  }
  inline FADd atanh(const FADd u) {
    return FADd(::atanh(u.value), u.grad/(1.0-u.value*u.value) );
  }


#endif
  

  inline FADd operator +(const double u,const FADd v)
  { return FADd(u+v.value, v.grad); }
  inline FADd operator -(const double u,const FADd v) 
  { return FADd(u-v.value, -v.grad); }
  inline FADd operator *(const double u,const FADd v) 
  { return FADd(u*v.value,  v.grad*u); }
  inline FADd operator /(const double &u,const FADd &v) 
  { return FADd(u/v.value, -u*v.grad/v.value/v.value); }

  inline float realToFloat(double v) { return float(v) ; }
  inline float realToFloat(FADd v) { return float(v.value) ; }
  inline double realToDouble(double v) { return v ; }
  inline double realToDouble(FADd v) { return v.value ; }

  inline int signbit(FADd v) { return signbit(v.value) ; }
}

#endif

//#############################################################################
//#
//# Copyright 2017-2019, Mississippi State University
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
  
  //---------------------- second derivative type

  class FAD2d {
  public:

    double value;
    double grad;
    double grad2;
    /* MPGCOMMENT [05-12-2017 13:49] ---> constructor */
#define gradcheck() //if(grad!=0) debugger_() ; //fatal(grad!=0.0) ;
    template< class T1>
      FAD2d(T1 a0): value(a0), grad(0.0),grad2(0.0) { gradcheck();}
  FAD2d(double a0, double b0, double c0): value(a0), grad(b0), grad2(c0) { gradcheck();}
  FAD2d(const FAD2d &u) : value(u.value), grad(u.grad), grad2(u.grad2) { gradcheck();}

  FAD2d() : value(0.0), grad(0.0),grad2(0.0) { gradcheck();}

    FAD2d& operator =(const FAD2d &u) {
      value = u.value;
      grad = u.grad;
      grad2 = u.grad2 ;
      return *this;
    }
    FAD2d& operator =(const int &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0 ;
      return *this;
    }
    FAD2d& operator =(const float &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }
    FAD2d& operator =(const double &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }
    FAD2d& operator =(const long double &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }


    FAD2d operator +(const FAD2d &v) const{
      return FAD2d(value+v.value, grad+v.grad, grad2+v.grad2);
    }
    FAD2d operator +() const{
      return FAD2d(value, grad, grad2);
    }
    FAD2d operator +(const int &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const float &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const double &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const long double &v) const{
      return FAD2d(double(value+v), grad, grad2);
    }

    FAD2d& operator +=(const FAD2d &u) {
      value+=u.value;
      grad+=u.grad;
      grad2+=u.grad2 ;
      return *this;
    }
    FAD2d& operator +=(const int &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const float &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const double &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const long double &u) {
      value+=u;
      return *this;
    }

    FAD2d operator -(const FAD2d &v) const{
      return FAD2d(value-v.value, grad-v.grad,grad2-v.grad2);
    }
    FAD2d operator -() const {
      return FAD2d(-value, -grad, -grad2);
    }
    FAD2d operator -(const int &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const float &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const double &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const long double &v) const{
      return FAD2d(double(value-v), grad, grad2);
    }

    FAD2d& operator -=(const FAD2d &u) {
      value-=u.value;
      grad-=u.grad;
      grad2 -= u.grad2 ;
      return *this;
    }
    FAD2d& operator -=(const int &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const float &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const double &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const long double &u) {
      value-=u;
      return *this;
    }

    FAD2d operator *(const FAD2d &v) const {
      return FAD2d(value*v.value,  grad*v.value + v.grad*value, 
		   grad2*v.value + value*v.grad2 + 2.*grad*v.grad);
    }
    FAD2d operator *(const int &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const float &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const double &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const long double &v) const{
      return FAD2d(double(value*v),  double(grad*v), double(grad2*v));
    }
    //    template<class T> FAD2d& operator *=(const T &u) {
    //      value*=(double)u;
    //      return *this;
    //    }
    FAD2d& operator *=(const FAD2d &u) {
      grad2 = grad2*u.value + value*u.grad2 + 2.*grad*u.grad ;
      grad = grad*u.value + u.grad * value;
      value*=u.value;
      return *this;
    }
    FAD2d& operator *=(const int &u) {
      value*=u;
      grad *=u;
      grad2 *= u ;
      return *this;
    }
    FAD2d& operator *=(const float &u) {
      value*=u;
      grad *=u;
      grad2 *= u;
      return *this;
    }
    FAD2d& operator *=(const double &u) {
      value*=u;
      grad *=u;
      grad2 *= u ;
      return *this;
    }
    FAD2d& operator *=(const long double &u) {
      value*=u;
      grad *=u;
      grad2 *= u;
      return *this;
    }

    FAD2d operator /(const FAD2d &v) const{
      double d = v.value ;
      double d2 = d*d ;
      double d3 = d*d2 ;
      return FAD2d(value/d, 
		   (grad*v.value - value*v.grad)/d2,
		   grad2/d - grad*v.grad/d2 +
		   2.*v.grad*v.grad*value/d3 -
		   (v.grad*grad+v.grad2*value)/d2) ;
    }
    FAD2d operator /(const int &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const float &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const double &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const long double &v) const{
      return FAD2d(double(value/v), double(grad/v), double(grad2/v));
    }

    FAD2d& operator /=(const FAD2d &u) {
      double d = u.value ;
      double d2 = d*d ;
      double d3 = d*d2 ;
      grad2 = (grad2/d - grad*u.grad/d2 +2.*u.grad*u.grad*value/d3 -
	       (u.grad*grad+u.grad2*value)/d2) ;
      grad = (grad*u.value - value * u.grad) / d2 ;
      value/=u.value;
      return *this;
    }
    FAD2d& operator /=(const int &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }
    FAD2d& operator /=(const float &u) {
      value/=u;
      grad /=u;
      grad2 /=u ;
      return *this;
    }
    FAD2d& operator /=(const double &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }
    FAD2d& operator /=(const long double &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }


    bool operator ==(const FAD2d &u)const  {
      return ((value==u.value)?true:false);
    }
    bool operator !=(const FAD2d &u) const {
      return ((value!=u.value)?true:false);
    }
    bool operator >(const FAD2d &u) const {
      return ((value>u.value)?true:false);
    }
    bool operator <(const FAD2d &u) const {
      return ((value<u.value)?true:false);
    }
    bool operator >=(const FAD2d &u) const {
      return ((value>=u.value)?true:false);
    }
    bool operator <=(const FAD2d &u) const {
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

  };

  inline std::ostream &operator <<(std::ostream& stream, const FAD2d &u) {
    stream << u.value << '^' << u.grad << "^^" << u.grad2 ;
    return stream;
  }
  inline std::istream &operator >> (std::istream& stream, FAD2d &u) {
    u.grad2 = 0 ;
    u.grad = 0 ;
    stream >> u.value;
    if(stream.peek() == '^') {
      stream.get() ;
      stream >> u.grad ;
      if(stream.peek() == '^') {
	stream.get() ;
	if(stream.peek() == '^') 
	  stream.get() ;
	stream >> u.grad2 ;
      }
    }

    return stream;
  }

#ifdef NO_CMATH
  inline double ceil(const FAD2d u) {
    return ::ceil(u.value) ;
  }
  inline double floor(const FAD2d u) {
    return ::floor(u.value) ;
  }
  inline FAD2d sin(const FAD2d u) {
    double su = ::sin(u.value) ;
    double cu = ::cos(u.value) ;
    return FAD2d(su, u.grad*cu, cu*u.grad2 - u.grad*u.grad*su);
  }
  inline FAD2d cos(const FAD2d u) {
    double su = ::sin(u.value) ;
    double cu = ::cos(u.value) ;
    return FAD2d(cu, -u.grad*su, -su*u.grad2 - u.grad*u.grad*cu);
  }
  inline FAD2d tan(const FAD2d u) {
    double tanu = ::tan(u.value) ;
    double secu = 1./::cos(u.value) ;
    double secu2 = secu*secu ;
    return FAD2d(tanu, u.grad*secu2,
		 u.grad2*secu2+2.*u.grad*u.grad*secu2*tanu) ;
  }


  inline FAD2d exp(const FAD2d u) {
    double eu = ::exp(u.value) ;
    return FAD2d(eu, u.grad*eu, eu*u.grad2+u.grad*u.grad*eu);
  }
  inline FAD2d log(const FAD2d u) {
    return FAD2d(::log(u.value), u.grad/u.value, 
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value));
  }
  inline FAD2d log10(const FAD2d u) {
    return FAD2d(::log(u.value), u.grad/u.value,
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value))/::log(10.0);
  }

  inline FAD2d fabs(FAD2d u) {
    return FAD2d(::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  inline FAD2d abs(FAD2d u) {
    return FAD2d(::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FAD2d pow(const FAD2d u, const int k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const float k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const double k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const long double k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d sqrt(const FAD2d u) {
    double su = ::sqrt(u.value) ;
    return FAD2d(su,0.5*u.grad/max(su,1e-30),
		 0.5*u.grad2/(max(su,1e-30)) - 0.25*u.grad*u.grad/(max(su*su*su,1e-30))) ;
  }
  inline FAD2d pow(const FAD2d k, const FAD2d u) {
    double kpu = ::pow(k.value,u.value) ;
    double kpu1 = ::pow(k.value,u.value-1.) ;
    double kpu2 = ::pow(k.value,u.value-2.) ;
    double lxu = ::log(k.value) ;
    return FAD2d(kpu,kpu1*k.grad*u.value + kpu*lxu*u.grad,
		 u.value*(k.grad*(kpu2*(k.grad*(u.value - 1)) + 
				  kpu1*(lxu*u.grad)) + kpu1*k.grad2) +
		 2*(kpu1*(k.grad*u.grad)) + 
		 lxu*(u.grad*(kpu1*(k.grad*u.value) + 
			      kpu*(lxu*u.grad)) + kpu*u.grad2)) ;

  }
  inline FAD2d pow(const int k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const float k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const double k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const long double k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }

  inline FAD2d sinh(const FAD2d u) {
    double expu = ::exp(u.value) ;
    double expmu = ::exp(-u.value) ;
    return FAD2d(::sinh(u.value),
		 0.5*u.grad*(expu + expmu),
		 0.5*(expmu*(u.grad2-u.grad*u.grad)+expu*(u.grad2+u.grad*u.grad))) ;
  }
  inline FAD2d cosh(const FAD2d u) {
    double expu = ::exp(u.value) ;
    double expmu = ::exp(-u.value) ;
    return FAD2d(::cosh(u.value),
		 0.5*u.grad*(expu - expmu),
		 0.5*(expmu*(u.grad*u.grad-u.grad2)+expu*(u.grad2+u.grad*u.grad))) ;
  }

  inline FAD2d tanh(const FAD2d u) {
    double ex = ::exp(min(u.value,350.0)) ;
    double exm = ::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    double rex = dex/sex ;
    double rex2 = rex*rex ;
    return FAD2d(::tanh(u.value),u.grad*(1.-rex2),
		 u.grad*u.grad*2.*(rex2-1.)*rex+
		 u.grad2*(1.-rex2)) ;
  }
  inline FAD2d asin(const FAD2d u) {
    double rsrt = 1./::sqrt(1.-u.value*u.value) ;
    return FAD2d(::asin(u.value), u.grad*rsrt,
		 u.grad2*rsrt + u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d acos(const FAD2d u) {
    double rsrt = 1./::sqrt(1.-u.value*u.value) ;
    return FAD2d(::acos(u.value), -u.grad*rsrt,
		 -u.grad2*rsrt - u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d atan(const FAD2d u) {
    double rval = 1./(1.+u.value*u.value) ;
    return FAD2d(::atan(u.value), u.grad*rval,
		 u.grad2*rval-2.*u.grad*u.grad*u.value*rval*rval);
  }
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // This will not work in general
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  inline FAD2d atan2(const FAD2d u, const FAD2d v) {  
    return atan(u/v);
  }

  inline FAD2d asinh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = ::sqrt(uv2+1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::asinh(uv), u.grad*(uv/strm+1.)/uvstrm,
		 ((1. / strm - uv2 / ::pow(strm, 3)) / uvstrm 
		  - ::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;
  }
  inline FAD2d acosh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = ::sqrt(uv*uv-1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::acosh(uv), u.grad*(uv/strm +1.)/uvstrm,
		 ((1 / strm - uv2 / ::pow(strm, 3)) / uvstrm 
		  - ::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad 
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;	 
  }
  inline FAD2d atanh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double factor = .5*(1./(1.+uv)+1./(1.-uv)) ;
    return FAD2d(::atanh(uv), .5*u.grad*factor,
		 0.5*(1./ (1. - 2.*uv + uv2) - 1./ (1. + 2.*uv + uv2))*u.grad*u.grad + u.grad2*factor) ;
  }


#else
  inline double ceil(const FAD2d u) {
    return std::ceil(u.value) ;
  }
  inline double floor(const FAD2d u) {
    return std::floor(u.value) ;
  }
  inline FAD2d sin(const FAD2d u) {
    double su = std::sin(u.value) ;
    double cu = std::cos(u.value) ;
    return FAD2d(su, u.grad*cu, cu*u.grad2 - u.grad*u.grad*su);
  }
  inline FAD2d cos(const FAD2d u) {
    double su = std::sin(u.value) ;
    double cu = std::cos(u.value) ;
    return FAD2d(cu, -u.grad*su, -su*u.grad2 - u.grad*u.grad*cu);
  }
  inline FAD2d tan(const FAD2d u) {
    double tanu = std::tan(u.value) ;
    double secu = 1./std::cos(u.value) ;
    double secu2 = secu*secu ;
    return FAD2d(tanu, u.grad*secu2,
		 u.grad2*secu2+2.*u.grad*u.grad*secu2*tanu) ;
  }

  inline FAD2d exp(const FAD2d u) {
    double eu = std::exp(u.value) ;
    return FAD2d(eu, u.grad*eu, eu*u.grad2+u.grad*u.grad*eu);
  }
  inline FAD2d log(const FAD2d u) {
    return FAD2d(std::log(u.value), u.grad/u.value, 
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value));
  }
  inline FAD2d log10(const FAD2d u) {
    return FAD2d(std::log(u.value), u.grad/u.value,
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value))/std::log(10.0);
  }

  inline FAD2d fabs(FAD2d u) {
    return FAD2d(std::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  inline FAD2d abs(FAD2d u) {
    return FAD2d(std::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FAD2d pow(const FAD2d u, const int k) {
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k - 2)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const float k) {
    return FAD2d(std::pow(u.value, double(k)), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k-2.0)*std::pow(u.grad, 2)) + std::pow(u.value, k-1.0)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const double k) {
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k-2.0)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const long double ki) {
    double k = double(ki) ;
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k - 2)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d sqrt(const FAD2d u) {

    double su = std::sqrt(u.value) ;
    return FAD2d(su,0.5*u.grad/max(su,1e-30),
		 0.5*u.grad2/(max(su,1e-30)) - 0.25*u.grad*u.grad/(max(su*su*su,1e-30))) ;
  }
  inline FAD2d pow(const FAD2d k, const FAD2d u) {
    double kpu = std::pow(k.value,u.value) ;
    double kpu1 = std::pow(k.value,u.value-1.) ;
    double kpu2 = std::pow(k.value,u.value-2.) ;
    double lxu = std::log(k.value) ;
    return FAD2d(kpu,kpu1*k.grad*u.value + kpu*lxu*u.grad,
		 u.value*(k.grad*(kpu2*(k.grad*(u.value - 1)) + 
				  kpu1*(lxu*u.grad)) + kpu1*k.grad2) +
		 2*(kpu1*(k.grad*u.grad)) + 
		 lxu*(u.grad*(kpu1*(k.grad*u.value) + 
			      kpu*(lxu*u.grad)) + kpu*u.grad2)) ;

  }
  inline FAD2d pow(const int k, const FAD2d u) {
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const float k, const FAD2d u) {
    double kpu = std::pow(double(k),u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const double k, const FAD2d u) {
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const long double ki, const FAD2d u) {
    double k = double(ki) ;
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }

  inline FAD2d sinh(const FAD2d u) {
    double expu = std::exp(u.value) ;
    double expmu = std::exp(-u.value) ;
    return FAD2d(std::sinh(u.value),
		 0.5*u.grad*(expu + expmu),
		 0.5*(expmu*(u.grad2-u.grad*u.grad)+expu*(u.grad2+u.grad*u.grad))) ;
  }
  inline FAD2d cosh(const FAD2d u) {
    double expu = std::exp(u.value) ;
    double expmu = std::exp(-u.value) ;
    return FAD2d(std::cosh(u.value),
		 0.5*u.grad*(expu - expmu),
		 0.5*(expmu*(u.grad*u.grad-u.grad2)+
		      expu*(u.grad2+u.grad*u.grad))) ;
  }

  inline FAD2d tanh(const FAD2d u) {
    double ex = std::exp(min(u.value,350.0)) ;
    double exm = std::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    double rex = dex/sex ;
    double rex2 = rex*rex ;
    return FAD2d(std::tanh(u.value),u.grad*(1.-rex2),
		 u.grad*u.grad*2.*(rex2-1.)*rex+
		 u.grad2*(1.-rex2)) ;
  }
  inline FAD2d asin(const FAD2d u) {
    double rsrt = 1./std::sqrt(1.-u.value*u.value) ;
    return FAD2d(std::asin(u.value), u.grad*rsrt,
		 u.grad2*rsrt + u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d acos(const FAD2d u) {
    double rsrt = 1./std::sqrt(1.-u.value*u.value) ;
    return FAD2d(std::acos(u.value), -u.grad*rsrt,
		 -u.grad2*rsrt - u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d atan(const FAD2d u) {
    double rval = 1./(1.+u.value*u.value) ;
    return FAD2d(std::atan(u.value), u.grad*rval,
		 u.grad2*rval-2.*u.grad*u.grad*u.value*rval*rval);
  }
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // This will not work in general
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  inline FAD2d atan2(const FAD2d u, const FAD2d v) {  
    return atan(u/v);
  }

  inline FAD2d asinh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = std::sqrt(uv2+1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::asinh(uv), u.grad*(uv/strm+1.)/uvstrm,
		 ((1. / strm - uv2 / std::pow(strm, 3)) / uvstrm 
		  - std::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;
  }
  inline FAD2d acosh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = std::sqrt(uv*uv-1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::acosh(uv), u.grad*(uv/strm +1.)/uvstrm,
		 ((1 / strm - uv2 / std::pow(strm, 3)) / uvstrm 
		  - std::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad 
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;	 
  }
  inline FAD2d atanh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double factor = .5*(1./(1.+uv)+1./(1.-uv)) ;
    return FAD2d(::atanh(uv), .5*u.grad*factor,
		 0.5*(1./ (1. - 2.*uv + uv2) -
		      1./ (1. + 2.*uv + uv2))*u.grad*u.grad + u.grad2*factor) ;
  }


#endif
  

  inline FAD2d operator +(const double u,const FAD2d v)
  { return FAD2d(u+v.value, v.grad, v.grad2); }
  inline FAD2d operator -(const double u,const FAD2d v) 
  { return FAD2d(u-v.value, -v.grad, -v.grad2); }
  inline FAD2d operator *(const double u,const FAD2d v) 
    { return FAD2d(u*v.value,  v.grad*u, v.grad2*u ); }
  inline FAD2d operator /(const double &u,const FAD2d &v) 
  { return FAD2d(u/v.value, -u*v.grad/v.value/v.value,
		 2.*u*v.grad*v.grad/(v.value*v.value*v.value) -
		 u*v.grad2/(v.value*v.value)); }

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
  FADd(const FAD2d &u): value(u.value), grad(u.grad) { gradcheck();}
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
  inline double ceil(const FADd u) {
    return ::ceil(u.value) ;
  }
  inline double floor(const FADd u) {
    return ::floor(u.value) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(::sin(u.value), u.grad*::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(::cos(u.value), -u.grad*::sin(u.value));
  }
  inline FADd tan(const FADd u) {
    return FADd(::tan(u.value), u.grad/pow(::cos(u.value),2)) ;
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
    return FADd(::sinh(u.value),
		0.5*u.grad*(::exp(u.value) + ::exp(-(u.value)))) ;
  }
  inline FADd cosh(const FADd u) {
    return FADd(::cosh(u.value),
		0.5*u.grad*(::exp(u.value) - ::exp(-(u.value)))) ;
  }
  inline FADd tanh(const FADd u) {
    double ex = ::exp(min(u.value,350.0)) ;
    double exm = ::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    return FADd(::tanh(u.value),u.grad*(1.-dex*dex/(sex*sex))) ;
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
    double strm = ::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::asinh(u.value), u.grad*(u.value/strm+1.)/uvstrm) ;
  }
  inline FADd acosh(const FADd u) {
    double strm = ::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::acosh(u.value), u.grad*(u.value/strm +1.)/uvstrm) ;
  }
  inline FADd atanh(const FADd u) {
    const double uv = u.value ;
    return FADd(::atanh(uv), .5*u.grad*(1./(1.+uv)+1./(1.-uv))) ;
  }


#else
  inline double ceil(const FADd u) {
    return std::ceil(u.value) ;
  }
  inline double floor(const FADd u) {
    return std::floor(u.value) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(std::sin(u.value), u.grad*std::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(std::cos(u.value), -u.grad*std::sin(u.value));
  }
  inline FADd tan(const FADd u) {
    return FADd(std::tan(u.value), u.grad/pow(std::cos(u.value),2)) ;
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
    return FADd(std::pow(u.value, double(k)), 
		(double)k * std::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const double k) {
    return FADd(std::pow(u.value, k), 
		k * std::pow(u.value, k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const long double k) {
    return FADd(std::pow(u.value,  double(k)), 
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
    double kpu = std::pow(double(k),u.value) ;
    return FADd(kpu,kpu*std::log(double(k))*u.grad) ;
  }
  inline FADd pow(const double k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(k)*u.grad) ;
  }
  inline FADd pow(const long double k, const FADd u) {
    double kpu = std::pow(double(k),u.value) ;
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
    double strm = std::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::asinh(u.value), u.grad*(u.value/strm+1.)/uvstrm) ;
  }
  inline FADd acosh(const FADd u) {
    double strm = std::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::acosh(u.value), u.grad*(u.value/strm +1.)/uvstrm) ;
  }
  inline FADd atanh(const FADd u) {
    const double uv = u.value ;
    return FADd(::atanh(uv), .5*u.grad*(1./(1.+uv)+1./(1.-uv))) ;
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

  //----------------helpers
  inline float realToFloat(double v) { return float(v) ; }
  inline double realToDouble(double v) { return v ; } 

  inline float realToFloat(const FADd &v) { return float(v.value) ; }
  inline double realToDouble(const FADd &v) { return v.value ; }

  inline float realToFloat(const FAD2d &v) { return float(v.value) ; }
  inline double realToDouble(const FAD2d &v) { return v.value ; }

}

#endif

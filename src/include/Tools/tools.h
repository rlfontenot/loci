#ifndef TOOLS_H
#define TOOLS_H

#ifdef NO_CMATH
#include <math.h>
#else
#include <cmath>

// trigonemetric
using std::acos ;
using std::asin ;
using std::atan ;
using std::atan2 ;
using std::cos ;
using std::sin ;
using std::tan ;
// hyperbolic
using std::sinh ;
using std::cosh ;
using std::tanh ;
// exponetial and logrithmic
using std::exp ;
using std::log ;
using std::log10 ;
using std::sqrt ;
using std::pow ;
using std::fabs ;
// misc
using std::ceil ;
using std::floor ;
using std::fmod ;
#endif
#ifdef NO_CSTDLIB
#include <stdlib.h>
inline float abs(float f1) { return f1>0?f1:-f1; }
inline double abs(double f1) { return f1>0?f1:-f1; }
inline long double abs(long double f1) { return f1>0?f1:-f1; }
#else
#include <cstdlib>
using std::abs ;
#endif

#include <algorithm>
using std::max ;
using std::min ;

#include <utility>
#include <iosfwd>

template <class T> inline T square(const T &a)   { return a*a ; }

template <class T> inline T sign(T a, T b) { return (b>=0)?abs(a):-abs(a) ; }


#endif

#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <iosfwd>

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
// misc
using std::max ;
using std::min ;
using std::abs ;
using std::fabs ;
using std::ceil ;
using std::floor ;
using std::fmod ;

template <class T> inline T square(const T &a)   { return a*a ; }

template <class T> inline T sign(T a, T b) { return (b>=0)?abs(a):-abs(a) ; }


#endif

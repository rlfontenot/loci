#ifndef TOOLS_STREAM_H
#define TOOLS_STREAM_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#ifdef GXX_FIXES
#include <g++-fixes/istream>
#include <g++-fixes/ostream>
#include <g++-fixes/sstream>
#else
#include <istream>
#include <ostream>
#include <sstream>
#endif
#include <iostream>
#include <fstream>
#include <string>

using std::istream ;
using std::ostream ;
using std::endl ;
using std::cin ;
using std::cout ;
using std::cerr ;
using std::ios ;
using std::ofstream ;
using std::ifstream ;

using std::istringstream ;
using std::ostringstream ;
using std::string ;
using std::char_traits ;

  
#endif

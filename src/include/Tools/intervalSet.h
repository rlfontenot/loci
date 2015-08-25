//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#ifndef INTERVALSET_H
#define INTERVALSET_H 
#define ENTITY64

#include <Tools/intervalSet_def.h>
#include <Tools/intervalSet_impl.h>
namespace Loci {
  typedef  genIntervalSet<int_type> intervalSet;
  typedef genSequence<int_type> sequence;
  extern const genIntervalSet<int_type> EMPTY ;
  extern const int UNIVERSE_MIN ;
  extern const int UNIVERSE_MAX ;
  typedef std::pair<int_type, int_type> interval;
}
#endif

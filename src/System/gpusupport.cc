//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include "loci_globs.h"
#include "sched_tools.h"
#include "dist_tools.h"
#include "param_rule.h"
#include "thread.h"
#include <Tools/except.h>
#include <constraint.h>
#include <new>
using std::bad_alloc ;
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ios ;
using std::ofstream ;
using std::istream ;
using std::ostream ;

///////////////////////////////////
#include <sstream>
#include <algorithm>
#include <sys/time.h> // for gettimeofday function
///////////////////////////////////

namespace Loci {
  rule_db rename_gpu_containers(fact_db  &facts,
				const rule_db &rdb)  {
    return rdb ;
  }
}

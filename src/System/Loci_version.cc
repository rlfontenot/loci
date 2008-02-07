//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#include <Loci_version.h>

namespace Loci {
  
  const char *revision_name = "$Name:  $" ;

  std::string version() {
    const char *p = revision_name;
    while(*p!=':' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    while(*p!=' ' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    std::string rn ;
    while(*p!='$' &&  *p!=' ' && *p!='\0') 
      rn += *p++ ;

    rn += " Compiled On " ;
    rn += __DATE__ ;
    rn += " " ;
    rn += __TIME__ ;
    return rn ;
  }
}

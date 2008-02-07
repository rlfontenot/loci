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
#ifndef LOCI_XDR_HEADER
#define LOCI_XDR_HEADER
#include "Config/conf.h"

#include <rpc/rpc.h>
#include <rpc/xdr.h>
#ifdef NO_XDR_CPP_PROTOTYPES 
extern
"C" {
  extern bool_t xdr_int(XDR *xdrs, int *ip) ;
  extern bool_t xdr_double(XDR *xdrs, double *dp) ;
  extern bool_t xdr_char(XDR *xdrs, char *cp) ;
  extern void xdrstdio_create (XDR *xdrs, FILE *file, enum xdr_op xop) ;
  //  extern void xdr_destroy(XDR *xdrs) ;
  // HACK
#undef xdr_destroy
#define xdr_destroy(X)

}
#endif

#endif

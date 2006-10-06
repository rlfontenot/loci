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
  extern void xdrstdio_create (XDR *xdrs, FILE *file, enum xdr_op xop) ;
  //  extern void xdr_destroy(XDR *xdrs) ;
  // HACK
#undef xdr_destroy
#define xdr_destroy(X)

}
#endif

#endif

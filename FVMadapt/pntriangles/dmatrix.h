#ifndef DMATRIX_H
#define DMATRIX_H
#include <vector>

class dmatrix {
  int rstart, cstart ;
  int nr, nc ;
  std::vector<double> data ;
 public:
  dmatrix(int nrl, int nrh, int ncl, int nch) {
    nr = nrh-nrl+1 ;
    nc = nch-ncl+1 ;
    rstart = nrl ;
    cstart = ncl ;
    std::vector<double> tmp(nr*nc) ;
    data.swap(tmp) ;
  }
  double *operator[](int r) { double *col =  &data[(r-rstart)*nc] ; return col-rstart ; }
} ;

class dvector {
  int start ;
  std::vector<double> data ;
 public:
  dvector(int nl, int nh) {
    std::vector<double> tmp(nh-nl+1) ;
    start = nl ;
    data.swap(tmp) ;
  }
  double &operator[](int i) { return data[i-start]; }
} ;

void svdcmp(dmatrix &a, int m, int n, dvector &w, dmatrix &v);

#endif

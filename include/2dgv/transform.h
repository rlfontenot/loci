#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "grid.h"


class transform2d {
  double mat[2][3] ;
public:
  void identity() {
    mat[0][0]=1.;mat[0][1]=0.;mat[0][2]=0.;
    mat[1][0]=0.;mat[1][1]=1.;mat[1][2]=0.;
  }    
  transform2d() {
    identity() ;
  }
  positions apply_transform(const positions &p) {
    positions pt ;
    pt.x = mat[0][0]*p.x+mat[0][1]*p.y+mat[0][2] ;
    pt.y = mat[1][0]*p.x+mat[1][1]*p.y+mat[1][2] ;
    return pt ;
  }
  void scale_transform(double scale) {
    for(int i=0;i<2;++i)
      for(int j=0;j<3;++j)
        mat[i][j] *= scale ;
  }
  void scale_transform(double scalex, double scaley) {
    for(int i=0;i<3;++i) {
      mat[0][i] *= scalex ;
      mat[1][i] *= scaley ;
    }
  }
  void pan_transform(const positions &p) {
    mat[0][2] += p.x ;
    mat[1][2] += p.y ;
  }
  transform2d inverse_transform() {
    transform2d tr ;
    double r = 1.0/(mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]) ;
    tr.mat[0][0] = mat[1][1]*r ;
    tr.mat[0][1] = -mat[0][1]*r ;
    tr.mat[0][2] = (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1])*r ;
    tr.mat[1][0] = - mat[1][0]*r ;
    tr.mat[1][1] = mat[0][0]*r ;
    tr.mat[1][2] = (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2])*r ;
    return tr ;
  }
} ;

#endif
 

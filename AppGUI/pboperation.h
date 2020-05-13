#ifndef PBOPERATION_H
#define PBOPERATION_H
#include "defines.h"
typedef vector3d<double> vect3d;
vect3d get_wireframe_center(const vector<vect3d> &pos,
                            const vector<int>  &trias);

vect3d get_average_normal(const vector<vect3d> &pos,
                          const vector<int>  &trias);

bool angleBetween(positions3d v1,positions3d v2, double& angle, vect3d& axis) ;


#endif

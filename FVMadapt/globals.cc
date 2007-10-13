#include "globals.h"

double Globals::tolerance = 1e-10;
double  Globals::fold = 90.0*3.1415926/180.0; //90 degree
int Globals::levels = 1;
int Globals::factor = 2;
vect3d Globals::split = vect3d(0.0, 0.0, 1.0);
vect3d Globals::nosplit = vect3d(0.0, 0.0, 1.0);

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
#include "globals.h"

double Globals::tolerance = 1e-10;
double  Globals::fold = 90.0*3.1415926/180.0; //90 degree
int Globals::levels = 1;
int Globals::factor = 2;
vect3d Globals::split = vect3d(0.0, 0.0, 1.0);
vect3d Globals::nosplit = vect3d(0.0, 0.0, 1.0);

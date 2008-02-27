/*-*- C++ -*-*/
#ifndef SCITYPES_H
#define SCITYPES_H

#include <Loci_types.h>

const double EPSILON = 1e-30 ;

typedef double real;

typedef Loci::vector3d<unsigned int> vect3u ;
typedef Loci::vector3d<double> vect3d ;
typedef Loci::vector3d<vect3d> tens3d ;

// Inner product of two 3D vectors is the dot product
// dot(vect3d, vect3d) is already defined
// double dot(const vect3d lhs, const vect3d rhs);

// Inner product of a 3D tensor of rank 2 and a 3D vector
vect3d inner(const tens3d & lhs, const vect3d & rhs);

// Inner product of a 3D vector and a 3D tensor of rank 2
vect3d inner(const vect3d & v, const tens3d & t);

// Outer product of two 3d vectors
tens3d outer(const vect3d & lhs, const vect3d & rhs);

// This operation is similar to 3x3 Matrix-Matrix multiplication
tens3d product(const tens3d & t1, const tens3d & t2);

// Compute the determinant of the 3x3 Matrix
double determinant(const tens3d & M);

// Compute the transpose of the 3x3 Matrix
tens3d transpose(const tens3d & M);

// Compute the inverse of the 3x3 Matrix
tens3d inverse(const tens3d & M);

#endif

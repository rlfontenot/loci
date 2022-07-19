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
#ifndef SIMPLE_PARTITION_LONG_H
#define SIMPLE_PARTITION_LONG_H
#include <Tools/intervalSet.h>

using std::vector;

namespace Loci{
  //return the start number of each processor  
template<class T> vector<T> g_simple_partition_vec(T mn, T mx, int p) {
    vector<T> nums(p+1) ;
    T n = mx-mn+1 ;
    T dn = n/p ; // divisor
    T rn = n%p ; // remainder
    T start = mn ;
    nums[0] = start ;
    for(int i=0;i<p;++i) {
      start += dn+((i<rn)?1:0) ;
      nums[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return nums ;
}
  
  //return num of entities of each processor
  template<class T> vector<T> g_simple_partition_dist(T mn, T mx, int p) {
    vector<T> nums(p) ;
    T n = mx-mn+1 ;
    T dn = n/p ; // divisor
    T rn = n%p ; // remainder
  
    for(int i=0;i<p;++i) {
      nums[i]= dn+((i<rn)?1:0) ;
    }
    return nums ;
  }

  //return the entitiySet of each processor
template<class T> vector<genIntervalSet<T> > g_simple_partition(T mn, T mx, MPI_Comm comm) {
  int p = 1 ;
  MPI_Comm_size(comm,&p) ;
  vector<T> pl = g_simple_partition_vec(mn,mx,p) ;
  vector<genIntervalSet<T> > ptn(p) ;
  for(int i=0;i<p;++i)
    ptn[i] = std::pair<T, T>(pl[i],pl[i+1]-1) ;
  return ptn ;
}

//return the entitiySet of each processor
template<class T> vector<genIntervalSet<T> > g_simple_partition(T num_entities, MPI_Comm comm) {
  T mn = 0;
  T mx = num_entities -1;
  return g_simple_partition(mn,mx,comm) ;
}  

}
#endif

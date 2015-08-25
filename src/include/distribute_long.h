//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#ifndef DISTRIBUTE_LONG_H
#define DISTRIBUTE_LONG_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <typeinfo>
#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <mpi.h>



namespace Loci {

  
  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  
  template<class T> inline bool g_spec_ival_compare(const std::pair<T, T> &i1,
                                                  const std::pair<T, T> &i2) {
    if(i1.first < i2.first)
      return true ;
    if(i1.first == i2.first && i1.second > i2.second)
      return true ;
    return false ;
  }
 //  template<class T> T g_GLOBAL_OR(T b, MPI_Comm comm=MPI_COMM_WORLD) {
//     T result ;
//     if(sizeof(T)==sizeof(int)){
//       MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LOR, comm) ;
//     }else{
//       MPI_Allreduce(&b, &result, 1, MPI_LONG, MPI_LOR, comm) ;
//     }
//     return result ;
//   }
  
//   template<class T> T  g_GLOBAL_AND(T b, MPI_Comm comm=MPI_COMM_WORLD) {
//     T result ;
//     if(sizeof(T)==sizeof(int)){ 
//       MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LAND, comm) ;
//     }else{
//       MPI_Allreduce(&b, &result, 1, MPI_LONG, MPI_LAND, comm) ;
//     }
//     return result ;
//   }
  
  template<class T> T g_GLOBAL_MAX(T b, MPI_Comm comm=MPI_COMM_WORLD) {
    T result ;
    if(sizeof(T)==sizeof(int)){
      MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MAX, comm) ;
    }else{
      MPI_Allreduce(&b, &result, 1, MPI_GENTITY_TYPE, MPI_MAX, comm) ;
     }
    return result ;
  }
  
  template<class T> T g_GLOBAL_MIN(T b, MPI_Comm comm=MPI_COMM_WORLD) {
    T result ;
    if(sizeof(T)==sizeof(int)){
      MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MIN, comm) ;
    }else{
      MPI_Allreduce(&b, &result, 1, MPI_GENTITY_TYPE, MPI_MIN, comm) ;
    }
    return result ;
  }
    
  // Collect largest interval of entitySet from all processors
  template<class T> genIntervalSet<T>  g_collectLargest(const  genIntervalSet<T>&e) {
    const int p = MPI_processes ;
    // Else we compute set
    //First get largest interval
    std::pair<T, T> ivl_large(1,-1) ;
    
    const long esz = e.num_intervals() ;
    for(long i=0;i<esz;++i)
      if((ivl_large.second-ivl_large.first) < (e[i].second-e[i].first))
        ivl_large = e[i] ;
    
    std::vector<std::pair<T, T> > ivl_large_p(p) ;
    if(sizeof(T)==sizeof(int)){
      MPI_Allgather(&ivl_large,2,MPI_INT,&(ivl_large_p[0]),2,MPI_INT,
                    MPI_COMM_WORLD) ;
    }else{
      MPI_Allgather(&ivl_large,2,MPI_GENTITY_TYPE,&(ivl_large_p[0]),2,MPI_GENTITY_TYPE,
                    MPI_COMM_WORLD) ; 
    }
    genIntervalSet<T> lset ;
    for(int i=0;i<p;++i)
      if(ivl_large_p[i].first <= ivl_large_p[i].second)
        lset += ivl_large_p[i] ;
    
    return lset ;
  }

  // Return union of all entitySets from all processors
  template<class T> genIntervalSet<T>  g_all_gather_entitySet(const  genIntervalSet<T> &e) {
    const int p = MPI_processes ;
    if(p == 1)
      return e ;
    
    T send_count = 2*e.num_intervals() ; //must be int*
    std::vector<T> recv_count(p) ;//must be int*
    if(sizeof(T)==sizeof(int)){
      MPI_Allgather(&send_count,1,MPI_INT,&recv_count[0],1,MPI_INT,
                    MPI_COMM_WORLD) ;
    }else{
      MPI_Allgather(&send_count,1,MPI_GENTITY_TYPE,&recv_count[0],1,MPI_GENTITY_TYPE,
                    MPI_COMM_WORLD) ;
    }
    std::vector<T> recv_disp(p) ;//must be int*
    recv_disp[0] = 0 ;
    for(int i=1;i<p;++i)
      recv_disp[i] = recv_disp[i-1]+recv_count[i-1] ;
    T tot = recv_disp[p-1]+recv_count[p-1] ;
    if(tot == 0)
      return  genIntervalSet<T>::EMPTY ;
    std::vector<std::pair<T, T> > ivl_list(tot/2) ;
    std::vector<std::pair<T, T> > snd_list(send_count/2) ;
    for(long i=0;i<send_count/2;++i)
      snd_list[i] = e[i] ;
     if(sizeof(T)==sizeof(int)){
       MPI_Allgatherv(&(snd_list[0]),send_count,MPI_INT,
                      &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_INT,
                      MPI_COMM_WORLD) ;
     }else{
       MPI_Allgatherv(&(snd_list[0]),send_count,MPI_GENTITY_TYPE,
                      &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_GENTITY_TYPE,
                      MPI_COMM_WORLD) ; 
     }
     std::sort(ivl_list.begin(),ivl_list.end(),g_spec_ival_compare<T>) ;
     genIntervalSet<T> tmp = ivl_list[0] ;
     for(size_t i=1;i<ivl_list.size();++i)
       tmp += ivl_list[i] ;
     return tmp ;
  }
  
  template<class T> genIntervalSet<T>  g_all_collect_entitySet(const genIntervalSet<T> &e) {
    const int p = MPI_processes ;
    // no operation for single processor
    if(p == 1)
      return e ;
    

    // Check to see if the result should be EMPTY or UNIVERSE
    int code = 0 ;
    if(e != genIntervalSet<T>::EMPTY)
      code = 1 ;
    if(e == ~(genIntervalSet<T>::EMPTY))
      code = 2 ;

    code = g_GLOBAL_MAX(code,MPI_COMM_WORLD) ;
    
    if(code == 0) // All empty, so return empty
      return genIntervalSet<T>::EMPTY ;
    if(code == 2) // at least one UNIVERSE, so return UNIVERSE
      return ~genIntervalSet<T>::EMPTY ;
    
    // First collect largest intervals, most of the time this will get
    // us the final set.  Compute remainder to determine if we have
    // more work to do.  Repeat collecting the largest interval for a
    // fixed number of tries, then finally just gather the remainder
#ifdef VERBOSE
    stopWatch s ;
    s.start() ;
#endif
    genIntervalSet<T> lset = g_collectLargest(e) ;
    genIntervalSet<T> rem = e-lset ;
    for(int i=0;i<4;++i) {
      T remsz = rem.num_intervals() ;
      T szmx = g_GLOBAL_MAX(remsz,MPI_COMM_WORLD) ;
      if(szmx == 0) {
#ifdef VERBOSE
        debugout << "time to get lset = " << s.stop() << endl ;
#endif
        return lset ;
      }
      lset += g_collectLargest(rem) ;
      rem -= lset ;
    }
#ifdef VERBOSE
    debugout << "time to get lset = " << s.stop() << endl ;
    s.start() ;
    
    debugout << "e="<< e.num_intervals() << ",rem=" << rem.num_intervals()
             << ",lset=" << lset.num_intervals() << endl ;
#endif
    genIntervalSet<T> remtot = g_all_gather_entitySet(rem) ;
#ifdef VERBOSE
    debugout << "time to gather rem = " << s.stop() << endl ;
#endif
    return lset + remtot ;
  } ;

 
}

#endif
 

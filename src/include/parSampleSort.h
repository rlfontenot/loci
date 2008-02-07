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
#ifndef PAR_SAMPLE_SORT_H
#define PAR_SAMPLE_SORT_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <mpi.h>

#include <vector>
#include <algorithm>

namespace Loci {
  // With user supplied comparison operator

  // Utility routine for sample sort
  template <class T,class Cmp>
  void parGetSplitters(std::vector<T> &splitters,
                       const std::vector<T> &input,
                       Cmp cmp,
                       MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    splitters = std::vector<T>(p-1) ;
    std::vector<T> allsplits(p*(p-1)) ;

    int nlocal = input.size() ;
    if(nlocal < p) {
      std::cerr << "sample sort needs at least p elements per processor"
                << std::endl ;
    }
    for(int i=1;i<p;++i) 
      splitters[i-1] = input[(i*nlocal)/p] ;

    int tsz = sizeof(T) ;
    MPI_Allgather(&splitters[0],(p-1)*tsz,MPI_BYTE,
                  &allsplits[0],(p-1)*tsz,MPI_BYTE,comm) ;
    
    std::sort(allsplits.begin(),allsplits.end(),cmp) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }


  template <class T, class Cmp>
  void parSplitSort(std::vector<T> &list, std::vector<T> &splitters,
                    Cmp cmp, MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;

    int s=0 ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || cmp(list[i] , splitters[s]) ) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !cmp(list[i],splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(T) ;

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    std::vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(T) ;

    std::vector<T> sorted_pnts(result_size) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    list.swap(sorted_pnts) ;
    std::sort(list.begin(),list.end(),cmp) ;
    return ;
  }

  template <class T, class Cmp>
  void parSampleSort(std::vector<T> &list, Cmp cmp, MPI_Comm comm) {
    // First sort list locally
    std::sort(list.begin(),list.end(),cmp) ;

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    if(p == 1) // if serial run, we are finished
      return ;

    std::vector<T> splitters ;
    parGetSplitters(splitters,list,cmp,comm) ;

    parSplitSort(list,splitters,cmp,comm) ;
  }


  // Now same code with default operator

  // Utility routine for sample sort
  template <class T> void parGetSplitters(std::vector<T> &splitters,
                                          const std::vector<T> &input,
                                          MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    splitters = std::vector<T>(p-1) ;
    std::vector<T> allsplits(p*(p-1)) ;

    int nlocal = input.size() ;
    if(nlocal < p) {
      std::cerr << "sample sort needs at least p elements per processor"
                << std::endl ;
    }
    for(int i=1;i<p;++i) 
      splitters[i-1] = input[(i*nlocal)/p] ;

    int tsz = sizeof(T) ;
    MPI_Allgather(&splitters[0],(p-1)*tsz,MPI_BYTE,
                  &allsplits[0],(p-1)*tsz,MPI_BYTE,comm) ;
    
    std::sort(allsplits.begin(),allsplits.end()) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }


  template <class T> void parSplitSort(std::vector<T> &list,
                                       std::vector<T> &splitters,
                                       MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;

    int s=0 ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || (list[i] < splitters[s]) ) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !(list[i] < splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(T) ;

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    std::vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(T) ;

    std::vector<T> sorted_pnts(result_size) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    list.swap(sorted_pnts) ;
    std::sort(list.begin(),list.end()) ;
    return ;
  }

  template <class T> void parSampleSort(std::vector<T> &list, MPI_Comm comm) {
    // First sort list locally
    std::sort(list.begin(),list.end()) ;

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    if(p == 1) // if serial run, we are finished
      return ;

    std::vector<T> splitters ;
    parGetSplitters(splitters,list,comm) ;

    parSplitSort(list,splitters,comm) ;
  }

  template <class T> void balanceDistribution(std::vector<T> &list,
                                              MPI_Comm comm) {
    
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    if(p == 1) // if serial run, we are finished
      return ;

    int sz = list.size() ;
    std::vector<int> sizes(p) ;
    MPI_Allgather(&sz,1,MPI_INT,&sizes[0],1,MPI_INT,comm) ;
    int tot_size = 0  ;
    for(int i=0;i<p;++i) {
      tot_size += sizes[i] ;
    }
    std::vector<int> bsizes(p) ;
    int psz = tot_size/p ;
    int rem = tot_size%p ;
    for(int i=0;i<p;++i)
      bsizes[i] = psz + (i<rem?1:0) ;
    int r=0 ;
    MPI_Comm_rank(comm,&r) ;
    std::vector<T> recv(bsizes[r]) ;
    std::vector<int> sdist(p+1),rdist(p+1) ;
    sdist[0] = 0 ;
    rdist[0] = 0 ;
    for(int i=0;i<p;++i) {
      sdist[i+1] = sizes[i]+sdist[i] ;
      rdist[i+1] = bsizes[i]+rdist[i] ;
    }
    FATAL(sdist[p] != rdist[p]) ;
    // Perform Irecvs
    std::vector<MPI_Request> requests(p) ;
    int req = 0 ;
    int i1 = rdist[r] ;
    int i2 = rdist[r+1]-1 ;
    for(int i=0;i<p;++i) {
      if(sdist[i]<=i2 && sdist[i+1]-1>=i1) { // intersection
        int li = std::max(i1,sdist[i]) ;
        int ri = std::min(i2,sdist[i+1]-1) ;
        int len = ri-li+1 ;
        FATAL(len <= 0) ;
        int s2 = li-rdist[r] ;
        if(i == r) { // local copy
          int s1 = li-sdist[r] ;
          FATAL(s1 < 0 || s2 < 0) ;
          for(int j=0;j<len;++j)
            recv[s2+j] = list[s1+j] ;
        } else {
          MPI_Irecv(&recv[s2],len*sizeof(T),MPI_BYTE,i,55,comm,
                    &requests[req++]) ;
        }
      }
    }

    // Perform sends
    i1 = sdist[r] ;
    i2 = sdist[r+1]-1 ;
    for(int i=0;i<p;++i) {
      if(i != r && rdist[i]<=i2 && rdist[i+1]-1>=i1) { // intersection
        int li = std::max(i1,rdist[i]) ;
        int ri = std::min(i2,rdist[i+1]-1) ;
        int len = ri-li+1 ;
        int s1 = li-sdist[r] ;
        FATAL(s1 < 0) ;
        FATAL(len <= 0) ;
        MPI_Send(&list[s1],len*sizeof(T),MPI_BYTE,i,55,comm) ;
      }
    }

    if(req > 0) {
      std::vector<MPI_Status> status(p) ;
      MPI_Waitall(req,&requests[0],&status[0]) ;
    }
    list.swap(recv) ;
  }
}

#endif

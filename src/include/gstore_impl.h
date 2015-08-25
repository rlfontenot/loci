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

#ifndef GSTORE_IMPL_H
#define GSTORE_IMPL_H 1

#include <Tools/except.h>
#include <parSampleSort.h>
namespace Loci {
  
  
  namespace {
    template<class T> inline  bool gstore_porder(const std::pair<GEntity, T> &i1, 
                                                 const std::pair<GEntity, T> &i2) {
      return i1.first < i2.first ;
    }
  }


  template<class T> void GStore<T>::local_sort(){
    std::sort(Rep.begin(), Rep.end(), gstore_porder<T>);
  }

  template<class T> void GStore<T>::sort(){
    std::sort(Rep.begin(), Rep.end(), gstore_porder<T>);
    parSampleSort(Rep,gstore_porder,MPI_COMM_WORLD) ; 
  }

 

  

  //   typedef std::vector<std::pair<int,int> > protoMap;
  
  //   inline bool equiFF(const std::pair<int,int> &v1,
  //                      const std::pair<int,int> &v2) {
  //     return v1.first < v2.first ;
  //   }

  //   inline void equiJoinFF(protoMap &in1, protoMap &in2, protoMap &out) {
  //     std::sort(in1.begin(),in1.end(),equiFF) ;
  //     std::sort(in2.begin(),in2.end(),equiFF) ;
    
  //     int p = 0 ;
  //     MPI_Comm_size(MPI_COMM_WORLD,&p) ;

  //     // Sort inputs using same splitters (this will make sure that
  //     // data that needs to be on the same processor ends up on the
  //     // same processor
  //     if(p != 1) {
  //       std::vector<std::pair<int,int> > splitters ;
  //       parGetSplitters(splitters,in1,equiFF,MPI_COMM_WORLD) ;

  //       parSplitSort(in1,splitters,equiFF,MPI_COMM_WORLD) ;
  //       parSplitSort(in2,splitters,equiFF,MPI_COMM_WORLD) ;
  //     }

  //     // Find pairs where first entry are the same and create joined protomap
  //     out.clear() ;
  //     size_t j = 0 ;
  //     for(size_t i=0;i<in1.size();++i) {
  //       while(j<in2.size() && equiFF(in2[j],in1[i]))
  //         ++j ;
  //       size_t k=j ;
  //       while(k<in2.size() && in2[k].first == in1[i].first) {
  //         out.push_back(std::pair<int,int>(in1[i].second,in2[k].second)) ;
  //         k++ ;
  //       }
  //     }

  //     // Remove duplicates from protomap
  //     parSampleSort(out,equiFF,MPI_COMM_WORLD) ;
  //     std::sort(out.begin(),out.end()) ;
  //     out.erase(std::unique(out.begin(),out.end()),out.end()) ;
  //   }

                         
      
  
 

  
  //  //remap cr from negative ids to ref GEntity
  //   template<class T>  void  GStore<T>::remap(const std::map<T,T>& remap){
  //     typename gRep::iterator itr = Rep->begin();
  //     for(itr = Rep->begin(); itr != Rep->end(); itr++){
  //       T val = itr->second;
  //       typename std::map<T,T>::const_iterator mi = remap.find(val);
  //       if(mi != remap.end()){
  //         itr->second = mi->second;
  //       }
  //     }
  //   }
  
 
  
    
  
  template<class T> void GStore<T>::make_consistent(){
  }

  
  //*********************************************************************/

  template<class T> 
  std::ostream &GStore<T>::Print(std::ostream &s) const 
  {
    s << '{' << std::endl ;
    typename gRep::const_iterator itr;
    for( itr = Rep->begin();
         itr != Rep->end(); itr++){
      s << itr->first << " -> " << itr->second << std::endl ;
    }
   
    s << '}' << std::endl ;

    return s ;
  } 
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const GStore<T> &t)
  { return t.Print(s) ; }
}
#endif

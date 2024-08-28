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
#ifndef DISTRIBUTE_CONTAINER_H
#define DISTRIBUTE_CONTAINER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <vector>
#include <algorithm>

#include <DMap.h>
#include <store_rep.h>
#include <DMultiMap.h>
#include <fact_db.h>
#include <Map.h>
#include <multiMap.h>
#include <distribute_io.h>

namespace Loci {
  
  void distributed_inverseMap(multiMap &result,
                              std::vector<std::pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn) ;

  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) ;

  void distributed_inverseMap(dmultiMap &result, const dMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result, const Map &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const dmultiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const multiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn); 

  void distributed_inverseMap(dmultiMap &result,
                              const const_dMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_Map &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_dmultiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_multiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn); 

  inline void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }
  
  inline void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts,size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  inline void distributed_inverseMap(dmultiMap &result,
                                     const multiMap &input_map,
                                     const entitySet &input_image,
                                     const entitySet &input_preimage,
                                     fact_db &facts,
				     size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn);
  

  inline void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts,size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }


  typedef std::vector<std::pair<int,int> > protoMap;
  
  inline bool equiFF(const std::pair<int,int> &v1,
                     const std::pair<int,int> &v2) {
    return v1.first < v2.first ;
  }

  inline void equiJoinFF(protoMap &in1, protoMap &in2, protoMap &out) {
    std::sort(in1.begin(),in1.end(),equiFF) ;
    std::sort(in2.begin(),in2.end(),equiFF) ;
    
    int p = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD,&p) ;

    // Sort inputs using same splitters (this will make sure that
    // data that needs to be on the same processor ends up on the
    // same processor
    if(p != 1) {
      std::vector<std::pair<int,int> > splitters ;
      parGetSplitters(splitters,in1,equiFF,MPI_COMM_WORLD) ;

      parSplitSort(in1,splitters,equiFF,MPI_COMM_WORLD) ;
      parSplitSort(in2,splitters,equiFF,MPI_COMM_WORLD) ;
    }

    // Find pairs where first entry are the same and create joined protomap
    out.clear() ;
    size_t j = 0 ;
    for(size_t i=0;i<in1.size();++i) {
      while(j<in2.size() && equiFF(in2[j],in1[i]))
        ++j ;
      size_t k=j ;
      while(k<in2.size() && in2[k].first == in1[i].first) {
        out.push_back(std::pair<int,int>(in1[i].second,in2[k].second)) ;
        k++ ;
      }
    }

    // Remove duplicates from protomap
    parSampleSort(out,equiFF,MPI_COMM_WORLD) ;
    std::sort(out.begin(),out.end()) ;
    out.erase(std::unique(out.begin(),out.end()),out.end()) ;
  }

  inline void removeIdentity(protoMap &mp) {
      // Remove self references
    protoMap::iterator ii,ij ;
    for(ii=ij=mp.begin();ij!=mp.end();++ij) {
      if(ij->first != ij->second) {
        *ii = *ij ;
        ii++ ;
      }
    }
    mp.erase(ii,mp.end()) ;
  }
  inline void addToProtoMap(multiMap &m, protoMap &mp) {
    entitySet dom = m.domain() ;
    FORALL(dom,e) {
      int sz = m[e].size() ;
      for(int i=0;i<sz;++i)
        mp.push_back(std::pair<int,int>(e,m[e][i])) ;
    } ENDFORALL ;
  }

  inline void addToProtoMap(Map &m, protoMap &mp) {
    entitySet dom = m.domain() ;
    FORALL(dom,e) {
      mp.push_back(std::pair<int,int>(e,m[e])) ;
    } ENDFORALL ;
  }

}

#endif

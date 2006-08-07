#ifndef DISTRIBUTE_CONTAINER_H
#define DISTRIBUTE_CONTAINER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <vector>
#include <DMap.h>
#include <store_rep.h>
#include <DMultiMap.h>
#include <fact_db.h>
#include <Map.h>
#include <multiMap.h>

namespace Loci {
  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) ;

  extern fact_db *exec_current_fact_db ;
  storeRepP distribute_store(storeRepP &sp, fact_db &facts) ;
  inline storeRepP distribute_store(storeRepP &sp)
    { return distribute_store(sp,*exec_current_fact_db) ; }

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

  inline void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }
  
  inline void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  inline void distributed_inverseMap(dmultiMap &result,
                                     const multiMap &input_map,
                                     const entitySet &input_image,
                                     const entitySet &input_preimage,
                                     fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  storeRepP collect_global_store(storeRepP &sp) ;

  storeRepP collect_store(storeRepP &sp, fact_db &facts) ;
  extern storeRepP collect_store(storeRepP &sp, fact_db &facts) ;
  
  inline storeRepP collect_store(storeRepP &sp) 
  { return collect_store(sp,*exec_current_fact_db) ;}

  void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn);
  

  inline void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }


}

#endif

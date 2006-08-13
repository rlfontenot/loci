#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

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

#include <Map.h>
#include <DMap.h>
#include <fact_db.h>
#include <store_rep.h>

namespace Loci {

  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  
  void Init(int* argc, char*** argv) ;
  void Finalize() ; 
  void Abort() ;
  
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<dMap> send_global_map(Map &attrib_data, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  void fill_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  storeRepP send_clone_non(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) ;
  
  entitySet collect_entitySet(entitySet e, fact_db &facts) ;

  extern fact_db *exec_current_fact_db ;

  inline entitySet collect_entitySet(entitySet e)
    { return collect_entitySet(e,*exec_current_fact_db) ; }

  entitySet all_collect_entitySet(const entitySet &e) ;

  inline entitySet all_collect_entitySet(entitySet localset,fact_db &facts) {
    if(facts.is_distributed_start())
      return Loci::all_collect_entitySet(localset) ;
    return localset ;
  }
  
  std::vector<entitySet> all_collect_vectors(entitySet &e) ;
  int GLOBAL_OR(int b) ;
  int GLOBAL_AND(int b) ;
  int GLOBAL_MAX(int b) ;
  int GLOBAL_MIN(int b) ;
  
  // We've added these back as they seem to be used
  // in the fuel cell program//from distribute.h //////
  Map distribute_global_map(Map &m, const std::vector<entitySet> &vset) ;
  Map distribute_gmap(Map &m, const std::vector<entitySet> &vset) ;

}

#endif
 

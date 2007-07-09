#ifndef PNN_H
#define PNN_H
#include "kd_tree.h"
#include <mpi.h>

namespace Loci {
  // Find nearest neighbor in parallel.
  // Inputs:
  // target_pnts:  An array of that are searched for a closest match
  // target_ids:   An array of integer identifiers of target_pnts
  // search_pnts:  An array of search points (for each of these points we are
  //               searching for the nearest neighbor from target_pnts
  // comm:         MPI communicator for processors participating in NN search
  // rebalance:    Useful for randomly scattered data, do not use if inputs
  //               have good spatial locality, this option is not optimal for
  //               surface/volume searches, but can be good for volume/volume
  //               searches.
  // Outputs:
  // closest:      An array for each search_pnts giving the closest target_pnts
  //               integer id (as passed in as target_ids)
  void parallelNearestNeighbors(const std::vector<kdTree::coord3d> &target_pnts,
                                const std::vector<int> &target_ids,
                                const std::vector<kdTree::coord3d> &search_pnts,
                                std::vector<int> &closest,
                                MPI_Comm comm,bool rebalance=false) ;
}


#endif



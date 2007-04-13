#ifndef PNN_H
#define PNN_H
#include <Loci>
#include "kd_tree.h"
#include <mpi.h>

namespace Loci {
  void parallelNearestNeighbors(const std::vector<kdTree::coord3d> &target_pnts,
                                const std::vector<int> &target_ids,
                                const std::vector<kdTree::coord3d> &search_pnts,
                                std::vector<int> &closest,
                                MPI_Comm comm) ;
}


#endif



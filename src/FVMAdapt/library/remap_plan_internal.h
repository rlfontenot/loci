//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//#############################################################################
#ifndef FVMADAPT_REMAP_PLAN_INTERNAL_H
#define FVMADAPT_REMAP_PLAN_INTERNAL_H

#include <Loci.h>
#include <FVMAdapt/remap_plan.h>

namespace Loci {
  namespace detail {

    /// Assemble and collectively validate the target-owned cell remap slices.
    bool buildDistributedCellRemapPlan(
      CPTR<AMRRemapPlan>& plan,
      AMRRemapReport& report,
      gatherCommSchedule& sourceGather,
      const std::vector<std::pair<int,int> >& targetSource,
      const store<double>& sourceVolume,
      const_store<vector3d<double> >& sourceCenter,
      const store<int>& sourceTargetCount,
      dataPartitionP sourcePartition,
      const entitySet& localTargetCells,
      const Map& targetLocalToGlobal,
      const const_store<vector3d<double> >& targetCenter,
      const const_store<double>& targetVolume,
      const store<int>& refinedSource,
      const multiStore<int>& refinedSourceToTarget,
      const store<double>& gatheredTargetVolume,
      const store<vector3d<double> >& gatheredTargetCenter,
      MPI_Comm comm) ;
  }
}

#endif

//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# This program is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//#############################################################################

#include <Loci.h>
#include <FVMAdapt2/gridInterface.h>

#include <doctest.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Loci ;

namespace {

  void setup_refinement_mapping(AMRrefinementMapping& mapping) {
    const entitySet sourceCells = interval(0,4) ;
    const entitySet targetCells = interval(0,1) ;
    const entitySet oneSource = interval(0,0) ;

    mapping.geom_cells_local = targetCells ;
    mapping.geom_cells_global = targetCells ;
    mapping.l2g.allocate(targetCells) ;
    mapping.l2g[0] = 0 ;
    mapping.l2g[1] = 1 ;

    store<int> stencilSize ;
    stencilSize.allocate(oneSource) ;
    stencilSize[0] = 4 ;

    mapping.parent.allocate(oneSource) ;
    mapping.parent[0] = 0 ;
    store<int> targetCount ;
    targetCount.allocate(oneSource) ;
    targetCount[0] = 2 ;
    mapping.parent2child_l.allocate(targetCount) ;
    mapping.parent2child_l[0][0] = 0 ;
    mapping.parent2child_l[0][1] = 1 ;
    mapping.refinedCells = targetCells ;

    mapping.gradCellStencil.allocate(stencilSize) ;
    mapping.gradCellStencil[0][0] = 0 ;
    mapping.gradCellStencil[0][1] = 1 ;
    mapping.gradCellStencil[0][2] = 2 ;
    mapping.gradCellStencil[0][3] = 3 ;
    mapping.stencilWeights.allocate(stencilSize) ;
    mapping.stencilWeights[0][0] = vector3d<double>(0.5,0.0,0.0) ;
    mapping.stencilWeights[0][1] = vector3d<double>(-0.5,0.0,0.0) ;
    mapping.stencilWeights[0][2] = vector3d<double>(0.0,0.5,0.0) ;
    mapping.stencilWeights[0][3] = vector3d<double>(0.0,-0.5,0.0) ;
    mapping.grad_dvs.allocate(stencilSize) ;
    mapping.grad_dvs[0][0] = vector3d<double>(1.0,0.0,0.0) ;
    mapping.grad_dvs[0][1] = vector3d<double>(-1.0,0.0,0.0) ;
    mapping.grad_dvs[0][2] = vector3d<double>(0.0,1.0,0.0) ;
    mapping.grad_dvs[0][3] = vector3d<double>(0.0,-1.0,0.0) ;
    mapping.child_dvs.allocate(targetCount) ;
    // The two offsets have a zero volume-weighted average for fractions
    // 0.25 and 0.75, matching the conservative refinement reconstruction.
    mapping.child_dvs[0][0] = vector3d<double>(-0.75,0.0,0.0) ;
    mapping.child_dvs[0][1] = vector3d<double>(0.25,0.0,0.0) ;

    std::vector<entitySet> sourcePartition(1,sourceCells) ;
    std::vector<entitySet> targetPartition(1,targetCells) ;
    mapping.gradientComm.generateSchedule(
      interval(1,4),createPartition(sourcePartition,MPI_COMM_WORLD)) ;
    mapping.refineCellComm.generateSchedule(
      targetCells,createPartition(targetPartition,MPI_COMM_WORLD)) ;
    mapping.directMapComm.generateSchedule(
      EMPTY,createPartition(sourcePartition,MPI_COMM_WORLD)) ;
  }

  void setup_coarsening_mapping(AMRrefinementMapping& mapping) {
    const entitySet sourceCells = interval(0,1) ;
    const entitySet targetCells = interval(0,0) ;

    mapping.geom_cells_local = targetCells ;
    mapping.geom_cells_global = targetCells ;
    mapping.l2g.allocate(targetCells) ;
    mapping.l2g[0] = 0 ;

    store<int> noEntries ;
    noEntries.allocate(EMPTY) ;
    mapping.parent.allocate(EMPTY) ;
    mapping.parent2child_l.allocate(noEntries) ;
    mapping.gradCellStencil.allocate(noEntries) ;
    mapping.stencilWeights.allocate(noEntries) ;
    mapping.grad_dvs.allocate(noEntries) ;
    mapping.child_dvs.allocate(noEntries) ;
    mapping.refinedCells = EMPTY ;

    std::vector<entitySet> sourcePartition(1,sourceCells) ;
    std::vector<entitySet> targetPartition(1,targetCells) ;
    mapping.gradientComm.generateSchedule(
      EMPTY,createPartition(sourcePartition,MPI_COMM_WORLD)) ;
    mapping.refineCellComm.generateSchedule(
      EMPTY,createPartition(targetPartition,MPI_COMM_WORLD)) ;
    mapping.directMapComm.generateSchedule(
      sourceCells,createPartition(sourcePartition,MPI_COMM_WORLD)) ;

    mapping.c2lp.push_back(std::make_pair(0,0)) ;
    mapping.c2lp.push_back(std::make_pair(0,1)) ;
    mapping.directWeights.push_back(0.25) ;
    mapping.directWeights.push_back(0.75) ;
  }

  void allocate_state(storeVec<double>& state, const entitySet& domain) {
    state.setVecSize(4) ;
    state.allocate(domain) ;
  }

  std::vector<AMRcomponentGroup> closure_group() {
    std::vector<AMRcomponentGroup> groups(1) ;
    groups[0].components = {0,1,2} ;
    return groups ;
  }

} // namespace


/// A shared limiter should preserve an existing component identity on every
/// refined target while retaining conservation, bounds, and legacy behavior.
TEST_CASE("coupled AMR components preserve linear closure during refinement") {
  REQUIRE(MPI_processes == 1) ;
  AMRrefinementMapping mapping ;
  setup_refinement_mapping(mapping) ;

  storeVec<double> source ;
  allocate_state(source,interval(0,4)) ;
  source[0][0] = 1.0 ;
  source[0][1] = 0.4 ;
  source[0][2] = 0.6 ;
  source[0][3] = 10.0 ;
  source[1][0] = 2.0 ;
  source[1][1] = 1.2 ;
  source[1][2] = 0.8 ;
  source[1][3] = 12.0 ;
  source[2][0] = 0.8 ;
  source[2][1] = 0.3 ;
  source[2][2] = 0.5 ;
  source[2][3] = 9.0 ;
  source[3][0] = 3.0 ;
  source[3][1] = 1.2 ;
  source[3][2] = 1.8 ;
  source[3][3] = 10.0 ;
  source[4][0] = 0.9 ;
  source[4][1] = 0.4 ;
  source[4][2] = 0.5 ;
  source[4][3] = 10.0 ;

  storeVec<double> independent ;
  storeVec<double> coupled ;
  allocate_state(independent,interval(0,1)) ;
  allocate_state(coupled,interval(0,1)) ;
  mapping.interpolateData(source,independent) ;
  mapping.interpolateData(source,coupled,closure_group()) ;

  for(int component=0;component<4;++component) {
    const double coupledAverage =
      0.25*coupled[0][component] + 0.75*coupled[1][component] ;
    CHECK(coupledAverage == doctest::Approx(source[0][component])) ;
  }

  for(int target=0;target<2;++target) {
    CAPTURE(target) ;
    CHECK(std::abs(independent[target][0]-independent[target][1]-
                   independent[target][2]) > 1.0e-6) ;
    CHECK(coupled[target][0] ==
          doctest::Approx(coupled[target][1]+coupled[target][2])) ;

    for(int component=0;component<3;++component) {
      double lower = source[0][component] ;
      double upper = source[0][component] ;
      for(int stencilCell=1;stencilCell<5;++stencilCell) {
        lower = std::min(lower,source[stencilCell][component]) ;
        upper = std::max(upper,source[stencilCell][component]) ;
      }
      CHECK(coupled[target][component] >= lower) ;
      CHECK(coupled[target][component] <= upper) ;
    }

    CHECK(coupled[target][3] == doctest::Approx(independent[target][3])) ;
  }
}


/// Derefinement should apply the same geometric weights to every component,
/// preserving the same identity without requiring a limiter.
TEST_CASE("coupled AMR components preserve linear closure during coarsening") {
  REQUIRE(MPI_processes == 1) ;
  AMRrefinementMapping mapping ;
  setup_coarsening_mapping(mapping) ;

  storeVec<double> source ;
  allocate_state(source,interval(0,1)) ;
  source[0][0] = 1.0 ;
  source[0][1] = 0.25 ;
  source[0][2] = 0.75 ;
  source[0][3] = 4.0 ;
  source[1][0] = 3.0 ;
  source[1][1] = 2.0 ;
  source[1][2] = 1.0 ;
  source[1][3] = 8.0 ;

  storeVec<double> target ;
  allocate_state(target,interval(0,0)) ;
  mapping.interpolateData(source,target,closure_group()) ;

  CHECK(target[0][0] == doctest::Approx(2.5)) ;
  CHECK(target[0][1] == doctest::Approx(1.5625)) ;
  CHECK(target[0][2] == doctest::Approx(0.9375)) ;
  CHECK(target[0][0] == doctest::Approx(target[0][1]+target[0][2])) ;
  CHECK(target[0][3] == doctest::Approx(7.0)) ;
}

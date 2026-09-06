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
#include <FVMAdapt2/remap_plan.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <iostream>
#include <numeric>
#include <vector>

using Loci::AMRCellContribution ;
using Loci::AMRCellGeometry ;
using Loci::AMRRemapPlan ;
using Loci::AMRRemapReport ;
using Loci::CPTR ;
using Loci::vector3d ;

namespace {

  // These geometry fixtures use matching transition indices and persistent IDs.
  std::vector<AMRCellGeometry> source_geometry() {
    return std::vector<AMRCellGeometry>{
      AMRCellGeometry(0,0,1.0,vector3d<double>(0.0,0.0,0.0)),
      AMRCellGeometry(1,1,0.25,vector3d<double>(1.0,0.0,0.0)),
      AMRCellGeometry(2,2,0.75,vector3d<double>(2.0,0.0,0.0)),
      AMRCellGeometry(3,3,2.0,vector3d<double>(3.0,0.0,0.0))
    } ;
  }

  std::vector<AMRCellGeometry> target_geometry() {
    return std::vector<AMRCellGeometry>{
      AMRCellGeometry(10,10,0.25,vector3d<double>(-0.75,0.0,0.0)),
      AMRCellGeometry(11,11,0.75,vector3d<double>(0.25,0.0,0.0)),
      AMRCellGeometry(12,12,1.0,vector3d<double>(1.75,0.0,0.0)),
      AMRCellGeometry(13,13,2.0,vector3d<double>(3.0,0.0,0.0))
    } ;
  }

  std::vector<AMRCellContribution> valid_contributions() {
    return std::vector<AMRCellContribution>{
      AMRCellContribution(0,10,0,10,0.25,vector3d<double>(-0.75,0.0,0.0)),
      AMRCellContribution(0,11,0,11,0.75,vector3d<double>(0.25,0.0,0.0)),
      AMRCellContribution(1,12,1,12,0.25,vector3d<double>(1.0,0.0,0.0)),
      AMRCellContribution(2,12,2,12,0.75,vector3d<double>(2.0,0.0,0.0)),
      AMRCellContribution(3,13,3,13,2.0,vector3d<double>(3.0,0.0,0.0))
    } ;
  }

  bool has_plan(const CPTR<AMRRemapPlan>& plan) {
    return plan != static_cast<AMRRemapPlan*>(0) ;
  }

  bool has_plan(const Loci::const_CPTR<AMRRemapPlan>& plan) {
    return plan != static_cast<const AMRRemapPlan*>(0) ;
  }

  Loci::fact_db::distribute_infoP
  make_serial_distribution(const Loci::entitySet& owned) {
    Loci::fact_db::distribute_infoP dist =
      new Loci::fact_db::distribute_info ;
    dist->myid = 0 ;
    dist->isDistributed = 1 ;
    dist->my_entities = owned ;
    dist->comp_entities = owned ;
    dist->copy_total_size = 0 ;
    dist->xmit_total_size = 0 ;

    Loci::Map l2g ;
    Loci::Map l2f ;
    Loci::store<unsigned char> keyDomain ;
    l2g.allocate(owned) ;
    l2f.allocate(owned) ;
    keyDomain.allocate(owned) ;
    FORALL(owned,cell) {
      l2g[cell] = cell ;
      l2f[cell] = cell ;
      keyDomain[cell] = 0 ;
    } ENDFORALL ;
    dist->l2g = l2g.Rep() ;
    dist->l2f = l2f.Rep() ;
    dist->key_domain = keyDomain.Rep() ;
    return dist ;
  }

  void setup_runtime_mapping(Loci::AMRinterpolation& interpolation) {
    const Loci::entitySet sources = Loci::interval(0,3) ;
    const Loci::entitySet targets = Loci::interval(10,13) ;

    Loci::store<Loci::CellId> sourceIds, targetIds ;
    sourceIds.allocate(sources) ;
    targetIds.allocate(targets) ;
    FORALL(sources,cell) { sourceIds[cell] = cell ; } ENDFORALL ;
    FORALL(targets,cell) { targetIds[cell] = cell ; } ENDFORALL ;

    Loci::store<std::pair<int,int> > targetSource ;
    targetSource.allocate(Loci::interval(0,4)) ;
    targetSource[0] = std::make_pair(10,0) ;
    targetSource[1] = std::make_pair(11,0) ;
    targetSource[2] = std::make_pair(12,1) ;
    targetSource[3] = std::make_pair(12,2) ;
    targetSource[4] = std::make_pair(13,3) ;

    Loci::store<double> sourceVolume ;
    Loci::store<vector3d<double> > sourceCenter ;
    sourceVolume.allocate(sources) ;
    sourceCenter.allocate(sources) ;
    const std::vector<AMRCellGeometry> sourceData = source_geometry() ;
    for(size_t source=0;source<sourceData.size();++source) {
      sourceVolume[sourceData[source].cell] = sourceData[source].volume ;
      sourceCenter[sourceData[source].cell] = sourceData[source].centroid ;
    }

    Loci::store<int> stencilSizes ;
    stencilSizes.allocate(sources) ;
    FORALL(sources,source) {
      stencilSizes[source] = source == 0 ? 2 : 0 ;
    } ENDFORALL ;
    Loci::multiStore<int> gradientCells ;
    Loci::multiStore<vector3d<double> > gradientDeltas ;
    gradientCells.allocate(stencilSizes) ;
    gradientDeltas.allocate(stencilSizes) ;
    gradientCells[0][0] = 1 ;
    gradientCells[0][1] = 2 ;
    gradientDeltas[0][0] = vector3d<double>(1.0,0.0,0.0) ;
    gradientDeltas[0][1] = vector3d<double>(2.0,0.0,0.0) ;

    Loci::store<double> targetVolume ;
    Loci::store<vector3d<double> > targetCenter ;
    targetVolume.allocate(targets) ;
    targetCenter.allocate(targets) ;
    const std::vector<AMRCellGeometry> targetData = target_geometry() ;
    for(size_t target=0;target<targetData.size();++target) {
      targetVolume[targetData[target].cell] = targetData[target].volume ;
      targetCenter[targetData[target].cell] = targetData[target].centroid ;
    }

    Loci::fact_db facts ;
    Loci::constraint geomCells ;
    geomCells = targets ;
    facts.create_fact("geom_cells",geomCells) ;
    facts.put_distribute_info(make_serial_distribution(targets)) ;

    Loci::const_store<vector3d<double> > sourceCenterConst(
      sourceCenter.Rep()) ;
    Loci::const_store<vector3d<double> > targetCenterConst(
      targetCenter.Rep()) ;
    Loci::const_store<double> targetVolumeConst(targetVolume.Rep()) ;
    Loci::const_store<Loci::CellId> sourceIdsConst(sourceIds.Rep()) ;
    Loci::const_store<Loci::CellId> targetIdsConst(targetIds.Rep()) ;
    interpolation.setupRefinementMapping(
      targetSource,sourceVolume,sourceCenterConst,sourceIdsConst,gradientCells,
      gradientDeltas,targetCenterConst,targetVolumeConst,targetIdsConst,facts) ;
  }

} // namespace


/// One contribution model should describe retained, refined, and derefined
/// cells while preserving constant, linear, and extensive cell data.
TEST_CASE("cell remap contributions preserve geometric transfer behavior") {
  AMRRemapReport report ;
  CPTR<AMRRemapPlan> plan = AMRRemapPlan::createCellPlan(
    source_geometry(),target_geometry(),valid_contributions(),report) ;

  REQUIRE(has_plan(plan)) ;
  REQUIRE(report.valid) ;
  CHECK(report.maximumSourceVolumeError == doctest::Approx(0.0)) ;
  CHECK(report.maximumTargetVolumeError == doctest::Approx(0.0)) ;

  Loci::amr_cell_transition::value transition ;
  REQUIRE(plan->cellTransition(10,transition)) ;
  CHECK(transition == Loci::amr_cell_transition::refined) ;
  REQUIRE(plan->cellTransition(12,transition)) ;
  CHECK(transition == Loci::amr_cell_transition::derefined) ;
  REQUIRE(plan->cellTransition(13,transition)) ;
  CHECK(transition == Loci::amr_cell_transition::retained) ;

  std::vector<double> targetValues ;
  REQUIRE(plan->remapCellAverages(std::vector<double>(4,7.0),targetValues)) ;
  CHECK(targetValues == std::vector<double>(4,7.0)) ;

  std::vector<double> sourceValues{0.0,1.0,2.0,3.0} ;
  std::vector<vector3d<double> > gradients(
    4,vector3d<double>(1.0,0.0,0.0)) ;
  REQUIRE(plan->remapCellAverages(sourceValues,gradients,targetValues)) ;
  const std::vector<AMRCellGeometry> targets = plan->targetCellGeometry() ;
  REQUIRE(targetValues.size() == targets.size()) ;
  for(size_t target=0;target<targets.size();++target)
    CHECK(targetValues[target] == doctest::Approx(targets[target].centroid.x)) ;

  std::vector<double> targetIntegrals ;
  const std::vector<double> sourceIntegrals{2.0,0.5,1.5,4.0} ;
  REQUIRE(plan->remapCellIntegrals(sourceIntegrals,targetIntegrals)) ;
  CHECK(std::accumulate(targetIntegrals.begin(),targetIntegrals.end(),0.0) ==
        doctest::Approx(8.0)) ;

  SUBCASE("a target-owned partition retains global source cardinality") {
    const std::vector<AMRCellGeometry> localSources{
      AMRCellGeometry(0,0,1.0,vector3d<double>(0.0,0.0,0.0))
    } ;
    const std::vector<AMRCellGeometry> localTargets{
      AMRCellGeometry(10,10,0.25,vector3d<double>(-0.75,0.0,0.0))
    } ;
    const std::vector<AMRCellContribution> localContributions{
      AMRCellContribution(0,10,0,10,0.25,vector3d<double>(-0.75,0.0,0.0))
    } ;
    CPTR<AMRRemapPlan> localPlan = AMRRemapPlan::createCellPlanPartition(
      localSources,std::vector<size_t>(1,2),localTargets,
      localContributions,report) ;
    REQUIRE(has_plan(localPlan)) ;
    CHECK_FALSE(report.sourceCoverageChecked) ;
    REQUIRE(localPlan->cellTransition(10,transition)) ;
    CHECK(transition == Loci::amr_cell_transition::refined) ;
  }

  SUBCASE("the runtime interpolation setup publishes the same relation") {
    REQUIRE(Loci::MPI_processes == 1) ;
    Loci::AMRinterpolation interpolation ;
    setup_runtime_mapping(interpolation) ;
    const Loci::const_CPTR<AMRRemapPlan> runtimePlan =
      interpolation.getRemapPlan() ;
    REQUIRE(has_plan(runtimePlan)) ;
    CHECK(interpolation.getRemapReport().valid) ;
    CHECK(interpolation.getRemapReport().sourceCoverageChecked) ;
    CHECK(interpolation.getRemapReport().sourceCells == 4) ;
    CHECK(interpolation.getRemapReport().targetCells == 4) ;
    CHECK(interpolation.getRemapReport().contributions == 5) ;
    REQUIRE(runtimePlan->cellTransition(12,transition)) ;
    CHECK(transition == Loci::amr_cell_transition::derefined) ;

    Loci::store<double> sourceAverages ;
    sourceAverages.allocate(Loci::interval(0,3)) ;
    FORALL(sourceAverages.domain(),source) {
      sourceAverages[source] = 7.0 ;
    } ENDFORALL ;
    Loci::store<double> targetAverages ;
    REQUIRE(interpolation.remapCellAverages(sourceAverages,targetAverages)) ;
    FORALL(targetAverages.domain(),target) {
      CHECK(targetAverages[target] == doctest::Approx(7.0)) ;
    } ENDFORALL ;

    Loci::store<double> sourceIntegrals ;
    sourceIntegrals.allocate(Loci::interval(0,3)) ;
    sourceIntegrals[0] = 2.0 ;
    sourceIntegrals[1] = 0.5 ;
    sourceIntegrals[2] = 1.5 ;
    sourceIntegrals[3] = 4.0 ;
    Loci::store<double> targetIntegrals ;
    REQUIRE(interpolation.remapCellIntegrals(
      sourceIntegrals,targetIntegrals)) ;
    CHECK(targetIntegrals[10] == doctest::Approx(0.5)) ;
    CHECK(targetIntegrals[11] == doctest::Approx(1.5)) ;
    CHECK(targetIntegrals[12] == doctest::Approx(2.0)) ;
    CHECK(targetIntegrals[13] == doctest::Approx(4.0)) ;
  }
}


/// Invalid coverage, duplicate contributions, and many-to-many relations must
/// be rejected before a solver can use the remap plan.
TEST_CASE("cell remap plans reject ambiguous or incomplete relations") {
  AMRRemapReport report ;

  SUBCASE("incomplete target coverage") {
    std::vector<AMRCellContribution> contributions = valid_contributions() ;
    contributions[0].overlapVolume = 0.20 ;
    CHECK_FALSE(has_plan(AMRRemapPlan::createCellPlan(
      source_geometry(),target_geometry(),contributions,report))) ;
    CHECK_FALSE(report.valid) ;
    CHECK(report.missingSourceCells > 0) ;
    CHECK(report.missingTargetCells > 0) ;
  }

  SUBCASE("duplicate source-target pair") {
    std::vector<AMRCellContribution> contributions = valid_contributions() ;
    contributions.push_back(contributions[0]) ;
    CHECK_FALSE(has_plan(AMRRemapPlan::createCellPlan(
      source_geometry(),target_geometry(),contributions,report))) ;
    CHECK(report.duplicateContributions == 1) ;
  }

  SUBCASE("many-to-many relation") {
    const std::vector<AMRCellGeometry> sources{
      AMRCellGeometry(0,0,1.0,vector3d<double>(0.5,0.0,0.0)),
      AMRCellGeometry(1,1,1.0,vector3d<double>(1.5,0.0,0.0))
    } ;
    const std::vector<AMRCellGeometry> targets{
      AMRCellGeometry(10,10,1.0,vector3d<double>(0.5,0.0,0.0)),
      AMRCellGeometry(11,11,1.0,vector3d<double>(1.5,0.0,0.0))
    } ;
    const std::vector<AMRCellContribution> contributions{
      AMRCellContribution(0,10,0,10,0.5,vector3d<double>(0.25,0.0,0.0)),
      AMRCellContribution(0,11,0,11,0.5,vector3d<double>(0.75,0.0,0.0)),
      AMRCellContribution(1,10,1,10,0.5,vector3d<double>(1.25,0.0,0.0)),
      AMRCellContribution(1,11,1,11,0.5,vector3d<double>(1.75,0.0,0.0))
    } ;
    CHECK_FALSE(has_plan(AMRRemapPlan::createCellPlan(
      sources,targets,contributions,report))) ;
    CHECK(report.unsupportedRelations == 4) ;
  }
}


int main(int argc, char **argv) {
  Loci::Init(&argc,&argv) ;

  doctest::Context context ;
  context.applyCommandLine(argc,argv) ;
  const int result = context.run() ;

  Loci::Finalize() ;
  if(result == 0)
    std::cout << "SUCCESS!" << std::endl ;
  return result ;
}

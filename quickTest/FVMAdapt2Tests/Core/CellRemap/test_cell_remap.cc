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
#include <FVMAdapt2/remap_plan.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <numeric>
#include <vector>

using namespace Loci ;

namespace {

  const CellId sourceRetained = (CellId(1) << 48) + 1 ;
  const CellId sourceRefined = (CellId(1) << 48) + 2 ;
  const CellId sourceDerefinedLeft = (CellId(1) << 48) + 3 ;
  const CellId sourceDerefinedRight = (CellId(1) << 48) + 4 ;

  const CellId targetRetained = (CellId(1) << 52) + 1 ;
  const CellId targetRefinedLeft = (CellId(1) << 52) + 2 ;
  const CellId targetRefinedRight = (CellId(1) << 52) + 3 ;
  const CellId targetDerefined = (CellId(1) << 52) + 4 ;

  CPTR<AMRRemapPlan> make_transition_plan(AMRRemapReport& report) {
    const std::vector<AMRCellGeometry> sourceGeometry{
          AMRCellGeometry(
                10, sourceRetained, 2.0, vector3d<double>(0.0, 0.0, 0.0)),
          AMRCellGeometry(
                20, sourceRefined, 2.0, vector3d<double>(2.0, 0.0, 0.0)),
          AMRCellGeometry(
                30, sourceDerefinedLeft, 1.0, vector3d<double>(4.5, 0.0, 0.0)),
          AMRCellGeometry(40, sourceDerefinedRight, 1.0,
                vector3d<double>(5.5, 0.0, 0.0))} ;
    const std::vector<AMRCellGeometry> targetGeometry{
          AMRCellGeometry(
                110, targetRetained, 2.0, vector3d<double>(0.0, 0.0, 0.0)),
          AMRCellGeometry(
                120, targetRefinedLeft, 1.0, vector3d<double>(1.5, 0.0, 0.0)),
          AMRCellGeometry(
                130, targetRefinedRight, 1.0, vector3d<double>(2.5, 0.0, 0.0)),
          AMRCellGeometry(
                140, targetDerefined, 2.0, vector3d<double>(5.0, 0.0, 0.0))} ;
    const std::vector<AMRCellContribution> contributions{
          AMRCellContribution(10, 110, sourceRetained, targetRetained, 2.0,
                vector3d<double>(0.0, 0.0, 0.0)),
          AMRCellContribution(20, 120, sourceRefined, targetRefinedLeft, 1.0,
                vector3d<double>(1.5, 0.0, 0.0)),
          AMRCellContribution(20, 130, sourceRefined, targetRefinedRight, 1.0,
                vector3d<double>(2.5, 0.0, 0.0)),
          AMRCellContribution(30, 140, sourceDerefinedLeft, targetDerefined,
                1.0, vector3d<double>(4.5, 0.0, 0.0)),
          AMRCellContribution(40, 140, sourceDerefinedRight, targetDerefined,
                1.0, vector3d<double>(5.5, 0.0, 0.0))} ;
    return AMRRemapPlan::createCellPlan(
          sourceGeometry, targetGeometry, contributions, report) ;
  }

  bool has_plan(const CPTR<AMRRemapPlan>& plan) {
    return plan != static_cast<AMRRemapPlan*>(0) ;
  }

  void check_transition(const AMRRemapPlan& plan, CellId target,
        amr_cell_transition::value expected) {
    amr_cell_transition::value actual = amr_cell_transition::retained ;
    REQUIRE(plan.cellTransition(target, actual)) ;
    CHECK(actual == expected) ;
  }

} // namespace

/// Persistent identities select retained, refined, and derefined cell
/// relations without depending on transition-local entity numbers.
TEST_CASE("cell remap exposes transition relations by persistent identity") {
  AMRRemapReport report ;
  CPTR<AMRRemapPlan> plan = make_transition_plan(report) ;
  REQUIRE(has_plan(plan)) ;
  REQUIRE(report.valid) ;

  check_transition(*plan, targetRetained, amr_cell_transition::retained) ;
  check_transition(*plan, targetRefinedLeft, amr_cell_transition::refined) ;
  check_transition(*plan, targetRefinedRight, amr_cell_transition::refined) ;
  check_transition(*plan, targetDerefined, amr_cell_transition::derefined) ;

  size_t begin = 0 ;
  size_t end = 0 ;
  REQUIRE(plan->cellContributions(targetDerefined, begin, end)) ;
  REQUIRE(end - begin == 2) ;
  CHECK(plan->cellContributions()[begin].sourceCellId == sourceDerefinedLeft) ;
  CHECK(plan->cellContributions()[begin + 1].sourceCellId ==
        sourceDerefinedRight) ;

  std::vector<CellId> targets ;
  plan->targetCells(sourceRefined, targets) ;
  CHECK(targets == std::vector<CellId>{targetRefinedLeft, targetRefinedRight}) ;
  REQUIRE(plan->sourceCellGeometry(sourceRefined) != 0) ;
  CHECK(plan->sourceCellGeometry(sourceRefined)->cell == 20) ;
  REQUIRE(plan->targetCellGeometry(targetDerefined) != 0) ;
  CHECK(plan->targetCellGeometry(targetDerefined)->cell == 140) ;
}

/// Cell averages and extensive integrals follow the same overlap relation and
/// preserve the global integral across refinement and derefinement.
TEST_CASE("cell remap conserves average and integral transfers") {
  AMRRemapReport report ;
  CPTR<AMRRemapPlan> plan = make_transition_plan(report) ;
  REQUIRE(has_plan(plan)) ;

  const std::vector<double> sourceAverages{3.0, 5.0, 7.0, 11.0} ;
  std::vector<double> targetAverages ;
  REQUIRE(plan->remapCellAverages(sourceAverages, targetAverages)) ;
  REQUIRE(targetAverages.size() == 4) ;
  CHECK(targetAverages[0] == doctest::Approx(3.0)) ;
  CHECK(targetAverages[1] == doctest::Approx(5.0)) ;
  CHECK(targetAverages[2] == doctest::Approx(5.0)) ;
  CHECK(targetAverages[3] == doctest::Approx(9.0)) ;

  const std::vector<double> targetVolumes{2.0, 1.0, 1.0, 2.0} ;
  double targetAverageIntegral = 0.0 ;
  for (size_t target = 0; target < targetAverages.size(); ++target)
    targetAverageIntegral += targetAverages[target] * targetVolumes[target] ;
  CHECK(targetAverageIntegral == doctest::Approx(34.0)) ;

  const std::vector<double> sourceIntegrals{6.0, 10.0, 7.0, 11.0} ;
  std::vector<double> targetIntegrals ;
  REQUIRE(plan->remapCellIntegrals(sourceIntegrals, targetIntegrals)) ;
  REQUIRE(targetIntegrals.size() == 4) ;
  CHECK(targetIntegrals[0] == doctest::Approx(6.0)) ;
  CHECK(targetIntegrals[1] == doctest::Approx(5.0)) ;
  CHECK(targetIntegrals[2] == doctest::Approx(5.0)) ;
  CHECK(targetIntegrals[3] == doctest::Approx(18.0)) ;
  CHECK(std::accumulate(targetIntegrals.begin(), targetIntegrals.end(), 0.0) ==
        doctest::Approx(std::accumulate(
              sourceIntegrals.begin(), sourceIntegrals.end(), 0.0))) ;
}

/// A persistent cell identity must name one source cell; duplicate identities
/// are rejected even when transition-local entity numbers are distinct.
TEST_CASE("cell remap rejects duplicate persistent source identities") {
  const CellId duplicateSource = (CellId(1) << 56) + 1 ;
  const CellId firstTarget = (CellId(1) << 57) + 1 ;
  const CellId secondTarget = (CellId(1) << 57) + 2 ;
  const std::vector<AMRCellGeometry> sourceGeometry{
        AMRCellGeometry(
              1, duplicateSource, 1.0, vector3d<double>(0.0, 0.0, 0.0)),
        AMRCellGeometry(
              2, duplicateSource, 1.0, vector3d<double>(1.0, 0.0, 0.0))} ;
  const std::vector<AMRCellGeometry> targetGeometry{
        AMRCellGeometry(3, firstTarget, 1.0, vector3d<double>(0.0, 0.0, 0.0)),
        AMRCellGeometry(4, secondTarget, 1.0, vector3d<double>(1.0, 0.0, 0.0))} ;
  const std::vector<AMRCellContribution> contributions{
        AMRCellContribution(1, 3, duplicateSource, firstTarget, 1.0,
              vector3d<double>(0.0, 0.0, 0.0)),
        AMRCellContribution(2, 4, duplicateSource, secondTarget, 1.0,
              vector3d<double>(1.0, 0.0, 0.0))} ;

  AMRRemapReport report ;
  CPTR<AMRRemapPlan> plan = AMRRemapPlan::createCellPlan(
        sourceGeometry, targetGeometry, contributions, report) ;
  CHECK_FALSE(has_plan(plan)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.invalidIdentities == 1) ;
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv) ;
  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;
  Loci::Finalize() ;
  return result ;
}

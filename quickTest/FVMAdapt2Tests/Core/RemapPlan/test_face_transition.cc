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
#include "mesh_state.h"
#include <FVMAdapt2/face.h>
#include <FVMAdapt2/prism.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <algorithm>
#include <map>
#include <queue>
#include <utility>
#include <vector>

using namespace Loci ;

namespace {

  typedef std::vector<vector3d<double>> Polygon ;

  typedef std::pair<Prism*, std::vector<int>> PrismPathNode ;
  typedef std::pair<Face*, std::vector<int>> FacePathNode ;

  Polygon square(double xmin, double xmax, double ymin, double ymax,
        bool reversed = false) {
    Polygon result{vector3d<double>(xmin, ymin, 0.0),
          vector3d<double>(xmax, ymin, 0.0), vector3d<double>(xmax, ymax, 0.0),
          vector3d<double>(xmin, ymax, 0.0)} ;
    if (reversed)
      std::reverse(result.begin(), result.end()) ;
    return result ;
  }

  Polygon interior_square() {
    return Polygon{vector3d<double>(0.5, 0.0, 0.0),
          vector3d<double>(0.5, 0.0, 1.0), vector3d<double>(0.5, 1.0, 1.0),
          vector3d<double>(0.5, 1.0, 0.0)} ;
  }

  Polygon interior_patch(double ymin, double ymax, double zmin, double zmax) {
    return Polygon{vector3d<double>(0.5, ymin, zmin),
          vector3d<double>(0.5, ymin, zmax), vector3d<double>(0.5, ymax, zmax),
          vector3d<double>(0.5, ymax, zmin)} ;
  }

  Polygon child_interior_square() {
    return Polygon{vector3d<double>(0.25, 0.0, 0.0),
          vector3d<double>(0.25, 0.0, 1.0), vector3d<double>(0.25, 1.0, 1.0),
          vector3d<double>(0.25, 1.0, 0.0)} ;
  }

  bool prism_replay_paths(
        const std::vector<char>& plan, std::vector<std::vector<int>>& paths) {
    Prism root ;
    const int leafCount = root.empty_resplit(plan) ;
    if (leafCount < 1)
      return false ;
    paths.assign(size_t(leafCount), std::vector<int>()) ;
    std::vector<bool> assigned(size_t(leafCount), false) ;
    std::queue<PrismPathNode> pending ;
    pending.push(PrismPathNode(&root, std::vector<int>())) ;
    while (!pending.empty()) {
      Prism* cell = pending.front().first ;
      const std::vector<int> path = pending.front().second ;
      pending.pop() ;
      const int children = cell->numChildren() ;
      if (children == 0) {
        const int leaf = cell->getCellIndex() - 1 ;
        if (leaf < 0 || leaf >= leafCount || assigned[size_t(leaf)])
          return false ;
        paths[size_t(leaf)] = path ;
        assigned[size_t(leaf)] = true ;
        continue ;
      }
      if (children < 0)
        return false ;
      const int splitCode = int(cell->getMySplitCode()) ;
      for (int child = 0; child < children; ++child) {
        std::vector<int> childPath = path ;
        childPath.push_back(splitCode) ;
        childPath.push_back(child) ;
        pending.push(PrismPathNode(cell->getChildCell(child), childPath)) ;
      }
    }
    return std::find(assigned.begin(), assigned.end(), false) == assigned.end() ;
  }

  bool general_face_replay_paths(const std::vector<char>& plan, int edgeCount,
        std::vector<std::vector<int>>& paths) {
    Face root(edgeCount) ;
    std::vector<Face*> leaves ;
    root.empty_resplit(plan, leaves) ;
    std::map<Face*, std::vector<int>> leafPaths ;
    std::queue<FacePathNode> pending ;
    pending.push(FacePathNode(&root, std::vector<int>())) ;
    while (!pending.empty()) {
      Face* face = pending.front().first ;
      const std::vector<int> path = pending.front().second ;
      pending.pop() ;
      if (face->child == 0) {
        leafPaths[face] = path ;
        continue ;
      }
      for (int child = 0; child < face->numEdge; ++child) {
        std::vector<int> childPath = path ;
        childPath.push_back(1) ;
        childPath.push_back(child) ;
        pending.push(FacePathNode(face->child[child], childPath)) ;
      }
    }
    paths.clear() ;
    for (size_t leaf = 0; leaf < leaves.size(); ++leaf) {
      const std::map<Face*, std::vector<int>>::const_iterator path =
            leafPaths.find(leaves[leaf]) ;
      if (path == leafPaths.end())
        return false ;
      paths.push_back(path->second) ;
    }
    return paths.size() == leafPaths.size() ;
  }

  bool cell_leaf_paths(const std::vector<char>& plan,
        cell_topology::value topology, std::vector<std::vector<int>>& paths) {
    std::vector<int> encoded ;
    const bool encodedPlan =
          topology == cell_topology::hex
                ? detail::encodeHexCellLeafPaths(plan, encoded)
                : (topology == cell_topology::prism
                              ? detail::encodePrismCellLeafPaths(plan, encoded)
                              : false) ;
    return encodedPlan && detail::decodeLeafPaths(encoded, paths) ;
  }

  RootCellState root_state_from_plan(int root, const std::vector<char>& plan,
        cell_topology::value topology = cell_topology::hex) {
    std::vector<std::vector<int>> paths ;
    cell_leaf_paths(plan, topology, paths) ;
    return RootCellState(root, paths, topology) ;
  }

  CPTR<FaceState> coarse_state(FaceTransitionReport& report) {
    const std::vector<FaceIdentity> identities{FaceIdentity(
          10, FaceKey(face_origin::base_face, 10, std::vector<int>()))} ;
    const std::vector<Polygon> polygons{square(0.0, 1.0, 0.0, 1.0)} ;
    const std::vector<RootCellState> roots{
          root_state_from_plan(42, std::vector<char>())} ;
    return FaceState::create(identities, polygons, roots, report) ;
  }

  CPTR<FaceState> refined_state(FaceTransitionReport& report) {
    const std::vector<FaceIdentity> identities{
          FaceIdentity(100, FaceKey(face_origin::base_face, 10, {3, 0})),
          FaceIdentity(101, FaceKey(face_origin::base_face, 10, {3, 1})),
          FaceIdentity(102, FaceKey(face_origin::base_face, 10, {3, 2})),
          FaceIdentity(103, FaceKey(face_origin::base_face, 10, {3, 3})),
          FaceIdentity(
                104, FaceKey(face_origin::cell_interior, 42, {7, 0}, {7, 1}))} ;
    const std::vector<Polygon> polygons{square(0.0, 0.5, 0.0, 0.5),
          square(0.5, 1.0, 0.0, 0.5, true), square(0.0, 0.5, 0.5, 1.0),
          square(0.5, 1.0, 0.5, 1.0), interior_square()} ;
    const std::vector<RootCellState> roots{
          root_state_from_plan(42, std::vector<char>{7, 0})} ;
    return FaceState::create(identities, polygons, roots, report) ;
  }

  CPTR<FaceState> first_level_state(FaceTransitionReport& report) {
    const std::vector<FaceIdentity> identities{
          FaceIdentity(100, FaceKey(face_origin::base_face, 10, {3, 0})),
          FaceIdentity(
                104, FaceKey(face_origin::cell_interior, 42, {7, 0}, {7, 1}))} ;
    const std::vector<Polygon> polygons{
          square(0.0, 1.0, 0.0, 1.0), interior_square()} ;
    const std::vector<RootCellState> roots{
          root_state_from_plan(42, std::vector<char>{7, 0})} ;
    return FaceState::create(identities, polygons, roots, report) ;
  }

  CPTR<FaceState> second_level_state(FaceTransitionReport& report) {
    const std::vector<FaceIdentity> identities{
          FaceIdentity(200, FaceKey(face_origin::base_face, 10, {3, 0})),
          FaceIdentity(204,
                FaceKey(face_origin::cell_interior, 42, {7, 0, 7, 0}, {7, 1})),
          FaceIdentity(205,
                FaceKey(face_origin::cell_interior, 42, {7, 0, 7, 1}, {7, 1})),
          FaceIdentity(206,
                FaceKey(face_origin::cell_interior, 42, {7, 0, 7, 2}, {7, 1})),
          FaceIdentity(207,
                FaceKey(face_origin::cell_interior, 42, {7, 0, 7, 3}, {7, 1})),
          FaceIdentity(208, FaceKey(face_origin::cell_interior, 42,
                                  {7, 0, 7, 0}, {7, 0, 7, 1}))} ;
    const std::vector<Polygon> polygons{square(0.0, 1.0, 0.0, 1.0),
          interior_patch(0.0, 0.5, 0.0, 0.5),
          interior_patch(0.5, 1.0, 0.0, 0.5),
          interior_patch(0.0, 0.5, 0.5, 1.0),
          interior_patch(0.5, 1.0, 0.5, 1.0), child_interior_square()} ;
    const std::vector<RootCellState> roots{
          root_state_from_plan(42, std::vector<char>{7, 7})} ;
    return FaceState::create(identities, polygons, roots, report) ;
  }

  CPTR<FaceState> prism_state(const std::vector<char>& plan,
        const std::vector<FaceIdentity>& identities,
        const std::vector<Polygon>& polygons, FaceTransitionReport& report) {
    const std::vector<RootCellState> roots{
          root_state_from_plan(52, plan, cell_topology::prism)} ;
    return FaceState::create(identities, polygons, roots, report) ;
  }

  bool has_state(const CPTR<FaceState>& state) {
    return state != static_cast<FaceState*>(0) ;
  }

  bool has_transition(const CPTR<FaceRemap>& transition) {
    return transition != static_cast<FaceRemap*>(0) ;
  }

  FaceId face_id(face_origin::value origin, int root,
        const std::vector<int>& first,
        const std::vector<int>& second = std::vector<int>()) {
    return persistentFaceId(FaceKey(origin, root, first, second)) ;
  }

} // namespace

/// A first refinement must associate each inherited face piece with its source
/// and identify a new cell-interior face without inventing a source face.
TEST_CASE("face transition reports split and created faces") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = coarse_state(stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = refined_state(stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.remap.contributions == 4) ;
  CHECK(report.remap.createdFaces == 1) ;
  CHECK(report.remap.removedFaces == 0) ;

  const_CPTR<FaceRemap> remap = transition ;
  std::vector<FaceId> targetFaces ;
  remap->targetFaces(face_id(face_origin::base_face, 10, {}), targetFaces) ;
  std::vector<FaceId> expectedTargetFaces{
        face_id(face_origin::base_face, 10, {3, 0}),
        face_id(face_origin::base_face, 10, {3, 1}),
        face_id(face_origin::base_face, 10, {3, 2}),
        face_id(face_origin::base_face, 10, {3, 3})} ;
  std::sort(expectedTargetFaces.begin(), expectedTargetFaces.end()) ;
  CHECK(targetFaces == expectedTargetFaces) ;

  CellId sourceCell = 0 ;
  REQUIRE(remap->isCreatedFace(
        face_id(face_origin::cell_interior, 42, {7, 0}, {7, 1}), sourceCell)) ;
  CHECK(sourceCell == persistentCellId(42, {})) ;
  remap->targetFaces(999, targetFaces) ;
  CHECK(targetFaces.empty()) ;

  std::vector<double> targetValues ;
  std::vector<unsigned char> mapped ;
  REQUIRE(remap->remapFaceAverages(
        std::vector<double>(1, 5.0), targetValues, mapped)) ;
  REQUIRE(targetValues.size() == 5) ;
  const FaceId createdFace =
        face_id(face_origin::cell_interior, 42, {7, 0}, {7, 1}) ;
  const FaceId reversedFace = face_id(face_origin::base_face, 10, {3, 1}) ;
  size_t reversedIndex = targetValues.size() ;
  for (size_t target = 0; target < remap->targetFaceGeometry().size();
        ++target) {
    const FaceId face = remap->targetFaceGeometry()[target].face ;
    if (face == createdFace) {
      CHECK(targetValues[target] == doctest::Approx(0.0)) ;
      CHECK(mapped[target] == 0) ;
    } else {
      CHECK(targetValues[target] == doctest::Approx(5.0)) ;
      CHECK(mapped[target] == 1) ;
    }
    if (face == reversedFace)
      reversedIndex = target ;
  }

  // Face averages are unsigned by default even when mesh orientation reverses.
  REQUIRE(reversedIndex < targetValues.size()) ;
  bool foundReversedContribution = false ;
  for (size_t contribution = 0; contribution < remap->overlaps().size();
        ++contribution) {
    if (remap->overlaps()[contribution].target == reversedFace) {
      CHECK(remap->overlaps()[contribution].orientation == -1) ;
      foundReversedContribution = true ;
    }
  }
  REQUIRE(foundReversedContribution) ;
  REQUIRE(remap->remapFaceAverages(
        std::vector<double>(1, 5.0), targetValues, mapped, true)) ;
  CHECK(targetValues[reversedIndex] == doctest::Approx(-5.0)) ;
}

/// Full derefinement must merge inherited face pieces and identify internal
/// source faces that disappear into the accepted coarse cell.
TEST_CASE("face transition reports coarsened and removed faces") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = refined_state(stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = coarse_state(stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.remap.contributions == 4) ;
  CHECK(report.remap.createdFaces == 0) ;
  CHECK(report.remap.removedFaces == 1) ;

  const_CPTR<FaceRemap> remap = transition ;
  for (int child = 0; child < 4; ++child) {
    std::vector<FaceId> targetFaces ;
    remap->targetFaces(
          face_id(face_origin::base_face, 10, {3, child}), targetFaces) ;
    CHECK(targetFaces ==
          std::vector<FaceId>{face_id(face_origin::base_face, 10, {})}) ;
  }

  std::vector<FaceId> targetFaces{face_id(face_origin::base_face, 10, {})} ;
  const FaceId removedFace =
        face_id(face_origin::cell_interior, 42, {7, 0}, {7, 1}) ;
  remap->targetFaces(removedFace, targetFaces) ;
  CHECK(targetFaces.empty()) ;

  CellId targetCell = 0 ;
  REQUIRE(remap->isRemovedFace(removedFace, targetCell)) ;
  CHECK(targetCell == persistentCellId(42, {})) ;
}

/// Refining one existing leaf must retain inherited faces and associate each
/// newly created internal face with the source leaf cell that contained it.
TEST_CASE("face transition supports incremental refinement") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = first_level_state(stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = second_level_state(stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.unsupportedRootPlans == 0) ;
  CHECK(report.remap.contributions == 5) ;
  CHECK(report.remap.createdFaces == 1) ;
  CHECK(report.remap.removedFaces == 0) ;

  const_CPTR<FaceRemap> remap = transition ;
  std::vector<FaceId> targetFaces ;
  remap->targetFaces(
        face_id(face_origin::cell_interior, 42, {7, 0}, {7, 1}), targetFaces) ;
  std::vector<FaceId> expectedTargetFaces{
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 0}, {7, 1}),
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 1}, {7, 1}),
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 2}, {7, 1}),
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 3}, {7, 1})} ;
  std::sort(expectedTargetFaces.begin(), expectedTargetFaces.end()) ;
  CHECK(targetFaces == expectedTargetFaces) ;

  CellId sourceCell = 0 ;
  REQUIRE(remap->isCreatedFace(
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 0}, {7, 0, 7, 1}),
        sourceCell)) ;
  CHECK(sourceCell == persistentCellId(42, {7, 0})) ;
}

/// Derefining one branch of a nonempty plan must merge inherited faces and
/// associate disappearing internal faces with the accepted containing cell.
TEST_CASE("face transition supports incremental derefinement") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = second_level_state(stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = first_level_state(stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.unsupportedRootPlans == 0) ;
  CHECK(report.remap.contributions == 5) ;
  CHECK(report.remap.createdFaces == 0) ;
  CHECK(report.remap.removedFaces == 1) ;

  const_CPTR<FaceRemap> remap = transition ;
  CellId targetCell = 0 ;
  REQUIRE(remap->isRemovedFace(
        face_id(face_origin::cell_interior, 42, {7, 0, 7, 0}, {7, 0, 7, 1}),
        targetCell)) ;
  CHECK(targetCell == persistentCellId(42, {7, 0})) ;
}

/// Prism split codes must use the prism tree when locating the source leaf
/// that contains a newly created internal face.
TEST_CASE("face transition locates created faces in prism leaves") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = prism_state(std::vector<char>{2},
        std::vector<FaceIdentity>(), std::vector<Polygon>(), stateReport) ;
  REQUIRE(has_state(source)) ;
  const FaceKey createdKey(face_origin::cell_interior, 52,
        std::vector<int>{2, 0, 2, 2}, std::vector<int>{2, 0, 2, 3}) ;
  CPTR<FaceState> target = prism_state(std::vector<char>{2, 2},
        std::vector<FaceIdentity>{FaceIdentity(500, createdKey)},
        std::vector<Polygon>{interior_square()}, stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.remap.createdFaces == 1) ;
  CHECK(report.remap.removedFaces == 0) ;

  CellId sourceCell = 0 ;
  REQUIRE(transition->isCreatedFace(persistentFaceId(createdKey), sourceCell)) ;
  CHECK(sourceCell == persistentCellId(52, {2, 0})) ;
}

/// General-cell transitions use the same normalized leaf paths and can locate
/// a newly created face without interpreting a cell-plan byte stream.
TEST_CASE("face transition locates created faces in general-cell leaves") {
  FaceTransitionReport stateReport ;
  const std::vector<RootCellState> sourceRoots{
        RootCellState(62, std::vector<std::vector<int>>{std::vector<int>()},
              cell_topology::general)} ;
  CPTR<FaceState> source = FaceState::create(std::vector<FaceIdentity>(),
        std::vector<Polygon>(), sourceRoots, stateReport) ;
  REQUIRE(has_state(source)) ;

  const FaceKey createdKey(face_origin::cell_interior, 62,
        std::vector<int>{1, 0}, std::vector<int>{1, 1}) ;
  const std::vector<RootCellState> targetRoots{RootCellState(
        62, {{1, 0}, {1, 1}, {1, 2}, {1, 3}}, cell_topology::general)} ;
  CPTR<FaceState> target = FaceState::create(
        std::vector<FaceIdentity>{FaceIdentity(600, createdKey)},
        std::vector<Polygon>{interior_square()}, targetRoots, stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.remap.createdFaces == 1) ;
  CellId sourceCell = 0 ;
  REQUIRE(transition->isCreatedFace(persistentFaceId(createdKey), sourceCell)) ;
  CHECK(sourceCell == persistentCellId(62, {})) ;
}

/// Refined prism faces may retain a point on a straight edge; roundoff at that
/// point must not make an otherwise convex face appear geometrically invalid.
TEST_CASE("face state accepts collinear vertices from prism refinement") {
  const std::vector<FaceIdentity> identities{
        FaceIdentity(
              500, FaceKey(face_origin::cell_interior, 52, {3, 1}, {3, 2})),
        FaceIdentity(
              501, FaceKey(face_origin::cell_interior, 52, {3, 1}, {3, 4}))} ;
  const std::vector<Polygon> polygons{
        Polygon{
              vector3d<double>(0.35355339059327373, 0.35355339059327373, 0.25),
              vector3d<double>(0.35355339059327379, 0.35355339059327379, 0.5),
              vector3d<double>(0.5, 0.0, 0.5), vector3d<double>(0.5, 0.0, 0.0),
              vector3d<double>(0.35355339059327373, 0.35355339059327373, 0.0)},
        Polygon{vector3d<double>(0.0, 0.0, 0.5),
              vector3d<double>(0.5, 0.0, 0.5),
              vector3d<double>(0.35355339059327379, 0.35355339059327379, 0.5),
              vector3d<double>(0.17677669529663689, 0.42677669529663687, 0.5),
              vector3d<double>(0.0, 0.5, 0.5)}} ;

  FaceTransitionReport report ;
  CPTR<FaceState> state =
        prism_state(std::vector<char>{3}, identities, polygons, report) ;
  REQUIRE(has_state(state)) ;
  CHECK(report.valid) ;
  CHECK(report.invalidPolygons == 0) ;
}

/// Face ancestry follows refinement paths even when unrelated faces happen to
/// occupy the same geometric region.
TEST_CASE("face transition narrows matching by leaf path") {
  const std::vector<FaceIdentity> sourceIdentities{
        FaceIdentity(10, FaceKey(face_origin::base_face, 10, {3, 0})),
        FaceIdentity(11, FaceKey(face_origin::base_face, 10, {3, 1}))} ;
  const std::vector<Polygon> sourcePolygons{
        square(0.0, 1.0, 0.0, 1.0), square(0.0, 1.0, 0.0, 1.0)} ;
  const std::vector<FaceIdentity> targetIdentities{
        FaceIdentity(20, FaceKey(face_origin::base_face, 10, {3, 0, 3, 2})),
        FaceIdentity(21, FaceKey(face_origin::base_face, 10, {3, 0, 3, 3})),
        FaceIdentity(22, FaceKey(face_origin::base_face, 10, {3, 1}))} ;
  const std::vector<Polygon> targetPolygons{square(0.0, 0.5, 0.0, 1.0),
        square(0.5, 1.0, 0.0, 1.0), square(0.0, 1.0, 0.0, 1.0)} ;

  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = FaceState::create(sourceIdentities, sourcePolygons,
        std::vector<RootCellState>(), stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = FaceState::create(targetIdentities, targetPolygons,
        std::vector<RootCellState>(), stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CPTR<FaceRemap> transition = buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
  REQUIRE(has_transition(transition)) ;
  REQUIRE(report.valid) ;
  CHECK(report.remap.contributions == 3) ;

  std::vector<FaceId> targets ;
  transition->targetFaces(face_id(face_origin::base_face, 10, {3, 0}), targets) ;
  std::vector<FaceId> expected{
        face_id(face_origin::base_face, 10, {3, 0, 3, 2}),
        face_id(face_origin::base_face, 10, {3, 0, 3, 3})} ;
  std::sort(expected.begin(), expected.end()) ;
  CHECK(targets == expected) ;
  transition->targetFaces(face_id(face_origin::base_face, 10, {3, 1}), targets) ;
  CHECK(targets ==
        std::vector<FaceId>{face_id(face_origin::base_face, 10, {3, 1})}) ;
}

/// Geometric validation and overlap measurement must use the face's own scale
/// and must not change when the same mesh is translated.
TEST_CASE("face transition is scale and translation independent") {
  const double scales[2] = {1.0e-12, 1.0} ;
  const double offsets[2] = {0.0, 1.0e12} ;
  for (int geometryCase = 0; geometryCase < 2; ++geometryCase) {
    const double scale = scales[geometryCase] ;
    const double offset = offsets[geometryCase] ;
    const std::vector<FaceIdentity> sourceIdentities{
          FaceIdentity(10, FaceKey(face_origin::base_face, 10, {}))} ;
    const std::vector<FaceIdentity> targetIdentities{
          FaceIdentity(20, FaceKey(face_origin::base_face, 10, {1, 0})),
          FaceIdentity(21, FaceKey(face_origin::base_face, 10, {1, 1}))} ;
    const std::vector<Polygon> sourcePolygons{
          square(offset, offset + scale, offset, offset + scale)} ;
    const std::vector<Polygon> targetPolygons{
          square(offset, offset + 0.5 * scale, offset, offset + scale),
          square(offset + 0.5 * scale, offset + scale, offset, offset + scale)} ;

    FaceTransitionReport stateReport ;
    CPTR<FaceState> source = FaceState::create(sourceIdentities, sourcePolygons,
          std::vector<RootCellState>(), stateReport) ;
    REQUIRE(has_state(source)) ;
    CPTR<FaceState> target = FaceState::create(targetIdentities, targetPolygons,
          std::vector<RootCellState>(), stateReport) ;
    REQUIRE(has_state(target)) ;

    FaceTransitionReport report ;
    CPTR<FaceRemap> transition = buildFaceRemap(
          const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report) ;
    REQUIRE(has_transition(transition)) ;
    REQUIRE(report.valid) ;
    CHECK(report.remap.contributions == 2) ;
    for (size_t contribution = 0; contribution < transition->overlaps().size();
          ++contribution)
      CHECK(transition->overlaps()[contribution].area ==
            doctest::Approx(0.5 * scale * scale)) ;
  }
}

/// Coverage validation is relative to face area, including for meshes whose
/// characteristic length is much smaller than one.
TEST_CASE("face remap detects relative area loss on a small face") {
  const double area = 1.0e-24 ;
  const vector3d<double> centroid(0.5e-12, 0.5e-12, 0.0) ;
  const std::vector<FaceGeometry> source{FaceGeometry(10, area, centroid)} ;
  const std::vector<FaceGeometry> target{FaceGeometry(20, area, centroid)} ;
  const std::vector<FaceOverlap> contributions{
        FaceOverlap(10, 20, 0.5 * area, centroid, 1)} ;

  FaceRemapReport report ;
  CPTR<FaceRemap> remap = FaceRemap::create(source, target, contributions,
        std::vector<CreatedFace>(), std::vector<RemovedFace>(), report) ;
  CHECK(remap == static_cast<FaceRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.missingSourceFaces == 1) ;
  CHECK(report.missingTargetFaces == 1) ;
}

/// A root cannot change cell topology between two accepted mesh states even
/// when the compact split-code values happen to be legal for both trees.
TEST_CASE("face transition rejects a root topology change") {
  FaceTransitionReport stateReport ;
  const std::vector<FaceIdentity> identities{
        FaceIdentity(10, FaceKey(face_origin::base_face, 10, {}))} ;
  const std::vector<Polygon> polygons{square(0.0, 1.0, 0.0, 1.0)} ;
  CPTR<FaceState> source = FaceState::create(identities, polygons,
        {root_state_from_plan(42, std::vector<char>{1})}, stateReport) ;
  CPTR<FaceState> target = FaceState::create(identities, polygons,
        {root_state_from_plan(42, std::vector<char>{1}, cell_topology::prism)},
        stateReport) ;
  REQUIRE(has_state(source)) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CHECK_FALSE(has_transition(buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report))) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == face_transition_status::invalid_identity) ;
  CHECK(report.invalidIdentities == 1) ;
}

/// Leaf-path validation must interpret split codes according to the declared
/// root topology rather than accepting the hexahedral code set everywhere.
TEST_CASE("face state validates leaf paths against root topology") {
  FaceTransitionReport report ;
  CPTR<FaceState> hexState =
        FaceState::create(std::vector<FaceIdentity>(), std::vector<Polygon>(),
              {root_state_from_plan(42, std::vector<char>{7})}, report) ;
  REQUIRE(has_state(hexState)) ;
  REQUIRE(report.valid) ;

  CPTR<FaceState> prismState = prism_state(std::vector<char>{3},
        std::vector<FaceIdentity>(), std::vector<Polygon>(), report) ;
  REQUIRE(has_state(prismState)) ;
  REQUIRE(report.valid) ;

  CPTR<FaceState> invalidPrism = prism_state(std::vector<char>{7},
        std::vector<FaceIdentity>(), std::vector<Polygon>(), report) ;
  CHECK_FALSE(has_state(invalidPrism)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == face_transition_status::invalid_identity) ;
}

/// A face state accepts only canonical keys whose cell-interior paths name two
/// distinct accepted leaves of the identified root cell.
TEST_CASE("face state validates canonical face keys") {
  const std::vector<RootCellState> roots{
        root_state_from_plan(42, std::vector<char>{1})} ;
  const std::vector<FaceKey> invalidKeys{
        FaceKey(face_origin::base_face, -1, {}),
        FaceKey(face_origin::base_face, 10, {1, 0}, {1, 1}),
        FaceKey(face_origin::cell_interior, 42, {1, 1}, {1, 0}),
        FaceKey(face_origin::cell_interior, 42, {1, 0}, {1, 0}),
        FaceKey(face_origin::cell_interior, 43, {1, 0}, {1, 1}),
        FaceKey(face_origin::cell_interior, 42, {}, {1, 1}),
        FaceKey(static_cast<face_origin::value>(99), 10, {})} ;

  for (size_t key = 0; key < invalidKeys.size(); ++key) {
    FaceTransitionReport report ;
    CPTR<FaceState> state =
          FaceState::create({FaceIdentity(10, invalidKeys[key])},
                {interior_square()}, roots, report) ;
    CHECK_FALSE(has_state(state)) ;
    CHECK_FALSE(report.valid) ;
    CHECK(report.status == face_transition_status::invalid_identity) ;
    CHECK(report.invalidIdentities == 1) ;
  }
}

/// Replacing a retained split direction gives child paths different geometric
/// meanings, so the transition must fail closed instead of guessing ancestry.
TEST_CASE("face transition rejects incompatible incremental plan changes") {
  FaceTransitionReport stateReport ;
  const std::vector<FaceIdentity> identity{
        FaceIdentity(10, FaceKey(face_origin::base_face, 10, {1, 0}))} ;
  const std::vector<Polygon> polygon{square(0.0, 1.0, 0.0, 1.0)} ;
  CPTR<FaceState> source = FaceState::create(identity, polygon,
        {root_state_from_plan(42, std::vector<char>{1, 0})}, stateReport) ;
  CPTR<FaceState> target = FaceState::create(identity, polygon,
        {root_state_from_plan(42, std::vector<char>{2, 0})}, stateReport) ;
  REQUIRE(has_state(source)) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CHECK_FALSE(has_transition(buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report))) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == face_transition_status::unsupported_plan_change) ;
  CHECK(report.unsupportedRootPlans == 1) ;
}

/// A retained general-cell split must keep the same number of children because
/// otherwise an unchanged child ordinal no longer identifies the same region.
TEST_CASE("face transition rejects incompatible general-cell arity") {
  FaceTransitionReport stateReport ;
  CPTR<FaceState> source = FaceState::create({}, {},
        {RootCellState(
              62, {{1, 0}, {1, 1}, {1, 2}, {1, 3}}, cell_topology::general)},
        stateReport) ;
  REQUIRE(has_state(source)) ;
  CPTR<FaceState> target = FaceState::create({}, {},
        {RootCellState(62, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}},
              cell_topology::general)},
        stateReport) ;
  REQUIRE(has_state(target)) ;

  FaceTransitionReport report ;
  CHECK_FALSE(has_transition(buildFaceRemap(
        const_CPTR<FaceState>(source), const_CPTR<FaceState>(target), report))) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == face_transition_status::unsupported_plan_change) ;
  CHECK(report.unsupportedRootPlans == 1) ;
}

/// Persistent identities include the split code as well as the child ordinal,
/// because child zero has a different meaning under different split modes.
TEST_CASE("persistent identities distinguish split modes") {
  const FaceId firstFace = face_id(face_origin::base_face, 10, {1, 0}) ;
  const FaceId secondFace = face_id(face_origin::base_face, 10, {2, 0}) ;
  const CellId firstCell = persistentCellId(42, {1, 0}) ;
  const CellId secondCell = persistentCellId(42, {2, 0}) ;

  CHECK(firstFace != 0) ;
  CHECK(secondFace != 0) ;
  CHECK(firstFace != secondFace) ;
  CHECK(firstCell != 0) ;
  CHECK(secondCell != 0) ;
  CHECK(firstCell != secondCell) ;
  CHECK(face_id(face_origin::base_face, 10, {1}) == 0) ;
  CHECK(face_id(face_origin::base_face, 10, {0, 0}) == 0) ;
  CHECK(face_id(face_origin::base_face, 10, {1, -1}) == 0) ;
  CHECK(face_id(face_origin::base_face, -1, {}) == 0) ;
  CHECK(persistentFaceId(FaceKey(face_origin::base_face, 10, {}, {1, 0})) == 0) ;
  CHECK(persistentFaceId(
              FaceKey(face_origin::cell_interior, 42, {1, 1}, {1, 0})) == 0) ;
  CHECK(persistentFaceId(
              FaceKey(static_cast<face_origin::value>(99), 10, {})) == 0) ;
  CHECK(persistentCellId(42, {1}) == 0) ;
  CHECK(persistentCellId(-1, {}) == 0) ;
}

/// Scheduler facts encode each breadth-first leaf as a step count followed by
/// (split code, child ordinal) pairs, without retaining tree pointers.
TEST_CASE("accepted plan leaf paths round trip") {
  std::vector<int> encoded ;
  std::vector<std::vector<int>> paths ;

  REQUIRE(detail::encodeQuadFaceLeafPaths(std::vector<char>{3, 0}, encoded)) ;
  REQUIRE(detail::decodeLeafPaths(encoded, paths)) ;
  CHECK(paths == std::vector<std::vector<int>>{{3, 0}, {3, 1}, {3, 2}, {3, 3}}) ;
  std::vector<int> reencoded ;
  REQUIRE(detail::encodeLeafPaths(paths, reencoded)) ;
  CHECK(reencoded == encoded) ;

  REQUIRE(detail::encodeHexCellLeafPaths(std::vector<char>{7, 0}, encoded)) ;
  REQUIRE(detail::decodeLeafPaths(encoded, paths)) ;
  REQUIRE(paths.size() == 8) ;
  for (size_t child = 0; child < paths.size(); ++child)
    CHECK(paths[child] == std::vector<int>{7, int(child)}) ;

  CHECK_FALSE(detail::decodeLeafPaths(std::vector<int>{1, 1, 7}, paths)) ;
}

/// Prism paths must follow the varying child counts used by Prism::empty_split,
/// including the four-sided children produced by transverse refinement.
TEST_CASE("prism cell leaf paths match core replay") {
  const std::vector<std::vector<char>> plans{std::vector<char>(),
        std::vector<char>{1}, std::vector<char>{2}, std::vector<char>{3},
        std::vector<char>{1, 2}, std::vector<char>{2, 2},
        std::vector<char>{3, 3}} ;
  for (size_t plan = 0; plan < plans.size(); ++plan) {
    std::vector<int> encoded ;
    std::vector<std::vector<int>> encodedPaths ;
    std::vector<std::vector<int>> replayPaths ;
    REQUIRE(detail::encodePrismCellLeafPaths(plans[plan], encoded)) ;
    REQUIRE(detail::decodeLeafPaths(encoded, encodedPaths)) ;
    REQUIRE(prism_replay_paths(plans[plan], replayPaths)) ;
    CHECK(encodedPaths == replayPaths) ;
  }

  std::vector<int> encoded ;
  CHECK_FALSE(detail::encodePrismCellLeafPaths(std::vector<char>{4}, encoded)) ;
  CHECK_FALSE(
        detail::encodePrismCellLeafPaths(std::vector<char>{0, 2}, encoded)) ;
}

/// Prism end faces begin as triangles, while every child after the first
/// split is a quadrilateral; encoded paths must retain that change in arity.
TEST_CASE("general face leaf paths match core replay") {
  const std::vector<std::vector<char>> plans{
        std::vector<char>(), std::vector<char>{1}, std::vector<char>{1, 1}} ;
  for (size_t plan = 0; plan < plans.size(); ++plan) {
    std::vector<int> encoded ;
    std::vector<std::vector<int>> encodedPaths ;
    std::vector<std::vector<int>> replayPaths ;
    REQUIRE(detail::encodeGeneralFaceLeafPaths(plans[plan], 3, encoded)) ;
    REQUIRE(detail::decodeLeafPaths(encoded, encodedPaths)) ;
    REQUIRE(general_face_replay_paths(plans[plan], 3, replayPaths)) ;
    CHECK(encodedPaths == replayPaths) ;
  }

  std::vector<int> encoded ;
  CHECK_FALSE(
        detail::encodeGeneralFaceLeafPaths(std::vector<char>{1}, 2, encoded)) ;
  CHECK_FALSE(
        detail::encodeGeneralFaceLeafPaths(std::vector<char>{2}, 3, encoded)) ;
  CHECK_FALSE(detail::encodeGeneralFaceLeafPaths(
        std::vector<char>{0, 1}, 3, encoded)) ;
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv) ;
  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;
  Loci::Finalize() ;
  return result ;
}

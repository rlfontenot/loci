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
#include <FVMAdapt2/node_edge.h>
#include <FVMAdapt2/node_transition.h>
#include "mesh_state.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <map>
#include <vector>

using namespace Loci ;

namespace {

  NodeConstruction base_node(long long fileNumber, double x) {
    return NodeConstruction::baseNode(
          fileNumber, vector3d<double>(x, 0.0, 0.0)) ;
  }

  NodeConstruction midpoint(
        const NodeConstruction& first, const NodeConstruction& second) {
    return NodeConstruction::constructed(node_construction::edge,
          0.5 * (first.position + second.position),
          std::vector<NodeParent>{
                NodeParent(first.node, 0.5), NodeParent(second.node, 0.5)}) ;
  }

  FineNodeConstruction transport_record(const NodeConstruction& construction,
        const std::map<NodeId, int>& currentNumbers) {
    FineNodeConstruction result ;
    result.node = construction.node ;
    result.kind = construction.kind ;
    result.baseFileNumber = construction.baseFileNumber ;
    result.parentCount = int(construction.parents.size()) ;
    for (size_t parent = 0; parent < construction.parents.size(); ++parent) {
      result.parentIds[parent] = construction.parents[parent].node ;
      const std::map<NodeId, int>::const_iterator number =
            currentNumbers.find(construction.parents[parent].node) ;
      result.parentNodeNumbers[parent] =
            number == currentNumbers.end() ? -1 : number->second ;
      result.parentWeights[parent] = construction.parents[parent].weight ;
    }
    return result ;
  }

  CPTR<NodeRemap> build_serial_remap(
        const std::vector<NodeConstruction>& previous,
        const std::vector<NodeConstruction>& current,
        NodeTransitionReport& report, double relativeTolerance = 1.0e-10) {
    REQUIRE(MPI_processes == 1) ;

    std::vector<entitySet> previousPartition(1) ;
    std::vector<entitySet> currentPartition(1) ;
    if (!previous.empty())
      previousPartition[0] = interval(0, int(previous.size()) - 1) ;
    if (!current.empty())
      currentPartition[0] = interval(0, int(current.size()) - 1) ;

    store<NodeId> previousIds ;
    store<vector3d<double>> previousPositions ;
    previousIds.allocate(previousPartition[0]) ;
    previousPositions.allocate(previousPartition[0]) ;
    for (size_t node = 0; node < previous.size(); ++node) {
      previousIds[int(node)] = previous[node].node ;
      previousPositions[int(node)] = previous[node].position ;
    }

    std::map<NodeId, int> currentNumbers ;
    for (size_t node = 0; node < current.size(); ++node)
      currentNumbers[current[node].node] = int(node) ;
    store<NodeId> currentIds ;
    store<FineNodeConstruction> currentConstructions ;
    store<vector3d<double>> currentPositions ;
    currentIds.allocate(currentPartition[0]) ;
    currentConstructions.allocate(currentPartition[0]) ;
    currentPositions.allocate(currentPartition[0]) ;
    for (size_t node = 0; node < current.size(); ++node) {
      currentIds[int(node)] = current[node].node ;
      currentConstructions[int(node)] =
            transport_record(current[node], currentNumbers) ;
      currentPositions[int(node)] = current[node].position ;
    }

    return detail::buildNodeRemap(previousIds, previousPositions,
          previousPartition, currentIds, currentConstructions, currentPositions,
          currentPartition, report, relativeTolerance) ;
  }

  std::map<NodeId, double> sources(const NodeRemap& remap, NodeId target) {
    size_t begin = 0 ;
    size_t end = 0 ;
    std::map<NodeId, double> result ;
    if (!remap.nodeContributions(target, begin, end))
      return result ;
    for (size_t entry = begin; entry < end; ++entry) {
      const NodeContribution& contribution = remap.nodeContributions()[entry] ;
      result[contribution.sourceNode] = contribution.weight ;
    }
    return result ;
  }

  void check_position_reproduction(const NodeRemap& remap) {
    std::vector<vector3d<double>> sourcePositions ;
    for (size_t node = 0; node < remap.sourceNodeGeometry().size(); ++node)
      sourcePositions.push_back(remap.sourceNodeGeometry()[node].position) ;
    std::vector<vector3d<double>> targetPositions ;
    REQUIRE(remap.interpolateNodeData(sourcePositions, targetPositions)) ;
    REQUIRE(targetPositions.size() == remap.targetNodeGeometry().size()) ;
    for (size_t node = 0; node < targetPositions.size(); ++node) {
      const vector3d<double>& expected =
            remap.targetNodeGeometry()[node].position ;
      CHECK(targetPositions[node].x == doctest::Approx(expected.x)) ;
      CHECK(targetPositions[node].y == doctest::Approx(expected.y)) ;
      CHECK(targetPositions[node].z == doctest::Approx(expected.z)) ;
    }
  }

  FineNodeConstruction transport_record(const NodeConstruction& construction,
        const std::vector<int>& parentNumbers) {
    FineNodeConstruction result ;
    result.node = construction.node ;
    result.kind = construction.kind ;
    result.baseFileNumber = construction.baseFileNumber ;
    result.parentCount = int(construction.parents.size()) ;
    REQUIRE(parentNumbers.size() == construction.parents.size()) ;
    for (size_t parent = 0; parent < construction.parents.size(); ++parent) {
      result.parentIds[parent] = construction.parents[parent].node ;
      result.parentNodeNumbers[parent] = parentNumbers[parent] ;
      result.parentWeights[parent] = construction.parents[parent].weight ;
    }
    return result ;
  }
}

/// Refinement helpers that only carry geometry still construct the correct
/// point, but do not claim ancestry that their parent nodes cannot prove.
TEST_CASE("geometry-only nodes preserve refinement geometry") {
  Node left(vector3d<double>(0.0, 2.0, 4.0)) ;
  Node right(vector3d<double>(2.0, 4.0, 6.0)) ;
  Node* middle = Node::constructed(node_construction::edge,
        std::vector<Node*>{&left, &right}, std::vector<double>{0.5, 0.5}) ;

  REQUIRE(middle != 0) ;
  CHECK(middle->p.x == doctest::Approx(1.0)) ;
  CHECK(middle->p.y == doctest::Approx(3.0)) ;
  CHECK(middle->p.z == doctest::Approx(5.0)) ;
  FineNodeConstruction construction ;
  CHECK_FALSE(middle->fineConstruction(construction)) ;
  delete middle ;
}

/// Nodes already present in the previous mesh retain their identity and
/// contribute to themselves with unit weight.
TEST_CASE("node remap retains existing nodes by identity") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(10, 0.0) ;
  const NodeConstruction right = base_node(11, 2.0) ;
  NodeTransitionReport report ;
  CPTR<NodeRemap> remap =
        build_serial_remap(std::vector<NodeConstruction>{left, right},
              std::vector<NodeConstruction>{right, left}, report) ;

  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  REQUIRE(report.valid) ;
  CHECK(report.retainedNodes == 2) ;
  CHECK(report.createdNodes == 0) ;
  CHECK(report.contributions == 2) ;
  CHECK(sources(*remap, left.node) ==
        std::map<NodeId, double>{{left.node, 1.0}}) ;
  CHECK(sources(*remap, right.node) ==
        std::map<NodeId, double>{{right.node, 1.0}}) ;
  check_position_reproduction(*remap) ;
}

/// A first edge split reconstructs its midpoint from the two previous
/// endpoints and exposes the same weights used to construct its geometry.
TEST_CASE("node remap exposes midpoint ancestry") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(20, 0.0) ;
  const NodeConstruction right = base_node(21, 2.0) ;
  const NodeConstruction middle = midpoint(left, right) ;
  CHECK(middle.node == midpoint(right, left).node) ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap =
        build_serial_remap(std::vector<NodeConstruction>{left, right},
              std::vector<NodeConstruction>{left, middle, right}, report) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  REQUIRE(report.valid) ;
  CHECK(report.retainedNodes == 2) ;
  CHECK(report.createdNodes == 1) ;
  CHECK(sources(*remap, middle.node) ==
        std::map<NodeId, double>{{left.node, 0.5}, {right.node, 0.5}}) ;
  check_position_reproduction(*remap) ;
}

/// Recursion expands nodes created in the current cycle, but stops at a fine
/// node already present in the immediately previous accepted mesh.
TEST_CASE("node remap stops ancestry at previous mesh nodes") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(30, 0.0) ;
  const NodeConstruction right = base_node(31, 2.0) ;
  const NodeConstruction middle = midpoint(left, right) ;
  const NodeConstruction quarter = midpoint(left, middle) ;
  const std::vector<NodeConstruction> deep{left, quarter, middle, right} ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> direct = build_serial_remap(
        std::vector<NodeConstruction>{left, right}, deep, report) ;
  REQUIRE(direct != static_cast<NodeRemap*>(0)) ;
  CHECK(sources(*direct, quarter.node) ==
        std::map<NodeId, double>{{left.node, 0.75}, {right.node, 0.25}}) ;

  CPTR<NodeRemap> secondCycle = build_serial_remap(
        std::vector<NodeConstruction>{left, middle, right}, deep, report) ;
  REQUIRE(secondCycle != static_cast<NodeRemap*>(0)) ;
  const std::map<NodeId, double> secondCycleSources =
        sources(*secondCycle, quarter.node) ;
  CHECK(secondCycleSources ==
        std::map<NodeId, double>{{left.node, 0.5}, {middle.node, 0.5}}) ;
  CHECK(secondCycleSources.count(right.node) == 0) ;
  check_position_reproduction(*secondCycle) ;
}

/// Derefinement drops fine nodes from the target mesh; recreating the same
/// topology later restores the same persistent identity and fresh ancestry.
TEST_CASE("node identity survives derefine and recreate") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(40, 0.0) ;
  const NodeConstruction right = base_node(41, 2.0) ;
  const NodeConstruction middle = midpoint(left, right) ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> derefine =
        build_serial_remap(std::vector<NodeConstruction>{left, middle, right},
              std::vector<NodeConstruction>{left, right}, report) ;
  REQUIRE(derefine != static_cast<NodeRemap*>(0)) ;
  CHECK(report.retainedNodes == 2) ;
  CHECK(report.createdNodes == 0) ;
  size_t begin = 0 ;
  size_t end = 0 ;
  CHECK_FALSE(derefine->nodeContributions(middle.node, begin, end)) ;

  const NodeConstruction recreated = midpoint(left, right) ;
  CHECK(recreated.node == middle.node) ;
  CPTR<NodeRemap> refine =
        build_serial_remap(std::vector<NodeConstruction>{left, right},
              std::vector<NodeConstruction>{left, recreated, right}, report) ;
  REQUIRE(refine != static_cast<NodeRemap*>(0)) ;
  CHECK(report.createdNodes == 1) ;
  CHECK(sources(*refine, recreated.node) ==
        std::map<NodeId, double>{{left.node, 0.5}, {right.node, 0.5}}) ;
}

/// A current base node missing from the previous mesh cannot be assigned a
/// source value, so the node remap fails closed.
TEST_CASE("node remap reports a missing previous base node") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(60, 0.0) ;
  const NodeConstruction right = base_node(61, 2.0) ;
  NodeTransitionReport report ;
  CPTR<NodeRemap> remap =
        build_serial_remap(std::vector<NodeConstruction>{left},
              std::vector<NodeConstruction>{left, right}, report) ;
  CHECK(remap == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::missing_source_node) ;
  CHECK(report.missingSourceNodes == 1) ;
}

/// Malformed construction weights are rejected before they can become an
/// incomplete solver-facing remap.
TEST_CASE("node remap rejects malformed construction weights") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(65, 0.0) ;
  const NodeConstruction right = base_node(66, 2.0) ;
  NodeConstruction middle = midpoint(left, right) ;
  middle.parents[0].weight = 0.0 ;
  middle.parents[1].weight = 1.0 ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap =
        build_serial_remap(std::vector<NodeConstruction>{left, right},
              std::vector<NodeConstruction>{left, middle, right}, report) ;
  CHECK(remap == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::missing_state) ;
}

/// The node builder rejects construction weights that do not reproduce the
/// target node's geometry.
TEST_CASE("node remap validates geometric reproduction") {
  if (MPI_processes != 1)
    return ;
  const NodeConstruction left = base_node(70, 0.0) ;
  const NodeConstruction right = base_node(71, 2.0) ;
  NodeConstruction middle = midpoint(left, right) ;
  middle.position.x = 1.25 ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap =
        build_serial_remap(std::vector<NodeConstruction>{left, right},
              std::vector<NodeConstruction>{left, middle, right}, report) ;
  CHECK(remap == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::inconsistent_positions) ;
  CHECK(report.inconsistentPositions == 1) ;
  CHECK(report.maximumPositionError == doctest::Approx(0.25)) ;
}

/// Empty accepted partitions produce an empty remap, which is required when
/// a mesh has fewer locally owned nodes than MPI ranks.
TEST_CASE("node remap accepts an empty mesh") {
  if (MPI_processes != 1)
    return ;
  NodeTransitionReport report ;
  CPTR<NodeRemap> remap = build_serial_remap(std::vector<NodeConstruction>(),
        std::vector<NodeConstruction>(), report) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  CHECK(report.valid) ;
  CHECK(report.sourceNodes == 0) ;
  CHECK(report.targetNodes == 0) ;
  CHECK(remap->nodeContributions().empty()) ;
}

/// Node IDs are 64-bit values. Distinct IDs with matching low 32 bits must
/// remain distinct while node origins are validated.
TEST_CASE("node remap preserves full width node identities") {
  const NodeId first = 17 ;
  const NodeId second = first + (NodeId(1) << 32) ;
  const std::vector<NodeGeometry> geometry{
        NodeGeometry(first, vector3d<double>(0.0, 0.0, 0.0)),
        NodeGeometry(second, vector3d<double>(1.0, 0.0, 0.0))} ;
  const std::vector<NodeContribution> contributions{
        NodeContribution(first, first, 1.0),
        NodeContribution(second, second, 1.0)} ;
  const std::vector<NodeOrigin> origins{
        NodeOrigin(first, node_origin::base_node),
        NodeOrigin(second, node_origin::base_node)} ;

  NodeRemapReport report ;
  CPTR<NodeRemap> remap =
        NodeRemap::create(geometry, geometry, contributions, origins, report) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  CHECK(report.valid) ;
  CHECK(report.invalidOrigins == 0) ;
}

/// The node assembler resolves a midpoint parent owned by another MPI rank
/// and accepts ranks that own no nodes.
TEST_CASE("node remap resolves remote parents and empty ranks") {
  const NodeConstruction left = base_node(80, 0.0) ;
  const NodeConstruction right = base_node(81, 2.0) ;
  const NodeConstruction middle = midpoint(left, right) ;

  std::vector<entitySet> previousPartition(MPI_processes) ;
  std::vector<entitySet> currentPartition(MPI_processes) ;
  if (MPI_processes == 1) {
    previousPartition[0] = interval(0, 1) ;
    currentPartition[0] = interval(0, 2) ;
  } else {
    previousPartition[0] += 0 ;
    previousPartition[1] += 1 ;
    currentPartition[0] += 0 ;
    currentPartition[0] += 2 ;
    currentPartition[1] += 1 ;
  }

  store<NodeId> previousIds ;
  store<vector3d<double>> previousPositions ;
  previousIds.allocate(previousPartition[MPI_rank]) ;
  previousPositions.allocate(previousPartition[MPI_rank]) ;
  if (previousPartition[MPI_rank].inSet(0)) {
    previousIds[0] = left.node ;
    previousPositions[0] = left.position ;
  }
  if (previousPartition[MPI_rank].inSet(1)) {
    previousIds[1] = right.node ;
    previousPositions[1] = right.position ;
  }

  store<NodeId> currentIds ;
  store<FineNodeConstruction> constructions ;
  store<vector3d<double>> currentPositions ;
  currentIds.allocate(currentPartition[MPI_rank]) ;
  constructions.allocate(currentPartition[MPI_rank]) ;
  currentPositions.allocate(currentPartition[MPI_rank]) ;
  if (currentPartition[MPI_rank].inSet(0)) {
    currentIds[0] = left.node ;
    currentPositions[0] = left.position ;
    constructions[0] = transport_record(left, std::vector<int>()) ;
  }
  if (currentPartition[MPI_rank].inSet(1)) {
    currentIds[1] = right.node ;
    currentPositions[1] = right.position ;
    constructions[1] = transport_record(right, std::vector<int>()) ;
  }
  if (currentPartition[MPI_rank].inSet(2)) {
    currentIds[2] = middle.node ;
    currentPositions[2] = middle.position ;
    constructions[2] = transport_record(middle, std::vector<int>{0, 1}) ;
  }

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap = detail::buildNodeRemap(previousIds, previousPositions,
        previousPartition, currentIds, constructions, currentPositions,
        currentPartition, report) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  CHECK(report.valid) ;
  CHECK(remap->targetNodeGeometry().size() == currentIds.domain().size()) ;
  if (currentIds.domain().inSet(2)) {
    CHECK(sources(*remap, middle.node) ==
          std::map<NodeId, double>{{left.node, 0.5}, {right.node, 0.5}}) ;
  }
}

/// A rank with incomplete node inputs makes the collective build fail
/// cleanly instead of leaving the other ranks inside an MPI exchange.
TEST_CASE("node remap rejects rank-local missing state collectively") {
  std::vector<entitySet> previousPartition(MPI_processes) ;
  std::vector<entitySet> currentPartition(MPI_processes) ;
  previousPartition[0] += 0 ;
  currentPartition[0] += 0 ;

  store<NodeId> previousIds ;
  store<vector3d<double>> previousPositions ;
  store<NodeId> currentIds ;
  store<FineNodeConstruction> constructions ;
  store<vector3d<double>> currentPositions ;
  const entitySet localNodes =
        MPI_rank == 0 ? EMPTY : previousPartition[MPI_rank] ;
  previousIds.allocate(localNodes) ;
  previousPositions.allocate(localNodes) ;
  currentIds.allocate(currentPartition[MPI_rank]) ;
  constructions.allocate(currentPartition[MPI_rank]) ;
  currentPositions.allocate(currentPartition[MPI_rank]) ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap = detail::buildNodeRemap(previousIds, previousPositions,
        previousPartition, currentIds, constructions, currentPositions,
        currentPartition, report) ;
  CHECK(remap == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::missing_state) ;
}

/// A rank with incomplete redistribution inputs makes every rank stop before
/// the first MPI exchange rather than returning on that rank alone.
TEST_CASE("node remap rejects rank-local redistribution state collectively") {
  NodeRemapReport remapReport ;
  CPTR<NodeRemap> remap = NodeRemap::create(std::vector<NodeGeometry>(),
        std::vector<NodeGeometry>(), std::vector<NodeContribution>(),
        std::vector<NodeOrigin>(), remapReport) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;

  std::vector<entitySet> nodePartition(MPI_processes) ;
  if (MPI_rank == 0)
    nodePartition.pop_back() ;
  store<NodeId> generatedNodeIds ;
  generatedNodeIds.allocate(EMPTY) ;

  NodeTransitionReport report ;
  CPTR<NodeRemap> installed = detail::redistributeNodeRemap(
        remap, nodePartition, generatedNodeIds, report) ;
  CHECK(installed == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::missing_state) ;
}

/// An unchanged mesh has one unit contribution per node, regardless of the
/// order in which the caller supplies its records.
TEST_CASE("node remap retains an unordered mesh by identity") {
  const int count = 10000 ;
  std::vector<NodeGeometry> geometry ;
  std::vector<NodeContribution> contributions ;
  std::vector<NodeOrigin> origins ;
  for (int node = 1; node <= count; ++node) {
    geometry.push_back(NodeGeometry(node, vector3d<double>(node, 0, 0))) ;
    const NodeId reverse = count + 1 - node ;
    contributions.push_back(NodeContribution(reverse, reverse, 1.0)) ;
    origins.push_back(NodeOrigin(reverse, node_origin::base_node)) ;
  }

  NodeRemapReport report ;
  CPTR<NodeRemap> remap =
        NodeRemap::create(geometry, geometry, contributions, origins, report) ;
  REQUIRE(remap != static_cast<NodeRemap*>(0)) ;
  CHECK(report.valid) ;
  for (const NodeGeometry& node : geometry) {
    size_t begin = 0, end = 0 ;
    REQUIRE(remap->nodeContributions(node.node, begin, end)) ;
    REQUIRE(end - begin == 1) ;
    CHECK(remap->nodeContributions()[begin].sourceNode == node.node) ;
    CHECK(remap->nodeContributions()[begin].weight == 1.0) ;
  }
}

/// A malformed parent count on one rank must be rejected before indexing the
/// parent arrays, and ranks with no nodes must leave the build with it.
TEST_CASE("node remap rejects invalid parent counts collectively") {
  int parentCount = -1 ;
  SUBCASE("negative parent count") {}
  SUBCASE("too many parents") {
    parentCount = FineNodeConstruction::maximumParents + 1 ;
  }

  std::vector<entitySet> partition(MPI_processes) ;
  partition[0] += 0 ;
  store<NodeId> ids ;
  store<vector3d<double>> positions ;
  store<FineNodeConstruction> constructions ;
  ids.allocate(partition[MPI_rank]) ;
  positions.allocate(partition[MPI_rank]) ;
  constructions.allocate(partition[MPI_rank]) ;
  if (MPI_rank == 0) {
    const NodeConstruction node = base_node(90, 0.0) ;
    ids[0] = node.node ;
    positions[0] = node.position ;
    constructions[0] = transport_record(node, std::vector<int>()) ;
    constructions[0].parentCount = parentCount ;
  }

  NodeTransitionReport report ;
  CPTR<NodeRemap> remap = detail::buildNodeRemap(ids, positions, partition, ids,
        constructions, positions, partition, report) ;
  CHECK(remap == static_cast<NodeRemap*>(0)) ;
  CHECK_FALSE(report.valid) ;
  CHECK(report.status == node_transition_status::missing_state) ;
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv) ;
  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;
  Loci::Finalize() ;
  return result ;
}

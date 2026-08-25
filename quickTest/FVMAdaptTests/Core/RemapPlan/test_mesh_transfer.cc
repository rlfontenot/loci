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
#include <FVMAdapt/mesh_transfer.h>

#include <doctest.h>

#include <numeric>
#include <vector>

using namespace Loci ;

namespace {

  bool has_face_remap(const CPTR<AMRFaceRemap>& remap) {
    return remap != static_cast<AMRFaceRemap*>(0) ;
  }

  bool has_node_remap(const CPTR<AMRNodeRemap>& remap) {
    return remap != static_cast<AMRNodeRemap*>(0) ;
  }

  std::vector<AMRFaceGeometry> source_faces() {
    return std::vector<AMRFaceGeometry>{
      AMRFaceGeometry(10,2.0,vector3d<double>(0.0,0.0,0.0)),
      AMRFaceGeometry(20,1.0,vector3d<double>(2.5,0.0,0.0)),
      AMRFaceGeometry(21,1.0,vector3d<double>(3.5,0.0,0.0)),
      AMRFaceGeometry(30,1.0,vector3d<double>(5.0,0.0,0.0)),
      AMRFaceGeometry(40,1.0,vector3d<double>(7.0,0.0,0.0))
    } ;
  }

  std::vector<AMRFaceGeometry> target_faces() {
    return std::vector<AMRFaceGeometry>{
      AMRFaceGeometry(100,1.0,vector3d<double>(-0.5,0.0,0.0)),
      AMRFaceGeometry(101,1.0,vector3d<double>(0.5,0.0,0.0)),
      AMRFaceGeometry(102,2.0,vector3d<double>(3.0,0.0,0.0)),
      AMRFaceGeometry(103,1.0,vector3d<double>(5.0,0.0,0.0)),
      AMRFaceGeometry(104,1.0,vector3d<double>(8.0,0.0,0.0))
    } ;
  }

  std::vector<AMRFaceContribution> face_contributions() {
    return std::vector<AMRFaceContribution>{
      AMRFaceContribution(10,100,1.0,
                          vector3d<double>(-0.5,0.0,0.0),1),
      AMRFaceContribution(10,101,1.0,
                          vector3d<double>(0.5,0.0,0.0),1),
      AMRFaceContribution(20,102,1.0,
                          vector3d<double>(2.5,0.0,0.0),1),
      AMRFaceContribution(21,102,1.0,
                          vector3d<double>(3.5,0.0,0.0),-1),
      AMRFaceContribution(30,103,1.0,
                          vector3d<double>(5.0,0.0,0.0),-1)
    } ;
  }

  CPTR<AMRFaceRemap> make_face_remap(AMRFaceRemapReport& report) {
    return AMRFaceRemap::create(
      source_faces(),target_faces(),face_contributions(),
      std::vector<AMRCreatedFace>{AMRCreatedFace(104,60)},
      std::vector<AMRRemovedFace>{AMRRemovedFace(40,70)},report) ;
  }

  std::vector<AMRNodeGeometry> source_nodes() {
    return std::vector<AMRNodeGeometry>{
      AMRNodeGeometry(0,vector3d<double>(0.0,0.0,0.0)),
      AMRNodeGeometry(1,vector3d<double>(2.0,0.0,0.0)),
      AMRNodeGeometry(2,vector3d<double>(0.0,2.0,0.0)),
      AMRNodeGeometry(3,vector3d<double>(0.0,0.0,2.0))
    } ;
  }

  std::vector<AMRNodeGeometry> target_nodes() {
    return std::vector<AMRNodeGeometry>{
      AMRNodeGeometry(10,vector3d<double>(0.0,0.0,0.0)),
      AMRNodeGeometry(11,vector3d<double>(1.0,0.0,0.0)),
      AMRNodeGeometry(12,vector3d<double>(2.0/3.0,2.0/3.0,0.0)),
      AMRNodeGeometry(13,vector3d<double>(0.5,0.5,0.5))
    } ;
  }

  std::vector<AMRNodeContribution> node_contributions() {
    return std::vector<AMRNodeContribution>{
      AMRNodeContribution(0,10,1.0),
      AMRNodeContribution(0,11,0.5),
      AMRNodeContribution(1,11,0.5),
      AMRNodeContribution(0,12,1.0/3.0),
      AMRNodeContribution(1,12,1.0/3.0),
      AMRNodeContribution(2,12,1.0/3.0),
      AMRNodeContribution(0,13,0.25),
      AMRNodeContribution(1,13,0.25),
      AMRNodeContribution(2,13,0.25),
      AMRNodeContribution(3,13,0.25)
    } ;
  }

  std::vector<AMRNodeOrigin> node_origins() {
    return std::vector<AMRNodeOrigin>{
      AMRNodeOrigin(10,amr_node_origin::retained,0),
      AMRNodeOrigin(11,amr_node_origin::edge,100),
      AMRNodeOrigin(12,amr_node_origin::face,200),
      AMRNodeOrigin(13,amr_node_origin::cell,300)
    } ;
  }

} // namespace


/// Face ancestry should cover splitting, coarsening, reversed orientation,
/// and faces created or removed inside cells without guessing an old face ID.
TEST_CASE("face remap distinguishes inherited and cell-internal faces") {
  AMRFaceRemapReport report ;
  CPTR<AMRFaceRemap> remap = make_face_remap(report) ;
  REQUIRE(has_face_remap(remap)) ;
  REQUIRE(report.valid) ;
  CHECK(report.maximumSourceAreaError == doctest::Approx(0.0)) ;
  CHECK(report.maximumTargetAreaError == doctest::Approx(0.0)) ;

  int cell = -1 ;
  REQUIRE(remap->isCreatedFace(104,cell)) ;
  CHECK(cell == 60) ;
  REQUIRE(remap->isRemovedFace(40,cell)) ;
  CHECK(cell == 70) ;

  const std::vector<double> sourceAverages{4.0,2.0,6.0,8.0,10.0} ;
  std::vector<double> targetAverages ;
  std::vector<unsigned char> mapped ;
  REQUIRE(remap->remapFaceAverages(
    sourceAverages,targetAverages,mapped,true)) ;
  CHECK(targetAverages[0] == doctest::Approx(4.0)) ;
  CHECK(targetAverages[1] == doctest::Approx(4.0)) ;
  CHECK(targetAverages[2] == doctest::Approx(-2.0)) ;
  CHECK(targetAverages[3] == doctest::Approx(-8.0)) ;
  CHECK(mapped[4] == 0) ;

  REQUIRE(remap->remapFaceAverages(
    sourceAverages,targetAverages,mapped,false)) ;
  CHECK(targetAverages[2] == doctest::Approx(4.0)) ;
  CHECK(targetAverages[3] == doctest::Approx(8.0)) ;

  const std::vector<double> sourceIntegrals{8.0,2.0,6.0,8.0,10.0} ;
  std::vector<double> targetIntegrals ;
  REQUIRE(remap->remapFaceIntegrals(
    sourceIntegrals,targetIntegrals,mapped,true)) ;
  CHECK(targetIntegrals[0] == doctest::Approx(4.0)) ;
  CHECK(targetIntegrals[1] == doctest::Approx(4.0)) ;
  CHECK(targetIntegrals[2] == doctest::Approx(-4.0)) ;
  CHECK(targetIntegrals[3] == doctest::Approx(-8.0)) ;
  CHECK(mapped[4] == 0) ;

  SUBCASE("crossing anisotropic partitions may be many-to-many") {
    const std::vector<AMRFaceGeometry> sources{
      AMRFaceGeometry(0,1.0,vector3d<double>(0.25,0.5,0.0)),
      AMRFaceGeometry(1,1.0,vector3d<double>(0.75,0.5,0.0))
    } ;
    const std::vector<AMRFaceGeometry> targets{
      AMRFaceGeometry(10,1.0,vector3d<double>(0.5,0.25,0.0)),
      AMRFaceGeometry(11,1.0,vector3d<double>(0.5,0.75,0.0))
    } ;
    const std::vector<AMRFaceContribution> contributions{
      AMRFaceContribution(0,10,0.5,vector3d<double>(0.25,0.25,0.0),1),
      AMRFaceContribution(0,11,0.5,vector3d<double>(0.25,0.75,0.0),1),
      AMRFaceContribution(1,10,0.5,vector3d<double>(0.75,0.25,0.0),1),
      AMRFaceContribution(1,11,0.5,vector3d<double>(0.75,0.75,0.0),1)
    } ;
    CPTR<AMRFaceRemap> crossing = AMRFaceRemap::create(
      sources,targets,contributions,std::vector<AMRCreatedFace>(),
      std::vector<AMRRemovedFace>(),report) ;
    REQUIRE(has_face_remap(crossing)) ;
    REQUIRE(crossing->remapFaceAverages(
      std::vector<double>{2.0,6.0},targetAverages,mapped,false)) ;
    CHECK(targetAverages == std::vector<double>{4.0,4.0}) ;
  }
}


/// Incomplete face coverage and ambiguous contribution relations must be
/// rejected before a solver transfers face history or fluxes.
TEST_CASE("face remap rejects incomplete or ambiguous ancestry") {
  AMRFaceRemapReport report ;

  SUBCASE("a new internal face must be identified explicitly") {
    CHECK_FALSE(has_face_remap(AMRFaceRemap::create(
      source_faces(),target_faces(),face_contributions(),
      std::vector<AMRCreatedFace>(),
      std::vector<AMRRemovedFace>{AMRRemovedFace(40,70)},report))) ;
    CHECK(report.missingTargetFaces == 1) ;
  }

  SUBCASE("orientation is a sign") {
    std::vector<AMRFaceContribution> contributions = face_contributions() ;
    contributions[0].orientation = 0 ;
    CHECK_FALSE(has_face_remap(AMRFaceRemap::create(
      source_faces(),target_faces(),contributions,
      std::vector<AMRCreatedFace>{AMRCreatedFace(104,60)},
      std::vector<AMRRemovedFace>{AMRRemovedFace(40,70)},report))) ;
    CHECK(report.invalidGeometry > 0) ;
  }
}


/// Node contributor weights should reproduce geometry and any linear nodal
/// field for retained, edge-created, face-created, and cell-created nodes.
TEST_CASE("node remap reproduces geometry and linear nodal data") {
  AMRNodeRemapReport report ;
  CPTR<AMRNodeRemap> remap = AMRNodeRemap::create(
    source_nodes(),target_nodes(),node_contributions(),node_origins(),report) ;
  REQUIRE(has_node_remap(remap)) ;
  REQUIRE(report.valid) ;
  CHECK(report.maximumWeightError == doctest::Approx(0.0)) ;
  CHECK(report.maximumPositionError == doctest::Approx(0.0)) ;

  AMRNodeOrigin origin ;
  REQUIRE(remap->nodeOrigin(11,origin)) ;
  CHECK(origin.kind == amr_node_origin::edge) ;
  CHECK(origin.sourceEntity == 100) ;

  const std::vector<AMRNodeGeometry> sources = source_nodes() ;
  const std::vector<AMRNodeGeometry> targets = target_nodes() ;
  std::vector<double> sourceValues ;
  for(size_t source=0;source<sources.size();++source) {
    const vector3d<double>& point = sources[source].position ;
    sourceValues.push_back(1.0+2.0*point.x+3.0*point.y+4.0*point.z) ;
  }
  std::vector<double> targetValues ;
  REQUIRE(remap->interpolateNodeData(sourceValues,targetValues)) ;
  for(size_t target=0;target<targets.size();++target) {
    const vector3d<double>& point = targets[target].position ;
    CHECK(targetValues[target] == doctest::Approx(
      1.0+2.0*point.x+3.0*point.y+4.0*point.z)) ;
  }

  std::vector<vector3d<double> > sourceVectors ;
  for(size_t source=0;source<sources.size();++source)
    sourceVectors.push_back(sources[source].position) ;
  std::vector<vector3d<double> > targetVectors ;
  REQUIRE(remap->interpolateNodeData(sourceVectors,targetVectors)) ;
  for(size_t target=0;target<targets.size();++target)
    CHECK(norm(targetVectors[target]-targets[target].position) ==
          doctest::Approx(0.0)) ;
}


/// Node plans must reject weights or retained-node identities that do not
/// reproduce the stated target geometry.
TEST_CASE("node remap rejects inconsistent provenance") {
  AMRNodeRemapReport report ;

  SUBCASE("weights do not reproduce the target node") {
    std::vector<AMRNodeContribution> contributions = node_contributions() ;
    contributions[1].weight = 0.25 ;
    CHECK_FALSE(has_node_remap(AMRNodeRemap::create(
      source_nodes(),target_nodes(),contributions,node_origins(),report))) ;
    CHECK(report.inconsistentWeights > 0) ;
    CHECK(report.inconsistentPositions > 0) ;
  }

  SUBCASE("retained origin names the actual source node") {
    std::vector<AMRNodeOrigin> origins = node_origins() ;
    origins[0].sourceEntity = 1 ;
    CHECK_FALSE(has_node_remap(AMRNodeRemap::create(
      source_nodes(),target_nodes(),node_contributions(),origins,report))) ;
    CHECK(report.invalidOrigins > 0) ;
  }

  SUBCASE("origin kind names a supported source entity") {
    std::vector<AMRNodeOrigin> origins = node_origins() ;
    origins[1].kind = static_cast<amr_node_origin::value>(99) ;
    CHECK_FALSE(has_node_remap(AMRNodeRemap::create(
      source_nodes(),target_nodes(),node_contributions(),origins,report))) ;
    CHECK(report.invalidOrigins > 0) ;
  }
}

//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//# This file is part of the Loci Framework, distributed under the LGPL.
//#
//#############################################################################

#include <Loci.h>
#include <LociGridReaders.h>
#include <FVMAdapt2/dataxferDB.h>
#include <FVMAdapt2/gridInterface.h>
#include <FVMAdapt2/mesh_transfer.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <cmath>
#include <map>
#include <set>

using namespace Loci ;

namespace {
  // Make every rank leave a failed stage together before the next collective.
  bool all_pass(bool passed) {
    int local = passed ? 1 : 0 ;
    int global = 0 ;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
    return global != 0 ;
  }

  std::set<FaceId> global_ids(const store<FaceId>& ids) {
    std::vector<FaceId> local ;
    FORALL(ids.domain(),entity) {
      local.push_back(ids[entity]) ;
    } ENDFORALL ;
    int count = int(local.size()) ;
    std::vector<int> counts(MPI_processes), offsets(MPI_processes,0) ;
    MPI_Allgather(&count,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD) ;
    for(int rank=1;rank<MPI_processes;++rank)
      offsets[rank] = offsets[rank-1]+counts[rank-1] ;
    std::vector<FaceId> global(offsets.back()+counts.back()) ;
    MPI_Allgatherv(local.data(),count,MPI_LONG_LONG,global.data(),
                   counts.data(),offsets.data(),MPI_LONG_LONG,MPI_COMM_WORLD) ;
    return std::set<FaceId>(global.begin(),global.end()) ;
  }

  // The fixture's four boundary planes are x=0, y=0, z=0, and x+y+z=1.
  int boundary_plane(const vector3d<double>& point) {
    if(std::abs(point.x) < 1.e-10) return 0 ;
    if(std::abs(point.y) < 1.e-10) return 1 ;
    if(std::abs(point.z) < 1.e-10) return 2 ;
    if(std::abs(point.x+point.y+point.z-1.0) < 1.e-10) return 3 ;
    return -1 ;
  }

  CPTR<FaceRemap> installed_remap(refinedGridData& grid, fact_db& facts) {
    INFO("face status = " << int(grid.faceTransitionReport.status)
         << ", invalid polygons = " << grid.faceTransitionReport.invalidPolygons) ;
    REQUIRE(all_pass(grid.faceTransitionReport.valid)) ;
    REQUIRE(all_pass(setupFVMGridFromContainer(facts,grid))) ;
    storeRepP rep = facts.get_variable("faceRemap") ;
    REQUIRE(all_pass(rep != 0)) ;
    blackbox<CPTR<FaceRemap> > remap ;
    remap = rep ;
    REQUIRE(all_pass(*remap != static_cast<FaceRemap*>(0))) ;
    store<FaceId> faceIds ;
    faceIds = facts.get_variable("faceId") ;
    constraint faces ;
    faces = facts.get_variable("faces") ;
    CHECK(faceIds.domain() == *faces) ;
    CHECK((*remap)->targetFaceGeometry().size() == faceIds.domain().size()) ;
    FORALL(*faces,face) {
      size_t begin = 0, end = 0 ;
      CellId origin = 0 ;
      CHECK(((*remap)->overlaps(faceIds[face],begin,end) && begin < end)
            != (*remap)->isCreatedFace(faceIds[face],origin)) ;
    } ENDFORALL ;
    return *remap ;
  }

  // Check the public overlap rows against the tetrahedron's geometry.
  void check_boundary_overlaps(const FaceRemap& remap,
                                   const std::set<FaceId>& oldFaces,
                                   bool refining) {
    std::map<FaceId,FaceGeometry> source ;
    for(const FaceGeometry& face : remap.sourceFaceGeometry())
      source[face.face] = face ;
    for(const FaceGeometry& face : remap.targetFaceGeometry()) {
      size_t begin = 0, end = 0 ;
      const bool found = remap.overlaps(face.face,begin,end) ;
      CHECK(found) ;
      if(!found) continue ;
      const int plane = boundary_plane(face.centroid) ;
      if(plane < 0) {
        CHECK(begin == end) ;
        continue ;
      }
      if(refining) CHECK(end-begin == 1) ;
      else CHECK(end-begin > 1) ;
      double area = 0.0 ;
      for(size_t row=begin;row<end;++row) {
        const FaceOverlap& contribution = remap.overlaps()[row] ;
        CHECK(oldFaces.count(contribution.source) == 1) ;
        CHECK(contribution.target == face.face) ;
        CHECK(contribution.orientation == 1) ;
        CHECK(source.count(contribution.source) == 1) ;
        if(source.count(contribution.source) != 1) continue ;
        const FaceGeometry& old = source.at(contribution.source) ;
        CHECK(boundary_plane(old.centroid) == plane) ;
        CHECK(contribution.area > 0.0) ;
        CHECK(contribution.area <= old.area*(1.0+1.e-10)) ;
        CHECK(boundary_plane(contribution.centroid) == plane) ;
        area += contribution.area ;
        if(MPI_processes == 1)
          std::cout << "face " << face.face << ": source=" << old.face
                    << ", area=" << contribution.area
                    << ", orientation=" << contribution.orientation
                    << ", sourceFraction=" << contribution.area/old.area
                    << '\n' ;
      }
      CHECK(area == doctest::Approx(face.area)) ;
    }
  }
}

/// Split one tetrahedron: boundary faces inherit old-face flux, while every
/// new internal face identifies the old cell in which it was created.
/// Coarsen it again: boundary contributions recombine and internal faces vanish.
TEST_CASE("general-cell faces expose old contributors and originating cells") {
  bool replay = false ;
  SUBCASE("retain the previous mesh state") {}
  SUBCASE("replay the previous mesh from its plan") { replay = true ; }

  DataXFER_DB.deleteItem("currentPlan") ;
  DataXFER_DB.deleteItem("nextPlan") ;
  DataXFER_DB.deleteItem("cellweights") ;
  DataXFER_DB.deleteItem("c2pglobal") ;

  fact_db original ;
  REQUIRE(all_pass(setupFVMGrid(original,"tet.vog"))) ;
  rule_db rules, consumerRules ;
  load_module("fvmadapt2",rules) ;
  load_module("fvm",consumerRules) ;
  REQUIRE(all_pass(makeQuery(consumerRules,original,
    "fileNumber(pos),fileNumber(face2node),fileNumber(geom_cells)"))) ;
  REQUIRE(all_pass(installBaseMeshIds(original))) ;
  store<FaceId> oldFaceIds ;
  oldFaceIds = original.get_variable("faceId") ;
  const std::set<FaceId> baseFaces = global_ids(oldFaceIds) ;
  store<CellId> oldCellIds ;
  oldCellIds = original.get_variable("cellId") ;
  const std::set<CellId> baseCells = global_ids(oldCellIds) ;
  REQUIRE(baseFaces.size() == 4) ;
  REQUIRE(baseCells.size() == 1) ;

  store<char> tags ;
  tags.allocate(oldCellIds.domain()) ;
  FORALL(tags.domain(),cell) { tags[cell] = 1 ; } ENDFORALL ;
  CPTR<refinedGridData> refined ;
  CPTR<MeshState> history ;
  onlineRefineMesh(refined,history,rules,2,0,tags.Rep(),"tet") ;
  REQUIRE(all_pass(refined != static_cast<refinedGridData*>(0))) ;
  fact_db refinedFacts ;
  const CPTR<FaceRemap> split = installed_remap(*refined,refinedFacts) ;
  CHECK(refined->faceTransitionReport.remap.createdFaces > 0) ;
  CHECK(refined->faceTransitionReport.remap.removedFaces == 0) ;
  check_boundary_overlaps(*split,baseFaces,true) ;
  for(const CreatedFace& face : split->createdFaces()) {
    CellId sourceCell = 0 ;
    CHECK(split->isCreatedFace(face.targetFace,sourceCell)) ;
    CHECK(sourceCell == *baseCells.begin()) ;
    if(MPI_processes == 1)
      std::cout << "created face " << face.targetFace << ": sourceCell="
                << sourceCell << '\n' ;
  }

  // A unit outward flux density has an integrated flux equal to face area.
  std::vector<double> oldFlux, newFlux ;
  std::vector<unsigned char> mapped ;
  for(const FaceGeometry& face : split->sourceFaceGeometry())
    oldFlux.push_back(face.area) ;
  REQUIRE(all_pass(split->remapFaceIntegrals(oldFlux,newFlux,mapped,true))) ;
  for(size_t face=0;face<newFlux.size();++face) {
    const FaceGeometry& target = split->targetFaceGeometry()[face] ;
    CHECK(bool(mapped[face]) == (boundary_plane(target.centroid) >= 0)) ;
    if(mapped[face]) CHECK(newFlux[face] == doctest::Approx(target.area)) ;
  }

  oldFaceIds = refinedFacts.get_variable("faceId") ;
  const std::set<FaceId> refinedFaces = global_ids(oldFaceIds) ;
  const size_t internalFaces = refined->faceTransitionReport.remap.createdFaces ;
  tags.allocate(refined->cellIds.domain()) ;
  // Tag 0 retains cells. Even a warped internal face must retain its identity
  // and receive its entire contribution from that same previous face.
  FORALL(tags.domain(),cell) { tags[cell] = 0 ; } ENDFORALL ;
  CPTR<refinedGridData> retained ;
  onlineRefineMesh(retained,history,rules,2,1,tags.Rep(),"tet") ;
  REQUIRE(all_pass(retained != static_cast<refinedGridData*>(0))) ;
  fact_db retainedFacts ;
  const CPTR<FaceRemap> unchanged = installed_remap(*retained,retainedFacts) ;
  CHECK(retained->faceTransitionReport.remap.createdFaces == 0) ;
  CHECK(retained->faceTransitionReport.remap.removedFaces == 0) ;
  for(const FaceGeometry& face : unchanged->targetFaceGeometry()) {
    size_t begin = 0, end = 0 ;
    CHECK(unchanged->overlaps(face.face,begin,end)) ;
    CHECK(end-begin == 1) ;
    if(end-begin != 1) continue ;
    const FaceOverlap& contribution = unchanged->overlaps()[begin] ;
    CHECK(contribution.source == face.face) ;
    CHECK(contribution.orientation == 1) ;
    CHECK(contribution.area == doctest::Approx(face.area)) ;
  }

  // The module's cell tags use 2 for derefinement.
  tags.allocate(retained->cellIds.domain()) ;
  FORALL(tags.domain(),cell) { tags[cell] = 2 ; } ENDFORALL ;
  if(replay) history = CPTR<MeshState>() ;
  CPTR<refinedGridData> coarsened ;
  onlineRefineMesh(coarsened,history,rules,2,2,tags.Rep(),"tet") ;
  REQUIRE(all_pass(coarsened != static_cast<refinedGridData*>(0))) ;
  fact_db coarsenedFacts ;
  const CPTR<FaceRemap> merged = installed_remap(*coarsened,coarsenedFacts) ;
  CHECK(coarsened->faceTransitionReport.remap.createdFaces == 0) ;
  CHECK(coarsened->faceTransitionReport.remap.removedFaces == internalFaces) ;
  check_boundary_overlaps(*merged,refinedFaces,false) ;
  oldFlux.clear() ;
  for(const FaceGeometry& face : merged->sourceFaceGeometry())
    oldFlux.push_back(face.area) ;
  REQUIRE(all_pass(merged->remapFaceIntegrals(oldFlux,newFlux,mapped,true))) ;
  for(size_t face=0;face<newFlux.size();++face) {
    CHECK(mapped[face] != 0) ;
    CHECK(newFlux[face] == doctest::Approx(merged->targetFaceGeometry()[face].area)) ;
  }
  for(const RemovedFace& face : merged->removedFaces())
    CHECK(face.targetCell == *baseCells.begin()) ;
  oldFaceIds = coarsenedFacts.get_variable("faceId") ;
  CHECK(global_ids(oldFaceIds) == baseFaces) ;
  oldCellIds = coarsenedFacts.get_variable("cellId") ;
  CHECK(global_ids(oldCellIds) == baseCells) ;
}

int main(int argc, char** argv) {
  Loci::Init(&argc,&argv) ;
  doctest::Context context(argc,argv) ;
  const int result = context.run() ;
  const bool passed = all_pass(result == 0) ;
  Loci::Finalize() ;
  return passed ? 0 : 1 ;
}

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
#include <LociGridReaders.h>
#include <FVMAdapt2/dataxferDB.h>
#include "mesh_state.h"
#include <FVMAdapt2/gridInterface.h>

#include <hdf5.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace Loci ;

namespace {
  CPTR<FaceState> internal_face_state(const CPTR<MeshState>& state) {
    if (state == static_cast<MeshState*>(0))
      return CPTR<FaceState>() ;
    return CPTR<FaceState>(state) ;
  }

  bool all_ranks_pass(bool localPassed) {
    int local = localPassed ? 1 : 0 ;
    int global = 0 ;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    return global != 0 ;
  }

  std::vector<FaceId> gather_ids(const std::vector<FaceId>& local) {
    std::vector<int> counts(MPI_processes, 0) ;
    std::vector<int> offsets(MPI_processes, 0) ;
    const int localCount = int(local.size()) ;
    MPI_Allgather(
          &localCount, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    for (int process = 1; process < MPI_processes; ++process)
      offsets[process] = offsets[process - 1] + counts[process - 1] ;
    const int globalCount = offsets.back() + counts.back() ;
    std::vector<FaceId> global(static_cast<size_t>(globalCount)) ;
    MPI_Allgatherv(local.empty() ? 0 : &local[0], localCount, MPI_LONG_LONG,
          global.empty() ? 0 : &global[0], &counts[0], &offsets[0],
          MPI_LONG_LONG, MPI_COMM_WORLD) ;
    std::sort(global.begin(), global.end()) ;
    return global ;
  }

  std::vector<FaceId> gather_store_ids(storeRepP idsRep) {
    const bool available = all_ranks_pass(idsRep != 0) ;
    if (!available)
      return std::vector<FaceId>() ;
    const_store<FaceId> ids(idsRep) ;
    std::vector<FaceId> local(ids.domain().size()) ;
    size_t entry = 0 ;
    FORALL(ids.domain(), entity) {
      local[entry++] = ids[entity] ;
    }
    ENDFORALL ;
    return gather_ids(local) ;
  }

  bool unique_ids(const std::vector<FaceId>& ids) {
    return std::adjacent_find(ids.begin(), ids.end()) == ids.end() ;
  }

  bool write_restart_plan(fact_db& sourceFacts, const std::string& filename,
        int& expectedGeneratedCells, bool refineFirstCell) {
    storeRepP cellsRep = sourceFacts.get_variable("geom_cells") ;
    if (!all_ranks_pass(cellsRep != 0))
      return false ;
    constraint cells ;
    cells = cellsRep ;
    int localCellCount = int((*cells).size()) ;
    int globalCellCount = 0 ;
    MPI_Allreduce(&localCellCount, &globalCellCount, 1, MPI_INT, MPI_SUM,
          MPI_COMM_WORLD) ;
    if (globalCellCount == 0)
      return false ;

    fact_db::distribute_infoP distribution = sourceFacts.get_distribute_info() ;
    dMap globalToFile ;
    const size_t cellKeySpace = cellsRep->getDomainKeySpace() ;
    if (distribution != 0 && cellKeySpace < distribution->g2fv.size())
      globalToFile = distribution->g2fv[cellKeySpace].Rep() ;
    const bool fileMapCovers =
          distribution != 0 && ((*cells) - globalToFile.domain()).size() == 0 ;
    if (!all_ranks_pass(distribution == 0 || fileMapCovers))
      return false ;

    int localFirstFile = std::numeric_limits<int>::max() ;
    FORALL(*cells, cell) {
      localFirstFile =
            std::min(localFirstFile, fileMapCovers ? globalToFile[cell] : cell) ;
    }
    ENDFORALL ;
    int firstFile = 0 ;
    MPI_Allreduce(
          &localFirstFile, &firstFile, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    expectedGeneratedCells = globalCellCount + (refineFirstCell ? 7 : 0) ;

    // Code 7 splits one hexahedral root in all three directions. Other roots,
    // and the general-cell case, retain an explicit empty plan.
    entitySet serializedDomain ;
    FORALL(*cells, cell) {
      serializedDomain += fileMapCovers ? globalToFile[cell] : cell ;
    }
    ENDFORALL ;
    store<std::vector<char>> serializedPlan ;
    serializedPlan.allocate(serializedDomain) ;
    FORALL(*cells, cell) {
      const int file = fileMapCovers ? globalToFile[cell] : cell ;
      serializedPlan[file] = std::vector<char>(
            1, refineFirstCell && file == firstFile ? char(7) : 'C') ;
    }
    ENDFORALL ;

    hid_t fileId = Loci::hdf5CreateFile(
          filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) ;
    hid_t groupId = 0 ;
    if (MPI_rank == 0)
      groupId = H5Gcreate(
            fileId, "cellPlan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
    int offset = 0 ;
    Loci::write_store(groupId, serializedPlan.Rep(), serializedDomain, offset,
          MPI_COMM_WORLD) ;
    if (MPI_rank == 0)
      H5Gclose(groupId) ;

    param<int> refineLevel ;
    *refineLevel = 1 ;
    Loci::writeContainer(fileId, "refineLevel", refineLevel.Rep(), sourceFacts) ;
    Loci::hdf5CloseFile(fileId) ;
    return true ;
  }

  bool generated_ids_are_complete(const CPTR<refinedGridData>& grid,
        std::vector<FaceId>& faceIds, std::vector<FaceId>& cellIds,
        int expectedGeneratedCells) {
    if (!all_ranks_pass(grid != static_cast<refinedGridData*>(0)))
      return false ;
    const bool localPartitionAvailable =
          grid->local_cells.size() == size_t(MPI_processes) ;
    if (!all_ranks_pass(localPartitionAvailable))
      return false ;

    const entitySet faces = grid->new_face2node.domain() ;
    const entitySet cells = grid->local_cells[MPI_rank] ;
    const CPTR<FaceState> faceState =
          internal_face_state(grid->transitionState) ;
    bool localValid = faceState != static_cast<FaceState*>(0) &&
                      grid->faceIds.domain() == faces &&
                      grid->cellIds.domain() == cells ;
    if (localValid) {
      const std::vector<FaceIdentity>& identities = faceState->faceIdentities() ;
      localValid = identities.size() == faces.size() ;
      for (size_t face = 0; face < identities.size(); ++face)
        localValid =
              localValid && faces.inSet(identities[face].face) &&
              grid->faceIds[identities[face].face] == identities[face].id ;
    }
    if (!all_ranks_pass(localValid))
      return false ;

    faceIds = gather_store_ids(grid->faceIds.Rep()) ;
    cellIds = gather_store_ids(grid->cellIds.Rep()) ;
    return !faceIds.empty() &&
           cellIds.size() == static_cast<size_t>(expectedGeneratedCells) &&
           unique_ids(faceIds) && unique_ids(cellIds) ;
  }

  bool installed_ids_are_complete(fact_db& facts,
        const std::vector<FaceId>& expectedFaces,
        const std::vector<FaceId>& expectedCells) {
    storeRepP facesRep = facts.get_variable("faces") ;
    storeRepP cellsRep = facts.get_variable("geom_cells") ;
    storeRepP faceIdsRep = facts.get_variable("faceId") ;
    storeRepP cellIdsRep = facts.get_variable("cellId") ;
    if (!all_ranks_pass(facesRep != 0 && cellsRep != 0 && faceIdsRep != 0 &&
                        cellIdsRep != 0))
      return false ;

    constraint faces ;
    constraint cells ;
    faces = facesRep ;
    cells = cellsRep ;
    const_store<FaceId> faceIds(faceIdsRep) ;
    const_store<CellId> cellIds(cellIdsRep) ;
    if (!all_ranks_pass(
              faceIds.domain() == *faces && cellIds.domain() == *cells))
      return false ;
    return gather_store_ids(faceIdsRep) == expectedFaces &&
           gather_store_ids(cellIdsRep) == expectedCells ;
  }

  std::vector<FaceId> transition_source_ids(const CPTR<FaceRemap>& remap) {
    const bool available = all_ranks_pass(remap != static_cast<FaceRemap*>(0)) ;
    if (!available)
      return std::vector<FaceId>() ;
    const std::vector<FaceGeometry>& geometry = remap->sourceFaceGeometry() ;
    std::vector<FaceId> local(geometry.size()) ;
    for (size_t face = 0; face < geometry.size(); ++face)
      local[face] = geometry[face].face ;
    return gather_ids(local) ;
  }

  int finish_test(bool passed, const std::string& planFile,
        const std::string& missingWeightFile, const char* description) {
    Loci::DataXFER_DB.deleteItem("currentPlan") ;
    Loci::DataXFER_DB.deleteItem("nextPlan") ;
    Loci::DataXFER_DB.deleteItem("cellweights") ;
    if (MPI_rank == 0) {
      std::remove(planFile.c_str()) ;
      std::remove(missingWeightFile.c_str()) ;
    }
    MPI_Barrier(MPI_COMM_WORLD) ;
    if (MPI_rank == 0)
      std::cout << description << (passed ? " passed" : " failed") << std::endl ;
    Loci::Finalize() ;
    return passed ? 0 : 1 ;
  }
}

int main(int argc, char* argv[]) {
  Loci::Init(&argc, &argv) ;

  const bool generalCell = argc == 3 && std::string(argv[1]) == "--general" ;
  bool passed = argc == 2 || generalCell ;
  if (!passed && MPI_rank == 0)
    std::cerr << "usage: test_plan_restart [--general] <case>" << std::endl ;
  const std::string caseName =
        passed ? argv[generalCell ? 2 : 1] : std::string() ;
  const std::string planFile = "plan_restart.h5" ;
  const std::string missingWeightFile = "missing_plan_restart_weights.h5" ;

  if (MPI_rank == 0) {
    std::remove(planFile.c_str()) ;
    std::remove(missingWeightFile.c_str()) ;
  }
  MPI_Barrier(MPI_COMM_WORLD) ;
  Loci::DataXFER_DB.deleteItem("currentPlan") ;
  Loci::DataXFER_DB.deleteItem("nextPlan") ;

  fact_db sourceFacts ;
  if (passed && !Loci::setupFVMGrid(sourceFacts, caseName + ".vog")) {
    if (MPI_rank == 0)
      std::cerr << "Unable to read the source grid" << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed) ;

  int expectedGeneratedCells = 0 ;
  if (passed && !write_restart_plan(sourceFacts, planFile,
                      expectedGeneratedCells, !generalCell)) {
    if (MPI_rank == 0)
      std::cerr << "Unable to write the nontrivial restart plan" << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed) ;

  rule_db refmeshRules ;
  if (passed)
    Loci::load_module("fvmadapt2", refmeshRules) ;

  int level = 0 ;
  CPTR<refinedGridData> restartedGrid ;
  if (passed)
    Loci::initializeGridFromPlan(restartedGrid, level, refmeshRules, caseName,
          missingWeightFile, planFile) ;

  // Face identities and remaps are independent of general-node support.
  if (generalCell)
    passed = all_ranks_pass(
          passed && restartedGrid != static_cast<refinedGridData*>(0) &&
          restartedGrid->nodeTransitionReport.status ==
                node_transition_status::unsupported_topology &&
          restartedGrid->nodeIds.domain() == EMPTY) ;

  std::vector<FaceId> generatedFaceIds ;
  std::vector<FaceId> generatedCellIds ;
  if (passed && !generated_ids_are_complete(restartedGrid, generatedFaceIds,
                      generatedCellIds, expectedGeneratedCells)) {
    if (MPI_rank == 0)
      std::cerr << "Plan restart did not produce complete persistent IDs"
                << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed && level == 1) ;

  fact_db installedFacts ;
  if (passed &&
        !Loci::setupFVMGridFromContainer(installedFacts, *restartedGrid, 0)) {
    if (MPI_rank == 0)
      std::cerr << "Unable to install the plan-restarted grid" << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed) ;

  if (passed && !installed_ids_are_complete(
                      installedFacts, generatedFaceIds, generatedCellIds)) {
    if (MPI_rank == 0)
      std::cerr << "Installed plan-restart IDs are incomplete or changed"
                << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed) ;

  store<char> deeperTags ;
  CPTR<MeshState> retainedTransitionState ;
  if (passed) {
    deeperTags.allocate(restartedGrid->cellIds.domain()) ;
    FORALL(deeperTags.domain(), cell) {
      deeperTags[cell] = 1 ;
    }
    ENDFORALL ;
    retainedTransitionState = restartedGrid->transitionState ;
    restartedGrid = CPTR<refinedGridData>() ;
  }

  CPTR<refinedGridData> nextGrid ;
  if (passed)
    Loci::onlineRefineMesh(nextGrid, retainedTransitionState, refmeshRules, 2,
          level, deeperTags.Rep(), caseName) ;

  if (passed) {
    const bool transitionValid =
          nextGrid != static_cast<refinedGridData*>(0) &&
          nextGrid->faceRemap != static_cast<FaceRemap*>(0) &&
          nextGrid->faceTransitionReport.valid &&
          nextGrid->faceTransitionReport.status ==
                face_transition_status::available &&
          nextGrid->faceTransitionReport.remap.valid ;
    passed = all_ranks_pass(transitionValid) ;
  }
  if (passed &&
        transition_source_ids(nextGrid->faceRemap) != generatedFaceIds) {
    if (MPI_rank == 0)
      std::cerr << "The next transition did not use the plan-restarted face "
                   "identities as its source"
                << std::endl ;
    passed = false ;
  }
  passed = all_ranks_pass(passed) ;

  return finish_test(
        passed, planFile, missingWeightFile, "plan restart identity test") ;
}

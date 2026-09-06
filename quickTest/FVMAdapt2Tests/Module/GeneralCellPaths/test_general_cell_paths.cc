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
#include <FVMAdapt2/diamondcell.h>
#include "mesh_state.h"
#include <FVMAdapt2/gridInterface.h>

#include "library/refinement_state_internal.h"

#include <iostream>
#include <string>
#include <vector>

using namespace Loci ;

namespace {

  bool all_ranks_pass(bool localPassed) {
    const int local = localPassed ? 1 : 0 ;
    int global = 0 ;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    return global != 0 ;
  }

  /// Seed tetrahedral root children without introducing geometric state.
  void seed_tetra_children(Cell& root) {
    root.numNode = 4 ;
    root.child = new DiamondCell*[root.numNode] ;
    for (int child = 0; child < root.numNode; ++child)
      root.child[child] = new DiamondCell(3) ;
  }

  bool replay_result(const std::vector<char>& plan, bool seedChildren,
        bool expectedValid, int expectedLeaves) {
    Cell root ;
    root.numNode = 4 ;
    if (seedChildren)
      seed_tetra_children(root) ;
    int leaves = -1 ;
    const bool valid = Loci::detail::replayGeneralCellPlan(&root, plan, leaves) ;
    return valid == expectedValid && (!valid || leaves == expectedLeaves) ;
  }

  /// Accept compact or zero-padded plans and reject unconsumed plan data.
  bool check_plan_validation() {
    return replay_result(std::vector<char>(), false, true, 1) &&
           replay_result(std::vector<char>(1, 1), true, true, 4) &&
           replay_result(std::vector<char>{1, 0, 0, 0, 0, 0}, true, true, 4) &&
           replay_result(std::vector<char>{1, 0, 1}, true, true, 11) &&
           replay_result(std::vector<char>(1, 0), false, false, 0) &&
           replay_result(std::vector<char>(1, 2), false, false, 0) &&
           replay_result(std::vector<char>{1, 0, 0, 0, 0, 1}, true, false, 0) ;
  }

  bool check_paths(rule_db& rules, const std::string& meshFile,
        const std::vector<char>& plan,
        const std::vector<std::vector<int>>& expectedPaths) {
    fact_db facts ;
    bool valid = Loci::setupFVMGrid(facts, meshFile) ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    Loci::createLowerUpper(facts) ;
    Loci::createEdgesPar(facts) ;
    Loci::parallelClassifyCell(facts) ;

    storeRepP cellsRep = facts.get_variable("geom_cells") ;
    storeRepP generalCellsRep = facts.get_variable("gnrlcells") ;
    valid = cellsRep != 0 && generalCellsRep != 0 ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    constraint cells ;
    constraint generalCells ;
    cells = cellsRep ;
    generalCells = generalCellsRep ;
    valid = *cells == *generalCells ;
    int localCellCount = int((*cells).size()) ;
    int globalCellCount = 0 ;
    MPI_Allreduce(&localCellCount, &globalCellCount, 1, MPI_INT, MPI_SUM,
          MPI_COMM_WORLD) ;
    valid = all_ranks_pass(valid) && globalCellCount == 1 ;
    if (!valid)
      return false ;

    store<std::vector<char>> acceptedPlan ;
    acceptedPlan.allocate(*cells) ;
    FORALL(*cells, cell) {
      acceptedPlan[cell] = plan ;
    }
    ENDFORALL ;
    facts.create_fact("balancedCellPlan", acceptedPlan) ;

    valid = Loci::makeQuery(rules, facts, "balanced_num_fine_cells") ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    storeRepP leafCountRep = facts.get_variable("balanced_num_fine_cells") ;
    valid = leafCountRep != 0 ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    const_store<int> leafCount(leafCountRep) ;
    valid = ((*cells) - leafCount.domain()) == EMPTY ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    FORALL(*cells, cell) {
      if (leafCount[cell] != int(expectedPaths.size()))
        valid = false ;
    }
    ENDFORALL ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    valid = Loci::makeQuery(rules, facts, "cellLeafPaths") ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    storeRepP encodedPathsRep = facts.get_variable("cellLeafPaths") ;
    valid = encodedPathsRep != 0 ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    const_store<std::vector<int>> encodedPaths(encodedPathsRep) ;
    valid = ((*cells) - encodedPaths.domain()) == EMPTY ;
    valid = all_ranks_pass(valid) ;
    if (!valid)
      return false ;

    FORALL(*cells, cell) {
      std::vector<std::vector<int>> paths ;
      if (!Loci::detail::decodeLeafPaths(encodedPaths[cell], paths) ||
            paths != expectedPaths)
        valid = false ;
    }
    ENDFORALL ;
    return all_ranks_pass(valid) ;
  }

}

int main(int argc, char* argv[]) {
  Loci::Init(&argc, &argv) ;

  bool passed = argc == 2 ;
  if (!passed && MPI_rank == 0)
    std::cerr << "usage: test_general_cell_paths <case>" << std::endl ;

  if (passed)
    passed = check_plan_validation() ;

  rule_db rules ;
  if (passed)
    Loci::load_module("fvmadapt2", rules) ;

  const std::vector<std::vector<int>> rootPath(1) ;
  if (passed)
    passed = check_paths(
          rules, std::string(argv[1]) + ".vog", std::vector<char>(), rootPath) ;

  std::vector<std::vector<int>> splitPaths(4) ;
  for (int child = 0; child < 4; ++child) {
    splitPaths[child].push_back(1) ;
    splitPaths[child].push_back(child) ;
  }
  if (passed)
    passed = check_paths(rules, std::string(argv[1]) + ".vog",
          std::vector<char>(1, 1), splitPaths) ;

  std::vector<std::vector<int>> nestedPaths{{1, 0}, {1, 2}, {1, 3}} ;
  for (int child = 0; child < 8; ++child)
    nestedPaths.push_back(std::vector<int>{1, 1, 1, child}) ;
  if (passed)
    passed = check_paths(rules, std::string(argv[1]) + ".vog",
          std::vector<char>{1, 0, 1}, nestedPaths) ;

  passed = all_ranks_pass(passed) ;
  if (MPI_rank == 0)
    std::cout << "FVMAdapt2 general-cell leaf paths: "
              << (passed ? "PASSED" : "FAILED") << std::endl ;

  Loci::Finalize() ;
  return passed ? 0 : 1 ;
}

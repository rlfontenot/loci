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
#ifndef FVMADAPT2_REMAP_PLAN_H
#define FVMADAPT2_REMAP_PLAN_H

#include <FVMAdapt2/mesh_transfer.h>

#include <Tools/basic_types.h>
#include <Tools/cptr.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace Loci {

  /// Geometry for one cell in either the source or target index space.
  ///
  /// reconstructionPoint is the point about which a source-cell average is
  /// reconstructed. It normally equals centroid; a refined source may use the
  /// volume-weighted centroid of its target contributions to preserve its
  /// integral exactly under linear reconstruction.
  struct AMRCellGeometry {
    /// Transition-local index used by the existing interpolation schedules.
    int cell ;
    /// Persistent identity independent of entity numbering and MPI ownership.
    CellId cellId ;
    double volume ;
    vector3d<double> centroid ;
    vector3d<double> reconstructionPoint ;

    AMRCellGeometry()
        : cell(0), cellId(0), volume(0.0), centroid(0.0, 0.0, 0.0),
          reconstructionPoint(0.0, 0.0, 0.0) {}
    AMRCellGeometry(int cellIndex, CellId persistentId, double cellVolume,
          const vector3d<double>& cellCentroid)
        : cell(cellIndex), cellId(persistentId), volume(cellVolume),
          centroid(cellCentroid), reconstructionPoint(cellCentroid) {}
    AMRCellGeometry(int cellIndex, CellId persistentId, double cellVolume,
          const vector3d<double>& cellCentroid,
          const vector3d<double>& sourceReconstructionPoint)
        : cell(cellIndex), cellId(persistentId), volume(cellVolume),
          centroid(cellCentroid),
          reconstructionPoint(sourceReconstructionPoint) {}
  } ;

  /// Geometric contribution from one source cell to one target cell.
  struct AMRCellContribution {
    /// Transition-local source and target indices.
    int sourceCell ;
    int targetCell ;
    /// Persistent source and target identities for public consumers.
    CellId sourceCellId ;
    CellId targetCellId ;
    double overlapVolume ;
    vector3d<double> overlapCentroid ;

    AMRCellContribution()
        : sourceCell(0), targetCell(0), sourceCellId(0), targetCellId(0),
          overlapVolume(0.0), overlapCentroid(0.0, 0.0, 0.0) {}
    AMRCellContribution(int sourceIndex, int targetIndex, CellId sourceId,
          CellId targetId, double volume, const vector3d<double>& centroid)
        : sourceCell(sourceIndex), targetCell(targetIndex),
          sourceCellId(sourceId), targetCellId(targetId), overlapVolume(volume),
          overlapCentroid(centroid) {}
  } ;

  namespace amr_cell_transition {
    enum value {
      retained,
      refined,
      derefined
    } ;
  }

  /// Structural and geometric checks performed while building a remap plan.
  struct AMRRemapReport {
    bool valid ;
    size_t sourceCells ;
    size_t targetCells ;
    size_t contributions ;
    size_t invalidGeometry ;
    size_t invalidIdentities ;
    size_t duplicateContributions ;
    size_t missingSourceCells ;
    size_t missingTargetCells ;
    size_t inconsistentSourceMoments ;
    size_t inconsistentTargetMoments ;
    size_t unsupportedRelations ;
    bool sourceCoverageChecked ;
    double maximumSourceVolumeError ;
    double maximumTargetVolumeError ;
    double maximumSourceCentroidError ;
    double maximumTargetCentroidError ;

    AMRRemapReport() ;
  } ;

  /// Immutable geometric relation between the cells of two adaptation states.
  ///
  /// Source and target cell numbers occupy distinct index spaces even when
  /// their integer values happen to match. Contributions are stored in target
  /// order and contain unnormalized overlap measures.
  class AMRRemapPlan : public CPTR_type {
  public:
    static CPTR<AMRRemapPlan>
    createCellPlan(const std::vector<AMRCellGeometry>& sourceGeometry,
                   const std::vector<AMRCellGeometry>& targetGeometry,
                   const std::vector<AMRCellContribution>& contributions,
                   AMRRemapReport& report,
                   double relativeTolerance = 1.0e-10) ;

    /// Build the target-owned portion of a distributed cell remap.
    ///
    /// sourceTargetCounts gives the global number of targets receiving data
    /// from each entry in sourceGeometry. Source-volume coverage must be
    /// checked collectively by the distributed builder; target coverage is
    /// checked here from the complete target-owned contribution ranges.
    static CPTR<AMRRemapPlan>
    createCellPlanPartition(
      const std::vector<AMRCellGeometry>& sourceGeometry,
      const std::vector<size_t>& sourceTargetCounts,
      const std::vector<AMRCellGeometry>& targetGeometry,
      const std::vector<AMRCellContribution>& contributions,
      AMRRemapReport& report,
      double relativeTolerance = 1.0e-10) ;

    const std::vector<AMRCellGeometry>& sourceCellGeometry() const {
      return sourceGeometry_ ;
    }
    const std::vector<AMRCellGeometry>& targetCellGeometry() const {
      return targetGeometry_ ;
    }
    const std::vector<AMRCellContribution>& cellContributions() const {
      return contributions_ ;
    }

    /// Return the contribution range for a target cell.
    bool cellContributions(int targetCell,
                           size_t& begin, size_t& end) const ;

    /// Return the contribution range for a persistent target-cell identity.
    bool cellContributions(
          CellId targetCellId, size_t& begin, size_t& end) const ;

    /// Find source or target geometry by persistent identity.
    const AMRCellGeometry* sourceCellGeometry(CellId sourceCellId) const ;
    const AMRCellGeometry* targetCellGeometry(CellId targetCellId) const ;

    /// Return all target cells that receive a contribution from a source cell.
    void targetCells(int sourceCell, std::vector<int>& targets) const ;

    /// Return locally owned targets receiving data from a persistent source.
    void targetCells(CellId sourceCellId, std::vector<CellId>& targets) const ;

    /// Classify a target cell from source/target relation cardinality.
    bool cellTransition(int targetCell,
                        amr_cell_transition::value& transition) const ;

    /// Classify a target cell selected by its persistent identity.
    bool cellTransition(
          CellId targetCellId, amr_cell_transition::value& transition) const ;

    /// Transfer source-cell averages with piecewise-constant reconstruction.
    bool remapCellAverages(const std::vector<double>& sourceValues,
                           std::vector<double>& targetValues) const ;

    /// Transfer source-cell averages using supplied source gradients.
    bool remapCellAverages(
      const std::vector<double>& sourceValues,
      const std::vector<vector3d<double> >& sourceGradients,
      std::vector<double>& targetValues) const ;

    /// Transfer extensive cell integrals according to geometric overlap.
    bool remapCellIntegrals(const std::vector<double>& sourceIntegrals,
                            std::vector<double>& targetIntegrals) const ;

  private:
    AMRRemapPlan() {}

    static CPTR<AMRRemapPlan>
    createCellPlanImpl(const std::vector<AMRCellGeometry>& sourceGeometry,
                       const std::vector<size_t>* sourceTargetCounts,
                       const std::vector<AMRCellGeometry>& targetGeometry,
                       const std::vector<AMRCellContribution>& contributions,
                       bool checkSourceCoverage,
                       AMRRemapReport& report,
                       double relativeTolerance) ;

    std::vector<AMRCellGeometry> sourceGeometry_ ;
    std::vector<AMRCellGeometry> targetGeometry_ ;
    std::vector<AMRCellContribution> contributions_ ;
    std::vector<size_t> targetOffsets_ ;
    std::vector<size_t> sourceDegrees_ ;
    std::vector<size_t> targetDegrees_ ;
    std::vector<std::pair<CellId, size_t>> sourceIdentityIndex_ ;
    std::vector<std::pair<CellId, size_t>> targetIdentityIndex_ ;
  } ;
}

#endif

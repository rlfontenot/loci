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

#include <FVMAdapt2/remap_plan.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace Loci {

  AMRRemapReport::AMRRemapReport()
    : valid(false), sourceCells(0), targetCells(0), contributions(0),
      invalidGeometry(0), duplicateContributions(0), missingSourceCells(0),
      missingTargetCells(0), inconsistentSourceMoments(0),
      inconsistentTargetMoments(0), unsupportedRelations(0),
      sourceCoverageChecked(false),
      maximumSourceVolumeError(0.0), maximumTargetVolumeError(0.0),
      maximumSourceCentroidError(0.0), maximumTargetCentroidError(0.0) {}

  namespace {
    bool finiteVector(const vector3d<double>& value) {
      return std::isfinite(value.x) && std::isfinite(value.y) &&
        std::isfinite(value.z) ;
    }

    bool geometryOrder(const AMRCellGeometry& left,
                       const AMRCellGeometry& right) {
      return left.cell < right.cell ;
    }

    bool contributionOrder(const AMRCellContribution& left,
                           const AMRCellContribution& right) {
      if(left.targetCell != right.targetCell)
        return left.targetCell < right.targetCell ;
      return left.sourceCell < right.sourceCell ;
    }

    bool withinTolerance(double actual, double expected,
                         double relativeTolerance, double& error) {
      error = std::abs(actual-expected) ;
      const double scale = std::max(1.0,std::abs(expected)) ;
      return error <= relativeTolerance*scale ;
    }

    bool withinTolerance(const vector3d<double>& actual,
                         const vector3d<double>& expected,
                         double relativeTolerance, double& error) {
      error = norm(actual-expected) ;
      const double scale = std::max(1.0,norm(expected)) ;
      return error <= relativeTolerance*scale ;
    }

    int geometryIndex(const std::vector<AMRCellGeometry>& geometry,
                      int cell) {
      AMRCellGeometry key ;
      key.cell = cell ;
      std::vector<AMRCellGeometry>::const_iterator location =
        std::lower_bound(geometry.begin(),geometry.end(),key,geometryOrder) ;
      if(location == geometry.end() || location->cell != cell)
        return -1 ;
      return int(location-geometry.begin()) ;
    }
  }

  CPTR<AMRRemapPlan> AMRRemapPlan::
  createCellPlan(const std::vector<AMRCellGeometry>& sourceGeometry,
                 const std::vector<AMRCellGeometry>& targetGeometry,
                 const std::vector<AMRCellContribution>& contributions,
                 AMRRemapReport& report,
                 double relativeTolerance) {
    return createCellPlanImpl(sourceGeometry,0,targetGeometry,contributions,
                              true,report,relativeTolerance) ;
  }

  CPTR<AMRRemapPlan> AMRRemapPlan::
  createCellPlanPartition(
    const std::vector<AMRCellGeometry>& sourceGeometry,
    const std::vector<size_t>& sourceTargetCounts,
    const std::vector<AMRCellGeometry>& targetGeometry,
    const std::vector<AMRCellContribution>& contributions,
    AMRRemapReport& report,
    double relativeTolerance) {
    return createCellPlanImpl(sourceGeometry,&sourceTargetCounts,
                              targetGeometry,contributions,false,report,
                              relativeTolerance) ;
  }

  CPTR<AMRRemapPlan> AMRRemapPlan::
  createCellPlanImpl(
    const std::vector<AMRCellGeometry>& sourceGeometry,
    const std::vector<size_t>* sourceTargetCounts,
    const std::vector<AMRCellGeometry>& targetGeometry,
    const std::vector<AMRCellContribution>& contributions,
    bool checkSourceCoverage,
    AMRRemapReport& report,
    double relativeTolerance) {
    report = AMRRemapReport() ;
    report.sourceCells = sourceGeometry.size() ;
    report.targetCells = targetGeometry.size() ;
    report.contributions = contributions.size() ;
    report.sourceCoverageChecked = checkSourceCoverage ;

    if(relativeTolerance < 0.0 || !std::isfinite(relativeTolerance)) {
      report.invalidGeometry++ ;
      return CPTR<AMRRemapPlan>() ;
    }
    if(sourceTargetCounts != 0 &&
       sourceTargetCounts->size() != sourceGeometry.size()) {
      report.invalidGeometry++ ;
      return CPTR<AMRRemapPlan>() ;
    }

    CPTR<AMRRemapPlan> plan = new AMRRemapPlan ;
    plan->sourceGeometry_ = sourceGeometry ;
    plan->targetGeometry_ = targetGeometry ;
    plan->contributions_ = contributions ;
    std::sort(plan->sourceGeometry_.begin(),plan->sourceGeometry_.end(),
              geometryOrder) ;
    std::sort(plan->targetGeometry_.begin(),plan->targetGeometry_.end(),
              geometryOrder) ;
    std::sort(plan->contributions_.begin(),plan->contributions_.end(),
              contributionOrder) ;

    for(size_t i=0;i<plan->sourceGeometry_.size();++i) {
      const AMRCellGeometry& geometry = plan->sourceGeometry_[i] ;
      if(!std::isfinite(geometry.volume) || geometry.volume <= 0.0 ||
         !finiteVector(geometry.centroid) ||
         !finiteVector(geometry.reconstructionPoint))
        report.invalidGeometry++ ;
      if(i != 0 && plan->sourceGeometry_[i-1].cell == geometry.cell)
        report.invalidGeometry++ ;
    }
    for(size_t i=0;i<plan->targetGeometry_.size();++i) {
      const AMRCellGeometry& geometry = plan->targetGeometry_[i] ;
      if(!std::isfinite(geometry.volume) || geometry.volume <= 0.0 ||
         !finiteVector(geometry.centroid) ||
         !finiteVector(geometry.reconstructionPoint))
        report.invalidGeometry++ ;
      if(i != 0 && plan->targetGeometry_[i-1].cell == geometry.cell)
        report.invalidGeometry++ ;
    }

    plan->sourceDegrees_.assign(plan->sourceGeometry_.size(),0) ;
    plan->targetDegrees_.assign(plan->targetGeometry_.size(),0) ;
    std::vector<double> sourceCoverage(plan->sourceGeometry_.size(),0.0) ;
    std::vector<double> targetCoverage(plan->targetGeometry_.size(),0.0) ;
    std::vector<vector3d<double> > sourceMoments(
      plan->sourceGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    std::vector<vector3d<double> > targetMoments(
      plan->targetGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    std::set<std::pair<int,int> > uniqueContributions ;

    for(size_t i=0;i<plan->contributions_.size();++i) {
      const AMRCellContribution& contribution = plan->contributions_[i] ;
      const int source = geometryIndex(plan->sourceGeometry_,
                                       contribution.sourceCell) ;
      const int target = geometryIndex(plan->targetGeometry_,
                                       contribution.targetCell) ;
      if(source < 0 || target < 0 ||
         !std::isfinite(contribution.overlapVolume) ||
         contribution.overlapVolume <= 0.0 ||
         !finiteVector(contribution.overlapCentroid)) {
        report.invalidGeometry++ ;
        continue ;
      }
      if(!uniqueContributions.insert(
           std::make_pair(contribution.sourceCell,
                          contribution.targetCell)).second) {
        report.duplicateContributions++ ;
        continue ;
      }
      plan->sourceDegrees_[source]++ ;
      plan->targetDegrees_[target]++ ;
      sourceCoverage[source] += contribution.overlapVolume ;
      targetCoverage[target] += contribution.overlapVolume ;
      sourceMoments[source] +=
        contribution.overlapVolume*contribution.overlapCentroid ;
      targetMoments[target] +=
        contribution.overlapVolume*contribution.overlapCentroid ;
    }

    if(sourceTargetCounts != 0) {
      // Geometry was sorted above, so associate supplied counts by cell id.
      std::map<int,size_t> targetCountBySource ;
      for(size_t source=0;source<sourceGeometry.size();++source)
        targetCountBySource[sourceGeometry[source].cell] =
          (*sourceTargetCounts)[source] ;
      for(size_t source=0;source<plan->sourceGeometry_.size();++source) {
        const std::map<int,size_t>::const_iterator count =
          targetCountBySource.find(plan->sourceGeometry_[source].cell) ;
        if(count == targetCountBySource.end() || count->second == 0)
          report.invalidGeometry++ ;
        else
          plan->sourceDegrees_[source] = count->second ;
      }
    }

    for(size_t source=0;source<plan->sourceGeometry_.size();++source) {
      double error = 0.0 ;
      const bool missing = plan->sourceDegrees_[source] == 0 ||
        (checkSourceCoverage &&
         !withinTolerance(sourceCoverage[source],
                          plan->sourceGeometry_[source].volume,
                          relativeTolerance,error)) ;
      if(missing)
        report.missingSourceCells++ ;
      report.maximumSourceVolumeError =
        std::max(report.maximumSourceVolumeError,error) ;
      if(checkSourceCoverage && sourceCoverage[source] > 0.0) {
        const vector3d<double> centroid =
          sourceMoments[source]/sourceCoverage[source] ;
        if(!withinTolerance(
             centroid,plan->sourceGeometry_[source].reconstructionPoint,
                            relativeTolerance,error))
          report.inconsistentSourceMoments++ ;
        report.maximumSourceCentroidError =
          std::max(report.maximumSourceCentroidError,error) ;
      }
    }
    for(size_t target=0;target<plan->targetGeometry_.size();++target) {
      double error = 0.0 ;
      const bool missing = plan->targetDegrees_[target] == 0 ||
        !withinTolerance(targetCoverage[target],
                         plan->targetGeometry_[target].volume,
                         relativeTolerance,error) ;
      if(missing)
        report.missingTargetCells++ ;
      report.maximumTargetVolumeError =
        std::max(report.maximumTargetVolumeError,error) ;
      if(targetCoverage[target] > 0.0) {
        const vector3d<double> centroid =
          targetMoments[target]/targetCoverage[target] ;
        if(!withinTolerance(centroid,plan->targetGeometry_[target].centroid,
                            relativeTolerance,error))
          report.inconsistentTargetMoments++ ;
        report.maximumTargetCentroidError =
          std::max(report.maximumTargetCentroidError,error) ;
      }
    }

    for(size_t i=0;i<plan->contributions_.size();++i) {
      const int source = geometryIndex(plan->sourceGeometry_,
                                       plan->contributions_[i].sourceCell) ;
      const int target = geometryIndex(plan->targetGeometry_,
                                       plan->contributions_[i].targetCell) ;
      if(source >= 0 && target >= 0 &&
         plan->sourceDegrees_[source] > 1 &&
         plan->targetDegrees_[target] > 1)
        report.unsupportedRelations++ ;
    }

    plan->targetOffsets_.assign(plan->targetGeometry_.size()+1,0) ;
    size_t contribution = 0 ;
    for(size_t target=0;target<plan->targetGeometry_.size();++target) {
      plan->targetOffsets_[target] = contribution ;
      while(contribution < plan->contributions_.size() &&
            plan->contributions_[contribution].targetCell ==
            plan->targetGeometry_[target].cell)
        ++contribution ;
    }
    plan->targetOffsets_[plan->targetGeometry_.size()] = contribution ;

    report.valid = report.invalidGeometry == 0 &&
      report.duplicateContributions == 0 &&
      report.missingSourceCells == 0 &&
      report.missingTargetCells == 0 &&
      report.inconsistentSourceMoments == 0 &&
      report.inconsistentTargetMoments == 0 &&
      report.unsupportedRelations == 0 ;
    if(!report.valid)
      return CPTR<AMRRemapPlan>() ;
    return plan ;
  }

  bool AMRRemapPlan::cellContributions(int targetCell,
                                       size_t& begin, size_t& end) const {
    const int target = geometryIndex(targetGeometry_,targetCell) ;
    if(target < 0)
      return false ;
    begin = targetOffsets_[target] ;
    end = targetOffsets_[target+1] ;
    return true ;
  }

  void AMRRemapPlan::targetCells(int sourceCell,
                                 std::vector<int>& targets) const {
    targets.clear() ;
    for(size_t i=0;i<contributions_.size();++i)
      if(contributions_[i].sourceCell == sourceCell)
        targets.push_back(contributions_[i].targetCell) ;
  }

  bool AMRRemapPlan::cellTransition(
    int targetCell, amr_cell_transition::value& transition) const {
    const int target = geometryIndex(targetGeometry_,targetCell) ;
    if(target < 0 || targetDegrees_[target] == 0)
      return false ;
    if(targetDegrees_[target] > 1) {
      transition = amr_cell_transition::derefined ;
      return true ;
    }
    const size_t contribution = targetOffsets_[target] ;
    const int source = geometryIndex(sourceGeometry_,
                                     contributions_[contribution].sourceCell) ;
    if(source < 0)
      return false ;
    transition = sourceDegrees_[source] > 1 ?
      amr_cell_transition::refined : amr_cell_transition::retained ;
    return true ;
  }

  bool AMRRemapPlan::
  remapCellAverages(const std::vector<double>& sourceValues,
                    std::vector<double>& targetValues) const {
    std::vector<vector3d<double> > gradients(
      sourceGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    return remapCellAverages(sourceValues,gradients,targetValues) ;
  }

  bool AMRRemapPlan::
  remapCellAverages(
    const std::vector<double>& sourceValues,
    const std::vector<vector3d<double> >& sourceGradients,
    std::vector<double>& targetValues) const {
    if(sourceValues.size() != sourceGeometry_.size() ||
       sourceGradients.size() != sourceGeometry_.size())
      return false ;
    targetValues.assign(targetGeometry_.size(),0.0) ;
    for(size_t target=0;target<targetGeometry_.size();++target) {
      double integral = 0.0 ;
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const AMRCellContribution& contribution = contributions_[entry] ;
        const int source = geometryIndex(sourceGeometry_,
                                         contribution.sourceCell) ;
        if(source < 0)
          return false ;
        const vector3d<double> displacement =
          contribution.overlapCentroid-
          sourceGeometry_[source].reconstructionPoint ;
        const double reconstructed = sourceValues[source] +
          dot(sourceGradients[source],displacement) ;
        integral += contribution.overlapVolume*reconstructed ;
      }
      targetValues[target] = integral/targetGeometry_[target].volume ;
    }
    return true ;
  }

  bool AMRRemapPlan::
  remapCellIntegrals(const std::vector<double>& sourceIntegrals,
                     std::vector<double>& targetIntegrals) const {
    if(sourceIntegrals.size() != sourceGeometry_.size())
      return false ;
    targetIntegrals.assign(targetGeometry_.size(),0.0) ;
    for(size_t target=0;target<targetGeometry_.size();++target) {
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const AMRCellContribution& contribution = contributions_[entry] ;
        const int source = geometryIndex(sourceGeometry_,
                                         contribution.sourceCell) ;
        if(source < 0)
          return false ;
        targetIntegrals[target] += sourceIntegrals[source]*
          contribution.overlapVolume/sourceGeometry_[source].volume ;
      }
    }
    return true ;
  }
}

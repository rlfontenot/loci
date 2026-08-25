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

#include <FVMAdapt/mesh_transfer.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace Loci {

  AMRFaceRemapReport::AMRFaceRemapReport()
    : valid(false), sourceFaces(0), targetFaces(0), contributions(0),
      createdFaces(0), removedFaces(0), invalidGeometry(0),
      duplicateContributions(0), missingSourceFaces(0),
      missingTargetFaces(0), inconsistentSourceMoments(0),
      inconsistentTargetMoments(0),
      maximumSourceAreaError(0.0), maximumTargetAreaError(0.0),
      maximumSourceCentroidError(0.0), maximumTargetCentroidError(0.0) {}

  AMRNodeRemapReport::AMRNodeRemapReport()
    : valid(false), sourceNodes(0), targetNodes(0), contributions(0),
      invalidGeometry(0), duplicateContributions(0), missingTargetNodes(0),
      invalidOrigins(0), inconsistentWeights(0), inconsistentPositions(0),
      maximumWeightError(0.0), maximumPositionError(0.0) {}

  namespace {
    bool finiteVector(const vector3d<double>& value) {
      return std::isfinite(value.x) && std::isfinite(value.y) &&
        std::isfinite(value.z) ;
    }

    bool faceGeometryOrder(const AMRFaceGeometry& left,
                           const AMRFaceGeometry& right) {
      return left.face < right.face ;
    }

    bool faceContributionOrder(const AMRFaceContribution& left,
                               const AMRFaceContribution& right) {
      if(left.targetFace != right.targetFace)
        return left.targetFace < right.targetFace ;
      return left.sourceFace < right.sourceFace ;
    }

    bool createdFaceOrder(const AMRCreatedFace& left,
                          const AMRCreatedFace& right) {
      return left.targetFace < right.targetFace ;
    }

    bool removedFaceOrder(const AMRRemovedFace& left,
                          const AMRRemovedFace& right) {
      return left.sourceFace < right.sourceFace ;
    }

    int faceGeometryIndex(const std::vector<AMRFaceGeometry>& geometry,
                          int face) {
      AMRFaceGeometry key ;
      key.face = face ;
      std::vector<AMRFaceGeometry>::const_iterator location =
        std::lower_bound(geometry.begin(),geometry.end(),key,
                         faceGeometryOrder) ;
      if(location == geometry.end() || location->face != face)
        return -1 ;
      return int(location-geometry.begin()) ;
    }

    bool nodeGeometryOrder(const AMRNodeGeometry& left,
                           const AMRNodeGeometry& right) {
      return left.node < right.node ;
    }

    bool nodeContributionOrder(const AMRNodeContribution& left,
                               const AMRNodeContribution& right) {
      if(left.targetNode != right.targetNode)
        return left.targetNode < right.targetNode ;
      return left.sourceNode < right.sourceNode ;
    }

    bool nodeOriginOrder(const AMRNodeOrigin& left,
                         const AMRNodeOrigin& right) {
      return left.targetNode < right.targetNode ;
    }

    int nodeGeometryIndex(const std::vector<AMRNodeGeometry>& geometry,
                          int node) {
      AMRNodeGeometry key ;
      key.node = node ;
      std::vector<AMRNodeGeometry>::const_iterator location =
        std::lower_bound(geometry.begin(),geometry.end(),key,
                         nodeGeometryOrder) ;
      if(location == geometry.end() || location->node != node)
        return -1 ;
      return int(location-geometry.begin()) ;
    }

    bool withinTolerance(double actual, double expected,
                         double relativeTolerance, double& error) {
      error = std::abs(actual-expected) ;
      return error <= relativeTolerance*std::max(1.0,std::abs(expected)) ;
    }

    bool withinTolerance(const vector3d<double>& actual,
                         const vector3d<double>& expected,
                         double relativeTolerance, double& error) {
      error = norm(actual-expected) ;
      return error <= relativeTolerance*std::max(1.0,norm(expected)) ;
    }
  }

  CPTR<AMRFaceRemap> AMRFaceRemap::
  create(const std::vector<AMRFaceGeometry>& sourceGeometry,
         const std::vector<AMRFaceGeometry>& targetGeometry,
         const std::vector<AMRFaceContribution>& contributions,
         const std::vector<AMRCreatedFace>& createdFaces,
         const std::vector<AMRRemovedFace>& removedFaces,
         AMRFaceRemapReport& report,
         double relativeTolerance) {
    report = AMRFaceRemapReport() ;
    report.sourceFaces = sourceGeometry.size() ;
    report.targetFaces = targetGeometry.size() ;
    report.contributions = contributions.size() ;
    report.createdFaces = createdFaces.size() ;
    report.removedFaces = removedFaces.size() ;
    if(relativeTolerance < 0.0 || !std::isfinite(relativeTolerance)) {
      report.invalidGeometry++ ;
      return CPTR<AMRFaceRemap>() ;
    }

    CPTR<AMRFaceRemap> remap = new AMRFaceRemap ;
    remap->sourceGeometry_ = sourceGeometry ;
    remap->targetGeometry_ = targetGeometry ;
    remap->contributions_ = contributions ;
    remap->createdFaces_ = createdFaces ;
    remap->removedFaces_ = removedFaces ;
    std::sort(remap->sourceGeometry_.begin(),remap->sourceGeometry_.end(),
              faceGeometryOrder) ;
    std::sort(remap->targetGeometry_.begin(),remap->targetGeometry_.end(),
              faceGeometryOrder) ;
    std::sort(remap->contributions_.begin(),remap->contributions_.end(),
              faceContributionOrder) ;
    std::sort(remap->createdFaces_.begin(),remap->createdFaces_.end(),
              createdFaceOrder) ;
    std::sort(remap->removedFaces_.begin(),remap->removedFaces_.end(),
              removedFaceOrder) ;

    for(size_t i=0;i<remap->sourceGeometry_.size();++i) {
      const AMRFaceGeometry& geometry = remap->sourceGeometry_[i] ;
      if(!std::isfinite(geometry.area) || geometry.area <= 0.0 ||
         !finiteVector(geometry.centroid))
        report.invalidGeometry++ ;
      if(i != 0 && remap->sourceGeometry_[i-1].face == geometry.face)
        report.invalidGeometry++ ;
    }
    for(size_t i=0;i<remap->targetGeometry_.size();++i) {
      const AMRFaceGeometry& geometry = remap->targetGeometry_[i] ;
      if(!std::isfinite(geometry.area) || geometry.area <= 0.0 ||
         !finiteVector(geometry.centroid))
        report.invalidGeometry++ ;
      if(i != 0 && remap->targetGeometry_[i-1].face == geometry.face)
        report.invalidGeometry++ ;
    }

    std::vector<size_t> sourceDegrees(remap->sourceGeometry_.size(),0) ;
    std::vector<size_t> targetDegrees(remap->targetGeometry_.size(),0) ;
    std::vector<double> sourceCoverage(remap->sourceGeometry_.size(),0.0) ;
    std::vector<double> targetCoverage(remap->targetGeometry_.size(),0.0) ;
    std::vector<vector3d<double> > sourceMoments(
      remap->sourceGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    std::vector<vector3d<double> > targetMoments(
      remap->targetGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    std::set<std::pair<int,int> > uniqueContributions ;
    for(size_t i=0;i<remap->contributions_.size();++i) {
      const AMRFaceContribution& contribution = remap->contributions_[i] ;
      const int source = faceGeometryIndex(remap->sourceGeometry_,
                                           contribution.sourceFace) ;
      const int target = faceGeometryIndex(remap->targetGeometry_,
                                           contribution.targetFace) ;
      if(source < 0 || target < 0 ||
         !std::isfinite(contribution.overlapArea) ||
         contribution.overlapArea <= 0.0 ||
         !finiteVector(contribution.overlapCentroid) ||
         (contribution.orientation != 1 && contribution.orientation != -1)) {
        report.invalidGeometry++ ;
        continue ;
      }
      if(!uniqueContributions.insert(
           std::make_pair(contribution.sourceFace,
                          contribution.targetFace)).second) {
        report.duplicateContributions++ ;
        continue ;
      }
      sourceDegrees[source]++ ;
      targetDegrees[target]++ ;
      sourceCoverage[source] += contribution.overlapArea ;
      targetCoverage[target] += contribution.overlapArea ;
      sourceMoments[source] +=
        contribution.overlapArea*contribution.overlapCentroid ;
      targetMoments[target] +=
        contribution.overlapArea*contribution.overlapCentroid ;
    }

    std::set<int> createdSet ;
    for(size_t i=0;i<remap->createdFaces_.size();++i) {
      const int target = faceGeometryIndex(remap->targetGeometry_,
                                           remap->createdFaces_[i].targetFace) ;
      if(target < 0 || !createdSet.insert(
           remap->createdFaces_[i].targetFace).second ||
         (target >= 0 && targetDegrees[target] != 0))
        report.invalidGeometry++ ;
    }
    std::set<int> removedSet ;
    for(size_t i=0;i<remap->removedFaces_.size();++i) {
      const int source = faceGeometryIndex(remap->sourceGeometry_,
                                           remap->removedFaces_[i].sourceFace) ;
      if(source < 0 || !removedSet.insert(
           remap->removedFaces_[i].sourceFace).second ||
         (source >= 0 && sourceDegrees[source] != 0))
        report.invalidGeometry++ ;
    }

    for(size_t source=0;source<remap->sourceGeometry_.size();++source) {
      if(sourceDegrees[source] == 0) {
        if(removedSet.count(remap->sourceGeometry_[source].face) == 0)
          report.missingSourceFaces++ ;
        continue ;
      }
      double error = 0.0 ;
      if(!withinTolerance(sourceCoverage[source],
                          remap->sourceGeometry_[source].area,
                          relativeTolerance,error))
        report.missingSourceFaces++ ;
      report.maximumSourceAreaError =
        std::max(report.maximumSourceAreaError,error) ;
      const vector3d<double> centroid =
        sourceMoments[source]/sourceCoverage[source] ;
      if(!withinTolerance(centroid,remap->sourceGeometry_[source].centroid,
                          relativeTolerance,error))
        report.inconsistentSourceMoments++ ;
      report.maximumSourceCentroidError =
        std::max(report.maximumSourceCentroidError,error) ;
    }
    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      if(targetDegrees[target] == 0) {
        if(createdSet.count(remap->targetGeometry_[target].face) == 0)
          report.missingTargetFaces++ ;
        continue ;
      }
      double error = 0.0 ;
      if(!withinTolerance(targetCoverage[target],
                          remap->targetGeometry_[target].area,
                          relativeTolerance,error))
        report.missingTargetFaces++ ;
      report.maximumTargetAreaError =
        std::max(report.maximumTargetAreaError,error) ;
      const vector3d<double> centroid =
        targetMoments[target]/targetCoverage[target] ;
      if(!withinTolerance(centroid,remap->targetGeometry_[target].centroid,
                          relativeTolerance,error))
        report.inconsistentTargetMoments++ ;
      report.maximumTargetCentroidError =
        std::max(report.maximumTargetCentroidError,error) ;
    }

    remap->targetOffsets_.assign(remap->targetGeometry_.size()+1,0) ;
    size_t contribution = 0 ;
    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      remap->targetOffsets_[target] = contribution ;
      while(contribution < remap->contributions_.size() &&
            remap->contributions_[contribution].targetFace ==
            remap->targetGeometry_[target].face)
        ++contribution ;
    }
    remap->targetOffsets_[remap->targetGeometry_.size()] = contribution ;

    report.valid = report.invalidGeometry == 0 &&
      report.duplicateContributions == 0 && report.missingSourceFaces == 0 &&
      report.missingTargetFaces == 0 &&
      report.inconsistentSourceMoments == 0 &&
      report.inconsistentTargetMoments == 0 ;
    if(!report.valid)
      return CPTR<AMRFaceRemap>() ;
    return remap ;
  }

  bool AMRFaceRemap::faceContributions(int targetFace,
                                       size_t& begin, size_t& end) const {
    const int target = faceGeometryIndex(targetGeometry_,targetFace) ;
    if(target < 0)
      return false ;
    begin = targetOffsets_[target] ;
    end = targetOffsets_[target+1] ;
    return true ;
  }

  bool AMRFaceRemap::isCreatedFace(int targetFace, int& sourceCell) const {
    AMRCreatedFace key ;
    key.targetFace = targetFace ;
    const std::vector<AMRCreatedFace>::const_iterator face =
      std::lower_bound(createdFaces_.begin(),createdFaces_.end(),key,
                       createdFaceOrder) ;
    if(face == createdFaces_.end() || face->targetFace != targetFace)
      return false ;
    sourceCell = face->sourceCell ;
    return true ;
  }

  bool AMRFaceRemap::isRemovedFace(int sourceFace, int& targetCell) const {
    AMRRemovedFace key ;
    key.sourceFace = sourceFace ;
    const std::vector<AMRRemovedFace>::const_iterator face =
      std::lower_bound(removedFaces_.begin(),removedFaces_.end(),key,
                       removedFaceOrder) ;
    if(face == removedFaces_.end() || face->sourceFace != sourceFace)
      return false ;
    targetCell = face->targetCell ;
    return true ;
  }

  bool AMRFaceRemap::
  remapFaceAverages(const std::vector<double>& sourceValues,
                    std::vector<double>& targetValues,
                    std::vector<unsigned char>& mapped,
                    bool orientValues) const {
    if(sourceValues.size() != sourceGeometry_.size())
      return false ;
    targetValues.assign(targetGeometry_.size(),0.0) ;
    mapped.assign(targetGeometry_.size(),0) ;
    for(size_t target=0;target<targetGeometry_.size();++target) {
      if(targetOffsets_[target] == targetOffsets_[target+1])
        continue ;
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const AMRFaceContribution& contribution = contributions_[entry] ;
        const int source = faceGeometryIndex(sourceGeometry_,
                                             contribution.sourceFace) ;
        if(source < 0)
          return false ;
        const int orientation = orientValues ? contribution.orientation : 1 ;
        targetValues[target] += orientation*contribution.overlapArea*
          sourceValues[source]/targetGeometry_[target].area ;
      }
      mapped[target] = 1 ;
    }
    return true ;
  }

  bool AMRFaceRemap::
  remapFaceIntegrals(const std::vector<double>& sourceIntegrals,
                     std::vector<double>& targetIntegrals,
                     std::vector<unsigned char>& mapped,
                     bool orientValues) const {
    if(sourceIntegrals.size() != sourceGeometry_.size())
      return false ;
    targetIntegrals.assign(targetGeometry_.size(),0.0) ;
    mapped.assign(targetGeometry_.size(),0) ;
    for(size_t target=0;target<targetGeometry_.size();++target) {
      if(targetOffsets_[target] == targetOffsets_[target+1])
        continue ;
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const AMRFaceContribution& contribution = contributions_[entry] ;
        const int source = faceGeometryIndex(sourceGeometry_,
                                             contribution.sourceFace) ;
        if(source < 0)
          return false ;
        const int orientation = orientValues ? contribution.orientation : 1 ;
        targetIntegrals[target] += orientation*sourceIntegrals[source]*
          contribution.overlapArea/sourceGeometry_[source].area ;
      }
      mapped[target] = 1 ;
    }
    return true ;
  }

  CPTR<AMRNodeRemap> AMRNodeRemap::
  create(const std::vector<AMRNodeGeometry>& sourceGeometry,
         const std::vector<AMRNodeGeometry>& targetGeometry,
         const std::vector<AMRNodeContribution>& contributions,
         const std::vector<AMRNodeOrigin>& origins,
         AMRNodeRemapReport& report,
         double relativeTolerance) {
    report = AMRNodeRemapReport() ;
    report.sourceNodes = sourceGeometry.size() ;
    report.targetNodes = targetGeometry.size() ;
    report.contributions = contributions.size() ;
    if(relativeTolerance < 0.0 || !std::isfinite(relativeTolerance)) {
      report.invalidGeometry++ ;
      return CPTR<AMRNodeRemap>() ;
    }

    CPTR<AMRNodeRemap> remap = new AMRNodeRemap ;
    remap->sourceGeometry_ = sourceGeometry ;
    remap->targetGeometry_ = targetGeometry ;
    remap->contributions_ = contributions ;
    remap->origins_ = origins ;
    std::sort(remap->sourceGeometry_.begin(),remap->sourceGeometry_.end(),
              nodeGeometryOrder) ;
    std::sort(remap->targetGeometry_.begin(),remap->targetGeometry_.end(),
              nodeGeometryOrder) ;
    std::sort(remap->contributions_.begin(),remap->contributions_.end(),
              nodeContributionOrder) ;
    std::sort(remap->origins_.begin(),remap->origins_.end(),nodeOriginOrder) ;

    for(size_t i=0;i<remap->sourceGeometry_.size();++i) {
      if(!finiteVector(remap->sourceGeometry_[i].position))
        report.invalidGeometry++ ;
      if(i != 0 && remap->sourceGeometry_[i-1].node ==
         remap->sourceGeometry_[i].node)
        report.invalidGeometry++ ;
    }
    for(size_t i=0;i<remap->targetGeometry_.size();++i) {
      if(!finiteVector(remap->targetGeometry_[i].position))
        report.invalidGeometry++ ;
      if(i != 0 && remap->targetGeometry_[i-1].node ==
         remap->targetGeometry_[i].node)
        report.invalidGeometry++ ;
    }

    std::vector<double> weightSum(remap->targetGeometry_.size(),0.0) ;
    std::vector<vector3d<double> > positions(
      remap->targetGeometry_.size(),vector3d<double>(0.0,0.0,0.0)) ;
    std::vector<size_t> targetDegrees(remap->targetGeometry_.size(),0) ;
    std::set<std::pair<int,int> > uniqueContributions ;
    for(size_t i=0;i<remap->contributions_.size();++i) {
      const AMRNodeContribution& contribution = remap->contributions_[i] ;
      const int source = nodeGeometryIndex(remap->sourceGeometry_,
                                           contribution.sourceNode) ;
      const int target = nodeGeometryIndex(remap->targetGeometry_,
                                           contribution.targetNode) ;
      if(source < 0 || target < 0 || !std::isfinite(contribution.weight) ||
         contribution.weight < -relativeTolerance) {
        report.invalidGeometry++ ;
        continue ;
      }
      if(!uniqueContributions.insert(
           std::make_pair(contribution.sourceNode,
                          contribution.targetNode)).second) {
        report.duplicateContributions++ ;
        continue ;
      }
      targetDegrees[target]++ ;
      weightSum[target] += contribution.weight ;
      positions[target] += contribution.weight*
        remap->sourceGeometry_[source].position ;
    }

    std::map<int,AMRNodeOrigin> originByTarget ;
    for(size_t i=0;i<remap->origins_.size();++i) {
      if(nodeGeometryIndex(remap->targetGeometry_,
                           remap->origins_[i].targetNode) < 0 ||
         remap->origins_[i].kind < amr_node_origin::retained ||
         remap->origins_[i].kind > amr_node_origin::cell ||
         !originByTarget.insert(std::make_pair(
           remap->origins_[i].targetNode,remap->origins_[i])).second)
        report.invalidOrigins++ ;
    }

    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      const int targetNode = remap->targetGeometry_[target].node ;
      if(targetDegrees[target] == 0)
        report.missingTargetNodes++ ;
      const std::map<int,AMRNodeOrigin>::const_iterator origin =
        originByTarget.find(targetNode) ;
      if(origin == originByTarget.end()) {
        report.invalidOrigins++ ;
      } else if(origin->second.kind == amr_node_origin::retained) {
        if(targetDegrees[target] != 1) {
          report.invalidOrigins++ ;
        } else {
          int retainedSource = 0 ;
          bool foundRetainedSource = false ;
          for(size_t entry=0;entry<remap->contributions_.size();++entry)
            if(remap->contributions_[entry].targetNode == targetNode) {
              retainedSource = remap->contributions_[entry].sourceNode ;
              foundRetainedSource = true ;
              break ;
            }
          if(!foundRetainedSource ||
             retainedSource != origin->second.sourceEntity)
            report.invalidOrigins++ ;
        }
      }
      double error = 0.0 ;
      if(!withinTolerance(weightSum[target],1.0,
                          relativeTolerance,error))
        report.inconsistentWeights++ ;
      report.maximumWeightError = std::max(report.maximumWeightError,error) ;
      if(weightSum[target] > 0.0 && !withinTolerance(
           positions[target]/weightSum[target],
           remap->targetGeometry_[target].position,
           relativeTolerance,error))
        report.inconsistentPositions++ ;
      report.maximumPositionError =
        std::max(report.maximumPositionError,error) ;
    }

    remap->targetOffsets_.assign(remap->targetGeometry_.size()+1,0) ;
    size_t contribution = 0 ;
    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      remap->targetOffsets_[target] = contribution ;
      while(contribution < remap->contributions_.size() &&
            remap->contributions_[contribution].targetNode ==
            remap->targetGeometry_[target].node)
        ++contribution ;
    }
    remap->targetOffsets_[remap->targetGeometry_.size()] = contribution ;

    report.valid = report.invalidGeometry == 0 &&
      report.duplicateContributions == 0 && report.missingTargetNodes == 0 &&
      report.invalidOrigins == 0 && report.inconsistentWeights == 0 &&
      report.inconsistentPositions == 0 ;
    if(!report.valid)
      return CPTR<AMRNodeRemap>() ;
    return remap ;
  }

  bool AMRNodeRemap::nodeContributions(int targetNode,
                                       size_t& begin, size_t& end) const {
    const int target = nodeGeometryIndex(targetGeometry_,targetNode) ;
    if(target < 0)
      return false ;
    begin = targetOffsets_[target] ;
    end = targetOffsets_[target+1] ;
    return true ;
  }

  bool AMRNodeRemap::nodeOrigin(int targetNode, AMRNodeOrigin& origin) const {
    AMRNodeOrigin key ;
    key.targetNode = targetNode ;
    const std::vector<AMRNodeOrigin>::const_iterator location =
      std::lower_bound(origins_.begin(),origins_.end(),key,nodeOriginOrder) ;
    if(location == origins_.end() || location->targetNode != targetNode)
      return false ;
    origin = *location ;
    return true ;
  }

  bool AMRNodeRemap::
  interpolateNodeData(const std::vector<double>& sourceValues,
                      std::vector<double>& targetValues) const {
    if(sourceValues.size() != sourceGeometry_.size())
      return false ;
    targetValues.assign(targetGeometry_.size(),0.0) ;
    for(size_t target=0;target<targetGeometry_.size();++target)
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const int source = nodeGeometryIndex(sourceGeometry_,
                                             contributions_[entry].sourceNode) ;
        if(source < 0)
          return false ;
        targetValues[target] +=
          contributions_[entry].weight*sourceValues[source] ;
      }
    return true ;
  }

  bool AMRNodeRemap::
  interpolateNodeData(const std::vector<vector3d<double> >& sourceValues,
                      std::vector<vector3d<double> >& targetValues) const {
    if(sourceValues.size() != sourceGeometry_.size())
      return false ;
    targetValues.assign(targetGeometry_.size(),
                        vector3d<double>(0.0,0.0,0.0)) ;
    for(size_t target=0;target<targetGeometry_.size();++target)
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const int source = nodeGeometryIndex(sourceGeometry_,
                                             contributions_[entry].sourceNode) ;
        if(source < 0)
          return false ;
        targetValues[target] +=
          contributions_[entry].weight*sourceValues[source] ;
      }
    return true ;
  }
}

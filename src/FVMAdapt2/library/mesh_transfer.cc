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

#include <FVMAdapt2/mesh_transfer.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace Loci {

  FaceRemapReport::FaceRemapReport()
      : valid(false), sourceFaces(0), targetFaces(0), contributions(0),
        createdFaces(0), removedFaces(0), invalidGeometry(0),
        duplicateContributions(0), missingSourceFaces(0), missingTargetFaces(0),
        inconsistentSourceMoments(0), inconsistentTargetMoments(0),
        maximumSourceAreaError(0.0), maximumTargetAreaError(0.0),
        maximumSourceCentroidError(0.0), maximumTargetCentroidError(0.0) {}

  NodeRemapReport::NodeRemapReport()
      : valid(false), sourceNodes(0), targetNodes(0), contributions(0),
        invalidGeometry(0), duplicateContributions(0), missingTargetNodes(0),
        invalidOrigins(0), inconsistentWeights(0), inconsistentPositions(0),
        maximumWeightError(0.0), maximumPositionError(0.0) {}

  namespace {
    bool finiteVector(const vector3d<double>& value) {
      return std::isfinite(value.x) && std::isfinite(value.y) &&
        std::isfinite(value.z) ;
    }

    bool faceGeometryOrder(
          const FaceGeometry& left, const FaceGeometry& right) {
      return left.face < right.face ;
    }

    bool faceOverlapOrder(const FaceOverlap& left, const FaceOverlap& right) {
      if (left.target != right.target)
        return left.target < right.target ;
      return left.source < right.source ;
    }

    bool createdFaceOrder(const CreatedFace& left, const CreatedFace& right) {
      return left.targetFace < right.targetFace ;
    }

    bool removedFaceOrder(const RemovedFace& left, const RemovedFace& right) {
      return left.sourceFace < right.sourceFace ;
    }

    int faceGeometryIndex(
          const std::vector<FaceGeometry>& geometry, FaceId face) {
      FaceGeometry key ;
      key.face = face ;
      std::vector<FaceGeometry>::const_iterator location = std::lower_bound(
            geometry.begin(), geometry.end(), key, faceGeometryOrder) ;
      if(location == geometry.end() || location->face != face)
        return -1 ;
      return int(location-geometry.begin()) ;
    }

    bool nodeGeometryOrder(
          const NodeGeometry& left, const NodeGeometry& right) {
      return left.node < right.node ;
    }

    bool nodeContributionOrder(
          const NodeContribution& left, const NodeContribution& right) {
      if(left.targetNode != right.targetNode)
        return left.targetNode < right.targetNode ;
      return left.sourceNode < right.sourceNode ;
    }

    bool nodeOriginOrder(const NodeOrigin& left, const NodeOrigin& right) {
      return left.targetNode < right.targetNode ;
    }

    int nodeGeometryIndex(
          const std::vector<NodeGeometry>& geometry, NodeId node) {
      NodeGeometry key ;
      key.node = node ;
      std::vector<NodeGeometry>::const_iterator location = std::lower_bound(
            geometry.begin(), geometry.end(), key, nodeGeometryOrder) ;
      if(location == geometry.end() || location->node != node)
        return -1 ;
      return int(location-geometry.begin()) ;
    }

    bool withinTolerance(double actual, double expected,
                         double relativeTolerance, double& error) {
      error = std::abs(actual-expected) ;
      const double scale = std::max(std::abs(actual), std::abs(expected)) ;
      return scale == 0.0 ? error == 0.0 : error <= relativeTolerance * scale ;
    }

    bool withinScale(double error, double scale, double relativeTolerance) {
      return scale == 0.0 ? error == 0.0 : error <= relativeTolerance * scale ;
    }

    bool withinTolerance(const vector3d<double>& actual,
                         const vector3d<double>& expected,
                         double relativeTolerance, double& error) {
      error = norm(actual-expected) ;
      return error <= relativeTolerance*std::max(1.0,norm(expected)) ;
    }
  }

  CPTR<FaceRemap> FaceRemap::create(
        const std::vector<FaceGeometry>& sourceGeometry,
        const std::vector<FaceGeometry>& targetGeometry,
        const std::vector<FaceOverlap>& contributions,
        const std::vector<CreatedFace>& createdFaces,
        const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
        double relativeTolerance) {
    return createImpl(sourceGeometry, targetGeometry, contributions,
          createdFaces, removedFaces, report, relativeTolerance, true) ;
  }

  CPTR<FaceRemap> FaceRemap::createTargetOwned(
        const std::vector<FaceGeometry>& sourceGeometry,
        const std::vector<FaceGeometry>& targetGeometry,
        const std::vector<FaceOverlap>& contributions,
        const std::vector<CreatedFace>& createdFaces,
        const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
        double relativeTolerance) {
    return createImpl(sourceGeometry, targetGeometry, contributions,
          createdFaces, removedFaces, report, relativeTolerance, false) ;
  }

  CPTR<FaceRemap> FaceRemap::createImpl(
        const std::vector<FaceGeometry>& sourceGeometry,
        const std::vector<FaceGeometry>& targetGeometry,
        const std::vector<FaceOverlap>& contributions,
        const std::vector<CreatedFace>& createdFaces,
        const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
        double relativeTolerance, bool validateSourceCoverage) {
    report = FaceRemapReport() ;
    report.sourceFaces = sourceGeometry.size() ;
    report.targetFaces = targetGeometry.size() ;
    report.contributions = contributions.size() ;
    report.createdFaces = createdFaces.size() ;
    report.removedFaces = removedFaces.size() ;
    if(relativeTolerance < 0.0 || !std::isfinite(relativeTolerance)) {
      report.invalidGeometry++ ;
      return CPTR<FaceRemap>() ;
    }

    CPTR<FaceRemap> remap = new FaceRemap ;
    remap->sourceGeometry_ = sourceGeometry ;
    remap->targetGeometry_ = targetGeometry ;
    remap->contributions_ = contributions ;
    remap->createdFaces_ = createdFaces ;
    remap->removedFaces_ = removedFaces ;
    std::sort(remap->sourceGeometry_.begin(),remap->sourceGeometry_.end(),
              faceGeometryOrder) ;
    std::sort(remap->targetGeometry_.begin(),remap->targetGeometry_.end(),
              faceGeometryOrder) ;
    std::sort(remap->contributions_.begin(), remap->contributions_.end(),
          faceOverlapOrder) ;
    std::sort(remap->createdFaces_.begin(),remap->createdFaces_.end(),
              createdFaceOrder) ;
    std::sort(remap->removedFaces_.begin(),remap->removedFaces_.end(),
              removedFaceOrder) ;

    for(size_t i=0;i<remap->sourceGeometry_.size();++i) {
      const FaceGeometry& geometry = remap->sourceGeometry_[i] ;
      if(!std::isfinite(geometry.area) || geometry.area <= 0.0 ||
         !finiteVector(geometry.centroid))
        report.invalidGeometry++ ;
      if(i != 0 && remap->sourceGeometry_[i-1].face == geometry.face)
        report.invalidGeometry++ ;
    }
    for(size_t i=0;i<remap->targetGeometry_.size();++i) {
      const FaceGeometry& geometry = remap->targetGeometry_[i] ;
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
    std::vector<double> sourceScales(remap->sourceGeometry_.size(), 0.0) ;
    std::vector<double> targetScales(remap->targetGeometry_.size(), 0.0) ;
    for (size_t source = 0; source < remap->sourceGeometry_.size(); ++source)
      sourceScales[source] =
            std::sqrt(std::max(0.0, remap->sourceGeometry_[source].area)) ;
    for (size_t target = 0; target < remap->targetGeometry_.size(); ++target)
      targetScales[target] =
            std::sqrt(std::max(0.0, remap->targetGeometry_[target].area)) ;
    std::set<std::pair<FaceId, FaceId>> uniqueContributions ;
    for(size_t i=0;i<remap->contributions_.size();++i) {
      const FaceOverlap& contribution = remap->contributions_[i] ;
      const int source =
            faceGeometryIndex(remap->sourceGeometry_, contribution.source) ;
      const int target =
            faceGeometryIndex(remap->targetGeometry_, contribution.target) ;
      if (source < 0 || target < 0 || !std::isfinite(contribution.area) ||
            contribution.area <= 0.0 || !finiteVector(contribution.centroid) ||
            (contribution.orientation != 1 && contribution.orientation != -1)) {
        report.invalidGeometry++ ;
        continue ;
      }
      if (!uniqueContributions
                  .insert(std::make_pair(
                        contribution.source, contribution.target))
                  .second) {
        report.duplicateContributions++ ;
        continue ;
      }
      sourceDegrees[source]++ ;
      targetDegrees[target]++ ;
      sourceCoverage[source] += contribution.area ;
      targetCoverage[target] += contribution.area ;
      sourceMoments[source] +=
            contribution.area *
            (contribution.centroid - remap->sourceGeometry_[source].centroid) ;
      targetMoments[target] +=
            contribution.area *
            (contribution.centroid - remap->targetGeometry_[target].centroid) ;
      sourceScales[source] = std::max(sourceScales[source],
            norm(contribution.centroid -
                  remap->sourceGeometry_[source].centroid)) ;
      targetScales[target] = std::max(targetScales[target],
            norm(contribution.centroid -
                  remap->targetGeometry_[target].centroid)) ;
    }

    std::set<FaceId> createdSet ;
    for(size_t i=0;i<remap->createdFaces_.size();++i) {
      const int target = faceGeometryIndex(remap->targetGeometry_,
                                           remap->createdFaces_[i].targetFace) ;
      if(target < 0 || !createdSet.insert(
           remap->createdFaces_[i].targetFace).second ||
         (target >= 0 && targetDegrees[target] != 0))
        report.invalidGeometry++ ;
    }
    std::set<FaceId> removedSet ;
    for(size_t i=0;i<remap->removedFaces_.size();++i) {
      const int source = faceGeometryIndex(remap->sourceGeometry_,
                                           remap->removedFaces_[i].sourceFace) ;
      if(source < 0 || !removedSet.insert(
           remap->removedFaces_[i].sourceFace).second ||
         (source >= 0 && sourceDegrees[source] != 0))
        report.invalidGeometry++ ;
    }

    if (validateSourceCoverage) {
      for (size_t source = 0; source < remap->sourceGeometry_.size();
            ++source) {
        if (sourceDegrees[source] == 0) {
          if (removedSet.count(remap->sourceGeometry_[source].face) == 0)
            report.missingSourceFaces++ ;
          continue ;
        }
        double error = 0.0 ;
        if (!withinTolerance(sourceCoverage[source],
                  remap->sourceGeometry_[source].area, relativeTolerance,
                  error))
          report.missingSourceFaces++ ;
        report.maximumSourceAreaError =
              std::max(report.maximumSourceAreaError, error) ;
        error = norm(sourceMoments[source] / sourceCoverage[source]) ;
        if (!withinScale(error, sourceScales[source], relativeTolerance))
          report.inconsistentSourceMoments++ ;
        report.maximumSourceCentroidError =
              std::max(report.maximumSourceCentroidError, error) ;
      }
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
      error = norm(targetMoments[target] / targetCoverage[target]) ;
      if (!withinScale(error, targetScales[target], relativeTolerance))
        report.inconsistentTargetMoments++ ;
      report.maximumTargetCentroidError =
        std::max(report.maximumTargetCentroidError,error) ;
    }

    remap->targetOffsets_.assign(remap->targetGeometry_.size()+1,0) ;
    size_t contribution = 0 ;
    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      remap->targetOffsets_[target] = contribution ;
      while (contribution < remap->contributions_.size() &&
             remap->contributions_[contribution].target ==
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
      return CPTR<FaceRemap>() ;

    remap->sourceOffsets_.assign(remap->sourceGeometry_.size() + 1, 0) ;
    for (size_t source = 0; source < remap->sourceGeometry_.size(); ++source)
      remap->sourceOffsets_[source + 1] =
            remap->sourceOffsets_[source] + sourceDegrees[source] ;
    remap->sourceTargets_.assign(remap->contributions_.size(), 0) ;
    std::vector<size_t> nextSourceOffset = remap->sourceOffsets_ ;
    for (size_t entry = 0; entry < remap->contributions_.size(); ++entry) {
      const FaceOverlap& contribution = remap->contributions_[entry] ;
      const int source =
            faceGeometryIndex(remap->sourceGeometry_, contribution.source) ;
      remap->sourceTargets_[nextSourceOffset[source]++] = contribution.target ;
    }
    return remap ;
  }

  bool FaceRemap::overlaps(
        FaceId targetFace, size_t& begin, size_t& end) const {
    const int target = faceGeometryIndex(targetGeometry_,targetFace) ;
    if(target < 0)
      return false ;
    begin = targetOffsets_[target] ;
    end = targetOffsets_[target+1] ;
    return true ;
  }

  void FaceRemap::targetFaces(
        FaceId sourceFace, std::vector<FaceId>& targets) const {
    targets.clear() ;
    const int source = faceGeometryIndex(sourceGeometry_, sourceFace) ;
    if (source < 0)
      return ;
    targets.insert(targets.end(),
          sourceTargets_.begin() + sourceOffsets_[source],
          sourceTargets_.begin() + sourceOffsets_[source + 1]) ;
  }

  bool FaceRemap::isCreatedFace(FaceId targetFace, CellId& sourceCell) const {
    CreatedFace key ;
    key.targetFace = targetFace ;
    const std::vector<CreatedFace>::const_iterator face = std::lower_bound(
          createdFaces_.begin(), createdFaces_.end(), key, createdFaceOrder) ;
    if(face == createdFaces_.end() || face->targetFace != targetFace)
      return false ;
    sourceCell = face->sourceCell ;
    return true ;
  }

  bool FaceRemap::isRemovedFace(FaceId sourceFace, CellId& targetCell) const {
    RemovedFace key ;
    key.sourceFace = sourceFace ;
    const std::vector<RemovedFace>::const_iterator face = std::lower_bound(
          removedFaces_.begin(), removedFaces_.end(), key, removedFaceOrder) ;
    if(face == removedFaces_.end() || face->sourceFace != sourceFace)
      return false ;
    targetCell = face->targetCell ;
    return true ;
  }

  bool FaceRemap::remapFaceAverages(const std::vector<double>& sourceValues,
        std::vector<double>& targetValues, std::vector<unsigned char>& mapped,
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
        const FaceOverlap& contribution = contributions_[entry] ;
        const int source =
              faceGeometryIndex(sourceGeometry_, contribution.source) ;
        if(source < 0)
          return false ;
        const int orientation = orientValues ? contribution.orientation : 1 ;
        targetValues[target] += orientation * contribution.area *
                                sourceValues[source] /
                                targetGeometry_[target].area ;
      }
      mapped[target] = 1 ;
    }
    return true ;
  }

  bool FaceRemap::remapFaceIntegrals(const std::vector<double>& sourceIntegrals,
        std::vector<double>& targetIntegrals,
        std::vector<unsigned char>& mapped, bool orientValues) const {
    if(sourceIntegrals.size() != sourceGeometry_.size())
      return false ;
    targetIntegrals.assign(targetGeometry_.size(),0.0) ;
    mapped.assign(targetGeometry_.size(),0) ;
    for(size_t target=0;target<targetGeometry_.size();++target) {
      if(targetOffsets_[target] == targetOffsets_[target+1])
        continue ;
      for(size_t entry=targetOffsets_[target];
          entry<targetOffsets_[target+1];++entry) {
        const FaceOverlap& contribution = contributions_[entry] ;
        const int source =
              faceGeometryIndex(sourceGeometry_, contribution.source) ;
        if(source < 0)
          return false ;
        const int orientation = orientValues ? contribution.orientation : 1 ;
        targetIntegrals[target] += orientation * sourceIntegrals[source] *
                                   contribution.area /
                                   sourceGeometry_[source].area ;
      }
      mapped[target] = 1 ;
    }
    return true ;
  }

  CPTR<NodeRemap> NodeRemap::create(
        const std::vector<NodeGeometry>& sourceGeometry,
        const std::vector<NodeGeometry>& targetGeometry,
        const std::vector<NodeContribution>& contributions,
        const std::vector<NodeOrigin>& origins, NodeRemapReport& report,
        double relativeTolerance) {
    report = NodeRemapReport() ;
    report.sourceNodes = sourceGeometry.size() ;
    report.targetNodes = targetGeometry.size() ;
    report.contributions = contributions.size() ;
    if(relativeTolerance < 0.0 || !std::isfinite(relativeTolerance)) {
      report.invalidGeometry++ ;
      return CPTR<NodeRemap>() ;
    }

    CPTR<NodeRemap> remap = new NodeRemap ;
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
    std::vector<NodeId> soleSource(remap->targetGeometry_.size(), 0) ;
    std::set<std::pair<NodeId, NodeId>> uniqueContributions ;
    for(size_t i=0;i<remap->contributions_.size();++i) {
      const NodeContribution& contribution = remap->contributions_[i] ;
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
      soleSource[target] = contribution.sourceNode ;
      positions[target] += contribution.weight*
        remap->sourceGeometry_[source].position ;
    }

    std::map<NodeId, NodeOrigin> originByTarget ;
    for(size_t i=0;i<remap->origins_.size();++i) {
      if (nodeGeometryIndex(
                remap->targetGeometry_, remap->origins_[i].targetNode) < 0 ||
            remap->origins_[i].kind < node_origin::base_node ||
            remap->origins_[i].kind > node_origin::cell ||
            !originByTarget
                   .insert(std::make_pair(
                         remap->origins_[i].targetNode, remap->origins_[i]))
                   .second)
        report.invalidOrigins++ ;
    }

    for(size_t target=0;target<remap->targetGeometry_.size();++target) {
      const NodeId targetNode = remap->targetGeometry_[target].node ;
      if(targetDegrees[target] == 0)
        report.missingTargetNodes++ ;
      const std::map<NodeId, NodeOrigin>::const_iterator origin =
            originByTarget.find(targetNode) ;
      if(origin == originByTarget.end()) {
        report.invalidOrigins++ ;
      } else if (origin->second.kind == node_origin::base_node) {
        if(targetDegrees[target] != 1) {
          report.invalidOrigins++ ;
        } else if (soleSource[target] != targetNode)
          report.invalidOrigins++ ;
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
      return CPTR<NodeRemap>() ;
    return remap ;
  }

  bool NodeRemap::nodeContributions(
        NodeId targetNode, size_t& begin, size_t& end) const {
    const int target = nodeGeometryIndex(targetGeometry_,targetNode) ;
    if(target < 0)
      return false ;
    begin = targetOffsets_[target] ;
    end = targetOffsets_[target+1] ;
    return true ;
  }

  bool NodeRemap::nodeOrigin(NodeId targetNode, NodeOrigin& origin) const {
    NodeOrigin key ;
    key.targetNode = targetNode ;
    const std::vector<NodeOrigin>::const_iterator location = std::lower_bound(
          origins_.begin(), origins_.end(), key, nodeOriginOrder) ;
    if(location == origins_.end() || location->targetNode != targetNode)
      return false ;
    origin = *location ;
    return true ;
  }

  bool NodeRemap::interpolateNodeData(const std::vector<double>& sourceValues,
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

  bool NodeRemap::interpolateNodeData(
        const std::vector<vector3d<double>>& sourceValues,
        std::vector<vector3d<double>>& targetValues) const {
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

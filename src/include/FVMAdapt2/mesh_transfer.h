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
#ifndef FVMADAPT2_MESH_TRANSFER_H
#define FVMADAPT2_MESH_TRANSFER_H

#include <Tools/basic_types.h>
#include <Tools/cptr.h>

#include <cstddef>
#include <vector>

namespace Loci {

  /// Geometry for one face in either the source or target mesh.
  struct AMRFaceGeometry {
    int face ;
    double area ;
    vector3d<double> centroid ;

    AMRFaceGeometry()
      : face(0), area(0.0), centroid(0.0,0.0,0.0) {}
    AMRFaceGeometry(int faceId, double faceArea,
                    const vector3d<double>& faceCentroid)
      : face(faceId), area(faceArea), centroid(faceCentroid) {}
  } ;

  /// Geometric contribution from one source face to one target face.
  ///
  /// orientation is +1 when their directed normals agree and -1 when they
  /// oppose. overlapArea is deliberately unnormalized so consumers can form
  /// either source or target fractions without losing geometric information.
  struct AMRFaceContribution {
    int sourceFace ;
    int targetFace ;
    double overlapArea ;
    vector3d<double> overlapCentroid ;
    int orientation ;

    AMRFaceContribution()
      : sourceFace(0), targetFace(0), overlapArea(0.0),
        overlapCentroid(0.0,0.0,0.0), orientation(1) {}
    AMRFaceContribution(int source, int target, double area,
                        const vector3d<double>& centroid,
                        int orientationSign)
      : sourceFace(source), targetFace(target), overlapArea(area),
        overlapCentroid(centroid), orientation(orientationSign) {}
  } ;

  /// A target face created inside one source cell during refinement.
  struct AMRCreatedFace {
    int targetFace ;
    int sourceCell ;

    AMRCreatedFace() : targetFace(0), sourceCell(0) {}
    AMRCreatedFace(int face, int cell)
      : targetFace(face), sourceCell(cell) {}
  } ;

  /// A source face removed inside one target cell during derefinement.
  struct AMRRemovedFace {
    int sourceFace ;
    int targetCell ;

    AMRRemovedFace() : sourceFace(0), targetCell(0) {}
    AMRRemovedFace(int face, int cell)
      : sourceFace(face), targetCell(cell) {}
  } ;

  struct AMRFaceRemapReport {
    bool valid ;
    size_t sourceFaces ;
    size_t targetFaces ;
    size_t contributions ;
    size_t createdFaces ;
    size_t removedFaces ;
    size_t invalidGeometry ;
    size_t duplicateContributions ;
    size_t missingSourceFaces ;
    size_t missingTargetFaces ;
    size_t inconsistentSourceMoments ;
    size_t inconsistentTargetMoments ;
    double maximumSourceAreaError ;
    double maximumTargetAreaError ;
    double maximumSourceCentroidError ;
    double maximumTargetCentroidError ;

    AMRFaceRemapReport() ;
  } ;

  /// Immutable geometric relation between source and target mesh faces.
  class AMRFaceRemap : public CPTR_type {
  public:
    static CPTR<AMRFaceRemap>
    create(const std::vector<AMRFaceGeometry>& sourceGeometry,
           const std::vector<AMRFaceGeometry>& targetGeometry,
           const std::vector<AMRFaceContribution>& contributions,
           const std::vector<AMRCreatedFace>& createdFaces,
           const std::vector<AMRRemovedFace>& removedFaces,
           AMRFaceRemapReport& report,
           double relativeTolerance = 1.0e-10) ;

    const std::vector<AMRFaceGeometry>& sourceFaceGeometry() const {
      return sourceGeometry_ ;
    }
    const std::vector<AMRFaceGeometry>& targetFaceGeometry() const {
      return targetGeometry_ ;
    }
    const std::vector<AMRFaceContribution>& faceContributions() const {
      return contributions_ ;
    }
    const std::vector<AMRCreatedFace>& createdFaces() const {
      return createdFaces_ ;
    }
    const std::vector<AMRRemovedFace>& removedFaces() const {
      return removedFaces_ ;
    }

    bool faceContributions(int targetFace,
                           size_t& begin, size_t& end) const ;
    bool isCreatedFace(int targetFace, int& sourceCell) const ;
    bool isRemovedFace(int sourceFace, int& targetCell) const ;

    /// Assemble source face averages on target faces.
    ///
    /// mapped is false for newly created internal faces. If orientValues is
    /// true, contributions are multiplied by their orientation signs.
    bool remapFaceAverages(const std::vector<double>& sourceValues,
                           std::vector<double>& targetValues,
                           std::vector<unsigned char>& mapped,
                           bool orientValues = true) const ;

    /// Distribute source face integrals to target faces by overlap area.
    bool remapFaceIntegrals(const std::vector<double>& sourceIntegrals,
                            std::vector<double>& targetIntegrals,
                            std::vector<unsigned char>& mapped,
                            bool orientValues = true) const ;

  private:
    AMRFaceRemap() {}

    std::vector<AMRFaceGeometry> sourceGeometry_ ;
    std::vector<AMRFaceGeometry> targetGeometry_ ;
    std::vector<AMRFaceContribution> contributions_ ;
    std::vector<AMRCreatedFace> createdFaces_ ;
    std::vector<AMRRemovedFace> removedFaces_ ;
    std::vector<size_t> targetOffsets_ ;
  } ;

  namespace amr_node_origin {
    enum value {
      retained,
      edge,
      face,
      cell
    } ;
  }

  /// Geometry for one node in either the source or target mesh.
  struct AMRNodeGeometry {
    int node ;
    vector3d<double> position ;

    AMRNodeGeometry() : node(0), position(0.0,0.0,0.0) {}
    AMRNodeGeometry(int nodeId, const vector3d<double>& nodePosition)
      : node(nodeId), position(nodePosition) {}
  } ;

  /// One source-node interpolation weight for a target node.
  struct AMRNodeContribution {
    int sourceNode ;
    int targetNode ;
    double weight ;

    AMRNodeContribution() : sourceNode(0), targetNode(0), weight(0.0) {}
    AMRNodeContribution(int source, int target, double nodeWeight)
      : sourceNode(source), targetNode(target), weight(nodeWeight) {}
  } ;

  /// Topological origin of a target node.
  ///
  /// sourceEntity uses the source mesh numbering for the indicated entity
  /// kind. For a retained node it is the corresponding source node.
  struct AMRNodeOrigin {
    int targetNode ;
    amr_node_origin::value kind ;
    int sourceEntity ;

    AMRNodeOrigin()
      : targetNode(0), kind(amr_node_origin::retained), sourceEntity(0) {}
    AMRNodeOrigin(int target, amr_node_origin::value originKind,
                  int originEntity)
      : targetNode(target), kind(originKind), sourceEntity(originEntity) {}
  } ;

  struct AMRNodeRemapReport {
    bool valid ;
    size_t sourceNodes ;
    size_t targetNodes ;
    size_t contributions ;
    size_t invalidGeometry ;
    size_t duplicateContributions ;
    size_t missingTargetNodes ;
    size_t invalidOrigins ;
    size_t inconsistentWeights ;
    size_t inconsistentPositions ;
    double maximumWeightError ;
    double maximumPositionError ;

    AMRNodeRemapReport() ;
  } ;

  /// Immutable source-node interpolation data for the target mesh nodes.
  class AMRNodeRemap : public CPTR_type {
  public:
    static CPTR<AMRNodeRemap>
    create(const std::vector<AMRNodeGeometry>& sourceGeometry,
           const std::vector<AMRNodeGeometry>& targetGeometry,
           const std::vector<AMRNodeContribution>& contributions,
           const std::vector<AMRNodeOrigin>& origins,
           AMRNodeRemapReport& report,
           double relativeTolerance = 1.0e-10) ;

    const std::vector<AMRNodeGeometry>& sourceNodeGeometry() const {
      return sourceGeometry_ ;
    }
    const std::vector<AMRNodeGeometry>& targetNodeGeometry() const {
      return targetGeometry_ ;
    }
    const std::vector<AMRNodeContribution>& nodeContributions() const {
      return contributions_ ;
    }
    const std::vector<AMRNodeOrigin>& nodeOrigins() const {
      return origins_ ;
    }

    bool nodeContributions(int targetNode,
                           size_t& begin, size_t& end) const ;
    bool nodeOrigin(int targetNode, AMRNodeOrigin& origin) const ;

    bool interpolateNodeData(const std::vector<double>& sourceValues,
                             std::vector<double>& targetValues) const ;
    bool interpolateNodeData(
      const std::vector<vector3d<double> >& sourceValues,
      std::vector<vector3d<double> >& targetValues) const ;

  private:
    AMRNodeRemap() {}

    std::vector<AMRNodeGeometry> sourceGeometry_ ;
    std::vector<AMRNodeGeometry> targetGeometry_ ;
    std::vector<AMRNodeContribution> contributions_ ;
    std::vector<AMRNodeOrigin> origins_ ;
    std::vector<size_t> targetOffsets_ ;
  } ;
}

#endif

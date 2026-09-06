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
#include <cstdint>
#include <vector>

namespace Loci {

  /// Persistent identities used by mesh-transition relations.
  ///
  /// These values are derived from canonical refinement-tree keys rather than
  /// from Loci entity numbers or the generated mesh's file order, both of
  /// which may change when that mesh is repartitioned.
  // Loci's long-long schema is the portable signed 64-bit store schema on
  // supported builds.  std::int64_t aliases long on LP64 systems, whose
  // legacy Loci schema is not a signed 64-bit atomic type.
  typedef long long FaceId ;
  typedef long long CellId ;
  typedef long long NodeId ;
  static_assert(sizeof(FaceId) == 8,
        "persistent AMR identities require 64-bit long long") ;
  static_assert(sizeof(CellId) == 8,
        "persistent AMR identities require 64-bit long long") ;
  static_assert(sizeof(NodeId) == 8,
        "persistent AMR identities require 64-bit long long") ;

  namespace node_transition_status {
    enum value {
      valid,
      invalid_tolerance,
      unsupported_topology,
      missing_state,
      missing_source_node,
      cyclic_construction,
      inconsistent_weights,
      inconsistent_positions,
      not_requested
    } ;
  }

  /// Validation result for one previous-to-current node transition.
  struct NodeTransitionReport {
    node_transition_status::value status ;
    bool valid ;
    size_t sourceNodes ;
    size_t targetNodes ;
    size_t retainedNodes ;
    size_t createdNodes ;
    size_t contributions ;
    size_t missingSourceNodes ;
    size_t cyclicConstructions ;
    size_t inconsistentWeights ;
    size_t inconsistentPositions ;
    double maximumWeightError ;
    double maximumPositionError ;

    NodeTransitionReport() ;
  } ;


  /// Geometry for one face in either the source or target mesh.
  struct FaceGeometry {
    FaceId face ;
    double area ;
    vector3d<double> centroid ;

    FaceGeometry() : face(0), area(0.0), centroid(0.0, 0.0, 0.0) {}
    FaceGeometry(
          FaceId faceId, double faceArea, const vector3d<double>& faceCentroid)
        : face(faceId), area(faceArea), centroid(faceCentroid) {}
  } ;

  /// Shared geometry between a face before and a face after adaptation.
  ///
  /// orientation is +1 when their directed normals agree and -1 when they
  /// oppose. area is deliberately unnormalized so consumers can form
  /// either source or target fractions without losing geometric information.
  struct FaceOverlap {
    FaceId source ;
    FaceId target ;
    double area ;
    vector3d<double> centroid ;
    int orientation ;

    FaceOverlap()
        : source(0), target(0), area(0.0), centroid(0.0, 0.0, 0.0),
          orientation(1) {}
    FaceOverlap(FaceId sourceId, FaceId targetId, double overlapArea,
          const vector3d<double>& overlapCentroid, int orientationSign)
        : source(sourceId), target(targetId), area(overlapArea),
          centroid(overlapCentroid), orientation(orientationSign) {}
  } ;

  /// A face created inside a cell of the immediately previous mesh.
  /// sourceCell is not necessarily a base-mesh root cell.
  struct CreatedFace {
    FaceId targetFace ;
    CellId sourceCell ;

    CreatedFace() : targetFace(0), sourceCell(0) {}
    CreatedFace(FaceId face, CellId cell)
        : targetFace(face), sourceCell(cell) {}
  } ;

  /// A previous-mesh face removed inside a resulting cell during derefinement.
  struct RemovedFace {
    FaceId sourceFace ;
    CellId targetCell ;

    RemovedFace() : sourceFace(0), targetCell(0) {}
    RemovedFace(FaceId face, CellId cell)
        : sourceFace(face), targetCell(cell) {}
  } ;

  struct FaceRemapReport {
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

    FaceRemapReport() ;
  } ;

  namespace face_transition_status {
    enum value {
      available,
      unsupported_restart,
      unsupported_topology,
      unsupported_plan_change,
      invalid_identity,
      invalid_geometry,
      inconsistent_relation
    } ;
  }

  /// Construction status for the live face transition.
  struct FaceTransitionReport {
    face_transition_status::value status ;
    bool valid ;
    size_t sourceFaces ;
    size_t targetFaces ;
    size_t unsupportedRootPlans ;
    size_t invalidIdentities ;
    size_t invalidPolygons ;
    FaceRemapReport remap ;

    FaceTransitionReport() ;
  } ;

  /// Opaque history retained between in-process mesh transitions.
  ///
  /// Solvers should carry this handle unchanged. The refinement-tree
  /// representation is private to FVMAdapt2.
  class MeshState : public CPTR_type {
  protected:
    MeshState() {}
  } ;

  /// Immutable geometric relation for locally owned target mesh faces.
  /// Source means the mesh before adaptation; target means the resulting mesh.
  /// Source geometry contains the subset required by this rank and is not
  /// necessarily locally owned on the source mesh.
  class FaceRemap : public CPTR_type {
  public:
    static CPTR<FaceRemap> create(
          const std::vector<FaceGeometry>& sourceGeometry,
          const std::vector<FaceGeometry>& targetGeometry,
          const std::vector<FaceOverlap>& contributions,
          const std::vector<CreatedFace>& createdFaces,
          const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
          double relativeTolerance = 1.0e-10) ;

    /// Build a target-owned view of an already globally validated relation.
    ///
    /// Source faces may be present only for contributions to locally owned
    /// targets, so source-side coverage is intentionally not revalidated.
    /// Target coverage and all local identities remain fully checked.
    static CPTR<FaceRemap> createTargetOwned(
          const std::vector<FaceGeometry>& sourceGeometry,
          const std::vector<FaceGeometry>& targetGeometry,
          const std::vector<FaceOverlap>& contributions,
          const std::vector<CreatedFace>& createdFaces,
          const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
          double relativeTolerance = 1.0e-10) ;

    const std::vector<FaceGeometry>& sourceFaceGeometry() const {
      return sourceGeometry_ ;
    }
    const std::vector<FaceGeometry>& targetFaceGeometry() const {
      return targetGeometry_ ;
    }
    const std::vector<FaceOverlap>& overlaps() const { return contributions_ ; }
    const std::vector<CreatedFace>& createdFaces() const {
      return createdFaces_ ;
    }
    const std::vector<RemovedFace>& removedFaces() const {
      return removedFaces_ ;
    }

    /// Rows for one target face; a created face has a valid empty row range.
    bool overlaps(FaceId targetFace, size_t& begin, size_t& end) const ;

    /// Return locally owned target faces that receive a contribution from a
    /// source face in this target-owned remap.
    void targetFaces(FaceId sourceFace, std::vector<FaceId>& targets) const ;

    bool isCreatedFace(FaceId targetFace, CellId& sourceCell) const ;
    bool isRemovedFace(FaceId sourceFace, CellId& targetCell) const ;

    /// Assemble source face averages on locally owned target faces.
    ///
    /// mapped is false for newly created internal faces. Set orientValues for
    /// directed quantities such as fluxes; ordinary scalar averages are
    /// assembled without an orientation sign by default. sourceValues follow
    /// sourceFaceGeometry() order; targetValues follow targetFaceGeometry().
    bool remapFaceAverages(const std::vector<double>& sourceValues,
          std::vector<double>& targetValues, std::vector<unsigned char>& mapped,
          bool orientValues = false) const ;

    /// Distribute source face integrals to locally owned target faces by
    /// overlap area. Input and output ordering matches sourceFaceGeometry()
    /// and targetFaceGeometry(), respectively.
    bool remapFaceIntegrals(const std::vector<double>& sourceIntegrals,
          std::vector<double>& targetIntegrals,
          std::vector<unsigned char>& mapped, bool orientValues = false) const ;

  private:
    FaceRemap() {}

    static CPTR<FaceRemap> createImpl(
          const std::vector<FaceGeometry>& sourceGeometry,
          const std::vector<FaceGeometry>& targetGeometry,
          const std::vector<FaceOverlap>& contributions,
          const std::vector<CreatedFace>& createdFaces,
          const std::vector<RemovedFace>& removedFaces, FaceRemapReport& report,
          double relativeTolerance, bool validateSourceCoverage) ;

    std::vector<FaceGeometry> sourceGeometry_ ;
    std::vector<FaceGeometry> targetGeometry_ ;
    std::vector<FaceOverlap> contributions_ ;
    std::vector<CreatedFace> createdFaces_ ;
    std::vector<RemovedFace> removedFaces_ ;
    std::vector<size_t> targetOffsets_ ;
    std::vector<size_t> sourceOffsets_ ;
    std::vector<FaceId> sourceTargets_ ;
  } ;

  namespace node_origin {
    enum value { base_node, edge, face, cell } ;
  }

  /// Geometry for one node in either the source or target mesh.
  struct NodeGeometry {
    NodeId node ;
    vector3d<double> position ;

    NodeGeometry() : node(0), position(0.0, 0.0, 0.0) {}
    NodeGeometry(NodeId nodeId, const vector3d<double>& nodePosition)
        : node(nodeId), position(nodePosition) {}
  } ;

  /// One source-node interpolation weight for a target node.
  struct NodeContribution {
    NodeId sourceNode ;
    NodeId targetNode ;
    double weight ;

    NodeContribution() : sourceNode(0), targetNode(0), weight(0.0) {}
    NodeContribution(NodeId source, NodeId target, double nodeWeight)
        : sourceNode(source), targetNode(target), weight(nodeWeight) {}
  } ;

  /// Topological origin of a target node.
  struct NodeOrigin {
    NodeId targetNode ;
    node_origin::value kind ;

    NodeOrigin() : targetNode(0), kind(node_origin::base_node) {}
    NodeOrigin(NodeId target, node_origin::value originKind)
        : targetNode(target), kind(originKind) {}
  } ;

  struct NodeRemapReport {
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

    NodeRemapReport() ;
  } ;

  /// Immutable source-node interpolation data for locally owned target nodes.
  /// Source geometry contains the subset needed by this rank and is not
  /// necessarily locally owned on the source mesh.
  class NodeRemap : public CPTR_type {
  public:
    static CPTR<NodeRemap> create(
          const std::vector<NodeGeometry>& sourceGeometry,
          const std::vector<NodeGeometry>& targetGeometry,
          const std::vector<NodeContribution>& contributions,
          const std::vector<NodeOrigin>& origins, NodeRemapReport& report,
          double relativeTolerance = 1.0e-10) ;

    const std::vector<NodeGeometry>& sourceNodeGeometry() const {
      return sourceGeometry_ ;
    }
    const std::vector<NodeGeometry>& targetNodeGeometry() const {
      return targetGeometry_ ;
    }
    const std::vector<NodeContribution>& nodeContributions() const {
      return contributions_ ;
    }
    const std::vector<NodeOrigin>& nodeOrigins() const { return origins_ ; }

    bool nodeContributions(NodeId targetNode, size_t& begin, size_t& end) const ;
    bool nodeOrigin(NodeId targetNode, NodeOrigin& origin) const ;

    /// Input follows sourceNodeGeometry() order; output follows
    /// targetNodeGeometry() order.
    bool interpolateNodeData(const std::vector<double>& sourceValues,
                             std::vector<double>& targetValues) const ;
    bool interpolateNodeData(
      const std::vector<vector3d<double> >& sourceValues,
      std::vector<vector3d<double> >& targetValues) const ;

  private:
    NodeRemap() {}

    std::vector<NodeGeometry> sourceGeometry_ ;
    std::vector<NodeGeometry> targetGeometry_ ;
    std::vector<NodeContribution> contributions_ ;
    std::vector<NodeOrigin> origins_ ;
    std::vector<size_t> targetOffsets_ ;
  } ;
}

#endif

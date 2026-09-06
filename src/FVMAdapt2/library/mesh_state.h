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
#ifndef FVMADAPT2_MESH_STATE_H
#define FVMADAPT2_MESH_STATE_H

#include <FVMAdapt2/defines.h>
#include <FVMAdapt2/node_transition.h>
#include <Loci.h>

#include <Tools/cptr.h>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace Loci {

  namespace face_origin {
    enum value { base_face, cell_interior } ;
  }

  /// Topology of the root cell that produced a refinement tree.
  namespace cell_topology {
    enum value { hex, prism, general } ;
  }

  /// Topological identity of a face independent of its current mesh number.
  ///
  /// A base-face key consists of the original face file number and its leaf
  /// path; its second path is empty. A cell-interior key consists of the
  /// original cell file number and two distinct accepted cell-leaf paths in
  /// sorted order. Each path stores its refinement steps as consecutive
  /// (split code, child ordinal) pairs.
  struct FaceKey {
    face_origin::value origin ;
    int root ;
    std::vector<int> firstPath ;
    std::vector<int> secondPath ;

    FaceKey() : origin(face_origin::base_face), root(0) {}
    FaceKey(face_origin::value faceOrigin, int rootEntity,
          const std::vector<int>& first,
          const std::vector<int>& second = std::vector<int>())
        : origin(faceOrigin), root(rootEntity), firstPath(first),
          secondPath(second) {}
  } ;

  bool operator<(const FaceKey& left, const FaceKey& right) ;
  bool operator==(const FaceKey& left, const FaceKey& right) ;

  /// Topological identity of a leaf cell independent of mesh numbering.
  /// path uses the same flattened (split code, child ordinal) representation.
  struct CellKey {
    int root ;
    std::vector<int> path ;

    CellKey() : root(0) {}
    CellKey(int rootCell, const std::vector<int>& leafPath)
        : root(rootCell), path(leafPath) {}
  } ;

  bool operator<(const CellKey& left, const CellKey& right) ;
  bool operator==(const CellKey& left, const CellKey& right) ;

  /// Deterministically encode a canonical face key as a persistent identity.
  ///
  /// The encoding is specified by FVMAdapt2 and is independent of process
  /// ownership, generated entity numbers, and VOG face ordering.
  FaceId persistentFaceId(const FaceKey& key) ;

  /// Deterministically identify a leaf cell from its base cell and tree path.
  CellId persistentCellId(const CellKey& key) ;
  CellId persistentCellId(int rootCell, const std::vector<int>& leafPath) ;

  /// Association between a persistent face key and one mesh face number.
  struct FaceIdentity {
    int face ;
    FaceId id ;
    FaceKey key ;

    FaceIdentity() : face(0), id(0) {}
    FaceIdentity(int faceNumber, const FaceKey& faceKey)
        : face(faceNumber), id(persistentFaceId(faceKey)), key(faceKey) {}
  } ;

  /// Association between a persistent cell key and one mesh cell number.
  struct CellIdentity {
    int cell ;
    CellId id ;
    CellKey key ;

    CellIdentity() : cell(0), id(0) {}
    CellIdentity(int cellNumber, const CellKey& cellKey)
        : cell(cellNumber), id(persistentCellId(cellKey)), key(cellKey) {}
  } ;

  /// Accepted leaf-cell topology for one original cell.
  ///
  /// leafPaths contains one canonical refinement path for every accepted
  /// leaf. It is independent of cell numbering and of the cell-plan encoding.
  struct RootCellState {
    int root ;
    cell_topology::value topology ;
    std::vector<std::vector<int>> leafPaths ;

    RootCellState() : root(0), topology(cell_topology::hex) {}
    RootCellState(int rootCell,
          const std::vector<std::vector<int>>& acceptedLeafPaths,
          cell_topology::value rootTopology = cell_topology::hex)
        : root(rootCell), topology(rootTopology), leafPaths(acceptedLeafPaths) {
    }
  } ;

  /// Immutable face identity and geometry for one accepted mesh state.
  ///
  /// Face polygons are retained internally so an in-process subsequent AMR
  /// cycle can construct ancestry without guessing from neighboring cells.
  class FaceState : public MeshState {
  public:
    static CPTR<FaceState> create(const std::vector<FaceIdentity>& identities,
          const std::vector<std::vector<vector3d<double>>>& polygons,
          const std::vector<RootCellState>& rootCells,
          FaceTransitionReport& report, double relativeTolerance = 1.0e-10) ;

    const std::vector<FaceIdentity>& faceIdentities() const {
      return identities_ ;
    }
    const std::vector<RootCellState>& rootCells() const { return rootCells_ ; }

    /// Own a snapshot: installing the grid consumes its temporary stores.
    void retainNodes(const store<NodeId>& ids,
        const store<vector3d<double>>& positions,
        const std::vector<entitySet>& partition) ;

    // Canonical node ownership is retained with the face history so the next
    // cycle need not rebuild the previous mesh just to recover its nodes.
    store<NodeId> nodeIds ;
    store<vector3d<double>> nodePositions ;
    std::vector<entitySet> nodePartition ;

  private:
    FaceState() {}

    std::vector<FaceIdentity> identities_ ;
    std::vector<std::vector<vector3d<double>>> polygons_ ;
    std::vector<RootCellState> rootCells_ ;

    friend CPTR<FaceRemap> buildFaceRemap(const_CPTR<FaceState>,
          const_CPTR<FaceState>, FaceTransitionReport&, double) ;
  } ;

  /// Construct the face relation between two accepted mesh states.
  CPTR<FaceRemap> buildFaceRemap(const_CPTR<FaceState> source,
        const_CPTR<FaceState> target, FaceTransitionReport& report,
        double relativeTolerance = 1.0e-10) ;

  namespace detail {
    // Exchange fixed-size records between old/new mesh owners. Persistent
    // IDs are not Loci entity numbers and cannot directly index a store.
    template <class T>
    std::vector<T> exchangeByDestination(
          const std::vector<std::vector<T>>& outgoing) {
      const int processes = MPI_processes ;
      std::vector<int> sendCounts(processes, 0) ;
      std::vector<int> receiveCounts(processes, 0) ;
      std::vector<int> sendOffsets(processes, 0) ;
      std::vector<int> receiveOffsets(processes, 0) ;
      for (int process = 0; process < processes; ++process)
        sendCounts[process] = int(outgoing[process].size()) ;
      MPI_Alltoall(sendCounts.data(), 1, MPI_INT, receiveCounts.data(), 1,
            MPI_INT, MPI_COMM_WORLD) ;
      for (int process = 1; process < processes; ++process) {
        sendOffsets[process] =
              sendOffsets[process - 1] + sendCounts[process - 1] ;
        receiveOffsets[process] =
              receiveOffsets[process - 1] + receiveCounts[process - 1] ;
      }
      const int sendSize =
            processes == 0 ? 0 : sendOffsets.back() + sendCounts.back() ;
      const int receiveSize =
            processes == 0 ? 0 : receiveOffsets.back() + receiveCounts.back() ;
      std::vector<T> sendRecords ;
      std::vector<T> receiveRecords ;
      sendRecords.resize(size_t(sendSize)) ;
      receiveRecords.resize(size_t(receiveSize)) ;
      for (int process = 0; process < processes; ++process)
        std::copy(outgoing[process].begin(), outgoing[process].end(),
              sendRecords.begin() + sendOffsets[process]) ;
      MPI_Datatype recordType ;
      MPI_Type_contiguous(int(sizeof(T)), MPI_BYTE, &recordType) ;
      MPI_Type_commit(&recordType) ;
      MPI_Alltoallv(sendRecords.empty() ? 0 : sendRecords.data(),
            sendCounts.data(), sendOffsets.data(), recordType,
            receiveRecords.empty() ? 0 : receiveRecords.data(),
            receiveCounts.data(), receiveOffsets.data(), recordType,
            MPI_COMM_WORLD) ;
      MPI_Type_free(&recordType) ;
      return receiveRecords ;
    }

    /// Encode normalized leaf paths for transport as a Loci fact.
    bool encodeLeafPaths(const std::vector<std::vector<int>>& paths,
          std::vector<int>& encodedPaths) ;

    /// Encode breadth-first leaf paths for a quadrilateral face plan.
    ///
    /// The Loci fact is [leaf count, step count, code, child, ...] for each
    /// leaf, so a path remains meaningful when a different split code uses
    /// the same child ordinal.
    bool encodeQuadFaceLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) ;

    /// Encode breadth-first leaf paths for a general polygonal face plan.
    bool encodeGeneralFaceLeafPaths(const std::vector<char>& plan,
          int initialEdgeCount, std::vector<int>& encodedPaths) ;

    /// Encode breadth-first leaf paths for a hexahedral cell plan.
    bool encodeHexCellLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) ;

    /// Encode breadth-first leaf paths for a triangular-prism cell plan.
    bool encodePrismCellLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) ;

    /// Decode the compact Loci-fact representation into flattened
    /// (split code, child ordinal) paths.
    bool decodeLeafPaths(const std::vector<int>& encodedPaths,
          std::vector<std::vector<int>>& paths) ;

    /// Validate persistent-ID collisions without replicating all face keys.
    bool validateFaceIdentityHashesDistributed(
          const std::vector<FaceIdentity>& identities,
          size_t& invalidIdentities) ;

    /// Validate persistent cell IDs without replicating all cell keys.
    bool validateCellIdentityHashesDistributed(
          const std::vector<CellIdentity>& identities,
          size_t& invalidIdentities) ;
    /// Original cell/face and FineFaces row bound to each assembled mesh face.
    struct FaceSources {
      Map cell ;
      Map face ;
      store<int> ordinal ;
    } ;

    bool collectOriginalFaceState(
          fact_db& facts, CPTR<FaceState>& state, FaceTransitionReport& report) ;

    bool collectAcceptedFaceState(int numNodes, int numFaces,
          const store<FineFaces>& fineFacesCell,
          const store<FineFaces>& fineFaces, const FaceSources& faceSources,
          const store<vector3d<double>>& positions, const multiMap& faceToNode,
          const std::vector<entitySet>& outputCellPartition, fact_db& facts,
          CPTR<FaceState>& state, store<FaceId>& faceIds,
          store<CellId>& cellIds, FaceTransitionReport& report) ;

    /// Deliver face contributions to the final owners of the generated mesh.
    CPTR<FaceRemap> redistributeFaceRemap(const CPTR<FaceRemap>& remap,
          const FaceTransitionReport& globalReport,
          const std::vector<entitySet>& facePartition,
          const std::vector<entitySet>& cellPartition,
          const store<FaceId>& faceIds, const store<CellId>& cellIds,
          FaceTransitionReport& localReport) ;

    /// Combine independently owned transition groups into one global report.
    void reduceFaceTransitionReport(
          const FaceTransitionReport& local, FaceTransitionReport& global) ;
    /// Build the local target-owned node remap for one accepted mesh change.
    /// Construction parents are cloned by canonical node number; previous
    /// mesh nodes are found by persistent ID through a distributed directory.
    CPTR<NodeRemap> buildNodeRemap(const store<NodeId>& previousIds,
          const store<vector3d<double>>& previousPositions,
          const std::vector<entitySet>& previousPartition,
          const store<NodeId>& currentIds,
          const store<FineNodeConstruction>& currentConstructions,
          const store<vector3d<double>>& currentPositions,
          const std::vector<entitySet>& currentPartition,
          NodeTransitionReport& report, double relativeTolerance = 1.0e-10) ;

    /// Combine local transition counts and errors into one collective report.
    void reduceNodeTransitionReport(
          const NodeTransitionReport& local, NodeTransitionReport& global) ;

    /// Move each target-owned remap row to the final owner of that mesh node.
    /// nodePartition uses the generated node numbering and generatedNodeIds
    /// supplies its persistent identity.
    CPTR<NodeRemap> redistributeNodeRemap(const CPTR<NodeRemap>& remap,
          const std::vector<entitySet>& nodePartition,
          const store<NodeId>& generatedNodeIds, NodeTransitionReport& report,
          double relativeTolerance = 1.0e-10) ;

  }
}

#endif

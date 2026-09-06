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

#include "mesh_state.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <set>

namespace Loci {
  namespace detail {
    namespace {
      template <class T>
      void appendBytes(std::vector<unsigned char>& bytes, const T& value) {
        const size_t offset = bytes.size() ;
        bytes.resize(offset + sizeof(T)) ;
        std::memcpy(&bytes[offset], &value, sizeof(T)) ;
      }

      template <class T>
      bool readBytes(
            const std::vector<unsigned char>& bytes, size_t& offset, T& value) {
        if (sizeof(T) > bytes.size() - offset)
          return false ;
        std::memcpy(&value, &bytes[offset], sizeof(T)) ;
        offset += sizeof(T) ;
        return true ;
      }

      void appendIdentity(
            std::vector<unsigned char>& bytes, const FaceIdentity& identity) {
        appendBytes(bytes, identity.id) ;
        const std::int32_t origin = std::int32_t(identity.key.origin) ;
        const std::int32_t root = std::int32_t(identity.key.root) ;
        const std::int32_t firstSize =
              std::int32_t(identity.key.firstPath.size()) ;
        const std::int32_t secondSize =
              std::int32_t(identity.key.secondPath.size()) ;
        appendBytes(bytes, origin) ;
        appendBytes(bytes, root) ;
        appendBytes(bytes, firstSize) ;
        for (size_t entry = 0; entry < identity.key.firstPath.size(); ++entry) {
          const std::int32_t pathEntry =
                std::int32_t(identity.key.firstPath[entry]) ;
          appendBytes(bytes, pathEntry) ;
        }
        appendBytes(bytes, secondSize) ;
        for (size_t entry = 0; entry < identity.key.secondPath.size();
              ++entry) {
          const std::int32_t pathEntry =
                std::int32_t(identity.key.secondPath[entry]) ;
          appendBytes(bytes, pathEntry) ;
        }
      }

      void appendCellIdentity(
            std::vector<unsigned char>& bytes, const CellIdentity& identity) {
        appendBytes(bytes, identity.id) ;
        const std::int32_t root = std::int32_t(identity.key.root) ;
        const std::int32_t pathSize = std::int32_t(identity.key.path.size()) ;
        appendBytes(bytes, root) ;
        appendBytes(bytes, pathSize) ;
        for (size_t entry = 0; entry < identity.key.path.size(); ++entry) {
          const std::int32_t pathEntry = std::int32_t(identity.key.path[entry]) ;
          appendBytes(bytes, pathEntry) ;
        }
      }

      void appendNodeIdentity(std::vector<unsigned char>& bytes, NodeId id,
            long long fileNumber) {
        appendBytes(bytes, id) ;
        appendBytes(bytes, fileNumber) ;
      }

      bool readIdentity(const std::vector<unsigned char>& bytes, size_t& offset,
            FaceId& id, FaceKey& key) {
        std::int32_t origin = 0 ;
        std::int32_t root = 0 ;
        std::int32_t firstSize = 0 ;
        std::int32_t secondSize = 0 ;
        if (!readBytes(bytes, offset, id) ||
              !readBytes(bytes, offset, origin) ||
              !readBytes(bytes, offset, root) ||
              !readBytes(bytes, offset, firstSize) || firstSize < 0 ||
              size_t(firstSize) >
                    (bytes.size() - offset) / sizeof(std::int32_t))
          return false ;
        key = FaceKey() ;
        key.origin = static_cast<face_origin::value>(origin) ;
        key.root = int(root) ;
        key.firstPath.resize(size_t(firstSize)) ;
        for (size_t entry = 0; entry < key.firstPath.size(); ++entry) {
          std::int32_t pathEntry = 0 ;
          if (!readBytes(bytes, offset, pathEntry))
            return false ;
          key.firstPath[entry] = int(pathEntry) ;
        }
        if (!readBytes(bytes, offset, secondSize) || secondSize < 0 ||
              size_t(secondSize) >
                    (bytes.size() - offset) / sizeof(std::int32_t))
          return false ;
        key.secondPath.resize(size_t(secondSize)) ;
        for (size_t entry = 0; entry < key.secondPath.size(); ++entry) {
          std::int32_t pathEntry = 0 ;
          if (!readBytes(bytes, offset, pathEntry))
            return false ;
          key.secondPath[entry] = int(pathEntry) ;
        }
        return id != 0 &&
               (origin == std::int32_t(face_origin::base_face) ||
                     origin == std::int32_t(face_origin::cell_interior)) &&
               persistentFaceId(key) == id ;
      }

      bool readCellIdentity(const std::vector<unsigned char>& bytes,
            size_t& offset, CellId& id, CellKey& key) {
        std::int32_t root = 0 ;
        std::int32_t pathSize = 0 ;
        if (!readBytes(bytes, offset, id) || !readBytes(bytes, offset, root) ||
              !readBytes(bytes, offset, pathSize) || pathSize < 0 ||
              size_t(pathSize) > (bytes.size() - offset) / sizeof(std::int32_t))
          return false ;
        key = CellKey() ;
        key.root = int(root) ;
        key.path.resize(size_t(pathSize)) ;
        for (size_t entry = 0; entry < key.path.size(); ++entry) {
          std::int32_t pathEntry = 0 ;
          if (!readBytes(bytes, offset, pathEntry))
            return false ;
          key.path[entry] = int(pathEntry) ;
        }
        return id != 0 && persistentCellId(key) == id ;
      }

      bool readNodeIdentity(const std::vector<unsigned char>& bytes,
            size_t& offset, NodeId& id, long long& fileNumber) {
        return readBytes(bytes, offset, id) &&
               readBytes(bytes, offset, fileNumber) && id != 0 &&
               fileNumber >= 0 && persistentBaseNodeId(fileNumber) == id ;
      }

      template <class FaceMap, class PositionStore>
      std::vector<vector3d<double>> facePolygon(int face,
            const FaceMap& faceToNode, const PositionStore& positions) {
        std::vector<vector3d<double>> polygon(faceToNode[face].size()) ;
        for (int node = 0; node < faceToNode[face].size(); ++node)
          polygon[node] = positions[faceToNode[face][node]] ;
        return polygon ;
      }

      template <class PositionStore, class FaceMap>
      bool expandPositionsForFaces(const PositionStore& positions,
            const FaceMap& faceToNode, const entitySet& ownedNodes,
            dstore<vector3d<double>>& expandedPositions) {
        entitySet owned = ownedNodes ;
        std::vector<entitySet> nodePartition =
              all_collect_vectors(owned, MPI_COMM_WORLD) ;
        entitySet requiredNodes =
              MapRepP(faceToNode.Rep())->image(faceToNode.domain()) ;
        FORALL(positions.domain(), node) {
          expandedPositions[node] = positions[node] ;
        }
        ENDFORALL ;
        storeRepP expandedRep = expandedPositions.Rep() ;
        entitySet requested = requiredNodes - positions.domain() ;
        fill_clone(expandedRep, requested, nodePartition) ;
        expandedPositions.setRep(expandedRep) ;
        return (requiredNodes - expandedPositions.domain()).size() == 0 ;
      }

      FaceKey interiorKey(
            int root, std::vector<int> leftPath, std::vector<int> rightPath) {
        if (rightPath < leftPath)
          leftPath.swap(rightPath) ;
        return FaceKey(face_origin::cell_interior, root, leftPath, rightPath) ;
      }

      bool rootCellTopology(Entity cell, const entitySet& hexCells,
            const entitySet& prismCells, const entitySet& generalCells,
            cell_topology::value& topology) {
        const bool isHex = hexCells.inSet(cell) ;
        const bool isPrism = prismCells.inSet(cell) ;
        const bool isGeneral = generalCells.inSet(cell) ;
        if (int(isHex) + int(isPrism) + int(isGeneral) != 1)
          return false ;
        topology =
              isHex ? cell_topology::hex
                    : (isPrism ? cell_topology::prism : cell_topology::general) ;
        return true ;
      }
    }

    bool validateFaceIdentityHashesDistributed(
          const std::vector<FaceIdentity>& identities,
          size_t& invalidIdentities) {
      invalidIdentities = 0 ;
      const int processCount = MPI_processes ;
      std::vector<std::vector<unsigned char>> perDestination(processCount) ;
      for (size_t entry = 0; entry < identities.size(); ++entry) {
        const std::uint64_t bits =
              static_cast<std::uint64_t>(identities[entry].id) ;
        const int destination = int(bits % std::uint64_t(processCount)) ;
        appendIdentity(perDestination[destination], identities[entry]) ;
      }

      const std::vector<unsigned char> receiveBytes =
            exchangeByDestination(perDestination) ;

      std::map<FaceId, FaceKey> byId ;
      size_t offset = 0 ;
      while (offset < receiveBytes.size()) {
        FaceId id = 0 ;
        FaceKey key ;
        if (!readIdentity(receiveBytes, offset, id, key)) {
          std::cerr << "rank " << MPI_rank
                    << " could not decode AMR face identity at byte " << offset
                    << " of " << receiveBytes.size() << std::endl ;
          invalidIdentities++ ;
          break ;
        }
        const std::map<FaceId, FaceKey>::const_iterator prior = byId.find(id) ;
        if (prior != byId.end()) {
          // A state owns each canonical face exactly once.  A repeated key is
          // therefore invalid just as surely as a true hash collision.
          std::cerr << "rank " << MPI_rank << " duplicate AMR face id " << id
                    << " prior(origin=" << int(prior->second.origin)
                    << ",root=" << prior->second.root
                    << ",first=" << prior->second.firstPath.size()
                    << ",second=" << prior->second.secondPath.size()
                    << ") current(origin=" << int(key.origin)
                    << ",root=" << key.root << ",first=" << key.firstPath.size()
                    << ",second=" << key.secondPath.size() << ")" << std::endl ;
          invalidIdentities++ ;
        } else {
          byId[id] = key ;
        }
      }
      unsigned long long localInvalid =
            static_cast<unsigned long long>(invalidIdentities) ;
      unsigned long long globalInvalid = 0 ;
      MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, MPI_COMM_WORLD) ;
      invalidIdentities = size_t(globalInvalid) ;
      return invalidIdentities == 0 ;
    }

    bool validateCellIdentityHashesDistributed(
          const std::vector<CellIdentity>& identities,
          size_t& invalidIdentities) {
      invalidIdentities = 0 ;
      const int processCount = MPI_processes ;
      std::vector<std::vector<unsigned char>> perDestination(processCount) ;
      for (size_t entry = 0; entry < identities.size(); ++entry) {
        const std::uint64_t bits =
              static_cast<std::uint64_t>(identities[entry].id) ;
        const int destination = int(bits % std::uint64_t(processCount)) ;
        appendCellIdentity(perDestination[destination], identities[entry]) ;
      }

      const std::vector<unsigned char> receiveBytes =
            exchangeByDestination(perDestination) ;

      std::map<CellId, CellKey> byId ;
      size_t offset = 0 ;
      while (offset < receiveBytes.size()) {
        CellId id = 0 ;
        CellKey key ;
        if (!readCellIdentity(receiveBytes, offset, id, key)) {
          std::cerr << "rank " << MPI_rank
                    << " could not decode AMR cell identity at byte " << offset
                    << " of " << receiveBytes.size() << std::endl ;
          invalidIdentities++ ;
          break ;
        }
        if (byId.find(id) != byId.end())
          invalidIdentities++ ;
        else
          byId[id] = key ;
      }
      unsigned long long localInvalid =
            static_cast<unsigned long long>(invalidIdentities) ;
      unsigned long long globalInvalid = 0 ;
      MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, MPI_COMM_WORLD) ;
      invalidIdentities = size_t(globalInvalid) ;
      return invalidIdentities == 0 ;
    }

    bool validateNodeIdentityHashesDistributed(
          const std::vector<std::pair<NodeId, long long>>& identities,
          size_t& invalidIdentities) {
      invalidIdentities = 0 ;
      const int processCount = MPI_processes ;
      std::vector<std::vector<unsigned char>> perDestination(processCount) ;
      for (size_t entry = 0; entry < identities.size(); ++entry) {
        const std::uint64_t bits =
              static_cast<std::uint64_t>(identities[entry].first) ;
        const int destination = int(bits % std::uint64_t(processCount)) ;
        appendNodeIdentity(perDestination[destination], identities[entry].first,
              identities[entry].second) ;
      }

      const std::vector<unsigned char> receiveBytes =
            exchangeByDestination(perDestination) ;

      std::map<NodeId, long long> byId ;
      size_t offset = 0 ;
      while (offset < receiveBytes.size()) {
        NodeId id = 0 ;
        long long fileNumber = -1 ;
        if (!readNodeIdentity(receiveBytes, offset, id, fileNumber)) {
          std::cerr << "rank " << MPI_rank
                    << " could not decode AMR node identity at byte " << offset
                    << " of " << receiveBytes.size() << std::endl ;
          invalidIdentities++ ;
          break ;
        }
        if (!byId.insert(std::make_pair(id, fileNumber)).second)
          invalidIdentities++ ;
      }
      unsigned long long localInvalid =
            static_cast<unsigned long long>(invalidIdentities) ;
      unsigned long long globalInvalid = 0 ;
      MPI_Allreduce(&localInvalid, &globalInvalid, 1, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, MPI_COMM_WORLD) ;
      invalidIdentities = size_t(globalInvalid) ;
      return invalidIdentities == 0 ;
    }

    void reduceFaceTransitionReport(
          const FaceTransitionReport& local, FaceTransitionReport& global) {
      global = FaceTransitionReport() ;
      int localStatus = int(local.status) ;
      int globalStatus = 0 ;
      int localValid = local.valid ? 1 : 0 ;
      int globalValid = 0 ;
      MPI_Allreduce(
            &localStatus, &globalStatus, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
      MPI_Allreduce(
            &localValid, &globalValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      global.status = static_cast<face_transition_status::value>(globalStatus) ;
      global.valid = globalValid != 0 ;

      unsigned long long localCounts[11] = {
            static_cast<unsigned long long>(local.sourceFaces),
            static_cast<unsigned long long>(local.targetFaces),
            static_cast<unsigned long long>(local.unsupportedRootPlans),
            static_cast<unsigned long long>(local.invalidIdentities),
            static_cast<unsigned long long>(local.invalidPolygons),
            static_cast<unsigned long long>(local.remap.contributions),
            static_cast<unsigned long long>(local.remap.createdFaces),
            static_cast<unsigned long long>(local.remap.removedFaces),
            static_cast<unsigned long long>(local.remap.duplicateContributions),
            static_cast<unsigned long long>(local.remap.missingSourceFaces),
            static_cast<unsigned long long>(local.remap.missingTargetFaces)} ;
      unsigned long long globalCounts[11] = {0} ;
      MPI_Allreduce(localCounts, globalCounts, 11, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, MPI_COMM_WORLD) ;
      global.sourceFaces = size_t(globalCounts[0]) ;
      global.targetFaces = size_t(globalCounts[1]) ;
      global.unsupportedRootPlans = size_t(globalCounts[2]) ;
      global.invalidIdentities = size_t(globalCounts[3]) ;
      global.invalidPolygons = size_t(globalCounts[4]) ;
      global.remap.sourceFaces = global.sourceFaces ;
      global.remap.targetFaces = global.targetFaces ;
      global.remap.contributions = size_t(globalCounts[5]) ;
      global.remap.createdFaces = size_t(globalCounts[6]) ;
      global.remap.removedFaces = size_t(globalCounts[7]) ;
      global.remap.duplicateContributions = size_t(globalCounts[8]) ;
      global.remap.missingSourceFaces = size_t(globalCounts[9]) ;
      global.remap.missingTargetFaces = size_t(globalCounts[10]) ;

      unsigned long long localMoreCounts[3] = {
            static_cast<unsigned long long>(local.remap.invalidGeometry),
            static_cast<unsigned long long>(
                  local.remap.inconsistentSourceMoments),
            static_cast<unsigned long long>(
                  local.remap.inconsistentTargetMoments)} ;
      unsigned long long globalMoreCounts[3] = {0} ;
      MPI_Allreduce(localMoreCounts, globalMoreCounts, 3,
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD) ;
      global.remap.invalidGeometry = size_t(globalMoreCounts[0]) ;
      global.remap.inconsistentSourceMoments = size_t(globalMoreCounts[1]) ;
      global.remap.inconsistentTargetMoments = size_t(globalMoreCounts[2]) ;

      double localErrors[4] = {local.remap.maximumSourceAreaError,
            local.remap.maximumTargetAreaError,
            local.remap.maximumSourceCentroidError,
            local.remap.maximumTargetCentroidError} ;
      double globalErrors[4] = {0.0, 0.0, 0.0, 0.0} ;
      MPI_Allreduce(
            localErrors, globalErrors, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) ;
      global.remap.maximumSourceAreaError = globalErrors[0] ;
      global.remap.maximumTargetAreaError = globalErrors[1] ;
      global.remap.maximumSourceCentroidError = globalErrors[2] ;
      global.remap.maximumTargetCentroidError = globalErrors[3] ;
      global.remap.valid = global.valid ;
    }

    bool collectOriginalFaceState(fact_db& facts, CPTR<FaceState>& state,
          FaceTransitionReport& report) {
      storeRepP facesRep = facts.get_variable("faces") ;
      storeRepP cellsRep = facts.get_variable("geom_cells") ;
      storeRepP hexCellsRep = facts.get_variable("hexcells") ;
      storeRepP prismCellsRep = facts.get_variable("prisms") ;
      storeRepP generalCellsRep = facts.get_variable("gnrlcells") ;
      storeRepP faceToNodeRep = facts.get_variable("face2node") ;
      storeRepP positionsRep = facts.get_variable("pos") ;
      storeRepP faceNumberRep = facts.get_variable("fileNumber(face2node)") ;
      storeRepP cellNumberRep = facts.get_variable("fileNumber(geom_cells)") ;
      int localFactsValid =
            facesRep != 0 && cellsRep != 0 && hexCellsRep != 0 &&
                        prismCellsRep != 0 && generalCellsRep != 0 &&
                        faceToNodeRep != 0 && positionsRep != 0 &&
                        faceNumberRep != 0 && cellNumberRep != 0
                  ? 1
                  : 0 ;
      int globalFactsValid = 0 ;
      MPI_Allreduce(&localFactsValid, &globalFactsValid, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalFactsValid == 0) {
        report = FaceTransitionReport() ;
        report.status = face_transition_status::invalid_identity ;
        report.invalidIdentities++ ;
        return false ;
      }

      constraint faces ;
      constraint cells ;
      constraint hexCells ;
      constraint prismCells ;
      constraint generalCells ;
      faces = facesRep ;
      cells = cellsRep ;
      hexCells = hexCellsRep ;
      prismCells = prismCellsRep ;
      generalCells = generalCellsRep ;
      const_multiMap faceToNode(faceToNodeRep) ;
      const_store<vector3d<double>> positions(positionsRep) ;
      const_store<int> faceNumber(faceNumberRep) ;
      const_store<int> cellNumber(cellNumberRep) ;
      int localFirstFace = std::numeric_limits<int>::max() ;
      int localFirstCell = std::numeric_limits<int>::max() ;
      FORALL(*faces, face) {
        localFirstFace = std::min(localFirstFace, faceNumber[face]) ;
      }
      ENDFORALL ;
      FORALL(*cells, cell) {
        localFirstCell = std::min(localFirstCell, cellNumber[cell]) ;
      }
      ENDFORALL ;
      int firstFace = localFirstFace ;
      int firstCell = localFirstCell ;
      MPI_Allreduce(
            &localFirstFace, &firstFace, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      MPI_Allreduce(
            &localFirstCell, &firstCell, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      entitySet ownedNodes = positions.domain() ;
      if (MPI_processes > 1) {
        const size_t nodeKeySpace = positionsRep->getDomainKeySpace() ;
        const std::vector<entitySet> nodePartition =
              facts.get_init_ptn(nodeKeySpace) ;
        const int localPartitionValid =
              nodePartition.size() == size_t(MPI_processes) ? 1 : 0 ;
        int globalPartitionValid = 0 ;
        MPI_Allreduce(&localPartitionValid, &globalPartitionValid, 1, MPI_INT,
              MPI_MIN, MPI_COMM_WORLD) ;
        if (globalPartitionValid == 0) {
          report = FaceTransitionReport() ;
          report.status = face_transition_status::invalid_geometry ;
          report.invalidPolygons++ ;
          return false ;
        }
        ownedNodes &= nodePartition[MPI_rank] ;
      }
      dstore<vector3d<double>> expandedPositions ;
      const bool positionsValid = expandPositionsForFaces(
            positions, faceToNode, ownedNodes, expandedPositions) ;
      int localPositionsValid = positionsValid ? 1 : 0 ;
      int globalPositionsValid = 0 ;
      MPI_Allreduce(&localPositionsValid, &globalPositionsValid, 1, MPI_INT,
            MPI_MIN, MPI_COMM_WORLD) ;
      if (globalPositionsValid == 0) {
        report = FaceTransitionReport() ;
        report.status = face_transition_status::invalid_geometry ;
        report.invalidPolygons++ ;
        return false ;
      }
      std::vector<FaceIdentity> identities ;
      std::vector<std::vector<vector3d<double>>> polygons ;
      FORALL(*faces, face) {
        const int fileFace = faceNumber[face] - firstFace ;
        identities.push_back(FaceIdentity(fileFace,
              FaceKey(face_origin::base_face, fileFace, std::vector<int>()))) ;
        polygons.push_back(facePolygon(face, faceToNode, expandedPositions)) ;
      }
      ENDFORALL ;

      std::vector<RootCellState> roots ;
      bool localTopologyValid = true ;
      FORALL(*cells, cell) {
        const int root = cellNumber[cell] - firstCell ;
        cell_topology::value topology = cell_topology::hex ;
        if (!rootCellTopology(
                  cell, *hexCells, *prismCells, *generalCells, topology)) {
          localTopologyValid = false ;
          continue ;
        }
        roots.push_back(
              RootCellState(root, std::vector<std::vector<int>>(1), topology)) ;
      }
      ENDFORALL ;
      if (positionsValid && localTopologyValid)
        state = FaceState::create(identities, polygons, roots, report) ;
      const int localStateValid =
            positionsValid && localTopologyValid &&
                        state != static_cast<FaceState*>(0)
                  ? 1
                  : 0 ;
      int globalStateValid = 0 ;
      MPI_Allreduce(&localStateValid, &globalStateValid, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalStateValid == 0) {
        state = CPTR<FaceState>() ;
        if (report.status == face_transition_status::available)
          report.status = face_transition_status::invalid_geometry ;
        report.valid = false ;
        return false ;
      }
      return true ;
    }

    bool collectAcceptedFaceState(int numNodes, int numFaces,
          const store<FineFaces>& fineFacesCell,
          const store<FineFaces>& fineFaces, const FaceSources& faceSources,
          const store<vector3d<double>>& positions, const multiMap& faceToNode,
          const std::vector<entitySet>& outputCellPartition, fact_db& facts,
          CPTR<FaceState>& state, store<FaceId>& faceIds,
          store<CellId>& cellIds, FaceTransitionReport& report) {
      storeRepP facesRep = facts.get_variable("faces") ;
      storeRepP cellsRep = facts.get_variable("geom_cells") ;
      storeRepP hexCellsRep = facts.get_variable("hexcells") ;
      storeRepP prismCellsRep = facts.get_variable("prisms") ;
      storeRepP generalCellsRep = facts.get_variable("gnrlcells") ;
      storeRepP faceNumberRep = facts.get_variable("fileNumber(face2node)") ;
      storeRepP rootNumberRep = facts.get_variable("planRootFileNumber") ;
      storeRepP cellOffsetRep = facts.get_variable("balanced_cell_offset") ;
      storeRepP cellPathsRep = facts.get_variable("cellLeafPaths") ;
      storeRepP facePathsRep = facts.get_variable("faceLeafPaths") ;
      int localFactsValid =
            facesRep != 0 && cellsRep != 0 && hexCellsRep != 0 &&
                        prismCellsRep != 0 && generalCellsRep != 0 &&
                        faceNumberRep != 0 && rootNumberRep != 0 &&
                        cellOffsetRep != 0 && cellPathsRep != 0 &&
                        facePathsRep != 0
                  ? 1
                  : 0 ;
      int globalFactsValid = 0 ;
      MPI_Allreduce(&localFactsValid, &globalFactsValid, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalFactsValid == 0) {
        if (MPI_rank == 0) {
          std::cerr << "Unable to collect accepted face state; missing facts:" ;
          if (facesRep == 0)
            std::cerr << " faces" ;
          if (cellsRep == 0)
            std::cerr << " geom_cells" ;
          if (hexCellsRep == 0)
            std::cerr << " hexcells" ;
          if (prismCellsRep == 0)
            std::cerr << " prisms" ;
          if (generalCellsRep == 0)
            std::cerr << " gnrlcells" ;
          if (faceNumberRep == 0)
            std::cerr << " fileNumber(face2node)" ;
          if (rootNumberRep == 0)
            std::cerr << " planRootFileNumber" ;
          if (cellOffsetRep == 0)
            std::cerr << " balanced_cell_offset" ;
          if (cellPathsRep == 0)
            std::cerr << " cellLeafPaths" ;
          if (facePathsRep == 0)
            std::cerr << " faceLeafPaths" ;
          std::cerr << std::endl ;
        }
        report = FaceTransitionReport() ;
        report.status = face_transition_status::invalid_identity ;
        report.invalidIdentities++ ;
        return false ;
      }

      constraint faces ;
      constraint cells ;
      constraint hexCells ;
      constraint prismCells ;
      constraint generalCells ;
      faces = facesRep ;
      cells = cellsRep ;
      hexCells = hexCellsRep ;
      prismCells = prismCellsRep ;
      generalCells = generalCellsRep ;
      const_store<int> faceNumber(faceNumberRep) ;
      const_store<int> rootNumber(rootNumberRep) ;
      const_store<int> cellOffset(cellOffsetRep) ;
      const_store<std::vector<int>> encodedCellPaths(cellPathsRep) ;
      const_store<std::vector<int>> encodedFacePaths(facePathsRep) ;
      int localFirstFace = std::numeric_limits<int>::max() ;
      int localFirstRoot = std::numeric_limits<int>::max() ;
      FORALL(*faces, face) {
        localFirstFace = std::min(localFirstFace, faceNumber[face]) ;
      }
      ENDFORALL ;
      FORALL(*cells, cell) {
        localFirstRoot = std::min(localFirstRoot, rootNumber[cell]) ;
      }
      ENDFORALL ;
      int firstFace = localFirstFace ;
      int firstRoot = localFirstRoot ;
      MPI_Allreduce(
            &localFirstFace, &firstFace, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      MPI_Allreduce(
            &localFirstRoot, &firstRoot, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      dstore<vector3d<double>> expandedPositions ;
      const bool positionsValid = expandPositionsForFaces(
            positions, faceToNode, positions.domain(), expandedPositions) ;
      int localPositionsValid = positionsValid ? 1 : 0 ;
      int globalPositionsValid = 0 ;
      MPI_Allreduce(&localPositionsValid, &globalPositionsValid, 1, MPI_INT,
            MPI_MIN, MPI_COMM_WORLD) ;
      if (globalPositionsValid == 0) {
        report = FaceTransitionReport() ;
        report.status = face_transition_status::invalid_geometry ;
        report.invalidPolygons++ ;
        return false ;
      }
      bool localValid = true ;

      std::map<int, std::vector<std::vector<int>>> cellPaths ;
      std::vector<RootCellState> roots ;
      std::vector<CellIdentity> cellIdentities ;
      entitySet outputCells ;
      FORALL(*cells, cell) {
        std::vector<std::vector<int>> paths ;
        if (!decodeLeafPaths(encodedCellPaths[cell], paths)) {
          localValid = false ;
          continue ;
        }
        cellPaths[cell] = paths ;
        const int generatedCellBase =
              MPI_processes == 1 ? numNodes + numFaces : numNodes ;
        const int firstCell = generatedCellBase + cellOffset[cell] ;
        const int root = rootNumber[cell] - firstRoot ;
        cell_topology::value topology = cell_topology::hex ;
        if (!rootCellTopology(
                  cell, *hexCells, *prismCells, *generalCells, topology)) {
          localValid = false ;
          continue ;
        }
        roots.push_back(RootCellState(root, paths, topology)) ;
        for (size_t leaf = 0; leaf < paths.size(); ++leaf) {
          const int outputCell = firstCell + int(leaf) ;
          outputCells += outputCell ;
          cellIdentities.push_back(
                CellIdentity(outputCell, CellKey(root, paths[leaf]))) ;
        }
      }
      ENDFORALL ;

      std::vector<FaceIdentity> identities ;
      std::vector<std::vector<vector3d<double>>> polygons ;
      entitySet outputFaces = faceToNode.domain() ;
      const entitySet cellFaces = faceSources.cell.domain() ;
      const entitySet originalFaces = faceSources.face.domain() ;
      localValid = localValid && (cellFaces + originalFaces) == outputFaces &&
                   (cellFaces & originalFaces) == EMPTY &&
                   faceSources.ordinal.domain() == outputFaces ;
      std::map<int, std::vector<std::vector<int>>> facePaths ;
      // Assembly supplies this association explicitly. Output face numbering
      // and iteration order need not match the original cells or faces.
      const entitySet boundFaces = outputFaces & (cellFaces + originalFaces) &
                                   faceSources.ordinal.domain() ;
      FORALL(boundFaces, outputFace) {
        const bool cellInterior = cellFaces.inSet(outputFace) ;
        const int entity = cellInterior ? faceSources.cell[outputFace]
                                        : faceSources.face[outputFace] ;
        const int ordinal = faceSources.ordinal[outputFace] ;
        if (ordinal < 0) {
          localValid = false ;
          continue ;
        }
        if (cellInterior && (*cells).inSet(entity) &&
              fineFacesCell.domain().inSet(entity)) {
          const std::vector<std::vector<int>>& paths = cellPaths[entity] ;
          if (size_t(ordinal) >= fineFacesCell[entity].size() ||
                fineFacesCell[entity][ordinal].size() < 5) {
            localValid = false ;
            continue ;
          }
          const int left =
                fineFacesCell[entity][ordinal][0] - cellOffset[entity] - 1 ;
          const int right =
                fineFacesCell[entity][ordinal][1] - cellOffset[entity] - 1 ;
          if (left < 0 || right < 0 || size_t(left) >= paths.size() ||
                size_t(right) >= paths.size()) {
            localValid = false ;
            continue ;
          }
          identities.push_back(FaceIdentity(
                outputFace, interiorKey(rootNumber[entity] - firstRoot,
                                  paths[left], paths[right]))) ;
        } else if (!cellInterior && (*faces).inSet(entity) &&
                   fineFaces.domain().inSet(entity)) {
          std::vector<std::vector<int>>& paths = facePaths[entity] ;
          if (paths.empty() &&
                !decodeLeafPaths(encodedFacePaths[entity], paths)) {
            localValid = false ;
            continue ;
          }
          if (paths.size() != fineFaces[entity].size() ||
                size_t(ordinal) >= paths.size()) {
            localValid = false ;
            continue ;
          }
          identities.push_back(FaceIdentity(outputFace,
                FaceKey(face_origin::base_face, faceNumber[entity] - firstFace,
                      paths[ordinal]))) ;
        } else {
          localValid = false ;
          continue ;
        }
        polygons.push_back(
              facePolygon(outputFace, faceToNode, expandedPositions)) ;
      }
      ENDFORALL ;

      if (localValid)
        state = FaceState::create(identities, polygons, roots, report) ;
      localValid = localValid && state != static_cast<FaceState*>(0) ;
      int localStateValid = localValid ? 1 : 0 ;
      int globalStateValid = 0 ;
      MPI_Allreduce(&localStateValid, &globalStateValid, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalStateValid == 0) {
        state = CPTR<FaceState>() ;
        report.valid = false ;
        if (report.status == face_transition_status::available)
          report.status = face_transition_status::invalid_geometry ;
        return false ;
      }
      faceIds.allocate(outputFaces) ;
      for (size_t face = 0; face < identities.size(); ++face)
        faceIds[identities[face].face] = identities[face].id ;
      size_t invalidCellIdentities = 0 ;
      if (!validateCellIdentityHashesDistributed(
                cellIdentities, invalidCellIdentities)) {
        report.valid = false ;
        report.status = face_transition_status::invalid_identity ;
        report.invalidIdentities += invalidCellIdentities ;
        state = CPTR<FaceState>() ;
        return false ;
      }
      const bool localPartitionValid =
            outputCellPartition.size() == size_t(MPI_processes) ;
      int localPartition = localPartitionValid ? 1 : 0 ;
      int globalPartition = 0 ;
      MPI_Allreduce(&localPartition, &globalPartition, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalPartition == 0) {
        report.valid = false ;
        report.status = face_transition_status::invalid_identity ;
        report.invalidIdentities++ ;
        state = CPTR<FaceState>() ;
        return false ;
      }
      dstore<CellId> availableCellIds ;
      for (size_t cell = 0; cell < cellIdentities.size(); ++cell)
        availableCellIds[cellIdentities[cell].cell] = cellIdentities[cell].id ;
      entitySet ownedOutputCells = outputCells ;
      std::vector<entitySet> initialCellPartition =
            all_collect_vectors(ownedOutputCells, MPI_COMM_WORLD) ;
      entitySet requiredCells = outputCellPartition[MPI_rank] ;
      storeRepP availableRep = availableCellIds.Rep() ;
      fill_clone(availableRep, requiredCells, initialCellPartition) ;
      availableCellIds.setRep(availableRep) ;
      const bool localCellIdsValid =
            (requiredCells - availableCellIds.domain()).size() == 0 ;
      int localIdsValid = localCellIdsValid ? 1 : 0 ;
      int globalIdsValid = 0 ;
      MPI_Allreduce(&localIdsValid, &globalIdsValid, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalIdsValid == 0) {
        report.valid = false ;
        report.status = face_transition_status::invalid_identity ;
        report.invalidIdentities++ ;
        state = CPTR<FaceState>() ;
        return false ;
      }
      cellIds.allocate(requiredCells) ;
      FORALL(requiredCells, cell) {
        cellIds[cell] = availableCellIds[cell] ;
      }
      ENDFORALL ;
      return true ;
    }
  }

  void FaceState::retainNodes(const store<NodeId>& ids,
      const store<vector3d<double>>& positions,
      const std::vector<entitySet>& partition) {
    nodeIds.allocate(ids.domain()) ;
    nodePositions.allocate(positions.domain()) ;
    FORALL(ids.domain(), node) {
      nodeIds[node] = ids[node] ;
      nodePositions[node] = positions[node] ;
    } ENDFORALL ;
    nodePartition = partition ;
  }

  bool installBaseMeshIds(fact_db& facts) {
    storeRepP facesRep = facts.get_variable("faces") ;
    storeRepP cellsRep = facts.get_variable("geom_cells") ;
    storeRepP positionsRep = facts.get_variable("pos") ;
    storeRepP faceNumberRep = facts.get_variable("fileNumber(face2node)") ;
    storeRepP cellNumberRep = facts.get_variable("fileNumber(geom_cells)") ;
    storeRepP nodeNumberRep = facts.get_variable("fileNumber(pos)") ;
    const int localPrerequisites = facesRep != 0 && cellsRep != 0 &&
                                               positionsRep != 0 &&
                                               nodeNumberRep != 0
                                         ? 1
                                         : 0 ;
    int globalPrerequisites = 0 ;
    MPI_Allreduce(&localPrerequisites, &globalPrerequisites, 1, MPI_INT,
          MPI_MIN, MPI_COMM_WORLD) ;
    if (globalPrerequisites == 0) {
      if (localPrerequisites == 0)
        cerr << "rank " << MPI_rank
             << " cannot install base AMR identities: missing faces, "
                "geom_cells, pos, or fileNumber(pos)"
             << endl ;
      return false ;
    }
    constraint faces ;
    constraint cells ;
    faces = facesRep ;
    cells = cellsRep ;
    const_store<vector3d<double>> positions(positionsRep) ;
    const_store<int> faceNumber ;
    const_store<int> cellNumber ;
    const_store<int> nodeNumber(nodeNumberRep) ;
    if (faceNumberRep != 0)
      faceNumber = faceNumberRep ;
    if (cellNumberRep != 0)
      cellNumber = cellNumberRep ;
    fact_db::distribute_infoP distribution = facts.get_distribute_info() ;
    Map localToFile ;
    dMap faceGlobalToFile ;
    dMap cellGlobalToFile ;
    if (distribution != 0)
      localToFile = distribution->l2f.Rep() ;
    const size_t faceKeySpace = facesRep->getDomainKeySpace() ;
    const size_t cellKeySpace = cellsRep->getDomainKeySpace() ;
    const size_t nodeKeySpace = positionsRep->getDomainKeySpace() ;
    if (distribution != 0 && faceKeySpace < distribution->g2fv.size())
      faceGlobalToFile = distribution->g2fv[faceKeySpace].Rep() ;
    if (distribution != 0 && cellKeySpace < distribution->g2fv.size())
      cellGlobalToFile = distribution->g2fv[cellKeySpace].Rep() ;
    const bool faceFactCovers =
          faceNumberRep != 0 && ((*faces) - faceNumber.domain()).size() == 0 ;
    const bool cellFactCovers =
          cellNumberRep != 0 && ((*cells) - cellNumber.domain()).size() == 0 ;
    const bool faceGlobalMapCovers =
          distribution != 0 &&
          ((*faces) - faceGlobalToFile.domain()).size() == 0 ;
    const bool cellGlobalMapCovers =
          distribution != 0 &&
          ((*cells) - cellGlobalToFile.domain()).size() == 0 ;
    const bool faceMapCovers =
          distribution != 0 && ((*faces) - localToFile.domain()).size() == 0 ;
    const bool cellMapCovers =
          distribution != 0 && ((*cells) - localToFile.domain()).size() == 0 ;
    const bool nodeFactCovers =
          (positions.domain() - nodeNumber.domain()).size() == 0 ;
    const bool localSourcesValid = (faceFactCovers || faceGlobalMapCovers ||
                                         faceMapCovers || distribution == 0) &&
                                   (cellFactCovers || cellGlobalMapCovers ||
                                         cellMapCovers || distribution == 0) &&
                                   nodeFactCovers ;
    int localValidSource = localSourcesValid ? 1 : 0 ;
    int globalValidSource = 0 ;
    MPI_Allreduce(&localValidSource, &globalValidSource, 1, MPI_INT, MPI_MIN,
          MPI_COMM_WORLD) ;
    if (globalValidSource == 0) {
      cerr << "rank " << MPI_rank
           << " cannot install base AMR identities: face domain=" << *faces
           << " cell domain=" << *cells
           << " fileNumber face domain=" << faceNumber.domain()
           << " fileNumber cell domain=" << cellNumber.domain()
           << " node domain=" << positions.domain()
           << " fileNumber node domain=" << nodeNumber.domain()
           << " face g2f domain=" << faceGlobalToFile.domain()
           << " cell g2f domain=" << cellGlobalToFile.domain()
           << " l2f domain=" << localToFile.domain() << endl ;
      return false ;
    }

    int localFirstFace = std::numeric_limits<int>::max() ;
    int localFirstCell = std::numeric_limits<int>::max() ;
    FORALL(*faces, face) {
      const int file =
            faceFactCovers
                  ? faceNumber[face]
                  : (faceGlobalMapCovers
                                ? faceGlobalToFile[face]
                                : (faceMapCovers ? localToFile[face] : face)) ;
      localFirstFace = std::min(localFirstFace, file) ;
    }
    ENDFORALL ;
    FORALL(*cells, cell) {
      const int file =
            cellFactCovers
                  ? cellNumber[cell]
                  : (cellGlobalMapCovers
                                ? cellGlobalToFile[cell]
                                : (cellMapCovers ? localToFile[cell] : cell)) ;
      localFirstCell = std::min(localFirstCell, file) ;
    }
    ENDFORALL ;
    int firstFace = localFirstFace ;
    int firstCell = localFirstCell ;
    MPI_Allreduce(
          &localFirstFace, &firstFace, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    MPI_Allreduce(
          &localFirstCell, &firstCell, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    if (firstFace == std::numeric_limits<int>::max() ||
          firstCell == std::numeric_limits<int>::max()) {
      cerr << "rank " << MPI_rank
           << " cannot install base AMR identities: empty global face or "
              "cell domain"
           << endl ;
      return false ;
    }

    entitySet ownedFaces = *faces ;
    entitySet ownedCells = *cells ;
    entitySet ownedNodes = positions.domain() ;
    if (MPI_processes > 1) {
      const std::vector<entitySet> facePartition =
            facts.get_init_ptn(faceKeySpace) ;
      const std::vector<entitySet> cellPartition =
            facts.get_init_ptn(cellKeySpace) ;
      const std::vector<entitySet> nodePartition =
            facts.get_init_ptn(nodeKeySpace) ;
      const int localPartitions =
            facePartition.size() == size_t(MPI_processes) &&
                        cellPartition.size() == size_t(MPI_processes) &&
                        nodePartition.size() == size_t(MPI_processes)
                  ? 1
                  : 0 ;
      int globalPartitions = 0 ;
      MPI_Allreduce(&localPartitions, &globalPartitions, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) ;
      if (globalPartitions == 0) {
        if (localPartitions == 0)
          cerr << "rank " << MPI_rank
               << " cannot install base AMR identities: missing key-space "
                  "ownership partitions"
               << endl ;
        return false ;
      }
      ownedFaces &= facePartition[MPI_rank] ;
      ownedCells &= cellPartition[MPI_rank] ;
      ownedNodes &= nodePartition[MPI_rank] ;
    }
    store<FaceId> faceIds ;
    store<CellId> cellIds ;
    store<NodeId> nodeIds ;
    faceIds.allocate(*faces) ;
    cellIds.allocate(*cells) ;
    nodeIds.allocate(positions.domain()) ;
    std::vector<FaceIdentity> faceIdentities ;
    std::vector<CellIdentity> cellIdentities ;
    std::vector<std::pair<NodeId, long long>> nodeIdentities ;
    faceIdentities.reserve((*faces).size()) ;
    cellIdentities.reserve((*cells).size()) ;
    nodeIdentities.reserve(ownedNodes.size()) ;
    FORALL(*faces, face) {
      const int file =
            faceFactCovers
                  ? faceNumber[face]
                  : (faceGlobalMapCovers
                                ? faceGlobalToFile[face]
                                : (faceMapCovers ? localToFile[face] : face)) ;
      const FaceKey key(
            face_origin::base_face, file - firstFace, std::vector<int>()) ;
      const FaceIdentity identity(face, key) ;
      if (ownedFaces.inSet(face))
        faceIdentities.push_back(identity) ;
      faceIds[face] = identity.id ;
    }
    ENDFORALL ;
    FORALL(*cells, cell) {
      const int file =
            cellFactCovers
                  ? cellNumber[cell]
                  : (cellGlobalMapCovers
                                ? cellGlobalToFile[cell]
                                : (cellMapCovers ? localToFile[cell] : cell)) ;
      const CellIdentity identity(
            cell, CellKey(file - firstCell, std::vector<int>())) ;
      if (ownedCells.inSet(cell))
        cellIdentities.push_back(identity) ;
      cellIds[cell] = identity.id ;
    }
    ENDFORALL ;
    FORALL(positions.domain(), node) {
      const long long file = nodeNumber[node] ;
      nodeIds[node] = persistentBaseNodeId(file) ;
      if (ownedNodes.inSet(node))
        nodeIdentities.push_back(std::make_pair(nodeIds[node], file)) ;
    }
    ENDFORALL ;
    size_t invalidFaceIdentities = 0 ;
    size_t invalidCellIdentities = 0 ;
    size_t invalidNodeIdentities = 0 ;
    if (!detail::validateFaceIdentityHashesDistributed(
              faceIdentities, invalidFaceIdentities) ||
          !detail::validateCellIdentityHashesDistributed(
                cellIdentities, invalidCellIdentities) ||
          !detail::validateNodeIdentityHashesDistributed(
                nodeIdentities, invalidNodeIdentities)) {
      cerr << "rank " << MPI_rank
           << " cannot install base AMR identities: " << invalidFaceIdentities
           << " invalid face identities and " << invalidCellIdentities
           << " invalid cell identities and " << invalidNodeIdentities
           << " invalid node identities" << endl ;
      return false ;
    }

    storeRepP existingFaceRep = facts.get_variable("faceId") ;
    storeRepP existingCellRep = facts.get_variable("cellId") ;
    storeRepP existingNodeRep = facts.get_variable("nodeId") ;
    const_store<FaceId> existingFaceIds ;
    const_store<CellId> existingCellIds ;
    const_store<NodeId> existingNodeIds ;
    if (existingFaceRep != 0)
      existingFaceIds = existingFaceRep ;
    if (existingCellRep != 0)
      existingCellIds = existingCellRep ;
    if (existingNodeRep != 0)
      existingNodeIds = existingNodeRep ;

    const bool localFaceComplete =
          existingFaceRep != 0 && existingFaceIds.domain() == faceIds.domain() ;
    const bool localCellComplete =
          existingCellRep != 0 && existingCellIds.domain() == cellIds.domain() ;
    const bool localNodeComplete =
          existingNodeRep != 0 && existingNodeIds.domain() == nodeIds.domain() ;
    const bool localFaceHasData =
          existingFaceRep != 0 && existingFaceIds.domain().size() != 0 ;
    const bool localCellHasData =
          existingCellRep != 0 && existingCellIds.domain().size() != 0 ;
    const bool localNodeHasData =
          existingNodeRep != 0 && existingNodeIds.domain().size() != 0 ;
    int localExisting[6] = {localFaceComplete ? 1 : 0,
          localCellComplete ? 1 : 0, localNodeComplete ? 1 : 0,
          localFaceHasData ? 1 : 0, localCellHasData ? 1 : 0,
          localNodeHasData ? 1 : 0} ;
    int allComplete[3] = {0, 0, 0} ;
    int anyData[3] = {0, 0, 0} ;
    MPI_Allreduce(
          localExisting, allComplete, 3, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    MPI_Allreduce(
          localExisting + 3, anyData, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

    // A rule declaration may leave a non-null store with no allocated
    // domain.  Populate that declaration, but do not overwrite a partially
    // installed distributed identity relation.
    if ((allComplete[0] == 0 && anyData[0] != 0) ||
          (allComplete[1] == 0 && anyData[1] != 0) ||
          (allComplete[2] == 0 && anyData[2] != 0)) {
      cerr << "rank " << MPI_rank
           << " cannot install base AMR identities over partially populated "
              "facts: face domain="
           << existingFaceIds.domain() << " expected=" << faceIds.domain()
           << " cell domain=" << existingCellIds.domain()
           << " expected=" << cellIds.domain()
           << " node domain=" << existingNodeIds.domain()
           << " expected=" << nodeIds.domain() << endl ;
      return false ;
    }

    bool localValid = true ;
    if (allComplete[0] != 0) {
      FORALL(faceIds.domain(), face) {
        if (existingFaceIds[face] != faceIds[face])
          localValid = false ;
      }
      ENDFORALL ;
    }
    if (allComplete[1] != 0) {
      FORALL(cellIds.domain(), cell) {
        if (existingCellIds[cell] != cellIds[cell])
          localValid = false ;
      }
      ENDFORALL ;
    }
    if (allComplete[2] != 0) {
      FORALL(nodeIds.domain(), node) {
        if (existingNodeIds[node] != nodeIds[node])
          localValid = false ;
      }
      ENDFORALL ;
    }
    int local = localValid ? 1 : 0 ;
    int global = 0 ;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
    if (global == 0) {
      cerr << "rank " << MPI_rank
           << " cannot install base AMR identities: existing values differ "
              "from the canonical base-grid identities"
           << endl ;
      return false ;
    }
    if (allComplete[0] == 0)
      facts.create_fact("faceId", faceIds) ;
    if (allComplete[1] == 0)
      facts.create_fact("cellId", cellIds) ;
    if (allComplete[2] == 0)
      facts.create_fact("nodeId", nodeIds) ;
    return true ;
  }

  namespace detail {
    struct PersistentCellOwner {
      CellId cell ;
      int owner ;
    } ;

    struct RemovedFaceWithGeometry {
      RemovedFace removed ;
      FaceGeometry source ;
    } ;

    CPTR<FaceRemap> redistributeFaceRemap(const CPTR<FaceRemap>& remap,
          const FaceTransitionReport& globalReport,
          const vector<entitySet>& facePtn, const vector<entitySet>& cellPtn,
          const store<FaceId>& generatedFaceIds,
          const store<CellId>& generatedCellIds,
          FaceTransitionReport& localReport) {
      localReport = globalReport ;
      if (remap == static_cast<FaceRemap*>(0))
        return CPTR<FaceRemap>() ;

      std::map<int, int> faceOwner ;
      for (size_t process = 0; process < facePtn.size(); ++process) {
        FORALL(facePtn[process], face) {
          if (!faceOwner.insert(std::make_pair(face, int(process))).second)
            faceOwner[face] = -1 ;
        }
        ENDFORALL ;
      }
      std::map<FaceId, int> targetOwner ;
      int localMissingOwner = 0 ;
      FORALL(generatedFaceIds.domain(), face) {
        const std::map<int, int>::const_iterator owner = faceOwner.find(face) ;
        if (owner == faceOwner.end() || owner->second < 0)
          localMissingOwner++ ;
        else
          targetOwner[generatedFaceIds[face]] = owner->second ;
      }
      ENDFORALL ;
      int globalMissingOwner = 0 ;
      MPI_Allreduce(&localMissingOwner, &globalMissingOwner, 1, MPI_INT,
            MPI_SUM, MPI_COMM_WORLD) ;
      if (globalMissingOwner != 0) {
        cerr << "rank " << MPI_rank << " transition has " << globalMissingOwner
             << " target faces without an owner" << endl ;
        localReport.valid = false ;
        localReport.status = face_transition_status::invalid_identity ;
        localReport.invalidIdentities += size_t(globalMissingOwner) ;
        return CPTR<FaceRemap>() ;
      }

      std::map<int, int> cellOwner ;
      for (size_t process = 0; process < cellPtn.size(); ++process) {
        FORALL(cellPtn[process], cell) {
          if (!cellOwner.insert(std::make_pair(cell, int(process))).second)
            cellOwner[cell] = -1 ;
        }
        ENDFORALL ;
      }
      std::vector<std::vector<PersistentCellOwner>> outgoingCellOwners(
            MPI_processes) ;
      int localInvalidCellOwners = 0 ;
      FORALL(generatedCellIds.domain(), cell) {
        const std::map<int, int>::const_iterator owner = cellOwner.find(cell) ;
        if (owner == cellOwner.end() || owner->second < 0) {
          localInvalidCellOwners++ ;
          continue ;
        }
        const CellId id = generatedCellIds[cell] ;
        const int lookupOwner = int(
              static_cast<std::uint64_t>(id) % std::uint64_t(MPI_processes)) ;
        PersistentCellOwner entry = {id, owner->second} ;
        outgoingCellOwners[lookupOwner].push_back(entry) ;
      }
      ENDFORALL ;
      const std::vector<PersistentCellOwner> localCellOwners =
            exchangeByDestination(outgoingCellOwners) ;
      std::map<CellId, int> targetCellOwner ;
      for (size_t entry = 0; entry < localCellOwners.size(); ++entry) {
        const std::map<CellId, int>::const_iterator existing =
              targetCellOwner.find(localCellOwners[entry].cell) ;
        if (existing != targetCellOwner.end()) {
          localInvalidCellOwners++ ;
          continue ;
        }
        targetCellOwner[localCellOwners[entry].cell] =
              localCellOwners[entry].owner ;
      }
      int globalInvalidCellOwners = 0 ;
      MPI_Allreduce(&localInvalidCellOwners, &globalInvalidCellOwners, 1,
            MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if (globalInvalidCellOwners != 0) {
        cerr << "rank " << MPI_rank << " transition has "
             << globalInvalidCellOwners
             << " ranks with an invalid target-cell ownership map" << endl ;
        localReport.valid = false ;
        localReport.status = face_transition_status::invalid_identity ;
        localReport.invalidIdentities += size_t(globalInvalidCellOwners) ;
        return CPTR<FaceRemap>() ;
      }

      std::map<FaceId, FaceGeometry> sourceById ;
      const std::vector<FaceGeometry>& sourceGeometry =
            remap->sourceFaceGeometry() ;
      for (size_t face = 0; face < sourceGeometry.size(); ++face)
        sourceById[sourceGeometry[face].face] = sourceGeometry[face] ;

      std::vector<std::vector<FaceGeometry>> outgoingSources(MPI_processes) ;
      std::vector<std::set<FaceId>> outgoingSourceIds(MPI_processes) ;
      std::vector<std::vector<FaceGeometry>> outgoingTargets(MPI_processes) ;
      std::vector<std::vector<FaceOverlap>> outgoingContributions(
            MPI_processes) ;
      std::vector<std::vector<CreatedFace>> outgoingCreated(MPI_processes) ;
      std::vector<std::vector<RemovedFace>> outgoingRemoved(MPI_processes) ;

      const std::vector<FaceGeometry>& targetGeometry =
            remap->targetFaceGeometry() ;
      for (size_t face = 0; face < targetGeometry.size(); ++face) {
        const std::map<FaceId, int>::const_iterator owner =
              targetOwner.find(targetGeometry[face].face) ;
        if (owner == targetOwner.end()) {
          localMissingOwner++ ;
          continue ;
        }
        outgoingTargets[owner->second].push_back(targetGeometry[face]) ;
      }
      const std::vector<FaceOverlap>& contributions = remap->overlaps() ;
      for (size_t entry = 0; entry < contributions.size(); ++entry) {
        const std::map<FaceId, int>::const_iterator owner =
              targetOwner.find(contributions[entry].target) ;
        const std::map<FaceId, FaceGeometry>::const_iterator source =
              sourceById.find(contributions[entry].source) ;
        if (owner == targetOwner.end() || source == sourceById.end()) {
          localMissingOwner++ ;
          continue ;
        }
        outgoingContributions[owner->second].push_back(contributions[entry]) ;
        if (outgoingSourceIds[owner->second].insert(source->first).second)
          outgoingSources[owner->second].push_back(source->second) ;
      }
      const std::vector<CreatedFace>& created = remap->createdFaces() ;
      for (size_t face = 0; face < created.size(); ++face) {
        const std::map<FaceId, int>::const_iterator owner =
              targetOwner.find(created[face].targetFace) ;
        if (owner == targetOwner.end()) {
          localMissingOwner++ ;
          continue ;
        }
        outgoingCreated[owner->second].push_back(created[face]) ;
      }
      std::vector<std::vector<RemovedFaceWithGeometry>> outgoingRemovedLookups(
            MPI_processes) ;
      const std::vector<RemovedFace>& removed = remap->removedFaces() ;
      for (size_t face = 0; face < removed.size(); ++face) {
        const std::map<FaceId, FaceGeometry>::const_iterator source =
              sourceById.find(removed[face].sourceFace) ;
        if (source == sourceById.end()) {
          localMissingOwner++ ;
          continue ;
        }
        const int lookupOwner =
              int(static_cast<std::uint64_t>(removed[face].targetCell) %
                    std::uint64_t(MPI_processes)) ;
        RemovedFaceWithGeometry lookup = {removed[face], source->second} ;
        outgoingRemovedLookups[lookupOwner].push_back(lookup) ;
      }
      const std::vector<RemovedFaceWithGeometry> localRemovedLookups =
            exchangeByDestination(outgoingRemovedLookups) ;
      for (size_t face = 0; face < localRemovedLookups.size(); ++face) {
        const std::map<CellId, int>::const_iterator owner =
              targetCellOwner.find(
                    localRemovedLookups[face].removed.targetCell) ;
        if (owner == targetCellOwner.end()) {
          localMissingOwner++ ;
          continue ;
        }
        outgoingRemoved[owner->second].push_back(
              localRemovedLookups[face].removed) ;
        const FaceGeometry& source = localRemovedLookups[face].source ;
        if (outgoingSourceIds[owner->second].insert(source.face).second)
          outgoingSources[owner->second].push_back(source) ;
      }
      MPI_Allreduce(&localMissingOwner, &globalMissingOwner, 1, MPI_INT,
            MPI_SUM, MPI_COMM_WORLD) ;
      if (globalMissingOwner != 0) {
        cerr << "rank " << MPI_rank << " transition has " << globalMissingOwner
             << " relation records without a target owner/source geometry"
             << endl ;
        localReport.valid = false ;
        localReport.status = face_transition_status::invalid_identity ;
        localReport.invalidIdentities += size_t(globalMissingOwner) ;
        return CPTR<FaceRemap>() ;
      }

      std::vector<FaceGeometry> localSources =
            exchangeByDestination(outgoingSources) ;
      const std::vector<FaceGeometry> localTargets =
            exchangeByDestination(outgoingTargets) ;
      const std::vector<FaceOverlap> localContributions =
            exchangeByDestination(outgoingContributions) ;
      const std::vector<CreatedFace> localCreated =
            exchangeByDestination(outgoingCreated) ;
      const std::vector<RemovedFace> localRemoved =
            exchangeByDestination(outgoingRemoved) ;

      std::sort(localSources.begin(), localSources.end(),
            [](const FaceGeometry& left, const FaceGeometry& right) {
              return left.face < right.face;
            }) ;
      int localConflictingSources = 0 ;
      for (size_t face = 1; face < localSources.size(); ++face) {
        const FaceGeometry& previous = localSources[face - 1] ;
        const FaceGeometry& current = localSources[face] ;
        if (previous.face == current.face &&
              (previous.area != current.area ||
                    previous.centroid.x != current.centroid.x ||
                    previous.centroid.y != current.centroid.y ||
                    previous.centroid.z != current.centroid.z))
          localConflictingSources++ ;
      }
      int globalConflictingSources = 0 ;
      MPI_Allreduce(&localConflictingSources, &globalConflictingSources, 1,
            MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if (globalConflictingSources != 0) {
        localReport.valid = false ;
        localReport.status = face_transition_status::invalid_geometry ;
        localReport.remap.invalidGeometry += size_t(globalConflictingSources) ;
        return CPTR<FaceRemap>() ;
      }
      localSources.erase(
            std::unique(localSources.begin(), localSources.end(),
                  [](const FaceGeometry& left, const FaceGeometry& right) {
                    return left.face == right.face;
                  }),
            localSources.end()) ;

      FaceRemapReport localRemapReport ;
      CPTR<FaceRemap> localRemap = FaceRemap::createTargetOwned(localSources,
            localTargets, localContributions, localCreated, localRemoved,
            localRemapReport) ;
      const bool localRemapValid =
            localRemap != static_cast<FaceRemap*>(0) && localRemapReport.valid ;
      const int localValid = localRemapValid ? 1 : 0 ;
      int globallyValid = 0 ;
      MPI_Allreduce(&localValid, &globallyValid, 1, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD) ;
      if (globallyValid == 0) {
        localReport.valid = false ;
        localReport.status = face_transition_status::inconsistent_relation ;
        localReport.remap = localRemapReport ;
        return CPTR<FaceRemap>() ;
      }
      return localRemap ;
    }
  }
}

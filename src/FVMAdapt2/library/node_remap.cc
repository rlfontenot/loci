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

#include <distribute.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <vector>

namespace Loci {
  namespace detail {
    namespace {

      int directoryOwner(NodeId node) {
        return int(
              static_cast<std::uint64_t>(node) % std::uint64_t(MPI_processes)) ;
      }

      struct NodeLookupRequest {
        NodeId node ;
        int requester ;
      } ;

      struct NodeLookupResponse {
        NodeId node ;
        double position[3] ;
        int found ;
      } ;

      struct CurrentNodeRecord {
        NodeId node ;
        int number ;
      } ;

      struct NodeGeometryRecord {
        NodeId node ;
        double position[3] ;
      } ;

      struct NodeContributionRecord {
        NodeId sourceNode ;
        NodeId targetNode ;
        double weight ;
      } ;

      struct NodeOriginRecord {
        NodeId targetNode ;
        int kind ;
      } ;

      bool finitePosition(const vector3d<double>& position) {
        return std::isfinite(position.x) && std::isfinite(position.y) &&
               std::isfinite(position.z) ;
      }

      bool validConstruction(const FineNodeConstruction& construction,
            size_t& errors, double relativeTolerance) {
        if (construction.node == 0 || construction.parentCount < 0 ||
              construction.parentCount > FineNodeConstruction::maximumParents) {
          ++errors ;
          return false ;
        }
        if (construction.kind == node_construction::base_node) {
          if (construction.baseFileNumber < 0 ||
                construction.parentCount != 0 ||
                construction.node !=
                      persistentBaseNodeId(construction.baseFileNumber)) {
            ++errors ;
            return false ;
          }
          return true ;
        }
        if (construction.kind < node_construction::edge ||
              construction.kind > node_construction::cell ||
              construction.baseFileNumber != -1 ||
              construction.parentCount == 0) {
          ++errors ;
          return false ;
        }
        std::vector<NodeId> parents ;
        std::set<NodeId> uniqueParents ;
        double weightSum = 0.0 ;
        for (int parent = 0; parent < construction.parentCount; ++parent) {
          const NodeId parentId = construction.parentIds[parent] ;
          if (parentId == 0 || parentId == construction.node ||
                construction.parentNodeNumbers[parent] < 0 ||
                !std::isfinite(construction.parentWeights[parent]) ||
                construction.parentWeights[parent] <= 0.0 ||
                !uniqueParents.insert(parentId).second) {
            ++errors ;
            return false ;
          }
          parents.push_back(parentId) ;
          weightSum += construction.parentWeights[parent] ;
        }
        if (std::abs(weightSum - 1.0) > relativeTolerance ||
              persistentConstructedNodeId(
                    node_construction::value(construction.kind), parents) !=
                    construction.node) {
          ++errors ;
          return false ;
        }
        return true ;
      }

      node_origin::value nodeOrigin(int construction) {
        switch (construction) {
        case node_construction::edge:
          return node_origin::edge ;
        case node_construction::face:
          return node_origin::face ;
        case node_construction::cell:
          return node_origin::cell ;
        case node_construction::base_node:
        case node_construction::invalid:
          return node_origin::base_node ;
        }
        return node_origin::base_node ;
      }

      typedef std::map<NodeId, double> NodeExpansion ;

      bool expandNode(NodeId node,
            const std::map<NodeId, NodeGeometry>& previous,
            const std::map<NodeId, FineNodeConstruction>& current,
            std::map<NodeId, NodeExpansion>& memo, std::set<NodeId>& active,
            std::set<NodeId>& missing, size_t& cycles,
            NodeExpansion& expansion) {
        if (previous.find(node) != previous.end()) {
          expansion.clear() ;
          expansion[node] = 1.0 ;
          return true ;
        }
        const std::map<NodeId, NodeExpansion>::const_iterator cached =
              memo.find(node) ;
        if (cached != memo.end()) {
          expansion = cached->second ;
          return true ;
        }
        if (!active.insert(node).second) {
          ++cycles ;
          return false ;
        }
        const std::map<NodeId, FineNodeConstruction>::const_iterator found =
              current.find(node) ;
        if (found == current.end() ||
              found->second.kind == node_construction::base_node) {
          missing.insert(node) ;
          active.erase(node) ;
          return false ;
        }
        NodeExpansion result ;
        bool complete = true ;
        for (int parent = 0; parent < found->second.parentCount; ++parent) {
          NodeExpansion parentExpansion ;
          if (!expandNode(found->second.parentIds[parent], previous, current,
                    memo, active, missing, cycles, parentExpansion)) {
            complete = false ;
            continue ;
          }
          for (NodeExpansion::const_iterator source = parentExpansion.begin();
                source != parentExpansion.end(); ++source)
            result[source->first] +=
                  found->second.parentWeights[parent] * source->second ;
        }
        active.erase(node) ;
        if (!complete)
          return false ;
        memo[node] = result ;
        expansion.swap(result) ;
        return true ;
      }

      bool samePosition(const vector3d<double>& left,
            const vector3d<double>& right, double relativeTolerance) {
        const double scale = std::max(1.0, std::max(norm(left), norm(right))) ;
        return norm(left - right) <= relativeTolerance * scale ;
      }

      bool collectiveInputsValid(bool localToleranceValid, bool localStateValid,
            NodeTransitionReport& report) {
        const int localValidity[2] = {
              localToleranceValid ? 1 : 0, localStateValid ? 1 : 0} ;
        int globalValidity[2] = {0, 0} ;
        MPI_Allreduce(localValidity, globalValidity, 2, MPI_INT, MPI_MIN,
              MPI_COMM_WORLD) ;
        if (globalValidity[0] == 0) {
          report.status = node_transition_status::invalid_tolerance ;
          report.valid = false ;
          return false ;
        }
        if (globalValidity[1] == 0) {
          report.status = node_transition_status::missing_state ;
          report.valid = false ;
          return false ;
        }
        return true ;
      }
    }

    void reduceNodeTransitionReport(
          const NodeTransitionReport& local, NodeTransitionReport& global) {
      unsigned long long localCounts[9] = {
            static_cast<unsigned long long>(local.sourceNodes),
            static_cast<unsigned long long>(local.targetNodes),
            static_cast<unsigned long long>(local.retainedNodes),
            static_cast<unsigned long long>(local.createdNodes),
            static_cast<unsigned long long>(local.contributions),
            static_cast<unsigned long long>(local.missingSourceNodes),
            static_cast<unsigned long long>(local.cyclicConstructions),
            static_cast<unsigned long long>(local.inconsistentWeights),
            static_cast<unsigned long long>(local.inconsistentPositions)} ;
      unsigned long long globalCounts[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0} ;
      MPI_Allreduce(localCounts, globalCounts, 9, MPI_UNSIGNED_LONG_LONG,
            MPI_SUM, MPI_COMM_WORLD) ;
      const double localErrors[2] = {
            local.maximumWeightError, local.maximumPositionError} ;
      double globalErrors[2] = {0.0, 0.0} ;
      MPI_Allreduce(
            localErrors, globalErrors, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) ;
      const int localStatus = local.valid ? int(node_transition_status::valid)
                                          : int(local.status) ;
      int globalStatus = 0 ;
      MPI_Allreduce(
            &localStatus, &globalStatus, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

      global = NodeTransitionReport() ;
      global.sourceNodes = size_t(globalCounts[0]) ;
      global.targetNodes = size_t(globalCounts[1]) ;
      global.retainedNodes = size_t(globalCounts[2]) ;
      global.createdNodes = size_t(globalCounts[3]) ;
      global.contributions = size_t(globalCounts[4]) ;
      global.missingSourceNodes = size_t(globalCounts[5]) ;
      global.cyclicConstructions = size_t(globalCounts[6]) ;
      global.inconsistentWeights = size_t(globalCounts[7]) ;
      global.inconsistentPositions = size_t(globalCounts[8]) ;
      global.maximumWeightError = globalErrors[0] ;
      global.maximumPositionError = globalErrors[1] ;
      global.status = node_transition_status::value(globalStatus) ;
      global.valid = globalStatus == int(node_transition_status::valid) ;
    }

    CPTR<NodeRemap> buildNodeRemap(const store<NodeId>& previousIds,
          const store<vector3d<double>>& previousPositions,
          const std::vector<entitySet>& previousPartition,
          const store<NodeId>& currentIds,
          const store<FineNodeConstruction>& currentConstructions,
          const store<vector3d<double>>& currentPositions,
          const std::vector<entitySet>& currentPartition,
          NodeTransitionReport& report, double relativeTolerance) {
      report = NodeTransitionReport() ;
      const bool localToleranceValid =
            relativeTolerance >= 0.0 && std::isfinite(relativeTolerance) ;
      bool localStateValid =
            previousPartition.size() == size_t(MPI_processes) &&
            currentPartition.size() == size_t(MPI_processes) ;
      if (localStateValid)
        localStateValid =
              previousIds.domain() == previousPositions.domain() &&
              previousIds.domain() == previousPartition[MPI_rank] &&
              currentIds.domain() == currentConstructions.domain() &&
              currentIds.domain() == currentPositions.domain() &&
              currentIds.domain() == currentPartition[MPI_rank] ;
      if (!collectiveInputsValid(localToleranceValid, localStateValid, report))
        return CPTR<NodeRemap>() ;

      std::vector<std::vector<NodeGeometryRecord>> outgoingPrevious(
            MPI_processes) ;
      size_t localInvalid = 0 ;
      FORALL(previousIds.domain(), node) {
        const NodeId id = previousIds[node] ;
        if (id == 0 || !finitePosition(previousPositions[node])) {
          ++localInvalid ;
          continue ;
        }
        NodeGeometryRecord record = {
              id, {previousPositions[node].x, previousPositions[node].y,
                        previousPositions[node].z}} ;
        outgoingPrevious[directoryOwner(id)].push_back(record) ;
      }
      ENDFORALL ;
      const std::vector<NodeGeometryRecord> directoryRecords =
            exchangeByDestination(outgoingPrevious) ;
      std::map<NodeId, NodeGeometry> previousDirectory ;
      for (size_t record = 0; record < directoryRecords.size(); ++record) {
        const NodeGeometryRecord& source = directoryRecords[record] ;
        const vector3d<double> position(
              source.position[0], source.position[1], source.position[2]) ;
        if (!previousDirectory
                    .insert(std::make_pair(
                          source.node, NodeGeometry(source.node, position)))
                    .second)
          ++localInvalid ;
      }

      std::vector<std::vector<CurrentNodeRecord>> outgoingCurrentIds(
            MPI_processes) ;
      FORALL(currentIds.domain(), node) {
        if (currentIds[node] == 0 ||
              currentIds[node] != currentConstructions[node].node ||
              !finitePosition(currentPositions[node]))
          ++localInvalid ;
        CurrentNodeRecord record = {currentIds[node], node} ;
        outgoingCurrentIds[directoryOwner(currentIds[node])].push_back(record) ;
      }
      ENDFORALL ;
      const std::vector<CurrentNodeRecord> currentIdRecords =
            exchangeByDestination(outgoingCurrentIds) ;
      std::set<NodeId> globalCurrentIds ;
      for (size_t record = 0; record < currentIdRecords.size(); ++record)
        if (!globalCurrentIds.insert(currentIdRecords[record].node).second)
          ++localInvalid ;

      dstore<FineNodeConstruction> constructionByNumber ;
      FORALL(currentConstructions.domain(), node) {
        constructionByNumber[node] = currentConstructions[node] ;
      }
      ENDFORALL ;
      entitySet inspected ;
      bool closureComplete = false ;
      while (true) {
        entitySet required ;
        const entitySet available = constructionByNumber.domain() ;
        const entitySet uninspected = available - inspected ;
        FORALL(uninspected, node) {
          const FineNodeConstruction& construction = constructionByNumber[node] ;
          if (!validConstruction(construction, localInvalid, relativeTolerance))
            continue ;
          for (int parent = 0; parent < construction.parentCount; ++parent)
            required += construction.parentNodeNumbers[parent] ;
        }
        ENDFORALL ;
        inspected += uninspected ;
        entitySet missing = required - constructionByNumber.domain() ;
        // Do not clone or follow parents until every new record is valid.
        const int localCounts[2] = {int(missing.size()), int(localInvalid)} ;
        int globalCounts[2] = {0, 0} ;
        MPI_Allreduce(
              localCounts, globalCounts, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
        if (globalCounts[1] != 0) {
          report.status = node_transition_status::missing_state ;
          report.missingSourceNodes = size_t(globalCounts[1]) ;
          return CPTR<NodeRemap>() ;
        }
        if (globalCounts[0] == 0) {
          closureComplete = true ;
          break ;
        }
        storeRepP expanded = constructionByNumber.Rep() ;
        std::vector<entitySet> ownerPartition = currentPartition ;
        fill_clone(expanded, missing, ownerPartition) ;
        constructionByNumber.setRep(expanded) ;
        const int localProgress =
              (constructionByNumber.domain() - available).size() ;
        int globalProgress = 0 ;
        MPI_Allreduce(&localProgress, &globalProgress, 1, MPI_INT, MPI_SUM,
              MPI_COMM_WORLD) ;
        if (globalProgress == 0)
          break ;
      }
      if (!closureComplete)
        ++localInvalid ;

      std::map<NodeId, FineNodeConstruction> currentById ;
      FORALL(constructionByNumber.domain(), node) {
        const FineNodeConstruction& construction = constructionByNumber[node] ;
        if (!currentById.insert(std::make_pair(construction.node, construction))
                    .second)
          ++localInvalid ;
        for (int parent = 0; parent < construction.parentCount; ++parent) {
          const int parentNumber = construction.parentNodeNumbers[parent] ;
          if (!constructionByNumber.domain().inSet(parentNumber) ||
                constructionByNumber[parentNumber].node !=
                      construction.parentIds[parent])
            ++localInvalid ;
        }
      }
      ENDFORALL ;

      int invalid = int(localInvalid) ;
      int globalInvalid = 0 ;
      MPI_Allreduce(
            &invalid, &globalInvalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if (globalInvalid != 0) {
        report.status = node_transition_status::missing_state ;
        report.missingSourceNodes = size_t(globalInvalid) ;
        return CPTR<NodeRemap>() ;
      }

      std::vector<std::vector<NodeLookupRequest>> outgoingRequests(
            MPI_processes) ;
      for (std::map<NodeId, FineNodeConstruction>::const_iterator node =
                  currentById.begin();
            node != currentById.end(); ++node) {
        NodeLookupRequest request = {node->first, MPI_rank} ;
        outgoingRequests[directoryOwner(node->first)].push_back(request) ;
      }
      const std::vector<NodeLookupRequest> requests =
            exchangeByDestination(outgoingRequests) ;
      std::vector<std::vector<NodeLookupResponse>> outgoingResponses(
            MPI_processes) ;
      for (size_t request = 0; request < requests.size(); ++request) {
        const std::map<NodeId, NodeGeometry>::const_iterator found =
              previousDirectory.find(requests[request].node) ;
        NodeLookupResponse response = {
              requests[request].node, {0.0, 0.0, 0.0}, 0} ;
        if (found != previousDirectory.end()) {
          response.position[0] = found->second.position.x ;
          response.position[1] = found->second.position.y ;
          response.position[2] = found->second.position.z ;
          response.found = 1 ;
        }
        outgoingResponses[requests[request].requester].push_back(response) ;
      }
      const std::vector<NodeLookupResponse> responses =
            exchangeByDestination(outgoingResponses) ;
      std::map<NodeId, NodeGeometry> previousForTargets ;
      for (size_t response = 0; response < responses.size(); ++response)
        if (responses[response].found != 0) {
          const vector3d<double> position(responses[response].position[0],
                responses[response].position[1],
                responses[response].position[2]) ;
          previousForTargets[responses[response].node] =
                NodeGeometry(responses[response].node, position) ;
        }

      std::vector<NodeGeometry> sourceGeometry ;
      std::vector<NodeGeometry> targetGeometry ;
      std::vector<NodeContribution> contributions ;
      std::vector<NodeOrigin> origins ;
      std::set<NodeId> usedSources ;
      std::set<NodeId> missingSources ;
      std::map<NodeId, NodeExpansion> memo ;
      std::set<NodeId> active ;
      size_t cycles = 0 ;
      FORALL(currentIds.domain(), node) {
        const NodeId target = currentIds[node] ;
        targetGeometry.push_back(NodeGeometry(target, currentPositions[node])) ;
        origins.push_back(
              NodeOrigin(target, nodeOrigin(currentConstructions[node].kind))) ;
        NodeExpansion expansion ;
        if (previousForTargets.find(target) != previousForTargets.end()) {
          ++report.retainedNodes ;
          expansion[target] = 1.0 ;
        } else {
          ++report.createdNodes ;
          if (!expandNode(target, previousForTargets, currentById, memo, active,
                    missingSources, cycles, expansion))
            continue ;
        }
        for (NodeExpansion::const_iterator source = expansion.begin();
              source != expansion.end(); ++source) {
          contributions.push_back(
                NodeContribution(source->first, target, source->second)) ;
          usedSources.insert(source->first) ;
        }
      }
      ENDFORALL ;
      for (std::set<NodeId>::const_iterator source = usedSources.begin();
            source != usedSources.end(); ++source) {
        const std::map<NodeId, NodeGeometry>::const_iterator geometry =
              previousForTargets.find(*source) ;
        if (geometry == previousForTargets.end())
          missingSources.insert(*source) ;
        else
          sourceGeometry.push_back(geometry->second) ;
      }

      // Each previous node is held by exactly one persistent-ID directory
      // owner, so this count reduces to the global number of source nodes
      // without double-counting sources referenced by several target ranks.
      report.sourceNodes = previousDirectory.size() ;
      report.targetNodes = targetGeometry.size() ;
      report.contributions = contributions.size() ;
      report.missingSourceNodes = missingSources.size() ;
      report.cyclicConstructions = cycles ;
      int localFailure = int(missingSources.size() + cycles) ;
      int globalFailure = 0 ;
      MPI_Allreduce(
            &localFailure, &globalFailure, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if (globalFailure != 0) {
        report.status = cycles != 0
                              ? node_transition_status::cyclic_construction
                              : node_transition_status::missing_source_node ;
        return CPTR<NodeRemap>() ;
      }

      NodeRemapReport remapReport ;
      CPTR<NodeRemap> remap = NodeRemap::create(sourceGeometry, targetGeometry,
            contributions, origins, remapReport, relativeTolerance) ;
      report.inconsistentWeights = remapReport.inconsistentWeights ;
      report.inconsistentPositions = remapReport.inconsistentPositions ;
      report.maximumWeightError = remapReport.maximumWeightError ;
      report.maximumPositionError = remapReport.maximumPositionError ;
      int localRemapFailure = remap == static_cast<NodeRemap*>(0) ? 1 : 0 ;
      int globalRemapFailure = 0 ;
      MPI_Allreduce(&localRemapFailure, &globalRemapFailure, 1, MPI_INT,
            MPI_SUM, MPI_COMM_WORLD) ;
      if (globalRemapFailure != 0) {
        report.status = remapReport.inconsistentWeights != 0
                              ? node_transition_status::inconsistent_weights
                              : (remapReport.inconsistentPositions != 0
                                            ? node_transition_status::
                                                    inconsistent_positions
                                            : node_transition_status::
                                                    missing_source_node) ;
        return CPTR<NodeRemap>() ;
      }
      report.status = node_transition_status::valid ;
      report.valid = true ;
      return remap ;
    }

    CPTR<NodeRemap> redistributeNodeRemap(const CPTR<NodeRemap>& remap,
          const std::vector<entitySet>& nodePartition,
          const store<NodeId>& generatedNodeIds, NodeTransitionReport& report,
          double relativeTolerance) {
      const bool localToleranceValid =
            relativeTolerance >= 0.0 && std::isfinite(relativeTolerance) ;
      const bool localStateValid =
            remap != static_cast<NodeRemap*>(0) &&
            nodePartition.size() == size_t(MPI_processes) ;
      if (!collectiveInputsValid(localToleranceValid, localStateValid, report))
        return CPTR<NodeRemap>() ;

      std::map<NodeId, int> targetOwners ;
      size_t localInvalid = 0 ;
      for (int process = 0; process < MPI_processes; ++process) {
        FORALL(nodePartition[process], node) {
          if (!generatedNodeIds.domain().inSet(node) ||
                generatedNodeIds[node] == 0 ||
                !targetOwners
                       .insert(std::make_pair(generatedNodeIds[node], process))
                       .second)
            ++localInvalid ;
        }
        ENDFORALL ;
      }

      const std::vector<NodeGeometry>& sourceGeometry =
            remap->sourceNodeGeometry() ;
      const std::vector<NodeGeometry>& targetGeometry =
            remap->targetNodeGeometry() ;
      const std::vector<NodeContribution>& contributions =
            remap->nodeContributions() ;
      std::map<NodeId, NodeGeometry> sourceById ;
      for (size_t source = 0; source < sourceGeometry.size(); ++source)
        if (!sourceById
                    .insert(std::make_pair(
                          sourceGeometry[source].node, sourceGeometry[source]))
                    .second)
          ++localInvalid ;

      std::vector<std::vector<NodeGeometryRecord>> outgoingSources(
            MPI_processes) ;
      std::vector<std::vector<NodeGeometryRecord>> outgoingTargets(
            MPI_processes) ;
      std::vector<std::vector<NodeContributionRecord>> outgoingContributions(
            MPI_processes) ;
      std::vector<std::vector<NodeOriginRecord>> outgoingOrigins(MPI_processes) ;
      for (size_t target = 0; target < targetGeometry.size(); ++target) {
        const NodeId targetId = targetGeometry[target].node ;
        const std::map<NodeId, int>::const_iterator owner =
              targetOwners.find(targetId) ;
        size_t begin = 0, end = 0 ;
        NodeOrigin origin ;
        if (owner == targetOwners.end() ||
              !remap->nodeContributions(targetId, begin, end) ||
              !remap->nodeOrigin(targetId, origin)) {
          ++localInvalid ;
          continue ;
        }
        const int destination = owner->second ;
        NodeGeometryRecord targetRecord = {
              targetId, {targetGeometry[target].position.x,
                              targetGeometry[target].position.y,
                              targetGeometry[target].position.z}} ;
        outgoingTargets[destination].push_back(targetRecord) ;
        NodeOriginRecord originRecord = {targetId, int(origin.kind)} ;
        outgoingOrigins[destination].push_back(originRecord) ;
        for (size_t entry = begin; entry < end; ++entry) {
          const NodeContribution& contribution = contributions[entry] ;
          const std::map<NodeId, NodeGeometry>::const_iterator source =
                sourceById.find(contribution.sourceNode) ;
          if (source == sourceById.end()) {
            ++localInvalid ;
            continue ;
          }
          NodeContributionRecord contributionRecord = {contribution.sourceNode,
                contribution.targetNode, contribution.weight} ;
          outgoingContributions[destination].push_back(contributionRecord) ;
          NodeGeometryRecord sourceRecord = {source->second.node,
                {source->second.position.x, source->second.position.y,
                      source->second.position.z}} ;
          outgoingSources[destination].push_back(sourceRecord) ;
        }
      }

      int invalid = int(localInvalid) ;
      int globalInvalid = 0 ;
      MPI_Allreduce(
            &invalid, &globalInvalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if (globalInvalid != 0) {
        report.status = node_transition_status::missing_state ;
        report.valid = false ;
        return CPTR<NodeRemap>() ;
      }

      const std::vector<NodeGeometryRecord> receivedSources =
            exchangeByDestination(outgoingSources) ;
      const std::vector<NodeGeometryRecord> receivedTargets =
            exchangeByDestination(outgoingTargets) ;
      const std::vector<NodeContributionRecord> receivedContributions =
            exchangeByDestination(outgoingContributions) ;
      const std::vector<NodeOriginRecord> receivedOrigins =
            exchangeByDestination(outgoingOrigins) ;

      std::map<NodeId, NodeGeometry> uniqueSources ;
      std::vector<NodeGeometry> installedTargets ;
      std::vector<NodeContribution> installedContributions ;
      std::vector<NodeOrigin> installedOrigins ;
      for (size_t source = 0; source < receivedSources.size(); ++source) {
        const vector3d<double> position(receivedSources[source].position[0],
              receivedSources[source].position[1],
              receivedSources[source].position[2]) ;
        const std::map<NodeId, NodeGeometry>::const_iterator prior =
              uniqueSources.find(receivedSources[source].node) ;
        if (prior == uniqueSources.end())
          uniqueSources[receivedSources[source].node] =
                NodeGeometry(receivedSources[source].node, position) ;
        else if (!samePosition(
                       prior->second.position, position, relativeTolerance))
          ++localInvalid ;
      }
      for (size_t target = 0; target < receivedTargets.size(); ++target) {
        const vector3d<double> position(receivedTargets[target].position[0],
              receivedTargets[target].position[1],
              receivedTargets[target].position[2]) ;
        installedTargets.push_back(
              NodeGeometry(receivedTargets[target].node, position)) ;
      }
      for (size_t entry = 0; entry < receivedContributions.size(); ++entry)
        installedContributions.push_back(
              NodeContribution(receivedContributions[entry].sourceNode,
                    receivedContributions[entry].targetNode,
                    receivedContributions[entry].weight)) ;
      for (size_t origin = 0; origin < receivedOrigins.size(); ++origin)
        installedOrigins.push_back(
              NodeOrigin(receivedOrigins[origin].targetNode,
                    node_origin::value(receivedOrigins[origin].kind))) ;
      std::vector<NodeGeometry> installedSources ;
      for (std::map<NodeId, NodeGeometry>::const_iterator source =
                  uniqueSources.begin();
            source != uniqueSources.end(); ++source)
        installedSources.push_back(source->second) ;

      NodeRemapReport remapReport ;
      CPTR<NodeRemap> installed = NodeRemap::create(installedSources,
            installedTargets, installedContributions, installedOrigins,
            remapReport, relativeTolerance) ;
      int localValid = localInvalid == 0 &&
                                   installed != static_cast<NodeRemap*>(0) &&
                                   remapReport.valid
                             ? 1
                             : 0 ;
      int globalValid = 0 ;
      MPI_Allreduce(
            &localValid, &globalValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ;
      if (globalValid == 0) {
        report.valid = false ;
        report.inconsistentWeights += remapReport.inconsistentWeights ;
        report.inconsistentPositions += remapReport.inconsistentPositions ;
        report.maximumWeightError = std::max(
              report.maximumWeightError, remapReport.maximumWeightError) ;
        report.maximumPositionError = std::max(
              report.maximumPositionError, remapReport.maximumPositionError) ;
        report.status =
              remapReport.inconsistentWeights != 0
                    ? node_transition_status::inconsistent_weights
                    : (remapReport.inconsistentPositions != 0
                                  ? node_transition_status::
                                          inconsistent_positions
                                  : node_transition_status::missing_state) ;
        return CPTR<NodeRemap>() ;
      }
      return installed ;
    }
  }
}

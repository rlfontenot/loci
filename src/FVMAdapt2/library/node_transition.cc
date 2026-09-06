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

#include <FVMAdapt2/node_transition.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>

namespace Loci {

  NodeTransitionReport::NodeTransitionReport()
      : status(node_transition_status::missing_state), valid(false),
        sourceNodes(0), targetNodes(0), retainedNodes(0), createdNodes(0),
        contributions(0), missingSourceNodes(0), cyclicConstructions(0),
        inconsistentWeights(0), inconsistentPositions(0),
        maximumWeightError(0.0), maximumPositionError(0.0) {}

  namespace {
    const std::uint64_t nodeIdOffset = UINT64_C(14695981039346656037) ;
    const std::uint64_t nodeIdPrime = UINT64_C(1099511628211) ;

    void appendHashByte(std::uint64_t& hash, unsigned char value) {
      hash ^= std::uint64_t(value) ;
      hash *= nodeIdPrime ;
    }

    void appendHashInt(std::uint64_t& hash, std::int64_t value) {
      const std::uint64_t bits = static_cast<std::uint64_t>(value) ;
      for (int byte = 0; byte < 8; ++byte)
        appendHashByte(hash, static_cast<unsigned char>(bits >> (8 * byte))) ;
    }

    NodeId signedHash(std::uint64_t hash) {
      NodeId result = 0 ;
      std::memcpy(&result, &hash, sizeof(result)) ;
      // Zero is reserved for malformed/default construction records.
      if (result == 0)
        result = std::numeric_limits<NodeId>::min() ;
      return result ;
    }

    std::vector<NodeId> parentIds(const std::vector<NodeParent>& parents) {
      std::vector<NodeId> result(parents.size()) ;
      for (size_t parent = 0; parent < parents.size(); ++parent)
        result[parent] = parents[parent].node ;
      return result ;
    }
  }

  NodeId persistentBaseNodeId(long long baseFileNumber) {
    if (baseFileNumber < 0)
      return 0 ;
    std::uint64_t hash = nodeIdOffset ;
    appendHashByte(hash, 0x4e) ; // 'N': node-key encoding version one
    appendHashByte(hash, 0x42) ; // 'B': original/base node
    appendHashInt(hash, static_cast<std::int64_t>(baseFileNumber)) ;
    return signedHash(hash) ;
  }

  NodeId persistentConstructedNodeId(node_construction::value kind,
        const std::vector<NodeId>& directParents) {
    if (kind < node_construction::edge || kind > node_construction::cell ||
          directParents.empty())
      return 0 ;
    std::vector<NodeId> parents = directParents ;
    std::sort(parents.begin(), parents.end()) ;
    std::uint64_t hash = nodeIdOffset ;
    appendHashByte(hash, 0x4e) ; // 'N': node-key encoding version one
    appendHashByte(hash, 0x47) ; // 'G': generated node
    appendHashInt(hash, static_cast<std::int64_t>(kind)) ;
    appendHashInt(hash, static_cast<std::int64_t>(parents.size())) ;
    for (size_t parent = 0; parent < parents.size(); ++parent)
      appendHashInt(hash, static_cast<std::int64_t>(parents[parent])) ;
    return signedHash(hash) ;
  }

  NodeConstruction NodeConstruction::baseNode(
        long long fileNumber, const vector3d<double>& nodePosition) {
    NodeConstruction result ;
    result.node = persistentBaseNodeId(fileNumber) ;
    result.position = nodePosition ;
    result.kind = node_construction::base_node ;
    result.baseFileNumber = fileNumber ;
    return result ;
  }

  NodeConstruction NodeConstruction::constructed(
        node_construction::value constructionKind,
        const vector3d<double>& nodePosition,
        const std::vector<NodeParent>& directParents) {
    NodeConstruction result ;
    result.position = nodePosition ;
    result.kind = constructionKind ;
    result.parents = directParents ;
    result.node = persistentConstructedNodeId(
          constructionKind, parentIds(directParents)) ;
    return result ;
  }
}

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
#ifndef FVMADAPT2_NODE_TRANSITION_H
#define FVMADAPT2_NODE_TRANSITION_H

#include <FVMAdapt2/mesh_transfer.h>

#include <cstddef>
#include <vector>

namespace Loci {

  namespace node_construction {
    /// These values are part of the persistent node-identity encoding.
    enum value { invalid = 0, base_node = 1, edge = 2, face = 3, cell = 4 } ;
  }

  /// Return the persistent identity of an original VOG node.
  NodeId persistentBaseNodeId(long long baseFileNumber) ;

  /// Return the persistent identity of a generated node.
  ///
  /// The identity depends on the construction kind and canonical parent-node
  /// identities, not entity numbering, coordinates, weights, or MPI layout.
  NodeId persistentConstructedNodeId(node_construction::value kind,
        const std::vector<NodeId>& directParents) ;

  /// One direct geometric parent of a generated node.
  struct NodeParent {
    NodeId node ;
    double weight ;

    NodeParent() : node(0), weight(0.0) {}
    NodeParent(NodeId parentNode, double parentWeight)
        : node(parentNode), weight(parentWeight) {}
  } ;

  /// A mesh node and the geometric parents used to construct it.
  struct NodeConstruction {
    NodeId node ;
    vector3d<double> position ;
    node_construction::value kind ;
    long long baseFileNumber ;
    std::vector<NodeParent> parents ;

    NodeConstruction()
        : node(0), position(0.0, 0.0, 0.0), kind(node_construction::invalid),
          baseFileNumber(-1) {}

    static NodeConstruction baseNode(
          long long fileNumber, const vector3d<double>& nodePosition) ;

    static NodeConstruction constructed(
          node_construction::value constructionKind,
          const vector3d<double>& nodePosition,
          const std::vector<NodeParent>& directParents) ;
  } ;


}

#endif

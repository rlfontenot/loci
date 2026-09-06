//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef FVMADAPT2_NODE_EDGE_H
#define FVMADAPT2_NODE_EDGE_H

#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <queue>
#include <cmath>
#include "defines.h"
#include "node_transition.h"
using std::queue;
//head--------> tail, edge point from head to tail
//edge2node[0]-------->edge2node[1]

//declaration of class Node, Edge and Face
class Node{
  struct ConstructionParent {
    Node* node ;
    double weight ;

    ConstructionParent(Node* parent, double parentWeight)
        : node(parent), weight(parentWeight) {}
  } ;

public:
  //constructors
  Node()
      : index(0), tag(0), constructionKind(Loci::node_construction::invalid),
        baseFileNumber(-1), persistentId(0) {}
  Node(vect3d& p0)
      : p(p0), index(0), tag(0),
        constructionKind(Loci::node_construction::invalid), baseFileNumber(-1),
        persistentId(0) {}
  Node(const Loci::vector3d<double>& p0)
      : p(p0), index(0), tag(0),
        constructionKind(Loci::node_construction::invalid), baseFileNumber(-1),
        persistentId(0) {}
  Node(vect3d& p0, int32 n)
      : p(p0), index(n), tag(0),
        constructionKind(Loci::node_construction::base_node), baseFileNumber(n),
        persistentId(Loci::persistentBaseNodeId(n)) {}
  Node(const Loci::vector3d<double>& p0, int32 n)
      : p(p0), index(n), tag(0),
        constructionKind(Loci::node_construction::base_node), baseFileNumber(n),
        persistentId(Loci::persistentBaseNodeId(n)) {}

  /// Construct a node from the exact geometric parents used by refinement.
  static Node* constructed(Loci::node_construction::value kind,
        const std::vector<Node*>& parents,
        const std::vector<double>& unnormalizedWeights) {
    if (parents.empty() || parents.size() != unnormalizedWeights.size())
      return new Node() ;
    double weightSum = 0.0 ;
    vect3d position(0.0, 0.0, 0.0) ;
    bool hasProvenance = kind >= Loci::node_construction::edge &&
                         kind <= Loci::node_construction::cell ;
    for (size_t parent = 0; parent < parents.size(); ++parent) {
      if (parents[parent] == 0 || !std::isfinite(unnormalizedWeights[parent]))
        return new Node() ;
      if (parents[parent]->persistentId == 0)
        hasProvenance = false ;
      weightSum += unnormalizedWeights[parent] ;
      position += unnormalizedWeights[parent] * parents[parent]->p ;
    }
    if (!std::isfinite(weightSum) || weightSum <= 0.0)
      return new Node() ;
    Node* node = new Node(position / weightSum) ;
    // Many legacy tree builders carry geometry only. Preserve their centroid
    // behavior while leaving ancestry unavailable for that construction path.
    if (!hasProvenance)
      return node ;
    node->constructionKind = kind ;
    std::vector<Loci::NodeId> parentIds ;
    parentIds.reserve(parents.size()) ;
    node->constructionParents.reserve(parents.size()) ;
    for (size_t parent = 0; parent < unnormalizedWeights.size(); ++parent) {
      const double weight = unnormalizedWeights[parent] / weightSum ;
      parentIds.push_back(parents[parent]->persistentId) ;
      node->constructionParents.push_back(
            ConstructionParent(parents[parent], weight)) ;
    }
    node->persistentId = Loci::persistentConstructedNodeId(kind, parentIds) ;
    return node ;
  }

  /// Return this node's persistent construction record.
  bool construction(Loci::NodeConstruction& result) const {
    if (constructionKind == Loci::node_construction::base_node) {
      result = Loci::NodeConstruction::baseNode(baseFileNumber, p) ;
      return result.node == persistentId && persistentId != 0 ;
    }
    if (constructionKind < Loci::node_construction::edge ||
          constructionKind > Loci::node_construction::cell ||
          persistentId == 0 || constructionParents.empty())
      return false ;
    std::vector<Loci::NodeParent> parents ;
    parents.reserve(constructionParents.size()) ;
    for (size_t parent = 0; parent < constructionParents.size(); ++parent) {
      if (constructionParents[parent].node == 0 ||
            constructionParents[parent].node->persistentId == 0)
        return false ;
      parents.push_back(
            Loci::NodeParent(constructionParents[parent].node->persistentId,
                  constructionParents[parent].weight)) ;
    }
    result = Loci::NodeConstruction::constructed(constructionKind, p, parents) ;
    return result.node == persistentId ;
  }

  /// Serialize construction provenance beside the node coordinate.
  bool fineConstruction(Loci::FineNodeConstruction& result) const {
    Loci::NodeConstruction constructionRecord ;
    if (!construction(constructionRecord) ||
          constructionRecord.parents.size() >
                size_t(Loci::FineNodeConstruction::maximumParents))
      return false ;
    result = Loci::FineNodeConstruction() ;
    result.node = constructionRecord.node ;
    result.kind = int(constructionRecord.kind) ;
    result.baseFileNumber = constructionRecord.baseFileNumber ;
    result.parentCount = int(constructionRecord.parents.size()) ;
    for (size_t parent = 0; parent < constructionRecord.parents.size();
          ++parent) {
      result.parentIds[parent] = constructionRecord.parents[parent].node ;
      result.parentNodeNumbers[parent] =
            int(constructionParents[parent].node->index) ;
      result.parentWeights[parent] = constructionRecord.parents[parent].weight ;
    }
    return true ;
  }

public:
  vect3d p; //coordiantes
  int32 index;//the index of node in input or output grid file, start with 1
  char tag; // 1 or 0, indicate the node need to be refined or not 

private:
  Loci::node_construction::value constructionKind ;
  long long baseFileNumber ;
  Loci::NodeId persistentId ;
  std::vector<ConstructionParent> constructionParents ;
};

class Edge{
public:
  //constructor
  Edge(Node* p0, Node* p1):head(p0), tail(p1), child(0), parent(0),level(0) {}
  Edge(Node* p0, Node* p1, int lev):head(p0), tail(p1),child(0),parent(0), level(lev){}
  Edge(Node* p0, Node* p1, int lev, Edge* p):head(p0), tail(p1),child(0),parent(p), level(lev){}
  //constructor
  Edge():child(0), parent(0),level(0){}
  //destructor
  ~Edge(){
    if(child!=0){
      if(child[1] != 0){
	delete child[1];
	child[1] = 0;
      }
      if(child[0] != 0){
	delete child[0];
	child[0] = 0;
      }
      delete[] child;
      child = 0;
    }
  }
  //calculate the middle point of the edge
  inline Node* centroid(){
    std::vector<Node*> parents(2) ;
    parents[0] = head ;
    parents[1] = tail ;
    return Node::constructed(
          Loci::node_construction::edge, parents, std::vector<double>(2, 0.5)) ;
  }

  inline double length(){
    return norm(head->p - tail->p);
  }
  
  //check if the edge contain the node,
  //if return 0, the node is the head
  //if return 1, the node is the tail
  //if return -1, the edge does't contain the node
  inline  int containNode(Node* theNode){
    if(head == theNode) return 0;
    if(tail == theNode) return 1;
    return -1;
  }
  
  //split function, split the edge in the middle
  inline void split(std::list<Node*>& node_list){
    if(child == 0){
      Node *center = centroid();
      node_list.push_back(center);
    
      child = new Edge*[2];
      child[0] = new Edge(head, center, level+1, this);
      child[1] = new Edge(center, tail, level+1, this);
    }
  }

  inline bool depth_greater_than_1(){
    if(child == 0) return false;
    if(child[0]->child != 0) return true;
    if(child[1]->child != 0) return true;
    return false;
  }


  inline double get_length(){
    return norm(tail->p - head->p);
  }



  
  void sort_leaves(std::list<Edge*>& leaves);
  int get_depth();
  //For QuadFace
  void resplit(const std::vector<char>& edgePlan, bool needReverse,
               std::list<Node*>& node_list);
  //for Face
  void resplit(const std::vector<char>& edgePlan,
               std::list<Node*>& node_list);
  
public:
  Node* head;  //edge2node[0]
  Node* tail;//edge2node[1]
  //head-------->tail the edge points from head to tail. the vector is tail->p - head->p
  Edge** child;
  Edge* parent;
  int level; //the level of tree structure
};


inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list){
   
  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  node_list.clear();
  
  for(std::list<Edge*>::iterator p = edge_list.begin(); p != edge_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  edge_list.clear();
  
}

inline void cleanup_list(std::list<Node*>& node_list){
                       
  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  node_list.clear();
}

#endif

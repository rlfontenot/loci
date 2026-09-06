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
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <utility>

namespace Loci {

  FaceTransitionReport::FaceTransitionReport()
      : status(face_transition_status::invalid_identity), valid(false),
        sourceFaces(0), targetFaces(0), unsupportedRootPlans(0),
        invalidIdentities(0), invalidPolygons(0) {}

  bool operator<(const FaceKey& left, const FaceKey& right) {
    if (left.origin != right.origin)
      return left.origin < right.origin ;
    if (left.root != right.root)
      return left.root < right.root ;
    if (left.firstPath != right.firstPath)
      return left.firstPath < right.firstPath ;
    return left.secondPath < right.secondPath ;
  }

  bool operator==(const FaceKey& left, const FaceKey& right) {
    return left.origin == right.origin && left.root == right.root &&
           left.firstPath == right.firstPath &&
           left.secondPath == right.secondPath ;
  }

  bool operator<(const CellKey& left, const CellKey& right) {
    if (left.root != right.root)
      return left.root < right.root ;
    return left.path < right.path ;
  }

  bool operator==(const CellKey& left, const CellKey& right) {
    return left.root == right.root && left.path == right.path ;
  }

  namespace {
    // FNV-1a over an explicitly byte-ordered key encoding.  This is a stable
    // FVMAdapt2 encoding, not an implementation-defined std::hash value.
    const std::uint64_t faceIdOffset = UINT64_C(14695981039346656037) ;
    const std::uint64_t faceIdPrime = UINT64_C(1099511628211) ;

    void appendHashByte(std::uint64_t& hash, unsigned char value) {
      hash ^= std::uint64_t(value) ;
      hash *= faceIdPrime ;
    }

    void appendHashInt(std::uint64_t& hash, std::int64_t value) {
      const std::uint64_t bits = static_cast<std::uint64_t>(value) ;
      for (int byte = 0; byte < 8; ++byte)
        appendHashByte(hash, static_cast<unsigned char>(bits >> (8 * byte))) ;
    }

    void appendHashPath(std::uint64_t& hash, const std::vector<int>& path) {
      appendHashInt(hash, static_cast<std::int64_t>(path.size() / 2)) ;
      for (size_t entry = 0; entry < path.size(); ++entry)
        appendHashInt(hash, static_cast<std::int64_t>(path[entry])) ;
    }

    bool validRefinementPath(const std::vector<int>& path) {
      if (path.size() % 2 != 0)
        return false ;
      for (size_t step = 0; step < path.size(); step += 2)
        if (path[step] <= 0 || path[step + 1] < 0)
          return false ;
      return true ;
    }

    bool validFaceKeyStructure(const FaceKey& key) {
      if (key.root < 0 || !validRefinementPath(key.firstPath) ||
            !validRefinementPath(key.secondPath))
        return false ;
      if (key.origin == face_origin::base_face)
        return key.secondPath.empty() ;
      return key.origin == face_origin::cell_interior &&
             key.firstPath < key.secondPath ;
    }

    FaceId signedHash(std::uint64_t hash) {
      FaceId result = 0 ;
      std::memcpy(&result, &hash, sizeof(result)) ;
      // Zero is reserved for an absent or malformed persistent identity.
      if (result == 0)
        result = std::numeric_limits<FaceId>::min() ;
      return result ;
    }
  }

  FaceId persistentFaceId(const FaceKey& key) {
    if (!validFaceKeyStructure(key))
      return 0 ;
    std::uint64_t hash = faceIdOffset ;
    appendHashByte(hash, 0x46) ; // 'F': face key
    appendHashByte(hash, 0x02) ; // encoding version two
    appendHashInt(hash, static_cast<std::int64_t>(key.origin)) ;
    appendHashInt(hash, static_cast<std::int64_t>(key.root)) ;
    appendHashPath(hash, key.firstPath) ;
    appendHashPath(hash, key.secondPath) ;
    return signedHash(hash) ;
  }

  CellId persistentCellId(const CellKey& key) {
    if (key.root < 0 || !validRefinementPath(key.path))
      return 0 ;
    std::uint64_t hash = faceIdOffset ;
    appendHashByte(hash, 0x43) ; // 'C': cell key
    appendHashByte(hash, 0x02) ; // encoding version two
    appendHashInt(hash, static_cast<std::int64_t>(key.root)) ;
    appendHashPath(hash, key.path) ;
    return signedHash(hash) ;
  }

  CellId persistentCellId(int rootCell, const std::vector<int>& leafPath) {
    return persistentCellId(CellKey(rootCell, leafPath)) ;
  }

  namespace {

    struct Point2 {
      double x ;
      double y ;

      Point2() : x(0.0), y(0.0) {}
      Point2(double px, double py) : x(px), y(py) {}
    } ;

    struct PolygonGeometry {
      bool valid ;
      bool planar ;
      int droppedAxis ;
      double scale ;
      double area ;
      vector3d<double> origin ;
      vector3d<double> centroid ;
      vector3d<double> normal ;
      std::vector<Point2> projected ;

      PolygonGeometry()
          : valid(false), planar(true), droppedAxis(0), scale(0.0), area(0.0),
            origin(0.0, 0.0, 0.0), centroid(0.0, 0.0, 0.0),
            normal(0.0, 0.0, 0.0) {}
    } ;

    bool finiteVector(const vector3d<double>& value) {
      return std::isfinite(value.x) && std::isfinite(value.y) &&
             std::isfinite(value.z) ;
    }

    Point2 projectPoint(const vector3d<double>& point, int droppedAxis) {
      if (droppedAxis == 0)
        return Point2(point.y, point.z) ;
      if (droppedAxis == 1)
        return Point2(point.z, point.x) ;
      return Point2(point.x, point.y) ;
    }

    vector3d<double> liftPoint(
          const Point2& point, const PolygonGeometry& plane) {
      vector3d<double> result(0.0, 0.0, 0.0) ;
      if (plane.droppedAxis == 0) {
        result.y = point.x ;
        result.z = point.y ;
        result.x = (-plane.normal.y * result.y - plane.normal.z * result.z) /
                   plane.normal.x ;
      } else if (plane.droppedAxis == 1) {
        result.z = point.x ;
        result.x = point.y ;
        result.y = (-plane.normal.z * result.z - plane.normal.x * result.x) /
                   plane.normal.y ;
      } else {
        result.x = point.x ;
        result.y = point.y ;
        result.z = (-plane.normal.x * result.x - plane.normal.y * result.y) /
                   plane.normal.z ;
      }
      return plane.origin + result ;
    }

    double cross2(
          const Point2& first, const Point2& second, const Point2& third) {
      return (second.x - first.x) * (third.y - first.y) -
             (second.y - first.y) * (third.x - first.x) ;
    }

    double signedArea2(const std::vector<Point2>& polygon) {
      double twiceArea = 0.0 ;
      if (polygon.size() < 3)
        return twiceArea ;
      for (size_t point = 1; point + 1 < polygon.size(); ++point)
        twiceArea += cross2(polygon[0], polygon[point], polygon[point + 1]) ;
      return twiceArea ;
    }

    bool convexPolygon(const std::vector<Point2>& polygon, double tolerance) {
      if (polygon.size() < 3)
        return false ;
      int sign = 0 ;
      for (size_t point = 0; point < polygon.size(); ++point) {
        const double turn =
              cross2(polygon[point], polygon[(point + 1) % polygon.size()],
                    polygon[(point + 2) % polygon.size()]) ;
        if (std::abs(turn) <= tolerance)
          continue ;
        const int currentSign = turn > 0.0 ? 1 : -1 ;
        if (sign != 0 && currentSign != sign)
          return false ;
        sign = currentSign ;
      }
      return sign != 0 ;
    }

    bool polygonGeometry(const std::vector<vector3d<double>>& vertices,
          double relativeTolerance, PolygonGeometry& geometry,
          bool allowWarped = false) {
      geometry = PolygonGeometry() ;
      if (vertices.size() < 3)
        return false ;
      for (size_t vertex = 0; vertex < vertices.size(); ++vertex)
        if (!finiteVector(vertices[vertex]))
          return false ;

      vector3d<double> areaVector(0.0, 0.0, 0.0) ;
      geometry.origin = vertices[0] ;
      const vector3d<double>& origin = geometry.origin ;
      for (size_t vertex = 1; vertex + 1 < vertices.size(); ++vertex)
        areaVector +=
              cross(vertices[vertex] - origin, vertices[vertex + 1] - origin) ;
      areaVector *= 0.5 ;
      geometry.area = norm(areaVector) ;
      if (!std::isfinite(geometry.area) || geometry.area <= 0.0)
        return false ;
      geometry.normal = areaVector / geometry.area ;

      double scale = 0.0 ;
      for (size_t vertex = 0; vertex < vertices.size(); ++vertex)
        scale = std::max(scale, norm(vertices[vertex] - origin)) ;
      if (!std::isfinite(scale) || scale <= 0.0)
        return false ;
      geometry.scale = scale ;
      const double tolerance = relativeTolerance * scale ;
      for (size_t vertex = 1; vertex < vertices.size(); ++vertex)
        if (std::abs(dot(geometry.normal, vertices[vertex] - origin)) >
              tolerance)
          geometry.planar = false ;
      if (!geometry.planar && !allowWarped)
        return false ;

      const double nx = std::abs(geometry.normal.x) ;
      const double ny = std::abs(geometry.normal.y) ;
      const double nz = std::abs(geometry.normal.z) ;
      geometry.droppedAxis = nx >= ny && nx >= nz ? 0 : (ny >= nz ? 1 : 2) ;
      geometry.projected.resize(vertices.size()) ;
      for (size_t vertex = 0; vertex < vertices.size(); ++vertex)
        geometry.projected[vertex] =
              projectPoint(vertices[vertex] - origin, geometry.droppedAxis) ;
      const double twiceProjectedArea = signedArea2(geometry.projected) ;
      const double areaTolerance = relativeTolerance * scale * scale ;
      if (std::abs(twiceProjectedArea) <= areaTolerance ||
            !convexPolygon(geometry.projected, areaTolerance))
        return false ;
      if (twiceProjectedArea < 0.0)
        std::reverse(geometry.projected.begin(), geometry.projected.end()) ;

      // A new or removed interior face needs an origin, not an overlap.
      // For a warped face, retain the vector-area magnitude and wireframe
      // center used by Loci's default face geometry. This is not a surface
      // intersection convention: changing warped faces cannot be clipped.
      if (!geometry.planar) {
        vector3d<double> centerOffset(0.0, 0.0, 0.0) ;
        double perimeter = 0.0 ;
        for (size_t vertex = 0; vertex < vertices.size(); ++vertex) {
          const vector3d<double> first = vertices[vertex] - origin ;
          const vector3d<double> second =
                vertices[(vertex + 1) % vertices.size()] - origin ;
          const double length = norm(second - first) ;
          centerOffset += 0.5 * length * (first + second) ;
          perimeter += length ;
        }
        if (!std::isfinite(perimeter) || perimeter <= 0.0)
          return false ;
        geometry.centroid = origin + centerOffset / perimeter ;
        geometry.valid = finiteVector(geometry.centroid) ;
        return geometry.valid ;
      }

      vector3d<double> centroidOffset(0.0, 0.0, 0.0) ;
      double triangleAreaSum = 0.0 ;
      for (size_t vertex = 1; vertex + 1 < vertices.size(); ++vertex) {
        const vector3d<double> first = vertices[vertex] - origin ;
        const vector3d<double> second = vertices[vertex + 1] - origin ;
        const vector3d<double> triangleVector = 0.5 * cross(first, second) ;
        const double triangleArea = dot(triangleVector, geometry.normal) ;
        centroidOffset += triangleArea * (first + second) / 3.0 ;
        triangleAreaSum += triangleArea ;
      }
      if (std::abs(triangleAreaSum) <= areaTolerance)
        return false ;
      geometry.centroid = origin + centroidOffset / triangleAreaSum ;
      geometry.valid = finiteVector(geometry.centroid) ;
      return geometry.valid ;
    }

    Point2 lineIntersection(const Point2& first, const Point2& second,
          const Point2& clipFirst, const Point2& clipSecond, double tolerance) {
      const double sx = second.x - first.x ;
      const double sy = second.y - first.y ;
      const double cx = clipSecond.x - clipFirst.x ;
      const double cy = clipSecond.y - clipFirst.y ;
      const double denominator = sx * cy - sy * cx ;
      if (std::abs(denominator) <= tolerance)
        return second ;
      const double t =
            ((clipFirst.x - first.x) * cy - (clipFirst.y - first.y) * cx) /
            denominator ;
      return Point2(first.x + t * sx, first.y + t * sy) ;
    }

    std::vector<Point2> clipPolygon(std::vector<Point2> subject,
          const std::vector<Point2>& clip, double tolerance) {
      for (size_t edge = 0; edge < clip.size() && !subject.empty(); ++edge) {
        const Point2& clipFirst = clip[edge] ;
        const Point2& clipSecond = clip[(edge + 1) % clip.size()] ;
        const std::vector<Point2> input = subject ;
        subject.clear() ;
        Point2 previous = input.back() ;
        bool previousInside =
              cross2(clipFirst, clipSecond, previous) >= -tolerance ;
        for (size_t point = 0; point < input.size(); ++point) {
          const Point2 current = input[point] ;
          const bool currentInside =
                cross2(clipFirst, clipSecond, current) >= -tolerance ;
          if (currentInside != previousInside)
            subject.push_back(lineIntersection(
                  previous, current, clipFirst, clipSecond, tolerance)) ;
          if (currentInside)
            subject.push_back(current) ;
          previous = current ;
          previousInside = currentInside ;
        }
      }
      return subject ;
    }

    bool overlapGeometry(const PolygonGeometry& source,
          const std::vector<vector3d<double>>& targetVertices,
          const PolygonGeometry& target, double relativeTolerance,
          double& overlapArea, vector3d<double>& overlapCentroid,
          int& orientation) {
      overlapArea = 0.0 ;
      overlapCentroid = vector3d<double>(0.0, 0.0, 0.0) ;
      if (!source.planar || !target.planar)
        return false ;
      orientation = dot(source.normal, target.normal) >= 0.0 ? 1 : -1 ;
      const double normalAgreement =
            std::abs(dot(source.normal, target.normal)) ;
      if (1.0 - normalAgreement > 10.0 * relativeTolerance)
        return true ;
      const double scale = std::max(source.scale, target.scale) ;
      for (size_t vertex = 0; vertex < targetVertices.size(); ++vertex)
        if (std::abs(
                  dot(source.normal, targetVertices[vertex] - source.origin)) >
              relativeTolerance * scale)
          return true ;

      std::vector<Point2> targetProjected(targetVertices.size()) ;
      for (size_t vertex = 0; vertex < targetVertices.size(); ++vertex)
        targetProjected[vertex] = projectPoint(
              targetVertices[vertex] - source.origin, source.droppedAxis) ;
      if (signedArea2(targetProjected) < 0.0)
        std::reverse(targetProjected.begin(), targetProjected.end()) ;
      std::vector<Point2> overlap = clipPolygon(source.projected,
            targetProjected, relativeTolerance * scale * scale) ;
      if (overlap.size() < 3)
        return true ;

      double twiceArea = signedArea2(overlap) ;
      if (twiceArea < 0.0) {
        std::reverse(overlap.begin(), overlap.end()) ;
        twiceArea = -twiceArea ;
      }
      const double projectedArea = 0.5 * twiceArea ;
      if (projectedArea <= relativeTolerance * scale * scale)
        return true ;

      double centroidX = 0.0 ;
      double centroidY = 0.0 ;
      double centroidDenominator = 0.0 ;
      for (size_t point = 0; point < overlap.size(); ++point) {
        const Point2& first = overlap[point] ;
        const Point2& second = overlap[(point + 1) % overlap.size()] ;
        const double edgeCross = first.x * second.y - second.x * first.y ;
        centroidX += (first.x + second.x) * edgeCross ;
        centroidY += (first.y + second.y) * edgeCross ;
        centroidDenominator += edgeCross ;
      }
      if (std::abs(centroidDenominator) <= relativeTolerance * scale * scale)
        return false ;
      const Point2 projectedCentroid(centroidX / (3.0 * centroidDenominator),
            centroidY / (3.0 * centroidDenominator)) ;
      const double normalProjection =
            source.droppedAxis == 0
                  ? std::abs(source.normal.x)
                  : (source.droppedAxis == 1 ? std::abs(source.normal.y)
                                             : std::abs(source.normal.z)) ;
      if (normalProjection <= std::numeric_limits<double>::epsilon())
        return false ;
      overlapArea = projectedArea / normalProjection ;
      overlapCentroid = liftPoint(projectedCentroid, source) ;
      return std::isfinite(overlapArea) && finiteVector(overlapCentroid) ;
    }

    // Exact retention needs no intersection, even for a warped face. Allow
    // the vertex loop to start elsewhere or reverse with cl/cr orientation.
    bool samePolygon(const std::vector<vector3d<double>>& source,
          const std::vector<vector3d<double>>& target, double tolerance) {
      if (source.size() != target.size() || source.empty())
        return false ;
      for (size_t start = 0; start < target.size(); ++start) {
        if (norm(source[0] - target[start]) > tolerance)
          continue ;
        bool forward = true, reverse = true ;
        for (size_t vertex = 1; vertex < source.size(); ++vertex) {
          forward =
                forward &&
                norm(source[vertex] -
                      target[(start + vertex) % target.size()]) <= tolerance ;
          reverse = reverse && norm(source[vertex] -
                                     target[(start + target.size() - vertex) %
                                            target.size()]) <= tolerance ;
        }
        if (forward || reverse)
          return true ;
      }
      return false ;
    }

    bool rootOrder(const RootCellState& left, const RootCellState& right) {
      return left.root < right.root ;
    }

    int rootStateIndex(const std::vector<RootCellState>& states, int root) {
      RootCellState key ;
      key.root = root ;
      const std::vector<RootCellState>::const_iterator location =
            std::lower_bound(states.begin(), states.end(), key, rootOrder) ;
      if (location == states.end() || location->root != root)
        return -1 ;
      return int(location - states.begin()) ;
    }

    int childrenForHexCode(int code) {
      if (code == 1 || code == 2 || code == 4)
        return 2 ;
      if (code == 3 || code == 5 || code == 6)
        return 4 ;
      if (code == 7)
        return 8 ;
      return code == 0 ? 0 : -1 ;
    }

    struct CellTreeTopology {
      std::vector<std::vector<int>> leaves ;
      std::map<std::vector<int>, int> splitCodes ;
      std::map<std::vector<int>, int> splitArities ;
    } ;

    struct PlanNode {
      std::vector<int> path ;
      int edgeCount ;

      PlanNode(const std::vector<int>& nodePath, int nodeEdgeCount)
          : path(nodePath), edgeCount(nodeEdgeCount) {}
    } ;

    int childrenForPrismCode(int code, int edgeCount) {
      if (code == 1)
        return 2 ;
      if (code == 2)
        return edgeCount ;
      if (code == 3)
        return 2 * edgeCount ;
      return code == 0 ? 0 : -1 ;
    }

    bool decodePrismPlan(
          const std::vector<char>& plan, CellTreeTopology& topology) {
      topology = CellTreeTopology() ;
      std::queue<PlanNode> pending ;
      pending.push(PlanNode(std::vector<int>(), 3)) ;
      size_t planIndex = 0 ;
      while (!pending.empty()) {
        const PlanNode node = pending.front() ;
        pending.pop() ;
        const int code = planIndex < plan.size() ? int(plan[planIndex++]) : 0 ;
        const int children = childrenForPrismCode(code, node.edgeCount) ;
        if (children < 0)
          return false ;
        if (children == 0) {
          topology.leaves.push_back(node.path) ;
          continue ;
        }
        topology.splitCodes[node.path] = code ;
        topology.splitArities[node.path] = children ;
        const int childEdgeCount = code == 1 ? node.edgeCount : 4 ;
        for (int child = 0; child < children; ++child) {
          std::vector<int> childPath = node.path ;
          childPath.push_back(code) ;
          childPath.push_back(child) ;
          pending.push(PlanNode(childPath, childEdgeCount)) ;
        }
      }
      for (; planIndex < plan.size(); ++planIndex)
        if (plan[planIndex] != 0)
          return false ;
      return true ;
    }

    bool analyzeCellTree(
          const RootCellState& state, CellTreeTopology& topology) {
      topology = CellTreeTopology() ;
      if (state.leafPaths.empty() ||
            (state.topology != cell_topology::hex &&
                  state.topology != cell_topology::prism &&
                  state.topology != cell_topology::general))
        return false ;

      std::set<std::vector<int>> leaves ;
      std::map<std::vector<int>, std::set<int>> children ;
      std::map<std::vector<int>, int> expectedChildCounts ;
      for (size_t leaf = 0; leaf < state.leafPaths.size(); ++leaf) {
        const std::vector<int>& path = state.leafPaths[leaf] ;
        if (!validRefinementPath(path) || !leaves.insert(path).second)
          return false ;

        std::vector<int> parent ;
        int prismEdgeCount = 3 ;
        for (size_t step = 0; step < path.size(); step += 2) {
          const int code = path[step] ;
          const int child = path[step + 1] ;
          int childCount = -1 ;
          if (state.topology == cell_topology::hex)
            childCount = childrenForHexCode(code) ;
          else if (state.topology == cell_topology::prism)
            childCount = childrenForPrismCode(code, prismEdgeCount) ;
          else if (code == 1)
            childCount = 0 ; // General-cell arity is supplied by the paths.
          if (childCount < 0 || (childCount != 0 && child >= childCount))
            return false ;

          const std::map<std::vector<int>, int>::const_iterator split =
                topology.splitCodes.find(parent) ;
          if (split != topology.splitCodes.end() && split->second != code)
            return false ;
          topology.splitCodes[parent] = code ;
          children[parent].insert(child) ;
          if (childCount != 0)
            expectedChildCounts[parent] = childCount ;

          parent.push_back(code) ;
          parent.push_back(child) ;
          if (state.topology == cell_topology::prism && code != 1)
            prismEdgeCount = 4 ;
        }
      }

      for (std::set<std::vector<int>>::const_iterator leaf = leaves.begin();
            leaf != leaves.end(); ++leaf) {
        std::vector<int> prefix ;
        for (size_t length = 0; length < leaf->size(); length += 2) {
          prefix.assign(leaf->begin(), leaf->begin() + length) ;
          if (leaves.find(prefix) != leaves.end())
            return false ;
        }
      }
      for (std::map<std::vector<int>, std::set<int>>::const_iterator split =
                  children.begin();
            split != children.end(); ++split) {
        const std::map<std::vector<int>, int>::const_iterator count =
              expectedChildCounts.find(split->first) ;
        if (count != expectedChildCounts.end()) {
          if (split->second.size() != size_t(count->second))
            return false ;
          topology.splitArities[split->first] = count->second ;
        } else {
          if (split->second.size() < 2 || *split->second.begin() != 0 ||
                size_t(*split->second.rbegin() + 1) != split->second.size())
            return false ;
          topology.splitArities[split->first] = int(split->second.size()) ;
        }
      }

      topology.leaves.assign(leaves.begin(), leaves.end()) ;
      return true ;
    }

    bool compatibleCellTrees(
          const RootCellState& sourceState, const RootCellState& targetState) {
      if (sourceState.topology != targetState.topology)
        return false ;
      CellTreeTopology source ;
      CellTreeTopology target ;
      if (!analyzeCellTree(sourceState, source) ||
            !analyzeCellTree(targetState, target))
        return false ;
      for (std::map<std::vector<int>, int>::const_iterator split =
                  source.splitCodes.begin();
            split != source.splitCodes.end(); ++split) {
        const std::map<std::vector<int>, int>::const_iterator targetSplit =
              target.splitCodes.find(split->first) ;
        if (targetSplit != target.splitCodes.end()) {
          const std::map<std::vector<int>, int>::const_iterator sourceArity =
                source.splitArities.find(split->first) ;
          const std::map<std::vector<int>, int>::const_iterator targetArity =
                target.splitArities.find(split->first) ;
          if (targetSplit->second != split->second ||
                sourceArity == source.splitArities.end() ||
                targetArity == target.splitArities.end() ||
                sourceArity->second != targetArity->second)
            return false ;
        }
      }
      return true ;
    }

    typedef std::map<int, std::set<std::vector<int>>> RootLeafIndex ;

    bool validFaceKey(const FaceKey& key, const RootLeafIndex& rootLeaves) {
      if (!validFaceKeyStructure(key))
        return false ;
      if (key.origin == face_origin::base_face)
        return true ;

      const RootLeafIndex::const_iterator root = rootLeaves.find(key.root) ;
      return root != rootLeaves.end() &&
             root->second.find(key.firstPath) != root->second.end() &&
             root->second.find(key.secondPath) != root->second.end() ;
    }

    bool pathPrefix(
          const std::vector<int>& prefix, const std::vector<int>& path) {
      return prefix.size() <= path.size() &&
             std::equal(prefix.begin(), prefix.end(), path.begin()) ;
    }

    bool relatedPaths(
          const std::vector<int>& first, const std::vector<int>& second) {
      return pathPrefix(first, second) || pathPrefix(second, first) ;
    }

    bool relatedFaces(const FaceKey& first, const FaceKey& second) {
      if (first.origin != second.origin || first.root != second.root)
        return false ;
      if (first.origin == face_origin::base_face)
        return relatedPaths(first.firstPath, second.firstPath) ;
      return (relatedPaths(first.firstPath, second.firstPath) &&
                   relatedPaths(first.secondPath, second.secondPath)) ||
             (relatedPaths(first.firstPath, second.secondPath) &&
                   relatedPaths(first.secondPath, second.firstPath)) ;
    }

    typedef std::map<std::vector<int>, std::vector<size_t>> FacePathIndex ;

    struct FaceGroupIndex {
      FacePathIndex firstPaths ;
      FacePathIndex secondPaths ;
    } ;

    void addRelatedFaces(const std::vector<int>& path,
          const FacePathIndex& index, std::set<size_t>& faces) {
      std::vector<int> prefix ;
      for (size_t length = 0; length < path.size(); length += 2) {
        prefix.assign(path.begin(), path.begin() + length) ;
        const FacePathIndex::const_iterator ancestor = index.find(prefix) ;
        if (ancestor != index.end())
          faces.insert(ancestor->second.begin(), ancestor->second.end()) ;
      }

      FacePathIndex::const_iterator descendant = index.lower_bound(path) ;
      while (descendant != index.end() && pathPrefix(path, descendant->first)) {
        faces.insert(descendant->second.begin(), descendant->second.end()) ;
        ++descendant ;
      }
    }

    bool containingLeafCell(const RootCellState& state,
          const FaceKey& interiorFace, CellId& cell) {
      int containingLeaf = -1 ;
      for (size_t leaf = 0; leaf < state.leafPaths.size(); ++leaf) {
        if (!pathPrefix(state.leafPaths[leaf], interiorFace.firstPath) ||
              !pathPrefix(state.leafPaths[leaf], interiorFace.secondPath))
          continue ;
        if (containingLeaf >= 0)
          return false ;
        containingLeaf = int(leaf) ;
      }
      if (containingLeaf < 0)
        return false ;
      cell = persistentCellId(state.root, state.leafPaths[containingLeaf]) ;
      return true ;
    }

    bool encodePaths(const std::vector<std::vector<int>>& leaves,
          std::vector<int>& encodedPaths) {
      encodedPaths.clear() ;
      encodedPaths.push_back(int(leaves.size())) ;
      for (size_t leaf = 0; leaf < leaves.size(); ++leaf) {
        if (!validRefinementPath(leaves[leaf])) {
          encodedPaths.clear() ;
          return false ;
        }
        encodedPaths.push_back(int(leaves[leaf].size() / 2)) ;
        encodedPaths.insert(
              encodedPaths.end(), leaves[leaf].begin(), leaves[leaf].end()) ;
      }
      return true ;
    }

    bool encodePlanLeafPaths(const std::vector<char>& plan, bool quadrilateral,
          std::vector<int>& encodedPaths) {
      std::queue<std::vector<int>> pending ;
      std::vector<std::vector<int>> leaves ;
      pending.push(std::vector<int>()) ;
      size_t planIndex = 0 ;
      while (!pending.empty()) {
        const std::vector<int> path = pending.front() ;
        pending.pop() ;
        const int code = planIndex < plan.size() ? int(plan[planIndex++]) : 0 ;
        const int children =
              quadrilateral ? (code == 0 ? 0
                                         : (code == 1 || code == 2
                                                       ? 2
                                                       : (code == 3 ? 4 : -1)))
                            : childrenForHexCode(code) ;
        if (children < 0)
          return false ;
        if (children == 0) {
          leaves.push_back(path) ;
          continue ;
        }
        for (int child = 0; child < children; ++child) {
          std::vector<int> childPath = path ;
          childPath.push_back(code) ;
          childPath.push_back(child) ;
          pending.push(childPath) ;
        }
      }
      for (; planIndex < plan.size(); ++planIndex)
        if (plan[planIndex] != 0)
          return false ;

      return encodePaths(leaves, encodedPaths) ;
    }

    bool encodeGeneralFacePaths(const std::vector<char>& plan,
          int initialEdgeCount, std::vector<int>& encodedPaths) {
      if (initialEdgeCount < 3)
        return false ;
      std::queue<PlanNode> pending ;
      std::vector<std::vector<int>> leaves ;
      pending.push(PlanNode(std::vector<int>(), initialEdgeCount)) ;
      size_t planIndex = 0 ;
      while (!pending.empty()) {
        const PlanNode node = pending.front() ;
        pending.pop() ;
        const int code = planIndex < plan.size() ? int(plan[planIndex++]) : 0 ;
        if (code == 0) {
          leaves.push_back(node.path) ;
          continue ;
        }
        if (code != 1)
          return false ;
        for (int child = 0; child < node.edgeCount; ++child) {
          std::vector<int> childPath = node.path ;
          childPath.push_back(code) ;
          childPath.push_back(child) ;
          pending.push(PlanNode(childPath, 4)) ;
        }
      }
      for (; planIndex < plan.size(); ++planIndex)
        if (plan[planIndex] != 0)
          return false ;
      return encodePaths(leaves, encodedPaths) ;
    }
  }

  CPTR<FaceState> FaceState::create(const std::vector<FaceIdentity>& identities,
        const std::vector<std::vector<vector3d<double>>>& polygons,
        const std::vector<RootCellState>& rootCells,
        FaceTransitionReport& report, double relativeTolerance) {
    report = FaceTransitionReport() ;
    report.targetFaces = identities.size() ;
    if (identities.size() != polygons.size() || relativeTolerance < 0.0 ||
          !std::isfinite(relativeTolerance)) {
      report.invalidIdentities++ ;
      report.status = face_transition_status::invalid_identity ;
      return CPTR<FaceState>() ;
    }

    std::set<int> roots ;
    RootLeafIndex rootLeaves ;
    for (size_t root = 0; root < rootCells.size(); ++root) {
      CellTreeTopology topology ;
      if (!roots.insert(rootCells[root].root).second ||
            rootCells[root].root < 0 ||
            !analyzeCellTree(rootCells[root], topology)) {
        report.invalidIdentities++ ;
        continue ;
      }
      rootLeaves[rootCells[root].root] = std::set<std::vector<int>>(
            rootCells[root].leafPaths.begin(), rootCells[root].leafPaths.end()) ;
    }
    std::set<int> faceNumbers ;
    std::set<FaceId> faceIds ;
    std::set<FaceKey> keys ;
    for (size_t face = 0; face < identities.size(); ++face) {
      if (!validFaceKey(identities[face].key, rootLeaves) ||
            identities[face].id == 0 ||
            identities[face].id != persistentFaceId(identities[face].key) ||
            !faceNumbers.insert(identities[face].face).second ||
            !faceIds.insert(identities[face].id).second ||
            !keys.insert(identities[face].key).second)
        report.invalidIdentities++ ;
      PolygonGeometry geometry ;
      if (!polygonGeometry(polygons[face], relativeTolerance, geometry,
                identities[face].key.origin == face_origin::cell_interior))
        report.invalidPolygons++ ;
    }
    if (report.invalidIdentities != 0 || report.invalidPolygons != 0) {
      report.status = report.invalidIdentities != 0
                            ? face_transition_status::invalid_identity
                            : face_transition_status::invalid_geometry ;
      return CPTR<FaceState>() ;
    }

    CPTR<FaceState> state = new FaceState ;
    state->identities_ = identities ;
    state->polygons_ = polygons ;
    state->rootCells_ = rootCells ;
    for (size_t root = 0; root < state->rootCells_.size(); ++root)
      std::sort(state->rootCells_[root].leafPaths.begin(),
            state->rootCells_[root].leafPaths.end()) ;
    std::sort(state->rootCells_.begin(), state->rootCells_.end(), rootOrder) ;
    report.valid = true ;
    report.status = face_transition_status::available ;
    return state ;
  }

  CPTR<FaceRemap> buildFaceRemap(const_CPTR<FaceState> source,
        const_CPTR<FaceState> target, FaceTransitionReport& report,
        double relativeTolerance) {
    report = FaceTransitionReport() ;
    if (source == static_cast<FaceState*>(0) ||
          target == static_cast<FaceState*>(0)) {
      report.invalidIdentities++ ;
      return CPTR<FaceRemap>() ;
    }
    report.sourceFaces = source->identities_.size() ;
    report.targetFaces = target->identities_.size() ;

    if (source->rootCells_.size() != target->rootCells_.size()) {
      report.invalidIdentities++ ;
      return CPTR<FaceRemap>() ;
    }
    for (size_t root = 0; root < source->rootCells_.size(); ++root) {
      const RootCellState& previous = source->rootCells_[root] ;
      const RootCellState& accepted = target->rootCells_[root] ;
      if (previous.root != accepted.root ||
            previous.topology != accepted.topology) {
        report.invalidIdentities++ ;
        continue ;
      }
      if (!compatibleCellTrees(previous, accepted))
        report.unsupportedRootPlans++ ;
    }
    if (report.invalidIdentities != 0) {
      report.status = face_transition_status::invalid_identity ;
      return CPTR<FaceRemap>() ;
    }
    if (report.unsupportedRootPlans != 0) {
      report.status = face_transition_status::unsupported_plan_change ;
      return CPTR<FaceRemap>() ;
    }

    std::vector<PolygonGeometry> sourceGeometry(source->polygons_.size()) ;
    std::vector<PolygonGeometry> targetGeometry(target->polygons_.size()) ;
    for (size_t face = 0; face < source->polygons_.size(); ++face)
      if (!polygonGeometry(source->polygons_[face], relativeTolerance,
                sourceGeometry[face],
                source->identities_[face].key.origin ==
                      face_origin::cell_interior))
        report.invalidPolygons++ ;
    for (size_t face = 0; face < target->polygons_.size(); ++face)
      if (!polygonGeometry(target->polygons_[face], relativeTolerance,
                targetGeometry[face],
                target->identities_[face].key.origin ==
                      face_origin::cell_interior))
        report.invalidPolygons++ ;
    if (report.invalidPolygons != 0) {
      report.status = face_transition_status::invalid_geometry ;
      return CPTR<FaceRemap>() ;
    }

    std::vector<FaceGeometry> sourceFaces ;
    std::vector<FaceGeometry> targetFaces ;
    for (size_t face = 0; face < source->identities_.size(); ++face)
      sourceFaces.push_back(FaceGeometry(source->identities_[face].id,
            sourceGeometry[face].area, sourceGeometry[face].centroid)) ;
    for (size_t face = 0; face < target->identities_.size(); ++face)
      targetFaces.push_back(FaceGeometry(target->identities_[face].id,
            targetGeometry[face].area, targetGeometry[face].centroid)) ;

    std::vector<FaceOverlap> contributions ;
    std::vector<size_t> sourceDegree(source->identities_.size(), 0) ;
    std::vector<size_t> targetDegree(target->identities_.size(), 0) ;
    typedef std::pair<int, int> FaceGroup ;
    std::map<FaceGroup, FaceGroupIndex> targetGroups ;
    for (size_t targetFace = 0; targetFace < target->identities_.size();
          ++targetFace) {
      const FaceKey& key = target->identities_[targetFace].key ;
      FaceGroupIndex& group =
            targetGroups[FaceGroup(static_cast<int>(key.origin), key.root)] ;
      group.firstPaths[key.firstPath].push_back(targetFace) ;
      if (key.origin == face_origin::cell_interior)
        group.secondPaths[key.secondPath].push_back(targetFace) ;
    }
    for (size_t sourceFace = 0; sourceFace < source->identities_.size();
          ++sourceFace) {
      const FaceKey& sourceKey = source->identities_[sourceFace].key ;
      const std::map<FaceGroup, FaceGroupIndex>::const_iterator group =
            targetGroups.find(FaceGroup(
                  static_cast<int>(sourceKey.origin), sourceKey.root)) ;
      if (group == targetGroups.end())
        continue ;
      std::set<size_t> candidates ;
      addRelatedFaces(
            sourceKey.firstPath, group->second.firstPaths, candidates) ;
      if (sourceKey.origin == face_origin::cell_interior)
        addRelatedFaces(
              sourceKey.firstPath, group->second.secondPaths, candidates) ;
      for (std::set<size_t>::const_iterator candidate = candidates.begin();
            candidate != candidates.end(); ++candidate) {
        const size_t targetFace = *candidate ;
        if (!relatedFaces(sourceKey, target->identities_[targetFace].key))
          continue ;
        double area = 0.0 ;
        vector3d<double> centroid ;
        int orientation = 1 ;
        if ((!sourceGeometry[sourceFace].planar ||
                  !targetGeometry[targetFace].planar) &&
              sourceKey == target->identities_[targetFace].key &&
              samePolygon(source->polygons_[sourceFace],
                    target->polygons_[targetFace],
                    relativeTolerance *
                          std::max(sourceGeometry[sourceFace].scale,
                                targetGeometry[targetFace].scale))) {
          area = targetGeometry[targetFace].area ;
          centroid = targetGeometry[targetFace].centroid ;
          orientation = dot(sourceGeometry[sourceFace].normal,
                              targetGeometry[targetFace].normal) >= 0.0
                              ? 1
                              : -1 ;
        } else if (!overlapGeometry(sourceGeometry[sourceFace],
                         target->polygons_[targetFace],
                         targetGeometry[targetFace], relativeTolerance, area,
                         centroid, orientation)) {
          report.invalidPolygons++ ;
          continue ;
        }
        if (area <= 0.0)
          continue ;
        contributions.push_back(FaceOverlap(source->identities_[sourceFace].id,
              target->identities_[targetFace].id, area, centroid, orientation)) ;
        sourceDegree[sourceFace]++ ;
        targetDegree[targetFace]++ ;
      }
    }
    if (report.invalidPolygons != 0) {
      report.status = face_transition_status::invalid_geometry ;
      return CPTR<FaceRemap>() ;
    }

    std::vector<CreatedFace> created ;
    std::vector<RemovedFace> removed ;
    for (size_t targetFace = 0; targetFace < target->identities_.size();
          ++targetFace) {
      if (targetDegree[targetFace] != 0)
        continue ;
      const FaceKey& key = target->identities_[targetFace].key ;
      const int previousRoot = rootStateIndex(source->rootCells_, key.root) ;
      CellId sourceCell = 0 ;
      if (key.origin != face_origin::cell_interior || previousRoot < 0 ||
            !containingLeafCell(
                  source->rootCells_[previousRoot], key, sourceCell)) {
        report.status = face_transition_status::inconsistent_relation ;
        return CPTR<FaceRemap>() ;
      }
      created.push_back(
            CreatedFace(target->identities_[targetFace].id, sourceCell)) ;
    }
    for (size_t sourceFace = 0; sourceFace < source->identities_.size();
          ++sourceFace) {
      if (sourceDegree[sourceFace] != 0)
        continue ;
      const FaceKey& key = source->identities_[sourceFace].key ;
      const int acceptedRoot = rootStateIndex(target->rootCells_, key.root) ;
      CellId targetCell = 0 ;
      if (key.origin != face_origin::cell_interior || acceptedRoot < 0 ||
            !containingLeafCell(
                  target->rootCells_[acceptedRoot], key, targetCell)) {
        report.status = face_transition_status::inconsistent_relation ;
        return CPTR<FaceRemap>() ;
      }
      removed.push_back(
            RemovedFace(source->identities_[sourceFace].id, targetCell)) ;
    }

    CPTR<FaceRemap> remap = FaceRemap::create(sourceFaces, targetFaces,
          contributions, created, removed, report.remap, relativeTolerance) ;
    if (remap == static_cast<FaceRemap*>(0)) {
      report.status = face_transition_status::inconsistent_relation ;
      return CPTR<FaceRemap>() ;
    }

    report.valid = true ;
    report.status = face_transition_status::available ;
    return remap ;
  }

  namespace detail {
    bool encodeLeafPaths(const std::vector<std::vector<int>>& paths,
          std::vector<int>& encodedPaths) {
      return encodePaths(paths, encodedPaths) ;
    }

    bool encodeQuadFaceLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) {
      return encodePlanLeafPaths(plan, true, encodedPaths) ;
    }

    bool encodeGeneralFaceLeafPaths(const std::vector<char>& plan,
          int initialEdgeCount, std::vector<int>& encodedPaths) {
      return encodeGeneralFacePaths(plan, initialEdgeCount, encodedPaths) ;
    }

    bool encodeHexCellLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) {
      return encodePlanLeafPaths(plan, false, encodedPaths) ;
    }

    bool encodePrismCellLeafPaths(
          const std::vector<char>& plan, std::vector<int>& encodedPaths) {
      CellTreeTopology topology ;
      if (!decodePrismPlan(plan, topology))
        return false ;
      return encodePaths(topology.leaves, encodedPaths) ;
    }

    bool decodeLeafPaths(const std::vector<int>& encodedPaths,
          std::vector<std::vector<int>>& paths) {
      paths.clear() ;
      if (encodedPaths.empty() || encodedPaths[0] < 0 ||
            size_t(encodedPaths[0]) > encodedPaths.size() - 1)
        return false ;
      size_t entry = 1 ;
      paths.resize(size_t(encodedPaths[0])) ;
      for (size_t path = 0; path < paths.size(); ++path) {
        if (entry >= encodedPaths.size() || encodedPaths[entry] < 0)
          return false ;
        const size_t steps = size_t(encodedPaths[entry++]) ;
        if (steps > (encodedPaths.size() - entry) / 2)
          return false ;
        const size_t length = 2 * steps ;
        paths[path].assign(encodedPaths.begin() + entry,
              encodedPaths.begin() + entry + length) ;
        entry += length ;
      }
      return entry == encodedPaths.size() ;
    }
  }
}

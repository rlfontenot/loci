//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#include "kd_tree.h"
#include <iostream>
using std::endl ;
using std::cout ;
using std::cerr ;

using std::vector ;
using std::pair ;
using std::min ;
using std::max ;
using std::swap ;
using std::sort ;
using std::numeric_limits ;
namespace Loci {
  namespace kdTree {
    // Constructor
    kd_tree::kd_tree(const vector<coord3d> &inpnts,
                     const vector<int> &ids) {
      // allocate space for pnts, copy input values into data structure
      pnts.reserve(inpnts.size()) ;

      for(int d=0;d<3;++d) {
        bbox.minc[d] =  std::numeric_limits<float>::max() ;
        bbox.minc[0] = -std::numeric_limits<float>::max() ;
      }

      for(size_t i=0;i<inpnts.size();++i) {
        coord_info ci ;
        ci.coords = inpnts[i] ;
        ci.id = ids[i] ;
        for(int d=0;d<3;++d) {
          bbox.minc[d] = min(bbox.minc[d],ci.coords[d]) ;
          bbox.maxc[d] = max(bbox.maxc[d],ci.coords[d]) ;
        }
        pnts.push_back(ci) ;
      }
      // allocate space for the splits
      vector<int> tmp(inpnts.size(),-1) ;
      splits.swap(tmp) ;

      // Recursively build the tree
      build_kd(0,pnts.size(),0) ;
    }
    // Constructor
    kd_tree::kd_tree(vector<coord_info> &inpnts) {
      // allocate space for pnts, copy input values into data structure
      pnts.swap(inpnts) ;

      // Compute bounding box
      for(int d=0;d<3;++d) {
        bbox.minc[d] =  std::numeric_limits<float>::max() ;
        bbox.minc[0] = -std::numeric_limits<float>::max() ;
      }
      for(size_t i=0;i<pnts.size();++i) {
        for(int d=0;d<3;++d) {
          bbox.minc[d] = min(bbox.minc[d],pnts[i].coords[d]) ;
          bbox.maxc[d] = max(bbox.maxc[d],pnts[i].coords[d]) ;
        }
      }

      // allocate space for the splits
      vector<int> tmp(pnts.size(),-1) ;
      splits.swap(tmp) ;

      // Recursively build the tree
      build_kd(0,pnts.size(),0) ;
    }

    namespace {
      // A tool used for the median computation
      inline bool pair_sort(const std::pair<float,int> &p1,
                            const std::pair<float,int> &p2) {
        return p1.first < p2.first ;
      }
    }

    // Select a pivot and move to begining of array.
    void kd_tree::select_pivot(int start, int end, int coord) {
      const int sz = (end-start) ;
      if(sz < 3) return ;
    

      // Choose a pivot by sampling 10000 points.  If size is less than 10000,
      // sample half
      int factor = std::max(sz/10000,2) ;
      int npnts = sz/factor ;
      int remain = sz%factor ;

      // select a sample and sort it
      std::vector<std::pair<float,int> > pv ;
      pv.reserve(npnts) ;
      for(int i=start+remain/2;i<end;i+=factor) 
        pv.push_back(std::pair<float,int>(pnts[i].coords[coord],i)) ;
      sort(pv.begin(),pv.end(),pair_sort) ;

      // Move pivot to first element
      int id = pv[(pv.size()-1)/2].second ;
      std::swap(pnts[start],pnts[id]) ;
    }


    // Recursively build the kd_tree.  start and end define the interval
    // of the array that will be inserted into the tree.  depth is the current
    // depth that starts with 0 at the root.  Coordinates that are split cycle
    // with the depth of the tree.
    void kd_tree::build_kd(int start, int end, int depth) {
      const int sz = end-start ;
      if(sz <= kd_bin_size) // Were at a leaf node, so return
        return ;

      // Current coordinate for bisection (0==x,1==y,2==z)
      const int coord = depth%3 ;

      // selects pivot and moves it to start
      select_pivot(start,end,coord) ;

      // Pivot is at start, move to pval for convenience
      double pval = pnts[start].coords[coord] ;

      // partition remaining array based on pivot value
      int p1 = start+1 ; // starting point for partition
      int p2 = end-1 ;   // last point to partition

      // Now partition into two parts: 1 part <= than pivot, other part > than pivot
      while(p1 <= p2) {
        // Note start+1 to p1 is part 1, while p2 to end is part 2
        if(pnts[p1].coords[coord] <= pval) // If less than pivot,
          p1++ ; // leave in first list by moving first pointer forward
        else { // otherwise move to the end of the list
          while(pnts[p2].coords[coord] > pval && p1 < p2)
            p2-- ;
          swap(pnts[p1],pnts[p2]) ;
          p2-- ; // move end of list pointer 
        }
      }
      // Record how remaining array is split
      int split = p1 ;
      splits[start] = split ;

#ifdef DEBUG
      // Verify split
      for(int i=start+1;i<split;++i)
        if(pnts[i].coords[coord] > pval)
          std::cerr << "split error" << std::endl ;
      for(int i=split;i<end;++i)
        if(pnts[i].coords[coord] <= pval)
          std::cerr << "split error" << std::endl ;
#endif
      
      // Construct left and right tree branches recursively
      build_kd(start+1, split, depth+1) ;
      build_kd(split,   end,   depth+1) ;
    }

    // Find depth of tree
    int kd_tree::depth_search(int start,int end,int depth) const {
      const int sz = end-start ;
      if(sz <= kd_bin_size) // Were done now, so return
        return depth ;
      const int split = splits[start] ;
      // The depth of the tree at any point is the maximum of the two branch depths
      return max(depth_search(start+1, split, depth+1),
                 depth_search(split,   end,   depth+1)) ;
    }


    // Here is a bit of code that checks to see if the current kd tree bounding box
    // contains the sphere formed by the current search space and the currently
    // known closest point radius.  It basically checks to see if the sphere
    // intersects any of the bounding box sides. (Note, we have a priority to
    // search the kd tree bounding box that contains the center, so this test
    // is only applied for bouding boxes for which the center is outside, but
    // the radius may project into the bounding box.
    bool kd_tree::insphere(const coord3d &v,double r, double proj_d,
                           bounds bnd, int coord) const {
      if(proj_d > r) // sphere not touching dividing plane, so test fails
        return false ;
      // Now check if circle in bounding plane is in bounds
      const int c1 = (coord+1)%3 ; 
      const int c2 = (coord+2)%3 ;
      // Transform to a two dimensional coordinate system that represents
      // the current cutting plane
      double x = v[c1] ;
      double y = v[c2] ;
      double xl = bnd.minc[c1] ;
      double xh = bnd.maxc[c1] ;
      double yl = bnd.minc[c2] ;
      double yh = bnd.maxc[c2] ;
      // The sphere projected onto this plane forms a circle centered at x,y
      // check to see if this center is within the bounding face, if so
      // the sphere intersects this bounding box
      if(xl <= x && x <= xh && yl <= y && y <= yh) // center of circle in bounds
        return true ;
      // Now compute the radius of the circle in the cut plane using
      // the (pathagorean theorem)
      double r2 = r-proj_d ;
      // Compute smallest Linf distance to bounding box corners

      // We will determine if this circle intersects the bounding box by just
      // comparing bounding coordinates
      double min_x = std::min((x-xl)*(x-xl),(x-xh)*(x-xh)) ;
      double min_y = std::min((y-yl)*(y-yl),(y-yh)*(y-yh)) ;
    
      double min_dist = std::min(min_x,min_y) ;

      // square circle intersection test
      if(min_dist > r2) { // If box containing circle is outside of bounding box
        return false ;    // then not in sphere
      }
      if(min_x <= r2) { // Check to see if circle intersects bbox x=const lines
        if(y >= yl && y <= yh)
          return true ;
        double r3 = r2-min_x ;
        if(min_y <= r3)
          return true ;
      }
      if(min_y <= r2) { // Check to see if circle intersects bbox y=const lines
        if(x >= xl && x <= xh)
          return true ;
        double r3 = r2-min_y ;
        if(min_x <= r3)
          return true ;
      }
      // Sphere does not intersect bbox
      return false ;
    }

    // Search for the closest corresponding point to the vector defined by v
    // rmin is the current best known estimate for the distance to the closest
    // point.  bnds is the current bounding box.
    int kd_tree::find_closest(int start,int end, int depth, const coord3d &v,
                              double &rmin,bounds bnds) const {
  
      const int sz = end-start ;

      // In base case perform linear search to find closest point
      if(sz <= kd_bin_size)  {
        if(sz == 0) // If size is zero, just return the root pivot
          return -1 ;
        int pt = start ;
        double rmint = dist2(v,pnts[pt].coords) ;
        for(int i=start+1;i<end;++i) {
          double rtmp = dist2(v,pnts[i].coords) ;
          if(rmint > rtmp) {
            rmint = rtmp ;
            pt = i ;
          }
        }
        // pt is the closest point in this group
        // Update rmin if this is closer than anything currently found
        if(rmint<rmin) {
          rmin = rmint ;
          return pt ;
        }
        return -1 ;
      }

      const int coord = depth%3 ;

      // otherwise search for which leaf contains our point
      int id = start ;

      // Check to see if the partitioning node is closer than current closest find
      double rtmp = dist2(v,pnts[id].coords) ;
      if(rtmp < rmin) 
        rmin = rtmp ;
      else
        id = -1 ;

      // Find the projected distance to the current cutting plane
      double proj_d = v[coord]-pnts[start].coords[coord] ;
      proj_d = proj_d*proj_d ; // square since we are comparing squared distances

      // is the point we are searching for on the left side of this partition?
      // we want to search in the containing partition first so that we
      // get a good estimate of the closest radius as early as possible
      const bool left = (v[coord] <= pnts[start].coords[coord]) ;
      // Compute the bounding box for the left and right kd-tree partitions
      bounds bnds_left = bnds ;
      bounds bnds_right = bnds ;
      bnds_left.maxc[coord] = pnts[start].coords[coord] ;
      bnds_right.minc[coord] = pnts[start].coords[coord] ;

      if(left) { // traverse left side first
        double proj_d2 = v[coord]-bnds_left.minc[coord] ;
        proj_d2 = proj_d2*proj_d2 ;

        // Check if sphere is in left side bounding box
        if(v[coord] >= bnds_left.minc[coord] ||
           insphere(v,rmin,proj_d2,bnds_left,coord)) {
          int id1= find_closest(start+1,splits[start], depth+1, v, rmin, bnds_left) ;
          // If a closer point is found, update id.
          if(id1 >= 0)
            id = id1 ;
        }
        // Now check the right side bounding box
        if(insphere(v,rmin,proj_d,bnds_right,coord)) {
          // If so, search this side for closest point
          int id2= find_closest(splits[start],end, depth+1, v, rmin, bnds_right) ;
          // If we find a closer point, update
          if(id2 >= 0)
            id = id2 ;
        }
      } else { // Search right side first
        // Do same thing as left side just search right side of tree first
        double proj_d2 = v[coord]-bnds_right.maxc[coord] ;
        proj_d2 = proj_d2*proj_d2 ;
        if(v[coord] <= bnds_right.maxc[coord] ||
           insphere(v,rmin,proj_d2,bnds_right,coord)) {
          int id1 = find_closest(splits[start],end, depth+1, v, rmin, bnds_right) ;
          if(id1 >= 0)
            id = id1 ;
        }
        if(insphere(v,rmin,proj_d,bnds_left,coord)) {
          int id2= find_closest(start+1,splits[start], depth+1, v, rmin, bnds_left) ;
          if(id2 >= 0)
            id = id2 ;
        }
      }
      // Return the id of the current closest point
      return id ;
    }

    inline bool pt_in_box(coord3d v,const kd_tree::bounds &b1) {
      return
        (v[0] >= b1.minc[0] && v[0] <= b1.maxc[0]) &&
        (v[1] >= b1.minc[1] && v[1] <= b1.maxc[1]) &&
        (v[2] >= b1.minc[2] && v[2] <= b1.maxc[2]) ;
    }
    inline bool pt_in_boxb(coord3d v,const kd_tree::bounds &b1) {
      return
        (v[0] >= b1.minc[0] && v[0] < b1.maxc[0]) &&
        (v[1] >= b1.minc[1] && v[1] < b1.maxc[1]) &&
        (v[2] >= b1.minc[2] && v[2] < b1.maxc[2]) ;
    }
    
    inline coord3d pselect(const kd_tree::bounds &b, int sel) {
      coord3d v ;
      v[0] = ((sel&1)==0)?b.maxc[0]:b.minc[0] ;
      v[1] = ((sel&2)==0)?b.maxc[1]:b.minc[1] ;
      v[2] = ((sel&4)==0)?b.maxc[2]:b.minc[2] ;
      return v ;
    }
    
    bool box_intersect(const kd_tree::bounds &b1, const kd_tree::bounds &b2) {
      bool test = false ;
      for(int i=0;i<8;++i)
        test = test || pt_in_box(pselect(b1,i),b2) || pt_in_box(pselect(b2,i),b1) ;
      return test ;
    }
    // recursively search for points in a bounding box, add them to the
    // vector passed in.
    void kd_tree::find_box(int start, int end, int depth,
                           vector<coord_info> &found_pts,
                           const bounds &box, bounds bnds) const {
      const int sz = end-start ;

      // Check if the entire box is in bounds then add all of the current
      // points...

      if(box.minc[0] >= bnds.minc[0] &&
         box.maxc[0] <= bnds.maxc[0] &&
         box.minc[1] >= bnds.minc[1] &&
         box.maxc[1] <= bnds.maxc[1] &&
         box.minc[2] >= bnds.minc[2] &&
         box.maxc[2] <= bnds.maxc[2]) {
        for(int i=start;i<end;++i) 
          found_pts.push_back(pnts[i]);

        return ;
      }
         

      // In base case perform linear search to find closest point
      if(sz <= kd_bin_size)  {
        if(sz == 0) // If size is zero, just return the root pivot
          return ;
        // search points and add them if they are in the box
        for(int i=start;i<end;++i)
          if(pt_in_box(pnts[i].coords,box))
            found_pts.push_back(pnts[i]);
        return ;
      }

      const int coord = depth%3 ;

      // otherwise search for which leaf contains our point
      int id = start ;

      if(pt_in_box(pnts[id].coords,box))
        found_pts.push_back(pnts[id]);
         
      // Compute the bounding box for the left and right kd-tree partitions
      bounds bnds_left = bnds ;
      bounds bnds_right = bnds ;
      bnds_left.maxc[coord] = pnts[start].coords[coord] ;
      bnds_right.minc[coord] = pnts[start].coords[coord] ;

      if(box.minc[coord] <= bnds_left.maxc[coord])
        find_box(start+1,splits[start], depth+1, found_pts, box, bnds_left) ;
      if(box.maxc[coord] >= bnds_right.minc[coord])
        find_box(splits[start],end, depth+1, found_pts, box, bnds_right) ;
    }

    // Search for the closest corresponding point to the vector defined by v
    // rmin is the current best known estimate for the distance to the closest
    // point.  bnds is the current bounding box.
    int kd_tree::find_closest_box(int start,int end, int depth,
                                  const coord3d &v, double &rmin,
                                  const bounds &box,
                                  bounds bnds) const {
  
      const int sz = end-start ;

      // In base case perform linear search to find closest point
      if(sz <= kd_bin_size)  {
        if(sz == 0) // If size is zero, just return the root pivot
          return -1 ;
        int pt = -1 ;
        // search points and add them if they are in the box
        for(int i=start;i<end;++i)
          if(pt_in_boxb(pnts[i].coords,box)) {
            double rtmp = dist2(v,pnts[i].coords) ;
            if(rtmp < rmin) {
              pt = i ;
              rmin = rtmp ;
            }
          }
        return pt ;
      }

      const int coord = depth%3 ;

      // otherwise search for which leaf contains our point
      int id = -1 ;
      // Check to see if the partitioning node is closer than current closest find
      if(pt_in_boxb(pnts[start].coords,box)) {
        double rtmp = dist2(v,pnts[start].coords) ;
        if(rtmp < rmin) {
          rmin = rtmp ;
          id = start ;
        }
      }

      // Find the projected distance to the current cutting plane
      double proj_d = v[coord]-pnts[start].coords[coord] ;
      proj_d = proj_d*proj_d ; // square since we are comparing squared distances

      // is the point we are searching for on the left side of this partition?
      // we want to search in the containing partition first so that we
      // get a good estimate of the closest radius as early as possible
      const bool left = (v[coord] <= pnts[start].coords[coord]) ;
      // Compute the bounding box for the left and right kd-tree partitions
      bounds bnds_left = bnds ;
      bounds bnds_right = bnds ;
      bnds_left.maxc[coord] = pnts[start].coords[coord] ;
      bnds_right.minc[coord] = pnts[start].coords[coord] ;

      if(left) { // traverse left side first
        double proj_d2 = v[coord]-bnds_left.minc[coord] ;
        proj_d2 = proj_d2*proj_d2 ;

        // Check if sphere is in left side bounding box
        if(box.minc[coord] <= bnds_left.maxc[coord]) {
          if(v[coord] >= bnds_left.minc[coord] ||
             insphere(v,rmin,proj_d2,bnds_left,coord)) {
            int id1= find_closest_box(start+1,splits[start], depth+1, v, rmin, box,bnds_left) ;
            // If a closer point is found, update id.
            if(id1 >= 0)
              id = id1 ;
          }
        }
        // Now check the right side bounding box
        if(box.maxc[coord] >= bnds_right.minc[coord]) {
          if(insphere(v,rmin,proj_d,bnds_right,coord)) {
            // If so, search this side for closest point
            int id2= find_closest_box(splits[start],end, depth+1, v, rmin, box,bnds_right) ;
            // If we find a closer point, update
            if(id2 >= 0)
              id = id2 ;
          }
        }
      } else { // Search right side first
        // Do same thing as left side just search right side of tree first
        double proj_d2 = v[coord]-bnds_right.maxc[coord] ;
        proj_d2 = proj_d2*proj_d2 ;
        if(box.maxc[coord] >= bnds_right.minc[coord]) {
          if(v[coord] <= bnds_right.maxc[coord] ||
             insphere(v,rmin,proj_d2,bnds_right,coord)) {
            int id1 = find_closest_box(splits[start],end, depth+1, v, rmin, box, bnds_right) ;
            if(id1 >= 0)
              id = id1 ;
          }
        }
        if(box.minc[coord] <= bnds_left.maxc[coord]) {
          if(insphere(v,rmin,proj_d,bnds_left,coord)) {
            int id2= find_closest_box(start+1,splits[start], depth+1, v, rmin, box, bnds_left) ;
            if(id2 >= 0)
              id = id2 ;
          }
        }
      }
      // Return the id of the current closest point
      return id ;
    }
 
  }
}

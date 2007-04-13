#ifndef KDTREE_H
#define KDTREE_H
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Loci {
  namespace kdTree {
    // 3-D vector class
    struct coord3d {
      double coords[3] ;
      double &operator[](size_t i) {return coords[i]; }
      double operator[](size_t i) const { return coords[i] ; }
      coord3d() {}
      coord3d(double x, double y, double z) {coords[0]=x;coords[1]=y,coords[2]=z;}
    } ;

    // Euclidean Distance Function

    // Compute the squared distance between two points
    inline double dist2(const coord3d &v1,const coord3d &v2) {
      const double vd[3] = {v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]} ;
      return vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2] ;
    }


    // Minimum partitioning size for kd_tree (e.g. the maximum number of nodes
    // found in a leaf node.
    const int kd_bin_size = 32 ;
    //const int kd_bin_size = 2 ;

    // The kd_tree is stored in an implicit tree datastructure.  Points that
    // become pivots are stored in the tree structure, only non-pivot points
    // are stored in the leaves.  The entire tree is stored in the two arrays
    // pnts and splits.  The first point in this array is the first pivot,
    // The two branches of the tree are stored in the remaining array. The left
    // half immediately follows the pivot, and the right half starts at the
    // array index stored in split.  Once the size of a sub array is smaller than
    // kd_bin_size, the remaining array becomes a leaf of the tree.
    class kd_tree {
    public:
      // Coordinate Info for kd-tree
      struct coord_info {
        coord3d coords ;
        int id ;
      } ;


      // This structure is used to store the bounds of the kd tree partition
      // the bounds of the tree partition are computed as part of the
      // recursive search procedure.
      struct bounds {
        double minc[3] ;
        double maxc[3] ;
        bounds() {
          minc[0] = -std::numeric_limits<float>::max() ;
          minc[1] = minc[0] ;
          minc[2] = minc[0] ;
          maxc[0] =  std::numeric_limits<float>::max() ;
          maxc[1] = maxc[0] ;
          maxc[2] = maxc[0] ;
        }
      } ;

    private:
      // These two arrays form an implicit kdtree where split gives the dividing
      // line between the two halves of the tree for the current pivot
      std::vector<coord_info> pnts ;
      std::vector<int> splits ;

      // Select a pivot and move to begining of array.
      void select_pivot(int start, int end, int coord) ;

      // This will test if a kd tree region contains a given sphere
      // defined by the center v, and the radius squared r.  proj_d
      // is the projected distance (squared) from the sphere center to
      // the current dividing plane and coord is the current dividing plane.
      bool insphere(const coord3d &v,double r, double proj_d,
                    bounds bnd, int coord) ;

      // Recursively build the kd tree
      // start points to the pivot element, end is one past the last element
      // in the array.  Depth is the current depth, which is zero at the root
      void build_kd(int start, int end, int depth) ;

      // Recursively search for the depth of the kd tree
      int depth_search(int start,int end,int depth) ;

      // recursively search for the closest point, refining minimum distance
      // and bounds as we go
      int find_closest(int start,int end, int depth, const coord3d &v,
                       double &rmin,bounds bnds) ;
  

    public:
      kd_tree(const std::vector<coord3d> &inpnts,
              const std::vector<int> &ids) ;
      kd_tree(std::vector<coord_info> &inpnts) ;

      // User routine to get depth of the kd_tree
      int depth() {
        return depth_search(0,pnts.size(),0) ;
      }
      // Search for the closest point using the kd_tree
      int find_closest(coord3d v, double &rmin) {
        const int sp = find_closest(0,pnts.size(),0,v,rmin,bounds()) ;
        if(sp < 0)
          return std::numeric_limits<int>::min() ;
        // Look up id of matched point
        return pnts[sp].id ;
      }
      int find_closest(coord3d v) {
        // Start with infinite closest point distance
        double rmin = std::numeric_limits<double>::max() ;
        return find_closest(v,rmin) ;
      }
    } ;
  }
}

#endif

#ifndef CUTPLANE_H
#define CUTPLANE_H

#include "grid.h"
#include "hdf5.h"

class CutPlane
{
public:
   void cut(cutplane_info &info, const LoadInfo& load_info,const positions3d& center, grid *figure);

private:
  bool registerFace(size_t faceNode[], int nNodes, int cellNum);
  void disambiguateFace(size_t faceNode[], int nNodes, int cellNum);
  void checkLoop(size_t start, size_t end);
  int sizeElementType(hid_t group_id, const char *element_name);
  void close();
  void insertEdges(edges &ed);
  void write_tets(int tets[], int ntets);
  void write_pyrm(int pyrm[], int npyrm);
  void write_prsm(int prsm[], int nprsm);
  void write_hexs(int hexs[], int nhexs);
  void write_general_cell(int nfaces[], int nnfaces,
			  int nsides[], 
			  int nodes[]);

  grid* fig;  // Points to the grid object to load with cutting plane
  vector<positions3d> nodes;  // Stores all 3d points in the grid
  vector<float> nodeVal;   // Stores scalar property values for each node
  vector<Array<size_t,5> > intersects ;  // List of all edges formed from 
                                         // plane intersection
  map<edges, size_t> nodeMap ;  // Maps old 3d node numbers to the new smaller 
                                // set of 2d nodes
  list<vector<size_t> > disambiguatedFaces ;  // Stores all unmatched 
            //previously disambiguated faces and corresponding fabricated nodes
  int cellCount;    // Tracks number of cells cut so far
};

#endif

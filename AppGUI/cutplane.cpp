/////////////////////////////////////////////////////////////////////////////
//  Filename: cutplane.cpp
//
//  Contains: The algorithm used to produce the 2d grid made from cutting a
//            3d grid with a plane.
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <cmath>
using std::signbit ;
#include <QString>
#include <QtOpenGL>
#include <algorithm>
#include <map>
#include <vector>
#include <list>
using std::map ;
using std::vector ;
using std::list ;

#include <QtDebug>

#include "hdf5.h"
#include "cutplane.h"
#include "grid.h"

// **** Call back functions for GLUtesselator ****

/////////////////////////////////////////////////////////////////////////////
//  Global Function
//
//  Callback function used by GLUtesselator. Organizes and places triangles
//  into the appropriate data structure.
/////////////////////////////////////////////////////////////////////////////

void cbVertex(void *vertex_data, void *user_data) {
  static void* data[3];
  static  int count = 0;

  data[count] = vertex_data;
  count = ++count % 3;

  if (count == 0) {
    grid* fig = (grid*) user_data;
    fig->triangle_list.push_back(triangles( *(int*)data[0], 
                                            *(int*)data[1], 
                                            *(int*)data[2] ));
  }
}

// void cbVertex(void *vertex_data, void *user_data) {
//   vector<int>* tri = (vector<int>*) user_data;
//   tri->push_back(*(int*)vertex_data);
// }


///////////////////////////////////////////////////////////////////
//  Global Function
//
//  Dummy function to force GLUtesselator into GL_TRIANLGES mode.
///////////////////////////////////////////////////////////////////

void cbEdgeFlag(GLboolean) {}

// **** CutPlane implementation ****

//////////////////////////////////////////////////////////////////////////////
//  public:
//    void cut(cutplane_info &info, grid *figure);
//
//  Reads in the grid and calls all functions necessary to generate a cutting
//  plane.
//////////////////////////////////////////////////////////////////////////////

void CutPlane::cut(cutplane_info &info, const LoadInfo &load_info,const positions3d& center, grid *figure)
{
 
  // Turn off hdf5 error reports
  H5Eset_auto(NULL, NULL);

  fig = figure;
  cellCount = 0;

  // Get list of nodes
  hsize_t npnts;

 
  QString posname = load_info.directory + "/output/grid_pos." + load_info.iteration + 
    '_' + load_info.casename ;
 
  hid_t file_id = H5Fopen(posname.toLocal8Bit(),
			  H5F_ACC_RDONLY, H5P_DEFAULT) ;

  hid_t dataset_id = H5Dopen(file_id, "/pos/data");
  hid_t dataspace_id = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(dataspace_id, &npnts, NULL);

  hid_t pos_tid = H5Tcreate(H5T_COMPOUND, sizeof(positions3d));
  
  // H5Tinsert(pos_tid, "x", HOFFSET(positions3d, x), H5T_IEEE_F64LE);
//   H5Tinsert(pos_tid, "y", HOFFSET(positions3d, y), H5T_IEEE_F64LE);
//   H5Tinsert(pos_tid, "z", HOFFSET(positions3d, z), H5T_IEEE_F64LE);

  H5Tinsert(pos_tid, "x", 0, H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "y", sizeof(double), H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "z", 2*sizeof(double), H5T_IEEE_F64LE);

  positions3d null3d;
  null3d.x = null3d.y = null3d.z = 0.0;
  nodes.assign(npnts, null3d);
  H5Dread(dataset_id, pos_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodes[0]);

  H5Tclose(pos_tid);
  H5Dclose(dataset_id);
  H5Fclose(file_id);

  // Transform nodes into cutting position
  affineMapping transMatrix;
  positions3d negTranslate = positions3d(-info.translate.x,
                                         -info.translate.y,
                                         -info.translate.z); 
  
  transMatrix.translate(center);
  transMatrix.rotateX(-info.rotate.x);
  transMatrix.rotateY(-info.rotate.y);
  transMatrix.rotateZ(-info.rotate.z);
  transMatrix.translate(negTranslate);  

  
  for(size_t i = 0; i < npnts; i++)
    nodes[i] = transMatrix.MapNode(nodes[i]);

  // Get connectivity infomation
  QString gridtopo = load_info.directory + "/output/" + load_info.casename +".topo" ;
  hid_t topo_id = H5Fopen(gridtopo.toLocal8Bit(),
			  H5F_ACC_RDONLY, H5P_DEFAULT) ;

  hid_t elg = H5Gopen(topo_id,"elements") ;

  int ntets = sizeElementType(elg,"tetrahedra") ;
  int nhexs = sizeElementType(elg,"hexahedra") ;
  int nprsm = sizeElementType(elg,"prism") ;
  int npyrm = sizeElementType(elg,"pyramid") ;
  int ngenc = sizeElementType(elg,"GeneralCellNfaces") ;

 

  // Get variable information
  QString filename = load_info.directory + "/output/" + load_info.variable + "_sca." + 
                     load_info.iteration + "_" + load_info.casename;
  hid_t scalar_id = H5Fopen(filename.toLocal8Bit(), 
			    H5F_ACC_RDONLY, H5P_DEFAULT);

  nodeVal.assign(npnts, 0.0);
  QString datasetName = "/" + load_info.variable + "/data";
  dataset_id = H5Dopen(scalar_id, datasetName.toLocal8Bit());
  H5Dread(dataset_id, H5T_IEEE_F32LE, 
	  H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodeVal[0]);
  H5Fclose(scalar_id);

  // Perform the cutting operations for ...
  if(ntets > 0) {  // Tetrahedra
    int *tets = new int[ntets * 4];

    hsize_t four[1] = {4};
    hid_t tets_tid = H5Tarray_create(H5T_NATIVE_INT, 1, four, NULL);

    dataset_id = H5Dopen(elg, "tetrahedra");
    H5Dread(dataset_id, tets_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, tets);

    H5Dclose(dataset_id);
    H5Tclose(tets_tid);
    write_tets(tets, ntets);
    delete [] tets;
  }
  if(npyrm > 0) {  // Pyramids
    int *pyrm = new int[npyrm * 5];

    hsize_t five[1] = {5};
    hid_t pyrm_tid = H5Tarray_create(H5T_NATIVE_INT, 1, five, NULL);

    dataset_id = H5Dopen(elg, "pyramid");
    H5Dread(dataset_id, pyrm_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, pyrm);

    H5Dclose(dataset_id);
    H5Tclose(pyrm_tid);
    write_pyrm(pyrm, npyrm);
    delete [] pyrm;
  }
  if(nprsm > 0) {  // Prisms
    int *prsm = new int[nprsm * 6];

    hsize_t six[1] = {6};
    hid_t prsm_tid = H5Tarray_create(H5T_NATIVE_INT, 1, six, NULL);

    dataset_id = H5Dopen(elg, "prism");
    H5Dread(dataset_id, prsm_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, prsm);

    H5Dclose(dataset_id);
    H5Tclose(prsm_tid);
    write_prsm(prsm, nprsm);
    delete [] prsm;
  }
  if(nhexs > 0) {  // Hexahera
    int *hexs = new int[nhexs * 8];

    hsize_t eight[1] = {8};
    hid_t hexs_tid = H5Tarray_create(H5T_NATIVE_INT, 1, eight, NULL);

    dataset_id = H5Dopen(elg, "hexahedra");
    H5Dread(dataset_id, hexs_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, hexs);

    H5Dclose(dataset_id);
    H5Tclose(hexs_tid);
    write_hexs(hexs, nhexs);
    delete [] hexs;
  }
  if(ngenc > 0) {  // General cells
    int *GeneralCellNfaces = new int[ngenc];
    dataset_id = H5Dopen(elg, "GeneralCellNfaces");
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	    GeneralCellNfaces);
    H5Dclose(dataset_id);

    int nside = sizeElementType(elg, "GeneralCellNsides");
    int *GeneralCellNsides = new int[nside];
    dataset_id = H5Dopen(elg, "GeneralCellNsides");
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	    GeneralCellNsides);
    H5Dclose(dataset_id);

    int nnodes = sizeElementType(elg, "GeneralCellNodes");
    int *GeneralCellNodes = new int[nnodes];
    dataset_id = H5Dopen(elg, "GeneralCellNodes");
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	    GeneralCellNodes);
    H5Dclose(dataset_id);

    write_general_cell(GeneralCellNfaces, ngenc,
                       GeneralCellNsides, 
                       GeneralCellNodes);
    delete [] GeneralCellNfaces;
    delete [] GeneralCellNsides;
    delete [] GeneralCellNodes;
  }

  H5Gclose(elg);
  H5Fclose(topo_id);
 
  // Perform closing operations
  close();
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    int sizeElementType(hid_t group_id, const char *element_name);
//
//  This function uses the HDF5 interface to find the number of given elements
//  in the given group.
//////////////////////////////////////////////////////////////////////////////

int CutPlane::sizeElementType(hid_t group_id, const char *element_name) {
  // Open dataset
  hid_t dataset = H5Dopen(group_id,element_name) ;
  if(dataset < 0) {
    H5Eclear() ;
    return 0 ;
  }

  // Get the dataspace
  hid_t dspace = H5Dget_space(dataset) ;

  // Retrieve the size of the dataspace
  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;

  H5Dclose(dataset);

  return int(size);
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void close();
//
//  This function takes all the cut edge information and compiles it into the
//  finished cutting plane cross section.
//////////////////////////////////////////////////////////////////////////////

void CutPlane::close() {
  if(intersects.size()==0) return;

  edges tempArr2 ;

  // Order each entry in intersects least to greatest
  vector<Array<size_t, 5> >::iterator it;
  for (it = intersects.begin() ; it != intersects.end() ; ++it) {
    if ( (*it)[0] > (*it)[2] || 
	 ((*it)[0] == (*it)[2] && (*it)[1] > (*it)[3]) ) {
      std::swap((*it)[0], (*it)[2]);
      std::swap((*it)[1], (*it)[3]);
    }
  }

  // Build a mapping of and interpolate all points for the new 2d grid
  for ( it = intersects.begin(); it != intersects.end(); ++it) {
    tempArr2.l = (*it)[0];
    tempArr2.r = (*it)[1];
    insertEdges(tempArr2);
    tempArr2.l = (*it)[2];
    tempArr2.r = (*it)[3];
    insertEdges(tempArr2);
  }

  // Make a list of indeces (used for the GLUtesselator callbacks)
  int *pntIndex;
  pntIndex = new int[fig->pos.size()];
  for (unsigned int i = 0; i < fig->pos.size(); ++i) {
    pntIndex[i] = i;
  }

  // Bubble sort intersects (takes very few passes)
  bool sorted = false;
  int pass = 0;
  while(!sorted) {
    pass++;
    sorted = true;
    vector<Array<size_t, 5> >::iterator scan, next;
    scan = intersects.begin();
    next = scan + 1;
    while(next != intersects.end()) {
      if ((*next)[4] < (*scan)[4]) {
	std::swap(*next, *scan);
        sorted = false;
      }

      scan++;
      next++;
    }
  }

  // Triangulate every cell formed by cutting plane
  unsigned int start = 0, end = 0;
  unsigned int ntris= 0, nquads =0;
  while (start < intersects.size()) {
    size_t currCell = intersects[start][4];

    // Find range of entries forming a single cell
    while (end <= intersects.size() && intersects[end][4] == currCell)
      end++;
    end--;
    
    // Make list of indeces for all edges forming the cell border
    list<int> edgeList;
    for (unsigned int i = start+1; i <= end; ++i)
      edgeList.push_back(i);
    
    edges firstNode(intersects[start][0], intersects[start][1]);
    edges scanNode(intersects[start][2], intersects[start][3]);

    vector<int> cellNodes;  // Keeps record of the cell's border traversal
    cellNodes.reserve((edgeList.size()+1) * 2);
    cellNodes.push_back(nodeMap[firstNode]);
    cellNodes.push_back(nodeMap[scanNode]);

    // Travel to every edge
    while(!edgeList.empty()) {

      // Find the next edge in the cell's border
      list<int>::iterator scan;
      for (scan = edgeList.begin(); scan != edgeList.end(); ++scan) {

	// When next edge is found, record it and remove it from the list
	if (intersects[*scan][0] == scanNode.l && 
	    intersects[*scan][1] == scanNode.r) {

	  scanNode.l = intersects[*scan][2];
	  scanNode.r = intersects[*scan][3];
	  
	  cellNodes.push_back(nodeMap[scanNode]);
	  edgeList.erase(scan);
	  break;

	} else if (intersects[*scan][2] == scanNode.l && 
		   intersects[*scan][3] == scanNode.r) {

	  scanNode.l = intersects[*scan][0];
	  scanNode.r = intersects[*scan][1];
	  
	  cellNodes.push_back(nodeMap[scanNode]);
	  edgeList.erase(scan);
	  break;
	}
      }
    }

    // If the first vertex is listed twice, remove the duplicate
    if (firstNode.l == scanNode.l && firstNode.r == scanNode.r)
      cellNodes.pop_back();

    // If the cell is not a triangle, tesselate it into triangles
    if (cellNodes.size() > 3) {
      nquads++;
     
      // Make GLUtesselator object
      GLUtesselator* myTess = gluNewTess();

      // Use the custom callback functions
      gluTessCallback(myTess, GLU_TESS_VERTEX_DATA,
		      (GLvoid (*) ()) &cbVertex);
      gluTessCallback(myTess, GLU_TESS_EDGE_FLAG,
		      (GLvoid (*) ()) &cbEdgeFlag);

      // Send the cell border to the tesselator
      gluTessBeginPolygon(myTess, fig);
      gluTessBeginContour(myTess);
      for (unsigned int i = 0; i < cellNodes.size(); ++i) {
	double point[3];
	point[0] = fig->pos[cellNodes[i]].x;
	point[1] = fig->pos[cellNodes[i]].y;
	point[2] = 0.0;
	gluTessVertex(myTess, point, &pntIndex[cellNodes[i]]);
      }
      gluTessEndContour(myTess);
      gluTessEndPolygon(myTess);
     
    }

    
    // If the cell was just already a triangle, push it to the grid's stack
    else if (cellNodes.size() == 3){
      ntris++;
      fig->triangle_list.push_back(triangles(cellNodes[0], 
					     cellNodes[1], 
					     cellNodes[2]));
    }else{
      qDebug()<<"not a polygon";
    }
    
    end++;
    start = end;
    
  } // End triangulation

  sort(intersects.begin(), intersects.end());

  // Build a unique list of edges
  vector<edges> exterior_edges;
  edges new_edge;
  int lastUnique = 0;
  fig->interior = 0;
  fig->edge_list.clear();
  for(size_t scan = 0; scan < intersects.size(); ++scan) {

    // If this edge is the boundary for two cells,
    if(intersects[scan][0] == intersects[scan + 1][0] && 
       intersects[scan][1] == intersects[scan + 1][1] && 
       intersects[scan][2] == intersects[scan + 1][2] && 
       intersects[scan][3] == intersects[scan + 1][3]) {

      // Record it as an inside edge.
      tempArr2.l = intersects[scan][0];
      tempArr2.r = intersects[scan][1];
      new_edge.l = nodeMap[tempArr2];
      tempArr2.l = intersects[scan][2];
      tempArr2.r = intersects[scan][3];
      new_edge.r = nodeMap[tempArr2];
      fig->edge_list.push_back(new_edge);

      scan++;
    } else {  // If this edge is the boundary for only one cell,
      // Record it as a boundary edge.
      tempArr2.l = intersects[scan][0];
      tempArr2.r = intersects[scan][1];
      new_edge.l = nodeMap[tempArr2];
      tempArr2.l = intersects[scan][2];
      tempArr2.r = intersects[scan][3];
      new_edge.r = nodeMap[tempArr2];
      exterior_edges.push_back(new_edge);
    }
    lastUnique++;
  }

  // Record number of internal edges
  fig->interior = fig->edge_list.size();

  // Append boundary edges to edge list
  fig->edge_list.insert(fig->edge_list.end(), 
			exterior_edges.begin(), exterior_edges.end());


  
  // Find the maximum and minimum x and y values in the grid
  positions MinPos = fig->pos[0];
  positions MaxPos = MinPos;

  for(size_t i=0;i<exterior_edges.size();++i) {
    MaxPos.x = qMax(MaxPos.x,fig->pos[exterior_edges[i].l].x) ;
    MaxPos.x = qMax(MaxPos.x,fig->pos[exterior_edges[i].r].x) ;
    MinPos.x = qMin(MinPos.x,fig->pos[exterior_edges[i].l].x) ;
    MinPos.x = qMin(MinPos.x,fig->pos[exterior_edges[i].r].x) ;
    MaxPos.y = qMax(MaxPos.y,fig->pos[exterior_edges[i].l].y) ;
    MaxPos.y = qMax(MaxPos.y,fig->pos[exterior_edges[i].r].y) ;
    MinPos.y = qMin(MinPos.y,fig->pos[exterior_edges[i].l].y) ;
    MinPos.y = qMin(MinPos.y,fig->pos[exterior_edges[i].r].y) ; 
  }
  fig->minview = fig->minpos = MinPos ;
  fig->maxview = fig->maxpos = MaxPos ;
  fig->size = qMax(MaxPos.x - MinPos.x, MaxPos.y - MinPos.y);


 
  // Find max/min of the scalar values as well
   fig->set_value_range();

  // Construct the contour curves
 fig->generate_contour_curves();

}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    bool registerFace(size_t faceNode[], int nNodes, int cellNum);
//
//  All the faces of cut grid cells get sent here. This function finds which
//  of the edges, if any, are cut and pushes them onto intersects.
//////////////////////////////////////////////////////////////////////////////

bool CutPlane::registerFace(size_t faceNode[], int nNodes, int cellNum) {
  Array<size_t, 5> faceCut;
  size_t thisNode, prevNode, cutsFound = 0;
  prevNode = faceNode[nNodes - 1];

  // Iterate over each edge of the face
  faceCut[4] = cellNum;
  for(int i = 0; i < nNodes; ++i) {
    thisNode = faceNode[i];

    // If the edge is cut, add it to intersects
    if (signbit(nodes[thisNode - 1].z * nodes[prevNode - 1].z)) {
      if (cutsFound < 2) {
	if (thisNode < prevNode) {
	  faceCut[0 + cutsFound*2] = thisNode ;
	  faceCut[1 + cutsFound*2] = prevNode ;
	} else {
	  faceCut[0 + cutsFound*2] = prevNode ;
	  faceCut[1 + cutsFound*2] = thisNode ;
	}
      }
      cutsFound++;
    }
    prevNode = thisNode;
  }

  if (cutsFound == 2) {
    intersects.push_back(faceCut);
  } else if (cutsFound > 2)
    // If more than two cut edges are found, this face needs disambiguating
    disambiguateFace(faceNode, nNodes, cellNum);

  // Report if the face was indeed cut (used for general cells)
  return (cutsFound > 0);
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void insertEdges(edges &ed);
//
//  Any edges sent to this function are added to the list of unique edges and
//  interpolated for the 2d grid. Duplicate edges are ignored.
//////////////////////////////////////////////////////////////////////////////

void CutPlane::insertEdges(edges &ed)
{
  // If this edge has been recorded already, ignore it
  if (nodeMap.find(ed) != nodeMap.end()) return;

  float t, x;
  size_t a, b;
  positions newPos;
  
  // Add edge to the mapping data structure
  nodeMap.insert(std::pair<edges,size_t>(ed, nodeMap.size()) );
  
  // Interpolate x and y values of new point and push to stack
  a = ed.l - 1;
  b = ed.r - 1;
  t = (nodes[b].z)/(nodes[b].z - nodes[a].z);
  
  newPos.x = t*nodes[a].x + (1.0 - t)*nodes[b].x;
  newPos.y = t*nodes[a].y + (1.0 - t)*nodes[b].y;
  fig->pos.push_back(newPos);

 
  x = t*nodeVal[a] + (1.0 - t)*nodeVal[b];
  fig->val.push_back(x);
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void disambiguateFace(size_t faceNode[], int nNodes, int cellNum);
//
//  If a face needs disambiguating (more than two edges are cut) this function
//  makes a new point in the center of the face and connects it with all the 
//  other points defining the face, forming triangles. These new triangle faces
//  are then re-sent to the registerFace() function.
//
//  It also checks if the same face has already been disambiguated for the
//  adjacent cell. In this case the point fabricated for the previous face
//  is used to prevent a duplicate being made.
//////////////////////////////////////////////////////////////////////////////

void CutPlane::disambiguateFace(size_t faceNode[], int nNodes, int cellNum) {
  bool needNewNode = true;
  size_t thisNode, prevNode, trigFace[3], centerNode;
  vector<size_t> tmpfaceNode(faceNode, faceNode + nNodes);
  list<vector<size_t> >::iterator face, remember;

  // Search for the given face in the list of previously disambiguated faces
  if (disambiguatedFaces.size() != 0) {
    sort(tmpfaceNode.begin(), tmpfaceNode.end());

    int i;
    for (face = disambiguatedFaces.begin(); 
	 face != disambiguatedFaces.end(); face++) {
      if ((*face).size() - 1 == static_cast<size_t>(nNodes)) {
	for (i = 0; i < nNodes; i++) {
	  if (tmpfaceNode[i] != (*face)[i+1])
	    break;
	}

	// If this face has already been disambiguated, reuse the previously
        //   fabricated point
	if (i == nNodes) {
	  centerNode = (*face)[0];
	  disambiguatedFaces.erase(face);
	  needNewNode = false;
	  break;
	}
      }
    }
  }
  
  if (needNewNode) {  // If the face has not been disambiguated
    positions3d newNode;
    newNode.x = newNode.y = newNode.z = 0;
    float tmp;

    // Use the average of the face's points to use as a new center point
    for (int i = 0; i < nNodes; ++i) {
      newNode.x += nodes[faceNode[i]-1].x;
      newNode.y += nodes[faceNode[i]-1].y;
      newNode.z += nodes[faceNode[i]-1].z;
    }
    newNode.x /= nNodes;
    newNode.y /= nNodes;
    newNode.z /= nNodes;
    nodes.push_back(newNode);  // Add the new point to the point list
    centerNode = nodes.size();

    // Also use the average scalar value as the new point's scalar value
    tmp = 0;
    for (int j = 0; j < nNodes; j++)
      tmp += nodeVal[faceNode[j]-1];
    tmp /= nNodes;
    nodeVal.push_back(tmp);

    // Add this face to the list of disambiguated faces
    vector<size_t> newFaceList;
    newFaceList.push_back(centerNode);
    for (int i = 0; i < nNodes; ++i) {
      newFaceList.push_back(faceNode[i]);
    }
    sort(newFaceList.begin()+1, newFaceList.end());
    disambiguatedFaces.push_back(newFaceList);
  }

  // Register the cuts on the disambiguated face
  trigFace[2] = centerNode;
  prevNode = faceNode[nNodes - 1];
  for (int i = 0; i < nNodes; ++i) {
    thisNode = faceNode[i];
    trigFace[0] = thisNode;
    trigFace[1] = prevNode;
    prevNode = thisNode;
    registerFace(trigFace, 3, cellNum);
  }
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void checkLoop(size_t start, size_t end);
//
//  This function traverses the edges formed by a cut cell to ensure that they
//  do indeed form a full loop. In the case that multiple loops exist, the
//  cell number for the following loops are incremented to ensure cell
//  number uniqueness.
//////////////////////////////////////////////////////////////////////////////

void CutPlane::checkLoop(size_t start, size_t end)
{
  size_t loopNumber = 0;

  // Form list of edge indeces
  list<size_t> edgeList;
  for (size_t i = start+1; i < end; i++)
    edgeList.push_back(i);

  edges firstNode, nextNode;
  firstNode.l = intersects[start][0];
  firstNode.r = intersects[start][1];
  nextNode.l = intersects[start][2];
  nextNode.r = intersects[start][3];

  // Traverse the entire list of edges
  list<size_t>::iterator iter;
  while (!edgeList.empty()) {
    for (iter = edgeList.begin(); iter != edgeList.end(); iter++) {
      if (intersects[*iter][0] == nextNode.l && 
	  intersects[*iter][1] == nextNode.r) {
	nextNode.l = intersects[*iter][2];
	nextNode.r = intersects[*iter][3];
	intersects[*iter][4] += loopNumber;
	edgeList.erase(iter);
	break;
      } else if (intersects[*iter][2] == nextNode.l && 
		 intersects[*iter][3] == nextNode.r) {
	nextNode.l = intersects[*iter][0];
	nextNode.r = intersects[*iter][1];
	intersects[*iter][4] += loopNumber;
	edgeList.erase(iter);
	break;
      }
    }

    // If a full loop has been made and there are more edges to traverse,
    if ((firstNode.l == nextNode.l && firstNode.r == nextNode.r) && 
	!edgeList.empty()) {
      // Increment the cell number and traverse the remaining edges.
      cellCount++;
      loopNumber++;
      firstNode.l = intersects[edgeList.front()][0];
      firstNode.r = intersects[edgeList.front()][1];
      nextNode.l = intersects[edgeList.front()][2];
      nextNode.r = intersects[edgeList.front()][3];

      intersects[edgeList.front()][4] += loopNumber;
      edgeList.erase(edgeList.begin());
    }
  }

  // If there is an incomplete loop, report the problem
  if (firstNode.l != nextNode.l || firstNode.r != nextNode.r)
    qDebug("** Problem cell (failed loop test) **");
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void write_tets(int tets[], int ntets);
//
//  This function checks if each tetrahedra is cut. Each cut tetrahedra then
//  has each of its faces sent to the registerFace() function and checkLoop()
//  is called with the cell's appropriate range in intersects[][].
//////////////////////////////////////////////////////////////////////////////

void CutPlane::write_tets(int tets[], int ntets) {
  bool isCut;
  int tetsCut = 0;
  size_t faceNodeIndex[3];
  size_t faceNode[3], intersectStart = intersects.size();

  for (int i = 0 ; i < ntets ; ++i) {
    isCut = false;
    for (int j = 1; j < 4; ++j) {
      if (signbit(nodes[tets[i*4 + 0] - 1].z * nodes[tets[i*4 + j] - 1].z)) {
	isCut = true;
	break;
      }
    }

    if(isCut) {
      tetsCut++ ;
      for (faceNodeIndex[2] = 2 ; faceNodeIndex[2] < 4 ; ++faceNodeIndex[2]) {
 	for (faceNodeIndex[1] = 1 ; faceNodeIndex[1] < faceNodeIndex[2] ; ++faceNodeIndex[1]) {
	  for (faceNodeIndex[0] = 0 ; faceNodeIndex[0] < faceNodeIndex[1] ; ++faceNodeIndex[0]) {
	    
	    for (int j = 0; j < 3; j++)
	      faceNode[j] = tets[i*4 + faceNodeIndex[j]];
	    registerFace(faceNode, 3, i+cellCount);
	    
	  }
	}
      }
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }

  cellCount += ntets;
 
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void write_pyrm(int pyrm[], int npyrm);
//
//  This function checks if each pyramid is cut. Each cut pyramid then has each
//  of its faces sent to the registerFace() function and checkLoop() is called
//  with the cell's appropriate range in intersects[][].
//////////////////////////////////////////////////////////////////////////////

void CutPlane::write_pyrm(int pyrm[], int npyrm) {
  bool isCut;
  int pyrmCut = 0, prevNode;
  size_t faceNode[4], intersectStart = intersects.size();

  for (int i = 0 ; i < npyrm ; ++i) {
    isCut = false;
    for (int j = 1; j < 5; ++j) {
      if (signbit(nodes[pyrm[i*5 + 0] - 1].z * nodes[pyrm[i*5 + j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 3;
      for (int thisNode = 0; thisNode < 4; ++thisNode) {
	faceNode[0] = pyrm[i*5 + prevNode];
	faceNode[1] = pyrm[i*5 + thisNode];
	faceNode[2] = pyrm[i*5 + 4];
	registerFace(faceNode, 3, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 4; ++j)
	faceNode[j] = pyrm[i*5 + j];
      registerFace(faceNode, 4, i+cellCount);
      pyrmCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }
  
  cellCount += npyrm;
  
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void write_prsm(int prsm[], int nprsm);
//
//  This function checks if each prism is cut. Each cut prism then has each of
//  its faces sent to the registerFace() function and checkLoop() is called
//  with the cell's appropriate range in intersects[][].
//////////////////////////////////////////////////////////////////////////////

void CutPlane::write_prsm(int prsm[], int nprsm) {
  bool isCut;
  int prsmCut = 0, prevNode;
  size_t faceNode[4], intersectStart = intersects.size();

  for (int i = 0 ; i < nprsm ; ++i) {
    isCut = false;
    for (int j = 1; j < 6; ++j) {
      if (signbit(nodes[prsm[i*6 + 0] - 1].z * nodes[prsm[i*6 + j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 2;
      for (int thisNode = 0; thisNode < 3; ++thisNode) {
	faceNode[0] = prsm[i*6 + prevNode];
	faceNode[1] = prsm[i*6 + thisNode];
	faceNode[2] = prsm[i*6 + thisNode + 3];
	faceNode[3] = prsm[i*6 + prevNode + 3];
	registerFace(faceNode, 4, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < 3; ++k) {
	  faceNode[k] = prsm[i*6 + (k + j*3)];
	}
	registerFace(faceNode, 3, i+cellCount);
      }
      prsmCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }

  }

  cellCount += nprsm;
  
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void write_hexs(int hexs[], int nhexs);
//
//  This function checks if each hexahedra is cut. Each cut hexahedra then
//  has each of its faces sent to the registerFace() function and checkLoop()
//  is called with the cell's appropriate range in intersects[][].
//////////////////////////////////////////////////////////////////////////////

void CutPlane::write_hexs(int hexs[], int nhexs) {
  bool isCut;
  int hexsCut = 0, prevNode;
  size_t faceNode[4], intersectStart = intersects.size();

  for (int i = 0 ; i < nhexs ; ++i) {
    isCut = false;
    for (int j = 1; j < 8; ++j) {
      if (signbit(nodes[hexs[i*8 + 0] - 1].z * nodes[hexs[i*8 + j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 3;
      for (int thisNode = 0; thisNode < 4; ++thisNode) {
	faceNode[0] = hexs[i*8 + prevNode];
	faceNode[1] = hexs[i*8 + thisNode];
	faceNode[2] = hexs[i*8 + thisNode + 4];
	faceNode[3] = hexs[i*8 + prevNode + 4];
	registerFace(faceNode, 4, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < 4; ++k) {
	  faceNode[k] = hexs[i*8 + (k + j*4)];
	}
	registerFace(faceNode, 4, i+cellCount);
      }
      hexsCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }

  cellCount += nhexs;
 
}


//////////////////////////////////////////////////////////////////////////////
//  private:
//    void write_general_cell(int nfaces[], int nnfaces,
//                            int nsides[], int nnsides,
//                            int nodes[], int nnodes);
//
//  This function checks if each general cell is cut. Each cut general cell
//  then has each of its faces sent to the registerFace() function and
//  checkLoop() is called with the cell's appropriate range in intersects[][].
//////////////////////////////////////////////////////////////////////////////

void CutPlane::write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[],
			     int nodes[]) {
  int genCut = 0, celloffset = 0, faceoffset = 0;
  bool isCut;
  size_t *faceNode, intersectStart = intersects.size();
  
  for (int cell = 0; cell < nnfaces; ++cell) {
    isCut = false;
    for (int face = 0; face < nfaces[cell]; ++face) {
      faceNode = new size_t[nsides[face + faceoffset]];
      for (int node = 0; node < nsides[face + faceoffset]; ++node) {
	faceNode[node] = nodes[node + celloffset];
      }
      if (registerFace(faceNode, nsides[face + faceoffset], cell + cellCount))
	isCut = true;
      delete [] faceNode;
      celloffset += nsides[face + faceoffset];
    }
    faceoffset += nfaces[cell];
    if (isCut) {
      genCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }
  
 
}


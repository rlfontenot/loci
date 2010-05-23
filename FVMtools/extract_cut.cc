

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
#include <Loci.h> 
#include <Tools/tools.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;
using std::list ;
using std::sort ;

#include "extract.h"

cuttingplane_topo_handler::cuttingplane_topo_handler(affineMapping &transformMatrix,
						     float xShift, float yShift, float zShift) {
  vector3d<float> transVec;
  transVec.x = xShift;
  transVec.y = yShift;
  transVec.z = zShift;

  transMatrix.translate(transVec);

  transMatrix.Combine(transformMatrix);
}

void cuttingplane_topo_handler::open(string casename, string iteration ,int npnts,
				     int ntets, int nprsm, int npyrm, int nhexs, int ngen,
				     const vector<string> &bc_names,
				     const vector<string> &variables,
				     const vector<int> &variable_types,
                                     double time) {
  strIter = iteration;
  numDisFaces = 0;
}

template <int T>  class ArrayLessThan {
public:
  bool operator()(const Array<int,T> &arr1, const Array<int,T> &arr2) {
    for (int i = 0; i < T; ++i) {
      if (arr1[i] == arr2[i]) continue;
      if (arr1[i] < arr2[i]) return true;
      break;
    }
    return false;
  }
} ;




void cuttingplane_topo_handler::close() {
  int temp1, temp2 ;
  Array<int,2> tempArr2 ;
  for (unsigned int i = 0 ; i < intersects.size() ; ++i) {
    tempArr2[0] = intersects[i][0];
    tempArr2[1] = intersects[i][1];
    nodeMap.push_back(tempArr2);
    tempArr2[0] = intersects[i][2];
    tempArr2[1] = intersects[i][3];
    nodeMap.push_back(tempArr2);

    cellMap[intersects[i][4]] = cellMap.size();

    if (intersects[i][0] > intersects[i][2]) {
      temp1 = intersects[i][0] ;
      temp2 = intersects[i][1] ;
      intersects[i][0] = intersects[i][2] ;
      intersects[i][1] = intersects[i][3] ;
      intersects[i][2] = temp1 ;
      intersects[i][3] = temp2 ;
    }
    else if ((intersects[i][0] == intersects[i][2]) && (intersects[i][1] > intersects[i][3])) {
      temp1 = intersects[i][0] ;
      temp2 = intersects[i][1] ;
      intersects[i][0] = intersects[i][2] ;
      intersects[i][1] = intersects[i][3] ;
      intersects[i][2] = temp1 ;
      intersects[i][3] = temp2 ;
    }
  }

  sort(intersects.begin(), intersects.end(), ArrayLessThan<5>());
  sort(nodeMap.begin(), nodeMap.end(), ArrayLessThan<2>());

  int lastUnique = 0;
  for (unsigned int scan = 0; scan < nodeMap.size(); ++scan) {
    if (!(nodeMap[scan][0] == nodeMap[lastUnique][0] && nodeMap[scan][1] == nodeMap[lastUnique][1])) {
      lastUnique++;
      nodeMap[lastUnique] = nodeMap[scan];
    }
  }
  nodeMap.resize(lastUnique+1);

  lastUnique = 0;
  for(unsigned int scan = 0; scan < intersects.size(); ++scan) {
    if(intersects[scan][0] == intersects[scan + 1][0] && intersects[scan][1] == intersects[scan + 1][1] && 
       intersects[scan][2] == intersects[scan + 1][2] && intersects[scan][3] == intersects[scan + 1][3]) {
      tempArr2[0] = intersects[scan][0];
      tempArr2[1] = intersects[scan][1];
      intersects[lastUnique][0] = searchNodeMap(tempArr2);
      tempArr2[0] = intersects[scan][2];
      tempArr2[1] = intersects[scan][3];
      intersects[lastUnique][1] = searchNodeMap(tempArr2);
      intersects[lastUnique][2] = cellMap[intersects[scan][4]] - 1;
      intersects[lastUnique][3] = cellMap[intersects[scan+1][4]] - 1;
      intersects[lastUnique][4] = 0;
      scan++;
    } else {
      tempArr2[0] = intersects[scan][0];
      tempArr2[1] = intersects[scan][1];
      intersects[lastUnique][0] = searchNodeMap(tempArr2);
      tempArr2[0] = intersects[scan][2];
      tempArr2[1] = intersects[scan][3];
      intersects[lastUnique][1] = searchNodeMap(tempArr2);
      intersects[lastUnique][2] = cellMap[intersects[scan][4]] - 1;
      intersects[lastUnique][3] = -1;
      intersects[lastUnique][4] = 0;
    }
    lastUnique++;
  }
  intersects.resize(lastUnique);

  vector<float> coords[2];
  float t, x;
  int a, b;
  for (unsigned int i = 0; i < nodeMap.size(); i++) {
    a = nodeMap[i][0] - 1;
    b = nodeMap[i][1] - 1;
    t = (nodes[b].z)/(nodes[b].z - nodes[a].z);
    coords[0].push_back( t*nodes[a].x + (1.0 - t)*nodes[b].x );
    coords[1].push_back( t*nodes[a].y + (1.0 - t)*nodes[b].y );
  }

  ofstream outFile;
  for (unsigned int varPass = 0; varPass < variables.size(); varPass++) {
    string out = variables[varPass] + "." + strIter;
    outFile.open(out.c_str());
    outFile << "general" << endl;
    outFile << nodeMap.size() << " 0" << endl;
    for (unsigned int i = 0; i < nodeMap.size(); i++)
      outFile << coords[0][i] << ' ' << coords[1][i] << endl;
    
    outFile << intersects.size() << " 0 " << cellMap.size() << " 0" << endl;
    for (unsigned int i = 0; i < intersects.size(); i++)
      outFile << intersects[i][0] << ' ' << intersects[i][1] << ' ' << intersects[i][2] << ' ' << intersects[i][3] << endl;

    for (unsigned int i = 0; i < nodeMap.size(); i++) {
      a = nodeMap[i][0] - 1;
      b = nodeMap[i][1] - 1;
      t = (nodes[b].z)/(nodes[b].z - nodes[a].z);
      x = t*nodeVal[varPass][a] + (1.0 - t)*nodeVal[varPass][b];
      outFile << x << endl;
    }
    outFile.close();
  }
}

template <int T>
bool cuttingplane_topo_handler::arrLessThan(Array<int,T> arr1, Array<int,T> arr2) {
  for (int i = 0; i < T; ++i) {
    if (arr1[i] == arr2[i]) continue;
    if (arr1[i] < arr2[i]) return true;
    break;
  }
  return false;
}

int cuttingplane_topo_handler::searchNodeMap(Array<int, 2> arr1) {
  int middle, front = 0, back = nodeMap.size() - 1;
  while(front + 1 < back) {
    middle = (front + back)/2;
    if (arrLessThan<2>(arr1, nodeMap[middle]))
      back = middle;
    else
      front = middle;
  }
  if (arrLessThan<2>(arr1, nodeMap[back]))
    return front;
  else 
    return back;
}

bool cuttingplane_topo_handler::registerFace(int faceNode[], int nNodes, int cellNum) {
  Array<int, 5> faceCut;
  int thisNode, prevNode, cutsFound = 0;
  prevNode = faceNode[nNodes - 1];

  faceCut[4] = cellNum;
  for(int i = 0; i < nNodes; ++i) {
    thisNode = faceNode[i];

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
  if (cutsFound == 2)
    intersects.push_back(faceCut);
  else if (cutsFound > 2)
    disambiguateFace(faceNode, nNodes, cellNum);

  return (cutsFound > 0);
}

void cuttingplane_topo_handler::disambiguateFace(int faceNode[], int nNodes, int cellNum) {
  bool needNewNode = true;
  int thisNode, prevNode, trigFace[3], centerNode;
  vector<int> tmpfaceNode(faceNode, faceNode + nNodes);
  list<vector<int> >::iterator face, remember;
  numDisFaces++;

  if (disambiguatedFaces.size() != 0) {
    sort(tmpfaceNode.begin(), tmpfaceNode.end());

    int i;
    for (face = disambiguatedFaces.begin(); face != disambiguatedFaces.end(); face++) {
      if ((*face).size() - 1 == static_cast<unsigned int>(nNodes)) {
	for (i = 0; i < nNodes; i++) {
	  if (tmpfaceNode[i] != (*face)[i+1])
	    break;
	}
	if (i == nNodes) {
	  centerNode = (*face)[0];
	  disambiguatedFaces.erase(face);
	  needNewNode = false;
	  break;
	}
      }
    }
  }
  
  if (needNewNode) {
    vector3d<float> newNode;
    newNode.x = newNode.y = newNode.z = 0;
    float tmp;

    for (int i = 0; i < nNodes; ++i) {
      newNode.x += nodes[faceNode[i]-1].x;
      newNode.y += nodes[faceNode[i]-1].y;
      newNode.z += nodes[faceNode[i]-1].z;
    }
    newNode.x /= nNodes;
    newNode.y /= nNodes;
    newNode.z /= nNodes;
    nodes.push_back(newNode);
    centerNode = nodes.size();

    for (unsigned int i = 0; i < nodeVal.size(); i++) {
      tmp = 0;
      for (int j = 0; j < nNodes; j++) {
	tmp += nodeVal[i][faceNode[j]-1];
      }
      tmp /= nNodes;
      nodeVal[i].push_back(tmp);
    }

    vector<int> newFaceList;
    newFaceList.push_back(centerNode);
    for (int i = 0; i < nNodes; ++i) {
      newFaceList.push_back(faceNode[i]);
    }
    sort(newFaceList.begin()+1, newFaceList.end());
    disambiguatedFaces.push_back(newFaceList);
  }

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

void cuttingplane_topo_handler::checkLoop(int start, int end) {
  bool loopIsGood = true;
  int loopNumber = 0;

  list<int> edges;
  for (int i = start+1; i < end; i++)
    edges.push_back(i);

  Array<int, 2> firstNode, nextNode;
  firstNode[0] = intersects[start][0];
  firstNode[1] = intersects[start][1];
  nextNode[0] = intersects[start][2];
  nextNode[1] = intersects[start][3];

  list<int>::iterator iter;
  while (!edges.empty() && loopIsGood) {
    for (iter = edges.begin(); iter != edges.end(); iter++) {
      if (intersects[*iter][0] == nextNode[0] && intersects[*iter][1] == nextNode[1]) {
	nextNode[0] = intersects[*iter][2];
	nextNode[1] = intersects[*iter][3];
	intersects[*iter][4] += loopNumber;
	edges.erase(iter);
	break;
      } else if (intersects[*iter][2] == nextNode[0] && intersects[*iter][3] == nextNode[1]) {
	nextNode[0] = intersects[*iter][0];
	nextNode[1] = intersects[*iter][1];
	intersects[*iter][4] += loopNumber;
	edges.erase(iter);
	break;
      }
    }
    if (iter == edges.end())
      loopIsGood = false;
    if ((firstNode[0] == nextNode[0] && firstNode[1] == nextNode[1]) && !edges.empty()) {
      cellCount++;
      loopNumber++;
      firstNode[0] = intersects[edges.front()][0];
      firstNode[1] = intersects[edges.front()][1];
      nextNode[0] = intersects[edges.front()][2];
      nextNode[1] = intersects[edges.front()][3];
      intersects[edges.front()][4] += loopNumber;
      edges.erase(edges.begin());
    }
  }
  if (firstNode[0] != nextNode[0] || firstNode[1] != nextNode[1])
    cerr << "** Problem cell: " << intersects[start][4] << " ** (failed loop test)" << endl;
}

void cuttingplane_topo_handler::create_mesh_positions(vector3d<float> pos[], int npnts) {
  for(int i = 0 ; i < npnts ; ++i) {
    nodes.push_back(pos[i]);       // Make a local copy of all node positions
  }
  
  for(int i = 0; i < npnts; i++) {
    nodes[i] = transMatrix.MapNode(nodes[i]);   // Transform topology into cutting position
  }

  cellCount = 0;
}

void cuttingplane_topo_handler::create_mesh_elements() {
}

void cuttingplane_topo_handler::write_tets(Array<int,4> tets[], int ntets) {
  bool isCut;
  int tetsCut = 0, intersectStart = intersects.size();
  int faceNodeIndex[3];
  int faceNode[3];

  for (int i = 0 ; i < ntets ; ++i) {
    isCut = false;
    for (int j = 1; j < 4; ++j) {
      if (signbit(nodes[tets[i][0] - 1].z * nodes[tets[i][j] - 1].z)) {
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
	      faceNode[j] = tets[i][faceNodeIndex[j]];
	    registerFace(faceNode, 3, i);

	  }
	}
      }
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }

  cellCount += ntets;
  cout << "Number of tets cut: " << tetsCut << endl ;
}

void cuttingplane_topo_handler::write_pyrm(Array<int,5> pyrm[], int npyrm) {  
  bool isCut;
  int pyrmCut = 0, prevNode, intersectStart = intersects.size();
  int faceNode[4];

  for (int i = 0 ; i < npyrm ; ++i) {
    isCut = false;
    for (int j = 1; j < 5; ++j) {
      if (signbit(nodes[pyrm[i][0] - 1].z * nodes[pyrm[i][j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 3;
      for (int thisNode = 0; thisNode < 4; ++thisNode) {
	faceNode[0] = pyrm[i][prevNode];
	faceNode[1] = pyrm[i][thisNode];
	faceNode[2] = pyrm[i][4];
	registerFace(faceNode, 3, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 4; ++j)
	faceNode[j] = pyrm[i][j];
      registerFace(faceNode, 4, i+cellCount);
      pyrmCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }
  
  cellCount += npyrm;
  cout << "Number of pyrm cut: " << pyrmCut << endl ;
}

void cuttingplane_topo_handler::write_prsm(Array<int,6> prsm[], int nprsm) {
  bool isCut;
  int prsmCut = 0, prevNode, intersectStart = intersects.size();
  int faceNode[4];

  for (int i = 0 ; i < nprsm ; ++i) {
    isCut = false;
    for (int j = 1; j < 6; ++j) {
      if (signbit(nodes[prsm[i][0] - 1].z * nodes[prsm[i][j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 2;
      for (int thisNode = 0; thisNode < 3; ++thisNode) {
	faceNode[0] = prsm[i][prevNode];
	faceNode[1] = prsm[i][thisNode];
	faceNode[2] = prsm[i][thisNode + 3];
	faceNode[3] = prsm[i][prevNode + 3];
	registerFace(faceNode, 4, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < 3; ++k) {
	  faceNode[k] = prsm[i][k + j*3];
	}
	registerFace(faceNode, 3, i+cellCount);
      }
      prsmCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }

  }

  cellCount += nprsm;
  cout << "Number of prsm cut: " << prsmCut << endl ;
}

void cuttingplane_topo_handler::write_hexs(Array<int,8> hexs[], int nhexs) {
  bool isCut;
  int hexsCut = 0, prevNode, intersectStart = intersects.size();
  int faceNode[4];

  for (int i = 0 ; i < nhexs ; ++i) {
    isCut = false;
    for (int j = 1; j < 8; ++j) {
      if (signbit(nodes[hexs[i][0] - 1].z * nodes[hexs[i][j] - 1].z)) {
	isCut = true;
	break;
      }
    }
    if (isCut) {
      prevNode = 3;
      for (int thisNode = 0; thisNode < 4; ++thisNode) {
	faceNode[0] = hexs[i][prevNode];
	faceNode[1] = hexs[i][thisNode];
	faceNode[2] = hexs[i][thisNode + 4];
	faceNode[3] = hexs[i][prevNode + 4];
	registerFace(faceNode, 4, i+cellCount);
	prevNode = thisNode;
      }
      for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < 4; ++k) {
	  faceNode[k] = hexs[i][k + j*4];
	}
	registerFace(faceNode, 4, i+cellCount);
      }
      hexsCut++;
      checkLoop(intersectStart, intersects.size());
      intersectStart = intersects.size();
    }
  }

  cellCount += nhexs;
  cout << "Number of hexs cut: " << hexsCut << endl ;
}

void cuttingplane_topo_handler::write_general_cell(int nfaces[], int nnfaces,
                                                   int nsides[], int nnsides,
						   int nodes[], int nnodes) {
  int genCut = 0, celloffset = 0, faceoffset = 0;
  bool isCut;
  int *faceNode, intersectStart = intersects.size();

  for (int cell = 0; cell < nnfaces; ++cell) {
    isCut = false;
    for (int face = 0; face < nfaces[cell]; ++face) {
      faceNode = new int[nsides[face + faceoffset]];
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

  cout << "Number of genc cut: " << genCut << endl;
}

void cuttingplane_topo_handler::close_mesh_elements() {
}

void cuttingplane_topo_handler::create_boundary_part(string name,int node_set[],int npnts) {
}

void cuttingplane_topo_handler::write_quads(Array<int,4> quads[], int quad_ids[],
					    int nquads) {
}

void cuttingplane_topo_handler::write_trias(Array<int,3> trias[], int tria_ids[],
					    int ntrias) {
}

void cuttingplane_topo_handler::write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
						   int nside_nodes[], int nside_nodes_size) {
}

void cuttingplane_topo_handler::close_boundary_part() {
}

void cuttingplane_topo_handler::create_nodal_output() {
}

void cuttingplane_topo_handler::output_nodal_scalar(float val[], int npnts, string valname) {
  vector<float> newVal;

  for (int i = 0; i < npnts; ++i) {
    newVal.push_back(val[i]);
  }

  nodeVal.push_back(newVal);
  variables.push_back(valname);
}

void cuttingplane_topo_handler::output_nodal_vector(vector3d<float> val[],
						    int npnts, string valname) {
}

void cuttingplane_topo_handler::close_nodal_output() {
}

void cuttingplane_topo_handler::output_boundary_scalar(float val[], int node_set[],
						       int nvals, string valname) {
}

void cuttingplane_topo_handler::output_boundary_vector(vector3d<float> val[], int node_set[],
						       int nvals, string valname) {
}

affineMapping::affineMapping() {
  for(int i=0;i<4;++i) {
    for(int j=0;j<4;++j)
      M[i][j] = 0 ;
    M[i][i] = 1 ;
  }
}
void affineMapping::Combine(affineMapping a) {
  float Mtmp[4][4] ;
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j) {
      float mtmp = 0 ;
      for(int k=0;k<4;++k)
	mtmp += a.M[i][k]*M[k][j] ;
      Mtmp[i][j] = mtmp ;
    }
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j)
      M[i][j] = Mtmp[i][j] ;
  }
void affineMapping::translate(vector3d<float> tv) {
  affineMapping tmp ;
  tmp.M[0][3] = tv.x ;
  tmp.M[1][3] = tv.y ;
  tmp.M[2][3] = tv.z ;
  Combine(tmp) ;
}
void affineMapping::rotateX(float theta) {
  float th = theta*2.*M_PI/360. ;
  float sth = sin(th) ;
  float cth = cos(th) ;
  affineMapping tmp ;
  
  tmp.M[1][1] =  cth ;
  tmp.M[1][2] =  sth ;
  tmp.M[2][1] = -sth ;
  tmp.M[2][2] =  cth ;
  Combine(tmp) ;
}
void affineMapping::rotateY(float theta) {
  float th = theta*2.*M_PI/360. ;
  float sth = sin(th) ;
  float cth = cos(th) ;
  affineMapping tmp ;
  
  tmp.M[0][0] =  cth ;
  tmp.M[0][2] = -sth ;
  tmp.M[2][0] =  sth ;
  tmp.M[2][2] =  cth ;
  Combine(tmp) ;
}
void affineMapping::rotateZ(float theta) {
  float th = theta*2.*M_PI/360. ;
  float sth = sin(th) ;
  float cth = cos(th) ;
  affineMapping tmp ;
  
  tmp.M[0][0] =  cth ;
  tmp.M[0][1] =  sth ;
  tmp.M[1][0] = -sth ;
  tmp.M[1][1] =  cth ;
  Combine(tmp) ;
}
vector3d<float> affineMapping::MapNode(vector3d<float> v) {
  float tmp[4] ;
  tmp[0] = v.x ;
  tmp[1] = v.y ;
  tmp[2] = v.z ;
  tmp[3] = 1. ;
  float res[4] ;
  for(int i=0;i<4;++i)
    res[i] = 0 ;
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j)
      res[i] += M[i][j]*tmp[j] ;
  vector3d<float> r(res[0],res[1],res[2]) ;
  return r ;
}

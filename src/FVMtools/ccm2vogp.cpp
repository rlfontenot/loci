//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//# Copyright 2021, ATA Engineering
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
#include <stdio.h>
#include <strings.h>
#include <iostream>
#include <cstdio>
#include <stdlib.h>	
#include <vector>
#include <map>
#include <string>
#include "ADF.h"
#include "vogtools.h"
#include <Loci.h>
#include <sstream>
#include <vector>
#include <sstream>
#include <numeric>
using std::cerr;
using std::cout;
using std::endl;
using std::istringstream;
using std::pair;
using std::vector;

// ADF docs say that in C, 1 greater than the length should be allocated to make
// room for the null termination character
#define ADF_NAME_LENGTH_C (ADF_NAME_LENGTH + 1)
#define ADF_DATA_TYPE_LENGTH_C (ADF_DATA_TYPE_LENGTH + 1)

typedef int int32;
typedef long int int64;

int64 readNodeis(double parentID, const char *name, int32 **value);
int64 readNodeBlockis(double parentID, const char *name, int32 **value,
                      const int64 &blockStart, const int64 &blockEnd);
void getMeshID(double processorID, double *verticesID, double *topoID);

void checkError(int err, double nodeID, string str) {
  if (err != -1) {
    cerr << str << " (error " << err << ")";
    char nodeName[ADF_NAME_LENGTH_C] = {'\0'};
    ADF_Get_Name(nodeID, nodeName, &err);
    cerr << " node(" << nodeName << ")" << endl;
    exit(1);
  }
}

inline int type2size(string type) {
  if (type == "MT") return 0;
  if (type == "R4" || type == "I4") return 4;
  if (type == "C1") return 1;
  if (type == "R8" || type == "I8") return 8;
  cerr << "Bad type: " << type << endl;
  return 0;
}

void getDimensions(double node, int64 *nDims, int64 *dims,
                   bool ignoreExtended = false) {
  int i;
  int32 nDims32, dims32[ADF_MAX_DIMENSIONS];
  int err;
  double nodeEx;

  {
    if (!nDims) return;
    *nDims = 0;
    ADF_Get_Node_ID(node, "ExtendedSize", &nodeEx, &err);
    if (err == -1 && !ignoreExtended) { /* This node is an extendend node */
      /* Check that this node is really a 1D array like it should be */
      i = 0;
      ADF_Get_Number_of_Dimensions(nodeEx, &i, &err);
      if (err != -1 || i != 1) return;

      /* Now find out how big it is */
      ADF_Get_Dimension_Values(nodeEx, dims32, &err);
      if (err != -1) return;
      *nDims = int64(dims32[0]);

      /* Now read in the number of dimensions */
      if (dims) {
        ADF_Read_All_Data(nodeEx, (char *)dims, &err);
        if (err != -1) return;
      }
    } else { /* Just a regular node;  call the regular ADF functions */
      ADF_Get_Number_of_Dimensions(node, &nDims32, &err);
      if (err != -1) return;
      *nDims = int64(nDims32);

      if (dims) {
        ADF_Get_Dimension_Values(node, dims32, &err);
        if (err == 32 /* No data */ || err == 27 /* Dimension is zero */) {
          *nDims = int64(0);
          dims[0] = int64(0);
        } else if (err != -1)
          return;

        for (i = 0; i < *nDims; ++i) dims[i] = int64(dims32[i]);
      }
    }

    return;
  } /* end CHECK_ERROR scope */
}

string getType(double nodeID) {
  int err;
  char data_type[ADF_DATA_TYPE_LENGTH_C] = {'\0'};
  ADF_Get_Data_Type(nodeID, data_type, &err);
  checkError(err, nodeID, "Error getting type");
  data_type[2] = '\0';
  return string(data_type);
}

// if it is a regular node this will just be the dimensions of the node;
// if it is an extended node, these will be the dimensions of the extended
// data
int64 getSize(double nodeID, bool ignoreExtended = false) {
  int64 num_dims, dim_values[ADF_MAX_DIMENSIONS];
  string type = getType(nodeID);
  if (type2size(type) == 0) return 0;
  getDimensions(nodeID, &num_dims, dim_values, ignoreExtended);
  if (num_dims == 0) return 0;
  int64 size = 1;
  for (int i = 0; i < num_dims; i++) size *= dim_values[i];
  return size;
}

void readAllData(double const nodeID, char *buffer) {
  int err;
  ADF_Read_All_Data(nodeID, buffer, &err);
  checkError(err, nodeID, "Error reading data");

  double tmpID;
  ADF_Get_Node_ID(nodeID, "ExtendedSize", &tmpID, &err);
  // 29 - Specified child is not a child of the specified parent
  if (err == 29) {
    return;
  }

  // is extended data

  // Read the rest of the extended data.
  int64 size = getSize(nodeID, true);
  string type = getType(nodeID);
  char *local_buffer = buffer + type2size(type) * size;

  int num_children;
  ADF_Number_of_Children(nodeID, &num_children, &err);
  checkError(err, nodeID, "Error getting number of children");

  int inum_ret;
  char *names = new char[(ADF_NAME_LENGTH_C + 1) * num_children];
  ADF_Children_Names(nodeID, 1, num_children, ADF_NAME_LENGTH_C, &inum_ret,
                     names, &err);
  checkError(err, nodeID, "Error getting children names");

  char *cp;
  char *end = names + (ADF_NAME_LENGTH_C + 1) * inum_ret;
  for (cp = names; cp != end; cp += (ADF_NAME_LENGTH_C + 1)) {
    static char const *extendedData = "ExtendedData-";
    if (strncmp(cp, extendedData, strlen(extendedData))) continue;

    double ID;
    ADF_Get_Node_ID(nodeID, cp, &ID, &err);
    checkError(err, nodeID,
               string("Error getting node ID ") + string(extendedData));

    ADF_Read_All_Data(ID, local_buffer, &err);
    checkError(err, nodeID, "Error reading data");

    size = getSize(ID, true);
    local_buffer += type2size(type) * size;
  }

  delete[] names;
}

void readBlockData(double const nodeID, char *buffer, const int64 &blockStart,
                   const int64 &blockEnd) {
  int64 totalNodeSize = getSize(nodeID);
  int64 nodeSizeNoExtended = getSize(nodeID, true);
  int64 readSize = 0;
  Loci::debugout << "total node size is " << totalNodeSize
                 << " and without extended data it is " << nodeSizeNoExtended
                 << endl;
  Loci::debugout << "reading block from " << blockStart << " to " << blockEnd
                 << endl;
  int err;
  if (blockStart <= nodeSizeNoExtended) {  // starting from parent block
    if (blockEnd <= nodeSizeNoExtended) {  // all data on parent block
      Loci::debugout << "all data on parent block, don't need extended" << endl;
      ADF_Read_Block_Data(nodeID, blockStart, blockEnd, buffer, &err);
      readSize = blockEnd - blockStart + 1;
    } else {  // some data on extended block
      Loci::debugout << "some data on extended block" << endl;
      ADF_Read_Block_Data(nodeID, blockStart, nodeSizeNoExtended, buffer, &err);
      readSize = nodeSizeNoExtended - blockStart + 1;
    }
    checkError(err, nodeID, "Error reading block data");
  }
  double tmpID;
  ADF_Get_Node_ID(nodeID, "ExtendedSize", &tmpID, &err);
  // 29 - Specified child is not a child of the specified parent
  if (err == 29 || blockEnd <= nodeSizeNoExtended) {
    return;
  }

  // is extended data -- read the rest of the extended data
  Loci::debugout << "Reading 'ExtendedSize' for node" << endl;
  string type = getType(nodeID);
  char *local_buffer = buffer + type2size(type) * readSize;

  int num_children;
  ADF_Number_of_Children(nodeID, &num_children, &err);
  Loci::debugout << "extended data number of children " << num_children << endl;
  checkError(err, nodeID, "Error getting number of children");

  int inum_ret;
  char *names = new char[(ADF_NAME_LENGTH_C + 1) * num_children];
  ADF_Children_Names(nodeID, 1, num_children, ADF_NAME_LENGTH_C, &inum_ret,
                     names, &err);
  checkError(err, nodeID, "Error getting children names");

  char *cp;
  char *end = names + (ADF_NAME_LENGTH_C + 1) * inum_ret;
  Loci::debugout << "read data of size " << readSize << " from parent" << endl;
  int64 runningTotal = nodeSizeNoExtended;
  Loci::debugout << "position in total node size is " << runningTotal << endl;
  for (cp = names; cp != end; cp += (ADF_NAME_LENGTH_C + 1)) {
    static char const *extendedData = "ExtendedData-";
    if (strncmp(cp, extendedData, strlen(extendedData))) continue;

    double ID;
    ADF_Get_Node_ID(nodeID, cp, &ID, &err);
    checkError(err, nodeID,
               string("Error getting node ID ") + string(extendedData));

    // adjust start if some data was read from parent
    Loci::debugout << "on child " << cp << endl;
    int64 childSize = getSize(ID);
    if (blockStart > runningTotal + childSize) {  // start on next child
      runningTotal += childSize;
      Loci::debugout
          << "adjusting position in total node size to advance to next child "
          << runningTotal << endl;
      continue;
    }

    int64 childStart = max(blockStart - runningTotal, 1L);
    int64 childEnd = min(blockEnd - runningTotal, childSize);
    Loci::debugout << "reading data for child " << cp << " from " << childStart
                   << " to " << childEnd << " and child is of size "
                   << childSize << endl;
    ADF_Read_Block_Data(ID, childStart, childEnd, local_buffer, &err);
    checkError(err, nodeID, "Error reading data");

    readSize = childEnd - childStart + 1;
    Loci::debugout << "read data of size " << readSize << " from child" << endl;
    local_buffer += type2size(type) * readSize;
    if (blockEnd == childEnd + runningTotal)  { // done reading chunk
      Loci::debugout << "done reading block" << endl;
      break;
    } else if (childEnd == childSize) {  // advance child
      Loci::debugout << "advancing to next child" << endl;
      runningTotal += childSize;
    }
  }

  delete[] names;
}

string readNodestr(double parentID, const char *name) {
  double childID;
  int err;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  if (err != -1) return string();  // node doesn't exist

  int64 size = getSize(childID);
  if (size == 0) return string();  // node exist and size is 0

  char value[ADF_FILENAME_LENGTH];
  readAllData(childID, &value[0]);
  value[size] = '\0';

  return string(value);
}

vector<float> readCoordinatesBlock32(double parentID, const char *name,
                                     const int64 &blockStart,
                                     const int64 &blockEnd) {
  double nodeID;
  int err;

  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID,
             string("Error in getting node ID of ") + string(name));

  int64 size = blockEnd - blockStart + 1;
  if (size == 0) return vector<float>();
  vector<float> data(size);
  readBlockData(nodeID, reinterpret_cast<char *>(&(*data.begin())), blockStart,
                blockEnd);
  return data;
}


int64 readNodeis(double parentID, const char *name, int32 **value) {
  double nodeID;
  int err;

  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID,
             string("Error in getting node ID of ") + string(name));

  int64 size = getSize(nodeID);
  if (size == 0) return 0;
  *value = new int32[size];
  readAllData(nodeID, (char *)(*value));
  return size;
}

void readFacesBlock32(const double &nodeID, const int64 &nodeStart,
                      const int64 &nodeEnd, vector<int32> &data) {
  int64 size = nodeEnd - nodeStart + 1;
  if (size == 0) return;
  readBlockData(nodeID, reinterpret_cast<char *>(&(*(data.begin()))), nodeStart,
                nodeEnd);
}

void readFacesConnBlock32(const double &nodeID, const int64 &numFaces,
                          const int64 &nodeStart, int64 &nodeEnd,
                          int64 &numFacesRead, vector<int32> &data,
                          const int faceSizeGuess = 7) {
  Loci::debugout << "readFacesConnBlock32()" << endl;
  const int64 nodeSize = getSize(nodeID);
  Loci::debugout << "node size is " << nodeSize << endl;
  nodeEnd = min(nodeStart + faceSizeGuess * numFaces - 1, nodeSize);
  int64 size = nodeEnd - nodeStart + 1;
  Loci::debugout << "node start, end, guess size is " << nodeStart << " "
                 << nodeEnd << " " << size << endl;
  data.resize(size);
  if (size == 0) return;
  readBlockData(nodeID, reinterpret_cast<char *>(&(*(data.begin()))), nodeStart,
                nodeEnd);
  Loci::debugout << "last data read is " << data.back() << endl;
  // find out how many faces we actually read
  numFacesRead = 0;
  int64 index = 0;
  bool readingFaces = true;
  while (readingFaces) {
    int32 faceSize = data[index];
  
    if (faceSize > 20) {
      Loci::debugout
          << "WARNING face size is " << faceSize
          << " is this a polyhedral mesh? Is a face of this size expected?"
          << endl;
    }

    if (index + faceSize < size) {  // have whole face
      numFacesRead++;
      index += faceSize + 1;  // advance to next face
    } else {  // have partial face
      readingFaces = false;
    }

    if (numFacesRead == numFaces) {  // read alloted faces, stop
      readingFaces = false;
    } else if (index >= static_cast<int64>(data.size() - 1)) {  // read all data
      readingFaces = false;
    }
  }

  if (index > static_cast<int64>(data.size())) {
    cerr << "ERROR index is bigger than data size read" << endl;
    Loci::Abort();
  }

  if (numFacesRead < numFaces && nodeSize > nodeEnd) {
    // our guess was too small, need to try again and read more
    // instead of copy/concatenating vectors, just reading from start again
    int64 numShort = numFaces - numFacesRead;
    Loci::debugout << "guess was too small by " << numShort
                   << " faces, try again" << endl;
    int newFaceSizeGuess = data.size() / numFacesRead + 2;
    Loci::debugout << "face size guess is now " << newFaceSizeGuess << endl;
    data.clear();
    readFacesConnBlock32(nodeID, numFaces, nodeStart, nodeEnd, numFacesRead,
                         data, newFaceSizeGuess);
  } else {
    // read too many, just the right number of faces, or read all faces on node
    // resize vector to ensure just the correct amount of data is kept
    Loci::debugout << "read too many, the right number, or all faces on node; "
                      "resizing from "
                   << data.size() << " to " << index << endl;
    Loci::debugout << "read " << numFacesRead << " faces, targeting "
                   << numFaces << " faces " << endl;

    data.resize(index);
    nodeEnd = nodeStart + index - 1;
    Loci::debugout << "now last data read is " << data.back() << endl;
    Loci::debugout << "nodeEnd is " << nodeEnd << endl;
  }
}

vector<double> readCoordinatesBlock64(double parentID, const char *name,
                                      const int64 &blockStart,
                                      const int64 &blockEnd) {
  double nodeID;
  int err;

  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID,
             string("Error in getting node ID of ") + string(name));

  int64 size = blockEnd - blockStart + 1;
  if (size == 0) return vector<double>();
  vector<double> data(size);
  readBlockData(nodeID, reinterpret_cast<char *>(&(*data.begin())), blockStart,
                blockEnd);
  return data;
}

int32 readNodei32(double parentID, const char *name) {
  double childID;
  int err;
  int32 value;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID,
             string("Error in  getting node ID of ") + string(name));
  ADF_Read_All_Data(childID, (char *)&value, &err);
  checkError(err, childID, "Error reading data");
  return value;
}

float readNodef(double parentID, const char *name) {
  double childID;
  int err;
  float value;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID,
             string("Error in getting  node ID of ") + string(name));
  ADF_Read_All_Data(childID, (char *)&value, &err);
  checkError(err, childID, "Error reading data");
  return value;
}

string getDataType(double parentID, const char *name) {
  double childID;
  int err;

  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID,
             string("Error in  getting node ID of ") + string(name));
  return getType(childID);
}

int getNumNodes(double vertices) {
  double nodeID;
  int err;
  int64 num_dims, dim_values[ADF_MAX_DIMENSIONS];

  ADF_Get_Node_ID(vertices, "Coordinates", &nodeID, &err);
  checkError(err, vertices, "Error in getting node ID of Coordinates");
  getDimensions(nodeID, &num_dims, dim_values);

  if (num_dims != 2) return 0;
  return dim_values[1];
}

int getNumFaces(double topoID) {
  int numFaces = readNodei32(topoID, "InternalFaces/NumFaces");
  int err;
  int numChildren;
  ADF_Number_of_Children(topoID, &numChildren, &err);
  checkError(err, topoID,
             "Cannot get number of children of the first child of States");

  char name[ADF_NAME_LENGTH_C] = {'\0'};
  int num_ret = 0;

  for (int start = 1; start <= numChildren; start++) {
    ADF_Children_Names(topoID, start, 1, ADF_NAME_LENGTH_C, &num_ret, name, &err);
    checkError(err, topoID, string("Cannot get node ID of ") + string(name));
    double bfaceID;
    if (string(name).substr(0, 14) == "BoundaryFaces-") {
      ADF_Get_Node_ID(topoID, name, &bfaceID, &err);
      checkError(err, topoID, string("Cannot get node ID of ") + string(name));
      numFaces += readNodei32(bfaceID, "NumFaces");
    }
  }

  return numFaces;
}

vector<std::pair<double, int> > getBoundaryFaceIDs(double topoID) {
  int err;
  int numChildren;
  ADF_Number_of_Children(topoID, &numChildren, &err);
  checkError(err, topoID,
             "Cannot get number of children of the first child of States");

  char name[ADF_NAME_LENGTH_C] = {'\0'};
  int num_ret = 0;

  int numBoundaries = readNodei32(topoID, "NumBoundaryTypes");

  vector<std::pair<double, int> > boundaryIDs;
  for (int start = 1; start <= numChildren; start++) {
    ADF_Children_Names(topoID, start, 1, ADF_NAME_LENGTH_C, &num_ret, name,
                       &err);
    checkError(err, topoID, string("Cannot get node ID of ") + string(name));
    double bfaceID;
    if (string(name).substr(0, 14) == "BoundaryFaces-") {
      ADF_Get_Node_ID(topoID, name, &bfaceID, &err);
      checkError(err, topoID, string("Cannot get node ID of ") + string(name));
      string index = string(name).substr(14);
      int btype = atoi(index.c_str());
      if (btype == 0) btype = numBoundaries;  // some index start with 0

      std::pair<double, int> bdata(bfaceID, btype);
      boundaryIDs.push_back(bdata);
    }
  }
  return boundaryIDs;
}

// Read in a block of nodal coordinates in a given processor-id
// Nodes reside on a given processor-id depending on the partition in STAR-CCM+
// A block of nodes will be read and sent to a given Loci processor
void readCoordinatesBlock(const vector<double> &processorIDs,
                          store<vector3d<double> > &pos, size_t &pid_start,
                          int64 &pos_start, int64 &node_start) {
  Loci::debugout << "reading coordinates block at proc-id " << pid_start
                 << endl;
  Loci::debugout << "pos_start is " << pos_start << endl;
  Loci::debugout << "node_start is " << node_start << endl;
  double vertices;
  getMeshID(processorIDs[pid_start], &vertices, NULL);
  string type = getDataType(vertices, "Coordinates");

  // nodes on processor-id
  const int64 numNodes = getNumNodes(vertices);
  Loci::debugout << "have " << numNodes << " vertices on proc-id" << endl;
  // num nodes left to read on processor-id, and in block
  const int64 leftToReadInProc = numNodes - pos_start;
  const int64 leftToReadInBlock = pos.domain()[0].second - pos_start + 1;
  // ending point on processor-id
  const int64 readSize = min(leftToReadInBlock, leftToReadInProc);
  Loci::debugout << "read size is " << readSize << endl;
  Loci::debugout << "left to read in proc " << leftToReadInProc << endl;
  Loci::debugout << "left to read in block " << leftToReadInBlock << endl;
  const int64 pos_end = pos_start + readSize;
  const int64 node_end = node_start + 3 * readSize - 1;

  float scale = readNodef(vertices, "ScaleFactor");

  if (type == "R4") {  // coordinates are 32 bit reals
    std::vector<float> verts =
        readCoordinatesBlock32(vertices, "Coordinates", node_start, node_end);
    for (int64 i = pos_start; i < pos_end; ++i) {
      int64 vid = i - pos_start;
      vector3d<double> p;
      p.x = verts[3 * vid] * scale;
      p.y = verts[3 * vid + 1] * scale;
      p.z = verts[3 * vid + 2] * scale;
      pos[size_t(i)] = p;
    }
  } else {  // coordinates are 64 bit reals
    std::vector<double> verts =
        readCoordinatesBlock64(vertices, "Coordinates", node_start, node_end);
    for (int64 i = pos_start; i < pos_end; ++i) {
      int64 vid = i - pos_start;
      vector3d<double> p;
      p.x = verts[3 * vid] * static_cast<double>(scale);
      p.y = verts[3 * vid + 1] * static_cast<double>(scale);
      p.z = verts[3 * vid + 2] * static_cast<double>(scale);
      pos[size_t(i)] = p;
    }
  }

  node_start = node_end + 1;
  pos_start = pos_end;
  if (leftToReadInBlock > leftToReadInProc) {
    // read all data in processor-id, but still need more data for block
    pid_start++;
    readCoordinatesBlock(processorIDs, pos, pid_start, pos_start, node_start);
  } else if (leftToReadInBlock == leftToReadInProc) {
    // read all data in processor-id, and filled block
    pid_start++;
  }
}

void determineBlockLocationCells(
    const int &processorID,
    const vector<vector<std::pair<double, int> > > &boundaryIDs,
    const int64 &blockSize, const int64 &facesRead, const double &nodeID,
    int64 &nodeStart, int64 &nodeEnd, int &boundaryIndex,
    bool &doneReadingBlock, bool &advanceProcessor, bool &advanceNode) {
  // nodeID/name are initialized to read internal faces first
  // block start/end are initialized to 0, or the start/end of previous block

  // determine if there is space in the current node
  const int64 nodeSize = getSize(nodeID);
  const bool onBoundary = boundaryIndex >= 0;
  const int64 leftToReadInNode = nodeSize - nodeEnd;
  // divide by 2 because we will augment the data on the boundary by 2x when
  // adding the boundary id
  const int64 leftToReadInBlock =
      onBoundary ? blockSize / 2 - facesRead : blockSize - 2 * facesRead;

  const int64 readSize = min(leftToReadInBlock, leftToReadInNode);

  Loci::debugout << "node size is " << nodeSize << endl;
  Loci::debugout << "left to read for current node " << leftToReadInNode
                 << endl;
  Loci::debugout << "left to read for current block " << leftToReadInBlock
                 << endl;
  Loci::debugout << "onBoundary, readSize: " << onBoundary << " " << readSize
                 << endl;
  nodeStart = nodeEnd + 1;
  nodeEnd += readSize;
  advanceProcessor = false;
  if (leftToReadInNode > leftToReadInBlock) {
    // read all data from current node, and stay on node, done with block
    doneReadingBlock = true;
    advanceNode = false;
  } else if (leftToReadInNode == leftToReadInBlock) {
    // read all data from current node, go to next node, done with block
    boundaryIndex++;
    if (boundaryIndex == static_cast<int>(boundaryIDs[processorID].size())) {
      advanceProcessor = true;
    }
    doneReadingBlock = true;
    advanceNode = true;
  } else {
    // read all data from current node, go to next node, not done with block
    boundaryIndex++;
    if (boundaryIndex == static_cast<int>(boundaryIDs[processorID].size())) {
      advanceProcessor = true;
    }
    doneReadingBlock = false;
    advanceNode = true;
  }
}

std::map<int, int> get_cell_g2l_map(const double &topoID, bool &isLocal) {
  std::map<int, int> cell_g2l_map;
  int err;
  double root;
  ADF_Get_Root_ID(topoID, &root, &err);
  checkError(err, topoID, "Error getting root ID");
  // first check if the node indexes and cell indexes are local or not
  isLocal = false;
  if (readNodestr(root, "/Meshes/Space") == "Local") isLocal = true;

  // if the cell indexes are not local read in the cell map data
  if (!isLocal) {
    int32 *cellMapData = 0;
    char nodeName[ADF_NAME_LENGTH_C] = {'\0'};
    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "Cells/MapId");
    int mapIndex = readNodei32(topoID, nodeName);

    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "/Maps/Map-%d/IdMap", mapIndex);
    int64 localNumCells = readNodeis(root, nodeName, &cellMapData);

    if (cellMapData == 0) {
      cerr << " Error reading vertex map " << endl;
      exit(1);
    }

    // build global2local map
    for (int i = 0; i < localNumCells; i++) {
      cell_g2l_map[cellMapData[i]] = i + 1;
    }
    if (cellMapData) delete[] cellMapData;
  }
  return cell_g2l_map;
}

void readFace2CellInChunks(
    const vector<double> &processorIDs,
    const vector<vector<std::pair<double, int> > > &boundaryIDs,
    const vector<entitySet> &local_faces, Map &cl, Map &cr) {
  Loci::debugout << "readFace2CellInChunks()" << endl;
  // get rank and number of processors
  const int R = Loci::MPI_rank;
  const int P = Loci::MPI_processes;

  vector<int32> faceCellData;
  if (R == 0) {
    cout << "Reading face-to-cell maps..." << endl;
    // for each Loci proc, read data chunk and send to corresponding proc
    int pidx = -1;  // processor id index
    int boundaryIndex = -1;  // start with internal faces
    int64 nodeStart = 0;
    int64 nodeEnd = 0;
    string name = "InternalFaces/Cells";
    bool advanceProcessor = true;
    bool advanceNode = false;
    double nodeID = 0;
    std::map<int, int> cell_g2l_map;

    for (int pp = 0; pp < P; ++pp) {
      int64 facesRead = 0;
      // 2x because each face has a left and right cell
      const int64 blockSize = 2 * local_faces[pp].size();
      // outer vector for mesh boundaries (interior treated as boundary)
      vector<vector<int32> > faceCellBlockBnd;
      Loci::debugout << "\nreading chunk for processor " << pp << " of size "
                     << blockSize << " for faces " << local_faces[pp] << endl;
      bool doneReadingBlock = false;
      bool isLocal = false;

      while (!doneReadingBlock) {
        Loci::debugout << "have read " << facesRead << " of "
             << local_faces[pp].size() << " faces" << endl;
        // get mesh node for current processor
        // if we have moved to a new processor, get new mesh id and map
        if (advanceProcessor) {
          pidx++;
          Loci::debugout << "advancing processor to index " << pidx << endl;
          // reset variables
          string name = "InternalFaces/Cells";
          boundaryIndex = -1;
          advanceProcessor = false;
          advanceNode = false;
          double topoID;
          getMeshID(processorIDs[pidx], NULL, &topoID);
          int err;
          ADF_Get_Node_ID(topoID, name.c_str(), &nodeID, &err);
          checkError(err, topoID,
                     string("Error in getting node ID of ") + name);
          cell_g2l_map = get_cell_g2l_map(topoID, isLocal);
          Loci::debugout << "cell global-to-local map is size "
                         << cell_g2l_map.size() << endl;
        }
        if (advanceNode) {
          Loci::debugout << "advancing node" << endl;
          int err;
          string name = "Cells";
          nodeStart = 0;
          nodeEnd = 0;
          Loci::debugout << "getting node ID for boundary at proc-index "
                         << pidx << " and boundary index " << boundaryIndex
                         << endl;
          ADF_Get_Node_ID(boundaryIDs[pidx][boundaryIndex].first, name.c_str(),
                          &nodeID, &err);
          checkError(err, boundaryIDs[pidx][boundaryIndex].first,
                     string("Error in getting node ID of ") + name);
          advanceNode = false;
        }

        bool onBoundary = boundaryIndex >= 0;
        int32 btype = 0;
        if (onBoundary) {  
          btype = boundaryIDs[pidx][boundaryIndex].second;
        }
        determineBlockLocationCells(
            pidx, boundaryIDs, blockSize, facesRead, nodeID, nodeStart, nodeEnd,
            boundaryIndex, doneReadingBlock, advanceProcessor, advanceNode);
        Loci::debugout << "DEBUG on processor-id " << pidx
                       << " and boundary id " << boundaryIndex
                       << " reading from node " << nodeID << " at " << nodeStart
                       << " to " << nodeEnd << endl;
        faceCellBlockBnd.push_back(vector<int32>(nodeEnd - nodeStart + 1));
        readFacesBlock32(nodeID, nodeStart, nodeEnd, faceCellBlockBnd.back());
        Loci::debugout << "done reading size " << faceCellBlockBnd.back().size()
                       << endl;
        // map cell data if not local
        if (!isLocal) {
          for (unsigned int i = 0; i < faceCellBlockBnd.back().size(); ++i) {
            if (cell_g2l_map.find(faceCellBlockBnd.back()[i]) ==
                cell_g2l_map.end()) {
              cerr << "ERROR: problem with cell_g2l_map" << endl;
              cerr << "face index, cell id: " << i << " "
                   << faceCellBlockBnd.back()[i] << endl;
              Loci::Abort();
            }
            faceCellBlockBnd.back()[i] =
                cell_g2l_map[faceCellBlockBnd.back()[i]];
          }
        }

        // if just read from boundary, rearrange cell data to include cr
        if (onBoundary) {
          vector<int32> faceCellsUpdate(2 * faceCellBlockBnd.back().size());
          for (unsigned int i = 0; i < faceCellBlockBnd.back().size(); ++i) {
            faceCellsUpdate[2 * i] = faceCellBlockBnd.back()[i];
            faceCellsUpdate[2 * i + 1] = -btype;
          }
          faceCellBlockBnd.back().swap(faceCellsUpdate);
        }
        facesRead += faceCellBlockBnd.back().size() / 2;
      }

      // send data to appropriate Loci processor
      int64 totalSize = 0;
      vector<int> starts(faceCellBlockBnd.size());
      for (unsigned int i = 0; i < faceCellBlockBnd.size(); ++i) {
        starts[i] = totalSize;
        totalSize += faceCellBlockBnd[i].size();
      }
      // check that total size matches expected size
      if (totalSize != blockSize) {
        cerr << "ERROR in reading face-cell data. Expected to read chunk size "
                "of "
             << blockSize << ", but read chunk size of " << totalSize << endl;
        Loci::Abort();
      }

      if (pp == 0) {
        Loci::debugout << "data is local for proc 0, concatenating" << endl;
        faceCellData.resize(blockSize);
        for (unsigned int i = 0; i < faceCellBlockBnd.size(); ++i) {
          std::copy(faceCellBlockBnd[i].begin(), faceCellBlockBnd[i].end(),
                    faceCellData.begin() + starts[i]);
        }
        Loci::debugout << "done with concatenation" << endl;
      } else {  // send to corresponding processor
        int numChunks = faceCellBlockBnd.size();
        Loci::debugout << "sending faceCell data in " << numChunks
                       << " chunks to processor " << pp << endl;
        MPI_Send(&numChunks, 1, MPI_INT, pp, 1, MPI_COMM_WORLD);
        MPI_Send(&(*starts.begin()), numChunks, MPI_INT, pp, 2, MPI_COMM_WORLD);
        for (unsigned int i = 0; i < faceCellBlockBnd.size(); ++i) {
          Loci::debugout << "sending faceCell data of size "
                         << faceCellBlockBnd[i].size() << " to rank " << pp
                         << endl;
          MPI_Send(&(*faceCellBlockBnd[i].begin()), faceCellBlockBnd[i].size(),
                   MPI_INT32_T, pp, i + 3, MPI_COMM_WORLD);
        }
      }
    }
  } else {  // not rank 0 - receive data
    MPI_Status status;
    faceCellData.resize(2 * local_faces[R].size());
    Loci::debugout << "processor " << R << " receiving faceCell data of size "
                   << faceCellData.size() << " for faces " << local_faces[R]
                   << endl;
    int numChunks = 0;
    MPI_Recv(&numChunks, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    vector<int> starts(numChunks);
    MPI_Recv(&(*starts.begin()), numChunks, MPI_INT, 0, 2, MPI_COMM_WORLD,
             &status);
    for (unsigned int i = 0; i < starts.size(); ++i) {
      int dataSize = faceCellData.size() - starts[i];
      if (i != starts.size() - 1) {  // not on last receive
        dataSize = starts[i + 1] - starts[i];
      }
      Loci::debugout << "receiving data of size " << dataSize << endl;
      MPI_Recv(&(*(faceCellData.begin() + starts[i])), dataSize, MPI_INT32_T, 0,
               i + 3, MPI_COMM_WORLD, &status);
    }
  }

  // put data into loci containers
  int cnt = 0;
  cl.allocate(local_faces[R]);
  cr.allocate(local_faces[R]);
  FORALL(local_faces[R], ei) {
    cl[ei] = faceCellData[cnt++];
    cr[ei] = faceCellData[cnt++];
  }ENDFORALL;

  if (R == 0) {
    cout << "  ...done reading face-to-cell maps" << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

std::map<int, int> get_node_g2l_map(const double &topoID, bool &isLocal) {
  std::map<int, int> node_g2l_map;
  int err;
  double root;
  ADF_Get_Root_ID(topoID, &root, &err);
  checkError(err, topoID, "Error getting root ID");
  // first check if the node indexes and cell indexes are local or not
  isLocal = false;
  if (readNodestr(root, "/Meshes/Space") == "Local") isLocal = true;

  // if the node indexes are not local, read in the vertices map data
  if (!isLocal) {
    int32 *vertexMapData = 0;
    char nodeName[ADF_NAME_LENGTH_C] = {'\0'};
    ADF_Get_Name(topoID, nodeName, &err);
    checkError(err, topoID, "Error getting node name");
    int index = atoi(string(nodeName).substr(18).c_str());
    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "/Meshes/Vertices-%d/MapId", index);
    int mapIndex = readNodei32(root, nodeName);

    // read vertices map data, assume FaceBasedTopology-i and  Vertices-i share
    // the same vertex map
    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "/Maps/Map-%d/IdMap", mapIndex);

    int64 localNumNode = readNodeis(root, nodeName, &vertexMapData);
    if (vertexMapData == 0) {
      cerr << " Error reading vertex map " << endl;
      exit(1);
    }

    // build global2local map
    for (int i = 0; i < localNumNode; i++) {
      node_g2l_map[vertexMapData[i]] = i + 1;
    }
    if (vertexMapData) delete[] vertexMapData;
  }
  return node_g2l_map;
}

// read in a block of facebased topology
void readFace2NodeInChunks(
    const vector<double> &processorIDs,
    const vector<vector<std::pair<double, int> > > &boundaryIDs,
    const vector<entitySet> &local_faces, multiMap &face2node) {
  Loci::debugout << "readFace2NodeInChunks()" << endl;
  // get rank and number of processors
  const int R = Loci::MPI_rank;
  const int P = Loci::MPI_processes;

  vector<int32> faceNodeData;
  if (R == 0) {
    cout << "Reading face connectivity..." << endl;
    // for each Loci proc, read data chunk and send to corresponding proc
    int pidx = -1;  // processor id index
    int boundaryIndex = -1;  // start with internal faces
    int64 nodeStart = 0;
    int64 nodeEnd = 0;
    string name = "InternalFaces/Vertices";
    bool advanceProcessor = true;
    bool advanceNode = false;
    double nodeID = 0;
    std::map<int, int> node_g2l_map;

    for (int pp = 0; pp < P; ++pp) {
      int64 facesRead = 0;
      const int64 facesInBlock = local_faces[pp].size();
      // outer vector for mesh boundaries (interior treated as boundary)
      vector<vector<int32> > faceNodeBlockBnd;
      Loci::debugout << "\nreading chunk for processor " << pp << " of size "
                     << facesInBlock << " for faces " << local_faces[pp]
                     << endl;
      bool isLocal = false;

      while (facesRead != facesInBlock) {
        Loci::debugout << "have read " << facesRead << " of " << facesInBlock
                       << " faces" << endl;
        // get mesh node for current processor
        // if we have moved to a new processor, get new mesh id and map
        if (advanceProcessor) {
          pidx++;
          Loci::debugout << "advancing processor to index " << pidx << endl;
          // reset variables
          string name = "InternalFaces/Vertices";
          boundaryIndex = -1;
          advanceProcessor = false;
          advanceNode = false;
          double topoID;
          getMeshID(processorIDs[pidx], NULL, &topoID);
          int err;
          ADF_Get_Node_ID(topoID, name.c_str(), &nodeID, &err);
          checkError(err, topoID,
                     string("Error in getting node ID of ") + name);
          node_g2l_map = get_node_g2l_map(topoID, isLocal);
          Loci::debugout << "node global-to-local map is size "
                         << node_g2l_map.size() << endl;
        }
        if (advanceNode) {
          Loci::debugout << "advancing node" << endl;
          int err;
          string name = "Vertices";
          nodeStart = 0;
          nodeEnd = 0;
          Loci::debugout << "getting node ID for boundary at proc-index "
                         << pidx << " and boundary index " << boundaryIndex
                         << endl;
          ADF_Get_Node_ID(boundaryIDs[pidx][boundaryIndex].first, name.c_str(),
                          &nodeID, &err);
          checkError(err, boundaryIDs[pidx][boundaryIndex].first,
                     string("Error in getting node ID of ") + name);
          advanceNode = false;
        }

        int64 facesToRead = facesInBlock - facesRead;
        nodeStart = nodeEnd + 1;
        Loci::debugout << "attempting to read " << facesToRead << " faces "
                       << endl;
        int64 facesReadInBlock = 0;
        faceNodeBlockBnd.push_back(vector<int32>());
        readFacesConnBlock32(nodeID, facesToRead, nodeStart, nodeEnd,
                             facesReadInBlock, faceNodeBlockBnd.back());
        Loci::debugout << "done reading size " << faceNodeBlockBnd.back().size()
                       << " which is " << facesReadInBlock << " faces" << endl;
        Loci::debugout << "after read node start, end " << nodeStart << ", "
                       << nodeEnd << endl;
        if (nodeEnd == getSize(nodeID)) {
          advanceNode = true;
          boundaryIndex++;
        }
        if (boundaryIndex == static_cast<int>(boundaryIDs[pidx].size())) {
          advanceProcessor = true;
        }

        // map node data if not local
        if (!isLocal) {
          // map vertices data
          unsigned int pointer = 0;
          while (pointer < faceNodeBlockBnd.back().size()) {
            int fs = faceNodeBlockBnd.back()[pointer];
            if (fs > 20) {
              Loci::debugout
                  << "WARNING have face size of " << fs
                  << " is this a polyhedral mesh? Is this face size expected?"
                  << endl;
              Loci::debugout << "pointer is " << pointer << " size is "
                             << faceNodeBlockBnd.back().size() << endl;
            }

            for (int i = 1; i <= faceNodeBlockBnd.back()[pointer]; ++i) {
              if (node_g2l_map.find(faceNodeBlockBnd.back()[pointer + i]) ==
                  node_g2l_map.end()) {
                cerr << "ERROR: problem with node_g2l_map" << endl;
                cerr << "connectivity index, pointer, node id: " << i << " "
                     << pointer << " " << faceNodeBlockBnd.back()[pointer + i]
                     << endl;
                Loci::Abort();
              }
              faceNodeBlockBnd.back()[pointer + i] =
                  node_g2l_map[faceNodeBlockBnd.back()[pointer + i]];
            }
            pointer += faceNodeBlockBnd.back()[pointer] + 1;
          }
        }

        Loci::debugout << "read " << facesReadInBlock << " faces in last block"
                       << endl;
        Loci::debugout << "for next read advanceNode: " << advanceNode
                       << " advance processor " << advanceProcessor << endl;
        facesRead += facesReadInBlock;
      }

      // send data to appropriate Loci processor
      int64 totalSize = 0;
      vector<int> starts(faceNodeBlockBnd.size());
      for (unsigned int i = 0; i < faceNodeBlockBnd.size(); ++i) {
        starts[i] = totalSize;
        totalSize += faceNodeBlockBnd[i].size();
      }

      if (pp == 0) {
        Loci::debugout << "data is local for proc 0, concatenating" << endl;
        faceNodeData.resize(totalSize);
        for (unsigned int i = 0; i < faceNodeBlockBnd.size(); ++i) {
          std::copy(faceNodeBlockBnd[i].begin(), faceNodeBlockBnd[i].end(),
                    faceNodeData.begin() + starts[i]);
        }
        Loci::debugout << "done with concatenation" << endl;
      } else {  // send to corresponding processor
        int numChunks = faceNodeBlockBnd.size();
        Loci::debugout << "sending faceNode data in " << numChunks
                       << " chunks to processor " << pp << endl;
        MPI_Send(&totalSize, 1, MPI_INT64_T, pp, 0, MPI_COMM_WORLD);
        MPI_Send(&numChunks, 1, MPI_INT, pp, 1, MPI_COMM_WORLD);
        MPI_Send(&(*starts.begin()), numChunks, MPI_INT, pp, 2, MPI_COMM_WORLD);
        for (unsigned int i = 0; i < faceNodeBlockBnd.size(); ++i) {
          Loci::debugout << "sending faceNode data of size "
                         << faceNodeBlockBnd[i].size() << " to rank " << pp
                         << endl;
          MPI_Send(&(*faceNodeBlockBnd[i].begin()), faceNodeBlockBnd[i].size(),
                   MPI_INT32_T, pp, i + 3, MPI_COMM_WORLD);
        }
      }
    }
  } else {  // not rank 0 - receive data
    MPI_Status status;
    int64 totalSize = 0;
    MPI_Recv(&totalSize, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD, &status);
    faceNodeData.resize(totalSize);
    Loci::debugout << "processor " << R << " receiving faceNode data of size "
                   << faceNodeData.size() << " for faces " << local_faces[R]
                   << endl;
    int numChunks = 0;
    MPI_Recv(&numChunks, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    vector<int> starts(numChunks);
    MPI_Recv(&(*starts.begin()), numChunks, MPI_INT, 0, 2, MPI_COMM_WORLD,
             &status);
    for (unsigned int i = 0; i < starts.size(); ++i) {
      int dataSize = faceNodeData.size() - starts[i];
      if (i != starts.size() - 1) {  // not on last receive
        dataSize = starts[i + 1] - starts[i];
      }
      Loci::debugout << "receiving data of size " << dataSize << endl;
      MPI_Recv(&(*(faceNodeData.begin() + starts[i])), dataSize, MPI_INT32_T, 0,
               i + 3, MPI_COMM_WORLD, &status);
    }
  }

  // put data into loci containers on all processors
  int cnt = 0;
  store<int> nodes_per_face;
  nodes_per_face.allocate(local_faces[R]);
  FORALL(local_faces[R], ei) {
    nodes_per_face[ei] = faceNodeData[cnt];
    cnt += nodes_per_face[ei] + 1;
  }ENDFORALL;

  cnt = 0;
  face2node.allocate(nodes_per_face);
  FORALL(local_faces[R], ei) {
    for (int ii = 0; ii < nodes_per_face[ei]; ++ii) {
      face2node[ei][ii] = faceNodeData[cnt + ii + 1] - 1;
    }
    cnt += nodes_per_face[ei] + 1;
  }ENDFORALL;

  MPI_Barrier(MPI_COMM_WORLD);
  if (R == 0) {
    cout << "  ...done reading face connectivity" << endl;
  }
}

// function to remove spaces in strings
void removeSpace(string &str) {
  str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
}

// Function the check that all boundary names are unique, and rename them
// if they are not.
void checkUnique(const vector<pair<int, string> > &surf_ids, int &bcIndex,
                 string &bcLabel) {
  // check to make sure bc labels are unique
  for (unsigned int ii = 0; ii < surf_ids.size(); ii++) {
    if (bcLabel == surf_ids[ii].second) {
      string newName = bcLabel + string("_1");
      cerr << "WARNING: Boundary ID " << surf_ids[ii].first
           << " has the same name as " << bcIndex << "!\nBoth are named "
           << bcLabel << ". Boundary ID " << bcIndex << " will be renamed to "
           << newName << endl;
      bcLabel = newName;
    }
  }
}

// This function gets the boundary condition names and IDs.
// These values are contained in the ProblemDescriptions node 
// (there is only one of these). Within the ProblemDescriptions node there 
// is one node for each set of problem description information of the form 
// ProblemDescription-#. Within each ProblemDescription-# node there is a 
// BoundaryRegion-# node to describe each boundary region. The 
// BoundaryRegion-# nodes have a Label field which contains the boundary
// label as named in STAR-CCM+; the # is used as the boundary ID
void getSurfaceNames(double &probDescID, vector<pair<int, string> > &surf_ids) {
  // probDescID -- ID for the ProblemDescriptions node
  // surf_ids -- (output) vector of boundary condtion IDs and names

  int err = 0;  // initialize error flag

  // get number of ProblemDescription children
  int numChildren;
  ADF_Number_of_Children(probDescID, &numChildren, &err);
  checkError(err, probDescID,
             "Cannot get number of children of the ProblemDescriptions");

  char name[ADF_NAME_LENGTH_C] = {'\0'};
  int num_ret = 0;

  // loop over all children to get to the boundary conditions held by each child
  // if there is more than one child the BCs for the first child will be stored,
  // then the BCs for the second child, and so on
  for (int start = 1; start <= numChildren; start++) {
    ADF_Children_Names(probDescID, start, 1, ADF_NAME_LENGTH_C, &num_ret, name,
                       &err);
    checkError(err, probDescID, "Cannot get node ID of " + string(name));

    // if child is ProblemDescription-# probe further for boundary condition
    // data
    double pDescNum;
    if (string(name).substr(0, 19) == "ProblemDescription-") {
      ADF_Get_Node_ID(probDescID, name, &pDescNum, &err);
      checkError(err, probDescID, "Cannot get node ID of " + string(name));

      // get number of ProblemDescription-# children
      int numDesc;
      ADF_Number_of_Children(pDescNum, &numDesc, &err);
      checkError(err, probDescID,
                 "Cannot get number of children of the ProblemDescriptions-#");

      // loop over each of the ProblemDescription-# childeren
      char cname[ADF_NAME_LENGTH_C] = {'\0'};
      int cnum_ret = 0;
      for (int ii = 1; ii <= numDesc; ii++) {
        ADF_Children_Names(pDescNum, ii, 1, ADF_NAME_LENGTH_C, &cnum_ret, cname,
                           &err);
        checkError(err, pDescNum, "Cannot get node ID of " + string(cname));

        // if child is BoundaryRegion-# probe further for boundary condtion data
        double bNum;
        if (string(cname).substr(0, 15) == "BoundaryRegion-") {
          ADF_Get_Node_ID(pDescNum, cname, &bNum, &err);
          checkError(err, pDescNum, "Cannot get node ID of " + string(cname));

          // get boundary condition index and name
          string index = string(cname).substr(15);
          int bcIndex = atoi(index.c_str());
          string bcLabel = readNodestr(bNum, "Label");

          removeSpace(bcLabel);
          checkUnique(surf_ids, bcIndex, bcLabel);

          surf_ids.push_back(pair<int, string>(bcIndex, bcLabel));
        }
      }
    }
  }
}

// for each processor, get its vertices and topology ID, might open another file
void getMeshID(double processorID, double *verticesID, double *topoID) {
  int numChildren;
  int err;

  double root;
  ADF_Number_of_Children(processorID, &numChildren, &err);
  checkError(err, processorID,
             "Cannot get number of children of the first child of States");

  // get vertices
  if (verticesID) {
    string verticesFileName = readNodestr(processorID, "VerticesFile");
    if (verticesFileName.size() > 0) {
      ADF_Database_Open(verticesFileName.c_str(), "READ_ONLY", "NATIVE", &root,
                        &err);
      checkError(err, root, "Error opening file");
    } else {
      ADF_Get_Root_ID(processorID, &root, &err);
      checkError(err, processorID, "Error get root ID");
    }

    int verticesIndex = readNodei32(processorID, "VerticesId");

    char nodeName[ADF_NAME_LENGTH_C] = {'\0'};
    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "/Meshes/Vertices-%d", verticesIndex);

    ADF_Get_Node_ID(root, nodeName, verticesID, &err);
    checkError(err, root, string("Error getting node ID ") + string(nodeName));
  }
  // get topoID
  if (topoID) {
    string topoFileName = readNodestr(processorID, "TopologyFile");
    if (topoFileName.size() > 0) {
      ADF_Database_Open(topoFileName.c_str(), "READ_ONLY", "NATIVE", &root,
                        &err);
      checkError(err, root, "Error opening file");
    } else {
      ADF_Get_Root_ID(processorID, &root, &err);
      checkError(err, processorID, "Error get root ID");
    }

    int topoIndex = readNodei32(processorID, "TopologyId");

    char nodeName[ADF_NAME_LENGTH_C] = {'\0'};
    bzero(nodeName, ADF_NAME_LENGTH_C);
    snprintf(nodeName, ADF_NAME_LENGTH_C, "/Meshes/FaceBasedTopology-%d",
             topoIndex);

    ADF_Get_Node_ID(root, nodeName, topoID, &err);
    checkError(err, root, string("Error getting node ID ") + string(nodeName));
  }
}


void readCoordinatesInChunks(const double &posScale,
                             const vector<double> processorIDs,
                             int64 &totalNumNode,
                             store<vector3d<double> > &pos) {
  Loci::debugout << "readCoordinatesInChunks()" << endl;
  // get rank and number of processors
  const int R = Loci::MPI_rank;
  const int P = Loci::MPI_processes;

  MPI_Bcast(&totalNumNode, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  vector<entitySet> local_nodes;
  local_nodes.resize(P);
  vector<int> localNodesSet = VOG::simplePartitionVec(0, totalNumNode - 1, P);
  for (int ii = 0; ii < P; ++ii) {
    local_nodes[ii] = interval(localNodesSet[ii], localNodesSet[ii + 1] - 1);
  }

  // read a block of nodes on rank 0 and send to corresponding processor
  pos.allocate(local_nodes[R]);
  Loci::debugout << "local nodes for this processor are " << local_nodes[R]
                 << endl;
  if (R == 0) {
    cout << "Reading nodal coordinate data..." << endl;
    size_t pid = 0;
    int64 pos_id = 0;
    int64 n_id = 1;  // ADF node data indices start at 1
    readCoordinatesBlock(processorIDs, pos, pid, pos_id, n_id);
    for (int i = 1; i < P; ++i) {
      store<vector3d<double> > pos_send;
      pos_send.allocate(local_nodes[i]);
      readCoordinatesBlock(processorIDs, pos_send, pid, pos_id, n_id);
      Loci::debugout << "Sending data of size " << 3 * pos_send.domain().size()
                     << " to processor " << i << endl;
      int first = pos_send.domain()[0].first;
      MPI_Send(&(pos_send[first]), 3 * pos_send.domain().size(),
               MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
  } else {  // receive node data from rank 0
    MPI_Status status;
    Loci::debugout << "Receiving data of size " << 3 * pos.domain().size()
                   << " from rank 0 at pos entity id " << pos.domain()[0].first
                   << endl;
    int first = pos.domain()[0].first;
    MPI_Recv(&(pos[first]), 3 * pos.domain().size(), MPI_DOUBLE, 0,
             0, MPI_COMM_WORLD, &status);
  }

  FORALL(pos.domain(), ei) { 
    pos[ei] *= posScale; 
  }ENDFORALL;

  MPI_Barrier(MPI_COMM_WORLD);
  if (R == 0) {
    cout << "  ...done reading nodal coordinate data" << endl;
  }
}

vector<entitySet> determineLocalFaces(const int64 &totalNumFace) {
  const int P = Loci::MPI_processes;
  vector<int> faceSetPartitions =
      VOG::simplePartitionVec(0, totalNumFace - 1, P);
  vector<entitySet> local_faces(P);
  for (int pp = 0; pp < P; ++pp) {
    entitySet local_bnd_faces =
        interval(faceSetPartitions[pp], faceSetPartitions[pp + 1] - 1);
    local_faces[pp] += local_bnd_faces;
  }
  return local_faces;
}

void readCCMInChunks(const double &posScale, const char *argv1,
                     store<vector3d<double> > &pos, Map &cl, Map &cr,
                     multiMap &face2node,
                     vector<pair<int, string> > &surf_ids) {
  // get rank
  const int R = Loci::MPI_rank;
  // initialize variables
  int64 totalNumNode = 0;
  int64 totalNumFace = 0;

  // open file on rank 0, get total number of nodes and faces
  double root;
  int err;
  vector<double> processorIDs;
  // outer vector of size processorIDs
  vector<vector<std::pair<double, int> > > boundaryIDs;
  if (R == 0) {
    // open file
    ADF_Database_Open(argv1, "READ_ONLY", "NATIVE", &root, &err);
    checkError(err, root, "Cannot open " + std::string(argv1) + " for reading");

    // get node states
    double statesID;
    ADF_Get_Node_ID(root, "States", &statesID, &err);
    checkError(err, root, "Cannot get node ID of States");

    // get first state
    int num_ret;
    char name[ADF_NAME_LENGTH_C] = {'\0'};
    double stateID;
    ADF_Children_Names(statesID, 1, 1, ADF_NAME_LENGTH_C, &num_ret, name, &err);
    checkError(err, statesID, "Cannot get children names of States");

    ADF_Get_Node_ID(statesID, name, &stateID, &err);
    checkError(err, stateID,
               "Cannot get the node ID of the first child of States");

    // get all processors
    int numChildren;
    ADF_Number_of_Children(stateID, &numChildren, &err);
    checkError(err, stateID,
               "Cannot get number of children of the first child of States");

    for (int start = 1; start <= numChildren; start++) {
      ADF_Children_Names(stateID, start, 1, ADF_NAME_LENGTH_C, &num_ret, name,
                         &err);
      checkError(err, stateID, "Cannot get node ID of ");

      double processorID;
      std::string nodeName(name);
      if (string(name).substr(0, 10) == "Processor-") {
        ADF_Get_Node_ID(stateID, name, &processorID, &err);
        checkError(err, stateID, string("Cannot get node ID of ") + nodeName);
        processorIDs.push_back(processorID);
      }
    }

    // go through all processors, compute total number of nodes and faces
    Loci::debugout << "ccm mesh has " << processorIDs.size() << " processor-ids"
                   << endl;
    for (size_t i = 0; i < processorIDs.size(); ++i) {
      double topoID = 0;
      double verticesID;
      getMeshID(processorIDs[i], &verticesID, &topoID);
      totalNumNode += getNumNodes(verticesID);
      totalNumFace += getNumFaces(topoID);
      boundaryIDs.push_back(getBoundaryFaceIDs(topoID));
    }
    cout << "ccm mesh has " << totalNumNode << " nodes and " << totalNumFace
         << " faces" << endl;

    // get node ProblemDescriptions
    double probDescID;
    ADF_Get_Node_ID(root, "ProblemDescriptions", &probDescID, &err);
    checkError(err, root, "Cannot get node ID of States");

    // get boundary condition names and indices
    getSurfaceNames(probDescID, surf_ids);
    cout << "ccm mesh has " << surf_ids.size() << " boundaries" << endl;
  }

  MPI_Bcast(&totalNumNode, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&totalNumFace, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  Loci::debugout << "Found " << totalNumNode << " nodes and " << totalNumFace
                 << " faces in ccm file" << endl;

  // read in node data in chunks and send to corresponding processors
  readCoordinatesInChunks(posScale, processorIDs, totalNumNode, pos);

  // determine allocation of local faces
  vector<entitySet> local_faces = determineLocalFaces(totalNumFace);

  // read in face2cell data in chunks and send to corresponding processors
  readFace2CellInChunks(processorIDs, boundaryIDs, local_faces, cl, cr);

  // read in face2node data in chunks and send to corresponding processors
  readFace2NodeInChunks(processorIDs, boundaryIDs, local_faces, face2node);

  // close file on rank 0
  if (R == 0) {
    ADF_Database_Close(root, &err);
  }
}


// facemap is not read, assume no duplicate faces between processors
int main(int argc, char *argv[]) {
  bool optimize = true;
  Loci::Init(&argc, &argv);  // Loci initialize

  string Lref = "NOSCALE";

  while (argc >= 2 && argv[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if (argc >= 3 && !strcmp(argv[1], "-Lref")) {
      Lref = argv[2];
      argc -= 2;
      argv += 2;
    } else if (argc >= 2 && !strcmp(argv[1], "-v")) {
      cout << "Loci version: " << Loci::version() << endl;
      if (argc == 2) {
        Loci::Finalize();
        exit(0);
      }
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-o")) {
      optimize = false;
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-in")) {
      Lref = "1 inch";
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-ft")) {
      Lref = "1 foot";
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-cm")) {
      Lref = "1 centimeter";
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-m")) {
      Lref = "1 meter";
      argc--;
      argv++;
    } else if (argc >= 2 && !strcmp(argv[1], "-mm")) {
      Lref = "1 millimeter";
      argc--;
      argv++;
    } else {
      cerr << "argument " << argv[1] << " is not understood." << endl;
      argc--;
      argv++;
    }
  }

  if (Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units."
         << endl;
    exit(-1);
  }

  if (Lref == "") Lref = "1 meter";

  if (!isdigit(Lref[0])) {
    Lref = string("1") + Lref;
  }

  Loci::UNIT_type tp;
  istringstream iss(Lref);
  iss >> tp;
  double posScale = tp.get_value_in("meter");

  if (argc < 2 || argc > 3) {
    cerr << "Usage: ccm2vog input.ccm[g] [output.vog]" << endl;
    exit(1);
  }

  string infile = string(argv[1]);
  string outfile;
  string case_name;

  // only accept "*.ccm" or "*.ccmg" file as input
  size_t p = infile.rfind('.');
  if (p != string::npos) {
    if (infile.substr(p, infile.size()) != ".ccm" &&
        infile.substr(p, infile.size()) != ".ccmg") {
      cerr << "Usage: ccm2vog input.ccm[g] [output.vog]" << endl;
      exit(1);
    }
    case_name = infile.substr(0, p);
  } else {
    cerr << "Usage: ccm2vog [options] input.ccm[g] [output.vog]" << endl
         << "options:" << endl
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -m  : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl;
    exit(1);
  }

  if (argc == 2)
    outfile = case_name + ".vog";
  else
    outfile = argv[2];

  if (string(argv[1]) == string("-o")) {
    optimize = false;
    argv++;
    argc--;
  }

  // read ccm file
  store<vector3d<double> > pos;
  multiMap face2node;
  Map cl, cr;
  vector<pair<int, string> > surf_ids;
  readCCMInChunks(posScale, argv[1], pos, cl, cr, face2node, surf_ids);

  // Loci takes over
  if (Loci::MPI_rank == 0) {
    cout << "converting to VOG format using " << Loci::MPI_processes
         << " processors" << endl;
  }

  if (Loci::MPI_rank == 0) cout << "coloring matrix" << endl;
  VOG::colorMatrix(pos, cl, cr, face2node);

  if (optimize) {
    if (Loci::MPI_rank == 0) cout << "optimizing mesh layout" << endl;
    VOG::optimizeMesh(pos, cl, cr, face2node);
  }

  if (Loci::MPI_rank == 0) cout << "writing VOG file" << endl;
  Loci::writeVOG(outfile, pos, cl, cr, face2node, surf_ids);

  Loci::Finalize();
  return 0;
}

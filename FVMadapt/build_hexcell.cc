// This file containing the functions that build a HexCell according to Loci data structures
//and refinementplans, when build a hexcell, each edge is defined as local edge in cell and resplit with
//needReverse. Each face is defined as local face in hexcell and resplit with orientcode


#include <vector>
#include <list>
#include <Loci.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "hexcell.h"


using std::cerr;
using std::endl;
using std::swap;
using std::cout;

//this function define the 6 faces of hexcell from loci data structures
Array<Entity, 6> collect_hex_faces( const Entity* lower, const Entity* upper,
                                    const Entity* boundary_map,  const Array<char, 6>& hex2face){
  
  // Collect the entity designation of all faces of the cell cc in the all_faces array
  Array<Entity, 6> faces;
  for (int f=0; f<6; f++) {
    switch (hex2face[f]/6) {
    case 0:
      faces[f] = lower[hex2face[f]%6];
      break;
    case 1:
      faces[f] = upper[hex2face[f]%6];
      break;
    case 2:
      faces[f] = boundary_map[hex2face[f]%6];
      break;
    default:
      cerr << " WARNING: illegal hex2face value" << endl;
      break;
    } 
  }
  return faces;
}


//this function define the 8 vertices of hexcell from loci data structures
Array<Entity, 8> collect_hex_vertices(const const_multiMap& face2node, const Array<Entity, 6>& faces,
                                  const Array<char, 8>& hex2node){
  Array<Entity, 8>  hex_vertices;
  for(int i = 0; i < 4; i++){
    hex_vertices[i] = face2node[faces[DOWN]][hex2node[i]];
  }
  for(int i = 4; i < 8; i++){
    hex_vertices[i] = face2node[faces[UP]][hex2node[i]];
  }
  
  std::swap(hex_vertices[1], hex_vertices[4]);
  std::swap(hex_vertices[3], hex_vertices[6]);
  return hex_vertices;
}




//collect the entity designation of all edges of the cell , and the direction of edges are stored in needReverse
Array<Entity, 12> collect_hex_edges(const Array<Entity, 6>& faces, const Array<Entity,8>& hex_vertices,
                                const const_multiMap& face2edge, const const_MapVec<2>& edge2node,
                                Array<bool, 12>& needReverse){
  entitySet all_edges;
  Array<Entity, 12> hex_edges;
  for(int fID = 0; fID < 6; fID++){
    for(int eID = 0; eID < 4; eID++){
      all_edges += face2edge[faces[fID]][eID];
    }
  }      
  if(all_edges.size() != 12){
    cerr << "WARNING: the number of edges is not 12" << endl;
    Loci:: Abort();
  }
  Entity head, tail;
  int node0[12] = {0, 1, 2, 3, 0, 1, 4, 5, 0, 2, 4, 6};
  int node1[12] = {4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7};
  
  for( int e = 0; e < 12; e++ ){
    head = hex_vertices[node0[e]];
    tail =  hex_vertices[node1[e]];
    entitySet::const_iterator ee;
    for(ee = all_edges.begin(); ee != all_edges.end(); ee++){
      if(edge2node[*ee][0] == head && edge2node[*ee][1] == tail){
        hex_edges[e] = *ee;
        needReverse[e] = false;
        break;
      }
      else if(edge2node[*ee][1] == head && edge2node[*ee][0] == tail){
        hex_edges[e] = *ee;
        needReverse[e] = true;
        break;
      }
      
    }
    if(ee ==  all_edges.end()){
      cerr << "WARNING: can not find edge" << endl;
      Loci::Abort();
    }
  }
  
  return hex_edges;
}






  
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,
                                               upper,
                                               boundary_map,
                                               hex2face);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                  face_entity,
                                                  hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                 node_entity,
                                                 face2edge,
                                                 edge2node,
                                                 edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;


    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}

// //serial version
// HexCell* build_hex_cell(const Entity* lower, int lower_size,
//                         const Entity* upper, int upper_size,
//                         const Entity* boundary_map, int boundary_map_size,
//                         const Array<char,6>& hex2face,
//                         const Array<char,8>& hex2node,
//                         const Array<char,6>& orientCode,
//                         const const_multiMap& face2node,
//                         const const_multiMap& face2edge,
//                         const const_MapVec<2>& edge2node,
//                         const const_store<vect3d>& pos,
//                         const const_store<std::vector<char> >& edgePlan,
//                         const const_store<std::vector<char> >& facePlan,
//                         const store<int>& node_offset,
//                         int offset_min,
//                         std::list<Node*>& bnode_list,
//                         std::list<Edge*>& edge_list,
//                         std::list<QuadFace*>& face_list){

//   Array<Entity, 6> face_entity = collect_hex_faces(lower,
//                                                upper,
//                                                boundary_map,
//                                                hex2face);
  
//   Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
//                                                   face_entity,
//                                                   hex2node);



//   Array<bool, 12> edge_reverse;
//   Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
//                                                  node_entity,
//                                                  face2edge,
//                                                  edge2node,
//                                                  edge_reverse);

 
//   //define each node and put it into node_list
//   std::map<Entity, Node*> n2n;
//   for(int i = 0; i < 8; i++){
//     Node* aNode = new Node(pos[node_entity[i]], node_entity[i]-offset_min+1);
//     bnode_list.push_back(aNode);
//     n2n[node_entity[i]] = aNode;
//   }
  
  
//   std::map<Entity, Edge*> e2e;
//   std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
//   for(int i = 0; i < 12; i++){
//     Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
//                             n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
//     edge_list.push_back(anEdge);
//     e2e[edge_entity[i]] = anEdge;


//     //resplit the edge
//     anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);

//      //index the node
//     int nindex = node_offset[edge_entity[i]];
//     for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
//       (*np)->index =  nindex++;
//     }
//     bnode_begin = --(bnode_list.end());
//   }
  
//   int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
//                   {1, 7, 3, 5}, {0, 6, 2, 4}};
  
//   //defines each face and put it into face_list
//   QuadFace** face = new QuadFace*[6];
//   for(int i  = 0; i < 6; i++){
//     face[i] = new QuadFace(4);
//     face_list.push_back(face[i]);
//     //define each edge
//     for(int j = 0; j < 4; j++){
//       face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
//     }
//     //resplit each face
//     face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);

//      //index the node
//     int nindex = node_offset[face_entity[i]];
//     for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
//       (*np)->index =  nindex++;
//     }
//     bnode_begin = --(bnode_list.end());
//   }
  
//   HexCell* aCell = new HexCell(face);
//   return aCell;
// }

HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                         const const_store<char>& posTag,
                        const const_store<std::vector<char> >& nodeTag,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,
                                               upper,
                                               boundary_map,
                                               hex2face);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                  face_entity,
                                                  hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                 node_entity,
                                                 face2edge,
                                                 edge2node,
                                                 edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;


    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag =  nodeTag[edge_entity[i]][nindex++];
      // checkNode(*np);
    }
    bnode_begin = --(bnode_list.end()); 
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);

    int   nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag = nodeTag[face_entity[i]][nindex++];
      //  checkNode(*np);
    } 
    
    bnode_begin= --(bnode_list.end());
    
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}

//parallel version
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        const const_store<int>& node_offset,
                        const Map& node_l2f,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,
                                               upper,
                                               boundary_map,
                                               hex2face);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                  face_entity,
                                                  hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                 node_entity,
                                                 face2edge,
                                                 edge2node,
                                                 edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]], node_l2f[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;


    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);

     //index the node
    int nindex = node_offset[edge_entity[i]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);

     //index the node
    int nindex = node_offset[face_entity[i]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}




// //build without restart(no edgeplan or faceplan)
// HexCell* build_hex_cell(const Entity* lower, int lower_size,
//                         const Entity* upper, int upper_size,
//                         const Entity* boundary_map, int boundary_map_size,
//                         const Array<char,6>& hex2face,
//                         const Array<char,8>& hex2node,
//                         const Array<char,6>& orientCode,
//                         const const_multiMap& face2node,
//                         const const_multiMap& face2edge,
//                         const const_MapVec<2>& edge2node,
//                         const const_store<vect3d>& pos,
//                         std::list<Node*>& bnode_list,
//                         std::list<Edge*>& edge_list,
//                         std::list<QuadFace*>& face_list){

//   Array<Entity, 6> face_entity = collect_hex_faces(lower,
//                                                upper,
//                                                boundary_map,
//                                                hex2face);
  
//   Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
//                                                   face_entity,
//                                                   hex2node);



//   Array<bool, 12> edge_reverse;
//   Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
//                                                  node_entity,
//                                                  face2edge,
//                                                  edge2node,
//                                                  edge_reverse);

 
//   //define each node and put it into node_list
//   std::map<Entity, Node*> n2n;
//   for(int i = 0; i < 8; i++){
//     Node* aNode = new Node(pos[node_entity[i]]);
//     bnode_list.push_back(aNode);
//     n2n[node_entity[i]] = aNode;
//   }
  
  
//   std::map<Entity, Edge*> e2e;
//   for(int i = 0; i < 12; i++){
//     Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
//                             n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
//     edge_list.push_back(anEdge);
//     e2e[edge_entity[i]] = anEdge;
    
//   }
  
//   int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
//                   {1, 7, 3, 5}, {0, 6, 2, 4}};
  
//   //defines each face and put it into face_list
//   QuadFace** face = new QuadFace*[6];
//   for(int i  = 0; i < 6; i++){
//     face[i] = new QuadFace(4);
//     face_list.push_back(face[i]);
//     //define each edge
//     for(int j = 0; j < 4; j++){
//       face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
//     }
//   }
  
//   HexCell* aCell = new HexCell(face);
//   return aCell;
// }
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                          const const_store<char>& posTag,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,
                                               upper,
                                               boundary_map,
                                               hex2face);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                  face_entity,
                                                  hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                 node_entity,
                                                 face2edge,
                                                 edge2node,
                                                 edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}

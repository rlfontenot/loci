///////////////////////////////////////////////////////////////////////////////////////////////
//                          quadface.h
//
// This function includes the declaration of class QuadFace.
// QuadFace is used for anisotropic refinement of hexcell and Prism
// Given a face with 4 edges, the direction of edges are:
//                3 -> 2
//                ^    ^
//                |    |
//                0 -> 1
// after split in both x and y directions,  the childID is :
//
//                1    3
//
//                0    2
//



#ifndef QUADFACE_H
#define QUADFACE_H
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <bitset>
#include "node_edge.h"
using std::ofstream;
using std::bitset;

struct Range2d;

//f2c orient functions are uesd when a quadface is built as in cell
// and the facePlan is for the face defined by face2node 
char orient_splitCode_f2c(char splitCode, char orientCode);
char orient_childID_f2c(char childID, char orientCode, char splitCode);
char orient_edgeID_f2c(char edgeID, char orientCode);

//c2f orient functions are used when the face is defined by face2node
//and facePlan is extracted from cell 
char orient_edgeID_c2f(char edgeID, char orientCode);



//for get_c1_hex and get_c1_prism, when a cell is empty_split, faceMap is created
// and face is split to generate leaves, by look-up the faceMap, c1 can be computed
//use Range2d to avoid full split cell and face.  
std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves);


//QuadFace,
class QuadFace{
public:
  //constructors 
  QuadFace(int numEdge):edge(new Edge*[numEdge]),child(0),childx(0), childy(0),code(char(0)){}
  //used for empty_split
  QuadFace():edge(0), child(0),childx(0), childy(0),code(char(0)){}
  QuadFace( Edge** e):edge(e),child(0), childx(0), childy(0),code(char(0)){}
 

  //destructor, it works this way without memory leakage
  ~QuadFace(){
    if(this != 0){
      switch(code){
      case 3:
       
        if(childx != 0){
          for(int i = 0; i < 2; i++){
            if(childx[i] != 0){
              delete childx[i];
              childx[i] = 0;
            }
          }
          delete[] childx;
          childx = 0;
        }
        if(childy != 0){
          for(int i = 0; i < 2; i++){
            if(childy[i] != 0) {
              delete[] childy[i]->childx;
              childy[i]->childx = 0;
              delete childy[i];
              childy[i] = 0;
            }
          }
          delete[] childy;
          childy = 0;
        }


        if(child!= 0){
          delete[] child;
          child = 0;
        }
        break;
      case 2:
        if(childx != 0){
          for(int i = 0; i < 2; i++){
            if(childx[i] != 0){
              delete childx[i];
               childx[i] = 0;
            }
          }
          delete[] childx;
          childx = 0;
        }
        break;
      case 1:
        if(childy != 0){
          for(int i = 0; i < 2; i++){
            if(childy[i] != 0){
              delete childy[i];
              childy[i] = 0;
            }
          }
          delete[] childy;
          childy = 0;
        }
        break;
      default:
        break;
      }
      if(edge != 0){
        delete [] edge;
        edge = 0;
      }
    }
  }

 
  
  
  //the center of the face, defined as the mass center of 4 nodes
  inline Node* QuadFace::centroid(){
    vect3d nodes[4];
    //get nodes
    for(int i = 0; i <2; i++){
      nodes[i] = edge[i]->head->p;
    }
    for(int i = 2; i <4; i++){
      nodes[i] = edge[i]->tail->p;
    }
    
    //calculate the mass center of nodes
    vect3d p = point_center(nodes, 4);
    return new Node(p);
  }
  

  //precondition: all its edges have been split
  //condition: edgecenter must be allocated and deallocated by caller
  inline void getEdgeCenter(Node** edgecenter)const{
    for(int i = 0; i <4; i++){
      edgecenter[i] = edge[i]->child[0]->tail;
    }
  }

  inline Node* getEdgeCenter(int edgeID)const{
    return edge[edgeID]->child[0]->tail;
  }
  
  
 
  
  inline Node* getNode(int nodeID)const{
    if(nodeID==0 || nodeID == 1){
      return edge[nodeID]->head;
    }
    else {
      return edge[nodeID]->tail;
    }
  }
  //only if code is 3
  inline Node* getCenter()const{
    return child[0]->edge[1]->tail; //unsafe version 
  }
  
  //get all the leaves of this, suppose this has been resplit
  void get_leaves(std::vector<QuadFace*>& leaves);

  // this has been built as defined in cell, and has been resplit with orientCode
  //get leaves that is in the same order as the face is bulit as defined by face2node and resplit without orientCode
  void get_leaves(const std::vector<char>& facePlan, char orientCode,std::vector<QuadFace*>& fine_faces); 

 //define face2node
  void set_f2n(std::list<int32>& f2n);

  int get_num_leaves()const;

  //used in building cells, quadface is built as defined in cell, and split with orientCode
  //all new nodes and edges are put into node_list and edge_list
 
  void split(char splitCode, char orientCode,
             std::list<Node*>& node_list,
             std::list<Edge*>& edge_list);

  void empty_split(char splitCode);

           
  //used in building cells, quadface is built as defined in cell, and split with orientCode
  //all new nodes and edges are put into node_list and edge_list 
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list);

  //used in building cells, quadface is built as defined in cell, and split with orientCode
  //all new nodes and edges are put into node_list and edge_list 
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<QuadFace*>& fine_faces);

  void empty_resplit(const std::vector<char>& facePlan, char orientCode,
                     std::vector<QuadFace*>& fine_faces);
  
public:
  
  Edge** edge;
  QuadFace** child; 
  QuadFace** childx; 
  QuadFace** childy;
  //split code, can change value during splitting
  //if code is 1, split only in y direction, only childy is defined. childx = child = 0
  //if code is 2, split only in x direction, 2 childx, childy= child = 0;
  //if coode is 3, split in both x and y direction, 2 childx, 2 childy, 4 child
  char code; 
};


QuadFace* build_quad_face( const Entity* face2node, 
                           const Entity* face2edge,
                           const const_multiMap& edge2node,
                           const const_store<vect3d>& pos,
                           const const_store<std::vector<char> >& edgePlan,
                           std::list<Node*>& bnode_list,
                           std::list<Edge*>& edge_list);


//parallel version
QuadFace* build_quad_face( const Entity* face2node, 
                           const Entity* face2edge,
                           const const_multiMap& edge2node,
                           const const_store<vect3d>& pos,
                           const const_store<std::vector<char> >& edgePlan,
                           const const_store<int>& node_offset,
                           const const_store<int>& node_l2f,
                           std::list<Node*>& bnode_list,
                           std::list<Edge*>& edge_list);

//this function is used in build_general_cell with quadface
QuadFace* build_tmp_quad_face( const Entity* face2node, 
                               const Entity* face2edge,
                               const const_multiMap& edge2node,
                               const const_store<std::vector<char> >& edgePlan,
                               std::list<Node*>& bnode_list,
                               std::list<Edge*>& edge_list);




//if the intersection of the leaves of f1 and the leaves of f2 is empty 
bool is_overlapped( QuadFace* f1,   QuadFace* f2);

//return the intersection of leaves of f1 and leaves of f2
std::vector<QuadFace*> overlap( QuadFace* f1,   QuadFace* f2);

//for serial version, write out .cog file
void  write_quad_inner_faces(const std::map<QuadFace*, NeibIndex>& faces,
                             int cell_offset, int& mxppf, ofstream& ofile);



//compile the facePlan according the tree structure of aQuadFace
//std::vector<char> make_faceplan( QuadFace* aFace);
inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list,
                         std::list<QuadFace*>& face_list){
  
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
  
   for(std::list<QuadFace*>::iterator p = face_list.begin();  p != face_list.end(); p++){
    if((*p) != 0){
    delete (*p);
    (*p) = 0;
    }
  }
  face_list.clear();
  
}

inline void cleanup_list(std::list<QuadFace*>& face_list){
  for(std::list<QuadFace*>::iterator p = face_list.begin();  p != face_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  face_list.clear();
  
}




void extract_quad_edge(const std::vector<char>&, std::vector<char>&, unsigned int);
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL);
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL,
                                  std::vector<char>& facePlanR, char orientCodeR);
#endif


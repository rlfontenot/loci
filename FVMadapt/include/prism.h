/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
//                               
//
//   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRISM_H
#define PRISM_H
#include <Loci.h>
#include <vector>
#include <bitset>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <utility> 
#include "hex_defines.h"
#include "quadface.h"
#include "face.h"


using std::cerr;
using std::endl;



class Prism
{
public:
  
  //constructor; 

  //this construct is used when gnrlface and quadface don't need to be defined.
  //it doesn't set nfold, usually after this constructor, setNfold() will be called
  Prism():nfold(0),mySplitCode(0),gnrlface(0),quadface(0), parentCell(0),
           childCell(0){faceOrient.reset();}
  
  Prism(int n):nfold(n),mySplitCode(0),gnrlface(new Face*[2]), quadface(new QuadFace*[n]), parentCell(0), childCell(0){
    faceOrient.reset();
  }
  
  //destructor
 ~Prism(){
    if(this != 0 ){
      if(childCell != 0){
        for(int i = 0; i < numChildren(); i++){
          if(childCell[i] != 0){
            delete  childCell[i];
            childCell[i] = 0;
          }
        }
        delete[] childCell;
        childCell = 0;
      }
      parentCell = 0;
    }
    if(gnrlface != 0){
      delete[] gnrlface;
      gnrlface = 0;
    }

    if(quadface != 0){
      delete[] quadface;
      quadface = 0;
    }
 }
  
  
  int32 getCellIndex() const {return cellIndex;}
  
  inline int getLevel(int d){
    switch(d){
    case 0: //xy direction
      return gnrlface[0]->edge[0]->level;
      break;
    case 1: //z direction
      return quadface[0]->edge[1]->level;
      break;
    default:
      cerr << "WARNING: illegal levelID" << endl;
      break;
    }
    return 0;
  }
  
  
  inline void setNfold(int d){
    nfold = d;
  }

  inline int getNfold() const{
    return nfold;
  }
  inline char getMySplitCode() const{
    return mySplitCode;
  }
  inline Prism* getChildCell(int i) const{
    return childCell[i];
  }
  
  inline int numChildren(){
    switch(mySplitCode){
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return nfold;
    case 3:
      return 2*nfold;
    default:
      cerr<< "WARNING: illegal split code" << endl;
      break;
    }
    return -1;
  }
  
  //used in build_prismcell.cc
  inline void setFace(int faceID, Face* aFace){
    gnrlface[faceID] = aFace;
  }
  
  inline void setFace(int faceID, QuadFace* aFace){
    quadface[faceID] = aFace;
  }
  
  double get_min_edge_length();
  //define if this is tagged for refinement
  bool get_tagged();
  int get_num_fine_faces();//for mxfpc
  //  void setParentCell( Prism* parent){parentCell = parent;};
  
  inline int whichChild(){
    if(parentCell == 0) return -1;
    for(int i = 0; i < parentCell->numChildren(); i++){
      if(this == parentCell->childCell[i]) return i;
    }
    return -1;
  }
  
    
  //this function splits diamondcell isotropially once
  // newly created nodes, edges and faces are put into the lists
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<QuadFace*>& quadface_list,
             std::list<Face*>& face_list);
  
  //this function splits diamondcell isotropically once
  //only define childCell
  //  void empty_split();
  // void empty_resplit(const std::vector<char>& cellPlan);

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& quadface_list,
               std::list<Face*>& face_list,
               std::vector<Prism*>& prism_cells);

  //used in make_prism_cellplan.cc
  void Prism::resplit( int level,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& quadface_list,
                       std::list<Face*>& face_list);
 
  void empty_split();
  void empty_resplit(const std::vector<char>& cellPlan);
  
  
  //this function check if aFace is one of the face
  //return -1: No
  //return i in [0, 6), Yes, face[i] == aFace
  //  inline int containFace(QuadFace* aFace){
  //  for(int i = 0; i <6; i++){
  //    if(face[i] == aFace) return i;
  //  }
  // return -1;
  // };
  
  //this function check if aCell is my sibling neighbor.
  //sibling means same size face, not necessarily has same parent
  bool  isSiblingNeighbor(const Prism* aCell, int dd, int &nf)const;

  //get real sibling neib
  Prism* Prism::getSiblingNeib(int dd, int& nf);

  //only use it when dd >=2
  inline int parentFace(int dd){
   if(parentCell->mySplitCode == 1) return dd;
    
    else if(dd == 2){
      return (whichChild()%(parentCell->nfold))+2;
    }
   else if(dd == 5){
     int childID = (whichChild()%parentCell->nfold);
     
     return (childID== 0)?
       (parentCell->nfold +1):(childID+1);
   }
   return -1;
  }
  
  //This function check if aCell is my neighbor is direction dd
  // edge neighbor is excluded. if two faces overlap area is not zero
  // their cells are neighbors
  bool  isNeighbor(const Prism* aCell, int dd, int& nf)const;
  
  //this function find  my face neighbor,
  //ff: in,  my faceID,
  //nf: out, the faceID of neibCell 
  Prism*  findNeighbor(int d, int& nf);

  bool getTagged();
  void setSplitCode(int split_mode);
 //after a cell is split, compose the cell plan according to the tree structure 
  std::vector<char> make_cellplan();
  
  //make a plan for isotropical split an original prism cell level levels 
  std::vector<char> make_cellplan(int level);
  
  inline std::vector<Edge*> get_edges(){
  
    std::vector<Edge*> edges(3*nfold);
    for(int i = 0; i < nfold; i++) edges[i] = gnrlface[0]->edge[i];
    for(int i = nfold; i< 2*nfold; i++) edges[i] = gnrlface[1]->edge[i-nfold];
    for(int i = 2*nfold; i < 3*nfold; i++){
      int j = i%nfold;
      edges[i] = quadface[j]->edge[faceOrient.test(j)?1:3];
    }

    return edges;
  }
  //if any edge is more than 1 levels down than my level, split myself, then balance each child 
  bool balance_cell(int split_mode,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<QuadFace*>& qface_list,
                    std::list<Face*>& gface_list);
  
  void rebalance_cells(int split_mode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& qface_list,
                       std::list<Face*>& gface_list);
  
  friend void set_prism_faces(const std::vector<Prism*>& cells,
                              std::map<QuadFace*, NeibIndex>& quadfaces,
                              std::map<Face*, NeibIndex>& faces);
  
  friend std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                         const std::vector<char>& facePlan,
                                         char orientCode,
                                         int faceID);
  friend std::vector<char>  extract_prism_face(const  std::vector<char>& cellPlan, int dd);
  
 void print();
private:
 
  int32 cellIndex; //the index of the cell, start with 1
  char nfold; //3 for normal prism, 4 for the children of prism when quadface is split
  
  //describe how the cell will be splitted. the order of dimensions is xyz, and 0 means no split;
  // 1 means split; for example, 011 means  y and z coordinations are splitted.
  char mySplitCode;
  Face** gnrlface;
  QuadFace** quadface;
  //orientCode is not necessary here
  
  Prism *parentCell; // the parent of the cell
  Prism **childCell;// an dynamic array of pointers to children cells

  
  
  std::bitset<6> faceMarked; //if the face in direction RIGHT, LEFT... has been checked
  //when a Prism is first created, gnrlface[0] points inward, gnrlface[1] points outward
  // thress quadfaces point outward. When it's split, gnrlface orientation won't change
  // quadface orientation of the children need to be defined because when a new face is created,
  //it can  not make both of two children sharing it keep their face orientation
  // so when nfold==4, if quadface[i] points inward, faceOrient[i] = 1; else faceOrient[i] = 0;
  std::bitset<4> faceOrient;
  
  //  char whichChild;
  //assignment and copying are prohibited
  void operator=(const Prism&);
  Prism(const Prism&);
  
private:
  //get all the leaves
  void get_leaves(std::vector<Prism*>& leaf_cell);

  //get 6 nodes
  inline void get_nodes(Node** node){
   
    for(int i = 0; i < nfold; i++){
      node[i] = (gnrlface[0]->needReverse[i])?(gnrlface[0]->edge[i]->tail):gnrlface[0]->edge[i]->head;
      node[i+nfold] = (gnrlface[1]->needReverse[i])?(gnrlface[1]->edge[i]->tail):gnrlface[1]->edge[i]->head;
    }
  }
  
  //get all the 4*3 edges
  inline void get_edges(Edge** edge){}
  
  
  //calculate the centroid of the Prism, it's defined as
  //the mean value of nodes
  
  inline Node* centroid(){
   
    Node* cellcenter = new Node();
    Node** vertices = new Node*[2*nfold];
    get_nodes(vertices);
    vect3d* nodes = new vect3d[2*nfold];
    for(int i = 0; i<2*nfold; i++){
      nodes[i] = vertices[i]->p;
    }
    cellcenter->p = point_center(nodes, 2*nfold);
    delete[] vertices;
    delete[] nodes;
    return cellcenter;
  }
  
  //get the facecenter
  //condition: facecenter will be allocated and deallocated by the caller
  //precondition: the faces have been splitted
  // inline void getFaceCenter(Node** facecenter){
 
  //  facecenter[0] = gnrlface[0]->child[0]->edge[2]->head;
  //  facecenter[1] = gnrlface[1]->child[0]->edge[2]->head;
  //  for(int i = 0; i <; i++){
  //   facecenter[i+2] = quadface[i]->child[0]->edge[2]->head;
  // }
  // }

  inline Node* getFaceCenter(int faceID){
    if(faceID < 2) return  gnrlface[faceID]->child[0]->edge[2]->head;
    return quadface[faceID-2]->getCenter();
  }

 
};

Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list);
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
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
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list);
// for no restart
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                         const const_store<char>& posTag,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list);


//serial version
// Prism* build_prism_cell(const Entity* lower, int lower_size,
//                         const Entity* upper, int upper_size,
//                         const Entity* boundary_map, int boundary_map_size,
//                         const Array<char,5>& prism2face,
//                         const Array<char,6>& prism2node,
//                         const Array<char,5>& orientCode,
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
//                         std::list<QuadFace*>& qface_list,
//                         std::list<Face*>& gface_list);

//parallel version
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
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
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list);



int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge);
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge);
#endif
  

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
////////////////////////////////////////////////////////////////////////////////////////
//                                get_fine_grid.cc                                    //
//                                by: Qiuhan Xue                                      //
//                                                                                    //
// This file generate inner_nodes  according refinement plans. The
//original grid is general element. Isotropic refinement is used. the rules should be
//usedAfter node_offset and cell_offset has been set up.
// FineFaces is defined as c1, c2,followed by face2node.
//
////////////////////////////////////////////////////////////////////////////////////////                             


#include <iostream>
#include <cstdlib>
#include <vector>
#include <Loci.h>
#include "diamondcell.h"

using std::cerr;
using std::endl;
using std::vector;
//get inner_nodes of general cells
class get_general_cell_nodes : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_multiMap edge2node;
  const_store<vect3d> pos;


  const_store<bool> is_quadface;
  store<Loci::FineNodes> inner_nodes;
  
  const_store<int> node_l2f;
public:
  get_general_cell_nodes(){
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);


    name_store("inner_nodes", inner_nodes);
    name_store("fileNumber(pos)", node_l2f);
    name_store("is_quadface", is_quadface);
    
    input("cellPlan");
    input("(lower, upper, boundary_map)->(facePlan,is_quadface)");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");
    input("(lower, upper, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
   
    output("inner_nodes");
     constraint("gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
     
      do_loop(seq, this);
    }
  }
  void calculate(Entity cc){
    
    // if(cell_offset[cc] > 44671 && cell_offset[cc] < 44735 ) cerr << "gnrl cell : " << cc << " offset: " << cell_offset[cc] << endl;
    
    if(cellPlan[cc].size() == 0){
      vector<vect3d>().swap(inner_nodes[cc]);
      return;
    }
    
    std::list<Node*> bnode_list; //boundary node
    std::list<Edge*> edge_list;
    std::list<Face*> face_list;
    std::list<Node*> node_list; //inner node


    //build a cell
    Cell* aCell =  build_general_cell(lower[cc].begin(), lower.num_elems(cc),
                                      upper[cc].begin(), upper.num_elems(cc),
                                      boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                      is_quadface,
                                      face2node,
                                      face2edge,
                                      edge2node,
                                      pos,
                                      edgePlan,
                                      facePlan,
                                      bnode_list,
                                      edge_list,
                                      face_list,
                                      node_l2f);
      
          
          
          
   
    std::vector<DiamondCell*> cells;
    //split the cell
    aCell->resplit( cellPlan[cc], 
                    node_list,
                    edge_list,
                    face_list,
                    cells);
      
    //put the node_list into inner_nodes, and index the nodes
    vector<vect3d>(node_list.size()).swap(inner_nodes[cc]); 
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++, nindex++){
      inner_nodes[cc][nindex]=(*np)->p;
    }
  
    
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    
    //aCell will clean up these
    cleanup_list(bnode_list, edge_list, face_list);
    cleanup_list(node_list);
  }
};
register_rule<get_general_cell_nodes> register_get_general_cell_nodes;  

//get inner nodes of edges
class get_edge_nodes : public pointwise_rule{
  const_store<std::vector<char> > edgePlan;
  const_store<vect3d> pos;
  const_multiMap edge2node;
  store<Loci::FineNodes> inner_nodes;
public:
  get_edge_nodes(){
    name_store("edgePlan", edgePlan);
    name_store("pos", pos); 
    name_store("edge2node", edge2node);
    name_store("inner_nodes", inner_nodes);
    input("edgePlan");
    input("edge2node->pos");
    output("inner_nodes");
   
  }
  virtual void compute(const sequence &seq){
        do_loop(seq, this);
  }
  void calculate(Entity e){
    if(edgePlan[e].size() != 0){
      std::list<Node*> node_list;
      
      Node* head = new Node(pos[edge2node[e][0]]);
      Node* tail = new Node(pos[edge2node[e][1]]);
      Edge* theEdge = new Edge(head, tail);
        
      theEdge->resplit(edgePlan[e], node_list);
      
      vector<vect3d>(node_list.size()).swap(inner_nodes[e]);
      int nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++, nindex++){
        inner_nodes[e][nindex] = (*np)->p;
      } 
      
      //clean up
      if(head != 0 ) {
        delete head;
        head = 0;
      }
      if(tail != 0){
        delete tail;
        tail = 0;
      }
      if(theEdge !=0){
        delete theEdge;
        theEdge = 0;
      }
      cleanup_list(node_list);
    }
    else{
      vector<vect3d>().swap(inner_nodes[e]);
      return;
    }
  }
};
register_rule<get_edge_nodes> register_get_edge_nodes;  

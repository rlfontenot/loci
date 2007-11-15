////////////////////////////////////////////////////////////////////////////////////////
//                                get_hex_fine_grid.cc                                    //
//                                by: Qiuhan Xue                                      //
//                                                                                    //
// This file generate inner_nodes  according refinement plans. The
//original grid is hexcell element. Anisotropic refinement is used. the rules should be
//used After node_offset and cell_offset has been set up.
// FineFaces is defined as c1, c2,followed by face2node.
//
////////////////////////////////////////////////////////////////////////////////////////                             


#include <iostream>
#include <cstdlib>
#include <vector>
#include <Loci.h>
#include "hexcell.h"

using std::cerr;
using std::endl;
using std::vector;
//get inner_nodes  of hex cells
class get_hex_cell_nodes : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
   const_store<Array<char,6> > orientCode;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  
 const_blackbox<Loci::storeRepP> node_remap;
  store<Loci::FineNodes> inner_nodes;
 


public:
  get_hex_cell_nodes(){
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("hexOrientCode", orientCode);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
    
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
  

    name_store("inner_nodes", inner_nodes);
    name_store("node_remap", node_remap);

    
    input("cellPlan");
    input("(hex2face, hex2node, hexOrientCode)");
    input("(lower, upper, boundary_map)->(facePlan)");
    input("(lower, upper, boundary_map)->face2edge->(edgePlan)");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
 input("node_remap");
    output("inner_nodes");
       constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
   
   
      do_loop(seq, this);
   
   
  }
  void calculate(Entity cc){
    
   
    
    if(cellPlan[cc].size() == 0){
      vector<vect3d>().swap(inner_nodes[cc]);
    
      return;
    }
    
    std::list<Node*> bnode_list; //boundary node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> face_list;
    std::list<Node*> node_list; //inner node
    

   //build a Cell
      HexCell* aCell =  build_hex_cell(lower[cc].begin(), lower.num_elems(cc),
                                       upper[cc].begin(), upper.num_elems(cc),
                                       boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                       hex2face[cc],
                                       hex2node[cc],
                                       orientCode[cc],
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       pos,
                                       edgePlan,
                                       facePlan,
                                       bnode_list,
                                       edge_list,
                                       face_list);
          
      std::vector<HexCell*> cells;
      
      aCell->resplit(cellPlan[cc], 
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
register_rule<get_hex_cell_nodes> register_get_hex_cell_nodes;  



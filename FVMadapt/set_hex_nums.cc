#include <iostream>
#include <cstdlib>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <Loci.h>
#include "hexcell.h"

using std::queue;
using std::cerr;
using std::cout;
using std::endl;
using std::set;
using std::map;

class set_hexcell_nums : public pointwise_rule{
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
  store<int> num_inner_nodes;
  store<int> num_fine_cells;
public:
  set_hexcell_nums(){
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
   
    name_store("num_inner_nodes", num_inner_nodes);
    name_store("num_fine_cells", num_fine_cells);

    input("cellPlan");
    input("(hex2face, hex2node, hexOrientCode)");
    input("(lower, upper, boundary_map)->facePlan");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node");
    input("(lower, upper, boundary_map)->face2edge");
    
    output("num_inner_nodes");
    output("num_fine_cells");
    constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
       
    if(cellPlan[cc].size() == 0){
      num_inner_nodes[cc] = 0;
      num_fine_cells[cc] = 1;
      return;
    }

    std::list<Node*> bnode_list; //boundary node
    std::list<Node*> node_list; //inner node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> face_list;
    
    HexCell* aCell = build_hex_cell(lower[cc].begin(), lower.num_elems(cc),
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
    aCell->resplit( cellPlan[cc],
                    node_list,
                    edge_list,
                    face_list,
                    cells);
  
  
  
    num_inner_nodes[cc] = node_list.size();
    num_fine_cells[cc] = cells.size();
    
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
register_rule<set_hexcell_nums> register_set_hexcell_nums;  



  
class set_quadface_nums : public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  store<int> num_inner_nodes;
  
public:
  set_quadface_nums(){
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("priority::quadface::num_inner_nodes", num_inner_nodes);
    input("facePlan");
    input("face2edge->(edgePlan, edge2node)");
    input("face2node->pos");
    output("priority::quadface::num_inner_nodes");
    constraint("quadfaces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    if(facePlan[f].size() == 0){
      num_inner_nodes[f] = 0;
      return;
    }
    std::list<Node*> bnode_list;
    std::list<Node*> node_list;
    std::list<Edge*> edge_list;
   
    
    QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                      face2edge[f].begin(),
                                      edge2node,
                                      pos,
                                      edgePlan,
                                      bnode_list,
                                      edge_list);
           
    aFace->resplit(facePlan[f],
                   char(0),
                   node_list,
                   edge_list);
    
    
    num_inner_nodes[f] = node_list.size();
    
    //clean up
    if(aFace != 0){
      delete aFace;
      aFace = 0;
    }
    
    cleanup_list(bnode_list, edge_list);
    cleanup_list(node_list);
    
  }
};
register_rule<set_quadface_nums> register_set_quadface_nums;






// class set_hexedge_nums : public pointwise_rule{
//   const_store<std::vector<char> > edgePlan;

//   store<int> num_inner_nodes;
  
// public:
//   set_hexedge_nums(){
//     name_store("edgePlan", edgePlan);
//     name_store("num_inner_nodes", num_inner_nodes);

//     input("edgePlan");


//     output("num_inner_nodes");
//     constraint("edge2node->pos");
//   }
//   virtual void compute(const sequence &seq){
//     do_loop(seq, this);
//     //  cerr << currentMem() << " After set_hexedge_nums " << Loci::MPI_rank << endl; 
//   }
//   void calculate(Entity e){

//     int num = 0;
//     if(edgePlan[e].size()!= 0){
//       for(unsigned int i = 0; i < edgePlan[e].size(); i++){  
//         if(edgePlan[e][i] == 1)num++;
//       }
//     }
//     num_inner_nodes[e] = num;
//   }
// };
// register_rule<set_hexedge_nums> register_set_hexedge_nums;  

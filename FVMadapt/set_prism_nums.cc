#include <iostream>
#include <cstdlib>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <Loci.h>
#include "prism.h"

using std::queue;
using std::cerr;
using std::cout;
using std::endl;
using std::set;
using std::map;

class set_prism_num_nodes : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > orientCode;
  const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  store<int> num_inner_nodes;

public:
  set_prism_num_nodes(){
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prismOrientCode", orientCode);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
   
    name_store("num_inner_nodes", num_inner_nodes);


    input("cellPlan");
    input("(prism2face, prism2node, prismOrientCode)");
    input("(lower, upper, boundary_map)->facePlan");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node");
    input("(lower, upper, boundary_map)->face2edge");
    
    output("num_inner_nodes");

    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
       
    if(cellPlan[cc].size() == 0){
      num_inner_nodes[cc] = 0;

      return;
    }

    std::list<Node*> bnode_list; //boundary node
    std::list<Node*> node_list; //inner node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> qface_list;
    std::list<Face*> gface_list;
    Prism* aCell = build_prism_cell(lower[cc].begin(), lower.num_elems(cc),
                                    upper[cc].begin(), upper.num_elems(cc),
                                    boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                    prism2face[cc],
                                    prism2node[cc],
                                    orientCode[cc],
                                    face2node,
                                    face2edge,
                                    edge2node,
                                    pos,
                                    edgePlan,
                                    facePlan,
                                    bnode_list,
                                    edge_list,
                                    qface_list,
                                    gface_list);
    
    std::vector<Prism*> cells;
    aCell->resplit( cellPlan[cc],
                    node_list,
                    edge_list,
                    qface_list,
                    gface_list,
                    cells);
  
  
  
    num_inner_nodes[cc] = node_list.size();

    
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    
    //aCell will clean up these
    cleanup_list(bnode_list, edge_list, gface_list);
    cleanup_list(qface_list);
    cleanup_list(node_list);
  }
};
register_rule<set_prism_num_nodes> register_set_prism_num_nodes;  


class set_prism_num_cells : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > orientCode;
  const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;

  store<int> num_fine_cells;
public:
  set_prism_num_cells(){
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prismOrientCode", orientCode);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
   

    name_store("num_fine_cells", num_fine_cells);

    input("cellPlan");
    input("(prism2face, prism2node, prismOrientCode)");
    input("(lower, upper, boundary_map)->facePlan");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node");
    input("(lower, upper, boundary_map)->face2edge");
    

    output("num_fine_cells");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
       
    if(cellPlan[cc].size() == 0){

      num_fine_cells[cc] = 1;
      return;
    }

    std::list<Node*> bnode_list; //boundary node
    std::list<Node*> node_list; //inner node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> qface_list;
    std::list<Face*> gface_list;
    Prism* aCell = build_prism_cell(lower[cc].begin(), lower.num_elems(cc),
                                    upper[cc].begin(), upper.num_elems(cc),
                                    boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                    prism2face[cc],
                                    prism2node[cc],
                                    orientCode[cc],
                                    face2node,
                                    face2edge,
                                    edge2node,
                                    pos,
                                    edgePlan,
                                    facePlan,
                                    bnode_list,
                                    edge_list,
                                    qface_list,
                                    gface_list);
    
    std::vector<Prism*> cells;
    aCell->resplit( cellPlan[cc],
                    node_list,
                    edge_list,
                    qface_list,
                    gface_list,
                    cells);
  
  
  

    num_fine_cells[cc] = cells.size();
    
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    
    //aCell will clean up these
    cleanup_list(bnode_list, edge_list, gface_list);
    cleanup_list(qface_list);
    cleanup_list(node_list);
  }
};
register_rule<set_prism_num_cells> register_set_prism_num_cells;  



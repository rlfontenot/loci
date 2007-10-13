#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <Loci.h>
#include <algorithm>
#include "diamondcell.h"

#include "globals.h"

using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;
using Loci::storeRepP;
//int currentMem(void);

//this rule make  a newCellPlan according to cellPlan 
//and nodeTag, posTag
class make_general_cellplan:public pointwise_rule{
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<bool> is_quadface;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<char>  posTag;
  const_store<std::vector<char> > nodeTag;
  const_store<bool> isIndivisible;
  store<std::vector<char> > newCellPlan;
  const_blackbox<storeRepP>  node_remap;
  Map node_l2f;
public:
  make_general_cellplan(){
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("posTag", posTag);
    name_store("nodeTag", nodeTag);
    name_store("isIndivisible", isIndivisible);
    name_store("priority::restart_par::newCellPlan", newCellPlan);
    name_store("node_remap", node_remap);
    name_store("is_quadface", is_quadface);
    input("(cellPlan,nodeTag) ");
    input("isIndivisible");
    input("(lower, upper, boundary_map) -> (facePlan,is_quadface, nodeTag)"); 
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->(edgePlan, nodeTag)");
    input("node_remap");
    output("priority::restart_par::newCellPlan");
    constraint("gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
   
  }
  void calculate(Entity cc){

    if(!isIndivisible[cc]){
      std::list<Node*> node_list;
      std::list<Edge*> edge_list;
      std::list<Face*> face_list;
      std::list<Node*> bnode_list;
      int nindex;

    
 

                                                 
      Cell* aCell = build_general_cell(lower[cc].begin(), lower.num_elems(cc),
                                       upper[cc].begin(), upper.num_elems(cc),
                                       boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                       is_quadface,
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       pos,
                                       edgePlan,
                                       facePlan,
                                       posTag,
                                       nodeTag,
                                       bnode_list,
                                       edge_list,
                                       face_list,
                                       node_l2f);

  
  
 
      //calculate min_edge_length of the cell
      double min_edge_length =aCell->get_min_edge_length();
      
      
      std::vector<DiamondCell*> cells;
      
      aCell->resplit( cellPlan[cc], 
                      node_list,
                      edge_list,
                      face_list,
                      cells);
      
      
      nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++){
        (*np)->tag = nodeTag[cc][nindex++];
      }
      
      
      
    
      // split the cell Globals::levels times
      int numCells = cells.size();
      
      if(numCells != 0){
        bool cell_split = false;
        for(int i = 0; i < numCells; i++){
        if(cells[i]->get_tagged()) {
          if((min_edge_length/double(1<<(cells[i]->getLevel()+Globals::levels))) > Globals::tolerance){
            cells[i]->resplit(Globals::levels,node_list, edge_list, face_list);
            cell_split = true;
          }
          else{
            int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0) - cells[i]->getLevel());
            if(split_level >= 1){
              cells[i]->resplit(split_level,node_list, edge_list, face_list);
              cell_split = true;
            }
          }
        }
        }
        if(cell_split)aCell->rebalance_cells(node_list, edge_list, face_list);
      }
      
      else{
        if(aCell->get_tagged()){
          
          
          int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
          if(split_level >= 1){
          aCell->split(node_list, edge_list, face_list);
          for(int i = 0; i < aCell->numNode; i++){
            aCell->child[i]->resplit(min(Globals::levels, split_level)-1, node_list, edge_list, face_list);
          }
          }
        }
      }
      //write new cellPlan
      newCellPlan[cc] = aCell->make_cellplan();
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(node_list, edge_list, face_list);
      cleanup_list(bnode_list);
      reduce_vector(newCellPlan[cc]);
    }
  }
};

register_rule<make_general_cellplan> register_make_general_cellplan;

class make_general_cellplan_norestart:public pointwise_rule{
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<char>  posTag;
  const_store<bool> isIndivisible;
  
  store<std::vector<char> > newCellPlan;
  const_blackbox<storeRepP> node_remap;
    Map node_l2f;
public:
  make_general_cellplan_norestart(){
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("posTag", posTag);
   
    name_store("isIndivisible", isIndivisible);
    name_store("restart_par::newCellPlan", newCellPlan);
    name_store("node_remap", node_remap);
    
  
    input("isIndivisible");
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("node_remap");
    output("restart_par::newCellPlan");
    constraint("gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }

    
  }
  void calculate(Entity cc){
    if(!isIndivisible[cc]){
    std::list<Node*> node_list;
    std::list<Edge*> edge_list;
    std::list<Face*> face_list;
    std::list<Node*> bnode_list;
  
    
    Cell* aCell = build_general_cell(lower[cc].begin(), lower.num_elems(cc),
                                     upper[cc].begin(), upper.num_elems(cc),
                                     boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                     face2node,
                                     face2edge,
                                     edge2node,
                                     pos,
                                     posTag,
                                     bnode_list,
                                     edge_list,
                                     face_list,
                                     node_l2f);

  
  
 
     //calculate min_edge_length of the cell
    double min_edge_length =aCell->get_min_edge_length();
    
       // split the cell Globals::levels times
    
   
   
    if(aCell->get_tagged()){
      int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
      
      if(split_level >= 1){
        aCell->split(node_list, edge_list, face_list);
        for(int i = 0; i < aCell->numNode; i++){
          aCell->child[i]->resplit(min(Globals::levels, split_level)-1, node_list, edge_list, face_list);
        }
      }
    }
      
    //write new cellPlan
    newCellPlan[cc] = aCell->make_cellplan();
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    cleanup_list(node_list, edge_list, face_list);
    cleanup_list(bnode_list);
    reduce_vector(newCellPlan[cc]);
    }
  }
};

register_rule<make_general_cellplan_norestart> register_make_general_cellplan_norestart;


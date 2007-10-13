#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <bitset>
#include "diamondcell.h"
#include <Loci.h>
#include <algorithm>
using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;


class merge_general_interior_tmpface:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_store<vect3d> pos;
  const_blackbox<Loci::storeRepP> node_remap;
   Map node_l2f;
public:
  merge_general_interior_tmpface(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);
    name_store("pos", pos);
    name_store("iface_remap", node_remap);
    input("iface_remap");
    input("(cl,cr)->tmpCellPlan");
    input("(cl, cr)->(lower, upper, boundary_map)->face2node->pos");
    input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("tmpFacePlan");
     constraint("interior_faces, (cl, cr)->gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
 
  }
  void calculate(Entity f){
   
   
    facePlan[f].clear();
   
    std::vector<char> facePlanL = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                                       boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cl[f]],
                                                       f,
                                                       node_l2f);
    
    std::vector<char> facePlanR = extract_general_face(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                                       upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                                       boundary_map[cr[f]].begin(), boundary_map.num_elems(cr[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cr[f]],
                                                       f,
                                                       node_l2f);
    
    
    
    facePlan[f] = merge_faceplan(facePlanL, facePlanR, face2node.num_elems(f));
    reduce_vector(facePlan[f]);     
   
  }
};

register_rule<merge_general_interior_tmpface> register_merge_general_interior_tmpface;



class merge_general_boundary_tmpface:public pointwise_rule{
  const_Map cl;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_store<vect3d> pos;//dummy
  const_blackbox<Loci::storeRepP> node_remap;
   Map node_l2f;
public:
  merge_general_boundary_tmpface(){
    name_store("cl", cl);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);
    name_store("pos", pos);
    name_store("bface_remap", node_remap);
    input("bface_remap");
    input("cl->tmpCellPlan");
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("tmpFacePlan");
    constraint("boundary_faces, cl->gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      node_l2f = *node_remap;
      do_loop(seq, this);
    }

  }
  void calculate(Entity f){
   
   
    facePlan[f] = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                       boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       cellPlan[cl[f]],
                                       f,
                                       node_l2f);  
    reduce_vector(facePlan[f]);
 
  }
  
};

register_rule<merge_general_boundary_tmpface> register_merge_general_boundary_tmpface;

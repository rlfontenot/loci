


#include <queue>
#include <vector>
#include <utility>
#include <Loci.h>
#include "hex_defines.h"
#include "defines.h"
#include "face.h"


using std::vector;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;


std::vector<char>   extract_general_face(const Entity* lower, int lower_size,
                                         const Entity* upper, int upper_size,
                                         const Entity* boundary_map, int boundary_map_size,
                                         const const_multiMap& face2node,
                                         const const_multiMap& face2edge,
                                         const const_MapVec<2>& edge2node,
                                         const std::vector<char>& cellPlan,
                                         Entity ff,
                                         const Map& node_remap
                                         );
std::vector<char>  extract_prism_face(const  std::vector<char>& cellPlan, int dd);
std::vector<char>  merge_tri_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1);
std::vector<char>  extract_hex_face(const  std::vector<char>& cellPlan,  DIRECTION dd);

std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL);

std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL,
                                  std::vector<char>& facePlanR, char orientCodeR);

std::vector<char>  merge_tri_face_pp(const  std::vector<char>& cellPlan1, int dd1, char orientCode1,
                                     const  std::vector<char>& cellPlan2, int dd2, char orientCode2);


std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes);

std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan);

  
    
  

class merge_face_gh:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fr;
  const_store<Array<char,6> > orientCode;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos; //dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_blackbox<Loci::storeRepP> node_remap;
  Map node_l2f;
public:
  merge_face_gh(){
    name_store("cl", cl);
    name_store("cr", cr);
     name_store("fr", fr);
    name_store("hexOrientCode", orientCode);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("iface_remap", node_remap);
    input("iface_remap");
    input("(cl,cr)->cellPlan");
    input("cr->hexOrientCode");
    input("fr");
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("cl->gnrlcells, cr->hexcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
    
  }
  void calculate(Entity f){
   
   
    
    facePlan[f].clear();
    std::vector<char> localPlanL = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                                       boundary_map[cl[f]].begin(),
                                                       boundary_map.num_elems(cl[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cl[f]],
                                                       f,
                                                       node_l2f);

    std::vector<char> facePlanL = transfer_plan_g2q(localPlanL);
    std::vector<char> facePlanR = extract_hex_face(cellPlan[cr[f]], DIRECTION(fr[f]));

    
    
    facePlan[f] = merge_quad_face(facePlanL, char(0),
                                  facePlanR,  orientCode[cr[f]][fr[f]]);  
     reduce_vector(facePlan[f]);
  }
};

register_rule<merge_face_gh> register_merge_face_gh;


class merge_face_hg:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fl;
  const_store<Array<char,6> > orientCode;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos; //dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_blackbox<Loci::storeRepP> node_remap;
  Map node_l2f;
public:
  merge_face_hg(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("fl", fl);
    name_store("hexOrientCode", orientCode);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("iface_remap", node_remap);
    input("iface_remap");
    input("(cl,cr)->cellPlan");
    input("cl->hexOrientCode");
    input("fl");
    input(" cr->(lower, upper, boundary_map)->face2node->pos");
    input(" cr->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("cr->gnrlcells, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
    
  }
  void calculate(Entity f){
   
   
    
    facePlan[f].clear();
    std::vector<char> localPlanR = extract_general_face(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                                       upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                                       boundary_map[cr[f]].begin(),
                                                       boundary_map.num_elems(cr[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cr[f]],
                                                       f,
                                                       node_l2f);
    std::vector<char> facePlanR = transfer_plan_g2q(localPlanR);
    std::vector<char> facePlanL = extract_hex_face(cellPlan[cl[f]], DIRECTION(fl[f]));
    
    facePlan[f] = merge_quad_face(facePlanR, char(0),
                                  facePlanL, orientCode[cl[f]][fl[f]]);  
    reduce_vector(facePlan[f]);
  }
};

register_rule<merge_face_hg> register_merge_face_hg;

class merge_face_gp:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fr;
  const_store<Array<char,5> > orientCode;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos; //dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_blackbox<Loci::storeRepP> node_remap;
   Map node_l2f;
public:
  merge_face_gp(){
    name_store("cl", cl);
    name_store("cr", cr);
     name_store("fr", fr);
    name_store("prismOrientCode", orientCode);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("iface_remap", node_remap);
    input("iface_remap");
    input("(cl,cr)->cellPlan");
     input("cr->prismOrientCode");
    input("fr");
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("cl->gnrlcells, cr->prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
    
  }
  void calculate(Entity f){
   
   
    
    facePlan[f].clear();
    std::vector<char> localPlanL = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                                       boundary_map[cl[f]].begin(),
                                                       boundary_map.num_elems(cl[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cl[f]],
                                                       f,
                                                       node_l2f);

   
    
    
    if(fr[f]< 2){
     std::vector<char> facePlanR = merge_tri_face_p(cellPlan[cr[f]], fr[f], orientCode[cr[f]][fr[f]]);
      facePlan[f] = merge_faceplan(localPlanL, facePlanR, face2node.num_elems(f));
    }
    else{
      std::vector<char> facePlanR =  extract_prism_face(cellPlan[cr[f]], fr[f]);
       std::vector<char> facePlanL = transfer_plan_g2q(localPlanL);
      facePlan[f] = merge_quad_face(facePlanL, char(0),
                                    facePlanR, orientCode[cr[f]][fr[f]]);
    }
    reduce_vector(facePlan[f]);
  }
};

register_rule<merge_face_gp> register_merge_face_gp;

class merge_face_pg:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fl;
  const_store<Array<char,5> > orientCode;
    const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos; //dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
  const_blackbox<Loci::storeRepP> node_remap;
   Map node_l2f;
public:
  merge_face_pg(){
    name_store("cl", cl);
    name_store("cr", cr);
     name_store("fl", fl);
    name_store("prismOrientCode", orientCode);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("iface_remap", node_remap);
    input("iface_remap");
    input("(cl,cr)->cellPlan");
     input("cl->prismOrientCode");
    input("fl");
    input("cr->(lower, upper, boundary_map)->face2node->pos");
    input("cr->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("cr->gnrlcells, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    node_l2f = *node_remap;
    do_loop(seq, this);
    }
    
  }
  void calculate(Entity f){
   
   
    
    facePlan[f].clear();
    std::vector<char> localPlanR = extract_general_face(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                                       upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                                       boundary_map[cr[f]].begin(),
                                                       boundary_map.num_elems(cr[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cr[f]],
                                                       f,
                                                       node_l2f);

   
    
    if(fl[f]< 2){
      std::vector<char> facePlanL = merge_tri_face_p(cellPlan[cl[f]], fl[f], orientCode[cl[f]][fl[f]]);
      facePlan[f] = merge_faceplan(facePlanL, localPlanR, face2node.num_elems(f));
    }
    else{
      std::vector<char> facePlanL =  extract_prism_face(cellPlan[cl[f]], fl[f]);
      std::vector<char> facePlanR = transfer_plan_g2q(localPlanR);

      facePlan[f] = merge_quad_face(facePlanR,char(0),
                                    facePlanL, orientCode[cl[f]][fl[f]]);
    }
    reduce_vector(facePlan[f]);
  }
};

register_rule<merge_face_pg> register_merge_face_pg;


class merge_face_hp:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fl;
  const_store<char> fr;
  const_store<Array<char,5> > orientCodeP;
  const_store<Array<char,6> > orientCodeH;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
public:
  merge_face_hp(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("fl", fl);
    name_store("fr", fr);
    name_store("prismOrientCode", orientCodeP);
    name_store("hexOrientCode", orientCodeH);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    input("(cl,cr)->cellPlan");
    input("cr->prismOrientCode");
    input("cl->hexOrientCode");
    input("(fl, fr)");
    output("facePlan");
    constraint("cl->hexcells, cr->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    std::vector<char> facePlanL = extract_hex_face(cellPlan[cl[f]], DIRECTION(fl[f]));
    std::vector<char> facePlanR = extract_prism_face(cellPlan[cr[f]], fr[f]);
    facePlan[f] = merge_quad_face(facePlanL, orientCodeH[cl[f]][fl[f]],
                                  facePlanR, orientCodeP[cr[f]][fr[f]]);
    reduce_vector(facePlan[f]);                                                                    
  }
};

register_rule<merge_face_hp> register_merge_face_hp;

class merge_face_ph:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_store<char> fl;
  const_store<char> fr;
  const_store<Array<char,5> > orientCodeP;
   const_store<Array<char,6> > orientCodeH;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
public:
  merge_face_ph(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("fl", fl);
    name_store("fr", fr);
    name_store("prismOrientCode", orientCodeP);
    name_store("hexOrientCode", orientCodeH);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    input("(cl,cr)->cellPlan");
    input("cl->prismOrientCode");
    input("cr->hexOrientCode");
    input("(fl, fr)");
    output("facePlan");
    constraint("cr->hexcells, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    std::vector<char> facePlanR = extract_hex_face(cellPlan[cr[f]], DIRECTION(fr[f]));
    std::vector<char> facePlanL = extract_prism_face(cellPlan[cl[f]], fl[f]);
    facePlan[f] = merge_quad_face(facePlanL, orientCodeP[cl[f]][fl[f]],
                                  facePlanR, orientCodeH[cr[f]][fr[f]]);
    reduce_vector(facePlan[f]);                                                                    
  }
};

register_rule<merge_face_ph> register_merge_face_ph;


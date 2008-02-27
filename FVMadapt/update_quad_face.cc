#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <bitset>
#include "hex_defines.h"
#include <Loci.h>

using std::list;
using std::queue;
using std::bitset;
using std::cerr;
using std::endl;
using std::cout;
void extract_quad_edge(const std::vector<char>&, std::vector<char>&, unsigned int);
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL);
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL,
                                  std::vector<char>& facePlanR, char orientCodeR);

std::vector<char>  extract_hex_face(const  std::vector<char>& cellPlan,  DIRECTION dd);

//this rule first extract faceplan from cl and cr, then merge them into one final plan 
//in the face's local coordinates system
//
//to merge the plans,  use a queue of ranges to simulate the tree splitting
//process, and put the inner eges in list xEdge and yEdges according to their directions.
//
//then sort and unique xEdge and yEdge. and then merge all the overlapped edges   
//
//The last step is compose a new refinement plan according to xEdges and yEdges. Search 
//the inner edges of a parent cell, the result implies the split code. repeat the same process for
//each child.

class merge_hex_interior_tmpface:public pointwise_rule{
  const_Map cl;
  const_Map cr;
 const_store<char> fl;
  const_store<char> fr;
  const_store<Array<char,6> > orientCode;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
 
public:
  merge_hex_interior_tmpface(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("fl", fl);
    name_store("fr", fr);
    name_store("hexOrientCode", orientCode);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);
   
    input("(cl,cr)->tmpCellPlan");
    input("(cl,cr)->hexOrientCode");
    input("(fl, fr)");
    output("tmpFacePlan");
     constraint("interior_faces, (cl, cr)->hexcells");
  }
  virtual void compute(const sequence &seq){
   
    do_loop(seq, this);
   
 
  }
  void calculate(Entity f){
     
    std::vector<char> facePlanL = extract_hex_face(cellPlan[cl[f]],DIRECTION(fl[f]));
    std::vector<char> facePlanR = extract_hex_face(cellPlan[cr[f]], DIRECTION(fr[f]));
    facePlan[f] = merge_quad_face(facePlanL, orientCode[cl[f]][fl[f]],
                                  facePlanR, orientCode[cr[f]][fr[f]]);
    reduce_vector(facePlan[f]);     
    
  }
};

register_rule<merge_hex_interior_tmpface> register_merge_hex_interior_tmpface;



class merge_hex_boundary_tmpface:public pointwise_rule{
  const_Map cl;
  const_store<char> fl;
   const_store<Array<char,6> > orientCode;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
public:
  merge_hex_boundary_tmpface(){
    name_store("cl", cl);
    name_store("fl", fl);
    name_store("hexOrientCode", orientCode);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);

    input("cl->tmpCellPlan");
    input("cl->hexOrientCode");
    input("fl");
    output("tmpFacePlan");
    constraint("boundary_faces, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    
    do_loop(seq, this);
    
  }
  void calculate(Entity f){
    std::vector<char> facePlanL = extract_hex_face( cellPlan[cl[f]], DIRECTION(fl[f]));
    facePlan[f] = merge_quad_face(facePlanL, orientCode[cl[f]][fl[f]]);
    reduce_vector(facePlan[f]);
 
  }
  
};

register_rule<merge_hex_boundary_tmpface> register_merge_hex_boundary_tmpface;

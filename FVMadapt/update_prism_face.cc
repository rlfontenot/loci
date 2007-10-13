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

void  extract_prism_face(const  std::vector<char>& cellPlan, std::vector<char>& facePlan, int dd);

std::vector<char>  merge_quad_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1);

std::vector<char>  merge_quad_face_pp(const  std::vector<char>& cellPlan1, int dd1, char orientCode1,
                                      const  std::vector<char>& cellPlan2, int dd2, char orientCode2);

std::vector<char>  merge_tri_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1);

std::vector<char>  merge_tri_face_pp(const  std::vector<char>& cellPlan1, int dd1, char orientCode1,
                                      const  std::vector<char>& cellPlan2, int dd2, char orientCode2);



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

class merge_prism_interior_tmpface:public pointwise_rule{
  const_Map cl;
  const_Map cr;
 const_store<char> fl;
  const_store<char> fr;
  const_store<Array<char,5> > orientCode;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
 
public:
  merge_prism_interior_tmpface(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("fl", fl);
    name_store("fr", fr);
    name_store("prismOrientCode", orientCode);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);
   
    input("(cl,cr)->tmpCellPlan");
    input("(cl,cr)->prismOrientCode");
    input("(fl, fr)");
    output("tmpFacePlan");
     constraint("(cl, cr)->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
    
 
  }
  void calculate(Entity f){
   if(fl[f]< 2)facePlan[f] = merge_tri_face_pp(cellPlan[cl[f]], fl[f], orientCode[cl[f]][fl[f]],
                                                cellPlan[cr[f]], fr[f], orientCode[cr[f]][fr[f]]);
    else facePlan[f] = merge_quad_face_pp(cellPlan[cl[f]], fl[f], orientCode[cl[f]][fl[f]],
                                          cellPlan[cr[f]], fr[f], orientCode[cr[f]][fr[f]]);
     reduce_vector(facePlan[f]);      
  
  }
};

register_rule<merge_prism_interior_tmpface> register_merge_prism_interior_tmpface;



class merge_prism_boundary_tmpface:public pointwise_rule{
  const_Map cl;
  const_store<char> fl;
   const_store<Array<char,5> > orientCode;
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;
public:
  merge_prism_boundary_tmpface(){
    name_store("cl", cl);
    name_store("fl", fl);
    name_store("prismOrientCode", orientCode);
    name_store("tmpCellPlan", cellPlan);
    name_store("tmpFacePlan", facePlan);

    input("cl->tmpCellPlan");
    input("cl->prismOrientCode");
    input("fl");
    output("tmpFacePlan");
    constraint("boundary_faces, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    if(fl[f]< 2) facePlan[f] = merge_tri_face_p(cellPlan[cl[f]], fl[f], orientCode[cl[f]][fl[f]]);
    
    else  facePlan[f] =  merge_quad_face_p(cellPlan[cl[f]], fl[f], orientCode[cl[f]][fl[f]]);   
    reduce_vector(facePlan[f]);
 
  }
  
};

register_rule<merge_prism_boundary_tmpface> register_merge_prism_boundary_tmpface;

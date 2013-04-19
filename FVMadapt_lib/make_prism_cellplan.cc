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
#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <Loci.h>
#include <algorithm>
#include "prism.h"

#ifdef USE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#endif

#include "globals.h"

using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;
using Loci::storeRepP;
//int currentMem(void);
#ifdef USE_LIBXML2
bool mark_node( xmlNode* root_element,
                std::list<Node*>::iterator begin_pnt,
                std::list<Node*>::iterator end_pnt);
#endif

//this rule make  a newCellPlan according to cellPlan 
//and nodeTag, posTag
class make_prism_cellplan:public pointwise_rule{
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<char>  posTag;
  const_store<std::vector<char> > nodeTag;
  const_store<bool> isIndivisible;
  const_param<int> split_mode_par;
  const_param<int> restart_tag_par;
  store<std::vector<char> > newCellPlan;

  const_store<int> node_l2f;
 
public:
  make_prism_cellplan(){
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
    name_store("posTag", posTag);
    name_store("nodeTag", nodeTag);
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);
    name_store("restart_tag_par", restart_tag_par);
 
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("split_mode_par");
    input("(cellPlan, nodeTag, prism2face, prism2node, prismOrientCode)");
    input("isIndivisible");
    input("(lower, upper, boundary_map) -> (fileNumber(face2node),facePlan, nodeTag)"); 
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->(edgePlan, nodeTag)");
    input("restart_tag_par");
    
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
  
    if(seq.size()!=0){
       
      do_loop(seq, this);
    }
  
     
  }
  void calculate(Entity cc){
    if(!isIndivisible[cc]){
      std::list<Node*> node_list;
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
      int nindex;

    
  

                                                 
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
                                      posTag,
                                      nodeTag,
                                      bnode_list,
                                      edge_list,
                                      qface_list,
                                      gface_list,
                                      node_l2f
                                      );
    
  
  
 
 
 
    
   
      std::vector<Prism*> cells;
    
      aCell->resplit( cellPlan[cc], 
                      node_list,
                      edge_list,
                      qface_list,
                      gface_list,
                      cells);
    
      nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++){
        (*np)->tag = nodeTag[cc][nindex++];
      }
      
    
      int numCells = cells.size();
      bool cell_split = false; //if any cell split, the whole big cell rebalance
      bool cell_merged = false;
      if(numCells != 0){//aCell is not a leaf
       
        //first if any cell need derefine
        std::set<Prism*> dparents;
        //mark the cell that will be eliminated
        for(int i = 0; i < numCells; i++){
          if(cells[i]->get_tagged() ==2){
            Prism* parent = cells[i]->getParentCell();
            if(parent!=0 && parent->needDerefine()){
              dparents.insert(parent);
              cells[i] =0;
            }
          }
        }
       
        //refine the other cells
        for(int i = 0; i < numCells; i++){
          if(cells[i] ==0) continue;//derefined cell
          if((cells[i]->get_tagged())==1) {
            cells[i]->setSplitCode(*split_mode_par, Globals::tolerance);
            if(cells[i]->getMySplitCode()!=0) {
              if(*split_mode_par == 2){
                double min_edge_length =cells[i]->get_min_edge_length();
                int split_level = Globals::levels;
                if(Globals::tolerance > 0.0) split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));  
                cells[i]->resplit(min(Globals::levels,split_level),node_list, edge_list, qface_list, gface_list);
                cell_split = true;
              }else{
                cell_split = true;
                cells[i]->split(node_list, edge_list, qface_list, gface_list);
              }
            }
          }
        }
        //derefine the cells
        for(std::set<Prism*>::const_iterator si = dparents.begin(); si!= dparents.end(); si++){
          (*si)->derefine();
          cell_merged = true;
        }
      
    }else{//aCell is a leaf
   
      aCell->setSplitCode(*split_mode_par, Globals::tolerance);
      if(aCell->getMySplitCode() != 0 ){
        if(*split_mode_par == 2){
          double min_edge_length =aCell->get_min_edge_length();
          int split_level = Globals::levels;
          if(Globals::tolerance > 0.0) split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
          aCell->resplit(min(Globals::levels, split_level), node_list, edge_list, qface_list, gface_list); 
        }
        else{
          aCell->split(node_list, edge_list, qface_list, gface_list);
        }
      }
    }
      
    //rebalance this cell
    if(cell_merged){
      //record current plan
      vector<char> tmpPlan = aCell->make_cellplan();
      //build a temp cell
      Prism* tmpCell = build_prism_cell(lower[cc].begin(), lower.num_elems(cc),
                                        upper[cc].begin(), upper.num_elems(cc),
                                        boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                        prism2face[cc],
                                        prism2node[cc],
                                        orientCode[cc],
                                        face2node,
                                        face2edge,
                                        edge2node,
                                        pos,
                                        bnode_list,
                                        edge_list,
                                        qface_list,
                                        gface_list,
                                        node_l2f);
        
      //resplit temp cell according to current plan 
      tmpCell->resplit( tmpPlan, 
                        node_list,
                        edge_list,
                        qface_list,
                        gface_list,
                        cells);
      //rebalance temp cell
      tmpCell->rebalance_cells(*split_mode_par,node_list, edge_list, qface_list, gface_list);
      //record balanced plan
      newCellPlan[cc] = tmpCell->make_cellplan();
      //clean up 
      delete tmpCell;
        
    }else if(cell_split){
      aCell->rebalance_cells(*split_mode_par,node_list, edge_list, qface_list, gface_list);
      newCellPlan[cc] = aCell->make_cellplan();
    }else{
      //write new cellPlan
      newCellPlan[cc] = aCell->make_cellplan();
    }
      
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    cleanup_list(node_list, edge_list, qface_list);
    cleanup_list(bnode_list);
    cleanup_list(gface_list);
    reduce_vector(newCellPlan[cc]);
  }
}
  };

register_rule<make_prism_cellplan> register_make_prism_cellplan;

class make_prism_cellplan_norestart:public pointwise_rule{
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<bool> isIndivisible;
  const_store<char>  posTag;
  const_param<int> split_mode_par;
  const_param<int> norestart_tag_par;
  store<std::vector<char> > newCellPlan;

  const_store<int > node_l2f;
public:
  make_prism_cellplan_norestart(){
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("posTag", posTag);
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);

    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);

    name_store("norestart_tag_par", norestart_tag_par);
    input("split_mode_par"); 
    input("isIndivisible, prism2face, prism2node, prismOrientCode");
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    input("norestart_tag_par");
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
       
      do_loop(seq, this);
    }
    
  }
  void calculate(Entity cc){
    if(!isIndivisible[cc]){   
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
     
    
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
                                      posTag,
                                      bnode_list,
                                      edge_list,
                                      qface_list,
                                      gface_list,
                                      node_l2f);


      aCell->setSplitCode(*split_mode_par, Globals::tolerance);
  
  
      
      //write new cellPlan
      if((aCell->getMySplitCode()) != 0  && *split_mode_par == 2){
        double min_edge_length =aCell->get_min_edge_length();
        int split_level =  Globals::levels;
        if(Globals::tolerance > 0.0) split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
        newCellPlan[cc] =  aCell->make_cellplan(min(Globals::levels, split_level));
      }
      else if(aCell->getMySplitCode() != 0   && *split_mode_par != 2){
        newCellPlan[cc].push_back(aCell->getMySplitCode());
      }
    
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, qface_list);
      cleanup_list(gface_list);
      reduce_vector(newCellPlan[cc]);
    }
  }
};

register_rule<make_prism_cellplan_norestart> register_make_prism_cellplan_norestart;

//this rule make  a newCellPlan according to cellPlan 
//and xml file
class make_prism_cellplan_xml:public pointwise_rule{
  const_param<std::string> xmlfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  
  const_store<bool> isIndivisible;
  const_param<int> split_mode_par;
  const_param<int> restart_xml_par;
  store<std::vector<char> > newCellPlan;
#ifdef USE_LIBXML2
  xmlDoc* doc ;
  xmlNode* root_element ;
#endif

  const_store<int>  node_l2f;
public:
  make_prism_cellplan_xml(){
    name_store("xmlfile_par", xmlfile_par);
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
  
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);
    name_store("restart_xml_par", restart_xml_par);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("restart_xml_par");
    input("(split_mode_par,xmlfile_par)");

    input("(cellPlan, prism2face, prism2node, prismOrientCode)");
    input("isIndivisible");
    input("(lower, upper, boundary_map) ->(fileNumber(face2node), facePlan)"); 
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");

    
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
#ifdef USE_LIBXML2
      doc = xmlReadFile((*xmlfile_par).c_str(), NULL, 0);
      if(doc==NULL){
        cerr << "WARNING: fail tp parse xml file" << endl;
        Loci::Abort();
      }
      root_element = xmlDocGetRootElement(doc);
      if(root_element == NULL) {
        cerr <<"WARNING: fail to parse xml file" << endl;
        Loci::Abort();
      }
      do_loop(seq, this);
      xmlFreeDoc(doc);
      xmlCleanupParser();
      xmlMemoryDump();
#else
      cerr << "libxml2 not available when adapt code compiled.  XML features not enabled!" << endl ;
      Loci::Abort() ;
#endif
    }
 
  }
  void calculate(Entity cc){
#ifdef USE_LIBXML2
    if(!isIndivisible[cc]){
   
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
      std::queue<Prism*> Q;

    
  

                                                 
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
                                      gface_list,
                                      node_l2f
                                      );
    
  
  
      if(!mark_node(root_element, bnode_list.begin(), bnode_list.end())){
        cerr << "WARNING: fail to mark nodes, please check xml file" << endl;
        Loci::Abort();
      }
 
 
 
    
   
      std::vector<Prism*> cells;
      std::list<Node*>::iterator former_pnt = bnode_list.end(); 
      aCell->resplit( cellPlan[cc], 
                      bnode_list,
                      edge_list,
                      qface_list,
                      gface_list,
                      cells);
      former_pnt--;
      if(!mark_node(root_element, bnode_list.begin(), bnode_list.end())){
        cerr << "WARNING: fail to mark nodes, please check xml file" << endl;
        Loci::Abort();
      }
    
    
      Prism* current;
      int numCells = cells.size();
    
      if(numCells != 0){
        for(int i = 0; i < numCells; i++)Q.push(cells[i]);
      }
      else{
        Q.push(aCell);
      }

      while(!Q.empty()){
        current =Q.front();
        //calculate min_edge_length of the cell
        

      
      
        if( current->getMySplitCode()==0)current->setSplitCode(*split_mode_par, Globals::tolerance);
        if(current->getMySplitCode()!=0){
          std::list<Node*>::iterator former_pnt = bnode_list.end();
          current->split(bnode_list, edge_list, qface_list, gface_list);
        

          former_pnt--;
          if(!mark_node(root_element, bnode_list.begin(), bnode_list.end())){
            cerr << "WARNING: fail to mark nodes, please check xml file" << endl;
            Loci::Abort();
          }
          for(int i = 0; i < current->numChildren(); i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(*split_mode_par, bnode_list, edge_list, qface_list, gface_list);
        
      //write new cellPlan
      newCellPlan[cc] = aCell->make_cellplan();
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, qface_list);
    
      cleanup_list(gface_list);
      reduce_vector(newCellPlan[cc]);
    
    
    }
#endif
  }
};

register_rule<make_prism_cellplan_xml> register_make_prism_cellplan_xml;

class make_prism_cellplan_xml_norestart:public pointwise_rule{
  const_param<std::string> xmlfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<bool> isIndivisible;
   
  const_param<int> split_mode_par;
  const_param<int> norestart_xml_par;
  store<std::vector<char> > newCellPlan;
#ifdef USE_LIBXML2
  xmlDoc* doc ;
  xmlNode* root_element ;
#endif

  const_store<int> node_l2f;
public:
  make_prism_cellplan_xml_norestart(){
    name_store("xmlfile_par", xmlfile_par);
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
     
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);
    name_store("norestart_xml_par", norestart_xml_par);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("norestart_xml_par");
    input("split_mode_par, xmlfile_par"); 
    input("isIndivisible, prism2face, prism2node, prismOrientCode");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
#ifdef USE_LIBXML2
      doc = xmlReadFile((*xmlfile_par).c_str(), NULL, 0);
      if(doc==NULL){
        cerr<< "WARNING: fail to parse xml file"<< endl;
        Loci::Abort();
      }
      root_element = xmlDocGetRootElement(doc);
      if(root_element==NULL) {
        cerr <<"WARNING: fail to  parse xml file" << endl;
        Loci::Abort();
      }
      do_loop(seq, this);
      xmlFreeDoc(doc);
      xmlCleanupParser();
      xmlMemoryDump();
#else
      cerr << "libxml2 not available when compiling adapt code.  XML features not enabled!" << endl ;
      Loci::Abort() ;
#endif
    }

  
  }
  void calculate(Entity cc){
#ifdef USE_LIBXML2
    if(!isIndivisible[cc]){   
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
      std::queue<Prism*> Q; 
    
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
                                      bnode_list,
                                      edge_list,
                                      qface_list,
                                      gface_list,
                                      node_l2f);
      if( !mark_node(root_element, bnode_list.begin(), bnode_list.end())){
        cerr << "WARNING: fail to mark nodes, please check xml file" << endl;
        Loci::Abort();
      }
      Q.push(aCell);
  
      Prism* current;
      while(!Q.empty()){
        current =Q.front();
     
     
        if(current->getMySplitCode()==0) current->setSplitCode(*split_mode_par, Globals::tolerance);
        if(current->getMySplitCode()!=0){
          std::list<Node*>::iterator former_pnt = bnode_list.end();
          former_pnt--;
          current->split(bnode_list, edge_list, qface_list, gface_list);
      
          if(!mark_node(root_element, bnode_list.begin(), bnode_list.end())){
            cerr << "WARNING: fail to mark nodes, please check xml file" << endl;
            Loci::Abort();
          }
          for(int i = 0; i < current->numChildren(); i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(*split_mode_par, bnode_list, edge_list, qface_list, gface_list);

           
      newCellPlan[cc] =  aCell->make_cellplan();
   
    
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, qface_list);
      cleanup_list(gface_list);
      reduce_vector(newCellPlan[cc]);
   
    }
#endif
  }
};

register_rule<make_prism_cellplan_xml_norestart> register_make_prism_cellplan_xml_norestart;


class make_prism_cellplan_par:public pointwise_rule{
  const_param<std::string> parfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  
  const_store<bool> isIndivisible;
  const_param<int> split_mode_par;
  const_param<int> restart_par_par;
  store<std::vector<char> > newCellPlan;
  vector<source_par> sources;
  const_store<int>  node_l2f;
public:
  make_prism_cellplan_par(){
    name_store("parfile_par", parfile_par);
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("edgePlan", edgePlan);
  
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);
    name_store("restart_par_par", restart_par_par);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("restart_par_par");
    input("(split_mode_par,parfile_par)");

    input("(cellPlan, prism2face, prism2node, prismOrientCode)");
    input("isIndivisible");
    input("(lower, upper, boundary_map) ->(fileNumber(face2node), facePlan)"); 
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");

    
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      readPar(*parfile_par, sources);
      if(sources.size()==0){
        cerr << "WARNING: fail to read par file" << endl;
        Loci::Abort();
      }
      do_loop(seq, this);  
    }
 
  }
  void calculate(Entity cc){

    if(!isIndivisible[cc]){
   
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
      std::queue<Prism*> Q;

    
  

                                                 
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
                                      gface_list,
                                      node_l2f
                                      );
    
  
  
    
 
 
 
    
   
      std::vector<Prism*> cells;
      aCell->resplit( cellPlan[cc], 
                      bnode_list,
                      edge_list,
                      qface_list,
                      gface_list,
                      cells);
   
    
      Prism* current;
      int numCells = cells.size();
    
      if(numCells != 0){
        for(int i = 0; i < numCells; i++)Q.push(cells[i]);
      }
      else{
        Q.push(aCell);
      }

      while(!Q.empty()){
        current =Q.front();
        if(current->get_tagged(sources)){
          current->split(bnode_list, edge_list, qface_list, gface_list);
          for(int i = 0; i < current->numChildren(); i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(*split_mode_par, bnode_list, edge_list, qface_list, gface_list);
        
      //write new cellPlan
      newCellPlan[cc] = aCell->make_cellplan();
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, qface_list);
    
      cleanup_list(gface_list);
      reduce_vector(newCellPlan[cc]);
    
    
    }
  }
};

register_rule<make_prism_cellplan_par> register_make_prism_cellplan_par;

class make_prism_cellplan_par_norestart:public pointwise_rule{
  const_param<std::string> parfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<bool> isIndivisible;
  
  const_param<int> split_mode_par;
  const_param<int> norestart_par_par;
  store<std::vector<char> > newCellPlan;

  vector<source_par> sources;
  const_store<int> node_l2f;
public:
  make_prism_cellplan_par_norestart(){
    name_store("parfile_par", parfile_par);
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("prismOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
     
    name_store("isIndivisible", isIndivisible);
    name_store("split_mode_par", split_mode_par);
    name_store("norestart_par_par", norestart_par_par);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("norestart_par_par");
    input("split_mode_par, parfile_par"); 
    input("isIndivisible, prism2face, prism2node, prismOrientCode");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    output("newCellPlan");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      readPar(*parfile_par, sources);
      if(sources.size()==0){
        cerr << "WARNING: fail to read par file" << endl;
        Loci::Abort();
      }
      do_loop(seq, this); 
    }
  }
  void calculate(Entity cc){
    if(!isIndivisible[cc]){   
      std::list<Edge*> edge_list;
      std::list<QuadFace*> qface_list;
      std::list<Face*> gface_list;
      std::list<Node*> bnode_list;
      std::queue<Prism*> Q; 
    
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
                                      bnode_list,
                                      edge_list,
                                      qface_list,
                                      gface_list,
                                      node_l2f);
   
      Q.push(aCell);
  
      Prism* current;
      while(!Q.empty()){
        current =Q.front();
     
     
        if(current->get_tagged(sources)){
          current->split(bnode_list, edge_list, qface_list, gface_list);
          for(int i = 0; i < current->numChildren(); i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(*split_mode_par, bnode_list, edge_list, qface_list, gface_list);
   
   
      newCellPlan[cc] =  aCell->make_cellplan();
     
    
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, qface_list);
      cleanup_list(gface_list);
      reduce_vector(newCellPlan[cc]);
   
    }

  }
};

register_rule<make_prism_cellplan_par_norestart> register_make_prism_cellplan_par_norestart;


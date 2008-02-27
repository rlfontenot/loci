#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <Loci.h>
#include <algorithm>
#include "prism.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include "globals.h"

using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;
using Loci::storeRepP;
//int currentMem(void);
void mark_node( xmlNode* root_element,
               std::list<Node*>::iterator begin_pnt,
               std::list<Node*>::iterator end_pnt);

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
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
   const_store<char>  posTag;
  const_store<std::vector<char> > nodeTag;
  const_store<bool> isIndivisible;
  const_param<int> split_mode_par;
  const_param<int> restart_no_xml_par;
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
    name_store("restart_no_xml_par", restart_no_xml_par);
 
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("split_mode_par");
    input("(cellPlan, nodeTag, prism2face, prism2node, prismOrientCode)");
    input("isIndivisible");
    input("(lower, upper, boundary_map) -> (fileNumber(face2node),facePlan, nodeTag)"); 
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->(edgePlan, nodeTag)");
    input("restart_no_xml_par");
    
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
   
    if(numCells != 0){
      bool cell_split = false; //if any cell split, the whole big cell rebalance
      for(int i = 0; i < numCells; i++){
   
        cells[i]->setSplitCode(*split_mode_par);
        if(cells[i]->getMySplitCode()!=0) {
          if(*split_mode_par == 2){
            double min_edge_length =cells[i]->get_min_edge_length();
            int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));  
            cells[i]->resplit(min(Globals::levels,split_level),node_list, edge_list, qface_list, gface_list);
            cell_split = true;
          }
          else{
            cell_split = true;
            cells[i]->split(node_list, edge_list, qface_list, gface_list);
          }
        }
      }
      if(cell_split)aCell->rebalance_cells(*split_mode_par, node_list, edge_list, qface_list, gface_list);
      
    }
    else{
   
      aCell->setSplitCode(*split_mode_par);
      if(aCell->getMySplitCode() != 0 ){
        if(*split_mode_par == 2){
          double min_edge_length =aCell->get_min_edge_length();
          int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
          aCell->resplit(min(Globals::levels, split_level), node_list, edge_list, qface_list, gface_list); 
        }
        else{
          aCell->split(node_list, edge_list, qface_list, gface_list);
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
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<bool> isIndivisible;
    const_store<char>  posTag;
   const_param<int> split_mode_par;
  const_param<int> no_restart_no_xml_par;
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

    name_store("no_restart_no_xml_par", no_restart_no_xml_par);
    input("split_mode_par"); 
    input("isIndivisible, prism2face, prism2node, prismOrientCode");
    input("(lower, upper, boundary_map)->face2node->(pos, posTag)");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    input("no_restart_no_xml_par");
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


    aCell->setSplitCode(*split_mode_par);
  
  
      
    //write new cellPlan
    if((aCell->getMySplitCode()) != 0  && *split_mode_par == 2){
      double min_edge_length =aCell->get_min_edge_length();
      int split_level = int(log(min_edge_length/Globals::tolerance)/log(2.0));
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
//and nodeTag, posTag
class make_prism_cellplan_xml:public pointwise_rule{
   const_param<std::string> xmlfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<Array<char, 5> > prism2face;
  const_store<Array<char, 6> > prism2node;
  const_store<Array<char, 5> > orientCode;
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  
  const_store<bool> isIndivisible;
  const_param<int> split_mode_par;
  const_param<int> restart_xml_par;
  store<std::vector<char> > newCellPlan;
   xmlDoc* doc ;
  xmlNode* root_element ;

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
      doc = xmlReadFile((*xmlfile_par).c_str(), NULL, 0);
      root_element = xmlDocGetRootElement(doc);
      
      do_loop(seq, this);
      xmlFreeDoc(doc);
      xmlCleanupParser();
      xmlMemoryDump();
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
    
  
  
  mark_node(root_element, bnode_list.begin(), bnode_list.end());
 
 
 
    
   
    std::vector<Prism*> cells;
    std::list<Node*>::iterator former_pnt = bnode_list.end(); 
    aCell->resplit( cellPlan[cc], 
                    bnode_list,
                    edge_list,
                    qface_list,
                    gface_list,
                    cells);
    former_pnt--;
    mark_node(root_element, bnode_list.begin(), bnode_list.end());
    
    
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
        

      
      
     if( current->getMySplitCode()==0)current->setSplitCode(*split_mode_par);
      if(current->getMySplitCode()!=0){
        std::list<Node*>::iterator former_pnt = bnode_list.end();
        current->split(bnode_list, edge_list, qface_list, gface_list);
        

        former_pnt--;
        mark_node(root_element, bnode_list.begin(), bnode_list.end());
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
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<bool> isIndivisible;
   
   const_param<int> split_mode_par;
  const_param<int> no_restart_xml_par;
  store<std::vector<char> > newCellPlan;
  xmlDoc* doc ;
  xmlNode* root_element ;

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
    name_store("no_restart_xml_par", no_restart_xml_par);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(face2node)", node_l2f);
    input("no_restart_xml_par");
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
      doc = xmlReadFile((*xmlfile_par).c_str(), NULL, 0);
      root_element = xmlDocGetRootElement(doc);
  
      do_loop(seq, this);
      xmlFreeDoc(doc);
      xmlCleanupParser();
      xmlMemoryDump();
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
    mark_node(root_element, bnode_list.begin(), bnode_list.end());
    Q.push(aCell);
  
   Prism* current;
   while(!Q.empty()){
     current =Q.front();
     
     
     if(current->getMySplitCode()==0) current->setSplitCode(*split_mode_par);
     if(current->getMySplitCode()!=0){
       std::list<Node*>::iterator former_pnt = bnode_list.end();
       former_pnt--;
       current->split(bnode_list, edge_list, qface_list, gface_list);
      
       mark_node(root_element, bnode_list.begin(), bnode_list.end());
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

register_rule<make_prism_cellplan_xml_norestart> register_make_prism_cellplan_xml_norestart;


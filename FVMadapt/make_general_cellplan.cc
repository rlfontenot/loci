#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <Loci.h>
#include <algorithm>
#include "diamondcell.h"
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
class make_general_cellplan:public pointwise_rule{
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<bool> is_quadface;
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<char>  posTag;
  const_store<std::vector<char> > nodeTag;
  const_store<bool> isIndivisible;
  const_param<int> restart_no_xml_par;
  store<std::vector<char> > newCellPlan;

  const_store<int> node_l2f;
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
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(pos)", node_l2f);
    name_store("is_quadface", is_quadface);
    name_store("restart_no_xml_par", restart_no_xml_par);
    input("restart_no_xml_par");
    input("(cellPlan,nodeTag) ");
    input("isIndivisible");
    input("(lower, upper, boundary_map) -> (facePlan,is_quadface, nodeTag)"); 
    input("(lower, upper, boundary_map)->face2node->(pos, posTag, fileNumber(pos))");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->(edgePlan, nodeTag)");

    output("newCellPlan");
    constraint("gnrlcells");
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
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<char>  posTag;
  const_store<bool> isIndivisible;
  const_param<int> no_restart_no_xml_par;
  store<std::vector<char> > newCellPlan;

  const_store<int> node_l2f;
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
    name_store("no_restart_no_xml_par", no_restart_no_xml_par);
    name_store("isIndivisible", isIndivisible);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(pos)", node_l2f);
    
  
    input("isIndivisible");
    input("(lower, upper, boundary_map)->face2node->(pos, posTag, fileNumber(pos))");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("no_restart_no_xml_par");
    output("newCellPlan");
    constraint("gnrlcells");
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



class make_general_cellplan_xml:public pointwise_rule{
  const_param<std::string> xmlfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_store<bool> is_quadface;
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<bool> isIndivisible;
  const_param<int> restart_xml_par;
  store<std::vector<char> > newCellPlan;

  const_store<int> node_l2f;
  xmlDoc* doc ;
  xmlNode* root_element ;
public:
  make_general_cellplan_xml(){
    name_store("xmlfile_par", xmlfile_par);
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
    name_store("restart_xml_par", restart_xml_par);
    name_store("isIndivisible", isIndivisible);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(pos)", node_l2f);
    name_store("is_quadface", is_quadface);
    input("xmlfile_par");
    input("restart_xml_par");
    input("cellPlan ");
    input("isIndivisible");
    input("(lower, upper, boundary_map) -> (facePlan,is_quadface)"); 
    input("(lower, upper, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edgePlan");
 
    output("newCellPlan");
    constraint("gnrlcells");
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
      std::list<Face*> face_list;
      std::list<Node*> bnode_list;
      std::queue<DiamondCell*> Q;

    
 

                                                 
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
                                       bnode_list,
                                       edge_list,
                                       face_list,
                                       node_l2f);

  
  
 


      
      
      
      mark_node(root_element, bnode_list.begin(), bnode_list.end());
     
      
      std::vector<DiamondCell*> cells;

      std::list<Node*>::iterator former_pnt = bnode_list.end();
      aCell->resplit( cellPlan[cc], 
                      bnode_list,
                      edge_list,
                      face_list,
                      cells);
     
      former_pnt--;
      mark_node(root_element, bnode_list.begin(), bnode_list.end());

      
     
      DiamondCell* current;
    
      // split the cell Globals::levels times
      int numCells = cells.size();

      if(numCells != 0){
        for(int i = 0; i < numCells; i++)Q.push(cells[i]);
      }
      else{
        double min_edge_length =aCell->get_min_edge_length();
        if(aCell->get_tagged() && min_edge_length > (2*Globals::tolerance)){
         

          former_pnt = bnode_list.end();
          aCell->split(bnode_list, edge_list, face_list);
          former_pnt--;
          mark_node(root_element, bnode_list.begin(), bnode_list.end());
         
          
          for(int i = 0; i < aCell->numNode; i++){
            Q.push(aCell->child[i]);
          }
        }
      }
      while(!Q.empty()){
        current =Q.front();
        //calculate min_edge_length of the cell
        double min_edge_length =current->get_min_edge_length();
        
        if(current->get_tagged() && min_edge_length > (2*Globals::tolerance)){

          former_pnt = bnode_list.end();
          current->split(bnode_list, edge_list, face_list);
          
          former_pnt--;
          mark_node(root_element, bnode_list.begin(), bnode_list.end());

          for(int i = 0; i < 2*current->getNfold()+2; i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(bnode_list, edge_list, face_list);
      
      
      //write new cellPlan
      newCellPlan[cc] = aCell->make_cellplan();
      //clean up
      if(aCell != 0){
        delete aCell;
        aCell = 0;
      }
      cleanup_list(bnode_list, edge_list, face_list);
     
      reduce_vector(newCellPlan[cc]);
    }
  }
};

register_rule<make_general_cellplan_xml> register_make_general_cellplan_xml;

class make_general_cellplan_xml_norestart:public pointwise_rule{
 const_param<std::string> xmlfile_par;
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_multiMap edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
   const_store<bool> isIndivisible;
  const_param<int> no_restart_xml_par;
  store<std::vector<char> > newCellPlan;

  const_store<int> node_l2f;
  xmlDoc* doc ;
  xmlNode* root_element ;
public:
  make_general_cellplan_xml_norestart(){
    name_store("xmlfile_par", xmlfile_par);
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("no_restart_xml_par", no_restart_xml_par);
   
    name_store("isIndivisible", isIndivisible);
    name_store("newCellPlan", newCellPlan);
    name_store("fileNumber(pos)", node_l2f);
    
    input("no_restart_xml_par");
    input("(isIndivisible,xmlfile_par)");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->face2node->fileNumber(pos)");
    output("newCellPlan");
    constraint("gnrlcells");
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
    std::list<Face*> face_list;
    std::list<Node*> bnode_list;
    std::queue<DiamondCell*> Q;
    
    Cell* aCell = build_general_cell(lower[cc].begin(), lower.num_elems(cc),
                                     upper[cc].begin(), upper.num_elems(cc),
                                     boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                     face2node,
                                     face2edge,
                                     edge2node,
                                     pos,
                                     bnode_list,
                                     edge_list,
                                     face_list,
                                     node_l2f);


   
    
    DiamondCell* current;
    
    mark_node(root_element,bnode_list.begin(), bnode_list.end());
    double min_edge_length =aCell->get_min_edge_length();


    if(aCell->get_tagged() && min_edge_length > (2*Globals::tolerance)){
       std::list<Node*>::iterator former_pnt = bnode_list.end();
       aCell->split(bnode_list, edge_list, face_list);
       former_pnt--;
       mark_node(root_element, bnode_list.begin(), bnode_list.end());
       for(int i = 0; i < aCell->numNode; i++){
         Q.push(aCell->child[i]);
       }
      
      while(!Q.empty()){
        current =Q.front();
        //calculate min_edge_length of the cell
        double min_edge_length =current->get_min_edge_length();
      
        if(current->get_tagged() && min_edge_length > (2*Globals::tolerance)){
          std::list<Node*>::iterator former_pnt = bnode_list.end();
          current->split(bnode_list, edge_list, face_list);
          
          former_pnt--;
          mark_node(root_element, bnode_list.begin(), bnode_list.end());
          for(int i = 0; i <2*current->getNfold()+2; i++){
            Q.push(current->getChildCell(i));
          }
        }
        Q.pop();
      }
      aCell->rebalance_cells(bnode_list, edge_list, face_list);
     
    }
    //write new cellPlan
    newCellPlan[cc] = aCell->make_cellplan();
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    cleanup_list(bnode_list, edge_list, face_list);
  
    reduce_vector(newCellPlan[cc]);
  
    }
  }
};

register_rule<make_general_cellplan_xml_norestart> register_make_general_cellplan_xml_norestart;


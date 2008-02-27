#include <iostream>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <vector>
#include <stack>
#include <Loci.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include "node_edge.h"

using std::stack;
using std::vector;
using std::cerr;
using std::endl;


// struct vect3d{
//   double x;
//   double y;
//   double z;
//   vect3d(double x0, double y0, double z0):x(x0), y(y0), z(z0){};
//   vect3d():x(0.0), y(0.0), z(0.0){};
// };
//code from Loci::vogmerge.cc
struct affineMapping {
  double M[4][4] ;
  //unit matrix
  affineMapping() {
    for(int i=0;i<4;++i) {
      for(int j=0;j<4;++j)
        M[i][j] = 0 ;
      M[i][i] = 1 ;
    }
  }
  
  void Combine(affineMapping a) {
    double Mtmp[4][4] ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) {
        double mtmp = 0 ;
        for(int k=0;k<4;++k)
          mtmp += a.M[i][k]*M[k][j] ;
        Mtmp[i][j] = mtmp ;
      }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        M[i][j] = Mtmp[i][j] ;
  }
  void translate(vect3d tv) {
    affineMapping tmp ;
    tmp.M[0][3] = tv.x ;
    tmp.M[1][3] = tv.y ;
    tmp.M[2][3] = tv.z ;
    Combine(tmp) ;
  }
  void scale(vect3d tv) {
    affineMapping tmp ;
    tmp.M[0][0] = tv.x ;
    tmp.M[1][1] = tv.y ;
    tmp.M[2][2] = tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  sth ;
    tmp.M[2][1] = -sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][2] = -sth ;
    tmp.M[2][0] =  sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][1] =  sth ;
    tmp.M[1][0] = -sth ;
    tmp.M[1][1] =  cth ;
    Combine(tmp) ;
  }
  vect3d Mapping(vect3d v) {
    double tmp[4] ;
    tmp[0] = v.x ;
    tmp[1] = v.y ;
    tmp[2] = v.z ;
    tmp[3] = 1. ;
    double res[4] ;
    for(int i=0;i<4;++i)
      res[i] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        res[i] += M[i][j]*tmp[j] ;
    vect3d r(res[0],res[1],res[2]) ;
    return r ;
  }
   
} ;





 std::vector<bool> process_shape(xmlNode* anode,  const std::vector<vect3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  
  xmlChar* the_shape = anode->children->content;
  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    vect3d p = pointSet[i];
    if(xmlStrEqual(the_shape, BAD_CAST("sphere")) ){
      result[i] = ((p.x*p.x + p.y*p.y + p.z*p.z) <= 1.0);
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("cone")) ){
      result[i] =  ((p.x*p.x + p.y*p.y ) <= p.z*p.z);
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("cylinder")) ){
      result[i] = ( (p.x*p.x + p.y*p.y ) <= 1.0);
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("box")) ){
      result[i] = (( p.x >= 0) && (p.x <= 1.0)
                   &&(p.y >= 0) && (p.y <= 1.0)
                   &&(p.z >=0) && (p.z <=1.0));
    } 
    
    else if(xmlStrEqual(the_shape, BAD_CAST("x+plane"))){
      result[i] = ( p.x >= 0) ;
    }
    
    else if(xmlStrEqual(the_shape, BAD_CAST("x-plane"))){
      result[i] = ( p.x < 0) ;
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("y+plane"))){
      result[i] =( p.y >= 0) ;
   }
    else if(xmlStrEqual(the_shape, BAD_CAST("y-plane"))){
      result[i] =( p.y < 0) ;
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("z+plane"))){
      result[i] = ( p.z >= 0) ;
    }
    else if(xmlStrEqual(the_shape, BAD_CAST("z-plane"))){
      result[i] =  (p.z < 0) ;
    }

  

    else{
      cerr << "WARNING: don't know the shape" << the_shape << endl;
     
    
    }
  }
  return result;
}

 vect3d  process_translate(xmlNode* anode){
  xmlNode* children = anode->children;
  xmlNode* cur_node =NULL;
  
  double x0, y0, z0;
  
  for(cur_node = children; cur_node; cur_node = cur_node->next){
    if(cur_node->type == 1){
      if(xmlStrEqual(cur_node->name, BAD_CAST("x0")) ){ 
        x0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
      else  if(xmlStrEqual(cur_node->name, BAD_CAST("y0")) ){ 
        y0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
      else  if(xmlStrEqual(cur_node->name, BAD_CAST("z0")) ){ 
        z0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
    }
    
  }
  
  return vect3d(x0, y0, z0);
  
}


 vect3d  process_scale(xmlNode* anode){
  xmlNode* children = anode->children;
  xmlNode* cur_node =NULL;
  
  double x0, y0, z0;
  
  for(cur_node = children; cur_node; cur_node = cur_node->next){
    if(cur_node->type == 1){
      if(xmlStrEqual(cur_node->name, BAD_CAST("x0")) ){ 
        x0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
      else  if(xmlStrEqual(cur_node->name, BAD_CAST("y0")) ){ 
        y0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
      else  if(xmlStrEqual(cur_node->name, BAD_CAST("z0")) ){ 
        z0 = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
    }
    
  }
  return vect3d(x0, y0, z0); 
  
}


 double  process_rotate(xmlNode* anode){
  xmlNode* children = anode->children;
  xmlNode* cur_node =NULL;
  
  double theta = 0.0;
  
  for(cur_node = children; cur_node; cur_node = cur_node->next){
    if(cur_node->type == 1){
      if(xmlStrEqual(cur_node->name, BAD_CAST("theta")) ){ 
        theta = xmlXPathCastNodeToNumber(cur_node->children);
       
      }
     
    }
    
  }
  return theta; 
}





 void  process_transform(xmlNode* anode,  std::vector<vect3d>& p){
  xmlNode* children = anode->children;
  xmlNode* cur_node =NULL;
  affineMapping aMatrix = affineMapping();
  
  for(cur_node = children; cur_node; cur_node = cur_node->next){
    if(cur_node->type == 1){
      if(xmlStrEqual(cur_node->name, BAD_CAST("translate")) ){ 
        aMatrix.translate( process_translate(cur_node));
      }
     else if(xmlStrEqual(cur_node->name, BAD_CAST("scale")) ){ 
        aMatrix.scale( process_scale(cur_node));
      }
      
      else  if(xmlStrEqual(cur_node->name, BAD_CAST("rotateX")) ){ 
        aMatrix.rotateX( process_rotate(cur_node));
      }

      else  if(xmlStrEqual(cur_node->name, BAD_CAST("rotateY")) ){ 
        aMatrix.rotateY( process_rotate(cur_node));
      }

      else  if(xmlStrEqual(cur_node->name, BAD_CAST("rotateZ")) ){ 
        aMatrix.rotateZ( process_rotate(cur_node));
      } 
    }
  }

  for(unsigned int i = 0; i< p.size(); i++){
    p[i] =   aMatrix.Mapping(p[i]);
  }
  
}

std::vector<bool>  process_object(xmlNode* anode,  std::vector<vect3d> p){
  
  xmlNode* children = anode->children;
  xmlNode* cur_node =NULL;
  for(cur_node = children; cur_node; cur_node = cur_node->next){
    if(cur_node->type == 1){
     
      if(xmlStrEqual(cur_node->name, BAD_CAST("transform")) ){
        process_transform(cur_node, p);
      }
      else if(xmlStrEqual(cur_node->name, BAD_CAST "shape") ){
        return process_shape(cur_node, p);
      }
    }
  }
  return std::vector<bool>();
}



 vector<bool> process_stack(stack<vector<bool> >& object_value, stack<char>& op_value){
  if(object_value.empty()){
    cerr << "WARNING: no object in process_stack()" << endl;
    return vector<bool>();
  }
  if(op_value.empty() && object_value.size() > 1){
    cerr << "WARNING: more than 1 object in process_stack()" << endl;
    return object_value.top();
  }
  
  if(op_value.empty()) return object_value.top(); // only one object

  
  char op  ;
  vector<bool> value1, value2 ;
  while(!op_value.empty()){
    op = op_value.top();
    op_value.pop();

    if(op == '-'){
      value1 = object_value.top();
      object_value.pop();
      for(unsigned int i = 0; i<value1.size(); i++)value1[i] = !value1[i]; 
      object_value.push(value1);
    }
    else if(op == '&'){
      
      value1 = object_value.top();
      object_value.pop();
      value2 = object_value.top();
      object_value.pop();

      for(unsigned int i = 0; i<value1.size(); i++)value1[i] = value1[i] && value2[i]; 
      object_value.push(value1);
    }
    else if(op == '|'){
      
      value1 = object_value.top();
      object_value.pop();
      value2 = object_value.top();
      object_value.pop();
       for(unsigned int i = 0; i<value1.size(); i++)value1[i] = value1[i] || value2[i]; 
      object_value.push(value1);
    }
    else{
      cerr << "WARNING:illegal op in process_stack()" << endl;

    }
  }

  if(object_value.size() > 1) cerr << "WARNING: too many objects in process_stack()"<< endl;
  return object_value.top();
}





 std::vector<bool> process_region(xmlNode* anode, std::vector<vect3d> p)
{
  stack<vector<bool> > object_value;
  stack<char> op_value;
 
  xmlNode* children = anode->children;
  
  
  xmlNode* cur_node =NULL;
   
  for(cur_node = children; cur_node; cur_node = cur_node->next){
   
    if(cur_node->type == 1){
 
      if(xmlStrEqual(cur_node->name, BAD_CAST "object") ){
        //make a copy of p, then process it
       
        object_value.push(process_object(cur_node, p));
      }
      else if(xmlStrEqual(cur_node->name, BAD_CAST "op") ){
       //don't use xmlNodeGetContent, otherwise have to free memory later   
        xmlChar* content = cur_node->children->content;
        if(xmlStrEqual(content, BAD_CAST "and"))op_value.push('&');
        else if(xmlStrEqual(content, BAD_CAST "or"))op_value.push('|');
        else if(xmlStrEqual(content, BAD_CAST "minus"))op_value.push('-');
        else
          {
            cerr << "WARNING: cann't recognize op" << content << endl;  
          }
        
      }
      //recursive function allow multi-level tree
      else if(xmlStrEqual(cur_node->name, BAD_CAST "region") ){
        //p is a copy
        object_value.push(process_region(cur_node, p));
      }
        
      
      else{
        cerr << "WARNING: cann't recognize the children of region" << cur_node->name <<endl;
      }
    }
  }

  
  return process_stack(object_value, op_value);
}
    




void mark_node( xmlNode* root_element, std::list<Node*>::iterator begin_pnt, std::list<Node*>::iterator end_pnt){
  std::list<Node*>::iterator current_pnt = begin_pnt;
  std::vector<vect3d> pointSet;
  
  for(current_pnt = begin_pnt; current_pnt != end_pnt; current_pnt++){
    pointSet.push_back((*current_pnt)->p);
  }
  std::vector<bool> result = process_region(root_element, pointSet);
  int index = 0;
  for(current_pnt = begin_pnt; current_pnt != end_pnt; current_pnt++, index++){
    if(result[index])(*current_pnt)->tag = 1;
    else (*current_pnt)->tag = 0;
  }
}

// int main(int argc, char** argv){
//   xmlDoc*doc = NULL;
//   xmlNode* root_element = NULL;
//   if(argc != 2) return 1;
  
//   LIBXML_TEST_VERSION;
//   doc = xmlReadFile(argv[1], NULL, 0);
//   if(doc == NULL){
//     cerr << " error: could not parse file " << argv[1] << endl;
//   }
  
//   root_element = xmlDocGetRootElement(doc);
//   cerr << root_element->name << endl;
//   std::vector<vect3d> pointSet;
//   pointSet.push_back(vect3d(1.0, 0.1, 0.2));
//   pointSet.push_back(vect3d(2.0, 0.1, 0.2));
//   vector<bool> result = process_region(root_element, pointSet);
//   cerr << "result: "<<result[0] << result[2] << endl;   
//   xmlFreeDoc(doc);
//   xmlCleanupParser();
//   xmlMemoryDump();
//   return 0;
// }
  

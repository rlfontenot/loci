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
#include <iostream>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <vector>
#include <stack>
#include "grid.h"
#include <QDomNode>


using std::stack;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;




 std::vector<bool> process_sphere(const QDomElement& anode,  const std::vector<positions3d>& pointSet){

   std::vector<bool> result(pointSet.size());
   if(anode.isNull())return result;

   double x0=0.0, y0=0.0, z0=0.0, r = 1.0;
  
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
     if(cur_node.tagName() == "x0"){
      if(cur_node.text()!="") x0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "y0"){
       if(cur_node.text()!="")  y0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "z0"){
        if(cur_node.text()!="") z0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "r"){
       if(cur_node.text()!="")  r = cur_node.text().toDouble();
       
     }else{
       cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
     }
   } 
     
   for(unsigned int i = 0; i < pointSet.size(); i++){
     positions3d p = pointSet[i];
     result[i] = (((p.x-x0)*(p.x-x0) + (p.y-y0)*(p.y-y0) + (p.z-z0)*(p.z-z0)) <= r*r);
    
   }
  
   return result;
 }


std::vector<bool> process_cone(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;

   double x0=0.0, y0=0.0, z0=0.0, r = 1.0, z1 = 0.0, z2 = 1.0;
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
     if(cur_node.tagName() == "x0"){
       if(cur_node.text()!="") x0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "y0"){
       if(cur_node.text()!="")  y0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "z0"){
        if(cur_node.text()!="") z0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "r"){
       if(cur_node.text()!="")  r = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "z1"){
       if(cur_node.text()!="")  z1 = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "z2"){
       if(cur_node.text()!="")  z2 = cur_node.text().toDouble();
       
     }else{
       cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
     }
   } 
     
   for(unsigned int i = 0; i < pointSet.size(); i++){
     positions3d p = pointSet[i];
     result[i] = (((p.x-x0)*(p.x-x0) + (p.y-y0)*(p.y-y0)) <= (r*r*(p.z-z0)*(p.z-z0)))&& (p.z >= z1) && (p.z <= z2) ;
     
   }
   return result;
}

std::vector<bool> process_cylinder(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  
  
  
  double x0=0.0, y0=0.0, r = 1.0, z1 = 0.0, z2 = 1.0;
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
     if(cur_node.tagName() == "x0"){
       if(cur_node.text()!="") x0 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "y0"){
       if(cur_node.text()!="")  y0 = cur_node.text().toDouble();
                        
     }else if(cur_node.tagName() == "r"){
       if(cur_node.text()!="")  r = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "z1"){
       if(cur_node.text()!="")  z1 = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "z2"){
       if(cur_node.text()!="")  z2 = cur_node.text().toDouble();
       
     }else{
       cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
     }
   } 
   
   
   for(unsigned int i = 0; i < pointSet.size(); i++){
     positions3d p = pointSet[i];
     result[i] = (((p.x-x0)*(p.x-x0) + (p.y-y0)*(p.y-y0)) <= r*r)&& (p.z >= z1) && (p.z <= z2) ;
     
   }
   return result;
}

std::vector<bool> process_box(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  
    
  double x1=0.0, x2=1.0, y1 = 0.0,y2 = 1.0, z1 = 0.0, z2 = 1.0;
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
     if(cur_node.tagName() == "x1"){
       if(cur_node.text()!="") x1 = cur_node.text().toDouble();
                
     }else if(cur_node.tagName() == "y1"){
       if(cur_node.text()!="")  y1 = cur_node.text().toDouble();
                        
     }else if(cur_node.tagName() == "z1"){
       if(cur_node.text()!="")  z1 = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "x2"){
       if(cur_node.text()!="")  x2 = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "y2"){
       if(cur_node.text()!="")  y2 = cur_node.text().toDouble();
       
     }else if(cur_node.tagName() == "z2"){
       if(cur_node.text()!="")  z2 = cur_node.text().toDouble();
       
     }else{
       cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
     }
   } 
   
    
   
   if(x1>x2 || y1>y2 ||z1>z2){
     cerr<<"WARNING: xml file, illegal value in the definition of box:" << endl;
     cerr <<"        x >= x1 && x <= x2 && y >= y1 && y <= y2 && z >= z1 && z <= z2 "<< endl;
     cerr <<"        no point is marked!" << endl;
     return result;
   }  
   for(unsigned int i = 0; i < pointSet.size(); i++){
     positions3d p = pointSet[i];
     result[i] = (p.x >= x1) && (p.x <= x2) && (p.y >= y1) && (p.y <= y2) && (p.z >= z1) && (p.z <=z2);
   }
   return result;
}

std::vector<bool> process_x_plus_plane(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  
  double x1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "x1"){
      if(cur_node.text()!="") x1 = cur_node.text().toDouble();
       
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 


  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.x >= x1);
    
  }
  return result;
}

std::vector<bool> process_x_minus_plane(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());

 if(anode.isNull())return result;
  
  double x1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "x1"){
      if(cur_node.text()!="") x1 = cur_node.text().toDouble();
       
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 



  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.x <= x1);
      
  }
  return result;
}

std::vector<bool> process_y_plus_plane(const QDomElement& anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  double y1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "y1"){
      if(cur_node.text()!="") y1 = cur_node.text().toDouble();
       
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 



  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.y >= y1);
      
  }
  return result;
}

std::vector<bool> process_y_minus_plane(const QDomElement anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  
  if(anode.isNull())return result;
  double y1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "y1"){
      if(cur_node.text()!="") y1 = cur_node.text().toDouble();
      
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 

  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.y <= y1);
    
  }
  return result;
}


std::vector<bool> process_z_plus_plane(const QDomElement anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  double z1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "z1"){
      if(cur_node.text()!="") z1 = cur_node.text().toDouble();
      
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 
  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.z >= z1);
      
  }
  return result;
}

std::vector<bool> process_z_minus_plane(const QDomElement anode,  const std::vector<positions3d>& pointSet){
  std::vector<bool> result(pointSet.size());
  if(anode.isNull())return result;
  double z1=0.0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "z1"){
      if(cur_node.text()!="") z1 = cur_node.text().toDouble();
      
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 
  
  
  for(unsigned int i = 0; i < pointSet.size(); i++){
    positions3d p = pointSet[i];
    result[i] = (p.z <= z1);
      
  }
  return result;
}


positions3d  process_translate(const QDomElement& anode){
   
  double x0=0, y0=0, z0=0;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "x0"){
      if(cur_node.text()!="") x0 = cur_node.text().toDouble();
      
    }else if(cur_node.tagName() == "y0"){
      if(cur_node.text()!="")  y0 = cur_node.text().toDouble();
      
    }else if(cur_node.tagName() == "z0"){
        if(cur_node.text()!="") z0 = cur_node.text().toDouble();
        
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 
       
 
  
  return positions3d(x0, y0, z0);
  
}

positions3d  process_scale(const QDomElement& anode){
   
  double x0=1, y0=1, z0=1;
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "x0"){
      if(cur_node.text()!="") x0 = cur_node.text().toDouble();
      
    }else if(cur_node.tagName() == "y0"){
      if(cur_node.text()!="")  y0 = cur_node.text().toDouble();
      
    }else if(cur_node.tagName() == "z0"){
        if(cur_node.text()!="") z0 = cur_node.text().toDouble();
        
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
  } 
       
 
  cerr << "scale: " << x0 << " " << y0 << " " << z0<< endl;  
  return positions3d(x0, y0, z0);
  
}
 


 double  process_rotate(const QDomElement& anode){
  
  
  double theta = 0.0;
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    if(cur_node.tagName() == "theta"){
      if(cur_node.text()!="") theta = cur_node.text().toDouble();
      
    }else{
      cerr<<"can not parse "<<cur_node.tagName().toStdString()<< endl;
    }
   } 
       
  
   return theta; 
 }





 void  process_transform(const QDomElement& anode,  std::vector<positions3d>& p){
   
   affineMapping aMatrix = affineMapping();
   for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){ 
     
     if(cur_node.tagName() == "translate"){ 
        aMatrix.translate( process_translate(cur_node));
     }else if(cur_node.tagName()=="scale"){ 
        aMatrix.scale( process_scale(cur_node));
      }else  if(cur_node.tagName() =="rotateX"){ 
        aMatrix.rotateX( process_rotate(cur_node));
      }else  if(cur_node.tagName()=="rotateY"){ 
        aMatrix.rotateY( process_rotate(cur_node));
      }else  if(cur_node.tagName()=="rotateZ"){ 
        aMatrix.rotateZ( process_rotate(cur_node));
      }else{
        cerr <<"can not parse "<<cur_node.tagName().toStdString()<< " in the children elements of 'transform'"<<endl;
        cerr <<"         the children elements can be: 'translate', 'scale', 'rotateX', 'rotateY' or 'rotateZ'" << endl;
        
      }
   }
 

   for(unsigned int i = 0; i< p.size(); i++){
     if(i == 1) cerr << "before: " << p[i].x << "  " << p[i].y << " " << p[i].z << endl;    
     p[i] =   aMatrix.MapNode(p[i]);
       if(i == 1) cerr << "after: " << p[i].x << "  " << p[i].y << " " << p[i].z << endl;    
   }
  
 }

std::vector<bool> process_shape(const QDomElement& anode, const std::vector<positions3d>& p){
  
  
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
    
    if(cur_node.tagName() =="sphere"){ 
      return process_sphere(cur_node, p);
    }else if(cur_node.tagName() =="cone"){ 
      return process_cone(cur_node, p);
    }else  if(cur_node.tagName() =="cylinder"){ 
      return process_cylinder(cur_node, p);
    }else  if(cur_node.tagName() =="box"){ 
      return process_box(cur_node, p);
    }else  if(cur_node.tagName() =="x_plus_plane"){ 
      return process_x_plus_plane(cur_node, p);
    }else  if(cur_node.tagName() =="x_minus_plane"){ 
        return process_x_minus_plane(cur_node, p);
    }else  if(cur_node.tagName() =="y_plus_plane"){ 
      return process_y_plus_plane(cur_node, p);
    }else  if(cur_node.tagName() =="y_minus_plane"){ 
      return process_y_minus_plane(cur_node, p);
    }else  if(cur_node.tagName() =="z_plus_plane"){ 
      return process_z_plus_plane(cur_node, p);
    } else  if(cur_node.tagName() =="z_minus_plane"){ 
      return process_z_minus_plane(cur_node, p);
    } else{
      cerr <<"WARNING: unknown element name: "<<cur_node.tagName().toStdString()<< " in the children elements of 'shape'"<<endl;
      cerr <<"         the children elements can be: 'sphere', 'cone', 'cylinder', 'box',"<<endl;
      cerr <<"                                      'x_plus_plane', 'x_minus_plane', 'y_plus_plane'"<<endl;
      cerr <<"                                      'y_minus_plane', 'z_plus_plane' or 'z_minus_plane'" << endl;
      
    }
  }

  return vector<bool>() ;
}




std::vector<bool>  process_object(const QDomElement& anode,  std::vector<positions3d> p){
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){     
    if(cur_node.tagName() =="transform"){
      process_transform(cur_node, p);
    } else if(cur_node.tagName() == "shape"){
      return process_shape(cur_node, p);
    }else{
      cerr <<"WARNING: unknown element name: "<<cur_node.tagName().toStdString()<< " in the children elements of 'object'"<<endl;
      cerr <<"         the children elements can be: 'transform' or 'shape'" << endl;
    
    }
  }
  return std::vector<bool>();
}



 vector<bool> process_stack(stack<vector<bool> >& object_value, stack<char>& op_value){
   if(object_value.empty()){
    cerr << "WARNING: no object in stack" << endl;
    return std::vector<bool>(); 
  }
  if(op_value.empty() && object_value.size() > 1){
    cerr << "WARNING: more than one objects left  in stack" << endl;
    return std::vector<bool>();  
  }
  
  if(op_value.empty()) return object_value.top(); // only one object
  
 
  char op  ;
  vector<bool> value1, value2 ;
  while(!op_value.empty()){
    op = op_value.top();
    op_value.pop();

    if(op == '!'){
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

    else if(op == '-'){
      
      value1 = object_value.top();
      object_value.pop();
      value2 = object_value.top();
      object_value.pop();
       for(unsigned int i = 0; i<value1.size(); i++)value1[i] = value1[i] &&(!value2[i]); 
      object_value.push(value1);
    } 
    else{
      cerr << "WARNING:illegal op in process_stack()" << endl;
      return std::vector<bool>();   
    }
  }

  if(object_value.size() > 1){
    cerr << "WARNING: too many objects in process_stack()"<< endl;
    return std::vector<bool>();  
  }
  return object_value.top();
 }





 std::vector<bool> process_region(const QDomElement& anode, std::vector<positions3d> p)
 {
   
  
  stack<vector<bool> > object_value;
  stack<char> op_value;
    
  for(QDomElement cur_node = anode.firstChildElement(); !cur_node.isNull(); cur_node = cur_node.nextSiblingElement()){
  
    if(cur_node.tagName()== "object"){
      //make a copy of p, then process it
      object_value.push(process_object(cur_node, p));
    }else if(cur_node.tagName()== "op"){
  
      if(cur_node.text()== "intersection")op_value.push('&');
      else if(cur_node.text()== "union")op_value.push('|');
      else if(cur_node.text()== "difference")op_value.push('-');
      else if(cur_node.text()== "complement")op_value.push('!');
      else{
        cerr <<"WARNING: unknown element name: "<<cur_node.tagName().toStdString()<< " in the children elements of 'op'"<<endl;
        cerr <<"         the children elements can be: 'intersection', 'union', 'difference', 'complement'"<<endl;
        return vector<bool>();
      }
        
    }else if(cur_node.tagName()== "region"){ //recursive function allow multi-level tree
    //p is a copy
      object_value.push(process_region(cur_node, p));
    }else{
      cerr <<"WARNING: unknown element name: "<<cur_node.tagName().toStdString()<< " in the children elements of 'region'"<<endl;
      cerr <<"         the children elements can be: 'object', 'op' or  'region'"<<endl;
      return vector<bool>();
    }
  }
  return process_stack(object_value, op_value);
}

    





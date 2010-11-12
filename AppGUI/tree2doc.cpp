#include "tree2doc.h"
#include <math.h>
#include "defines.h"
//#include "pboperation.h"
using std::vector;

#define PI 3.14159265358979323846264338327950


//order of rotation y first(heading), then z(atitude), then x(bank)
bool anglesBetween(positions3d v1,positions3d v2, double& heading, double& attitude, double& bank) {
  double norm1 = norm(v1);
  double norm2 = norm(v2);
  double angle = 0;
  positions3d axis = positions3d(0, 0, 1);
  
  if(norm1 > 1e-34 && norm2 > 1e-34){
        
    // turn vectors into unit vectors 
    positions3d n1 = v1/norm1;
    positions3d n2 = v2/norm2; 
    angle = acos(dot(n1,n2)); 
	
    if( fabs(angle) > 1e-5 && fabs(angle-PI)>1e-5){
      axis = cross(n1, n2);
    }else if(fabs(angle) <= 1e-5){
      // if no noticable rotation is available return zero rotation
      // this way we avoid Cross product artifacts 
      axis = positions3d(0, 0, 1);
      angle = 0;
    }else if(fabs(angle-PI) < 1e-5){
      // in this case there are 2 lines on the same axis 
      // there are an infinite number of normals 
      // in this case. Anyone of these normals will be 
      // a valid rotation (180 degrees). so I pick one axis that is perpendicular to n1
      angle = PI;
      if(fabs(fabs(n1.z)-1.0)<1e-5)axis=positions3d(1, 0, 0);
      else axis = positions3d(-1*n1.y, n1.x, 0);
    }
  }else{
    return false;
  }
  
  //from axis angle to Euler
  double s=sin(angle);
  double c=cos(angle);
  double t=1-c;
  //  if axis is not already normalized then uncomment this
  double magnitude = norm(axis);
  if (magnitude==0){
    cerr << " maginitude is zero" << endl;
    return false;
  }
  axis = axis/magnitude;
  double x = axis.x;
  double y = axis.y;
  double z = axis.z;
  
  if ((x*y*t + z*s) > 0.998) { // north pole singularity detected
    heading = 2*atan2(x*sin(angle/2.0),cos(angle/2.0));
    attitude = PI/2;
    bank = 0;
    return true;
  }
  if ((x*y*t + z*s) < -0.998) { // south pole singularity detected
    heading = -2*atan2(x*sin(angle/2),cos(angle/2));
    attitude = -PI/2;
    bank = 0;
    return true;
  }
  heading = atan2(y * s- x * z * t , 1 - (y*y+ z*z ) * t);
  attitude = asin(x * y * t + z * s) ;
  bank = atan2(x * s - y * z * t , 1 - (x*x + z*z) * t);
  
  return true;
}




QDomNode makeElement(QDomDocument& doc, const QTreeWidgetItem* item){

  if(item->text(0)=="translate"){
    double x0 = item->child(0)->text(1).toDouble();
    double y0 = item->child(1)->text(1).toDouble();
    double z0 = item->child(2)->text(1).toDouble();
    if(x0 == 0 && y0 == 0 && z0==0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      if(x0 != 0){
        QDomElement childelem = doc.createElement("x0");
        QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(y0 != 0){
        QDomElement childelem = doc.createElement("y0");
         QDomNode txtNode = doc.createTextNode(item->child(1)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(z0 != 0){
        QDomElement childelem = doc.createElement("z0");
        QDomNode txtNode = doc.createTextNode(item->child(2)->text(1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      return elem;
    }
  }else if(item->text(0)=="scale"){
    double x0 = item->child(0)->text(1).toDouble();
    double y0 = item->child(1)->text(1).toDouble();
    double z0 = item->child(2)->text(1).toDouble();
    if(x0 == 1 && y0 == 1 && z0==1)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      if(x0 != 1){
        QDomElement childelem = doc.createElement("x0");
         QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(y0 != 1){
        QDomElement childelem = doc.createElement("y0");
         QDomNode txtNode = doc.createTextNode(item->child(1)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(z0 != 1){
        QDomElement childelem = doc.createElement("z0");
        QDomNode txtNode = doc.createTextNode(item->child(2)->text(1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      return elem;
    }
  }else if(item->text(0)=="rotateX"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }else if(item->text(0)=="rotateY"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }else if(item->text(0)=="rotateZ"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }else if(item->text(0)=="cone"){
   
    double x0, y0, z0, r,z1, z2;
    //get values
    std::vector<double> para(8);
    for(int i = 0; i < 8; i++){
      para[i] = item->child(i)->text(1).toDouble();
    }
    double xx1= para[0];
    double yy1= para[1];
    double zz1= para[2];
    double xx2= para[3];
    double yy2= para[4];
    double zz2= para[5];
    double rr1= para[6];
    double rr2= para[7];

    if(zz2< zz1){
       xx2= para[0];
       yy2= para[1];
       zz2= para[2];
       xx1= para[3];
       yy1= para[4];
       zz1= para[5];
       rr2= para[6];
       rr1= para[7]; 
    }

    
    if(rr1==rr2) return QDomNode();
    if(xx1==xx2 &&yy1==yy2){
      if(zz1==zz2) return QDomNode();
      x0 = xx1;
      y0 = yy1;
      z0 = (rr1*zz2-rr2*zz1)/(rr1-rr2);
      r = (rr1-rr2)/(zz1-zz2);
      z1=zz1;
      z2=zz2;
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("x0");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      
      childelem = doc.createElement("y0");
      txtNode = doc.createTextNode(QString("%1").arg(y0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);

      childelem = doc.createElement("z0");
       txtNode = doc.createTextNode(QString("%1").arg(z0));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);

       childelem = doc.createElement("r");
       txtNode = doc.createTextNode(QString("%1").arg(r));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);

       childelem = doc.createElement("z1");
       txtNode = doc.createTextNode(QString("%1").arg(z1));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);

       childelem = doc.createElement("z2");
       txtNode = doc.createTextNode(QString("%1").arg(z2));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);
       return elem;
    }else{
      // z2>z1
      x0 = xx1;
      y0 = yy1;
      z1 = zz1;
      z2 = z1+sqrt((xx1-xx2)*(xx1-xx2)+(yy1-yy2)*(yy1-yy2)+(zz1-zz2)*(zz1-zz2));
      if(z1==z2) return QDomNode();
      z0 = (rr1*z2-rr2*z1)/(rr1-rr2);
      r = (rr1-rr2)/(z1-z2);

      positions3d p1 = positions3d(xx1, yy1, zz1);
      positions3d p2 = positions3d(xx2, yy2, zz2);
      positions3d pp2 = positions3d(xx1, yy1, z2);
      positions3d vfrom = pp2-p1;
      positions3d vto = p2-p1;
      double thetay=0, thetax=0, thetaz=0;
      if(anglesBetween(vfrom, vto, thetay, thetaz, thetax)){
        thetay = thetay*180/PI;
        thetaz = thetaz*180/PI;
        thetax = thetax*180/PI;
        QDomElement elem = doc.createElement(item->text(0));
        
        if(thetay !=0 || thetaz!=0 || thetax!=0){
          QDomElement tran_elem =  doc.createElement("transform");
            {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetay) >1e-3){
            QDomElement elem = doc.createElement("rotateY");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetay));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetaz) >1e-3){
            QDomElement elem = doc.createElement("rotateZ");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetaz));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetax) >1e-3){
            QDomElement elem = doc.createElement("rotateX");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetax));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
           {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(-1*p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          elem.appendChild(tran_elem);
        } 
         
        QDomElement childelem = doc.createElement("x0");
        QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      
        childelem = doc.createElement("y0");
        txtNode = doc.createTextNode(QString("%1").arg(y0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("z0");
        txtNode = doc.createTextNode(QString("%1").arg(z0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("r");
        txtNode = doc.createTextNode(QString("%1").arg(r));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("z1");
        txtNode = doc.createTextNode(QString("%1").arg(z1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("z2");
        txtNode = doc.createTextNode(QString("%1").arg(z2));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        return elem;
      
      }else{
        return QDomNode();
      }
      



    }
    return QDomNode();
  }else if(item->text(0)=="cylinder"){
    
    double x0, y0, r,z1, z2;
      //get values
      vector<double> para(7);
      for(int i = 0; i < 7; i++){
        para[i] = item->child(i)->text(1).toDouble();
      }
      double xx1= para[0];
      double yy1= para[1];
      double zz1= para[2];
      double xx2= para[3];
      double yy2= para[4];
      double zz2= para[5];
      double rr1= para[6];
   




    
     
      if(xx1==xx2 &&yy1==yy2){
       if(zz1 == zz2 ) return QDomNode();
        x0 = xx1;
        y0 = yy1;
        z1 = zz1;
        z2 = zz2;
        r = rr1;
        QDomElement elem = doc.createElement(item->text(0));
        QDomElement childelem = doc.createElement("x0");
        QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("y0");
        txtNode = doc.createTextNode(QString("%1").arg(y0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
            
        childelem = doc.createElement("r");
        txtNode = doc.createTextNode(QString("%1").arg(r));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("z1");
        txtNode = doc.createTextNode(QString("%1").arg(z1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
       childelem = doc.createElement("z2");
       txtNode = doc.createTextNode(QString("%1").arg(z2));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);
       return elem;
      }else{
         // z2>z1
      x0 = xx1;
      y0 = yy1;
      z1 = zz1;
      z2 = z1+sqrt((xx1-xx2)*(xx1-xx2)+(yy1-yy2)*(yy1-yy2)+(zz1-zz2)*(zz1-zz2));
      if(z1==z2) return QDomNode();
      r = rr1;

      positions3d p1 = positions3d(xx1, yy1, zz1);
      positions3d p2 = positions3d(xx2, yy2, zz2);
      positions3d pp2 = positions3d(xx1, yy1, z2);
      positions3d vfrom = pp2-p1;
      positions3d vto = p2-p1;
      double thetay=0, thetax=0, thetaz=0;
      if(anglesBetween(vfrom, vto, thetay, thetaz, thetax)){
        thetay = thetay*180/PI;
        thetaz = thetaz*180/PI;
        thetax = thetax*180/PI;
        QDomElement elem = doc.createElement(item->text(0));
        
        if(thetay !=0 || thetaz!=0 || thetax!=0){
          QDomElement tran_elem =  doc.createElement("transform");

          {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          
          if(fabs(thetay) >1e-3){
            QDomElement elem = doc.createElement("rotateY");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetay));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetaz) >1e-3){
            QDomElement elem = doc.createElement("rotateZ");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetaz));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetax) >1e-3){
            QDomElement elem = doc.createElement("rotateX");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetax));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
           {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(-1*p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          elem.appendChild(tran_elem);
        } 
        QDomElement childelem = doc.createElement("x0");
        QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("y0");
        txtNode = doc.createTextNode(QString("%1").arg(y0));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
            
        childelem = doc.createElement("r");
        txtNode = doc.createTextNode(QString("%1").arg(r));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
        childelem = doc.createElement("z1");
        txtNode = doc.createTextNode(QString("%1").arg(z1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        
       childelem = doc.createElement("z2");
       txtNode = doc.createTextNode(QString("%1").arg(z2));
       childelem.appendChild(txtNode);
       if(!childelem.isNull()) elem.appendChild(childelem);
       return elem; 

      }
      }
      return QDomNode();
  }else if(item->text(0)=="leftplane"){
   
    double x0;
    //get values
    vector<double> para(6);
    for(int i = 0; i < 6; i++){
      para[i] = item->child(i)->text(1).toDouble();
    }
    //point
    double xx1= para[0];
    double yy1= para[1];
    double zz1= para[2];
    //normal
    double xx2= para[3];
    double yy2= para[4];
    double zz2= para[5];
    
    if(xx2==0 &&yy2==0 &&zz2==1){
      
      x0 = zz1;
      QDomElement elem = doc.createElement("z_minus_plane");
      QDomElement childelem = doc.createElement("z1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else if(xx2==0 &&yy2==0 &&zz2==-1){
      
      x0 = zz1;
      QDomElement elem = doc.createElement("z_plus_plane");
      QDomElement childelem = doc.createElement("z1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else if(xx2==1 &&yy2==0 &&zz2==0){
      
      x0 = xx1;
      QDomElement elem = doc.createElement("x_minus_plane");
      QDomElement childelem = doc.createElement("x1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else if(xx2==-1 &&yy2==0 &&zz2==0){
      
      x0 = xx1;
      QDomElement elem = doc.createElement("x_plus_plane");
      QDomElement childelem = doc.createElement("x1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else if(xx2==0 &&yy2==1 &&zz2==0){
      
      x0 = yy1;
      QDomElement elem = doc.createElement("y_minus_plane");
      QDomElement childelem = doc.createElement("y1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else if(xx2==0 &&yy2==-1 &&zz2==0){
      
      x0 = yy1;
      QDomElement elem = doc.createElement("y_plus_plane");
      QDomElement childelem = doc.createElement("y1");
      QDomNode txtNode = doc.createTextNode(QString("%1").arg(x0));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }else{
      positions3d p1 = positions3d(xx1, yy1, zz1);
      positions3d vfrom = positions3d(0, 0, 1);
      positions3d vto = positions3d(xx2, yy2, zz2);
      double thetay=0, thetax=0, thetaz=0;
      if(anglesBetween(vfrom, vto, thetay, thetaz, thetax)){
        thetay = thetay*180/PI;
        thetaz = thetaz*180/PI;
        thetax = thetax*180/PI;
        QDomElement elem = doc.createElement("z_minus_plane");
        
        if(thetay !=0 || thetaz!=0 || thetax!=0){
          QDomElement tran_elem =  doc.createElement("transform");
          {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          
          if(fabs(thetay) >1e-3){
            QDomElement elem = doc.createElement("rotateY");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetay));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetaz) >1e-3){
            QDomElement elem = doc.createElement("rotateZ");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetaz));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          if(fabs(thetax) >1e-3){
            QDomElement elem = doc.createElement("rotateX");
            QDomElement childelem = doc.createElement("theta");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(thetax));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);
            if(!elem.isNull()) tran_elem.appendChild(elem);
          }

          {
            QDomElement elem = doc.createElement("translate");
            QDomElement childelem = doc.createElement("x0");
            QDomNode txtNode = doc.createTextNode(QString("%1").arg(-1*p1.x));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("y0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.y));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            childelem = doc.createElement("z0");
            txtNode = doc.createTextNode(QString("%1").arg(-1*p1.z));
            childelem.appendChild(txtNode);
            if(!childelem.isNull()) elem.appendChild(childelem);

            if(!elem.isNull()) tran_elem.appendChild(elem);
          }
          elem.appendChild(tran_elem);
        }
        
        QDomElement childelem = doc.createElement("z1");
        QDomNode txtNode = doc.createTextNode(QString("%1").arg(zz1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
        return elem;
      }
    }
      
    return QDomNode();
  }
  

    
  QDomElement elem = doc.createElement(item->text(0));
  if(item->text(1)!=""){
    QDomNode elt = doc.createTextNode(item->text(1));
    elem.appendChild(elt);
    // return elem;
  }
  if(item->childCount()!=0){
    for(int i = 0; i < item->childCount(); i++){
      QDomNode childElem = makeElement(doc, item->child(i));
      if(!childElem.isNull())elem.appendChild(childElem);
    }
    if(item->text(0)=="transform" && elem.firstChildElement().isNull()) return QDomNode();
   
  }
  return elem;
 
}


void graft_tree(QDomDocument& doc){
  QDomNodeList shapeList = doc.elementsByTagName("shape");
  for(int i = 0; i < shapeList.size(); i++){
    QDomElement  trans_elem = shapeList.at(i).firstChildElement().firstChildElement("transform");
    if(!trans_elem.isNull()){
      QDomNode tmpNode =shapeList.at(i).firstChildElement().removeChild(trans_elem);
      QDomNode parent = shapeList.at(i).parentNode();
      parent.insertBefore(tmpNode, parent.firstChildElement());
    }
  }

}

QDomDocument tree2dom(const QTreeWidgetItem* root){
  QDomDocument doc;
  QDomNode rootElem = makeElement(doc, root);
  doc.appendChild(rootElem);
  graft_tree(doc);
  return doc;
}




  

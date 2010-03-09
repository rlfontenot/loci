#include <QtDebug>

#include <set>
using std::set;
#include <cmath>
using std::signbit;

#include "hdf5.h"
#include "grid.h"
#include "glviewer.h"





//////////////////////////////////////////////////////////////////////////////
//  private:
//    void makeObjects();
//
//  This function calls all the other functions to remake all necessary display
//  lists.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::makeObjects()
{
  
 
  if (cpContourObject) glDeleteLists(cpContourObject, 1);

  if (mode == PLANE_AND_BOUND_MODE || mode == PLANE_ONLY_MODE) {
    if (gridObject) glDeleteLists(gridObject, 1);
    if (contourObject) glDeleteLists(contourObject, 1);
    if (borderObject) glDeleteLists(borderObject, 1);
    if (shadingObject) glDeleteLists(shadingObject, 1);
    // Make all cutting plane display lists
    shadingObject = makeShadingObject();
    gridObject = makeGridObject();
    contourObject = makeContourObject();
    borderObject = makeBorderObject();
  }

  else if (mode == BOUND_SELECT_MODE){
    if (!boundObjects.empty()) {
      for (size_t i = 0; i < boundObjects.size(); ++i)
        glDeleteLists(boundObjects[i], 1);
      boundObjects.clear();
    }
    makeBoundObjects();
  }
 
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    GLuint makeGridObject();
//
//  This function compiles the display list for the 2d cutting plane grid.
//////////////////////////////////////////////////////////////////////////////

GLuint GLViewer::makeGridObject()
{
  GLuint newList = glGenLists(1);
  int nedges = fig->interior;
  affineMapping transMatrix2;
  positions3d negCenter = positions3d(-centerx, -centery, -centerz);
  transMatrix2.translate(info.translate);
  transMatrix2.rotateZ(info.rotate.z);
  transMatrix2.rotateY(info.rotate.y);
  transMatrix2.rotateX(info.rotate.x);
  transMatrix2.translate(negCenter);
  
  vector<positions3d> tmpPos;
  for(size_t i =0; i <fig->pos.size();i++){
    positions3d aNode = positions3d(fig->pos[i].x, fig->pos[i].y, 0);
    tmpPos.push_back(transMatrix2.MapNode(aNode));
  }
  glNewList(newList, GL_COMPILE);
  glBegin(GL_LINES);
  glColor3f(0.75, 0.75, 0.75);
  for (int e = 0; e < nedges; ++e) {
    const edges &ed = fig->edge_list[e];
    glVertex3d(tmpPos[ed.l].x, tmpPos[ed.l].y, tmpPos[ed.l].z);
    glVertex3d(tmpPos[ed.r].x, tmpPos[ed.r].y, tmpPos[ed.r].z);
        
  }
  glEnd();
  glEndList();

  return newList;
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    GLuint makeContourObject();
//
//  This function compiles the display list for the cutting plane's contour
//  lines.
//////////////////////////////////////////////////////////////////////////////

GLuint GLViewer::makeContourObject()
{
  GLuint newList = glGenLists(1);
  int nsegs = fig->contour_curves.size();
  affineMapping transMatrix2;
  positions3d negCenter = positions3d(-centerx, -centery, -centerz);
  transMatrix2.translate(info.translate);
  transMatrix2.rotateZ(info.rotate.z);
  transMatrix2.rotateY(info.rotate.y);
  transMatrix2.rotateX(info.rotate.x);
  transMatrix2.translate(negCenter);

  vector<positions3d> tmpP1;
  vector<positions3d> tmpP2;
  
 for (int s = 0; s < nsegs; ++s) {
   const segments &seg = fig->contour_curves[s];
   positions3d p1 = positions3d(seg.p1.x, seg.p1.y, 0);
   positions3d p2 = positions3d(seg.p2.x, seg.p2.y, 0);
   p1 = transMatrix2.MapNode(p1);
   p2 = transMatrix2.MapNode(p2);
   tmpP1.push_back(p1);
   tmpP2.push_back(p2);
 }
    
  glNewList(newList, GL_COMPILE);
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  for (int s = 0; s < nsegs; ++s) {
    //  const segments &seg = fig->contour_curves[s];
    glVertex3d(tmpP1[s].x, tmpP1[s].y, tmpP1[s].z);
    glVertex3d(tmpP2[s].x, tmpP2[s].y, tmpP2[s].z);

  }
  glEnd();
  glEndList();

  return newList;
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    GLuint makeBorderObject();
//
//  This function compiles the display list for the cutting plane's border.
//////////////////////////////////////////////////////////////////////////////

GLuint GLViewer::makeBorderObject()
{
  GLuint newList = glGenLists(1);
  int nedges = fig->edge_list.size();

  affineMapping transMatrix2;
  positions3d negCenter = positions3d(-centerx, -centery, -centerz);
 transMatrix2.translate(info.translate);
  transMatrix2.rotateZ(info.rotate.z);
  transMatrix2.rotateY(info.rotate.y);
  transMatrix2.rotateX(info.rotate.x);
  transMatrix2.translate(negCenter);
 
  vector<positions3d> tmpP1;
  vector<positions3d> tmpP2;
  
  for (int e = fig->interior; e < nedges; ++e) {
    const edges &ed = fig->edge_list[e];
   positions3d p1 = positions3d(fig->pos[ed.l].x, fig->pos[ed.l].y, 0);
   positions3d p2 = positions3d(fig->pos[ed.r].x, fig->pos[ed.r].y, 0);
   p1 = transMatrix2.MapNode(p1);
   p2 = transMatrix2.MapNode(p2);
   tmpP1.push_back(p1);
   tmpP2.push_back(p2);
 }
  
  glNewList(newList, GL_COMPILE);
  glBegin(GL_LINES);
  glColor3f(0.0, 1.0, 0.0);
  for (int e = 0; e <(int)(nedges-fig->interior); ++e) {
    glVertex3d(tmpP1[e].x, tmpP1[e].y, tmpP1[e].z);
    glVertex3d(tmpP2[e].x, tmpP2[e].y, tmpP2[e].z);
    //glVertex3d(fig->pos[ed.l].x, fig->pos[ed.l].y, 0.0);
    //glVertex3d(fig->pos[ed.r].x, fig->pos[ed.r].y, 0.0);
  }
  glEnd();
  glEndList();

  return newList;
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    GLuint makeShadingObject();
//
//  This function compiles the display list for the cutting plane's color
//  shading.
//////////////////////////////////////////////////////////////////////////////

GLuint GLViewer::makeShadingObject()
{
  GLuint newList = glGenLists(1);
  int ntris = fig->triangle_list.size();

  affineMapping transMatrix2;
  positions3d negCenter = positions3d(-centerx, -centery, -centerz);
  transMatrix2.translate(info.translate);
  transMatrix2.rotateZ(info.rotate.z);
  transMatrix2.rotateY(info.rotate.y);
  transMatrix2.rotateX(info.rotate.x);
  transMatrix2.translate(negCenter);
 
  vector<positions3d> tmpPos;
  for(size_t i = 0; i < fig->pos.size(); i++){
    positions3d aNode = positions3d(fig->pos[i].x, fig->pos[i].y, 0);
    tmpPos.push_back(transMatrix2.MapNode(aNode));
  }
    
  glNewList(newList, GL_COMPILE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_TRIANGLES);
  for (int t = 0; t < ntris; ++t) {
    const triangles &tri = fig->triangle_list[t];
    
    positions3d p1 = shade(fig->val[tri.t1]);
    glColor3d(p1.x, p1.y, p1.z);
    glVertex3d(tmpPos[tri.t1].x, tmpPos[tri.t1].y, tmpPos[tri.t1].z);
    positions3d p2 = shade(fig->val[tri.t2]);
    glColor3d(p2.x, p2.y, p2.z);
    glVertex3d(tmpPos[tri.t2].x, tmpPos[tri.t2].y, tmpPos[tri.t2].z);
    positions3d p3 = shade(fig->val[tri.t3]);
    glColor3d(p3.x, p3.y, p3.z);
    glVertex3d(tmpPos[tri.t3].x, tmpPos[tri.t3].y, tmpPos[tri.t3].z);
  }
  glEnd();
  glEndList();
 
  return newList;
}



void GLViewer::drawBoundObject(int bid, QColor c){
  if(meshNodes.empty() ||mesh.empty()){
    qDebug("please load in boundary first");
    return;
  }
  qglColor(c);
  for (size_t i = 0; i < mesh[bid].size() / 3; ++i) {
    glBegin(GL_TRIANGLES);
    positions3d t0, t1, t2;
    t0 = meshNodes[mesh[bid][i*3 + 0]]; 
    t1 = meshNodes[mesh[bid][i*3 + 1]]; 
    t2 = meshNodes[mesh[bid][i*3 + 2]];
    glVertex3d(t0.x, t0.y, t0.z);
    glVertex3d(t1.x, t1.y, t1.z);
    glVertex3d(t2.x, t2.y, t2.z);
    glEnd();  
  }
}
 
  
 

  

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void makeBoundObjects();
//
//  This function compiles the display lists for each boundary and places them
//  in boundObjects.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::makeBoundShadingObject( int bid)
{
  if(meshNodes.empty() ||mesh.empty()||meshValues.empty()){
    qDebug("please use file menu load in boundary first");
    return;
  }
  
 
  // make display lists
  
  // qglColor(default_color[bid%12]);

  float r, g, b;
  r = default_color[bid%12].red();
  g = default_color[bid%12].green();
  b = default_color[bid%12].blue();
  // compile list
    GLuint newList = glGenLists(1);
    glNewList(newList, GL_COMPILE);
    glBegin(GL_LINE_LOOP);
    
    
    for (size_t i = 0; i < mesh[bid].size() / 3; ++i) {
      positions3d t0, t1, t2, u, v, norm;
      t0 = meshNodes[mesh[bid][i*3 + 0]]; 
      t1 = meshNodes[mesh[bid][i*3 + 1]]; 
      t2 = meshNodes[mesh[bid][i*3 + 2]];
      
      // Calculate the normals for lighting effects
      u.x = t2.x - t0.x; u.y = t2.y - t0.y; u.z = t2.z - t0.z;
      v.x = t1.x - t0.x; v.y = t1.y - t0.y; v.z = t1.z - t0.z;
      
      norm.x = u.y * v.z - u.z * v.y;
      norm.y = u.z * v.x - u.x * v.z;
      norm.z = u.x * v.y - u.y * v.x;
      double len = sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
      norm.x /= len; norm.y /= len; norm.z /= len;
      len = norm.x + norm.y + norm.z;
      len *= 0.24;
      len += 0.55;
      
      glColor3f(r*len, g*len, b*len);
      glVertex3d(t0.x, t0.y, t0.z);
      glColor3f(r*len, g*len, b*len);
      glVertex3d(t1.x, t1.y, t1.z);
      glColor3f(r*len, g*len, b*len);
      glVertex3d(t2.x, t2.y, t2.z);
    }
    
    glEnd();    
    glEndList();
    GLuint oldList = boundObjects[bid];
    boundObjects[bid] = newList;
    glDeleteLists(oldList, 1);
    
      
}

void GLViewer::makeBoundWireframeObject(int bid, QColor c)
{
    // compile list
  GLuint newList = glGenLists(1);
  glNewList(newList, GL_COMPILE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
  drawBoundObject(bid, c);   
  glEndList();
  
  GLuint oldList = boundObjects[bid];
  boundObjects[bid] = newList;
  glDeleteLists(oldList, 1);
}

void GLViewer::makeBoundFillObject(int bid, QColor c)
{
  // compile list
  GLuint newList = glGenLists(1);
  glNewList(newList, GL_COMPILE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);    
  drawBoundObject(bid, c); 
  glEndList();
  GLuint oldList = boundObjects[bid];
  boundObjects[bid] = newList;
  glDeleteLists(oldList, 1);
}



void GLViewer::makeBoundObjects()
{
  
  if(boundObjects.empty()){
      for (size_t bid = 0; bid < mesh.size(); ++bid) {
    
      
      // compile list
      GLuint newList = glGenLists(1);
      glNewList(newList, GL_COMPILE);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      drawBoundObject(bid, default_color[bid%12]);   
      glEndList();
      boundObjects.push_back(newList);
  
      }//for(bid...)
  }

}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    GLuint makeCPContour();
//
//  This function compiles the display list for the cut preview.
// cut preview cut the boundaries
//////////////////////////////////////////////////////////////////////////////

GLuint GLViewer::makeCPContour()
{

  vector<positions3d> contourLoop;

  // Make temporary node list
  vector<positions3d> tempNodes(meshNodes.size());

  // Set up transformation matrix
  positions3d negCenter = positions3d(-centerx, -centery, -centerz);
  
  positions3d center = positions3d(centerx, centery, centerz);
  positions3d negTranslate = positions3d(-previewInfo.translate.x,
                                         -previewInfo.translate.y,
                                         -previewInfo.translate.z); 
  
  
  affineMapping transMatrix;
  transMatrix.translate(center);
  transMatrix.rotateX(-previewInfo.rotate.x);
  transMatrix.rotateY(-previewInfo.rotate.y);
  transMatrix.rotateZ(-previewInfo.rotate.z);
  transMatrix.translate(negTranslate);  

  // Transform tempNodes
 
  for (int i = 0; i < (int)meshNodes.size(); ++i){
    tempNodes[i]= transMatrix.MapNode(meshNodes[i]);
  }
   
  double tol = 1e-34;
  // For each triangle in mesh
  for (size_t bid = 0; bid < mesh.size(); ++bid) {
    for(size_t i = 0; i < mesh[bid].size()/3; i++){
      positions3d tri[3];
      for (int j = 0; j < 3; ++j)
        tri[j] = tempNodes[ mesh[bid][i*3 + j] ];
      
      int cutsFound = 0;
      

      //if one node is on the z = 0 plane
      if(fabs(tri[0].z) < tol){
        if(fabs(tri[1].z) < tol){
          contourLoop.push_back(tri[0]);
          cutsFound++;
          contourLoop.push_back(tri[1]);
          cutsFound++;
          
        }
        if(fabs(tri[2].z) < tol){
          contourLoop.push_back(tri[2]);
          cutsFound++;
          contourLoop.push_back(tri[0]);
          cutsFound++;
        }
      }else if(fabs(tri[1].z) < tol){
        if(fabs(tri[2].z) < tol){
          contourLoop.push_back(tri[1]);
          cutsFound++;
          contourLoop.push_back(tri[2]);
          cutsFound++;
        }
      }

      //no node on z=0 plane
      if(fabs(tri[0].z) >= tol && fabs(tri[1].z) >= tol && fabs(tri[2].z) >=tol){ 
        // Try to cut each edge
        if((tri[0].z * tri[1].z)< 0) {
          double t;
          t = tri[0].z/(tri[0].z-tri[1].z);
                  
          positions3d newPnt;
          newPnt.x = (1-t)*tri[0].x + t*tri[1].x;
          newPnt.y = (1-t)*tri[0].y + t*tri[1].y;
          newPnt.z = 0.0;

         
          contourLoop.push_back(newPnt);
          cutsFound++;
        }      
        
        if((tri[1].z * tri[2].z)<0) {
        double t;
        t = tri[1].z/(tri[1].z-tri[2].z);
        positions3d newPnt;
        newPnt.x = (1-t)*tri[1].x + t*tri[2].x;
        newPnt.y = (1-t)*tri[1].y + t*tri[2].y;
        newPnt.z = 0.0;
       
        contourLoop.push_back(newPnt);
        cutsFound++;
        }
        
        if((tri[2].z * tri[0].z)<0) {
          double t;
          t = tri[2].z/(tri[2].z-tri[0].z);
          positions3d newPnt;
          newPnt.x = (1-t)*tri[2].x + t*tri[0].x;
          newPnt.y = (1-t)*tri[2].y + t*tri[0].y;
          newPnt.z = 0.0;
          
          contourLoop.push_back(newPnt);
          cutsFound++;
        }

      }

      if (cutsFound != 0 && cutsFound != 2) {
        qDebug() << "Major malfunction at tri #" << i;
        qDebug() << "Cuts: " << cutsFound;
        qDebug() << tri[0].z<< " "<< tri[1].z<<" " << tri[2].z;
        exit(0);
      }
    }
  }

  if(contourLoop.size()==0) return 0;
  
  //move the cut plane back
  affineMapping transMatrix2;
  transMatrix2.translate(previewInfo.translate);
  transMatrix2.rotateZ(previewInfo.rotate.z);
  transMatrix2.rotateY(previewInfo.rotate.y);
  transMatrix2.rotateX(previewInfo.rotate.x);
  transMatrix2.translate(negCenter);
 
  
  for (size_t i = 0; i < contourLoop.size(); ++i)
    contourLoop[i] = transMatrix2.MapNode(contourLoop[i]);
  
  
  // Compile display list
  GLuint newList = glGenLists(1);
  glNewList(newList, GL_COMPILE);
  glColor3f(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  for(size_t i = 0; i < contourLoop.size(); i++) 
    glVertex3d(contourLoop[i].x, contourLoop[i].y, contourLoop[i].z);
  glEnd();

  
  // glBegin(GL_LINE_LOOP);
  /*
  double len = size / 1.5;
  glVertex3d(len + center.x, len + center.y, 0.0);
  glVertex3d(len + center.x, -len + center.y, 0.0);
  glVertex3d(-len + center.x, -len + center.y, 0.0);
  glVertex3d(-len + center.x, len + center.y, 0.0);
  */
   //double len = size / 1.5;
    // glVertex3d(len, len, 0.0);
  // glVertex3d(len, -len, 0.0);
  // glVertex3d(-len, -len, 0.0);
  // glVertex3d(-len, len, 0.0);

  //glEnd();
  glEndList();

  return newList;
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void drawBoxes();
//
//  This function draws boxes around each boundary.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::drawBoxes()
{
  // int div = 4;
  //  float r, g, b;

  for (int i = 0; i < (int)objMinMax.size()/2; ++i) {
   
       // Draw box
    positions3d minpos, maxpos;
    minpos = objMinMax[2*i + 0];
    maxpos = objMinMax[2*i + 1];

    glBegin(GL_LINE_STRIP);
    qglColor(default_color[i%12]);
    glVertex3d(maxpos.x, maxpos.y, maxpos.z);
    glVertex3d(minpos.x, maxpos.y, maxpos.z);
    glVertex3d(minpos.x, minpos.y, maxpos.z);
    glVertex3d(maxpos.x, minpos.y, maxpos.z);
    glVertex3d(maxpos.x, maxpos.y, maxpos.z);
    glVertex3d(maxpos.x, maxpos.y, minpos.z);
    glVertex3d(minpos.x, maxpos.y, minpos.z);
    glVertex3d(minpos.x, minpos.y, minpos.z);
    glVertex3d(maxpos.x, minpos.y, minpos.z);
    glVertex3d(maxpos.x, maxpos.y, minpos.z);
    glEnd();
    
    glBegin(GL_LINES);
     qglColor(default_color[i%12]);
    glVertex3d(minpos.x, maxpos.y, maxpos.z);
    glVertex3d(minpos.x, maxpos.y, minpos.z);
    glVertex3d(minpos.x, minpos.y, maxpos.z);
    glVertex3d(minpos.x, minpos.y, minpos.z);
    glVertex3d(maxpos.x, minpos.y, maxpos.z);
    glVertex3d(maxpos.x, minpos.y, minpos.z);
    glEnd();
  }
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void drawCoordGrid();
//
//  This function draws Cartesian coordinates for a visual reference. Also
//  drawn is the center of camera (xPos, yPos) in red and the center of the
//  visible boundaries (center) in blue.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::drawCoordGrid()
{
 //  double length = size / 2.0;

//   // Draw grid
//   glBegin(GL_LINES);
//   glColor3f(0.0, 0.0, 0.0);
//   glVertex3d(length, 0.0, 0.0); glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(0.0, 0.0, length); glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(-length, 0.0, 0.0); glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(0.0, 0.0, -length); glVertex3d(0.0, 0.0, 0.0);

//   glColor3f(0.5, 0.5, 0.5);
//   for (int i = -10; i <= 10; ++i) {
//     if (i == 0) continue;
//     double scan = size * i/20.0;
//     glVertex3d(length, 0.0, scan); glVertex3d(-length, 0.0, scan);
//     glVertex3d(scan, 0.0, length); glVertex3d(scan, 0.0, -length);
//   }
//   glEnd();

//   // Draw view origin
//   glMatrixMode(GL_MODELVIEW);
//   glTranslatef(xPos, 0.0, yPos);

//   glColor3f(0.5, 0.2, 0.2);
//   glBegin(GL_LINES);
//   glVertex3d(rho*0.02, 0.0, 0.0); glVertex3d(-rho*0.02, 0.0, 0.0);
//   glVertex3d(0.0, rho*0.02, 0.0); glVertex3d(0.0, -rho*0.02, 0.0);
//   glVertex3d(0.0, 0.0, rho*0.02); glVertex3d(0.0, 0.0, -rho*0.02);
//   glEnd();

//   glLoadIdentity();
//   //  glMatrixMode(GL_PROJECTION);

//   // Draw center point
//   //  glMatrixMode(GL_MODELVIEW);
//   //  glLoadIdentity();
//   glTranslatef(center.x, center.y, center.z);

//   glColor3f(0.2, 0.2, 0.5);
//   glBegin(GL_LINES);
//   glVertex3d(rho*0.02, 0.0, 0.0); glVertex3d(-rho*0.02, 0.0, 0.0);
//   glVertex3d(0.0, rho*0.02, 0.0); glVertex3d(0.0, -rho*0.02, 0.0);
//   glVertex3d(0.0, 0.0, rho*0.02); glVertex3d(0.0, 0.0, -rho*0.02);
//   glEnd();

//   glLoadIdentity();
//   glMatrixMode(GL_PROJECTION);
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    double* shade(double value, double weight = 1.0);
//
//  This function takes a value and returns an array of doubles representing
//  the RGB color of that value under the current shading mode. The weight
//  passed is used to darken the color according to the shading of the node
//  due to the angle of the light source.
//////////////////////////////////////////////////////////////////////////////

positions3d GLViewer::shade(double value, double weight)
{
  double rgb[3];
  
  value = (value - min_val)/(max_val - min_val);
  // 0.0 <= value <= 1.0

  switch (shadeType) {
  case 1: // blue to red
    value *= 4;
    if (value < 1.0) {
      rgb[0] = 0.0;
      rgb[1] = value;
      rgb[2] = 1.0;
    } else if (value >= 1.0 && value < 2.0) {
      rgb[0] = 0.0;
      rgb[1] = 1.0;
      rgb[2] = -(value - 2.0);
    } else if (value >= 2.0 && value < 3.0) {
      rgb[0] = value - 2.0;
      rgb[1] = 1.0;
      rgb[2] = 0.0;
    } else { // (value >= 3.0)
      rgb[0] = 1.0;
      rgb[1] = -(value - 4.0);
      rgb[2] = 0.0;
    }
    break;
  case 2: // blackbody
    value *= 3;
    if (value < 1.0) {
      rgb[0] = value;
      rgb[1] = 0.0;
      rgb[2] = 0.0;
    } else if (value >= 1.0 && value < 2.0) {
      rgb[0] = 1.0;
      rgb[1] = value - 1;
      rgb[2] = 0.0;
    } else { // (value >= 2.0)
      rgb[0] = 1.0;
      rgb[1] = 1.0;
      rgb[2] = value - 2;
    }
    break;
  case 3: // blue to red (pressure)
    value *= 2;
    if (value < 1.0) {
      rgb[0] = value;
      rgb[1] = value;
      rgb[2] = 1.0;
    } else { // (value >= 1.0)
      rgb[0] = 1.0;
      rgb[1] = -(value - 2.0);
      rgb[2] = -(value - 2.0);
    }
    break;
  }

  rgb[0] *= weight;
  rgb[1] *= weight;
  rgb[2] *= weight;

  return positions3d(rgb[0], rgb[1], rgb[2]);
}

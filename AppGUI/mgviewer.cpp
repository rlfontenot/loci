#include <QtOpenGL>
#include <QString>
#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <GL/glut.h>

#include <map>
#include <vector>
#include <fstream>
#include <strings.h>
#include <stdio.h>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;
using std::set;
using std::map ;
using std::vector ;
using std::string ;
using std::cerr ;


#include "mgviewer.h"
#include "grid.h"
#include "hdf5.h"

#define PI 3.14159265358979323846264338327950
//the definition of affineMapping2 is copied from vogmerge.cc, should always
//be consistent with vogmerge.cc
















//////////////////////////////////////////////////////////////////////////////
//  Global Function
//
//  Places node index into the passed integar vector.
//////////////////////////////////////////////////////////////////////////////

void cbVertex2(void *vertex_data, void *user_data);
// {
//   vector<int>* tri = (vector<int>*) user_data;
//   tri->push_back(*(int*)vertex_data);
// }



///////////////////////////////////////////
//  public:
//    MGViewer(QWidget *parent = 0);
//
//  Initializes values.
///////////////////////////////////////////

MGViewer::MGViewer(QWidget *parent)
  : QGLWidget(parent)
{
  makeCurrent();
  centerx = centery = centerz=0;
  size = 0.1;
  currentObj = -1; //no select obj
  currentGrid = -1; //no grid loaded
  
  currentColor = default_color[0];
  tox= toy=toz = rox = roy = 0;
  
  scale = 1.0;
  min_val = max_val = 0.0;
  isFit = false;
  //  mode=BOUND_SELECT_MODE;
 
}

//////////////////////////////////////////////////////
//  public:
//    ~MGViewer();
//
//  Deletes display lists before widget destruction.
//////////////////////////////////////////////////////

MGViewer::~MGViewer()
{
  makeCurrent();
   if (boundObjects.size() > 0) {
    for (size_t i = 0; i < boundObjects.size(); ++i){
      for(size_t j = 0; j < boundObjects[i].size(); ++j)
        glDeleteLists(boundObjects[i][j], 1);
      boundObjects[i].clear();
    }
    boundObjects.clear();
   }
  
}








///////////////////////////////////////////
//  private:
//    void initializeGL();
//
//  Does some basic OpenGL configuration.
///////////////////////////////////////////

void MGViewer::initializeGL()
{
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
}
/////////////////////////////////////////////////////////////
//  protected:
//    void resizeGL(int height, int width);
//
//  Anytime the window is resized, this function handles it.
/////////////////////////////////////////////////////////////

void MGViewer::resizeGL(int width, int height)
{

  currentWidth = width;
  currentHeight = height;
  glViewport(0, 0, (GLint)width,(GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //change the parameters here also need change unproject parameters
  GLdouble near = 0.001*size;
  GLdouble far = size*30.0;
  gluPerspective(30, (GLdouble)width/height, near, far);
  glMatrixMode(GL_MODELVIEW);
 
 
}
////////////////////////////////////////////////////////////////////////////
//  protected:
//    void paintGL();
//
//  This function calls all the OpenGL commands and display lists to draw
//  whatever needs to be in the central widget.
////////////////////////////////////////////////////////////////////////////

void MGViewer::paintGL()
{
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
  gluLookAt(centerx, centery, 2.0*size+centerz, centerx, centery, centerz, 0.0, 1.0, 0.0); 
  glPushMatrix(); //without this, the rotation won't work well
  glScaled(scale,scale,scale);
  glTranslated(tox, toy, toz);
  glTranslated(centerx, centery, centerz);
  glRotated(rox, 1, 0, 0);
  glRotated(roy, 0, 1, 0);
  glRotated(roz, 0, 0, 1);
  glTranslated(-centerx, -centery, -centerz);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatr);
 
 
  
  for (size_t i = 0; i < boundObjects.size(); ++i)
    for(size_t j = 0; j <boundObjects[i].size(); ++j)
      if (objVisible[i][j])glCallList(boundObjects[i][j]);
      
  
  glPopMatrix();
  glFlush();
}
//select one boundary from current grid
void MGViewer::setCurrentColor(const IDColor& idColor){

  
  int gid = idColor.gridId;
  int bid = idColor.boundId;
  if(gid < 0 || bid < 0 ||gid >= (int)boundObjects.size() || bid >= (int)boundObjects[gid].size())return;
  if(currentGrid >= 0 && currentObj >=0){
    makeBoundWireframeObject(currentGrid, currentObj, currentColor);
  }
  currentGrid = gid;
  currentObj = bid;
  currentColor = idColor.color;
  makeBoundFillObject(gid, bid, idColor.color);
  updateGL();
  
}

void MGViewer::clearCurrent(){
  if(currentObj!=-1 && currentGrid!= -1){
    makeBoundWireframeObject(currentGrid, currentObj, currentColor);
    currentObj=-1;
    currentGrid = -1;
    updateGL();
  }
}

void MGViewer::reset(){
  
  tox= toy=toz=rox=roy =roz= 0;
  if(isFit){
    for(unsigned int i =0; i < objVisible.size(); i++){
      for(unsigned int j=0; j < objVisible[i].size(); j++){
      objVisible[i][j] = true;
      }
    }
    isFit = false;
    updateView();
  }
  updateGL();
}

void MGViewer::fit(){
  if(currentObj<0){
    QMessageBox::warning(window(), tr("fit"),
                         tr("please select a boundary first")
                         );
    return;
  }
  isFit = true;
  for( int i  = 0; i < (int)objVisible.size(); i++){
    for( int j =0; j < (int)objVisible[i].size(); j++){
      if(i==currentGrid && j == currentObj) objVisible[i][j] = true;
      else objVisible[i][j] = false;
    }
  }
  updateView();
  updateGL();
}




struct surface_info {
  vector<int> trias ;
  vector<int> quads ;
  vector<vector<int> > gen_faces ;
} ;

//first clear up, then load boundary
bool MGViewer::get_boundary(QString fileName) {
  if (boundObjects.size() > 0) {
    for (size_t i = 0; i < boundObjects.size(); ++i){
      for(size_t j = 0; j < boundObjects[i].size(); ++j)
        glDeleteLists(boundObjects[i][j], 1);
      boundObjects[i].clear();
    }
    boundObjects.clear();
  }
  meshNodes.clear();
  mesh.clear();
  objVisible.clear();
  objMinMax.clear();
  gridXform.clear();
  if(fileName.isEmpty()){
    updateGL();
    
    return false;
  }
  return load_boundary(fileName);
}


bool MGViewer::load_boundary(QString fileName) {

  //if .surf file doesn't exist or it's created before .vog file
  //run vog2surf
  QString surfFileName = fileName.section('.', 0, -2)+".surf";
  QString nameFileName =  fileName.section('.', 0, -2)+".names";
  QFileInfo surfInfo(surfFileName);
  QFileInfo vogInfo(fileName);
  QFileInfo nameInfo(nameFileName);
  
  if(!(surfInfo.exists()) || !(nameInfo.exists())|| surfInfo.created() < vogInfo.created()){
    QString command2 = "vog2surf " + fileName.section('.', 0, -2);
    int ret =  system(command2.toStdString().c_str());
    if(!WIFEXITED(ret))
      {
          if(WIFSIGNALED(ret))
            {
              QMessageBox::information(window(), "mainwindow",
                                       command2 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );

              return false;
            }
      }
    }
  //setup load info
  QStringList boundary_names;
  boundary_names << fileName;
  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
    QString casename = fileName.mid(first+1, last-first-1);
    QString directory = fileName.left(first);
    
    
    QFileInfo surfInfo2(surfFileName);
    
    if(!(surfInfo2.exists())) return false;
    //read in .surf file
        
    if(surfFileName.right(5) == ".surf"){
   
    //first read in meshNodes
    vector<vector<int> > tmpmesh;
    vector<vector3d<double> > pos;
    
    QFile file(surfFileName);  
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Application"),
                           tr("Cannot read file %1:\n%2.")
                           .arg(surfFileName)
                           .arg(file.errorString()));
      return false;
    }
    
    QTextStream in(&file); 
    int Read_Flag = 0; 
  
    
    int total_num_tris = 0, total_num_quads = 0, num_nodes = 0; 
    in >> total_num_tris >> total_num_quads >> num_nodes;
       int total_num_face = total_num_tris + total_num_quads;
    if(!total_num_face || !num_nodes){
      QMessageBox::warning(this, tr("Application"),
                           tr("Error in reading file %1\n")
                           .arg(fileName));
      return false;
    }
     

    pos.resize(num_nodes);
    
    QString text_line;
    const char* Text_Line;
    double dc0 = 0.0;

    text_line=in.readLine();
    for(int i = 0; i < num_nodes; i++){
      text_line= in.readLine();
      Text_Line = text_line.toStdString().c_str();
      vector3d<double> p;
      Read_Flag = sscanf (Text_Line, "%lf %lf %lf %lf %lf",
                          &p.x,
                          &p.y,
                          &p.z,
                          &dc0,
                          &dc0);
     
      pos[i] = p;
    }



    //if .name file exists, read in and read mesh   
    QString name_file = surfFileName.left(surfFileName.lastIndexOf('.'))+".names";
    QFile infile(name_file);
    if(infile.exists()){
      if (!infile.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(name_file)
                           .arg(infile.errorString()));
        return false;
      }
      
      QTextStream intext(&infile);
      vector<int> ids;
      while(!intext.atEnd()){
        QStringList id_name=intext.readLine().split(' ', QString::SkipEmptyParts);
        if(id_name.size() < 2){
          break;
        }
        ids.push_back(id_name[0].toInt());
        boundary_names << id_name[1];
        
      }
      std::map<int, int> id_map;
      for(size_t i = 0; i < ids.size(); i++){
        id_map[ids[i]] = i;
      }
      infile.close();
      //read in mesh
     
      tmpmesh.resize(ids.size());
        
    //input surf_list
      for(int i = 0; i < total_num_tris; i++){
        int p1=0, p2=0, p3=0, id=0, flag1 = 0, flag2=0;
        text_line= in.readLine();
        Text_Line = text_line.toStdString().c_str();
      
        Read_Flag = sscanf (Text_Line, "%i  %i  %i  %i  %i %i",
                            &p1 ,
                            &p2 ,
                            &p3 ,
                            &id ,
                            &flag1,
                          &flag2);
        
        tmpmesh[id_map[id]].push_back(p1-1);
        tmpmesh[id_map[id]].push_back(p2-1);
        tmpmesh[id_map[id]].push_back(p3-1);
      }
      //input surf_list
      for(int i = 0; i < total_num_quads; i++){
        int p1=0, p2=0, p3=0, p4 =0, id=0, flag1 = 0, flag2=0;
        text_line= in.readLine();
        Text_Line = text_line.toStdString().c_str();
        
        Read_Flag = sscanf (Text_Line, "%i  %i  %i  %i  %i %i %i",
                          &p1 ,
                          &p2 ,
                          &p3 ,
                          &p4 ,
                          &id ,
                          &flag1,
                          &flag2);
     
      tmpmesh[id_map[id]].push_back(p1-1);
      tmpmesh[id_map[id]].push_back(p2-1);
      tmpmesh[id_map[id]].push_back(p3-1);

      tmpmesh[id_map[id]].push_back(p1-1);
      tmpmesh[id_map[id]].push_back(p3-1);
      tmpmesh[id_map[id]].push_back(p4-1);
      
      }
      mesh.push_back(tmpmesh);
      meshNodes.push_back(pos);
    }
  }//end of read in .surf and .names
  
  
  vector<positions3d> tmpMinMax;
  vector<bool> tmpVisible;
  int last_grid = mesh.size() -1;
  if(last_grid < 0)return false;

  //  currentGrid = last_grid;
  for (size_t i = 0; i < mesh[last_grid].size(); ++i) {
    
    positions3d minpos = meshNodes[last_grid][mesh[last_grid][i][0]];
    positions3d maxpos = minpos;
    for(size_t j = 0 ; j < mesh[last_grid][i].size(); ++j){ 
      positions3d t0 = meshNodes[last_grid][mesh[last_grid][i][j]];
      // Find extrema
      minpos.x = qMin(minpos.x, t0.x);
      
      minpos.y = qMin(minpos.y, t0.y);
      
      minpos.z = qMin(minpos.z, t0.z);
      
      maxpos.x = qMax(maxpos.x, t0.x);
      
      maxpos.y = qMax(maxpos.y, t0.y);
     
      maxpos.z = qMax(maxpos.z, t0.z);
      
    }
    
    tmpMinMax.push_back(minpos);
    tmpMinMax.push_back(maxpos);
    tmpVisible.push_back(true);
  }
  objMinMax.push_back(tmpMinMax);
  objVisible.push_back(tmpVisible);
  gridXform.push_back(affineMapping2());
  
  updateView(); 
  emit gridLoaded(boundary_names);
  
  
  addBoundObjects();

  
  resizeGL(currentWidth, currentHeight);
  
  return true;
  
}


 
    














/////////////////////////////////////////////////////////////////////////////
//  public:
//    void setVisibility(int i, bool);
//
//  The main window sends the boundary visibility model here. This function
//  extracts the information from the model and leaves it unchanged.
/////////////////////////////////////////////////////////////////////////////

 void MGViewer::setCurrentVisibility(const IDVisibility& idVis)
{
  int gid = idVis.gridId;
  int bid = idVis.boundId;
  if(gid >= (int)objVisible.size() || bid >= (int)objVisible[gid].size())return;
    
 
  
    objVisible[gid][bid] = idVis.show;
    updateGL();
  }



//////////////////////////////////////////////////////////////////////////////
//  protected:
//    void mousePressEvent(QMouseEvent *event);
//
//  When the mouse is pressed, this function stores which pixel was pressed.
//////////////////////////////////////////////////////////////////////////////

void MGViewer::mousePressEvent(QMouseEvent *event)
{
  //glLoadIdentity();
  int x = event->x();
  int y = event->y();
  double xmove, ymove, zmove;
  GLint viewport[4];
  GLdouble  projectionMatrix[16];
  GLdouble  modelviewMatrix[16];
  
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatrix);
  
  gluUnProject(x, viewport[3]-y-1, 0.0, modelviewMatrix,
               projectionMatrix, viewport, &xmove, &ymove, &zmove);
  positions3d P1 = positions3d(xmove, ymove, zmove);

  gluUnProject(x, viewport[3]-y-1, 1.0, modelviewMatrix,
               projectionMatrix, viewport, &xmove, &ymove, &zmove);
  positions3d P2 = positions3d(xmove, ymove, zmove);
  
  double t = (centerz-P1.z)/(P2.z-P1.z);
  xmove = P1.x*(1.0-t)+P2.x*t;
  ymove = P1.y*(1.0-t)+P2.y*t;
  
  lastPos = positions(xmove, ymove);
}

///////////////////////////////////////////////////////////////////////////
//  protected:
//    void mouseMoveEvent(QMouseEvent *event);
//
//  This function handles changes to the camera when the mouse is dragged
//  across the screen
///////////////////////////////////////////////////////////////////////////

void MGViewer::mouseMoveEvent(QMouseEvent *event)
{
  //glLoadIdentity();
  int x = event->x();
  int y = event->y();
  double xmove, ymove, zmove;
  GLint viewport[4];
  GLdouble  projectionMatrix[16];
  GLdouble modelviewMatrix[16]; 
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);
  gluUnProject(x, viewport[3]-y-1, 0.0, modelviewMatrix,
               projectionMatrix, viewport, &xmove, &ymove, &zmove);
  positions3d P1 = positions3d(xmove, ymove, zmove);
  gluUnProject(x, viewport[3]-y-1, 1.0, modelviewMatrix,
               projectionMatrix, viewport, &xmove, &ymove, &zmove);
  positions3d P2 = positions3d(xmove, ymove, zmove);
  
  double t = (centerz-P1.z)/(P2.z-P1.z);
  xmove = P1.x*(1.0-t)+P2.x*t;
  ymove = P1.y*(1.0-t)+P2.y*t;

  // gluUnProject(x, viewport[3]-y-1, 0.93, modelviewMatrix,
  //            projectionMatrix, viewport, &xmove, &ymove, &zmove);

  
   double dx = xmove - lastPos.x;
   double dy = ymove - lastPos.y;

   double zbar = 2*size - toz;
   
   if (event->buttons() & Qt::RightButton) {//translate
     tox += dx/scale;
     toy += dy/scale;
   }
   else if (event->buttons() & Qt::MidButton) {//zoom
     if (zbar > 0.05*size) {
       toz += dx*zbar/size;
       scale = 1.0;
     }
     else if (scale >= 1.0) {
       scale += 0.1*dx/(scale*size);
     }
     else {
       toz = 2*size - zbar - 0.01*size;
     }
   } else if (event->buttons() & Qt::LeftButton) {//rotate
     if(fabs(dx) < fabs(dy)){
       if(dy<0) rox += 1;
       else rox -= 1;
     }else{
       if(dx>0) roy += 1;
       else roy -= 1;
     }
   }

  
   updateGL();
   lastPos = positions(xmove, ymove);
   
}



//////////////////////////////////////////////////////////////////////////////
//  protected:
//    void mouseDoubleClickEvent(QMouseEvent *event);
//
//  This function handles changes to the camera when a mouse button is double
//  clicked.
//////////////////////////////////////////////////////////////////////////////

void MGViewer::mouseDoubleClickEvent(QMouseEvent *event)
{

  if (event->buttons() & Qt::LeftButton) {
    int x = event->x();
    int y = event->y();
    double xmove, ymove, zmove;
    GLint viewport[4];
    GLdouble  projectionMatrix[16];
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);

    gluUnProject(x, viewport[3]-y-1, 0.0, modelviewMatr,
                 projectionMatrix, viewport, &xmove, &ymove, &zmove);
    positions3d P1 = positions3d(xmove, ymove, zmove);

    gluUnProject(x, viewport[3]-y-1, 1.0, modelviewMatr,
                 projectionMatrix, viewport, &xmove, &ymove, &zmove);
    positions3d P2 = positions3d(xmove, ymove, zmove);
    

    
    vector<std::pair<double, pair<int, int> > > picked;
    for(unsigned int gid = 0; gid<mesh.size(); gid++){
      for (unsigned int bid = 0; bid < mesh[gid].size(); ++bid) {
        if (objVisible[gid][bid]) {
          for(size_t i = 0; i <mesh[gid][bid].size()/3; i++){ 
              
        
            positions3d t1 =  gridXform[gid].Map( meshNodes[gid][mesh[gid][bid][i*3 + 0]]); 
            positions3d t2 = gridXform[gid].Map(meshNodes[gid][mesh[gid][bid][i*3 + 1]]); 
            positions3d  t3 = gridXform[gid].Map(meshNodes[gid][mesh[gid][bid][i*3 + 2]]);
           
            positions3d hitP;
            if(CheckLineTri(t1, t2, t3, P1, P2, hitP)){
              double winx, winy, winz;
              gluProject(hitP.x, hitP.y, hitP.z, modelviewMatr,
                         projectionMatrix, viewport, &winx, &winy, &winz);
              picked.push_back(std::make_pair(winz, make_pair(gid, bid)));
              break;
                             
       
            }
          }
        }
      }
    }
    
    if(picked.size() > 1) std::sort(picked.begin(), picked.end());
    if(picked.size() >0) {
    
      emit pickCurrent(IDOnly(picked[0].second)); 
    }
          
        
          
    
  } else if (event->buttons() & Qt::MidButton) {
    
    
  } else if (event->buttons() & Qt::RightButton) {
    //    tox = toy = toz =rox = roy = roz = 0;  
    
  }
  
  updateGL();
  
}

//////////////////////////////////////////////////////////////////////////////
//  protected:
//    void wheelEvent(QWheelEvent *event);
//
//  This function handles changes to the camera when the mouse wheel is spun.
//////////////////////////////////////////////////////////////////////////////

void MGViewer::wheelEvent(QWheelEvent *event)
{
  roz += event->delta() / 8.0;
  updateGL();
}



//////////////////////////////////////////////////////////////////////////////
//  private:
//    void updateView();
//
//  This function updates the camera position the fit all the visible
//  boundaries in the view at once.
//////////////////////////////////////////////////////////////////////////////

void MGViewer::updateView()
{
  // Calculate global min/max pos of visible boundaries
  positions3d minpos, maxpos;
  bool first = true;
  for(size_t k = 0; k <  objMinMax.size(); ++k){
    for (size_t i = 0; i < objMinMax[k].size()/2; ++i) {
      if (objVisible[k][i]) {
        if (first) {
          minpos =  gridXform[k].Map(objMinMax[k][2*i]);
          maxpos = gridXform[k].Map(objMinMax[k][2*i+1]);
          first = false;
        } else {
          positions3d tmpMin =  gridXform[k].Map(objMinMax[k][2*i + 0]);
          positions3d tmpMax =  gridXform[k].Map(objMinMax[k][2*i + 1]);
          
          minpos.x = qMin(minpos.x, tmpMin.x);
          minpos.y = qMin(minpos.y, tmpMin.y);
          minpos.z = qMin(minpos.z, tmpMin.z);
          
          maxpos.x = qMax(maxpos.x, tmpMax.x);
          maxpos.y = qMax(maxpos.y, tmpMax.y);
          maxpos.z = qMax(maxpos.z, tmpMax.z);
        }
      }
    }
  }

  // Calculate size of boundary region
  size = sqrt((maxpos.x - minpos.x)*(maxpos.x - minpos.x) + 
              (maxpos.y - minpos.y)*(maxpos.y - minpos.y) + 
              (maxpos.z - minpos.z)*(maxpos.z - minpos.z));
  
  // Calculate center of boundary region
  centerx = (maxpos.x + minpos.x) / 2.0;
  centery = (maxpos.y + minpos.y) / 2.0;
  centerz = (maxpos.z + minpos.z) / 2.0;
  tox = toy = toz = rox = roy = roz = 0.0;
  
}







 

void MGViewer::showBoundaries(){
 //  bool checked = true;
//   for(unsigned int j=0; j < objVisible[currentGrid].size(); j++){
//       checked = objVisible[currentGrid][j] && checked;
//     }

//     for(unsigned int i=0; i < objVisible[currentGrid].size(); i++){
      
//     objVisible[currentGrid][i] = !checked;
//   }

//   updateGL();
}

void MGViewer::drawBoundObject(int gid, int bid, QColor c){
  if(meshNodes.empty() ||mesh.empty()){
    qDebug("please load in boundary first");
    return;
  }
  qglColor(c);
  for (size_t i = 0; i < mesh[gid][bid].size() / 3; ++i) {
    glBegin(GL_TRIANGLES);
    positions3d t0, t1, t2;
    t0 =gridXform[gid].Map(meshNodes[gid][mesh[gid][bid][i*3 + 0]]); 
    t1 = gridXform[gid].Map(meshNodes[gid][mesh[gid][bid][i*3 + 1]]); 
    t2 = gridXform[gid].Map(meshNodes[gid][mesh[gid][bid][i*3 + 2]]);
    glVertex3d(t0.x, t0.y, t0.z);
    glVertex3d(t1.x, t1.y, t1.z);
    glVertex3d(t2.x, t2.y, t2.z);
    glEnd();  
  }
}
 
void MGViewer::makeBoundWireframeObject(int gid, int bid, QColor c)
{
    // compile list
  GLuint newList = glGenLists(1);
  glNewList(newList, GL_COMPILE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
  drawBoundObject(gid, bid, c);   
  glEndList();
  
  GLuint oldList = boundObjects[gid][bid];
  boundObjects[gid][bid] = newList;
  glDeleteLists(oldList, 1);
}

void MGViewer::makeBoundFillObject(int gid, int bid, QColor c)
{
  // compile list
  GLuint newList = glGenLists(1);
  glNewList(newList, GL_COMPILE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);    
  drawBoundObject(gid, bid, c); 
  glEndList();
  GLuint oldList = boundObjects[gid][bid];
  boundObjects[gid][bid] = newList;
  glDeleteLists(oldList, 1);
}



void MGViewer::addBoundObjects()
{
  
  if(boundObjects.size() < mesh.size()){
    int first = boundObjects.size();
    int last = mesh.size() -1;
    for(int gid = first; gid <= last; ++gid){
      vector<GLuint> tmpList;
      for (size_t bid = 0; bid < mesh[gid].size(); ++bid) {
       
        // compile list
        GLuint newList = glGenLists(1);
         glNewList(newList, GL_COMPILE);
         glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
         drawBoundObject(gid, bid, default_color[bid%12]);   
         glEndList();
         tmpList.push_back(newList);
      }
       boundObjects.push_back(tmpList);
  
    }
  }
}
void MGViewer::transGrid(const IDMatrix& tc){
  int gid = tc.gridId;
  
  //gridXform[gid].translate(tc.translate) ;
  //gridXform[gid].rotateX(tc.rotate.x);
  //gridXform[gid].rotateY(tc.rotate.y);
  //gridXform[gid].rotateZ(tc.rotate.z);
  //gridXform[gid].scale(tc.scale);

  gridXform[gid] = tc.matrix;
  
  for (size_t bid = 0; bid < mesh[gid].size(); ++bid) {
    // compile list
    GLuint newList = glGenLists(1);
    glNewList(newList, GL_COMPILE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
    drawBoundObject(gid, bid, default_color[bid%12]);   
    glEndList();
    GLuint oldList = boundObjects[gid][bid];
    boundObjects[gid][bid] = newList;
    glDeleteLists(oldList, 1);
   
    
  }
   updateView();
   updateGL();
}

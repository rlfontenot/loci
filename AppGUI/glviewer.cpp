#include <QtOpenGL>
#include <QString>
#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <GL/glut.h>
#include <QTreeWidgetItem>
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


#include "glviewer.h"
#include "grid.h"
#include "hdf5.h"

#define PI 3.14159265358979323846264338327950
#include <stdio.h>
#include <stdlib.h>

#ifndef CALLBACK
#define CALLBACK
#endif



void CALLBACK errorCallback(GLenum errorCode)
{
   const GLubyte *estring;
   estring = gluErrorString(errorCode);
   fprintf(stderr, "Quadric Error: %s\n", estring);
   exit(0);
}


//////////////////////////////////////////////////////////////////////////////
//  Global Function
//
//  Places node index into the passed integar vector.
//////////////////////////////////////////////////////////////////////////////

void cbVertex2(void *vertex_data, void *user_data)
{
  vector<int>* tri = (vector<int>*) user_data;
  tri->push_back(*(int*)vertex_data);
}



///////////////////////////////////////////
//  public:
//    GLViewer(QWidget *parent = 0);
//
//  Initializes values.
///////////////////////////////////////////

GLViewer::GLViewer(QWidget *parent)
  : QGLWidget(parent)
{
  makeCurrent();
  fig = 0;
  gridObject = 0;
  contourObject = 0;
  borderObject = 0;
  shadingObject = 0;
  cpContourObject = 0;
  shadingObject = 0;
  qobj = 0;
  extreme_value = 0;
  centerx = centery = centerz=0;
  size = 0.1;
  currentObj = -1; //no select obj
  currentColor = default_color[0];
  tox= toy=toz = rox = roy = 0;
  show_contours = true;
  show_preview =  show_shading = show_boundary_shading = show_grid = show_border = show_extrema = false;
  scale = 1.0;
  shadeType = 1;
  min_val = max_val = 0.0;
  isFit = false;
  mode=BOUND_SELECT_MODE;
   show_shapes = true;
  
}

//////////////////////////////////////////////////////
//  public:
//    ~GLViewer();
//
//  Deletes display lists before widget destruction.
//////////////////////////////////////////////////////

GLViewer::~GLViewer()
{
  makeCurrent();
  if (gridObject){
    glDeleteLists(gridObject, 1);
    gridObject = 0;
  }
  if (contourObject){
    glDeleteLists(contourObject, 1);
    contourObject = 0;
  }
  if (borderObject){
    glDeleteLists(borderObject, 1);
    borderObject = 0;
  }
  if (shadingObject){
    glDeleteLists(shadingObject, 1);
    shadingObject = 0;
  }
  if( cpContourObject){
    glDeleteLists(cpContourObject, 1);
    cpContourObject = 0;
  }
  if (boundObjects.size() > 0) {
    for (size_t i = 0; i < boundObjects.size(); ++i)
      if(boundObjects[i])glDeleteLists(boundObjects[i], 1);
    boundObjects.clear();
  }
  if(fig){
    delete fig;
    fig =0;
  }

  if(qobj){
    gluDeleteQuadric(qobj);
    qobj=0;
  }
}




void GLViewer::clearCurrent(){
  if(currentObj!=-1){
    makeBoundWireframeObject(currentObj, currentColor);
    currentObj=-1;
  
    updateGL();
  }
  // show_nodes = false;
  //show_shapes = false;
  // tags.clear();
  //clear shapes

  
}





///////////////////////////////////////////
//  private:
//    void initializeGL();
//
//  Does some basic OpenGL configuration.
///////////////////////////////////////////

void GLViewer::initializeGL()
{
  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glClearColor(0.5, 0.5, 0.5, 0.5);
  qobj = gluNewQuadric();
  gluQuadricCallback(qobj, GLU_ERROR, 
                    0);
  gluQuadricDrawStyle(qobj, GLU_LINE); /* smooth shaded */
}

  

/////////////////////////////////////////////////////////////
//  protected:
//    void resizeGL(int height, int width);
//
//  Anytime the window is resized, this function handles it.
/////////////////////////////////////////////////////////////

void GLViewer::resizeGL(int width, int height)
{
  currentWidth = width;
  currentHeight = height;
  glViewport(0, 0, (GLint)width,(GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //change the parameters here also need change unproject parameters
  GLdouble near = 0.001*size;
  GLdouble far = size*30.0;
  gluPerspective(30.0, (GLdouble)width/height, near, far);
  glMatrixMode(GL_MODELVIEW);
}
////////////////////////////////////////////////////////////////////////////
//  protected:
//    void paintGL();
//
//  This function calls all the OpenGL commands and display lists to draw
//  whatever needs to be in the central widget.
////////////////////////////////////////////////////////////////////////////


void GLViewer::drawSphere(const vector<double>& p){
 
  if(p.size() < 4)return;
  glPushMatrix(); 
  glTranslated(p[0], p[1], p[2]);
  gluSphere(qobj, p[3], 20, 20);
  glPopMatrix();
}

void GLViewer::drawCylinder(const vector<double>& p){
 
  if(p.size()<5)return;
   
  glPushMatrix(); 
  glTranslated(p[0], p[1], p[3]);
  gluCylinder(qobj, p[2], p[2], p[4]-p[3], 20, 20);
  glPopMatrix(); 
}

void GLViewer::drawCone(const vector<double>& p){
  if(p.size()<6)return;
  glPushMatrix();

  double x0 = p[0];
  double y0 = p[1];
  double z0 = p[2];
  double z1 = p[4];
  double ratio = p[3];
  double z2 = p[5];
  
  glTranslated(x0, y0, z1);
  gluCylinder(qobj, ratio*(z1-z0), ratio*(z2-z0), z2-z1, 20, 20);
  glPopMatrix(); 
}


void GLViewer::drawCube(const vector<double>& p){
   if(p.size()<6)return;
  glPushMatrix(); 
  double x1 = p[0];
  double y1 = p[1];
  double z1 = p[2];
  double x2 = p[3];
  double y2 = p[4];
  double z2 = p[5];
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  

  glVertex3d(x1, y1, z1);
  glVertex3d(x2, y1, z1);
  glVertex3d(x2, y1, z2);
  glVertex3d(x1, y1, z2);


  glVertex3d(x1, y2, z1);
  glVertex3d(x1, y2, z2);
  glVertex3d(x2, y2, z2);
  glVertex3d(x2, y2, z1);
    

  glVertex3d(x2, y1, z1);
  glVertex3d(x2, y2, z1);
  glVertex3d(x2, y2, z2);
  glVertex3d(x2, y1, z2);


  glVertex3d(x1, y1, z1);
  glVertex3d(x1, y1, z2);
  glVertex3d(x1, y2, z2);
  glVertex3d(x1, y2, z1);


  glVertex3d(x1, y1, z2);
  glVertex3d(x2, y1, z2);
  glVertex3d(x2, y2, z2);
  glVertex3d(x1, y2, z2);
      

  glVertex3d(x1, y1, z1);
  glVertex3d(x1, y2, z1);
  glVertex3d(x2, y2, z1);
  glVertex3d(x2, y1, z1);
  glEnd();
  glPopMatrix();   
}

void GLViewer::drawPxPlane(const vector<double>& p, double size){
  if(p.size()==0)return;
  double x = p[0];
  glPushMatrix(); 
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d(x, -1*size, -1*size);
  glVertex3d(x, -1*size, size);
  glVertex3d(x, size, size);
  glVertex3d(x, size, -1*size);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(x, 0, 0);
  glVertex3d(x+0.5*size, 0, 0);
  glVertex3d(x, -1*size, 0);
  glVertex3d(x, size, 0);
  glVertex3d(x, 0, -1*size);
  glVertex3d(x, 0, size);
  glEnd();
  glPopMatrix(); 
}
void GLViewer::drawNxPlane(const vector<double>& p, double size){
 
  if(p.size()==0)return;
  double x = p[0];
  glPushMatrix(); 
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d(x, -1*size, -1*size);
  glVertex3d(x, -1*size, size);
  glVertex3d(x, size, size);
  glVertex3d(x, size, -1*size);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(x, 0, 0);
  glVertex3d(x-0.5*size, 0, 0);
  glVertex3d(x, -1*size, 0);
  glVertex3d(x, size, 0);
  glVertex3d(x, 0, -1*size);
  glVertex3d(x, 0, size);
  glEnd();

  glPopMatrix();
  
}

void GLViewer::drawPyPlane(const vector<double>& p, double size){
  if(p.size()==0)return;
  double y = p[0];
  glPushMatrix(); 

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d( -1*size,y,  -1*size);
  glVertex3d( -1*size, y,size);
  glVertex3d( size, y,size);
  glVertex3d( size, y,-1*size);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0, y, 0);
  glVertex3d(0, y+0.5*size,0);
  glVertex3d( -1*size, y, 0);
  glVertex3d( size,y, 0);
  glVertex3d( 0, y,-1*size);
  glVertex3d( 0,y, size);
  glEnd();

  glPopMatrix();
  
}

void GLViewer::drawNyPlane(const vector<double>& p, double size){
  if(p.size()==0)return;
  double y = p[0];
  glPushMatrix();
  //  glLineStipple(1, 0xAAAA);
  //glEnable(GL_LINE_STIPPLE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d( -1*size,y,  -1*size);
  glVertex3d( -1*size, y,size);
  glVertex3d( size, y,size);
  glVertex3d( size, y,-1*size);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0, y, 0);
  glVertex3d(0, y-0.5*size,0);
  glVertex3d( -1*size, y, 0);
  glVertex3d( size,y, 0);
  glVertex3d( 0, y,-1*size);
  glVertex3d( 0,y, size);
  glEnd();
   glPopMatrix();
  
}

void GLViewer::drawPzPlane(const vector<double>& p, double size){
if(p.size()==0)return;
  double z = p[0];
  glPushMatrix();
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d( -1*size, -1*size,z);
  glVertex3d( -1*size, size,z);
  glVertex3d( size, size,z);
  glVertex3d( size, -1*size,z);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0,  0, z);
  glVertex3d(0, 0, z + 0.5*size);
  glVertex3d( -1*size, 0, z);
  glVertex3d( size,0, z);
  glVertex3d( 0, -1*size, z);
  glVertex3d( 0,size, z);
  glEnd();
   glPopMatrix();
  
}

void GLViewer::drawNzPlane(const vector<double>& p, double size){
  if(p.size()==0)return;
  double z = p[0];
  glPushMatrix();
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d( -1*size, -1*size,z);
  glVertex3d( -1*size, size,z);
  glVertex3d( size, size,z);
  glVertex3d( size, -1*size,z);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0,  0, z);
  glVertex3d(0, 0, z - 0.5*size);
   glVertex3d( -1*size, 0, z);
  glVertex3d( size,0, z);
  glVertex3d( 0, -1*size, z);
  glVertex3d( 0,size, z);
  glEnd();
  glPopMatrix();
  
}

void GLViewer::cleanDoc(){//signal not working
  doc = QDomDocument();
  updateGL();
 
}
void GLViewer::updateDoc(const QTreeWidgetItem* root){
  doc = tree2dom(root);
  updateGL();
}





void GLViewer::drawShapes(){
  show_shapes =  true;
 
  if(doc.isNull())return;
  
 
  QDomNodeList shapes = doc.elementsByTagName("object");
  
  if(shapes.isEmpty())return;


  glLineWidth(5);
  for( int i = 0; i <shapes.size(); i++){
    QDomElement elt = shapes.at(i).firstChildElement();
    if(elt.isNull())continue; //no child, do nothing
    glPushMatrix();
    double planesize = 1;
   
    for(; !elt.isNull(); elt = elt.nextSiblingElement()){
   
   
      if(elt.tagName()== "transform"){
      
        for( QDomElement elm = elt.firstChildElement(); !elm.isNull(); elm= elm.nextSiblingElement()){
          if(elm.tagName() =="translate"){
            double x0 =0 , y0 = 0, z0 = 0;
            for(QDomElement trans_elem = elm.firstChildElement();!trans_elem.isNull(); trans_elem = trans_elem.nextSiblingElement()){
              if(trans_elem.tagName()=="x0")x0 = trans_elem.text().toDouble();
              else  if(trans_elem.tagName()=="y0")y0 = trans_elem.text().toDouble();
              else  if(trans_elem.tagName()=="z0")z0 = trans_elem.text().toDouble();
              else{
              
                qDebug()<<tr("illegal child ") + trans_elem.tagName() + tr(" in 'translate'");
                                   
                return;
              }
            
            } 
             
            if(x0!=0 || y0 !=0 || z0 !=0)
              glTranslated(x0,
                           y0,
                           z0) ;
          }else if(elm.tagName() =="scale"){
            double x0 =1 , y0 = 1, z0 = 1;
            for(QDomElement trans_elem = elm.firstChildElement();!trans_elem.isNull(); trans_elem = trans_elem.nextSiblingElement()){
              if(trans_elem.tagName()=="x0")x0 = trans_elem.text().toDouble();
              else  if(trans_elem.tagName()=="y0")y0 = trans_elem.text().toDouble();
              else  if(trans_elem.tagName()=="z0")z0 = trans_elem.text().toDouble();
              else{
                
                qDebug()<<tr("illegal child ") + trans_elem.tagName() + tr(" in 'scale'");
                                   
                return;
              }
            
            } 
             
            
        
          
            if(x0!=1 || y0 !=1 || z0 !=1)
              glScaled(x0,
                       y0,
                       z0) ;
            planesize = planesize*x0;
          }else if(elm.tagName() =="rotateX"){

            double theta = 0;
            for(QDomElement trans_elem = elm.firstChildElement();!trans_elem.isNull(); trans_elem = trans_elem.nextSiblingElement()){
              if(trans_elem.tagName()=="theta")theta = trans_elem.text().toDouble();
              else{
                
                qDebug()<<tr("illegal child ") + trans_elem.tagName() + tr(" in 'rotateX'");
                
                return;
              }
              
          
            }
             
            if(theta != 0)glRotated(theta, 1.0, 0.0, 0.0);
          }else if(elm.tagName() =="rotateY"){

            double theta = 0;
            for(QDomElement trans_elem = elm.firstChildElement();!trans_elem.isNull(); trans_elem = trans_elem.nextSiblingElement()){
              if(trans_elem.tagName()=="theta")theta = trans_elem.text().toDouble();
              else{
                
                qDebug()<<tr("illegal child ") + trans_elem.tagName() + tr(" in 'rotateY'");
                
                return;
              }
              
          
            }
            
            if(theta != 0)glRotated(theta, 0.0, 1.0, 0.0);
            

          }else if(elm.tagName() =="rotateZ"){

            double theta = 0;
            for(QDomElement trans_elem = elm.firstChildElement();
                !trans_elem.isNull();
                trans_elem = trans_elem.nextSiblingElement()){
              if(trans_elem.tagName()=="theta")theta = trans_elem.text().toDouble();
              else{
                
                qDebug()<<tr("illegal child ") + trans_elem.tagName() + tr(" in 'rotateZ'");
                
                return;
              }
              
          
            }
            
            if(theta != 0)glRotated(theta, 0.0, 0.0, 1.0);
            
            
          }else{
             
              qDebug()<< tr("illegal child ") + elm.tagName() + tr(" in 'transform'");
              return;
              
          }
        }
      }else if(elt.tagName()=="shape"){
        QDomElement elm = elt.firstChildElement();
    
        vector<double> para;
        for(QDomElement obj_elem = elm.firstChildElement(); !obj_elem.isNull(); obj_elem = obj_elem.nextSiblingElement()){
            para.push_back(obj_elem.text().toDouble());
        }
    
       
        if(elm.tagName() =="sphere") 
          drawSphere(para);
        else if(elm.tagName() =="cone") 
          drawCone(para);
        else if(elm.tagName() =="cylinder") 
          drawCylinder(para); 
        else if(elm.tagName() =="box") 
          drawCube(para);
        else if(elm.tagName() =="x_plus_plane")
          drawPxPlane(para, planesize);
        else if(elm.tagName() =="x_minus_plane")
          drawNxPlane(para, planesize);
        else if(elm.tagName() =="y_plus_plane")
          drawPyPlane(para, planesize);
        else if(elm.tagName() =="y_minus_plane")
          drawNyPlane(para, planesize);
        else if(elm.tagName() =="z_plus_plane")
          drawPzPlane(para, planesize);
        else if(elm.tagName() =="z_minus_plane")
          drawNzPlane(para, planesize);
        
        

        glPopMatrix(); 
    
      }
    }
  }
  glLineWidth(1);
  drawMarkedNodes();
    
}
void GLViewer::markNodes(){
  if(doc.firstChildElement().isNull())return;
  tags.clear();
  QDomElement rootElement = doc.firstChildElement("region");
  tags = process_region(rootElement, meshNodes);
  updateGL();
}


unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

void GLViewer::markVolumeNodes(QString filename){
  //read in nodes
  hid_t input_fid ; 
  input_fid = H5Fopen(filename.toLocal8Bit(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(input_fid <= 0) {
    qDebug() << "unable to open file '" << filename << "'"<< endl ;
    return;
  }
  // read in positions
  hid_t fi = H5Gopen(input_fid,"file_info", H5P_DEFAULT) ;
  unsigned long numNodes = readAttributeLong(fi,"numNodes") ;

  
  H5Gclose(fi) ;
  
  hsize_t count = numNodes ;
  
  //#ifdef H5_INTERFACE_1_6_4
  hsize_t lstart = 0 ;
  //#else
    // hssize_t lstart = 0 ;
  //#endif
  
    // Read in pos data from file i
  vector<positions3d> pos_dat(numNodes) ;
  hid_t node_g = H5Gopen(input_fid,"node_info", H5P_DEFAULT) ;
  hid_t dataset = H5Dopen(node_g,"positions", H5P_DEFAULT) ;
  hid_t dspace = H5Dget_space(dataset) ;

  

  hid_t pos_tid = H5Tcreate(H5T_COMPOUND, sizeof(positions3d));
    
  H5Tinsert(pos_tid, "x", 0, H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "y", sizeof(double), H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "z", 2*sizeof(double), H5T_IEEE_F64LE);

  hsize_t stride = 1 ;
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count, NULL
) ;
  int rank = 1 ;
  hsize_t dimension = count ;
  hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
 
  hid_t err = H5Dread(dataset,pos_tid, memspace,dspace,H5P_DEFAULT,
		      &pos_dat[0]) ;
  if(err < 0) {
    qDebug() << "unable to read positions from '" << filename << "'" << endl ;
    return ;
  }
  H5Sclose(dspace) ;
  H5Dclose(dataset) ;
  H5Gclose(node_g) ;

  //mark the nodes
 
  if(doc.firstChildElement().isNull())return;

  
  
  QDomElement rootElement = doc.firstChildElement("region");
  vector<bool> vtags = process_region(rootElement, pos_dat);
  QString tagFileName =filename.section('.', 0, -2)+".tag";
  
  tagFileName = QFileDialog::getSaveFileName(this, tr("Save .tag File"),
                                             tagFileName,
                                             tr("tag Files (*.tag)"));
  if(tagFileName.section('.', -1, -1)!="tag") tagFileName +=".tag";
  



  QFileInfo tagInfo(tagFileName);
  if(tagInfo.exists()){
    QString command = "mv " + tagFileName+ " " + tagFileName+".bak";
    int ret =  system(command.toStdString().c_str());
    if(!WIFEXITED(ret))
      {
        if(WIFSIGNALED(ret))
          {
            QMessageBox::information(window(), "save tag file",
                                     command + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
            return;
          }
      }
  }
  if(tagFileName.isNull()){
    
     return ;
  }
  
  QFile file(tagFileName);
  if (!file.open(QFile::WriteOnly | QFile::Text)) {
    QMessageBox::warning(this, tr("save .vars file "),
                         tr("Cannot write file %1:\n%2.")
                         .arg(tagFileName)
                         .arg(file.errorString()));
    return;
  }

  
  QTextStream out(&file);
  for(unsigned int i = 0; i < vtags.size(); i++)
    out<< vtags[i]<<endl;
  file.close();
}

void GLViewer::drawMarkedNodes(){
    if(tags.size() != meshNodes.size()) return;
    if(tags.size()==0) return;
    glPointSize(5);
    glBegin(GL_POINTS);
    for(unsigned int i = 0; i < tags.size(); i++){
      if(tags[i]){
        glVertex3d(meshNodes[i].x, meshNodes[i].y, meshNodes[i].z);
      }
  }
     glEnd();
  glPointSize(1);
}

void GLViewer::drawExtremeNodes(double value){
  
  
  
  double start=min_val, end = max_val;
  

  if(value > mid_val){
    start = value;
    end  = max_val;
  }else{
    start = min_val;
    end = value;
  }
  
  if(extremeNodes.size() != extremeValues.size()) return;
  if(extremeNodes.size() == 0) return;
 
  glPointSize(5);
  glBegin(GL_POINTS);
     
  for(unsigned int i = 0; i < extremeNodes.size(); i++){
    if(extremeValues[i] >= start && extremeValues[i]<=end){
      positions3d c = shade(extremeValues[i]);
      glColor3d(c.x, c.y, c.z);
      glVertex3d(extremeNodes[i].x, extremeNodes[i].y, extremeNodes[i].z);
    }
  }
  
  glEnd();
  glPointSize(1);
}


  
void GLViewer::paintGL()
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


  glLineWidth(1);
  switch (mode) {
    
  case BOUND_SELECT_MODE: //after the grid is load, only draw boundaries
    for (size_t i = 0; i < boundObjects.size(); ++i)
      if (objVisible[i])glCallList(boundObjects[i]);
    if(show_preview){
      glEnable(GL_LINE_SMOOTH);
      if(cpContourObject)glCallList(cpContourObject);
      glDisable(GL_LINE_SMOOTH);
    }
    //drawBoxes();
    //drawCoordGrid();

    break;

  case PLANE_AND_BOUND_MODE:
  
    for (size_t i = 0; i < boundObjects.size(); ++i)
      if (objVisible[i])glCallList(boundObjects[i]);
    if (show_shading && shadingObject>0)glCallList(shadingObject);
    if (show_grid && gridObject > 0)glCallList(gridObject);
    glEnable(GL_LINE_SMOOTH);
    if(show_border && borderObject > 0)glCallList(borderObject);
    glDisable(GL_LINE_SMOOTH);
    if (show_contours && contourObject > 0)glCallList(contourObject);
    break;
     
 
  default:
    break;
  }
 
  glColor3f(0.0f, 0.0f, 0.0f);  
  //  if(show_nodes) drawMarkedNodes();
  if(show_extrema) drawExtremeNodes(extreme_value);
  if(show_shapes)drawShapes();
  glPopMatrix();

  glFlush();
  
}
void GLViewer::setCurrentObj(int i, QColor c){
  if( i>=0 && i<(int)boundObjects.size()){
    if(currentObj != -1) makeBoundWireframeObject(currentObj, currentColor);
      makeBoundFillObject(i, c);
      currentObj =i;
      currentColor = c;
      updateGL();
      
  }
}

void GLViewer::reset(){
  
  tox= toy=toz=rox=roy =roz= 0;
  if(isFit){
    for(unsigned int i =0; i < objVisible.size(); i++){
      objVisible[i] = true;
    }
    isFit = false;
    updateView();
  }
  updateGL();
}

void GLViewer::fit(){
  if(currentObj<0){
     QMessageBox::warning(window(), tr("fit"),
                      tr("please select a boundary first")
                       );
    return;
  }
  isFit = true;
  
  for( int i =0; i < (int)objVisible.size(); i++){
    if(i == currentObj) objVisible[i] = true;
    else objVisible[i] = false;
  }
  updateView();
  updateGL();
}




struct surface_info {
  vector<int> trias ;
  vector<int> quads ;
  vector<vector<int> > gen_faces ;
} ;


////////////////////////////////////////
//  public:
//    QSize sizeHint() const;
//
//  Requests a default size of 800x600.
////////////////////////////////////////

QSize GLViewer::sizeHint() const
{
  return QSize(600, 600);
}


void GLViewer::setLoadInfo(const LoadInfo& ld_info){
 
   loadInfo.casename = ld_info.casename;
   loadInfo.directory = ld_info.directory;
   loadInfo.iteration = ld_info.iteration;
   loadInfo.variable = ld_info.variable;
}

bool GLViewer::load_boundary(QString fileName,  QStringList& boundary_names) {
  reset();
  boundary_names.clear();  
  
  //assume this is vog file
  if(fileName.right(4)!=".vog"){
     QMessageBox::warning(window(), "load grid",
                         tr("the format of grid is not .vog" ));
     return false;
  }

  QString surfFileName = fileName.section('.', 0, -2)+".surface";
  QFileInfo surfInfo(surfFileName);
  QFileInfo vogInfo(fileName);
  
  if(!(surfInfo.exists())|| surfInfo.created() < vogInfo.created()){
    QString command2 = "vog2surf -surface " + surfFileName + " " + fileName.section('.', 0, -2);
    int ret =  system(command2.toStdString().c_str());
    if(!WIFEXITED(ret))
      {
        if(WIFSIGNALED(ret))
          {
            QMessageBox::information(window(), "load grid",
                                     command2 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );

              return false;
            }
      }
  }

  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
  loadInfo.casename = casename;
  loadInfo.directory = directory;

  if(surfFileName.right(8) ==".surface"){

    QFile file(surfFileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Application"),
                           tr("Cannot read file %1:\n%2.")
                           .arg(fileName)
                           .arg(file.errorString()));
      return false;
    }
    
    QTextStream in(&file); 

    mesh.clear();
    vector<vector3d<double> > pos;
    
    
    size_t  npos = 0;
    
  
    in >> npos ;
    pos.resize(npos);
    //input pos
    for(size_t i=0;i<npos;++i) {
      vector3d<double> p;
      in >> p.x >> p.y >> p.z ;
      pos[i] = p;
    }

  //input surf_list
  size_t nsurf = 0;
  in >> nsurf;
  in.readLine();
  boundary_names.clear();
  for(size_t i = 0; i < nsurf; i++){
    QString name=in.readLine();
    boundary_names << name;
  }
 
  
  if(nsurf<=0) return false;
  vector<surface_info> surf_list(nsurf);
  QList<int> bids;
  for(size_t i = 0; i < nsurf; i++){
    size_t ntris=0;
    size_t nquads = 0;
    size_t ngens = 0;
    int bid;
    in >>bid>>ntris >> nquads>>ngens ;
    bids<<bid;
    ntris *=3;
    nquads *=4;
    if(ntris>0) surf_list[i].trias.resize(ntris);
    if(nquads>0) surf_list[i].quads.resize(nquads);
    if(ngens>0)surf_list[i].gen_faces.resize(ngens);
    if(ntris>0){
      for(size_t j =0; j < ntris; j++){
        in>> surf_list[i].trias[j];
       
      }
    }
    

    if(nquads>0){
      for(size_t j=0;j<nquads;++j){
        in>>surf_list[i].quads[j];
      }
    }
    
    if(ngens >0){
      for(size_t j=0;j<ngens;++j){
        size_t nf;
        in >> nf ;
        vector<int> gen;
        if(nf >0){
          gen.resize(nf);
          for(size_t k=0;k<nf;++k)
            in>> gen[k] ;
          
        }
        surf_list[i].gen_faces[j] = gen;
      }
    }
  }

  // set up GLUtesselator
  GLUtesselator* myTess = gluNewTess();
  gluTessCallback(myTess, GLU_TESS_VERTEX_DATA,
                  (GLvoid (*) ()) &cbVertex2);
  gluTessCallback(myTess, GLU_TESS_EDGE_FLAG,
                  (GLvoid (*) ()) &cbEdgeFlag);
  
 // load node index vector
  vector<int>* pntIndex = new vector<int>;
  pntIndex->resize(pos.size()+1);
  for (unsigned int i = 0; i <= pos.size(); ++i)
    (*pntIndex)[i] = i;
  
  for(unsigned int id = 0; id <nsurf ; id++){
    size_t ntris=surf_list[id].trias.size()/3;
    size_t nquads = surf_list[id].quads.size()/4;
    //size_t ngens = surf_list[id].gen_faces.size();
    
    vector<int> vTri;
    for(size_t i = 0; i < ntris*3; i++)
      vTri.push_back(surf_list[id].trias[i]);
   // quads
    
    for (size_t j = 0; j < nquads; ++j) {
      gluTessBeginPolygon(myTess, &vTri);
      gluTessBeginContour(myTess);
      for (int k = 0; k < 4; ++k) {
        GLdouble point[3];
        point[0] = pos[surf_list[id].quads[j*4+k]-1].x;
        point[1] = pos[surf_list[id].quads[j*4+k]-1].y;
        point[2] = pos[surf_list[id].quads[j*4+k]-1].z;
        gluTessVertex(myTess, point, &(*pntIndex)[(surf_list[id].quads)[j*4+k]] );
      }
      gluTessEndContour(myTess);
      gluTessEndPolygon(myTess);
    }
       
    
    // gen cells

     
    for (size_t j = 0; j < surf_list[id].gen_faces.size(); ++j) {
      gluTessBeginPolygon(myTess, &vTri);
      gluTessBeginContour(myTess);
      for (size_t k = 0; k < surf_list[id].gen_faces[j].size(); ++k) {
        GLdouble point[3];
        point[0] = pos[surf_list[id].gen_faces[j][k]-1].x;
          point[1] = pos[surf_list[id].gen_faces[j][k]-1].y;
          point[2] = pos[surf_list[id].gen_faces[j][k]-1].z;
          gluTessVertex(myTess, point, &(*pntIndex)[(surf_list[id].gen_faces)[j][k]]);
      }
      gluTessEndContour(myTess);
      gluTessEndPolygon(myTess);
      
    }
    for(size_t ii = 0; ii < vTri.size(); ii++){
      vTri[ii] = vTri[ii]-1;
    }
    
      mesh.push_back(vTri);
 
  }//for(bid..)
  meshNodes.clear();
  meshNodes=vector<positions3d>(pos);
 
  gluDeleteTess(myTess);
  delete  pntIndex;
  meshMap.clear();
  meshMap.resize(npos);
  for (size_t  i = 0; i < npos; i++){
    in >>meshMap[i];
     }
  file.close();
 
  }else{}
 
  
  // Remap mesh to match nodes
  objMinMax.clear();
  objVisible.clear();
  
  for (size_t i = 0; i < mesh.size(); ++i) {
    
    positions3d minpos = meshNodes[mesh[i][0]];
    positions3d maxpos = minpos;
    for(size_t j = 0 ; j < mesh[i].size(); ++j){ 
      positions3d t0 = meshNodes[mesh[i][j]];
      // Find extrema
      minpos.x = qMin(minpos.x, t0.x);
      
      minpos.y = qMin(minpos.y, t0.y);
      
      minpos.z = qMin(minpos.z, t0.z);
      
      maxpos.x = qMax(maxpos.x, t0.x);
      
      maxpos.y = qMax(maxpos.y, t0.y);
     
      maxpos.z = qMax(maxpos.z, t0.z);
      
    }
    
    objMinMax.push_back(minpos);
    objMinMax.push_back(maxpos);
    objVisible.push_back(true);
  }
 
  
  //finish reading in all information
 
  updateView(); 
  extremeValues.clear();
  extremeNodes.clear();
  show_shapes = true;
  resizeGL(currentWidth, currentHeight);
  mode = BOUND_SELECT_MODE;
  currentObj=-1;
  makeObjects();
  updateGL();
  //  cleanDoc();
  return true;
}


 
bool GLViewer::load_boundary(QString fileName,  QStringList& boundary_names, QList<int>& bids) {
  reset();
  
  boundary_names.clear();
  bids.clear();
  //assume this is vog file
  if(fileName.right(4)!=".vog"){
     QMessageBox::warning(window(), "load grid",
                         tr("the format of grid is not .vog" ));
     return false;
  }

  QString surfFileName = fileName.section('.', 0, -2)+".surface";
  QFileInfo surfInfo(surfFileName);
  QFileInfo vogInfo(fileName);
  
  
  if(!(surfInfo.exists())|| surfInfo.created() < vogInfo.created()){
    QString command2 = "vog2surf -surface " + surfFileName + " " + fileName.section('.', 0, -2);
    int ret =  system(command2.toStdString().c_str());
    if(!WIFEXITED(ret))
      {
        if(WIFSIGNALED(ret))
          {
            QMessageBox::information(window(), "load grid",
                                     command2 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );

              return false;
            }
      }
  }

  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
  loadInfo.casename = casename;
  loadInfo.directory = directory;
 
  

  

  if(surfFileName.right(8) ==".surface"){

    QFile file(surfFileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Application"),
                           tr("Cannot read file %1:\n%2.")
                           .arg(fileName)
                           .arg(file.errorString()));
      return false;
    }
    
    QTextStream in(&file); 
  
    
    mesh.clear();
    vector<vector3d<double> > pos;
    
    
    size_t  npos = 0;
    
  
    in >> npos ;
    pos.resize(npos);
    //input pos
    for(size_t i=0;i<npos;++i) {
      vector3d<double> p;
      in >> p.x >> p.y >> p.z ;
      pos[i] = p;
    }
  
  //input surf_list
  size_t nsurf = 0;
  in >> nsurf;
  in.readLine();
  boundary_names.clear();
  for(size_t i = 0; i < nsurf; i++){
    QString name=in.readLine();
    boundary_names << name;
  }
 
  
  if(nsurf<=0) return false;
  vector<surface_info> surf_list(nsurf);
 
  for(size_t i = 0; i < nsurf; i++){
    size_t ntris=0;
    size_t nquads = 0;
    size_t ngens = 0;
    int bid;
    in >>bid>>ntris >> nquads>>ngens ;
    bids<<bid;
    ntris *=3;
    nquads *=4;
    if(ntris>0) surf_list[i].trias.resize(ntris);
    if(nquads>0) surf_list[i].quads.resize(nquads);
    if(ngens>0)surf_list[i].gen_faces.resize(ngens);
    if(ntris>0){
      for(size_t j =0; j < ntris; j++){
        in>> surf_list[i].trias[j];
       
      }
    }
    

    if(nquads>0){
      for(size_t j=0;j<nquads;++j){
        in>>surf_list[i].quads[j];
      }
    }
    
    if(ngens >0){
      for(size_t j=0;j<ngens;++j){
        size_t nf;
        in >> nf ;
        vector<int> gen;
        if(nf >0){
          gen.resize(nf);
          for(size_t k=0;k<nf;++k)
            in>> gen[k] ;
          
        }
        surf_list[i].gen_faces[j] = gen;
      }
    }
  }
 

  // set up GLUtesselator
  GLUtesselator* myTess = gluNewTess();
  gluTessCallback(myTess, GLU_TESS_VERTEX_DATA,
                  (GLvoid (*) ()) &cbVertex2);
  gluTessCallback(myTess, GLU_TESS_EDGE_FLAG,
                  (GLvoid (*) ()) &cbEdgeFlag);
  
 // load node index vector
  vector<int>* pntIndex = new vector<int>;
  pntIndex->resize(pos.size()+1);
  for (unsigned int i = 0; i <= pos.size(); ++i)
    (*pntIndex)[i] = i;
  
  for(unsigned int id = 0; id <nsurf ; id++){
    size_t ntris=surf_list[id].trias.size()/3;
    size_t nquads = surf_list[id].quads.size()/4;
    //size_t ngens = surf_list[id].gen_faces.size();
    
    vector<int> vTri;
    for(size_t i = 0; i < ntris*3; i++)
      vTri.push_back(surf_list[id].trias[i]);
   // quads
    
    for (size_t j = 0; j < nquads; ++j) {
      gluTessBeginPolygon(myTess, &vTri);
      gluTessBeginContour(myTess);
      for (int k = 0; k < 4; ++k) {
        GLdouble point[3];
        point[0] = pos[surf_list[id].quads[j*4+k]-1].x;
        point[1] = pos[surf_list[id].quads[j*4+k]-1].y;
        point[2] = pos[surf_list[id].quads[j*4+k]-1].z;
        gluTessVertex(myTess, point, &(*pntIndex)[(surf_list[id].quads)[j*4+k]] );
      }
      gluTessEndContour(myTess);
      gluTessEndPolygon(myTess);
    }
       
    
    // gen cells

     
    for (size_t j = 0; j < surf_list[id].gen_faces.size(); ++j) {
      gluTessBeginPolygon(myTess, &vTri);
      gluTessBeginContour(myTess);
      for (size_t k = 0; k < surf_list[id].gen_faces[j].size(); ++k) {
        GLdouble point[3];
        point[0] = pos[surf_list[id].gen_faces[j][k]-1].x;
          point[1] = pos[surf_list[id].gen_faces[j][k]-1].y;
          point[2] = pos[surf_list[id].gen_faces[j][k]-1].z;
          gluTessVertex(myTess, point, &(*pntIndex)[(surf_list[id].gen_faces)[j][k]]);
      }
      gluTessEndContour(myTess);
      gluTessEndPolygon(myTess);
      
    }
    for(size_t ii = 0; ii < vTri.size(); ii++){
      vTri[ii] = vTri[ii]-1;
    }
    
      mesh.push_back(vTri);
      
  }//for(bid..)
  meshNodes.clear();
  meshNodes=vector<positions3d>(pos);
 
  gluDeleteTess(myTess);
  delete  pntIndex;
  meshMap.clear();
  meshMap.resize(npos);
  for (size_t  i = 0; i < npos; i++){
    in >>meshMap[i];
     }
  file.close();
  
  }else{}
    
  
  // Remap mesh to match nodes
  objMinMax.clear();
  objVisible.clear();
  
  for (size_t i = 0; i < mesh.size(); ++i) {
    
    positions3d minpos = meshNodes[mesh[i][0]];
    positions3d maxpos = minpos;
    for(size_t j = 0 ; j < mesh[i].size(); ++j){ 
      positions3d t0 = meshNodes[mesh[i][j]];
      // Find extrema
      minpos.x = qMin(minpos.x, t0.x);
      
      minpos.y = qMin(minpos.y, t0.y);
      
      minpos.z = qMin(minpos.z, t0.z);
      
      maxpos.x = qMax(maxpos.x, t0.x);
      
      maxpos.y = qMax(maxpos.y, t0.y);
     
      maxpos.z = qMax(maxpos.z, t0.z);
      
    }
    
    objMinMax.push_back(minpos);
    objMinMax.push_back(maxpos);
    objVisible.push_back(true);
  }
 
  
  //finish reading in all information
  
  updateView(); 
 
  
 
  
  
  extremeValues.clear();
  extremeNodes.clear();
  
  show_shapes = true;
  resizeGL(currentWidth, currentHeight);
   mode = BOUND_SELECT_MODE;
   currentObj=-1;
   makeObjects();
  //  clearCurrent();
  cleanDoc();
  updateGL();
  return true;
}
    














/////////////////////////////////////////////////////////////////////////////
//  public:
//    void setVisibility(int i, bool);
//
//  The main window sends the boundary visibility model here. This function
//  extracts the information from the model and leaves it unchanged.
/////////////////////////////////////////////////////////////////////////////

void GLViewer::setVisibility(int i, bool show)
{
  
  if(i<0 || i>(int)objVisible.size()){
    qDebug() << "index out of range in GLViewer::setVisibility()";
  }else{
  
    objVisible[i] = show;
    updateGL();
  }
}



///////////////////////////////////////////////////////////////////////////
// 
//
//////////////////////////////////////////////////////////////////////////
void GLViewer::uncut(){
 
  for(unsigned int i=0; i < objVisible.size(); i++){
    objVisible[i] = true;
  }
  mode = BOUND_SELECT_MODE;
  show_preview = false;
  extremeValues.clear();
  extremeNodes.clear();
  extreme_value = 0;
  if(fig){
    delete fig;
    fig = 0;
  }

  if (cpContourObject){
    glDeleteLists(cpContourObject, 1);
    cpContourObject = 0;
  }
  
  if (mode == PLANE_AND_BOUND_MODE) {
    if (gridObject){
      glDeleteLists(gridObject, 1);
      gridObject = 0;
    }
    if (contourObject){
      glDeleteLists(contourObject, 1);
      contourObject = 0;
    }
    if (borderObject){
      glDeleteLists(borderObject, 1);
      borderObject = 0;
    }
    if (shadingObject){
      glDeleteLists(shadingObject, 1);
      shadingObject = 0;
    }
  }
  show_preview = show_contours = show_grid = show_shading = show_boundary_shading = show_border=show_extrema = false;
  makeObjects(); 
  updateGL();
}
//////////////////////////////////////////////////////////////////////////////
//  protected:
//    void mousePressEvent(QMouseEvent *event);
//
//  When the mouse is pressed, this function stores which pixel was pressed.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::mousePressEvent(QMouseEvent *event)
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

void GLViewer::mouseMoveEvent(QMouseEvent *event)
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
     double delta =(fabs(dx) > fabs(dy))?(-1*dx):(-1*dy);
     
     if (zbar > 0.05*size) {
       toz += delta*zbar/size;
       scale = 1.0;
     }
     else if (scale >= 1.0) {
       scale += 0.1*delta/(scale*size);
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


//
// Check for an intersection (HitPos) between a line(LP1,LP2) and a triangle face (TP1, TP2, TP3)
//
bool CheckLineTri( positions3d TP1, positions3d TP2, positions3d TP3, positions3d LP1, positions3d LP2, positions3d &HitPos)
  {
   positions3d Normal, IntersectPos;

   // Find Triangle Normal
   Normal = cross( TP2 - TP1, TP3 - TP1 );
   Normal = Normal/norm(Normal); // not really needed

   // Find distance from LP1 and LP2 to the plane defined by the triangle
   float Dist1 = dot((LP1-TP1), Normal );
   float Dist2 = dot((LP2-TP1), Normal );
   if ( (Dist1 * Dist2) >= 0.0f) return false;  // line doesn't cross the triangle.
   if ( Dist1 == Dist2) return false;// line and plane are parallel

   // Find point on the line that intersects with the plane
   IntersectPos = LP1 + (LP2-LP1) * ( -Dist1/(Dist2-Dist1) );

   // Find if the interesection point lies inside the triangle by testing it against all edges
   positions3d vTest;
   vTest = cross( Normal, TP2-TP1 );
   if ( dot( vTest, IntersectPos-TP1) < 0.0f ) return false;
   vTest = cross( Normal, TP3-TP2 );
   if ( dot( vTest, IntersectPos-TP2) < 0.0f ) return false;
   vTest = cross( Normal, TP1-TP3 );
   if ( dot( vTest, IntersectPos-TP1) < 0.0f ) return false;

   HitPos = IntersectPos;
   return true;
   }



//////////////////////////////////////////////////////////////////////////////
//  protected:
//    void mouseDoubleClickEvent(QMouseEvent *event);
//
//  This function handles changes to the camera when a mouse button is double
//  clicked.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::mouseDoubleClickEvent(QMouseEvent *event)
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
    

    
    vector<std::pair<double, unsigned int> > picked;
    for (unsigned int bid = 0; bid < mesh.size(); ++bid) {
      if (objVisible[bid]) {
        for(size_t i = 0; i <mesh[bid].size()/3; i++){ 
              
        
          positions3d t1 = meshNodes[mesh[bid][i*3 + 0]]; 
          positions3d t2 = meshNodes[mesh[bid][i*3 + 1]]; 
          positions3d  t3 = meshNodes[mesh[bid][i*3 + 2]];
          positions3d hitP;
          if(CheckLineTri(t1, t2, t3, P1, P2, hitP)){
            double winx, winy, winz;
             gluProject(hitP.x, hitP.y, hitP.z, modelviewMatr,
                 projectionMatrix, viewport, &winx, &winy, &winz);
            picked.push_back(std::make_pair(winz, bid));
            break;
                             
       
          }}}
    }
    if(picked.size() > 1) std::sort(picked.begin(), picked.end());
    if(picked.size() >0) {
     
      emit pickCurrent(int(picked[0].second)); 
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

void GLViewer::wheelEvent(QWheelEvent *event)
{
  roz += event->delta() / 8.0;
  updateGL();
}

//////////////////////////////////////////////////////////////////////////////
//  public slots:
//    void changeContours(int number);
//
//  This function is called by moving the contour slider. It calls the grid
//  object to change contour spacing and redraw the contours.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::changeContours(int number)
{
  if (fig) {
    fig->generate_contour_curves(number);
    if(contourObject)glDeleteLists(contourObject, 1);
    contourObject = makeContourObject();
    updateGL();
  }
}

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void updateView();
//
//  This function updates the camera position the fit all the visible
//  boundaries in the view at once.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::updateView()
{
  // Calculate global min/max pos of visible boundaries
  positions3d minpos, maxpos;
  bool first = true;
  for (size_t i = 0; i < objMinMax.size()/2; ++i) {
    if (objVisible[i]) {
      if (first) {
	minpos =  objMinMax[2*i];
        maxpos = objMinMax[2*i+1];
	first = false;
      } else {
	minpos.x = qMin(minpos.x, objMinMax[2*i + 0].x);
	minpos.y = qMin(minpos.y, objMinMax[2*i + 0].y);
	minpos.z = qMin(minpos.z, objMinMax[2*i + 0].z);
	
	maxpos.x = qMax(maxpos.x, objMinMax[2*i + 1].x);
	maxpos.y = qMax(maxpos.y, objMinMax[2*i + 1].y);
	maxpos.z = qMax(maxpos.z, objMinMax[2*i + 1].z);
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
  qDebug() << "size: " << size << " centerx: "<<centerx << "  centery: " << centery << " centerz: " << centerz ;
}




//////////////////////////////////////////////////////////////////////////////
//  public:
//    void previewCut(cutplane_info &Nfo);
//
//  Here is where the cut preview is made.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::previewCut(cutplane_info& Nfo)
{

  previewInfo = Nfo;
  mode = BOUND_SELECT_MODE;
  if (cpContourObject) {
    glDeleteLists(cpContourObject, 1);
    cpContourObject = 0;
  }
  cpContourObject = makeCPContour();
  show_preview = true;
  updateGL();

}


//////////////////////////////////////////////////////////////////////////////
//  public:
//    void cut(cutplane_info &Nfo);
//
//  This function takes the cutting plane info provided, initiates the cut,
//  and sets up the central widget to view the cut.
//
//  can not cut without previewCut first, because previewInfo not avaiable
//////////////////////////////////////////////////////////////////////////////

void GLViewer::cut()
{
  if(loadInfo.variable.isEmpty()||loadInfo.iteration.isEmpty()){
    QMessageBox::warning(this, tr("post_processing"),
                         tr("No scalar value, please run vogcheck first"));
   return;
  }
  if(cpContourObject==0)return;  
  
  info = previewInfo;
  if (fig){
    delete fig;
    fig = 0;
  }
  fig = new VolGrid;
  positions3d center = positions3d(centerx, centery, centerz);
   fig->cut(info, loadInfo, center);
  if(fig->triangle_list.size()==0) return;
  show_contours = true;
  show_grid = true;
  show_shading = true;
  show_border=true;
  show_boundary_shading = false;

  if (min_val == 0.0 && max_val == 0.0) {
    min_val = fig->min_val;
    max_val = fig->max_val;
  } else {
    min_val = qMin(min_val, fig->min_val);
    max_val = qMax(max_val, fig->max_val);
  }
  mode = PLANE_AND_BOUND_MODE;
  makeObjects();
  glViewport(0, 0, width(), height());
  updateGL();
  
}




void GLViewer::loadSca(){
  extremeValues.clear();
  extremeNodes.clear();
  
  // Get variable information
 hsize_t npnts;
 
 /* Turn off error handling */
 H5Eset_auto(H5E_DEFAULT,NULL,NULL);
 
 if(loadInfo.variable.isEmpty() || loadInfo.iteration.isEmpty()){
   QMessageBox::warning(this, tr("load scalar value"),
                        tr("no scalar value file exists, please run vogcheck first"));
   return;
 }
 
 QString posname = loadInfo.directory + "/output/grid_pos." + loadInfo.iteration + 
   '_' + loadInfo.casename ;
 
 hid_t file_id = H5Fopen(posname.toLocal8Bit(),
                         H5F_ACC_RDONLY, H5P_DEFAULT) ;

 if(file_id<0){
   QMessageBox::warning(this, tr("load scalar value"),
                        tr("unable to open ") + posname);
   return;
 }
 
 hid_t dataset_id = H5Dopen(file_id, "/pos/data", H5P_DEFAULT);
 hid_t dataspace_id = H5Dget_space(dataset_id);
 H5Sget_simple_extent_dims(dataspace_id, &npnts, NULL);
 hid_t pos_tid = H5Tcreate(H5T_COMPOUND, sizeof(positions3d));
 
 H5Tinsert(pos_tid, "x", 0, H5T_IEEE_F64LE);
 H5Tinsert(pos_tid, "y", sizeof(double), H5T_IEEE_F64LE);
 H5Tinsert(pos_tid, "z", 2*sizeof(double), H5T_IEEE_F64LE);

 vector<positions3d> nodePos;
 positions3d null3d;
 null3d.x = null3d.y = null3d.z = 0.0;
 nodePos.assign(npnts, null3d);
 H5Dread(dataset_id, pos_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodePos[0]);
 
 H5Tclose(pos_tid);
 H5Dclose(dataset_id);
 H5Fclose(file_id);



  
 QString filename = loadInfo.directory + "/output/" + loadInfo.variable + "_sca." + 
   loadInfo.iteration + "_" + loadInfo.casename;


 hid_t scalar_id = H5Fopen(filename.toLocal8Bit(), 
                           H5F_ACC_RDONLY, H5P_DEFAULT);
 
 if(scalar_id<0){
   QMessageBox::warning(this, tr("load scalar value"),
                        tr("unable to open ") + filename);
   return;
 }
 
 
 vector<float> nodeVal;
 nodeVal.assign(npnts, 0.0);
 QString datasetName = "/" + loadInfo.variable + "/data";
 dataset_id = H5Dopen(scalar_id, datasetName.toLocal8Bit(), H5P_DEFAULT);
 H5Dread(dataset_id, H5T_IEEE_F32LE, 
         H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodeVal[0]);
 H5Dclose(dataset_id);
 H5Fclose(scalar_id);
 

 {

   float minv = nodeVal[0];
   float maxv = nodeVal[0];
   
   for(size_t i = 0; i < npnts; i++){
     minv = min(minv, nodeVal[i]);
     maxv = max(maxv, nodeVal[i]);
   }
   min_val = minv;
   max_val = maxv;
 }
  
  if(max_val-min_val<1e-16) max_val = min_val = 0;
 
  
  if(meshMap.size()!= meshNodes.size()){
    QMessageBox::warning(this, tr("load scalar value"),
                         tr("mesh map doesn't exist"));
    return;
    
  }
  meshValues.resize(meshMap.size());
  for(unsigned int i = 0; i < meshValues.size(); i++){
    meshValues[i] = nodeVal[meshMap[i]];
  }

  //store 20% of nodes, grid_quality variable, top or bottom 20% nodes
  //other variables, top 10% and bottom 10% 
  {
    vector<pair<float, int> > sorted_nodeVal(npnts);
    for(unsigned int i = 0; i < npnts; i++){
      sorted_nodeVal[i] = pair<float, int>(nodeVal[i], i);
    }
    sort(sorted_nodeVal.begin(), sorted_nodeVal.end());

    mid_val = sorted_nodeVal[(int)(npnts/2)].first;
    
    int num_extreme_nodes = npnts /5;
    extremeNodes.resize(num_extreme_nodes);
    extremeValues.resize(num_extreme_nodes);
   
    if(loadInfo.variable=="cellVol"){
      for(int i = 0; i < num_extreme_nodes; i++){
        extremeValues[i] = nodeVal[sorted_nodeVal[i].second];
        extremeNodes[i] = nodePos[sorted_nodeVal[i].second];
      }
     
      
    }else if(loadInfo.variable=="cellFaceAngle"||
             loadInfo.variable=="cellShearTwist"||
             loadInfo.variable=="cellTwist"||
             loadInfo.variable=="nonconvex"||
             loadInfo.variable=="volumeRatio"){
      for(int i = 0; i < num_extreme_nodes; i++){
        extremeValues[i] = nodeVal[sorted_nodeVal[npnts-i-1].second];
        extremeNodes[i] = nodePos[sorted_nodeVal[npnts-i-1].second];
      }
     
    }else{
      int num_min_nodes = num_extreme_nodes/2;
      int num_max_nodes = num_extreme_nodes - num_min_nodes;
       for(int i = 0; i < num_min_nodes; i++){
         extremeValues[i] = nodeVal[sorted_nodeVal[i].second];
         extremeNodes[i] = nodePos[sorted_nodeVal[i].second];
       }
       for(int i = 1; i <= num_max_nodes; i++){
         extremeValues[num_extreme_nodes-i] = nodeVal[sorted_nodeVal[npnts-i].second];
         extremeNodes[num_extreme_nodes-i] = nodePos[sorted_nodeVal[npnts-i].second]; 
       } 
     
    }
    qDebug() << "min, max, medium of : "<<loadInfo.variable<<": "  <<min_val << "  " << max_val << "  " << mid_val;  
  }
 
  show_boundary_shading = true;
  mode = BOUND_SELECT_MODE;
  // extreme_value = mid_val;
  makeObjects();
  glViewport(0, 0, width(), height());
  updateGL();
}

void GLViewer::setExtrema(double value){
  extreme_value = value;
  show_extrema=true;
  // setShading(false);
  updateGL();
}

void GLViewer::setShading(bool b){
  if(show_boundary_shading !=b) toggleBoundaryShading();
}


double GLViewer::boundaryBoxSize(){
                                    
  return size;
}
  
//////////////////////////////////////////////////////////////////////////////
//  public slots:
//    void toggleContours();
//
//  This toggles the visiblility of the contours.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::toggleContours()
{
  show_contours = (show_contours)?false:true;
  updateGL();
}
void GLViewer::clearExtrema()
{
  show_extrema = false;
  updateGL();
}

void GLViewer::toggleBorder()
{
  show_border = (show_border)?false:true;
  updateGL();
}

void GLViewer::setShowShapes(bool b)
{
  show_shapes=b;
  updateGL();
}

void GLViewer::toggleBoundaryShading()
{
  show_boundary_shading = (show_boundary_shading)?false:true;
  makeObjects();
  updateGL();
}


// void GLViewer::toggleShowNodes()
// {
//   show_nodes = (show_nodes)?false:true;
//   updateGL();
// }


//////////////////////////////////////////////////////////////////////////////
//  public slots:
//    void toggleGrid();
//
//  This toggles the visibility of the grid.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::toggleGrid()
{
  show_grid = !show_grid;
  updateGL();
}

//////////////////////////////////////////////////////////////////////////////
//  public slots:
//    void toggleShading();
//
//  This toggles the visibility of the shading (coloring).
//////////////////////////////////////////////////////////////////////////////

void GLViewer::toggleShading()
{
  show_shading = (show_shading)?false:true;
 
  // makeObjects();
  updateGL();
}

//////////////////////////////////////////////////////////////////////////////
//  public slots:
//    setShadeType{1,2,3}()
//
//  These slots change the shading type and redraw the shading display list.
//  Since the signals and slots system in QT is very strict about type, this
//  the way I had to do it.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::setShadeType1()
{
  setShadeType(1);
}

void GLViewer::setShadeType2()
{
  setShadeType(2);
}

void GLViewer::setShadeType3()
{
  setShadeType(3);
}

void GLViewer::setShadeType(int type)
{
  this->show_shading = true;
  if (shadeType != type) {
    shadeType = type;
    makeObjects();
  }

 
  updateGL();
}


void GLViewer::showBoundaries(){
  bool checked = true;
  for(unsigned int i=0; i < objVisible.size(); i++){
    checked = objVisible[i] && checked;
  }
 

  for(unsigned int i=0; i < objVisible.size(); i++){
    objVisible[i] = !checked;
  }

  updateGL();
}
    

positions3d GLViewer::getTranslate(int b1, int b2){
  positions3d p1 = get_wireframe_center(meshNodes, mesh[b1]);
  positions3d p2 = get_wireframe_center(meshNodes, mesh[b2]);
 
  return p2-p1;
}

bool GLViewer::getRotate(int b1, int b2, double& angle, positions3d& axis, positions3d& center){
  positions3d p1 = get_average_normal(meshNodes, mesh[b1]);
  positions3d p2 = get_average_normal(meshNodes, mesh[b2]);
  
  positions3d p3 = get_wireframe_center(meshNodes, mesh[b1]);
  positions3d p4 = get_wireframe_center(meshNodes, mesh[b2]);

  p2  = -1.0*p2;
  bool result =  angleBetween(p1, p2, angle, axis);
  
  if(fabs(angle)> 1e-5){
    //round angle to 360/m, m is an integer
    double real_angle = angle*180/PI;
    int m= int(360.0/real_angle);
    real_angle = 360.0/m;
    angle=real_angle;
  }

  
  positions3d n= cross(p1, axis);
  
  
  double r = 0.5*norm(p4-p3);
  if(fabs(angle)> 1e-5) r= 0.5*norm(p4-p3)/sin(0.5*angle*PI/180);
  
  
  center = p3 + r*n; 
  
  
 return result;
  
  
}

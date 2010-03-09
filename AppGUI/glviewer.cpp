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


#include "glviewer.h"
#include "grid.h"
#include "hdf5.h"
#include "fvmadapt.h"
#define PI 3.14159265358979323846264338327950




#include <stdio.h>
#include <stdlib.h>

#ifndef CALLBACK
#define CALLBACK
#endif

GLuint startList;

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
  qobj = 0;
  //rgb = 0;
  centerx = centery = centerz=0;
  size = 0.1;
  currentObj = -1; //no select obj
  currentColor = default_color[0];
  tox= toy=toz = rox = roy = 0;
  show_contours = true;
  show_preview =  show_shading = show_grid = show_border = false;
  scale = 1.0;
  shadeType = 1;
  min_val = max_val = 0.0;
  isFit = false;
  mode=BOUND_SELECT_MODE;
  // rgb =0;
  adaptwindow = 0;
 show_shapes = true;
  show_nodes = false;
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
  if (gridObject)
    glDeleteLists(gridObject, 1);
  if (contourObject)
    glDeleteLists(contourObject, 1);
  if (borderObject)
    glDeleteLists(borderObject, 1);
  if (shadingObject)
    glDeleteLists(shadingObject, 1);

  if( cpContourObject) glDeleteLists(cpContourObject, 1);
  if (boundObjects.size() > 0) {
    for (size_t i = 0; i < boundObjects.size(); ++i)
      glDeleteLists(boundObjects[i], 1);
    boundObjects.clear();
  }
  if(fig){
    delete fig;
    fig =0;
  }

  gluDeleteQuadric(qobj);
}




void GLViewer::clearCurrent(){
  if(currentObj!=-1){
    makeBoundWireframeObject(currentObj, currentColor);
    currentObj=-1;
  
    updateGL();
  }
}





///////////////////////////////////////////
//  private:
//    void initializeGL();
//
//  Does some basic OpenGL configuration.
///////////////////////////////////////////

void GLViewer::initializeGL()
{
 //  glClearColor(1.0, 1.0, 1.0, 1.0);
//   glShadeModel(GL_SMOOTH);
//   glEnable(GL_DEPTH_TEST);
//   qobj = gluNewQuadric();
//   gluQuadricDrawStyle(qobj, GLU_FILL);
//   gluQuadricCallback(qobj, GLU_ERROR, 0);

 

  glClearColor(1.0, 1.0, 1.0, 1.0);
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


void GLViewer::drawSphere(const vector<double>& p){
 
  if(p.size() < 4)return;
  glPushMatrix(); 
  glTranslatef(p[0], p[1], p[2]);
  gluSphere(qobj, p[3], 20, 20);
  glPopMatrix();
  
}

void GLViewer::drawCylinder(const vector<double>& p){
 
  if(p.size()<5)return;
   
  glPushMatrix(); 
  glTranslatef(p[0], p[1], p[3]);
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
  
  glTranslatef(x0, y0, z1);
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
  glVertex3d(0,  0, 0);
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
  //glLineStipple(1, 0xAAAA);
  //  glEnable(GL_LINE_STIPPLE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex3d( -1*size, -1*size,z);
  glVertex3d( -1*size, size,z);
  glVertex3d( size, size,z);
  glVertex3d( size, -1*size,z);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(0,  0, 0);
  glVertex3d(0, 0, z - 0.5*size);
   glVertex3d( -1*size, 0, z);
  glVertex3d( size,0, z);
  glVertex3d( 0, -1*size, z);
  glVertex3d( 0,size, z);
  
  glEnd();
  //  glDisable(GL_LINE_STIPPLE);
  glPopMatrix();
  
}
void GLViewer::setAdaptWindow( QPointer<FVMAdapt> window){
  adaptwindow = window;
  // updateGL();
}
void GLViewer::adaptwindowClosed(){//signal not working
  adaptwindow = 0;
  updateGL();
}


void GLViewer::drawShapes(){
  show_shapes =  true;
  show_nodes = false;
  if(adaptwindow == 0)return;
  if(adaptwindow->tree == 0)return;
  QList<QTreeWidgetItem*> shapes = adaptwindow->tree->findItems(tr("object"), Qt::MatchRecursive);
  
  if(shapes.size() ==0)return;


  glLineWidth(5);
  for( int i = 0; i <shapes.size(); i++){

    if(shapes[i]->childCount() == 0)continue; //no child, do nothing
    glPushMatrix();
    double planesize = 1;
    QTreeWidgetItem* objItem = 0;
    for(int childId = 0; childId <shapes[i]->childCount(); childId++){
   
    QTreeWidgetItem* transItem = 0;
   
    if(shapes[i]->child(childId)->text(0)=="transform"){
      transItem =shapes[i]->child(childId);
    }else if(shapes[i]->child(childId)->text(0)=="shape")  objItem = shapes[i]->child(childId);
    
    if(transItem != 0){
      
      for( int j= 0; j < transItem->childCount(); j++){
        if(transItem->child(j)->text(0) =="translate"){
          double x0 =0 , y0 = 0, z0 = 0;
          for(int k = 0; k < transItem->child(j)->childCount(); k++){
            if(transItem->child(j)->child(k)->text(0)=="x0")x0 = transItem->child(j)->child(k)->text(1).toDouble();
            else  if(transItem->child(j)->child(k)->text(0)=="y0")y0 = transItem->child(j)->child(k)->text(1).toDouble();
            else  if(transItem->child(j)->child(k)->text(0)=="z0")z0 = transItem->child(j)->child(k)->text(1).toDouble();
            else{
              
              qDebug()<<tr("illegal child ") + transItem->child(j)->child(k)->text(0) + tr(" in 'translate'");
                                   
              return;
            }
            
          } 
             
        if(x0!=0 || y0 !=0 || z0 !=0)
          glTranslated(x0,
                       y0,
                       z0) ;
        }else if(transItem->child(j)->text(0) =="scale"){
          double x0 =1 , y0 = 1, z0 = 1;
          for(int k = 0; k < transItem->child(j)->childCount(); k++){
            if(transItem->child(j)->child(k)->text(0)=="x0")x0 = transItem->child(j)->child(k)->text(1).toDouble();
            else  if(transItem->child(j)->child(k)->text(0)=="y0")y0 = transItem->child(j)->child(k)->text(1).toDouble();
            else  if(transItem->child(j)->child(k)->text(0)=="z0")z0 = transItem->child(j)->child(k)->text(1).toDouble();
            else{
             
              qDebug() << tr("illegal child ") + transItem->child(j)->child(k)->text(0) + tr(" in 'scale'");
             
              return;
            }
            
          } 
          
        if(x0!=1 || y0 !=1 || z0 !=1)
          glScaled(x0,
                   y0,
                   z0) ;
         planesize = planesize*x0;
        }else if(transItem->child(j)->text(0) =="rotateX"){

          double theta = 0;
          for(int k = 0; k < transItem->child(j)->childCount(); k++){
            if(transItem->child(j)->child(k)->text(0)=="theta")theta = transItem->child(j)->child(k)->text(1).toDouble();
            else{
             
              qDebug()<< tr("illegal child ") + transItem->child(j)->child(k)->text(0) + tr(" in 'rotateX'");
             
              return;
            }
          }
             
          if(theta != 0)glRotatef(theta, 1, 0, 0);
        }else if(transItem->child(j)->text(0) =="rotateY"){

          double theta = 0;
          for(int k = 0; k < transItem->child(j)->childCount(); k++){
            if(transItem->child(j)->child(k)->text(0)=="theta")theta = transItem->child(j)->child(k)->text(1).toDouble();
            else{
           
              qDebug()<<tr("illegal child ") + transItem->child(j)->child(k)->text(0) + tr(" in 'rotateY'");
           
              return;
            }
          }
             
          if(theta != 0)glRotatef(theta, 0, 1, 0);
          }else if(transItem->child(j)->text(0) =="rotateZ"){

            double theta = 0;
            for(int k = 0; k < transItem->child(j)->childCount(); k++){
              if(transItem->child(j)->child(k)->text(0)=="theta")theta = transItem->child(j)->child(k)->text(1).toDouble();
              else{
               
                qDebug() << tr("illegal child ") + transItem->child(j)->child(k)->text(0) + tr(" in 'rotateZ'");
                return;
              }
            }
            
            if(theta != 0)glRotatef(theta, 0, 0, 1);
          }else{
             
            qDebug()<< tr("illegal child ") + transItem->child(j)->text(0) + tr(" in 'transform'");
            return;
                 
          }
      }
    }
    }
    if(objItem == 0 || objItem->childCount() == 0) {

     if(objItem) qDebug()<< objItem->text(0);
    
     qDebug()<<tr("no child 'shape' in 'object'");
     return;
    }
    
    

    objItem = objItem->child(0);
    vector<double> para;
    for(int j =0; j < objItem->childCount(); j++){
      para.push_back(objItem->child(j)->text(1).toDouble());
    }
    
    
    if(objItem->text(0) =="sphere") 
      drawSphere(para);
    else if(objItem->text(0) =="cone") 
      drawCone(para);
    else if(objItem->text(0) =="cylinder") 
      drawCylinder(para); 
    else if(objItem->text(0) =="box") 
      drawCube(para);
    else if(objItem->text(0) =="x+plane")
        drawPxPlane(para, planesize);
    else if(objItem->text(0) =="x-plane")
      drawNxPlane(para, planesize);
    else if(objItem->text(0) =="y+plane")
       drawPyPlane(para, planesize);
    else if(objItem->text(0) =="y-plane")
      drawNyPlane(para, planesize);
    else if(objItem->text(0) =="z+plane")
      drawPzPlane(para, planesize);
    else if(objItem->text(0) =="z-plane")
      drawNzPlane(para, planesize);
    
  

    glPopMatrix(); 
    
  }
   glLineWidth(1);
}
void GLViewer::markNodes(){
   show_shapes = false;
  show_nodes = true;
  if(adaptwindow==0) return;
  QDomDocument doc = adaptwindow->toDom();
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
  hid_t fi = H5Gopen(input_fid,"file_info") ;
   unsigned long numNodes = readAttributeLong(fi,"numNodes") ;
  
  H5Gclose(fi) ;
  
  hsize_t count = numNodes ;
  
#ifdef H5_INTERFACE_1_6_4
  hsize_t lstart = 0 ;
#else
  hssize_t lstart = 0 ;
#endif
  
    // Read in pos data from file i
  vector<positions3d> pos_dat(numNodes) ;
  hid_t node_g = H5Gopen(input_fid,"node_info") ;
  hid_t dataset = H5Dopen(node_g,"positions") ;
  hid_t dspace = H5Dget_space(dataset) ;

  

  hid_t pos_tid = H5Tcreate(H5T_COMPOUND, sizeof(positions3d));
    
  H5Tinsert(pos_tid, "x", 0, H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "y", sizeof(double), H5T_IEEE_F64LE);
  H5Tinsert(pos_tid, "z", 2*sizeof(double), H5T_IEEE_F64LE);

  hsize_t stride = 1 ;
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
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
  vector<bool> volumeTags(numNodes);

  QDomDocument doc = adaptwindow->toDom();
  if(doc.firstChildElement().isNull())return;
  
  QDomElement rootElement = doc.firstChildElement("region");
  vector<bool> vtags = process_region(rootElement, pos_dat);
  QString tagFileName =filename.section('.', 0, 0)+".tag";
  
  tagFileName = QFileDialog::getSaveFileName(this, tr("Save .tag File"),
                                             tagFileName,
                                             tr("tag Files (*.tag)"));
  
  



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

void GLViewer::drawNodes(){
    if(tags.size() != meshNodes.size()) return;
    if(tags.size()==0) return;
    glPointSize(5);
    glBegin(GL_POINTS);
    for(unsigned int i = 0; i < tags.size(); i++){
      if(tags[i]){
        glVertex3f(meshNodes[i].x, meshNodes[i].y, meshNodes[i].z);
      }
  }
     glEnd();
  glPointSize(1);
}
  
void GLViewer::paintGL()
{
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // qDebug() << "size in painGL: " << size << " centerx: "<<centerx << "  centery: " << centery << " centerz: " << centerz ;;
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
      glCallList(cpContourObject);
      glDisable(GL_LINE_SMOOTH);
    }
    //drawBoxes();
    //drawCoordGrid();

    break;

  case PLANE_AND_BOUND_MODE:
    
    for (size_t i = 0; i < boundObjects.size(); ++i)
      if (objVisible[i])glCallList(boundObjects[i]);
    if (show_shading)glCallList(shadingObject);
    if (show_grid)glCallList(gridObject);
    glEnable(GL_LINE_SMOOTH);
    if(show_border)glCallList(borderObject);
    glDisable(GL_LINE_SMOOTH);
    if (show_contours)glCallList(contourObject);
    break;
     
  case PLANE_ONLY_MODE:
    if (show_shading)
      glCallList(shadingObject);
    //    if (show_grid)
    // glCallList(gridObject);
    //glEnable(GL_LINE_SMOOTH);
    // glCallList(borderObject);
    //glDisable(GL_LINE_SMOOTH);
    //  if (show_contours)
    // glCallList(contourObject);
    
    break;
  default:
    break;
  }
 
  glColor3f(0.0, 0.0, 0.0);  
  if(show_nodes) drawNodes();
 
  if(show_shapes)drawShapes();
  glPopMatrix();
  glFlush();
  
  
}
void GLViewer::setCurrentObj(int i, QColor c){
  if( i>=0 && i<=(int)boundObjects.size()){
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




void GLViewer::setLoadInfo(const LoadInfo& ld_info){
  
   loadInfo.casename = ld_info.casename;
   loadInfo.directory = ld_info.directory;
   loadInfo.iteration = ld_info.iteration;
   loadInfo.variable = ld_info.variable;
}

bool GLViewer::load_boundary(QString fileName,  QStringList& boundary_names) {
  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
  loadInfo.casename = casename;
  loadInfo.directory = directory;
 
  

  QStringList format = fileName.split('.');
  if(format.size()!=2){
     QMessageBox::warning(this, tr("Application"),
                          fileName + tr(" has no postfix"));
     return false;
  }

  if(format[1] =="surface"){

    QFile file(fileName);
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
    int id;
    in >>id>>ntris >> nquads>>ngens ;
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
  file.close();
  gluDeleteTess(myTess);
  delete  pntIndex;
  }else if(format[1] == "surf"){
    //first read in meshNodes
    meshNodes.clear();
    QFile file(fileName);  
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Application"),
                           tr("Cannot read file %1:\n%2.")
                           .arg(fileName)
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
     

    meshNodes.resize(num_nodes);
    
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
      //  if(Read_Flag <3){
      //         QMessageBox::warning(this, tr("Application"),
      //                              tr("Error in reading file %1 nodes\n")
      //                              .arg(fileName));
      //         return;
      //       }
      meshNodes[i] = p;
    }



    //if .name file exists, read in and read mesh   
    QString name_file = format[0]+".names";
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
      boundary_names.clear();
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
    
      mesh.clear();
      mesh.resize(ids.size());
        
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
        //  if(Read_Flag <4){
        //         QMessageBox::warning(this, tr("Application"),
        //                              tr("Error in reading file %1 triangles\n")
        //                              .arg(fileName));
        //         return;
        //       } 
        mesh[id_map[id]].push_back(p1-1);
        mesh[id_map[id]].push_back(p2-1);
        mesh[id_map[id]].push_back(p3-1);
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
     //  if(Read_Flag <5){
//         QMessageBox::warning(this, tr("Application"),
//                              tr("Error in reading file %1 quads\n")
//                              .arg(fileName));
//         return;
//       } 
      mesh[id_map[id]].push_back(p1-1);
      mesh[id_map[id]].push_back(p2-1);
      mesh[id_map[id]].push_back(p3-1);

      mesh[id_map[id]].push_back(p1-1);
      mesh[id_map[id]].push_back(p3-1);
      mesh[id_map[id]].push_back(p4-1);
      
      }
    }else{

      mesh.clear();
      map<int, int> id_index_map;
        
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
        //  if(Read_Flag <4){
        //         QMessageBox::warning(this, tr("Application"),
        //                              tr("Error in reading file %1 triangles\n")
        //                              .arg(fileName));
        //         return;
        //       }
        if(id_index_map.find(id) == id_index_map.end()){
          mesh.push_back(vector<int>());
          id_index_map[id] = mesh.size() -1;
        }
        mesh[id_index_map[id]].push_back(p1-1);
        mesh[id_index_map[id]].push_back(p2-1);
        mesh[id_index_map[id]].push_back(p3-1);
        
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
        //  if(Read_Flag <5){
        //         QMessageBox::warning(this, tr("Application"),
        //                              tr("Error in reading file %1 quads\n")
        //                              .arg(fileName));
        //         return;
        //       }

        if(id_index_map.find(id) == id_index_map.end()){
          mesh.push_back(vector<int>());
          id_index_map[id] = mesh.size() -1;
        }
        mesh[id_index_map[id]].push_back(p1-1);
        mesh[id_index_map[id]].push_back(p2-1);
        mesh[id_index_map[id]].push_back(p3-1);
        
        mesh[id_index_map[id]].push_back(p1-1);
        mesh[id_index_map[id]].push_back(p3-1);
        mesh[id_index_map[id]].push_back(p4-1);
      
      }
      boundary_names.clear();
      for(map<int, int>::const_iterator p= id_index_map.begin(); p != id_index_map.end(); p++){
        boundary_names << QString("BC_%1").arg((*p).first);
       
      }
    }
    
    file.close();
  }
  
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
  mode = BOUND_SELECT_MODE;
  
  makeObjects();
  //  updateGL();
  tags.clear();
  if(adaptwindow)show_shapes = true;
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
// void GLViewer::uncut(){
//   showBoundaries(true);
//   mode = BOUND_SELECT_MODE;
//   show_preview = false;
 
//   updateGL();
// }
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
    //  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatr);//get the matrix in painGL so tranformation will take effect
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

void GLViewer::changeContours(int)
{
  if (fig) {
    //    fig->generate_contour_curves(number);
    glDeleteLists(contourObject, 1);
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
  if (cpContourObject) glDeleteLists(cpContourObject, 1);
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

  if(cpContourObject==0)return;  
  
  mode = PLANE_AND_BOUND_MODE;
  info = previewInfo;
  if (fig){
    delete fig;
    fig = 0;
  }
  fig = new grid;
 
  positions3d center = positions3d(centerx, centery, centerz);
 
  fig->cut(info, loadInfo, center);
  if(fig->triangle_list.size()==0) return;

 
  this->show_contours = true;
  this->show_grid = false;
  this->show_shading = false;


  if (min_val == 0.0 && max_val == 0.0) {
    min_val = fig->min_val;
    max_val = fig->max_val;
  } else {
    min_val = qMin(min_val, static_cast<float>(fig->min_val));
    max_val = qMax(max_val, static_cast<float>(fig->max_val));
  }

  makeObjects();

  glViewport(0, 0, width(), height());
  
  updateGL();

 
}

float GLViewer::boundaryBoxSize(){
                                    
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

void GLViewer::toggleBorder()
{
  show_border = (show_border)?false:true;
  updateGL();
}

void GLViewer::toggleShowShapes()
{
  show_shapes = (show_shapes)?false:true;
  updateGL();
}

void GLViewer::toggleShowNodes()
{
  show_nodes = (show_nodes)?false:true;
  updateGL();
}


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

//////////////////////////////////////////////////////////////////////////////
//  public:
//    void loadGrid();
//
//  This function opens a 2dgv grid from a file and views it.
//////////////////////////////////////////////////////////////////////////////

void GLViewer::loadGrid(const char* /*fileName*/)
{
  // View with plane only
  mode = PLANE_ONLY_MODE;

  // Reset grid object
  if (fig)
    delete fig;
  
  // Make new grid from file
  fig = new grid();
  //  fig->input(fileName);
    
  // Reset viewing parameters
  show_contours = true;
  show_grid = false;
  show_shading = false;

  min_val = fig->min_val;
  max_val = fig->max_val;
  
  // Make all display lists
  makeObjects();
  
  // Set up for viewing
  // glViewport(0, 0, width(), height());

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
    

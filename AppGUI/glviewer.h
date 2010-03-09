#ifndef GLVIEWER_H
#define GLVIEWER_H

#include <QGLWidget>
#include <vector>
#include <QDomElement>
#include <QPointer>
#include "grid.h"

class grid;
class cutplane_info;
class QStringList;
class FVMAdapt;
//class QStandardItemModel;

//BOUND_SELECT_MODE: before the cut plane is generated
//PLANE_AND_BOUND_MODE: after the cut plane is generated
//PLANE_ONLY_MODE:when only cut plane is read from a file
enum opMode {BOUND_SELECT_MODE,  PLANE_AND_BOUND_MODE, PLANE_ONLY_MODE};

std::vector<bool> process_region(QDomElement& anode, std::vector<positions3d> p);


//////////////////////////////////////////////////////////////////////////////
//  Global Function
//
//  Places node index into the passed integar vector.
//////////////////////////////////////////////////////////////////////////////

void cbVertex2(void *vertex_data, void *user_data);

void cbEdgeFlag(GLboolean);

// **** GLViewer implementation ****

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void makeObjects();
//
//  This function calls all the other functions to remake all necessary display
//  lists.
//////////////////////////////////////////////////////////////////////////////



class GLViewer : public QGLWidget
{
  Q_OBJECT

public:
  GLViewer(QWidget *parent = 0);
  ~GLViewer();
  void loadGrid(const char *fileName); //read in 2dgv file and view it
  float boundaryBoxSize(); // return the size of boudnary box, used in cutdialog
  bool load_boundary(QString filename,  QStringList& boundary_names); //read in .surf and .names file, or .surface file 
 public slots:
 //slots for cutplane display menu
  void toggleContours();
  void toggleGrid();
  void toggleShading();
  void toggleBorder();
  void toggleShowShapes();
  void toggleShowNodes();
  void setShadeType1();
  void setShadeType2();
  void setShadeType3();
  void changeContours(int number);
 
  void cut(); //read in volume grid and sca value, generate a cutplane
  // void uncut();
  void setCurrentObj(int i, QColor c); 
  void setVisibility(int , bool );
  void showBoundaries(); // slot for vis menu
  void clearCurrent();
  void previewCut(cutplane_info& Nfo);
  void reset();
  void fit();
  void setLoadInfo(const LoadInfo&);

  //for FVMadapt
  void drawShapes();
  void setAdaptWindow(QPointer<FVMAdapt>);
  void adaptwindowClosed();
  void markNodes();
  void markVolumeNodes(QString fileName);
  //  void refineGrids(QString );
signals:
  void pickCurrent(int);
protected:
  void initializeGL();
  void paintGL();
  void resizeGL(int height, int width);
   void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);

private:

  
  void updateView();
  void setShadeType(int type);
  void makeObjects();
  void drawBoxes();
  void drawCoordGrid();
  void triangulate(QString part, vector<int>* vTri);
  void makeBoundMesh();
  void makeBoundObjects();
  void makeBoundShadingObject(int bid);
  void drawBoundObject(int bid, QColor c);
  void makeBoundWireframeObject(int bid, QColor c);
  void makeBoundFillObject(int bid, QColor c);

  //  void draw();
  void drawSphere(const vector<double>& p);
  void drawCylinder(const vector<double>& p);
  void drawCone(const vector<double>& p);
  void drawCube(const vector<double>& p);
  void drawPxPlane(const vector<double>& p, double size);
  void drawNxPlane(const vector<double>& p, double size);
  void drawPyPlane(const vector<double>& p, double size);
  void drawNyPlane(const vector<double>& p, double size);
  void drawPzPlane(const vector<double>& p, double size);
  void drawNzPlane(const vector<double>& p, double size);
  void drawNodes();
 
  QPointer<FVMAdapt> adaptwindow ;
  vector<bool> tags;


  
  GLUquadricObj* qobj;//for FVMAdapt
  GLuint makeGridObject();
  GLuint gridObject;  // Holds grid display list

  GLuint makeContourObject();
  GLuint contourObject;  // Holds contour display list
  
  GLuint makeBorderObject();
  GLuint borderObject;  // Holds border display list

  GLuint makeShadingObject();
  GLuint shadingObject;  // Holds shading display list

  GLuint makeCPContour();
  GLuint cpContourObject;  // Holds cut preview display list
  vector<GLuint> boundObjects;  // Display lists of boundaries

  
  // Boundary mesh info
  vector<positions3d> meshNodes; 
  vector<float> meshValues;
  // map<int, int> meshMap;//from local index in meshNodes to global index
  vector<vector<int> > mesh;  // boundaries
 
 
 
 
  vector<bool> objVisible;  // Visibility of all boundaries
 
  vector<positions3d> objMinMax;  // Extrema of each boundary

  int currentObj;//current object
  QColor currentColor;// the color of current object

   double* rgb;  // Temporary used for display list creation
  double* shade(double value, double weight = 1.0);
  int shadeType;  // Stores which type of shading to use

  opMode mode;  // Which state the program is in
 
  grid *fig;  // Data structure that holds the cutting plane's topology

  positions lastPos;  // Last position visited by the mouse
  
  GLdouble size; //diagonal size of objects
  GLdouble centerx ,centery ,centerz, tox, toy, toz, rox, roy, roz;//for display
  GLdouble scale;
  
  bool isFit; // if this is in fit state
  
  
  
  int currentWidth, currentHeight;//viewport
  
  bool show_preview, show_contours, show_grid, show_shading, show_border, show_nodes, show_shapes;  // Visibility flags
  float min_val, max_val;  // Scalar value extrema over the whole grid
  
  cutplane_info info;  // The information for the current cutting plane
  cutplane_info previewInfo;  // Information for the cut being previewed
  LoadInfo loadInfo;

  GLdouble modelviewMatr[16];
  
};

#endif

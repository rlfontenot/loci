#ifndef GLVIEWER_H
#define GLVIEWER_H

#include <QGLWidget>
#include <vector>
#include <QDomElement>
#include <QPointer>
#include "grid.h"
#include "pboperation.h"

class VolGrid;
class cutplane_info;
class QStringList;
class QTreeWidgetItem;
//class QStandardItemModel;

//BOUND_SELECT_MODE: before the cut plane is generated
//PLANE_AND_BOUND_MODE: after the cut plane is generated
//PLANE_ONLY_MODE:when only cut plane is read from a file
enum opMode {BOUND_SELECT_MODE,  PLANE_AND_BOUND_MODE, PLANE_ONLY_MODE};

std::vector<bool> process_region(const QDomElement& anode, std::vector<positions3d> p);


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
  // void loadGrid(const char *fileName); //read in 2dgv file and view it
  double boundaryBoxSize(); // return the size of boudnary box, used in cutdialog
  bool load_boundary(QString filename,  QStringList& boundary_names); //read in .surface file
  bool load_boundary(QString filename,
                     QStringList& boundary_names,
                     QList<int>& bids); //for pb, boundary ids need to be read in
  double inline  get_min_val(){return min_val;};
  double inline  get_max_val(){return max_val;}; // Scalar value extrema over the whole grid                                  
                                                                                      
 public slots:
 //slots for cutplane display menu
  void toggleContours();
  void toggleGrid();
  void toggleShading();
  void toggleBoundaryShading();
  void setShading(bool); //show_exreme_nodes decide shading
  void toggleBorder();
  void toggleShowShapes();
  void clearExtrema();
  //  void toggleShowNodes();
  void setShadeType1();
  void setShadeType2();
  void setShadeType3();
  void changeContours(int number);
  void cut(); //read in volume grid and sca value, generate a cutplane
  void loadSca(); //read in sca value, shading the boundaries
  void setExtrema(double);
  void previewCut(cutplane_info& Nfo);
  void uncut();
  
  //for general display
  void setCurrentObj(int i, QColor c); 
  void setVisibility(int , bool );
  void showBoundaries(); // slot for vis menu
  void clearCurrent();
  void reset();
  void fit();
  void setLoadInfo(const LoadInfo&);

  //for FVMadapt
  void drawShapes();
  void cleanDoc();
  void updateDoc(const QTreeWidgetItem*);
  //void markNodes();
  //void markVolumeNodes(QString fileName);
  
  //fot pb

  positions3d getTranslate(int b1, int b2);
  bool getRotate(int b1, int b2, double& angle, positions3d& axis, positions3d& center);
  
  QSize sizeHint() const;
 
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
  void makeBoundShadingObjects();
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
  // void drawMarkedNodes();
  void drawExtremeNodes(double value);
 
  
  QDomDocument doc;//for FVMAdapt

  
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
  vector<double> meshValues;
  vector<int>  meshMap;//from local index in meshNodes to global index
  vector<vector<int> > mesh;  // boundaries
  vector<double> extremeValues;
  vector<positions3d> extremeNodes;//for scalar values display
   double extreme_value;
 
 
  vector<bool> objVisible;  // Visibility of all boundaries
 
  vector<positions3d> objMinMax;  // Extrema of each boundary

  int currentObj;//current object
  QColor currentColor;// the color of current object


  positions3d shade(double value, double weight = 1.0);
  int shadeType;  // Stores which type of shading to use

  opMode mode;  // Which state the program is in
 
  VolGrid *fig;  // Data structure that holds the cutting plane's topology
 

  

  positions lastPos;  // Last position visited by the mouse
  
  GLdouble size; //diagonal size of objects
  GLdouble centerx ,centery ,centerz, tox, toy, toz, rox, roy, roz;//for display
  GLdouble scale;
  
  bool isFit; // if this is in fit state
  
  
  
  int currentWidth, currentHeight;//viewport
  
  bool show_preview, show_contours, show_grid, show_shading,show_boundary_shading, show_border, show_extrema, show_shapes;  // Visibility flags
  double min_val, max_val;  // Scalar value extrema over the whole grid
  
  cutplane_info info;  // The information for the current cutting plane
  cutplane_info previewInfo;  // Information for the cut being previewed
  LoadInfo loadInfo;

  GLdouble modelviewMatr[16];
  
};

#endif

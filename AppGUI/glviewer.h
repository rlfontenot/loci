#ifndef GLVIEWER_H
#define GLVIEWER_H

#include <QGLWidget>
#include <vector>

#include "grid.h"

class grid;
class cutplane_info;
class QStringList;
//class QStandardItemModel;

//BOUND_SELECT_MODE: before the cut plane is generated
//PLANE_AND_BOUND_MODE: after the cut plane is generated
//PLANE_ONLY_MODE:when only cut plane is read from a file
enum opMode {BOUND_SELECT_MODE,  PLANE_AND_BOUND_MODE, PLANE_ONLY_MODE};




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
   void loadGrid(const char *fileName);
  float boundaryBoxSize();
  // void load_boundary(const LoadInfo& info,  QStringList& boundary_names);
  bool load_boundary(QString filename,  QStringList& boundary_names); 
public slots:
  void toggleContours();
  void toggleGrid();
  void toggleShading();
  void setShadeType1();
  void setShadeType2();
  void setShadeType3();
  void changeContours(int number);
  // void makeCut();
  void cut();
  void uncut();
  void setCurrentObj(int i, QColor c);
  void setVisibility(int , bool );
  void showBoundaries(bool);
  void clearCurrent();
  void previewCut(cutplane_info& Nfo);
  void reset();
  void fit();
  
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

  // double* rgb;  // Temporary used for display list creation
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
  
  bool show_preview, show_contours, show_grid, show_shading;  // Visibility flags
  float min_val, max_val;  // Scalar value extrema over the whole grid
  
  cutplane_info info;  // The information for the current cutting plane
  cutplane_info previewInfo;  // Information for the cut being previewed
  LoadInfo loadInfo;

  GLdouble modelviewMatr[16];
  
};

#endif

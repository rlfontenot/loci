#ifndef MGVIEWER_H
#define MGVIEWER_H

#include <QGLWidget>
#include <vector>
#include "defines.h"
#include "grid.h"

class grid;
class cutplane_info;
class QStringList;

bool CheckLineTri( positions3d TP1, positions3d TP2, positions3d TP3, positions3d LP1, positions3d LP2, positions3d &HitPos);

//////////////////////////////////////////////////////////////////////////////
//  Global Function
//
//  Places node index into the passed integar vector.
//////////////////////////////////////////////////////////////////////////////

void cbVertex2(void *vertex_data, void *user_data);

void cbEdgeFlag(GLboolean);

// **** MGViewer implementation ****

//////////////////////////////////////////////////////////////////////////////
//  private:
//    void makeObjects();
//
//  This function calls all the other functions to remake all necessary display
//  lists.
//////////////////////////////////////////////////////////////////////////////



class MGViewer : public QGLWidget
{
  Q_OBJECT

public:
  MGViewer(QWidget *parent = 0);
  ~MGViewer();
 
 
 
 public slots:
 
 void setCurrentColor(const IDColor&); 
  void setCurrentVisibility(const IDVisibility&);
  void showBoundaries(); // slot for vis menu
   void clearCurrent();
  void reset();
  void fit();
  bool load_boundary(QString filename); //read in .surf and .names file, or .surface file
  bool get_boundary(QString filename);
  void transGrid(const IDMatrix&);
  // void setLoadInfo(const LoadInfo&);
signals:
  void pickCurrent(const IDOnly&);
  void gridLoaded(const QStringList&);
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
  void makeObjects();
  void drawBoxes();

  void triangulate(QString part, vector<int>* vTri);

  void addBoundObjects();
  
  void drawBoundObject(int gid, int bid, QColor c);
  void makeBoundWireframeObject(int gid, int bid, QColor c);
  void makeBoundFillObject(int gid, int bid, QColor c);
  
  

    
  
  vector<vector<GLuint> > boundObjects;  // Display lists of boundaries
  
  // Boundary mesh info
  vector<vector<positions3d> > meshNodes; 
  //  vector<vector<float>  >meshValues;
  vector<vector<vector<int> > > mesh;  // boundaries
  vector<affineMapping2> gridXform ;
  
 
  vector<vector<bool> > objVisible;  // Visibility of all boundaries
 
  vector<vector<positions3d> > objMinMax;  // Extrema of each boundary
  
 
  positions lastPos;  // Last position visited by the mouse
  
  GLdouble size; //diagonal size of objects
  GLdouble centerx ,centery ,centerz, tox, toy, toz, rox, roy, roz;//for display
  GLdouble scale;
  
  bool isFit; // if this is in fit state
  
  int currentWidth, currentHeight;//viewport

  float min_val, max_val;  // Scalar value extrema over the whole grid

  GLdouble modelviewMatr[16];

  int currentGrid, currentObj;
  QColor currentColor;
  
};

#endif

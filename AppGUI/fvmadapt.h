#ifndef FVMADAPT_H
#define FVMADAPT_H

#include <QDialog>
#include <QDomDocument>
#include <QHideEvent>
#include "pages.h"
#include "vmergewindow.h"


class QLabel;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QCheckBox;
class QSlider;
class QListWidget;
class DomModel;
class QTreeWidgetItem;
class QStackedWidget;
class QTreeWidget;
class QSignalMapper;

  
enum BASIC_SHAPES{SPHERE, CONE, CYLINDER, BOX, LEFTPLANE};

QDomDocument tree2dom(const QTreeWidgetItem* root);

class Shape{
public:
  Shape(BASIC_SHAPES t, vector<double>& p):tp(t),para(p){};
  Shape(){};
  void reset();
public:
  BASIC_SHAPES tp;
  vector<double> para;
};

class ParaPage:public QGroupBox{
  Q_OBJECT
  public:
  ParaPage(Shape* s, QWidget* parent = 0);
                                         
public slots:
  void setValue(int);
  void showValue();
  void reset();
signals:
  void valueChanged();
public:
  
  vector<FloatSlider*> objs;
  Shape* shape;
  QSignalMapper* signalMapper;
};

struct IDMatrix2{
  int gridId;
  affineMapping matrix;
  IDMatrix2(int id, const affineMapping& m):gridId(id), matrix(m){}
};
  
class Transform: public QGroupBox
{
  Q_OBJECT
  public:
  Transform(QWidget *parent = 0);
  TranCoef value();
 
public slots:
  void setInfo();
  void clear();
  void setValue(positions3d* p);
signals:
  void tcChanged();
  
private:
  VectSlider* translate;
  VectSlider* rotateCenter;
  VectSlider* rotateAngle;
  VectSlider* scale;
};


class FVMAdapt: public QWidget
{
  Q_OBJECT

  public:
  FVMAdapt(QString fileName, QWidget *parent = 0);
  ~FVMAdapt();
  QTreeWidgetItem* getRoot();
 
  signals:
  void valueChanged(const QTreeWidgetItem*);
  void refineGrids();
                    
private slots:
  
  void helpClicked();
  void changePage(int);
  void addShape();
  void addTransform();
  void addRegion();
  void addOp();
  void removeNode();
  bool saveXml();
  void updateShape();
  void updateTransform();
  void showData(QTreeWidgetItem* item);  
  void done();
private:
  void createFlowBar();
  void createToolBar();
  void createTreeBar();
  QString filename;
  
  Transform* trans;
  QTreeWidget* tree;
  QTreeWidgetItem* root;
 
  QStackedWidget* paraPages;
  
  QGroupBox* flowbar;
  QGroupBox* toolbar;
  QGroupBox* treebar;
  vector<Shape*> defaultShapes;
};


#endif

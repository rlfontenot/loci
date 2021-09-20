
#ifndef FVMADAPT_H
#define FVMADAPT_H

#include <QDialog>
#include <QDomDocument>
#include <QHideEvent>
//#include "pages.h"
//#include "vmergewindow.h"
#include "glviewer.h"
#include "parapage.h"
#include "transform.h"
#include <QMainWindow>
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
class QTabWidget;
class QDockWidget; 
//enum BASIC_SHAPES{SPHERE, CONE, CYLINDER, BOX, LEFTPLANE};





class FVMAdapt : public QMainWindow
{
  Q_OBJECT

  public:
  FVMAdapt(QWidget *parent = 0);
  ~FVMAdapt();
  QTreeWidgetItem* getRoot();
 
signals:
  void valueChanged(const QTreeWidgetItem*);
private slots:
  
  void helpClicked();
  //void Viewer();
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
  void refineGrids();
  void loadGrid();
  void resizeTree();
  void validateTree();
private:
  void createFlowBar();
  void createToolBar();
  void buildTree();
  bool validateRegion(QTreeWidgetItem* item);
  bool validateObject(QTreeWidgetItem* item);
  // bool validateShape(QTreeWidgetItem* item);
  //void validateOp(QTreeWidgetItem* item);
  // void validateTransform(QTreeWidgetItem* item);
  //void validateTranslate(QTreeWidgetItem* item);
  
  
  QString filename;
  
  Transform* trans;
  QTreeWidget* tree;
  QTreeWidgetItem* root;
 
  QStackedWidget* paraPages;
  QTabWidget* tabWidget;

  QToolBar* toolbar;
  vector<Shape*> defaultShapes;

  GLViewer* viewer;
  QDockWidget* viewerDock;
 
 
};


#endif

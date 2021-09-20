
#ifndef VCUTWINDOW_H
#define VCUTWINDOW_H

#include <QDialog>
#include <QDomDocument>
#include <QHideEvent>
#include <QProcess>
#include "pages.h"
#include "glviewer.h"

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

class Transform;
class Shape;


class VCutWindow : public QMainWindow
{
  Q_OBJECT

  public:
  VCutWindow(QWidget *parent = 0);
  ~VCutWindow();
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
  void cutGrid();
  void afterCut(QString command, QProcess::ExitStatus status, QString directory);
  void loadGrid();
  void resizeTree();
  void validateTree();
  void addSphere();
private:
  void createFlowBar();
  void createToolBar();
  void buildTree();
  bool validateRegion(QTreeWidgetItem* item);
  bool validateObject(QTreeWidgetItem* item);
  
  QString inFilename;
  QString outFilename;
  Transform* trans;
  QTreeWidget* tree;
  QTreeWidgetItem* root;
 
  QStackedWidget* paraPages;
  QTabWidget* tabWidget;

  QToolBar* toolbar;
  vector<Shape*> defaultShapes;

  GLViewer* viewer;
  // QDockWidget* viewerDock;
  
  double x0, y0, z0, r;
 
};


#endif

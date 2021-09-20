#ifndef PBWINDOW_H
#define PBWINDOW_H

#include <QMainWindow>
#include "pages.h"
#include "glviewer.h"
#include <QMainWindow>
#include <QStandardItemModel>
#include <QTableView>
#include <QVBoxLayout>
#include <QList>
class QToolBar;
class QGroupBox;
class QCheckBox;
class QSlider;
class QSignalMapper;

  




class PbWindow : public QMainWindow
{
  Q_OBJECT

  public:
  PbWindow(QWidget *parent = 0);
 

 
signals:

private slots:
 
  void done();
  void loadGrid();
  void split();
  bool selectBoundary(const QStringList&);
 
  void showBoundary(QModelIndex, QModelIndex);
  void setCurrent(QModelIndex);
  void selectCurrent(int row);
  void updateBoundaryView(const QStringList&);
  void helpClicked();
  void setParameter();
private:
  void createFlowBar();
  void createToolBar();
  QString filename;
  QStringList bnames;
  QList<int> bids;
  QPointer<QStandardItemModel> modBoundaries;  // Boundary condition model
  QPointer<QTableView> boundaryView;  // Use for boundary Select
  QButtonGroup* boundaryGroup;
  QToolBar* toolbar;
  GLViewer* viewer;
  VectSpBox* translate;
  
  VectSpBox* rotateCenter;
  VectSpBox* rotateAxis;
  DoubleEdit* rotateAngle;
  QGroupBox* rotateGroup;
  QButtonGroup* buttonGroup;
  QStackedWidget* pagesWidget;
  
  int b1;
  int b2;
 
 
  QVBoxLayout* objLayout;
};


#endif

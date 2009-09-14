#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vector>
#include "grid.h"
#include "cutdialog.h"
#include <QDomDocument>

#include <QItemDelegate>
#include <QItemSelection>
#include <QTextEdit>
#include <QDomElement>
using std::vector;
class GLViewer;
class QDockWidget;
class QSlider;
class QTableView;
class QStandardItemModel;
 class BdCndWindow;
class QPushButton;
// Class that delegates for the boundary visibility model
class showDelegate: public QItemDelegate
{
   Q_OBJECT
public:
  QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem &option,
			const QModelIndex &index) const;
  showDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};

class colorDelegate: public QItemDelegate
{
  
  Q_OBJECT
public:
  bool  editorEvent(QEvent *event, QAbstractItemModel *model,
                    const QStyleOptionViewItem &option,
                    const QModelIndex &index);
  colorDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};





class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  
  MainWindow();
  QSize sizeHint() const;

  
 public slots:
 void newCase();
  void openCase();
  // void openCase(QString);
  void setGrid(QDomElement& theelem);
  void setBoundary(QDomElement& elem);
    bool selectBoundary();
  

  void cut();
  void resetSlider();
  void openSca();
  
  void openVog();
  bool saveVar();
  bool saveXml();
  bool saveImage();
  void aboutPreprocess();
  void aboutPostprocess();
  void showDisplayBar();
  void hideDisplayBar();
  void showVisBar();
  void hideVisBar();
  void showBoundary(QModelIndex, QModelIndex);
  void toggleViewer();
  
  void setCurrentObj(QModelIndex);
  void selectCurrent(int);
  void bdWindowClosed();
  void snapshot();

  void updateStatus(const QString&);
  void updateStatusTip(int);
  void clearAllStatus();
  void clearLastStatus();
  void changePage(int);
  signals:
  void setCurrent(QModelIndex);
  void stateChanged();
  void componentsChanged();
  
private:
   QString xmlpath;
  QString pngpath;


  void createActions();
  void createMenu();
  void createDisplayBar();
  void createVisBar();
  void createFlowBar();
  void createDockWindow();



  QToolBar* toolbar;//cutplane display
  QToolBar* visbar; //visualization
  QToolBar* tb; // file action
  QToolBar* flowbar;
  QMenu* viewMenu;
  QTextEdit* statusEdit;
  QGroupBox* statusWindow;

  GLViewer *viewer;  // Handle for central OpenGL widget
  QTableView* boundaryView;  // Use for boundary Select
  CutDialog* cutdialog;
  QStandardItemModel* modBoundaries;  // Boundary condition model
 
 
 
  QDockWidget *dock;
  QSlider* slider; // change number of contour

 
  QButtonGroup* flowbarButtons;
  QStackedWidget* central;
  
  
  
  BdCndWindow* bdWindow;


  
  QAction *cutAct;
  
  QAction *viewerAct; // toggle view of viewer
  QWidget* previousWidget;
  
  
  bool bdButtonDown; // for boundary select window reuse
  
  //the xml tree to store the  case information
  QDomDocument doc;
  bool isNewCase;
 
  
 };

 #endif


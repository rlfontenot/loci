#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <vector>
#include <QDomDocument>
#include <QItemDelegate>
#include <QItemSelection>
#include <QTextEdit>
#include <QDomElement>
#include <QPointer>
#include <QProcess>
using std::vector;
class GLViewer;
class QDockWidget;
class QSlider;
class QTableView;
class QStandardItemModel;
 class BdCndWindow;
class QPushButton;
class MGViewer;
class VMergeWindow;
class RefDialog;
class InitCndWindow;
#include "grid.h"
#include "cutdialog.h"
#include "vmergewindow.h"
#include "fvmadapt.h"

class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  
  MainWindow();
  QSize sizeHint() const;

  
public slots:

 
  
  void newCase();
  void openCase();
  void setGrid(QDomElement& theelem);
  void loadGrid(QString , QProcess::ExitStatus);
  void setBoundary(QDomElement& elem);
  bool selectBoundary();
  
  
  void cut();
  void check(const QString&);
  void vcheck();
  void resetSlider();
  void openSca();
  
  void openVog();
  bool saveVar();
  bool saveXml();
  bool saveImage();
  void aboutPreprocess();
  void aboutPostprocess();
 
  void showVisBar();
  void hideVisBar();
  void showBoundary(QModelIndex, QModelIndex);
  void toggleViewer();
  void selectBoundaryPressed();
  void updateConditionView();
  
  void setCurrentObj(QModelIndex);
  void selectCurrent(int);
  void bdWindowClosed();
  //visualization bar actions
  void snapshot();
  void reset();
  void fit();
  void clearCurrent();

  
  void updateStatus(const QString&);
  void updateStatusTip(int);
  void clearAllStatus();
  void clearLastStatus();
  void changePage(int);
  void toggleShowStatus();
  void vmClicked(); //vogmerge button clicked
  void vmClosed();
  void adaptClicked(); //FVMAdapt
 
  //void markVolumeNodes();
  void refineGrids();
  void showQuality(QString, QProcess::ExitStatus);
  
  signals:
  void setCurrent(QModelIndex);
  void stateChanged();
  void componentsChanged();
  void showStatus(const bool&);


 
private:
  QString xmlpath;
  QString pngpath;
  

  void createActions();
  void createMenu();
  
  void createVisBar();
  void createFlowBar();
  void createDockWindow();



 
  QPointer<QToolBar> visbar; //visualization
  QPointer<QToolBar> tb; // file action
  QPointer<QToolBar> flowbar;
  QPointer<QMenu> viewMenu;
  QPointer<QTextEdit> statusEdit;
  QPointer<QGroupBox> statusWindow;

 
  
  QPointer<QStandardItemModel> modBoundaries;  // Boundary condition model
  QPointer<QTableView> boundaryView;  // Use for boundary Select
  
  QPointer<GLViewer> viewer;  // Handle for central OpenGL widget
  QPointer<MGViewer> mgviewer;
  QPointer<VMergeWindow> vmwindow;

  QPointer<FVMAdapt> adaptwindow;
  
  QPointer<CutDialog> cutdialog;
 
 
 
 
  QPointer<QDockWidget> dock;
  QPointer<QDockWidget> bdock;

  

  
  QPointer<QButtonGroup> flowbarButtons;
  QPointer<QStackedWidget> central;
  
  
  
  QPointer<BdCndWindow> bdWindow;
  QPointer<InitCndWindow> initWindow;
  QPointer<RefDialog> refdialog;

  
 
  
  QAction *viewerAct; // toggle view of viewer
  QWidget* previousWidget;
  
  
  bool bdButtonDown; // for boundary select window reuse

  
  //the xml tree to store the  case information
  QDomDocument doc;
  bool isNewCase;
  bool displayStatus;

  bool waitForQualityFile;
  
 };

 #endif


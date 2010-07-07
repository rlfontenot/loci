#ifndef VARWINDOW_H
#define VARWINDOW_H
#include <QMainWindow>
#include <vector>
#include <QDomDocument>
#include <QItemDelegate>
#include <QItemSelection>
#include <QDomElement>
#include <QPointer>
#include <QProcess>
#include <QScrollArea>
using std::vector;
class GLViewer;
class QDockWidget;
class BdCndWindow;
class InitCndWindow;
class QMdiArea;
class QLabel;
class QButtonGroup;
class QStackedWidget;



class VarWindow : public QMainWindow
{
  Q_OBJECT
  public:
  VarWindow();
  QSize sizeHint() const;
public slots:
  void newCase();
  void openCase();
  void setGrid();
  void setBoundary();
  bool saveVar();
  bool saveXml();
  bool saveImage();
  void help(const QString&);
  void showVisBar();
  void hideVisBar();
 
 void check(const QString& fn);
  void showQuality(QString command, QProcess::ExitStatus status, QString directory);
  //visualization bar actions
  void snapshot();
  void updateStatus(const QString&);
  void updateStatusTip(int);
  
  void changePage(int);
  void toggleShowStatus();
  void updateModels();


  
signals:
  void setCurrent(QModelIndex);
  void stateChanged();
  void componentsChanged();
  void showStatus(const bool&);
  void caseChanged();

 
private:
  
 
  void createMenu();
  void createVisBar();
  void createFlowBar();
 
  /*!
    The tool bar for visualization. It include reset, fit, showBoundary, etc.  
  */
  QPointer<QToolBar> visbar; //visualization
   /*!
     The tool bar. It include file actions, dock widget toggle action, visualization bar toggle action,  etc.  
  */
  QPointer<QToolBar> tb; // file actions

  /*!
    The steps to generate a .var file.  
  */
  QPointer<QToolBar> flowbar;
  /*!
    The menu to turn on or off the dock widget or the tool bars.  
  */
  QPointer<QMenu> viewMenu;
   /*!
    The permanent message on status bar to show currrent status.  
  */
  QPointer<QLabel> statusLabel;

  /*!
    Graphics viewer
  */
  QPointer<GLViewer> viewer;  

  /*!
    The buttons in flow bar. This member data is needed when status tip is updated.  
  */
  QPointer<QButtonGroup> flowbarButtons;

   /*!
    Located in the central area of the window, contains all the user input interface.  
  */
  QPointer<QStackedWidget> central;
 /*!
   Located in the central area of the window, contains central.  
  */
  QPointer<QScrollArea> centralScrollArea;


  
   /*!
    The docked widget for \a viewer 
  */
  QPointer<QDockWidget> viewerDock;
  /*!
    Boundary condition window. This window might be recreated after a new grid in loaded.
  */ 
  QPointer<QWidget> bdWindow;
  /*!  
   the xml tree to store the  case information
  */  
  QDomDocument doc;
  /*!  
   If the status is displayed or not. This data is changed through the "show status" action on the tool bar, and signal showStatus() is emitted when its value changed.
   */  
  bool displayStatus;
  
};

#endif


#ifndef VMERGEWINDOW_H
#define VMERGEWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QModelIndex>
#include <QGroupBox>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QTableView>
#include <QList>
#include <QHBoxLayout>
#include <QMainWindow>
#include <utility>
#include "pages.h"
#include "defines.h"
#include "grid.h"
#include "progressdialog.h"
#include "mgviewer.h"
class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QCheckBox;
class QButtonGroup;
class QRadioButton;



class VMOption: public QGroupBox{
  Q_OBJECT
public:
  VMOption(int gridId, const QString &gridname, const QStringList &bcnames, QWidget *parent=0 );
 
  void setCurrentBid(int);
  signals:
  void tcChanged(const IDMatrix&);
  void  setCurrentColor(const IDColor&);
  void setCurrentVisibility(const IDVisibility&);
  void setCurrent(QModelIndex);
  
public slots:
  void setCurrentObj(QModelIndex);
  void setInfo();
  void showBoundary(QModelIndex, QModelIndex);
  void accept();
  void cancel();
  void previous();
  void next();
  void clear();
  void gRadioClicked();
public:
   QString currentText();
private:
  affineMapping2 currentM();
  void setValue(const TranCoef2&);
private:
  //translate
  VectSpBox* translateVec;
  //rotate angle
  VectSpBox* angleVec;
  //roate center
  VectSpBox* centerVec;
   //scale
  VectSpBox* scaleVec;
  //mirror option
  QButtonGroup* mButtonGroup;
  QRadioButton* radiog;
  
  QGroupBox* planeBox;
  //mirror plane origin
  VectSpBox* planeOriginVec;
  //mirror plane normal
  VectSpBox* planeNormalVec;

  DoubleEdit* tolEdit;
  QGroupBox* tolBox;
  
  //tag
  QLineEdit* tagEditor;

  QStandardItemModel* modBoundaries;  // Boundary condition model
  QTableView* boundaryView;  // Use for boundary Select
  

  QString gridName;
  int gridId;
  vector<TranCoef2> tc;
  int currentCoef;
  QString tag;
  QList<pair<QString, QString> > bdnames;
  QList<bool> needGlue;
 
};


class VMergeWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  VMergeWindow(QWidget* parent = 0);
 
 public slots:
  void vmClicked();
  void helpClicked();
  void loadGridClicked();
  void gridLoaded(const QStringList &);
  void changePage(int);
  void clearAll();
  void afterMerge(QString, QProcess::ExitStatus, QString);
  void selectCurrent(const IDOnly&);
  signals:
  void currentGridChanged(int);
  void tcChanged(const IDMatrix&);
  void  setCurrentColor(const IDColor&);
  void setCurrentVisibility(const IDVisibility&);
 
private:
  void clear();
  void createToolBar();
  void createFlowBar();
private:
  
  vector<vector<TranCoef2> > transcoef;
  //  QListWidget *typesWidget;
  QComboBox* typesWidget;
  QStackedWidget *pagesWidget;
  MGViewer* mgviewer;
  //  QHBoxLayout* visLayout;
 
};

#endif

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
public:
   QString currentText();
private:
  affineMapping2 currentM();
private:
  //translate
  DoubleEdit* xEditor1;
  DoubleEdit* yEditor1;
  DoubleEdit* zEditor1;
  //rotate angle
  DoubleEdit* xEditor2;
  DoubleEdit* yEditor2;
  DoubleEdit* zEditor2;
  //scale
  DoubleEdit* xEditor3;
  DoubleEdit* yEditor3;
  DoubleEdit* zEditor3;
  //roate center
  DoubleEdit* xEditor4;
  DoubleEdit* yEditor4;
  DoubleEdit* zEditor4;
  //tag
  QLineEdit* tagEditor;

  QStandardItemModel* modBoundaries;  // Boundary condition model
  QTableView* boundaryView;  // Use for boundary Select
  

  QString gridName;
  int gridId;
  vector<TranCoef> tc;
  int currentCoef;
  QString tag;
  QList<pair<QString, QString> > bdnames;
 
};


class VMergeWindow : public QWidget
{
  Q_OBJECT
  
public:
  VMergeWindow(QWidget* parent = 0);
 
 public slots:
 void vmClicked();
  void loadGridClicked();
  void gridLoaded(const QStringList &);
  void changePage(QListWidgetItem *, QListWidgetItem *);
  void clearAll();
  void afterMerge(QString, QProcess::ExitStatus);
  void selectCurrent(const IDOnly&);
  signals:
  void currentGridChanged(int);
  //  void loadGrid(QString);//add a grid
  // void getGrid(QString);//load in merged grid
  void tcChanged(const IDMatrix&);
  void  setCurrentColor(const IDColor&);
  void setCurrentVisibility(const IDVisibility&);
 
private:
  void clear();
  void createVisBar();
private:
  
  vector<vector<TranCoef> > transcoef;
  QListWidget *typesWidget;
  QStackedWidget *pagesWidget;
  MGViewer* mgviewer;
  QGroupBox* visbar;
 
};

#endif

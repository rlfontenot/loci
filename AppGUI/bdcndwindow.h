#ifndef BNCONDWINDOW_H
#define BNCONDWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QTextDocument>
#include <QModelIndex>
#include "pages.h"
#include "glviewer.h"
class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QTextEdit;
class QTableView;
class QComboBox;
class QStandardItemModel;
class BdCndWindow : public GeneralGroup
{
  Q_OBJECT
  
public:
  BdCndWindow(QDomElement& theelem, QPointer<GLViewer> theviewer,  QWidget* parent=0);
  
public slots:

  void setCurrent(QModelIndex);
  void selectCurrent(int);
 
private slots:
  bool selectBoundary();
  void selectBoundaryPressed();
  void setCurrent(int);
  
  void setCurrentObj(QModelIndex);
  void  updateConditionView();
  void updateBoundaryView();
  void showBoundary(QModelIndex, QModelIndex);
  void changePage(int);
  void checkStatus();
  void copyBdCnd(int previous, int current);
signals:
  void closed();
  
 
private:
  
 
  QStringList bdTypes;
  QComboBox *typesWidget;
  QStringList whatsThisList;
  QStringList toolTipList;
  QPointer<QStandardItemModel> modBoundaries;  // Boundary condition model
  QPointer<QTableView> boundaryView;  // Use for boundary Select
  QPointer<GLViewer> viewer;
  
  QStackedWidget *pagesWidget;
  QTextEdit* textEdit;
  QTextDocument text;
  QLabel* currentBdry;
};

#endif

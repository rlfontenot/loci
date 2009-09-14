#ifndef BNCONDWINDOW_H
#define BNCONDWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QTextDocument>
#include <QModelIndex>
#include "pages.h"
class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QTextEdit;
class QTableView;
class QComboBox;
class BdCndWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  BdCndWindow(QDomElement& theelem,  QDomElement& root, QStringList,  QTableView*, QWidget* parent=0);
  
public slots:

  void setCurrent(QModelIndex);
 
 private slots:
 void updateConditionView(const QString&);
  void updateConditionView();
  void changePage(int);
  void checkStatus();
  signals:
  void closed();

   
 
private:
  
  QStringList bdNames;
  QStringList bdTypes;
  QComboBox *typesWidget;
  QStringList whatsThisList;
  QStringList toolTipList;
 
  //QDomElement myroot;
   QDomElement cndNode;
  QStackedWidget *pagesWidget;
  QTextEdit* textEdit;
  QTextDocument text;
  QTableView* myTableView;  
};

#endif

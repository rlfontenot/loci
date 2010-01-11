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
  BdCndWindow(QDomElement& theelem,  QDomElement& root, QStringList,  QWidget* parent=0);
  
public slots:

void setCurrent(QModelIndex);
 
 
 private slots:

 void setCurrent(int);
  void changePage(int);
  void checkStatus();
  void copyBdCnd(int previous, int current);
  signals:
  void closed();
  void updateConditionView();
 
private:
  
  QStringList bdNames;
  QStringList bdTypes;
  QComboBox *typesWidget;
  QStringList whatsThisList;
  QStringList toolTipList;
 

  QDomElement cndNode;
  QStackedWidget *pagesWidget;
  QTextEdit* textEdit;
  QTextDocument text;
  QLabel* currentBdry;
};

#endif

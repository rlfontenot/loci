#ifndef INITCONDWINDOW_H
#define INITCONDWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QDomElement>
#include <QModelIndex>
#include "pages.h"
class QListWidget;
class QListWidgetItem;
class QStackedWidget;



class InitCndWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  InitCndWindow(QDomElement& elem, QDomElement& root, QWidget* parent = 0);
  
private slots:
void changePage(QListWidgetItem *current, QListWidgetItem *previous);
  void checkStatus();
  signals:

private:

  QListWidget *typesWidget;
  QString fileSelected;
  QStackedWidget *pagesWidget;
 
};

#endif

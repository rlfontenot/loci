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
class QTreeWidgetItem;


class InitCndWindow : public GeneralGroup
{
  Q_OBJECT
  
public:
  InitCndWindow(QDomElement& elem,  QWidget* parent = 0);
public slots:
  void setDirectory(QString s);
private slots:
  void changePage(int);
  void checkStatus();
signals:
  void valueChanged(const QTreeWidgetItem*);
  void directoryChanged(QString);
  void showShapes(bool);
private:

  QButtonGroup *typesWidget;
  QString fileSelected;
  QStackedWidget *pagesWidget;
 
};

#endif

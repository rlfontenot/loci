#ifndef SOLVERWINDOW_H
#define SOLVERWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QDomElement>
#include <QModelIndex>
#include "pages.h"
class QButtonGroup;
class QStackedWidget;



class SolverWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  SolverWindow(QDomElement& theelem, QDomElement& myroot, QWidget* parent = 0);
  
private slots:
void changePage(int id);
  void checkStatus();
  void updateState();// update buttons using pages's signal
  void updateShowStatus(const bool&);
  signals:
  
private:
  
  QButtonGroup *typesWidget;
   QStackedWidget *pagesWidget;
 
  
};

#endif

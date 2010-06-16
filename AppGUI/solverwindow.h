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
class QTabWidget;


class SolverWindow : public GeneralGroup
{
  Q_OBJECT
  
public:
  SolverWindow(QDomElement& theelem,QWidget* parent = 0);
  
private slots:
  void checkStatus();
  void changeState();// update buttons using pages's signal
  void updateShowStatus(const bool&);
  signals:
  
private:
  QTabWidget *pagesWidget;
 };

#endif

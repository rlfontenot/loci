#ifndef PHYSICSWINDOW_H
#define PHYSICSWINDOW_H

#include <QWidget>
#include <QDomDocument>
#include <QDomElement>
#include <QObject>
#include "pages.h"
#include "mdlwindow.h"

class QButtonGroup;
class QStackedWidget;
class QRadioButton;
class QSpinBox;

class PhysicsWindow : public GeneralGroup
{
  Q_OBJECT
  
  public:
  PhysicsWindow(QDomElement& theelem, QWidget* parent = 0);
signals:
    
private slots:
  void updateState(QString);
  void checkStatus();
private:

};

#endif

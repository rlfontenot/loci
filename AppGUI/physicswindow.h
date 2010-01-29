#ifndef PHYSICSWINDOW_H
#define PHYSICSWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QDomElement>
#include <QObject>
#include "pages.h"
#include "mdlwindow.h"

class QButtonGroup;
class QStackedWidget;
class QRadioButton;
class QSpinBox;













class PhysicsWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  PhysicsWindow(QDomElement& theelem, QDomElement& myroot, QWidget* parent = 0);
  // QStringList currentState();
  signals:
  //  void componentsChanged();
  void parentStateChanged(QString);
  void stateUpdated(const QStringList&);
private slots:
void updateState(QString);
  void checkStatus();
private:
  QStringList state; 
  QList<StackGroup2*> stacks;
 };

#endif

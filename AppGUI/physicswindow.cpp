#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "physicswindow.h"
#include "pages.h"
#include "getfile.h"

/*!
  \class PhysicsWindow
  
  \brief PhysicsWindow  allows the user to define the  physics state of the system and related variables. The user is be able to select a state, select a module or add a new model. When a model is selected, the components of the system is also changed.
  
  \mainclass
 
 
*/          
   

PhysicsWindow::PhysicsWindow(QDomElement& theelem, QWidget* parent)
  :GeneralGroup(theelem, parent)
{
  QStringList state;
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+ tr(" has no child")
                         );
    return;
  }
  int count=0;   
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) {   
    if(elem.attribute("element")=="stateStack"){
      StateStackGroup* stackGroup = new StateStackGroup(elem, this);
      state << stackGroup->myState();
      connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(updateState(QString)));
      connect(this, SIGNAL(stateChanged()), stackGroup, SLOT(changeState()));
      connect(this, SIGNAL(showStatus(const bool &)), stackGroup, SLOT(updateShowStatus(const bool &)));
      connect(stackGroup, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(stackGroup, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
      mainLayout->addWidget(stackGroup);
    }else{
      QMessageBox::warning(window(), myelem.tagName(),
                           tr("don't know how to handle it yet ")+ myelem.attribute("element") 
                           );
    }
  }
  
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  checkStatus();
  
  QString tmp = state.join(";");
  myelem.setAttribute("state", tmp);
  emit stateChanged();
}


/*!
  When child's state changed, update my state,
  and emit stateChanged()
 */ 
void PhysicsWindow::updateState(QString stat){
  //get original state
  QStringList state = myelem.attribute("state").split(';');
  //child state split into variableName and value
  QStringList list1=stat.split("=");
  bool Found = false;
  //find the variableName in the original state, update the new value
  //and emit stateChanged()
  for(int i = 0; i < state.size(); i++){
    if(state[i].contains(list1[0])){
      state[i] = stat;
      QString tmp = state.join(";");
      myelem.setAttribute("state", tmp);
      emit stateChanged();
      Found=true;
      break;
    }
   }
 
   if(!Found){
     QMessageBox::warning(window(), myelem.tagName(),
                          tr("illegal state ")+ stat 
                          );
     return;
   }
}





/*!
  Use child element's \a currentText and \a status update my \a currentText and \a status,
  emit updateStatusTip()
  
 */ 

 void PhysicsWindow::checkStatus(){
  int count = 0;
  int count_done = 0;
  QString text;
  for( QDomElement elt = myelem.firstChildElement(); !elt.isNull();
       elt = elt.nextSiblingElement()){
    if(elt.attribute("status")=="done")count_done++;
    count++;
    text += elt.attribute("currentText");
  }
   
  myelem.setAttribute("currentText", text);
  if(count==count_done) myelem.setAttribute("status", "done");
  else myelem.setAttribute("status", QString("%1  out of ").arg(count_done)+QString(" %1 finished").arg(count));
  emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
 }




        


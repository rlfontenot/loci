#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "physicswindow.h"
#include "pages.h"
#include "getfile.h"
#include "mdlwindow.h"













PhysicsWindow::PhysicsWindow(QDomElement& theelem, QDomElement& theroot, QWidget* parent)
  :GeneralWindow(theelem, theroot, parent)
{

  
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
      StackGroup2* stackGroup = new StackGroup2(elem, myroot, state, this);
      stacks << stackGroup;
      state << stackGroup->myState();
      connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(updateState(QString)));
      connect(this, SIGNAL(parentStateChanged(QString)), stackGroup, SLOT(parentStateChanged(QString)));
      connect(this, SIGNAL(stateUpdated(const QStringList&)), stackGroup, SLOT(setParentState(const QStringList&)));
       connect(this, SIGNAL(showStatus(const bool &)), stackGroup, SLOT(updateShowStatus(const bool &)));
      connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(checkStatus()));
      connect(stackGroup, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(stackGroup, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
      mainLayout->addWidget(stackGroup);
    }else{
      QMessageBox::warning(window(), myelem.tagName(),
                           tr("don't know how to handle it yet ")+ myelem.attribute("element") 
                           );
    }
  }
    emit stateUpdated(state);
    setLayout(mainLayout);
    setWindowTitle(myelem.attribute("title"));
    checkStatus();
}

void PhysicsWindow::updateState(QString stat){
  QStringList list1=stat.split("=");
  bool Found = false;
  for(int i = 0; i < state.size(); i++){
    if(state[i].contains(list1[0])){
      if(state[i]!=stat)emit parentStateChanged(stat);
      state[i] = stat;
      
      Found = true;
    }
  }
  if(Found){
    QString tmp = state.join(";");
    if(tmp != myelem.attribute("state")){
      myelem.setAttribute("state", tmp);
      if(!tmp.isEmpty())emit stateChanged();
    }
  }
  if(!Found){
    QMessageBox::warning(window(), myelem.tagName(),
                          tr("illegal state ")+ stat 
                         );
    return;
  }
  
  
}







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


//   QString tmp = state.join(";");
//   myelem.setAttribute("state", tmp);
//    if(!tmp.isEmpty())emit stateChanged();
  emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
 }



void StackGroup2::add(){
  
   QDomElement elm = myelem.firstChildElement();
   if(elm.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                              tr("stack elememt has no child")
                          );
     return;
   }
   for(int i = 0; i < current; i++)elm = elm.nextSiblingElement();
    if(elm.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                              tr("stack elememt has not enough children")
                          );
     return;
    }
  
    //   QDomElement elt_page = myroot.firstChildElement("models");
  
    // if(elt_page.isNull()){
//       QMessageBox::warning(window(), tr(".xml"),
//                            tr("can not find element 'models' ")
//                            );
//       return;
//     }
   
//   QDomElement elem= elt_page.firstChildElement(elm.attribute("define"));
  
//   if(elem.isNull()){
//     QMessageBox::warning(window(), tr(".xml"),
//                          tr("can not find the model ")+elm.attribute("define")
//                          );
//     return;
//   }

     
 //  if(elem.attribute("element")=="all"){
    
//     //   AllVBWindow* window= new AllVBWindow(elem, myroot);
    
//     ChemistryMdl* window = new ChemistryMdl;
//     window->show();
         
//   }else if(elem.attribute("element")=="CpWindow"){
//     SpeciesGroup* window=new SpeciesGroup();
//     window->show();
    
//   }

    if(elm.attribute("define")=="specified_ideal_gas_model"){

      QDomElement elt_page = myroot.firstChildElement("models");
      
      if(elt_page.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find element 'models' ")
                             );
        return;
      }
      
      QDomElement elem= elt_page.firstChildElement(elm.attribute("define"));
  
      if(elem.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find the model ")+elm.attribute("define")
                         );
        return;
      }
      
      
      if(elem.attribute("element")=="all"){
        
        AllVBWindow* window= new AllVBWindow(elem, myroot);
        window->show();
      }
      
    }else if(elm.attribute("define")=="gamma_model"){

      QDomElement elt_page = myroot.firstChildElement("models");
      
      if(elt_page.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find element 'models' ")
                             );
        return;
      }
      
      QDomElement elem= elt_page.firstChildElement(elm.attribute("define"));
  
      if(elem.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find the model ")+elm.attribute("define")
                             );
        return;
      }
      
      
      if(elem.attribute("element")=="all"){
        
        AllVBWindow* window= new AllVBWindow(elem, myroot);
        window->show();
      }
      

      
    }else if(elm.attribute("define")=="cp_model"){
      CpWindow* window = new CpWindow;
      window->show();
    }else if(elm.attribute("define")=="specified_mixture_model"){
       ChemistryMdl* window = new ChemistryMdl(true);
       window->show();
    }else if(elm.attribute("define")=="specified_chemistry_model"){
      ChemistryMdl* window = new ChemistryMdl;
      window->show();
    }
    
}
        




        


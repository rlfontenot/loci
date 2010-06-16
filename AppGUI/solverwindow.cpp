#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <QTabWidget>
#include "solverwindow.h"
#include "pages.h"
#include "getfile.h"

SolverWindow::SolverWindow(QDomElement& theelem, QWidget* parent):
  GeneralGroup(theelem, parent)
{
  pagesWidget = new QTabWidget;
  pagesWidget->setTabPosition(QTabWidget::West);
  int current=0;
  if(myelem.hasAttribute("current"))current=myelem.attribute("current").toInt();
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                        myelem.tagName()+ tr(" has no child")
                         );
    return;
  }
  int count=0;   
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) {   
                
    //this attribute is used for changeState of panel/page
    //with "index", the panel/page does hide/show
    elem.setAttribute("index", count);
    QWidget* newPage = 0;
    if(elem.attribute("element")=="panel"){
      newPage = new VarPanel(elem);
      qobject_cast<VarPanel*>(newPage)->setTitle("");
      qobject_cast<VarPanel*>(newPage)->setFlat(true);
      
      connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newPage, SLOT(updateShowStatus(const bool&)));
      //      connect(newPage, SIGNAL(stateChanged()), this, SLOT(updateState()));
    }else if (elem.attribute("element")=="page"){
      newPage = new VarPage(elem);
      connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
      //connect(newPage, SIGNAL(stateChanged()), this, SLOT(updateState()));
      connect(this, SIGNAL(showStatus(const bool&)), newPage, SLOT(updateShowStatus(const bool&)));
    }else{
      QMessageBox::warning(window(), ".xml",
                           elem.tagName()+ tr(": don't know how to handle it: ")
                           +elem.attribute("element")
                           );
    }
    pagesWidget->addTab(newPage, elem.tagName());
    pagesWidget->setToolTip(elem.attribute("toolTip"));
    pagesWidget->setWhatsThis(elem.attribute("whatsThis"));
  }//for elem
  
  pagesWidget->setCurrentIndex(current);
   QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(pagesWidget);
  mainLayout->addStretch(1);
  mainLayout->addSpacing(12);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  changeState();//make sure the state is checked
  checkStatus();
}

void SolverWindow::updateShowStatus(const bool& show){

  if(show){
    QDomElement elem = myelem.firstChildElement();
    
    if(elem.isNull()){
      QMessageBox::warning(window(), ".xml",
                        myelem.tagName()+ tr(" has no child")
                           );
      return;
    }
    
    int count=0;   
    for (; !elem.isNull(); elem = elem.nextSiblingElement()) {   
      if(elem.attribute("status")!="done")break;
      count++;
    }
    if(count<pagesWidget->count())pagesWidget->setCurrentIndex(count);
  }
  GeneralGroup::updateShowStatus(show);
}

void SolverWindow::changeState(){
  //first use this function to update the attribute "conditionSatisfied"
  GeneralGroup::changeState();
  //disable/able the tabs
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+ tr(" has no child")
                         );
    return;
  }
  int count=0;
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++){
    if(elem.hasAttribute("condition")){
      if(elem.attribute("conditionSatisfied")=="true"){
        pagesWidget->setTabEnabled(count, true);
      }else{
        pagesWidget->setTabEnabled(count, false);
      }
    }
  }
 
}
    
void SolverWindow::checkStatus(){
  
  QDomElement myroot = myelem.ownerDocument().documentElement();
   int count = 0;
   int count_done = 0;
   QString text;
   if(myelem.hasAttribute("comments")) text += myelem.attribute("comments");
   QString infix = "\n";
   for( QDomElement elt = myelem.firstChildElement(); !elt.isNull();
        elt = elt.nextSiblingElement()){
     if(elt.attribute("status")=="done")count_done++;
     count++;
     if(elt.hasAttribute("outputCondition")){
       if(conditionIsSatisfied(myroot, elt.attribute("outputCondition"))) text += elt.attribute("output")+ infix;
     }
     
     if((elt.hasAttribute("condition")
         && elt.attribute("conditionSatisfied")=="true")
        ||(!elt.hasAttribute("condition")))text +=elt.attribute("currentText") + infix; 
   }
   myelem.setAttribute("currentText", text);
   
   if(count==count_done) myelem.setAttribute("status", "done");
   else myelem.setAttribute("status", QString("%1  out of ").arg(count_done)+QString(" %1 finished").arg(count));
   emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
 }











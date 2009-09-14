#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "solverwindow.h"
#include "pages.h"
#include "getfile.h"

SolverWindow::SolverWindow(QDomElement& theelem, QDomElement& theroot, QWidget* parent):
  GeneralWindow(theelem, theroot, parent)
{
 
  typesWidget = new QButtonGroup(this);
  pagesWidget = new QStackedWidget;

  QVBoxLayout *mainLayout = new QVBoxLayout;

  
  QGroupBox* buttonGroupBox=new QGroupBox(myelem.attribute("label"));
  QVBoxLayout* buttonGroupLayout = new QVBoxLayout;
  buttonGroupBox->setLayout(buttonGroupLayout);
  
  
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
    
    QPushButton *bdCndiButton = new QPushButton(elem.tagName());
    typesWidget->addButton( bdCndiButton, count);
    buttonGroupLayout->addWidget(bdCndiButton);
    // bdCndiButton->setTextAlignment(Qt::AlignHCenter);
    bdCndiButton->setCheckable(true);
    
    //this attribute is used for changeState of panel/page
    //with "index", the panel/page does hide/show
    elem.setAttribute("index", count);
    
    if(count==current)bdCndiButton->setChecked(true);
    
    
    bdCndiButton->setToolTip(elem.attribute("toolTip"));
    bdCndiButton->setWhatsThis(elem.attribute("whatsThis"));
    QWidget* newPage = 0;
    if(elem.attribute("element")=="panel"){
      newPage = new OptionPage(elem, myroot);
      connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
      connect(newPage, SIGNAL(stateChanged()), this, SLOT(updateState()));
    }else if (elem.attribute("element")=="page"){
      newPage = new Page(elem, myroot);
      connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
      connect(newPage, SIGNAL(stateChanged()), this, SLOT(updateState()));
      
    }else{
      QMessageBox::warning(window(), ".xml",
                           elem.tagName()+ tr(": don't know how to handle it: ")
                           +elem.attribute("element")
                           );
    }
    pagesWidget->addWidget(newPage);
   
    
  }//for elem
  
  pagesWidget->setCurrentIndex(current);
  
  connect(typesWidget,
          SIGNAL(buttonClicked(int)),
          this, SLOT(changePage(int)));

  QHBoxLayout *horizontalLayout = new QHBoxLayout;
  horizontalLayout->addWidget(buttonGroupBox);
  horizontalLayout->addWidget(pagesWidget);
  

  mainLayout->addLayout(horizontalLayout);
  mainLayout->addStretch(1);
  mainLayout->addSpacing(12);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  updateState();
  checkStatus();
}

void SolverWindow::changePage(int id)
{
  pagesWidget->setCurrentIndex(id);
  checkStatus();
}

void SolverWindow::updateState(){
  
  
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
        typesWidget->button(count)->show();
      }else{
        typesWidget->button(count)->hide();
      }
    }
  }
}
    
 void SolverWindow::checkStatus(){
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











#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "initcndwindow.h"
#include "pages.h"
#include "getfile.h"

InitCndWindow::InitCndWindow(QDomElement& theelem, QDomElement& theroot, QWidget* parent): GeneralWindow(theelem, theroot, parent)
{
  typesWidget = new QListWidget;
  pagesWidget = new QStackedWidget;
  
  QHBoxLayout* aLayout = new QHBoxLayout;
  QLabel* aLabel = new QLabel(myelem.attribute("label"));
  aLayout->addWidget(aLabel);
  aLayout->addWidget(typesWidget);
    
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+tr(" has no child")
                         );
    return;
  }
  
  for (; !elem.isNull(); elem = elem.nextSiblingElement()) {   
      
    QListWidgetItem *bdCndiButton = new QListWidgetItem(typesWidget);
    bdCndiButton->setText(elem.tagName());
    bdCndiButton->setTextAlignment(Qt::AlignHCenter);
    bdCndiButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    bdCndiButton->setToolTip(elem.attribute("whatsThis"));
    bdCndiButton->setStatusTip(elem.attribute("whatsThis"));
    bdCndiButton->setWhatsThis(elem.attribute("whatsThis"));
    if(elem.hasAttribute("action") && elem.attribute("action")=="find file"){
      FindFileWindow* getFileWindow = new FindFileWindow(elem, fileSelected);
      pagesWidget->addWidget(getFileWindow);
    }else if(elem.hasAttribute("element")&&elem.attribute("element")=="panel"){
      OptionPage* bdCndPage = new OptionPage(elem,myroot);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      pagesWidget->addWidget(bdCndPage);
    }else if(elem.hasAttribute("element")&&elem.attribute("element")=="page"){
      Page* bdCndPage = new Page(elem, myroot);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      pagesWidget->addWidget(bdCndPage);
    }else{
      QMessageBox::warning(window(), elem.tagName(),
                           tr(" don't know how to handle it yet: ")
                           +elem.attribute("element")
                           );
      
    }
                    
  }
  if(myelem.hasAttribute("current")){
    typesWidget->setCurrentRow(myelem.attribute("current").toInt());
    pagesWidget->setCurrentIndex(myelem.attribute("current").toInt());
  }else{
    typesWidget->setCurrentRow(0);
    pagesWidget->setCurrentIndex(0);
    myelem.setAttribute("current", 0);
  }
  connect(typesWidget,
          SIGNAL(currentItemChanged(QListWidgetItem *, QListWidgetItem *)),
          this, SLOT(changePage(QListWidgetItem *, QListWidgetItem*)));
  
  
  
  QVBoxLayout *horizontalLayout = new QVBoxLayout;
  
  horizontalLayout->addLayout(aLayout);
  horizontalLayout->addWidget(pagesWidget);
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(horizontalLayout);
  mainLayout->addStretch(1);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  checkStatus();
}

void InitCndWindow::changePage(QListWidgetItem *current, QListWidgetItem* previous)
{
  if (!current) current = previous;
  myelem.setAttribute("current",typesWidget->row(current));
  pagesWidget->setCurrentIndex(typesWidget->row(current));
  checkStatus();
 }

 void InitCndWindow::checkStatus(){
   int current= myelem.attribute("current").toInt();
   QDomElement elt = myelem.firstChildElement();
   for(int i=0; i < current; i++) elt = elt.nextSiblingElement();
   if(elt.attribute("element")=="panel"||elt.attribute("element")=="page"){
     myelem.setAttribute("status", elt.attribute("status"));
     myelem.setAttribute("currentText", elt.attribute("currentText"));
     
     }else{
       myelem.setAttribute("status", "done");
     }
   emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
    
 }










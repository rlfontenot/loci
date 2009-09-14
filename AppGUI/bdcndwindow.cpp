#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "bdcndwindow.h"
#include "pages.h"

BdCndWindow::BdCndWindow(QDomElement& theelem, QDomElement& newroot,
                         QStringList names, QTableView* boundaryView,
                         QWidget* parent):GeneralWindow(theelem, newroot, parent){
 
  
  bdNames = names;
  typesWidget = new QComboBox;
  pagesWidget = new QStackedWidget;
  bdTypes.clear();
  myTableView =boundaryView;
  

  cndNode = myroot.firstChildElement("boundary_conditions");
  if(cndNode.isNull()){
    QMessageBox::warning(window(), ".xml",
                         tr("can not find element 'boundary_conditions'")
                         );
    return;
  }
  QDomElement elt_name = cndNode.firstChildElement();
  //if the node boundary_conditions has no child, create the children
  if(elt_name.isNull()){
    for(int i =0; i < bdNames.size(); i++){
      QDomElement  aNode = cndNode.ownerDocument().createElement(bdNames[i]);
      cndNode.appendChild(aNode);
    }
  }

       
      
  //now all boundary nodes are ready in cndNode   

  

   
  QDomElement elt = myelem.firstChildElement();
  if(elt.isNull()){
    QMessageBox::warning(window(), ".xml",
                         tr("'boundaryPage' has no child")
                         );
    return;
  }
 
  
  for (; !elt.isNull(); elt = elt.nextSiblingElement()) {
    bdTypes << elt.tagName();
    whatsThisList << elt.attribute("whatsThis");
    toolTipList <<elt.attribute("toolTip");
    QWidget* bdCndPage = 0;
    if(elt.attribute("element")=="panel") bdCndPage = new OptionPage( elt, myroot);
    else{
      QMessageBox::warning(window(), ".xml",
                           tr("Don't know how to handle it yet ")+ elt.tagName() + " " +elt.attribute("element")
                           );
      return;
    }
      
    pagesWidget->addWidget(bdCndPage);
    connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateConditionView(const QString&)));
    connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
    connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged())); 
  }
  
  typesWidget->addItems(bdTypes);
 
  connect(typesWidget,
          SIGNAL(currentIndexChanged(int)),
          this, SLOT(changePage(int)));


  if(!(whatsThisList[typesWidget->currentIndex()].isNull()))typesWidget->setWhatsThis(whatsThisList[typesWidget->currentIndex()]);
  if(!(toolTipList[typesWidget->currentIndex()].isNull()))typesWidget->setToolTip(toolTipList[typesWidget->currentIndex()]);
     
  QGroupBox *tableBox = new QGroupBox("please select one boundary:");
  QHBoxLayout *tableLayout  = new QHBoxLayout;
  tableLayout->addWidget(myTableView);
  tableBox->setLayout(tableLayout);
  tableBox->setFlat(true);
  
 
  QHBoxLayout *typesLayout  = new QHBoxLayout;
  QLabel *label = new QLabel(myelem.attribute("label"));
  QFont font;
  font.setBold(true);
  label->setFont(font);
  typesLayout->addWidget(label);
  typesLayout->addWidget(typesWidget);
 
  updateConditionView();  
    
  QVBoxLayout *mainLayout = new QVBoxLayout;

  mainLayout->addWidget(tableBox);
  mainLayout->addSpacing(10);
  mainLayout->addLayout(typesLayout);
 
  mainLayout->addSpacing(10);
  mainLayout->addWidget(pagesWidget);
  
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));


  
  
}
void BdCndWindow::updateConditionView(const QString& text){
  QModelIndex index = myTableView->currentIndex();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  if(index.row() >=0){
    QModelIndex index2 = model->index(index.row(), 4);
    model->setData(index2, text);
    myTableView->resizeRowsToContents();
    myTableView->resizeColumnsToContents();
  }else{
    QMessageBox::warning(window(), "boundary condition setup",
                         tr("Please select a boundary first")
                         );
    return;
  }

  QDomElement  elem = cndNode.firstChildElement();
  for(int i =0; i < index.row(); i++)elem = elem.nextSiblingElement();
  QModelIndex index3 = model->index(index.row(), 3);
  QString tmp = elem.firstChildElement().attribute("status");
  model->setData(index3, tmp);
  myTableView->resizeRowsToContents();
  myTableView->resizeColumnsToContents();
  checkStatus();
  
}

void BdCndWindow::updateConditionView(){
  QModelIndex index = myTableView->currentIndex();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  
  QDomElement  elem = cndNode.firstChildElement();
  int i =0;
  for(; !elem.isNull(); elem = elem.nextSiblingElement(), i++){
    
    if(!elem.firstChildElement().isNull()){ 
  
      QModelIndex index2 = model->index(i, 4);
      model->setData(index2, elem.firstChildElement().attribute("currentText"));
       QModelIndex index3 = model->index(i, 3);
      model->setData(index3, elem.firstChildElement().attribute("status"));
      myTableView->resizeRowsToContents();
      myTableView->resizeColumnsToContents();
    }
  }
}

 void BdCndWindow::checkStatus(){
  int count = 0;
  int count_done = 0;
  QString infix = ", \n";
  QString prefix = ": < \n " ;
  QString postfix = " > \n";
  QString text = cndNode.tagName() + prefix;
  
  for( QDomElement elt = cndNode.firstChildElement(); !elt.isNull();
       elt = elt.nextSiblingElement()){
    if(elt.firstChildElement().attribute("status")=="done")count_done++;
    text += elt.tagName() + "=" +elt.firstChildElement().attribute("currentText");
    if(elt != cndNode.lastChildElement())text += infix; 
    count++;
  }
  text += postfix;
  myelem.setAttribute("currentText", text);
  
  if(count==count_done) myelem.setAttribute("status", "done");
  else myelem.setAttribute("status", QString("%1  out of ").arg(count_done)+QString(" %1 finished").arg(count));
  emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
}
 

void BdCndWindow::changePage(int index)
{
  pagesWidget->setCurrentIndex(index);
  //replace the child node of current boundary

  //get current boundary node
  int bdCurrent = cndNode.attribute("currentIndex").toInt();
  QDomElement elt = cndNode.firstChildElement();
   //
   for(int i = 0; i< bdCurrent; i++){
     elt = elt.nextSiblingElement();
   }
   if(elt.isNull()){
     QMessageBox::warning(window(), "bounary condition",
                          tr("index %d out of range").arg(bdCurrent)
                          );
     return;
   }

   //get current boundary types node
   QDomElement pageNode = myelem.firstChildElement();
   for(int i=0; i < index; i++)pageNode = pageNode.nextSiblingElement();
   if(pageNode.isNull()){
     QMessageBox::warning(window(), "bounary condition",
                          tr("index %d out of range").arg(index)
                          );
     return;
   }

   //copy the current type Node to the child of current boundary node
   QDomNode newNode = pageNode.cloneNode(true);
   if(newNode.isNull()){
     QMessageBox::warning(window(), "bounary condition",
                          tr("cloneNode failed")
                          );
     return;
   }

   QDomNode tmpNode;
   if(elt.firstChildElement().isNull()){
     tmpNode = elt.appendChild(newNode);
     if(tmpNode.isNull()){
     QMessageBox::warning(window(), "main xml file",
                          tr("replace ")+elt.tagName()+tr(" node failed")
                          );
     return ;
   }  
   } else{
     QDomElement aNode = elt.firstChildElement();
     if(aNode.tagName() != newNode.toElement().tagName()){
       tmpNode = elt.replaceChild( newNode, aNode);
       if(tmpNode.isNull()){
         QMessageBox::warning(window(), "main xml file",
                              tr("replace ")+elt.tagName()+tr(" node failed")
                              );
         return ;
       } 
     }
   }
  
   
   QDomElement typeNode = elt.firstChildElement();
   //   //replace the stackedWidget to current set
   OptionPage* bdCndPage=0;
   if(typeNode.attribute("element")=="panel") bdCndPage = new OptionPage( typeNode, myroot);
   else{
     QMessageBox::warning(window(), "main xml file",
                          tr("Don't know how to handle it yet: ")+typeNode.tagName() + " " + typeNode.attribute("element")
                          );
     return ;
   }
   pagesWidget->removeWidget(pagesWidget->currentWidget());
   pagesWidget->insertWidget(index, bdCndPage);
   pagesWidget->setCurrentIndex(index);
   connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateConditionView(const QString&)));
   connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
   connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
   
   QString tmp =  bdCndPage->currentText();
   updateConditionView(tmp);

  
    if(!(whatsThisList[index].isNull()))typesWidget->setWhatsThis(whatsThisList[index]);
   else typesWidget->setWhatsThis("");
    if(!(toolTipList[index].isNull()))typesWidget->setToolTip(toolTipList[index]);
    else typesWidget->setToolTip("");
}

void  BdCndWindow::setCurrent(QModelIndex index){

  int bdCurrent=index.row();//current boundary index
  cndNode.setAttribute("currentIndex", bdCurrent);
 
  
  //go to current boundary
  QDomElement elt = cndNode.firstChildElement();
  for(int i = 0; i< bdCurrent; i++){
    elt = elt.nextSiblingElement();
  }
  if(elt.isNull()){
    QMessageBox::warning(window(), "bounary condition",
                         tr("index %d out of range").arg(bdCurrent)
                         );
    return;
  }
 
  //if currrent boudnary condition is not set, assign value 'notset'
  if(elt.firstChildElement().isNull()){
    
    typesWidget->setCurrentIndex(0);

  }else{//otherwise
    //find the index of boudary type
    QDomElement pageNode = elt.firstChildElement();
    int ind = typesWidget->findText(pageNode.tagName());
    if(ind ==-1){
      QMessageBox::warning(window(), "boundary condition setup",
                           tr("invalid boundary type ")+pageNode.tagName()
                           );
      return;
    }
    typesWidget->setCurrentIndex(ind);
        
  
    //replace the stackedWidget to current set
    QWidget* newPage=0;
    if(pageNode.attribute("element")=="panel") newPage = new OptionPage(pageNode, myroot);
    else{
      QMessageBox::warning(window(), "xml",
                           tr("Don't know how to handle it yet: ")+pageNode.tagName()+" "+ pageNode.attribute("element")
                           );
      return;
    }
    QWidget* oldPage = pagesWidget->currentWidget();
    pagesWidget->removeWidget(oldPage);
    pagesWidget->insertWidget(ind, newPage);
    pagesWidget->setCurrentWidget(newPage);
    
    connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateConditionView(const QString&)));
     connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
    connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
  }
  
}

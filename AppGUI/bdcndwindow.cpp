#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "bdcndwindow.h"
#include "pages.h"

BdCndWindow::BdCndWindow(QDomElement& theelem, QDomElement& newroot,
                         QStringList names, 
                         QWidget* parent):GeneralWindow(theelem, newroot, parent){
 
  
  bdNames = names;
  typesWidget = new QComboBox;
  pagesWidget = new QStackedWidget;
  bdTypes.clear();
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
    cndNode.setAttribute("currentIndex", "0");
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
  }
  
  typesWidget->addItems(bdTypes);
 
  connect(typesWidget,
          SIGNAL(currentIndexChanged(int)),
          this, SLOT(changePage(int)));
  

  if(!(whatsThisList[typesWidget->currentIndex()].isNull()))typesWidget->setWhatsThis(whatsThisList[typesWidget->currentIndex()]);
  if(!(toolTipList[typesWidget->currentIndex()].isNull()))typesWidget->setToolTip(toolTipList[typesWidget->currentIndex()]);
  
  QGroupBox *tableBox = new QGroupBox("current boundary:");
  QHBoxLayout *tableLayout  = new QHBoxLayout;
  currentBdry = new QLabel("");
  tableLayout->addWidget(currentBdry);
  tableBox->setLayout(tableLayout);
  tableBox->setFlat(true);
  
 
  QHBoxLayout *typesLayout  = new QHBoxLayout;
  QLabel *label = new QLabel(myelem.attribute("label"));
  QFont font;
  font.setBold(true);
  label->setFont(font);
  typesLayout->addWidget(label);
  typesLayout->addWidget(typesWidget);

  int current = cndNode.attribute("currentIndex").toInt();
  setCurrent(current);
  
  
  emit updateConditionView();  
    
  QVBoxLayout *mainLayout = new QVBoxLayout;

  mainLayout->addWidget(tableBox);
  mainLayout->addLayout(typesLayout);
  mainLayout->addWidget(pagesWidget);
 
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  
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
 
//change boundary type
void BdCndWindow::changePage(int index)
{
     //get current boundary node
    int bdCurrent = cndNode.attribute("currentIndex").toInt();
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
  
    //nonset node can not be deep cloned
     if(newNode.isNull()){
         QMessageBox::warning(window(), "bounary condition",
                              tr("cloneNode failed ")+pageNode.tagName()
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
   }else{
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
   //replace the stackedWidget to current set

   while(pagesWidget->count() != 0){
     QWidget* oldPage = pagesWidget->currentWidget();
     pagesWidget->removeWidget(oldPage);
     if(oldPage){
       delete oldPage;
       oldPage = 0;
     }
   }
   
   if(typeNode.attribute("element")=="panel"){
     OptionPage* newPage = new OptionPage( typeNode, myroot);
     pagesWidget->insertWidget(0, newPage);
     pagesWidget->setCurrentWidget(newPage);
     connect(newPage, SIGNAL(textChanged(const QString&)), this, SIGNAL(updateConditionView()));
     connect(this, SIGNAL(updateConditionView()), this, SLOT(checkStatus())); 
     connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
     connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
     connect(this, SIGNAL(showStatus(const bool&)), newPage, SLOT(updateShowStatus(const bool&)));
   }else{
     QMessageBox::warning(window(), "main xml file",
                          tr("Don't know how to handle it yet: ")+typeNode.tagName() + " " + typeNode.attribute("element")
                          );
     return ;
   }

  
   
   
     
   //QString tmp =  bdCndPage->currentText();
   //updateConditionView(tmp);
   emit updateConditionView();
  
    if(!(whatsThisList[index].isNull()))typesWidget->setWhatsThis(whatsThisList[index]);
    else typesWidget->setWhatsThis("");
    if(!(toolTipList[index].isNull()))typesWidget->setToolTip(toolTipList[index]);
    else typesWidget->setToolTip("");

}
//change boundary 


void  BdCndWindow::setCurrent(int bdCurrent){
 
  int previous  =  cndNode.attribute("currentIndex").toInt();
  cndNode.setAttribute("currentIndex", bdCurrent);
  if(bdCurrent >=0 )currentBdry->setText(bdNames[bdCurrent]);
  if(previous != bdCurrent){
  
   
  
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
    if(elt.firstChildElement().isNull() || elt.firstChildElement().tagName()=="notset"){
      copyBdCnd(previous, bdCurrent);
      
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
    }
  }
  emit updateConditionView();
}


void BdCndWindow::copyBdCnd(int previous, int current){
  QDomElement pn = cndNode.firstChildElement();
  for(int i = 0; i< previous; i++){
      pn = pn.nextSiblingElement();
  }
  if(pn.isNull()){
    QMessageBox::warning(window(), "bounary condition",
                           tr("index %d out of range").arg(previous)
                         );
    return;
    } 
  
  QDomElement cn = cndNode.firstChildElement();
  for(int i = 0; i< current; i++){
    cn = cn.nextSiblingElement();
  }
  if(cn.isNull()){
    QMessageBox::warning(window(), "bounary condition",
                           tr("index %d out of range").arg(current)
                         );
      return;
    } 

  //copy the current type Node to the child of current boundary node

  if(pn.firstChildElement().isNull())return;
  QDomNode newNode = pn.firstChildElement().cloneNode(true);
  
  if(newNode.isNull()){
    QMessageBox::warning(window(), "bounary condition, copyNdCnd",
                         tr("cloneNode failed") + pn.firstChildElement().tagName()
                         );
    return;
  }
   
  QDomNode tmpNode;
  if(cn.firstChildElement().isNull()){
    tmpNode = cn.appendChild(newNode);
    if(tmpNode.isNull()){
      QMessageBox::warning(window(), "main xml file",
                           tr("replace ")+cn.tagName()+tr(" node failed")
                           );
      return ;
    }  
  }else{
    QDomElement aNode = cn.firstChildElement();
     if(aNode.tagName() =="notset"){
       tmpNode = cn.replaceChild( newNode, aNode);
       if(tmpNode.isNull()){
         QMessageBox::warning(window(), "main xml file",
                              tr("replace ")+cn.tagName()+tr(" node failed")
                              );
         return ;
       } 
     }
  }
   
   while(pagesWidget->count() != 0){
     QWidget* oldPage = pagesWidget->currentWidget();
     pagesWidget->removeWidget(oldPage);
     if(oldPage){
       delete oldPage;
       oldPage = 0;
     }
   }
   
   if(cn.firstChildElement().attribute("element")=="panel"){
     QDomElement elm = cn.firstChildElement();
     OptionPage* newPage = new OptionPage(elm, myroot);
     pagesWidget->insertWidget(0, newPage);
     pagesWidget->setCurrentWidget(newPage);
     connect(newPage, SIGNAL(textChanged(const QString&)), this, SIGNAL(updateConditionView()));
     connect(this, SIGNAL(updateConditionView()), this, SLOT(checkStatus())); 
     connect(this, SIGNAL(stateChanged()), newPage, SLOT(changeState()));
     connect(this, SIGNAL(componentsChanged()), newPage, SIGNAL(componentsChanged()));
     connect(this, SIGNAL(showStatus(const bool&)), newPage, SLOT(updateShowStatus(const bool&)));
   }else{
     QMessageBox::warning(window(), "main xml file",
                          tr("Don't know how to handle it yet: ")+cn.firstChildElement().tagName() + " " + cn.firstChildElement().attribute("element")
                          );

   }
}

                            
void  BdCndWindow::setCurrent(QModelIndex index){

  int bdCurrent=index.row();//current boundary index
  // cndNode.setAttribute("currentIndex", bdCurrent);
//   if(bdCurrent >=0 )currentBdry->setText(bdNames[bdCurrent]);
  
//   //go to current boundary
//   QDomElement elt = cndNode.firstChildElement();
//   for(int i = 0; i< bdCurrent; i++){
//     elt = elt.nextSiblingElement();
//   }
//   if(elt.isNull()){
//     QMessageBox::warning(window(), "bounary condition",
//                          tr("index %d out of range").arg(bdCurrent)
//                          );
//     return;
//   }
 
//   //if currrent boudnary condition is not set, assign value 'notset'
//   if(elt.firstChildElement().isNull()){
//     typesWidget->setCurrentIndex(0);
//   }else{//otherwise
//     //find the index of boudary type
//     QDomElement pageNode = elt.firstChildElement();
//     int ind = typesWidget->findText(pageNode.tagName());
//     if(ind ==-1){
//       QMessageBox::warning(window(), "boundary condition setup",
//                            tr("invalid boundary type ")+pageNode.tagName()
//                            );
//       return;
//     }
//     typesWidget->setCurrentIndex(ind);
        
  
     
//   }
  setCurrent(bdCurrent);
  
  
}



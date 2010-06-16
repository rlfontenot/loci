#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "bdcndwindow.h"
#include "pages.h"

BdCndWindow::BdCndWindow(QDomElement& theelem,
                         QPointer<GLViewer> theviewer,
                         QWidget* parent):GeneralGroup(theelem, parent),viewer(theviewer){
  
  

  QDomElement myroot = myelem.ownerDocument().documentElement();
  typesWidget = new QComboBox;
  pagesWidget = new QStackedWidget;
  bdTypes.clear();
 
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
  
  QDomElement cndNode = myelem.ownerDocument().documentElement().firstChildElement("boundary_conditions");
  int current = cndNode.attribute("currentIndex").toInt();
  setCurrent(current);
  
  
  updateConditionView();  
    
  QVBoxLayout *mainLayout = new QVBoxLayout;
  if(selectBoundary()) mainLayout->addWidget(boundaryView);
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
  QDomElement cndNode = myelem.ownerDocument().documentElement().firstChildElement("boundary_conditions");
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
  QDomElement cndNode = myelem.ownerDocument().documentElement().firstChildElement("boundary_conditions");
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
     VarPanel* newPage = new VarPanel( typeNode);
     pagesWidget->insertWidget(0, newPage);
     pagesWidget->setCurrentWidget(newPage);
     connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateConditionView()));
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
   updateConditionView();
  
    if(!(whatsThisList[index].isNull()))typesWidget->setWhatsThis(whatsThisList[index]);
    else typesWidget->setWhatsThis("");
    if(!(toolTipList[index].isNull()))typesWidget->setToolTip(toolTipList[index]);
    else typesWidget->setToolTip("");

}
//change boundary 


void BdCndWindow::selectCurrent(int row){
  
  QColor value= modBoundaries->item(row, 0)->background().color();
  QModelIndex index =qobject_cast<const QAbstractItemModel*>(modBoundaries)->index(row, 1);
  qobject_cast<QAbstractItemView*>( boundaryView)->setCurrentIndex(index);
  emit setCurrent(index);
    
  viewer->setCurrentObj(row, value);
  
}


void BdCndWindow::showBoundary(QModelIndex top, QModelIndex ){
 
  if(top.column() ==2){//visibility item
    QString value = top.data(Qt::EditRole).toString();
    if(value == "show"){
      viewer->setVisibility(top.row(),true);
    }else{
      viewer->setVisibility(top.row(),false);
    }
  }
  else if(top.column() ==0) {//color item
   
    QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), top.column())->background().color();
    
    viewer->setCurrentObj(top.row(), value);  
  }
  

}

void BdCndWindow::setCurrentObj(QModelIndex top){
  QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), 0)->background().color();
  viewer->setCurrentObj(top.row(), value);
}

void BdCndWindow::updateBoundaryView(){
  //clean up the data
  if(modBoundaries){
    delete modBoundaries;
    modBoundaries = 0;
  }
    
  if(boundaryView){
    delete boundaryView;
    boundaryView = 0;
  }
      // update status
  updateStatusTip(myelem.attribute("buttonIndex").toInt());
  

  if(boundaryView == 0){
    selectBoundary();
    updateConditionView();
  }
}
void BdCndWindow::selectBoundaryPressed(){

  if(boundaryView == 0){
    selectBoundary();
  }
}
      
void BdCndWindow::updateConditionView(){
  if(boundaryView==0)return;
  QModelIndex index = boundaryView->currentIndex();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  QDomElement myroot = myelem.ownerDocument().documentElement();
  QDomElement cndNode = myroot.firstChildElement("boundary_conditions");
  if(cndNode.isNull()){
    QMessageBox::warning(window(), ".xml",
                         tr("can not find element 'boundary_conditions'")
                         );
    return;
  }
  QDomElement  elem = cndNode.firstChildElement();
  int i =0;
  for(; !elem.isNull(); elem = elem.nextSiblingElement(), i++){
    
    if(!elem.firstChildElement().isNull()){ 
      
      QModelIndex index2 = model->index(i, 4);
      model->setData(index2, elem.firstChildElement().attribute("currentText"));
      QModelIndex index3 = model->index(i, 3);
      model->setData(index3, elem.firstChildElement().attribute("status"));
      boundaryView->resizeRowsToContents();
      boundaryView->resizeColumnsToContents();
    }
  }
  checkStatus();
}  

bool BdCndWindow::selectBoundary(){
  
  QDomElement theroot = myelem.ownerDocument().documentElement();
  theroot = theroot.firstChildElement("mainWindow");
  QDomElement elem_bdname = theroot.firstChildElement("gridSetup");
  if(!elem_bdname.hasAttribute("boundary_names")){
    
    QMessageBox::warning(window(), tr("main, boundary.xml"),
                         tr("no boundary names, please use Grid Setup to load grid first")
                         );
    return false;
  }
  QStringList bdnames = elem_bdname.attribute("boundary_names").split(",");
  
    
  // Get boundary names from topo file
  if(bdnames.empty()){
    QMessageBox::warning(this, tr("select Boundary"),
                         tr("please  load boundaries first"));
    return false;
  }
  
  if(modBoundaries){
    delete modBoundaries;
    modBoundaries = 0;
  }
  if(boundaryView) {
    delete boundaryView;
    boundaryView = 0;
  }
  // Load information into data model
  modBoundaries = new QStandardItemModel(bdnames.size(), 5, this);
    
  modBoundaries->setHeaderData(0, Qt::Horizontal, QObject::tr("color"));
  modBoundaries->setHeaderData(1, Qt::Horizontal, QObject::tr("boundary name"));
  modBoundaries->setHeaderData(2, Qt::Horizontal, QObject::tr("show/hide"));
  modBoundaries->setHeaderData(3, Qt::Horizontal, QObject::tr("setup status"));
  modBoundaries->setHeaderData(4, Qt::Horizontal, QObject::tr("boundary conditions"));

    
    
  theroot = myelem.ownerDocument().documentElement();
  theroot = theroot.firstChildElement("boundary_conditions");
  for (int i = 0; i < bdnames.size(); ++i) {
    QColor newColor = default_color[i%12];
    QStandardItem* colorItem = new QStandardItem("");
    QStandardItem* nameItem = new QStandardItem(bdnames[i]);
    QStandardItem* showItem = new QStandardItem("show");
    QStandardItem* statusItem = new QStandardItem("not setup");
    QStandardItem* conditionItem = new QStandardItem("");
    colorItem->setBackground(QBrush(newColor));
      
    nameItem->setFlags(Qt::ItemIsSelectable | 
                       Qt::ItemIsUserCheckable | 
                       Qt::ItemIsEnabled);
      
    modBoundaries->setItem(i, 0, colorItem);
    modBoundaries->setItem(i, 1, nameItem);
    modBoundaries->setItem(i, 2, showItem);
    modBoundaries->setItem(i, 3, statusItem);
    modBoundaries->setItem(i, 4, conditionItem); 
  }
    
    

  QItemSelectionModel *selections = new QItemSelectionModel(modBoundaries, this);

 

  // Construct and show dock widget
  showDelegate* delBoundaries = new showDelegate(this);
  colorDelegate* delColor = new colorDelegate(this);
 

  boundaryView = new QTableView(this);

  boundaryView->setModel(modBoundaries);
  boundaryView->setSelectionModel(selections);
  boundaryView->setSelectionMode(QAbstractItemView::SingleSelection);
  boundaryView->setItemDelegateForColumn(2,delBoundaries);
  boundaryView->setItemDelegateForColumn(0,delColor);
 
  
  boundaryView->setColumnWidth(0, 20);
  boundaryView->setWordWrap(false);
  boundaryView->resizeRowsToContents();
  boundaryView->resizeColumnsToContents();
 
  
  connect(modBoundaries, SIGNAL(dataChanged( const QModelIndex&, const QModelIndex&)),
          this, SLOT(showBoundary(QModelIndex, QModelIndex)));
  connect(boundaryView, SIGNAL(clicked( const QModelIndex&)),
          this, SLOT(setCurrent(QModelIndex)));
  connect(boundaryView, SIGNAL(clicked( const QModelIndex&)),
          this, SLOT(setCurrentObj(QModelIndex)));

  connect(viewer, SIGNAL(pickCurrent(int)), this, SLOT(selectCurrent(int)));
 

  theroot = myelem.ownerDocument().documentElement();
  theroot = theroot.firstChildElement("boundary_conditions");
  int currentIndex = theroot.attribute("currentIndex").toInt();
  selectCurrent(currentIndex);


  return true;
}








void  BdCndWindow::setCurrent(int bdCurrent){
  QDomElement cndNode = myelem.ownerDocument().documentElement().firstChildElement("boundary_conditions"); 
  int previous  =  cndNode.attribute("currentIndex").toInt();
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
    
    currentBdry->setText(elt.tagName());
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

    updateConditionView();
}


void BdCndWindow::copyBdCnd(int previous, int current){
  QDomElement cndNode = myelem.ownerDocument().documentElement().firstChildElement("boundary_conditions");
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
     VarPanel* newPage = new VarPanel(elm);
     pagesWidget->insertWidget(0, newPage);
     pagesWidget->setCurrentWidget(newPage);
     connect(newPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateConditionView()));
    
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



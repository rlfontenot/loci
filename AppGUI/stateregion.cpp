#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <stdlib.h>
#include <unistd.h>
#include "stateregion.h"
#include "pages.h"
using namespace std;

void text2tree(QStringList &regions, QTreeWidget* tree){
  
  if(regions.size()==0) return;
  QTreeWidgetItem* root = tree->topLevelItem(0);
  
  for(int i = 0; i < regions.size(); i++){
    QString shape = regions[i].section("(",0, 0);
    QString par = regions[i].section("(",1, -1).section(",composition=",0, 0);
    QString composition = regions[i].section("(",1, -1).section(",composition=",-1, -1).section(")", 0, 0);

    //build a node
    QTreeWidgetItem* objItem = new QTreeWidgetItem();
    objItem->setText(0, "object");
    objItem->setText(1, composition);
    objItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled|Qt::ItemIsDragEnabled);
    QTreeWidgetItem* newItem = new QTreeWidgetItem(objItem);
    newItem->setText(0, "shape");
    newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    QTreeWidgetItem* valueNode = new QTreeWidgetItem(newItem);
    valueNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    root->addChild(objItem);
    
    if(shape=="inSphere"){
      QString r= par.section("radius=", -1, -1).section(",", 0, 0);
      QString center = par.section("center=[", -1, -1).section("]", 0, 0);
      QStringList xyz = center.split(",", QString::SkipEmptyParts);
      if(xyz.size()!=3){
        qDebug() << "can not get center of the sphere";
        return ;
      }
          
      valueNode->setText(0, "sphere");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "x0");
      x0Node->setText(1, xyz[0]);
           
      QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
      y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      y0Node->setText(0, "y0");
      y0Node->setText(1, xyz[1]);
      QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
      z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      z0Node->setText(0, "z0");
      z0Node->setText(1, xyz[2]);
      QTreeWidgetItem* rNode = new QTreeWidgetItem(valueNode);
      rNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      rNode->setText(0, "r");
      rNode->setText(1, r);
    }else if(shape=="inCylinder"){
      QString r= par.section("radius=", -1, -1).section(",", 0, 0);
      QString p1 = par.section("p1=[", -1, -1).section("]", 0, 0);
      QStringList xyz1 = p1.split(",", QString::SkipEmptyParts);
      if(xyz1.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      QString p2 = par.section("p2=[", -1, -1).section("]", 0, 0);
      QStringList xyz2 = p2.split(",", QString::SkipEmptyParts);
      if(xyz2.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      valueNode->setText(0, "cylinder");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x1");
        x0Node->setText(1, xyz1[0]);
           
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y1");
        y0Node->setText(1, xyz1[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z1");
        z0Node->setText(1, xyz1[2]);
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x2");
        x0Node->setText(1, xyz2[0]);
             
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y2");
        y0Node->setText(1, xyz2[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z2");
        z0Node->setText(1, xyz2[2]);
      }
           
           
      QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
      z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      z1Node->setText(0, "r");
      z1Node->setText(1, r);
           
    }else if(shape=="inCone"){
      QString r1= par.section("r1=", -1, -1).section(",", 0, 0);
      QString r2= par.section("r2=", -1, -1).section(",", 0, 0);
      QString p1 = par.section("p1=[", -1, -1).section("]", 0, 0);
      QStringList xyz1 = p1.split(",", QString::SkipEmptyParts);
      if(xyz1.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      QString p2 = par.section("p2=[", -1, -1).section("]", 0, 0);
      QStringList xyz2 = p2.split(",", QString::SkipEmptyParts);
      if(xyz2.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      valueNode->setText(0, "cone");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x1");
        x0Node->setText(1, xyz1[0]);
           
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y1");
        y0Node->setText(1, xyz1[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z1");
        z0Node->setText(1, xyz1[2]);
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x2");
        x0Node->setText(1, xyz2[0]);
             
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y2");
        y0Node->setText(1, xyz2[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z2");
        z0Node->setText(1, xyz2[2]);
      }
           
           
      QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
      z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      z1Node->setText(0, "r1");
      z1Node->setText(1, r1);
           
      z1Node = new QTreeWidgetItem(valueNode);
      z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      z1Node->setText(0, "r2");
      z1Node->setText(1, r2);  
       
    }else if(shape=="inBox"){
       
      QString p1 = par.section("p1=[", -1, -1).section("]", 0, 0);
      QStringList xyz1 = p1.split(",", QString::SkipEmptyParts);
      if(xyz1.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      QString p2 = par.section("p2=[", -1, -1).section("]", 0, 0);
      QStringList xyz2 = p2.split(",", QString::SkipEmptyParts);
      if(xyz2.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      valueNode->setText(0, "box");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x1");
        x0Node->setText(1, xyz1[0]);
           
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y1");
        y0Node->setText(1, xyz1[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z1");
        z0Node->setText(1, xyz1[2]);
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x2");
        x0Node->setText(1, xyz2[0]);
             
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y2");
        y0Node->setText(1, xyz2[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z2");
        z0Node->setText(1, xyz2[2]);
      }
    }else if(shape=="leftPlane"){
        
      QString p1 = par.section("point=[", -1, -1).section("]", 0, 0);
      QStringList xyz1 = p1.split(",", QString::SkipEmptyParts);
      if(xyz1.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      QString p2 = par.section("normal=[", -1, -1).section("]", 0, 0);
      QStringList xyz2 = p2.split(",", QString::SkipEmptyParts);
      if(xyz2.size()!=3){
        qDebug() << "can not get p1 of the cylinder";
        return ;
      }
      valueNode->setText(0, "leftplane");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "pointx");
        x0Node->setText(1, xyz1[0]);
           
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "pointy");
        y0Node->setText(1, xyz1[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "pointz");
        z0Node->setText(1, xyz1[2]);
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "normalx");
        x0Node->setText(1, xyz2[0]);
             
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "normaly");
        y0Node->setText(1, xyz2[1]);
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "normalz");
        z0Node->setText(1, xyz2[2]);
      }
           
    }
    
  }
  tree->expandAll();
  tree->resizeColumnToContents(0);
  
}
  

void RegionWindow::addStateClicked(){
  bool ok;
  int count = stateList->count();
  QString text = QInputDialog::getText(this, tr("Please specify the name of state function"),
                                       tr("State name:"), QLineEdit::Normal,
                                       tr("state%1").arg(count), &ok);
  if (ok && !text.isEmpty()){
    //add a node
    QDomElement pageNode = myelem.firstChildElement("state").firstChildElement();
    QDomNode newNode = pageNode.cloneNode(true);
  
    //nonset node can not be deep cloned
    if(newNode.isNull()){
      QMessageBox::warning(window(), "State Window",
                           tr("cloneNode failed ")+pageNode.tagName()
                           );
      return;
    }
    newNode.toElement().setTagName(text);
     
    QDomNode elt;
    elt = myelem.firstChildElement("state").appendChild(newNode);
    if(elt.isNull()){
      QMessageBox::warning(window(), "main xml file",
                           tr("append ")+newNode.toElement().tagName()+tr(" node failed")
                           );
      return ;
    }  
    //add a page
    QDomElement eltelem = elt.toElement();
    VarPanel* thegroup = new VarPanel(eltelem, this);
    connect(this, SIGNAL(componentsChanged()), thegroup, SIGNAL(componentsChanged()));
    connect(this, SIGNAL(showStatus(const bool &)), thegroup, SLOT(updateShowStatus(const bool &))); 
    connect(this, SIGNAL(stateChanged()), thegroup, SLOT(changeState()));
    statePages->addWidget(thegroup);
    QListWidgetItem *bdCndiButton = new QListWidgetItem(stateList);
    bdCndiButton->setText(thegroup->currentText());
    bdCndiButton->setTextAlignment(Qt::AlignHCenter);
    bdCndiButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    connect(thegroup, SIGNAL(textChanged(const QString&)), this, SLOT(setItemText(const QString&)));
    connect(thegroup, SIGNAL(textChanged(const QString&)), this, SLOT(updateText()));
    stateList->setCurrentItem(bdCndiButton);
  }
  stateFunctions << text; 
}


void RegionWindow::setItemText(const QString& text){
  QListWidgetItem* item = stateList->currentItem();
  if(item)item->setText(text);
}





RegionWindow::RegionWindow(  QDomElement& elem, QWidget *parent):GeneralGroup(elem, parent){
  if(myelem.hasAttribute("whatsThis"))setWhatsThis(myelem.attribute("whatsThis"));
  if(myelem.hasAttribute("toolTip"))setToolTip(myelem.attribute("toolTip"));
  if(myelem.hasAttribute("statusTip"))setStatusTip(myelem.attribute("statusTip"));
  
  //set up state window
  stateGroup = new QGroupBox("define all state functions", this);
  {
     
    stateList = new QListWidget;
    statePages = new QStackedWidget;
     
    QDomElement elt= myelem.firstChildElement("state").firstChildElement();
    for(; !elt.isNull(); elt=elt.nextSiblingElement()){
       
      VarPanel* thegroup = new VarPanel(elt, this);
      connect(this, SIGNAL(componentsChanged()), thegroup, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool &)), thegroup, SLOT(updateShowStatus(const bool &))); 
      connect(this, SIGNAL(stateChanged()), thegroup, SLOT(changeState()));
      connect(thegroup, SIGNAL(textChanged(const QString&)), this, SLOT(updateText()));
      statePages->addWidget(thegroup);
      QListWidgetItem *bdCndiButton = new QListWidgetItem(stateList);
      bdCndiButton->setText(thegroup->currentText());
      bdCndiButton->setTextAlignment(Qt::AlignHCenter);
      bdCndiButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
      connect(thegroup, SIGNAL(textChanged(const QString&)), this, SLOT(setItemText(const QString&)));
      stateList->setCurrentItem(bdCndiButton);
      stateFunctions<<elt.tagName();
    } 

  
 
  
    QPushButton *addStateButton = new QPushButton(tr("&Add  a State"));
    QPushButton *nextButton = new QPushButton(tr("next->" ));
    connect(addStateButton, SIGNAL(clicked()), this, SLOT(addStateClicked()));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(next()));
   
   
  
    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(addStateButton);
    buttonsLayout->addWidget(nextButton);
   
    connect(stateList,
            SIGNAL(currentRowChanged(int)),
            statePages, SLOT(setCurrentIndex(int)));
   

    QVBoxLayout* stateLayout = new QVBoxLayout;
    stateLayout->addWidget(stateList);
    stateLayout->addWidget(statePages);
    stateLayout->addLayout(buttonsLayout);
    stateGroup->setLayout(stateLayout);
  }
   
  //set up region window
   
  regionGroup = new QGroupBox("region window", this);
  {
    
    QPushButton *addRegionButton = new QPushButton(tr("&Add  a Region"));
    QPushButton *assignStateButton = new QPushButton(tr("Assign a State Function"));
    QPushButton *previousButton = new QPushButton(tr("<-Previous" ));
    connect(addRegionButton, SIGNAL(clicked()), this, SLOT(addShape()));
    connect(assignStateButton, SIGNAL(clicked()), this, SLOT(assignState()));
                                                     
    connect(previousButton, SIGNAL(clicked()), this, SLOT(previous()));
    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(addRegionButton);
    buttonsLayout->addWidget(assignStateButton);
    
    tree = new QTreeWidget;
    tree->setColumnCount(2);
    QStringList header;
    header << "name" << "value";
    tree->setHeaderLabels(header);
    QTreeWidgetItem*  root =  new QTreeWidgetItem(tree);
    root->setText(0, "region");
    tree->addTopLevelItem(root);
    tree->expandAll();
    tree->setSelectionMode(QAbstractItemView::SingleSelection);
    tree->viewport()->setAcceptDrops(true);
    tree->setDragEnabled(true);
    tree->setAcceptDrops(true);
    tree->setDropIndicatorShown(true);
    tree->setDragDropMode(QAbstractItemView::InternalMove);

    if(myelem.hasAttribute("currentText")){
      QString infix=",&#xa;";
      if(myelem.hasAttribute("infix")) infix = myelem.attribute("infix");
      
      
      
      QString regionPrefix="regions = [";
      if(myelem.hasAttribute("regionPrefix")) regionPrefix = myelem.attribute("regionPrefix");
      QString regionPostfix = "]";
      if(myelem.hasAttribute("regionPostfix")) regionPostfix = myelem.attribute("regionPostfix");
      
      QString currentText = myelem.attribute("currentText");
      currentText = currentText.section(regionPrefix, -1, -1, QString::SectionSkipEmpty);
      currentText = currentText.section(regionPostfix, 0, -2, QString::SectionSkipEmpty);
      QStringList regions = currentText.split(infix, QString::SkipEmptyParts); 
      if(regions.size()!=0) text2tree(regions, tree);
    }
    
    
    defaultShapes.push_back(0);
    //SPHERE
    {
      vector<double> vdefault(4);
      vdefault[0] = 0;
      vdefault[1] = 0;
      vdefault[2] = 0;
      vdefault[3] = 1;
      
      Shape* shape = new Shape(SPHERE, vdefault);
      defaultShapes.push_back(shape);
    }
    //CONE
    
    {
      vector<double> vdefault(8);
      vdefault[0] = 0;
      vdefault[1] = 0;
      vdefault[2] = 0;
      vdefault[3] = 0;
      vdefault[4] = 0;
      vdefault[5] = 1;
      vdefault[6] = 0;
      vdefault[7] = 1;
      Shape* shape = new Shape(CONE, vdefault);
      defaultShapes.push_back(shape);
    }
    //CYLINDER
    {
      vector<double> vdefault(7);
      vdefault[0] = 0;
      vdefault[1] = 0;
      vdefault[2] = 0;
      vdefault[3] = 0;
      vdefault[4] = 0;
      vdefault[5] = 1;
      vdefault[6] = 1;
      
      Shape* shape = new Shape(CYLINDER, vdefault);
      defaultShapes.push_back(shape);
    }
    
    //BOX
    {
      vector<double> vdefault(6);
      vdefault[0] = 0;
      vdefault[1] = 0;
      vdefault[2] = 0;
      vdefault[3] = 1;
      vdefault[4] = 1;
      vdefault[5] = 1;
      Shape* shape = new Shape(BOX, vdefault);
      defaultShapes.push_back(shape);
    }
    //PLANES
    {
      vector<double> vdefault(6);
      vdefault[0] = 0;
      vdefault[1] = 0;
      vdefault[2] = 0;
      vdefault[3] = 0;
      vdefault[4] = 0;
      vdefault[5] = 1;
      Shape* shape = new Shape(LEFTPLANE, vdefault);
      defaultShapes.push_back(shape);
    }


    regionPages = new QStackedWidget;

    for(unsigned int i =0; i< defaultShapes.size(); i++){
      ParaPage *paraPage = new ParaPage(defaultShapes[i]);
      regionPages->addWidget(paraPage);
      connect(paraPage, SIGNAL(valueChanged()), this, SLOT(updateShape()));     
    }

    connect(tree, SIGNAL(itemClicked(QTreeWidgetItem*, int)),
            this, SLOT(showData(QTreeWidgetItem*)));
   
    QHBoxLayout* objLayout = new QHBoxLayout; 
    objLayout->addWidget(regionPages);
    objLayout->addWidget(tree,2);
    
    QVBoxLayout *regionLayout = new QVBoxLayout;
    regionLayout->addLayout(objLayout);
    regionLayout->addStretch(1);
    regionLayout->addLayout(buttonsLayout);
    regionGroup->setLayout(regionLayout);
  }





  textEdit = new QTextEdit;
  
  myPages = new QStackedWidget(this);
  myPages->addWidget(stateGroup);
  myPages->addWidget(regionGroup);
  myPages->setCurrentIndex(0);
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(myPages);
  
  mainLayout->addWidget(textEdit);
  
  setLayout(mainLayout);
  emit valueChanged(getRoot());
  updateText();
}

void RegionWindow::next(){
  myPages->setCurrentIndex(1);
}
void RegionWindow::previous(){
  myPages->setCurrentIndex(0);
}

QString tree2str(const QTreeWidget* tree, QString prefix, QString postfix, QString infix){
  QString text = prefix;
  
  QList<QTreeWidgetItem*> nodelist= tree->findItems("object", Qt::MatchRecursive);
  
  for(int i = 0; i < nodelist.size(); i++){
    QTreeWidgetItem* item = nodelist[i]->child(0)->child(0);
    if(item){
      QString shape = item->text(0);
      if(shape== "sphere"){
        text += "inSphere(radius="+item->child(3)->text(1)+",";
        text += "center=[" +item->child(0)->text(1)+","+
          item->child(1)->text(1)+","+item->child(2)->text(1)+"],";
        text += "composition=" +nodelist[i]->text(1)+")"+infix;
       
        
      }else if(shape=="cone"){
        text += "inCone(r1="+item->child(6)->text(1)+",";
        text += "r2="+item->child(7)->text(1)+",";
        text += "p1=[" +item->child(0)->text(1)+","+
          item->child(1)->text(1)+","+item->child(2)->text(1)+"],";
        text += "p2=[" +item->child(3)->text(1)+","+
          item->child(4)->text(1)+","+item->child(5)->text(1)+"],";
        text += "composition=" +nodelist[i]->text(1)+")"+infix;
          
        
        
      }else if(shape=="cylinder"){
        text += "inCylinder(radius="+item->child(6)->text(1)+",";
        text += "p1=[" +item->child(0)->text(1)+","+
          item->child(1)->text(1)+","+item->child(2)->text(1)+"],";
        text += "p2=[" +item->child(3)->text(1)+","+
          item->child(4)->text(1)+","+item->child(5)->text(1)+"],";
        text += "composition=" +nodelist[i]->text(1)+")"+infix;

        
        
      }else if(shape == "box"){
        text += "inBox(";
        text += "p1=[" +item->child(0)->text(1)+","+
          item->child(1)->text(1)+","+item->child(2)->text(1)+"],";
        text += "p2=[" +item->child(3)->text(1)+","+
          item->child(4)->text(1)+","+item->child(5)->text(1)+"],";
        text += "composition=" +nodelist[i]->text(1)+")"+infix;
        
      }else if(shape=="leftplane"){
        text += "leftPlane(";
        
        text += "point=[" +item->child(0)->text(1)+","+
          item->child(1)->text(1)+","+item->child(2)->text(1)+"],";
        text += "normal=[" +item->child(3)->text(1)+","+
          item->child(4)->text(1)+","+item->child(5)->text(1)+"],";
        text += "composition=" +nodelist[i]->text(1)+")"+infix;  
          
      }
    }
  }
  text.remove(text.size()-infix.size(), infix.size());
  text += postfix;
  return text;
}
  

QString RegionWindow::currentText(){
  QString text;
  
  if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";
  QString prefix ="";
  if(myelem.hasAttribute("prefix")) prefix = myelem.attribute("prefix");
  QString postfix = "";
  if(myelem.hasAttribute("postfix")) postfix = myelem.attribute("postfix");
  QString infix=",&#xa;";
  if(myelem.hasAttribute("infix")) infix = myelem.attribute("infix");
  else if(prefix=="" && postfix=="") infix = "\n";

  
  QString regionPrefix="regions = [";
  if(myelem.hasAttribute("regionPrefix")) regionPrefix = myelem.attribute("regionPrefix");
  QString regionPostfix = "]";
  if(myelem.hasAttribute("regionPostfix")) regionPostfix = myelem.attribute("regionPostfix");
  
  text = myelem.tagName()+prefix;
  
  bool done = true;
  for( QDomElement elt = myelem.firstChildElement("state").firstChildElement();
       !elt.isNull(); elt=elt.nextSiblingElement()){
    text += elt.attribute("currentText")+infix;
    if(elt.attribute("status")!="done") done=false;
  }
  if(done){
    QList<QTreeWidgetItem*> nodelist= tree->findItems("object", Qt::MatchRecursive);
    for(int i = 0; i < nodelist.size(); i++){
      QTreeWidgetItem* item = nodelist[i];
      if(item->text(1).isEmpty()){
        done = false;
        break;
      }
    }
  }
  
  text += tree2str(tree, regionPrefix, regionPostfix, infix);
  
  
  text += postfix;
  myelem.setAttribute("currentText", text);

  //check status
  if(done) myelem.setAttribute("status", "done");
  else myelem.removeAttribute("status");
  
  return text;
}
  
  
void RegionWindow::assignState(){
  
  if(tree->currentItem() == 0 || tree->currentItem()->text(0)=="region")
    {     
      QMessageBox::warning(window(),
                           tr("assign a state function"),
                           tr("Please select a 'object' node on the tree  first"));
      
      return;
    }

  if(tree->currentItem()->text(0) != "object"){
    QTreeWidgetItem* item = tree->currentItem();
    while(item->text(0)!="object"){
      item = item->parent();
    }
    tree->setCurrentItem(item);
  }
  bool ok;
  QString item = QInputDialog::getItem(this, tr("select a state function:"),
                                       tr("states:"), stateFunctions, 0, false, &ok);
  
  if (ok && !item.isEmpty()){
    
   
    tree->currentItem()->setText(1, item);
  }
  updateText();
}
  

void RegionWindow::showData(QTreeWidgetItem* item ){

  if(item->text(0)=="shape"){
    QStringList items;
    items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
      tr("leftplane") ;
  
    QString tp = item->child(0)->text(0);

    int index = items.indexOf(tp);
    index += 1;
   
    regionPages->setCurrentIndex(index);
    vector<double> para(item->child(0)->childCount());
    for(int i = 0; i< item->child(0)->childCount(); i++){
      para[i] = item->child(0)->child(i)->text(1).toDouble();
    }
    Shape* ashape = qobject_cast<ParaPage*>(regionPages->currentWidget())->shape;
    if(para.size() != ashape->para.size()) return;
    for(int i = 0; i< item->child(0)->childCount(); i++){
      ashape->para[i] = para[i];
    }
    qobject_cast<ParaPage*>(regionPages->currentWidget())->showValue();
  }

}

  
      

void RegionWindow::updateShape(){
  if(tree->currentItem() == 0 ||
     tree->currentItem()->text(0) != "shape")
    {     
      QMessageBox::warning(window(),
                           tr("update shape"),
                           tr("Please select a 'shape' node on the tree  first"));
      
      return;
    }
  QTreeWidgetItem* valueNode = tree->currentItem()->child(0);
  
  
  Shape* ashape = qobject_cast<ParaPage*>(regionPages->currentWidget())->shape;
  if(valueNode->childCount() !=(int) ashape->para.size())return;
  for(int i = 0; i < valueNode->childCount(); i++)valueNode->child(i)->setText(1,QString("%1").arg(ashape->para[i]));
  
  emit valueChanged(getRoot());
  updateText();

}



void RegionWindow::addShape(){
  QStringList items;
  items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
    tr("leftplane");

  bool ok;
  QString item = QInputDialog::getItem(this, tr("please select a geometric primitive"),
                                       tr("shapes:"), items, 0, false, &ok);
     
  if (ok && !item.isEmpty()){

    QTreeWidgetItem* objItem = new QTreeWidgetItem();
    objItem->setText(0, "object");
    objItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled|Qt::ItemIsDragEnabled);
    QTreeWidgetItem* newItem = new QTreeWidgetItem(objItem);
    newItem->setText(0, "shape");
    newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    QTreeWidgetItem* valueNode = new QTreeWidgetItem(newItem);
    valueNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    if(tree->currentItem()!=0 && tree->currentItem()->text(0) =="region"){
      tree->currentItem()->addChild(objItem);
         
    }else{
      getRoot()->addChild(objItem);
    }
    tree->setCurrentItem(newItem);
    int index = items.indexOf(item);
    index += 1;
    regionPages->setCurrentIndex(index);
    Shape* ashape = qobject_cast<ParaPage*>(regionPages->currentWidget())->shape;
    switch(ashape->tp){
    case  SPHERE:
      {
        valueNode->setText(0, "sphere");
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x0");
        x0Node->setText(1, QString("%1").arg(ashape->para[0]));
           
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y0");
        y0Node->setText(1, QString("%1").arg(ashape->para[1]));
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z0");
        z0Node->setText(1, QString("%1").arg(ashape->para[2]));
        QTreeWidgetItem* rNode = new QTreeWidgetItem(valueNode);
        rNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        rNode->setText(0, "r");
        rNode->setText(1, QString("%1").arg(ashape->para[3]));
      }
         
      break;
    
    case CONE:
      {
        valueNode->setText(0, "cone");
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "x1");
          x0Node->setText(1, QString("%1").arg(ashape->para[0]));
           
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "y1");
          y0Node->setText(1, QString("%1").arg(ashape->para[1]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "z1");
          z0Node->setText(1, QString("%1").arg(ashape->para[2]));
        }
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "x2");
          x0Node->setText(1, QString("%1").arg(ashape->para[3]));
           
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "y2");
          y0Node->setText(1, QString("%1").arg(ashape->para[4]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "z2");
          z0Node->setText(1, QString("%1").arg(ashape->para[5]));
        }






           
        QTreeWidgetItem* rNode = new QTreeWidgetItem(valueNode);
        rNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        rNode->setText(0, "r1");
        rNode->setText(1, QString("%1").arg(ashape->para[6]));
           
        QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
        z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z1Node->setText(0, "r2");
        z1Node->setText(1, QString("%1").arg(ashape->para[7]));
      
      }





      
      break;
         
    case CYLINDER:
      {
        valueNode->setText(0, "cylinder");
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "x1");
          x0Node->setText(1, QString("%1").arg(ashape->para[0]));
           
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "y1");
          y0Node->setText(1, QString("%1").arg(ashape->para[1]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "z1");
          z0Node->setText(1, QString("%1").arg(ashape->para[2]));
        }
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "x2");
          x0Node->setText(1, QString("%1").arg(ashape->para[3]));
             
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "y2");
          y0Node->setText(1, QString("%1").arg(ashape->para[4]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "z2");
          z0Node->setText(1, QString("%1").arg(ashape->para[5]));
        }
           
           
        QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
        z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z1Node->setText(0, "r");
        z1Node->setText(1, QString("%1").arg(ashape->para[6]));
      
      }
         
      break;
         
    case BOX:
      {
        valueNode->setText(0, "box");
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x1");
        x0Node->setText(1, QString("%1").arg(ashape->para[0]));
      
        QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
        y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y0Node->setText(0, "y1");
        y0Node->setText(1, QString("%1").arg(ashape->para[1]));
        QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
        z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z0Node->setText(0, "z1");
        z0Node->setText(1, QString("%1").arg(ashape->para[2]));
      
        QTreeWidgetItem* x2Node = new QTreeWidgetItem(valueNode);
        x2Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x2Node->setText(0, "x2");
        x2Node->setText(1, QString("%1").arg(ashape->para[3]));
      
        QTreeWidgetItem* y2Node = new QTreeWidgetItem(valueNode);
        y2Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        y2Node->setText(0, "y2");
        y2Node->setText(1, QString("%1").arg(ashape->para[4]));

        QTreeWidgetItem* z2Node = new QTreeWidgetItem(valueNode);
        z2Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        z2Node->setText(0, "z2");
        z2Node->setText(1, QString("%1").arg(ashape->para[5]));
      }



    
      break;
    case LEFTPLANE:
      {
        valueNode->setText(0, "leftplane");
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "pointx");
          x0Node->setText(1, QString("%1").arg(ashape->para[0]));
        
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "pointy");
          y0Node->setText(1, QString("%1").arg(ashape->para[1]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "pointz");
          z0Node->setText(1, QString("%1").arg(ashape->para[2]));
        }
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "normalx");
          x0Node->setText(1, QString("%1").arg(ashape->para[3]));
        
          QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
          y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          y0Node->setText(0, "normaly");
          y0Node->setText(1, QString("%1").arg(ashape->para[4]));
          QTreeWidgetItem* z0Node = new QTreeWidgetItem(valueNode);
          z0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          z0Node->setText(0, "normalz");
          z0Node->setText(1, QString("%1").arg(ashape->para[5]));
        }
      
      
      }
      break;
    
    default:
      break;
    }

  }
  tree->expandAll();
  tree->resizeColumnToContents(0);
  emit valueChanged(getRoot());
  updateText();    
}


QTreeWidgetItem* RegionWindow::getRoot(){
  return tree->topLevelItem(0);
}

void RegionWindow::updateText(){
  textEdit->setText(currentText());
  emit textChanged();
}

void RegionWindow::updateShowStatus(const bool& show){
  QList<QTreeWidgetItem*> nodelist= tree->findItems("object", Qt::MatchRecursive);
  for(int i = 0; i < nodelist.size(); i++){
    QTreeWidgetItem* item = nodelist[i];
    item->setBackground(1,item->background(0));
  }
 
  
  if(show){//parent says show
   
    for(int i = 0; i < nodelist.size(); i++){
      QTreeWidgetItem* item = nodelist[i];
      if(item->text(1).isEmpty()){
        QBrush brush(QColor(255, 0, 0));
        item->setBackground(1,brush);
      }
    }
  }

  
  GeneralGroup::updateShowStatus(show);  
  
}

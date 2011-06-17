/////////////////////////////////////////////////
//  Filename: vcutwindow.cpp
//
//  Contains: Implementation of VCutWindow class
/////////////////////////////////////////////////

#include <QFileDialog>
#include <QLabel>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QtDebug>
#include <QGridLayout>
#include <QCheckBox>
#include <QMessageBox>
#include <QButtonGroup>
#include <QTreeView>
#include <QListWidget>
#include <QStackedWidget>
#include <QTreeWidget>
#include <QSignalMapper>
#include <QInputDialog>
#include <QScrollArea>
#include <QMainWindow>
#include <QToolBar>
#include <QDockWidget>
#include <QTabWidget>
#include <QProcess>
#include "vcutwindow.h"
#include "grid.h"
#include "pages.h"
#include "helpwindow.h"
#include "parapage.h"
#include "transform.h"
#include "tree2doc.h"
#include "progressdialog.h"


#define PI 3.14159265358979323846264338327950





void VCutWindow::createToolBar(){
  //  int spacing =0;  
  toolbar = addToolBar(tr("tree&vis"));
 //  QGroupBox* treebar = new QGroupBox("build a  tree");
//   QHBoxLayout* barLayout = new QHBoxLayout;

//   QPushButton* addShapeButton = new QPushButton(tr("Add\nShape"), this);
//   barLayout->addWidget(addShapeButton);
//   connect(addShapeButton, SIGNAL(clicked()), this, SLOT(addShape()));

//   barLayout->addSpacing(spacing);
//   QPushButton* addTransButton = new QPushButton(tr("Add\nTransform"), this);
//   barLayout->addWidget(addTransButton);
//   connect(addTransButton, SIGNAL(clicked()), this, SLOT(addTransform()));
  
//   barLayout->addSpacing(spacing);
//   QPushButton* addOpButton = new QPushButton(tr("Add\nOperator"), this);
//   barLayout->addWidget(addOpButton);
//   connect(addOpButton, SIGNAL(clicked()), this, SLOT(addOp()));
  
  

//   barLayout->addSpacing(spacing);
//   QPushButton* removeNodeButton = new QPushButton(tr("Remove\nNode"), this);
//   barLayout->addWidget(removeNodeButton);
//   connect(removeNodeButton, SIGNAL(clicked()), this, SLOT(removeNode()));
  
 


//   barLayout->addSpacing(spacing);
//   QPushButton* addRegionButton = new QPushButton(tr("Add\nRegion "), this);
//   barLayout->addWidget(addRegionButton);
//   connect(addRegionButton, SIGNAL(clicked()), this, SLOT(addRegion()));
  
//   barLayout->addSpacing(spacing);
//   QPushButton* helpButton = new QPushButton(tr("Help"), this);
//   // barLayout->addWidget(helpButton);
//   connect(helpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
 
//   treebar->setLayout(barLayout);
//   toolbar->addWidget(treebar);
//   treebar->hide();
//   addToolBarBreak();

 //  QGroupBox* viewbar = new QGroupBox("tree view");
//   QHBoxLayout* viewLayout = new QHBoxLayout;

//   QPushButton* expandButton = new QPushButton(tr("Expand"), this);
//   viewLayout->addWidget(expandButton);
//   connect(expandButton, SIGNAL(clicked()), tree, SLOT(expandAll()));
  
//   QPushButton* collapseButton = new QPushButton(tr("Collapse"), this);
//   viewLayout->addWidget(collapseButton);
//   connect(collapseButton, SIGNAL(clicked()), tree, SLOT(collapseAll()));
  
//   QPushButton* resizeButton = new QPushButton(tr("Resize"), this);
//   viewLayout->addWidget(resizeButton);
//   connect(resizeButton, SIGNAL(clicked()), this, SLOT(resizeTree()));
  
//   viewbar->setLayout(viewLayout);
  //  toolbar->addWidget(viewbar);
  //viewbar->hide(); 

  QGroupBox* visbar = new QGroupBox(tr("Visualization"));
  QHBoxLayout* visLayout = new QHBoxLayout;
   
 
  QPushButton *clearBoundaryAct = new QPushButton(tr("Clear"), this);
  visLayout->addWidget(clearBoundaryAct);
  connect(clearBoundaryAct, SIGNAL(clicked()),
          viewer, SLOT(clearCurrent())); 
   
  QPushButton* resetAct = new QPushButton(tr("Reset"), this);
  visLayout->addWidget(resetAct);
  connect(resetAct, SIGNAL(clicked()),
          viewer, SLOT(reset()));
 
  QPushButton* fitAct = new QPushButton(tr("Fit"), this);
  visLayout->addWidget(fitAct);
  connect(fitAct, SIGNAL(clicked()),
          viewer, SLOT(fit()));
  visbar->setLayout(visLayout);

 
  toolbar->addWidget(visbar);

  // toolbar->addAction(viewerDock->toggleViewAction());
  //  toolbar->addWidget(helpButton);
}
  
void VCutWindow::resizeTree(){
  tree->resizeColumnToContents(0);
  tree->resizeColumnToContents(1);
}

  
void VCutWindow::createFlowBar(){
  int spacing =20;  
  //create flowbar
 
  
  QGroupBox* flowbar = new QGroupBox("flow bar");
  QVBoxLayout* barLayout = new QVBoxLayout;
  

  barLayout->addSpacing(spacing);
  
  QPushButton* loadButton = new QPushButton(tr("Load\nGrid"), this);
  barLayout->addWidget(loadButton);
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadGrid()));

  barLayout->addSpacing(spacing);
  QPushButton* buildButton = new QPushButton(tr("Use Toolbar\nDefine a Region"), this);
  barLayout->addWidget(buildButton);
  buildButton->setDisabled(true);
  buildButton->hide();

  barLayout->addSpacing(spacing);
  QPushButton* modifyButton = new QPushButton(tr("Use Drag&Drop\nModify the Region\n(Optional)"), this);
  barLayout->addWidget(modifyButton);
  modifyButton->setDisabled(true);
  modifyButton->hide();
    
  barLayout->addSpacing(spacing);
  QPushButton* validateButton = new QPushButton(tr("Validate\nTree"), this);
  barLayout->addWidget(validateButton);
  connect(validateButton, SIGNAL(clicked()), this, SLOT(validateTree()));
  validateButton->hide();
  
  barLayout->addSpacing(spacing);
  QPushButton* saveButton = new QPushButton(tr("Save\nXml File"), this);
  barLayout->addWidget(saveButton);
  connect(saveButton, SIGNAL(clicked()), this, SLOT(saveXml()));
  saveButton->hide();
  
  barLayout->addSpacing(spacing);
  QPushButton* refineButton = new QPushButton(tr("Cut\n"), this);
  barLayout->addWidget(refineButton);
  connect(refineButton, SIGNAL(clicked()), this, SLOT(cutGrid()));
  

  
  barLayout->addSpacing(5*spacing);
  QPushButton* doneButton = new QPushButton(tr("Done\n"), this);
  barLayout->addWidget(doneButton);
  connect(doneButton, SIGNAL(clicked()), this, SLOT(done()));

  barLayout->addStretch(10);
  flowbar->setLayout(barLayout);
  
  
  QToolBar* flowToolBar = new QToolBar;
  addToolBar(Qt::LeftToolBarArea,flowToolBar );
  flowToolBar->addWidget(flowbar);
}

void VCutWindow::cutGrid(){
  
  QTreeWidgetItem* root = tree->topLevelItem(0);
  QTreeWidgetItem* sphereItem = root->child(0)->child(0)->child(0);
  
  if(sphereItem->text(0)=="sphere"){
  
    x0 = sphereItem->child(0)->text(1).toDouble();
    y0 =  sphereItem->child(1)->text(1).toDouble();
    z0 =  sphereItem->child(2)->text(1).toDouble();
    r =  sphereItem->child(3)->text(1).toDouble();
  }
  
  
  outFilename  = inFilename.section('.', 0, -2)+"_cut.vog";
  
  
  QString tmpOutFilename = QFileDialog::getSaveFileName(this, tr("Result Vog File"),
                                             outFilename,
                                             tr("Volume Grid files (*.vog)"));
  if(tmpOutFilename.isEmpty())return;

  outFilename = tmpOutFilename;
  if(outFilename.section('.', -1, -1)!="vog")outFilename+=".vog";
  
 
  if(!(outFilename.section('/', -1, -1).section('.',0,0).isEmpty())){
    QString command = QString("vogcut");
    command += " -g " + inFilename + " -center " +  QString("%1 ").arg(x0)+
      QString(" %1 ").arg(y0) + QString(" %1 ").arg(z0)
      + " -r " + QString(" %1 ").arg(r);
      
      command += " -o " + outFilename;
      
 
      ProgressDialog* progress = new ProgressDialog(command, QString(),false);
      progress->show();
      connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus, QString)), this, SLOT(afterCut(QString, QProcess::ExitStatus, QString)));
    }
}
void VCutWindow::afterCut(QString command, QProcess::ExitStatus status, QString directory){
  if(status==QProcess::NormalExit){
    
    QStringList bnames;  
    viewer->load_boundary(outFilename, bnames);
  }else{
    qDebug()<<"vogcut failed";
  }
}

  void VCutWindow::loadGrid(){
  
  
  QString tmpInFilename =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                 QDir::currentPath(),
                                 tr("vog Files (*.vog)"));
  if(tmpInFilename.isEmpty())return;
  inFilename = tmpInFilename;
  QStringList bnames;
  viewer->load_boundary(inFilename, bnames);
  }
void VCutWindow::done(){

  close();
}
void VCutWindow::buildTree(){
 //build tree
  tree = new QTreeWidget;
  tree->setColumnCount(2);
  QStringList header;
  header << "name" << "value";
  tree->setHeaderLabels(header);
  root =  new QTreeWidgetItem(tree);
  root->setText(0, "region");
 
  
  tree->addTopLevelItem(root);
  tree->expandAll();
  tree->setSelectionMode(QAbstractItemView::SingleSelection);
  tree->viewport()->setAcceptDrops(true);
  tree->setDragEnabled(true);
  tree->setAcceptDrops(true);
  tree->setDropIndicatorShown(true);
  tree->setDragDropMode(QAbstractItemView::InternalMove);

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


   paraPages = new QStackedWidget;

   for(unsigned int i =0; i< defaultShapes.size(); i++){
     ParaPage *paraPage = new ParaPage(defaultShapes[i]);
     paraPages->addWidget(paraPage);
     //connect(paraPage, SIGNAL(valueChanged()), this, SIGNAL(valueChanged()));
     connect(paraPage, SIGNAL(valueChanged()), this, SLOT(updateShape()));     
   }

}


  

VCutWindow::VCutWindow( QWidget *parent):QMainWindow(parent){
  QWidget::setAttribute(Qt::WA_DeleteOnClose, true);

  //QScrollArea* centralScrollArea = new QScrollArea;
  tabWidget = new QTabWidget;
 
  QGroupBox* central = new QGroupBox;
  central->setFlat(true);
  QHBoxLayout* objLayout = new QHBoxLayout;

  
  buildTree();

 
  
 
  trans = new Transform( this);
  connect(trans, SIGNAL(tcChanged()), this, SLOT(updateTransform()));
  connect(tree, SIGNAL(itemClicked(QTreeWidgetItem*, int)),
          this, SLOT(showData(QTreeWidgetItem*)));

  tabWidget->addTab(paraPages, tr("Geometry"));
  tabWidget->addTab(trans, tr("Transform"));
  
  viewer = new GLViewer();
 //  viewerDock  = new QDockWidget("Graphics Viewer", this); 
//   viewerDock->setAllowedAreas(Qt::RightDockWidgetArea );
//   viewerDock->setWidget(viewer);
//   addDockWidget(Qt::RightDockWidgetArea, viewerDock);
  // viewerDock->setFloating(true);
  connect(this, SIGNAL(valueChanged(const QTreeWidgetItem*)), viewer, SLOT(updateDoc(const QTreeWidgetItem*))); 
  
   objLayout->addWidget(tabWidget);
   objLayout->addWidget(tree);
   tree->hide();
   objLayout->addWidget(viewer,2);
   
   central->setLayout(objLayout);
  
   
   
   createFlowBar();
   createToolBar();
   setCentralWidget(central);
    
   setWindowTitle(tr("VCutWindow"));
   setMinimumSize(1000, 700);
   addSphere();
   x0 = y0 = z0 = 0.0;
   r = 1;
}


      

void VCutWindow::showData(QTreeWidgetItem* item ){

  if(item->text(0)=="shape"){
    tabWidget->setCurrentWidget(paraPages);
    QStringList items;
    items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
      tr("leftplane") ;
  
    QString tp = item->child(0)->text(0);

     int index = items.indexOf(tp);
     index += 1;
     changePage(index);

     vector<double> para(item->child(0)->childCount());
     for(int i = 0; i< item->child(0)->childCount(); i++){
       para[i] = item->child(0)->child(i)->text(1).toDouble();
     }
     Shape* ashape = qobject_cast<ParaPage*>(paraPages->currentWidget())->shape;
     if(para.size() != ashape->para.size()) return;
     for(int i = 0; i< item->child(0)->childCount(); i++){
       ashape->para[i] = para[i];
     }
     qobject_cast<ParaPage*>(paraPages->currentWidget())->showValue();
  }else if(item->text(0)=="transform"){
    tabWidget->setCurrentWidget(trans);
    positions3d p[4];
    for(int i =0; i < 2; i++){
      p[i].x = item->child(i)->child(0)->text(1).toDouble();
      p[i].y = item->child(i)->child(1)->text(1).toDouble();
      p[i].z = item->child(i)->child(2)->text(1).toDouble();
    }
    p[2].x = item->child(2)->child(0)->text(1).toDouble();
    p[2].y = item->child(3)->child(0)->text(1).toDouble();
    p[2].z = item->child(4)->child(0)->text(1).toDouble();

    p[3].x = item->child(6)->child(0)->text(1).toDouble();
    p[3].y = item->child(6)->child(1)->text(1).toDouble();
    p[3].z = item->child(6)->child(2)->text(1).toDouble();
    
    trans->setValue(p);
  }

}

  
      

void VCutWindow::updateShape(){
  if(tree->currentItem() == 0 ||
     tree->currentItem()->text(0) != "shape")
    {     
      QMessageBox::warning(window(),
                           tr("update shape"),
                           tr("Please select a 'shape' node on the tree  first"));
      
      return;
    }
  QTreeWidgetItem* valueNode = tree->currentItem()->child(0);
  
  
  Shape* ashape = qobject_cast<ParaPage*>(paraPages->currentWidget())->shape;
  if(valueNode->childCount() !=(int) ashape->para.size())return;
  for(int i = 0; i < valueNode->childCount(); i++)valueNode->child(i)->setText(1,QString("%1").arg(ashape->para[i]));
  
  emit valueChanged(root);

}



  
void VCutWindow::changePage(int i){
  
  paraPages->setCurrentIndex(i);
  // emit valueChanged(root); ?
}


void VCutWindow::addTransform(){
  tabWidget->setCurrentWidget(trans);
  if(tree->currentItem() == 0 ||
     ( tree->currentItem()->text(0) != "object"
       && tree->currentItem()->text(0) != "transform"
       && tree->currentItem()->text(0) != "shape")){
    QMessageBox::warning(window(),
                         tr("add transform"),
                         tr("Please select an object on the tree  first"));
    
    return;
  }
  
  QTreeWidgetItem* transformNode =  new QTreeWidgetItem();
  transformNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );

  transformNode->setText(0, "transform");

 
  if( tree->currentItem()->text(0) =="object"){
    tree->currentItem()->insertChild(0,transformNode);
  }else if(tree->currentItem()->text(0) =="shape" ||
           tree->currentItem()->text(0) =="transform" ){
    tree->currentItem()->parent()->insertChild(0,transformNode);
  }  
  tree->setCurrentItem(transformNode);

 
  TranCoef tc = trans->value();
  
  //translate
  {
    QTreeWidgetItem* translateNode =  new QTreeWidgetItem(transformNode);
    translateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        translateNode->setText(0, "translate");
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "x0");
          x0Node->setText(1, QString("%1").arg(tc.translate.x));
        }
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "y0");
          x0Node->setText(1, QString("%1").arg(tc.translate.y));
        }
        {
          QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
          x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "z0");
          x0Node->setText(1, QString("%1").arg(tc.translate.z));
        }
    }
    
  {
    {     
      QTreeWidgetItem* translateNode =  new QTreeWidgetItem(transformNode);
      translateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      translateNode->setText(0, "translate");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x0");
        x0Node->setText(1, QString("%1").arg(tc.rotateCenter.x));
      }
      {      
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "y0");
          x0Node->setText(1, QString("%1").arg(tc.rotateCenter.y));
        }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "z0");
          x0Node->setText(1, QString("%1").arg(tc.rotateCenter.z));
        }
      }

    {
      QTreeWidgetItem* rotateNode =  new QTreeWidgetItem(transformNode);
      rotateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        rotateNode->setText(0, "rotateX");
        QTreeWidgetItem* thetaNode =  new QTreeWidgetItem(rotateNode);
        thetaNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        thetaNode->setText(0, "theta");
        thetaNode->setText(1,  QString("%1").arg(tc.rotateAngle.x));
    }
    {
        QTreeWidgetItem* rotateNode =  new QTreeWidgetItem(transformNode);
        rotateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        rotateNode->setText(0, "rotateY");
        QTreeWidgetItem* thetaNode =  new QTreeWidgetItem(rotateNode);
        thetaNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        thetaNode->setText(0, "theta");
        thetaNode->setText(1,  QString("%1").arg(tc.rotateAngle.y));
      }
    {
        QTreeWidgetItem* rotateNode =  new QTreeWidgetItem(transformNode);
        rotateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        rotateNode->setText(0, "rotateZ");
        QTreeWidgetItem* thetaNode =  new QTreeWidgetItem(rotateNode);
        thetaNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        thetaNode->setText(0, "theta");
        thetaNode->setText(1,  QString("%1").arg(tc.rotateAngle.z));
      } 
      
      
      

    {
        
      QTreeWidgetItem* translateNode =  new QTreeWidgetItem(transformNode);
      translateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      translateNode->setText(0, "translate");
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "x0");
        x0Node->setText(1, QString("%1").arg(-tc.rotateCenter.x));
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "y0");
        x0Node->setText(1, QString("%1").arg(-tc.rotateCenter.y));
      }
      {
        QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
        x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
        x0Node->setText(0, "z0");
        x0Node->setText(1, QString("%1").arg(-tc.rotateCenter.z));
      }
    }
    } 
    
  {
    QTreeWidgetItem* translateNode =  new QTreeWidgetItem(transformNode);
    translateNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    translateNode->setText(0, "scale");
    {
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "x0");
      x0Node->setText(1, QString("%1").arg(tc.scale.x));
           }
    {
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "y0");
             x0Node->setText(1, QString("%1").arg(tc.scale.y));
        }
    {
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(translateNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
          x0Node->setText(0, "z0");
          x0Node->setText(1, QString("%1").arg(tc.scale.z));
    }
  }
    

  trans->clear();
  
  tree->expandItem(transformNode);
  tree->resizeColumnToContents(0);
  emit valueChanged(root);
}
void VCutWindow::updateTransform(){

  if(tree->currentItem() == 0 ||
     tree->currentItem()->text(0) != "transform")
    {     
      QMessageBox::warning(window(),
                           tr("update transform"),
                           tr("Please select an transform on the tree  first"));
      
      return;
    }
  
  QTreeWidgetItem* transformNode = tree->currentItem();
  QTreeWidgetItem* translateNode, x0Node;
  
  TranCoef tc = trans->value();
  
  
 //translate
  {
        translateNode =  transformNode->child(0);
        translateNode->child(0)->setText(1, QString("%1").arg(tc.translate.x));
        translateNode->child(1)->setText(1, QString("%1").arg(tc.translate.y));
        translateNode->child(2)->setText(1, QString("%1").arg(tc.translate.z));
        
  }
    
  
  {     
    translateNode =  transformNode->child(1);
    translateNode->child(0)->setText(1, QString("%1").arg(tc.rotateCenter.x));
    translateNode->child(1)->setText(1, QString("%1").arg(tc.rotateCenter.y));
    translateNode->child(2)->setText(1, QString("%1").arg(tc.rotateCenter.z));
  }
  {
    translateNode =  transformNode->child(2);
    translateNode->child(0)->setText(1, QString("%1").arg(tc.rotateAngle.x));
  }
  {
    translateNode =  transformNode->child(3);
    translateNode->child(0)->setText(1, QString("%1").arg(tc.rotateAngle.y));
  }
  {
    translateNode =  transformNode->child(4);
    translateNode->child(0)->setText(1, QString("%1").arg(tc.rotateAngle.z));
  }
   {     
    translateNode =  transformNode->child(5);
    translateNode->child(0)->setText(1, QString("%1").arg(-tc.rotateCenter.x));
    translateNode->child(1)->setText(1, QString("%1").arg(-tc.rotateCenter.y));
    translateNode->child(2)->setText(1, QString("%1").arg(-tc.rotateCenter.z));
  }      
      
  {     
    translateNode =  transformNode->child(6);
    translateNode->child(0)->setText(1, QString("%1").arg(tc.scale.x));
    translateNode->child(1)->setText(1, QString("%1").arg(tc.scale.y));
    translateNode->child(2)->setText(1, QString("%1").arg(tc.scale.z));
   }      
   
   tree->expandItem(transformNode);
   tree->resizeColumnToContents(0);
   emit valueChanged(root);
}



void VCutWindow::addOp(){
  if(tree->currentItem()==0
     ||(tree->currentItem()->text(0) !="object"
        && tree->currentItem()->text(0) !="region"
        && tree->currentItem()->text(1) !="complement")){
    QMessageBox::warning(window(),
                         tr("add operator"),
                         tr("Please select an 'object' node or a 'region' node first"));
    return;
  }
  if(tree->currentItem() == root){
     QMessageBox::warning(window(),
                         tr("add operator"),
                         tr("operator is not allowed to insert in front of root of the tree"));
     return;
  }
  
  QStringList items;
  items << tr("union") << tr("intersection") << tr("difference") << tr("complement");
  
  bool ok;
  QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                       tr("perators:"), items, 0, false, &ok);
  if (ok && !item.isEmpty()){
    if(tree->currentItem()!=0 &&
       (tree->currentItem()->text(0) =="object"||
        tree->currentItem()->text(0) == "region" ||
        tree->currentItem()->text(1) == "complement")){
      
      
      QTreeWidgetItem*  p=tree->currentItem()->parent();
      QTreeWidgetItem* newItem = new QTreeWidgetItem();
       newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled|Qt::ItemIsDragEnabled );

      newItem->setText(0, "op");
      newItem->setText(1, item);
      int index = p->indexOfChild(tree->currentItem());
      p->insertChild(index, newItem);
    }else{
      QMessageBox::warning(window(),
                           tr("add operator"),
                           tr("Please select an object on the tree  first"));
    }
  }
  

}

void VCutWindow::addSphere(){
  tabWidget->setCurrentWidget(paraPages);
  QString item = "sphere";
  
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
    root->addChild(objItem);
  }
  tree->setCurrentItem(newItem);
  int index = 0;
  index += 1;
  changePage(index);
  Shape* ashape = qobject_cast<ParaPage*>(paraPages->currentWidget())->shape;
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
  default:
    break;
  }
 emit valueChanged(root);
}
  

void VCutWindow::addShape(){
  tabWidget->setCurrentWidget(paraPages);
  QStringList items;
  items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
    tr("leftplane");

     bool ok;
     QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
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
         root->addChild(objItem);
       }
       tree->setCurrentItem(newItem);
       int index = items.indexOf(item);
       index += 1;
       changePage(index);
       Shape* ashape = qobject_cast<ParaPage*>(paraPages->currentWidget())->shape;
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
     emit valueChanged(root);
}



void VCutWindow::addRegion(){

  if(tree->currentItem()!=0 &&
     tree->currentItem()->text(0) =="region"){
    QTreeWidgetItem* newItem = new QTreeWidgetItem(tree->currentItem());
    newItem->setText(0, "region");
  }else{
    QTreeWidgetItem* newItem = new QTreeWidgetItem();
    newItem->setText(0, "region");
    root->addChild(newItem);
  }
  tree->expandAll();
  tree->resizeColumnToContents(0);
}


void VCutWindow::removeNode(){
  if(tree->currentItem()!=0 &&
     tree->currentItem()!= root){
    QTreeWidgetItem* selectedItem =  tree->currentItem();
    if(selectedItem->text(0)=="object"
       || selectedItem->text(0)=="op"
       || selectedItem->text(0)=="region"
       || selectedItem->text(0)=="transform"){
      QTreeWidgetItem* parent = selectedItem->parent();  
      parent->removeChild(selectedItem);
      delete selectedItem;
    }else{
      QMessageBox::warning(this, tr("Application"),
                           tr("remove this node may cause illegal tree structure"));
                          
    }   
    
  }
  emit valueChanged(root);
}
   

  
bool VCutWindow::saveXml(){
  QDomDocument doc = tree2dom(getRoot());
  if(doc.isNull())return false;
  QString fileName = inFilename.left(inFilename.lastIndexOf('.'))+".xml";

    fileName = QFileDialog::getSaveFileName(this, tr("Save .xml File"),
                                                 fileName,
                                                  tr("xml Files (*.xml)"));
    if(fileName.isEmpty())return false;
    
    if(fileName.section('.', -1, -1) != "xml") fileName = fileName + ".xml";
 
 
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Application"),
                           tr("Cannot write file %1:\n%2.")
                           .arg(fileName)
                           .arg(file.errorString()));
      return false;
  }
    
  QTextStream out(&file);

  doc.save(out, 2, QDomNode::EncodingFromDocument);
  
  //viewer->markVolumeNodes(filename);
  //viewer->markNodes();
  return true;
}


VCutWindow::~VCutWindow(){
   for(unsigned int i = 0 ; i< defaultShapes.size(); i++){
    if(defaultShapes[i]){
      delete defaultShapes[i];
      defaultShapes[i] = 0;
    }
  }
  defaultShapes.clear();
 
}



bool VCutWindow::validateRegion( QTreeWidgetItem* item){
  bool result = true;
  //root of the tree is 'region'
  if(item->text(0)!="region"){
    tree->setCurrentItem(item);
    int button = QMessageBox::question(this, tr("Error in region:"),
                                       tr("should be 'region'\n")+
                                       tr("Do you want to continue?"),
                                       QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok);
    result = false;
    if(button == QMessageBox::No)return false;
  }
  int numChild = item->childCount();

  //the children of region can only be 'op', 'object' or 'region'
  for(int i = 0;i < numChild; i++){
    QTreeWidgetItem* theChild = item->child(i); 
    if(theChild->text(0)!="region" &&
       theChild->text(0)!="object" &&
       theChild->text(0)!="op"){
      tree->setCurrentItem(item);
      int button = QMessageBox::question(this, tr("Error in region:"),
                                         tr("children of region can only be op, object or region.\n")+
                                         tr("Illegal child of region: ")+ theChild->text(0)+
                                         tr("\nDo you want to continue?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
      result = false;
      if(button == QMessageBox::No)return false;
    }
  }
    
  //between children
  for(int i = 0;i < numChild; i++){
    QTreeWidgetItem* theChild = item->child(i); 
    
    //'region' and 'object' should be followed by binary operator 
      if(theChild->text(0)=="region" ||
         theChild->text(0)=="object"){
        if(i< (numChild-1)){
          QTreeWidgetItem* nextChild = item->child(i+1);
          if(nextChild->text(0) != "op"){
       
            tree->setCurrentItem(theChild);
            int button = QMessageBox::question(this, tr("Error in region:"),
                                               theChild->text(0)+tr(" should be followed by op\n")+
                                               tr("Do you want to continue?"),
                                               QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
            result = false;
            if(button == QMessageBox::No)return false;
          }
          if(nextChild->text(0) == "op" && nextChild->text(1) =="complement"){
       
            tree->setCurrentItem(theChild);
            int button = QMessageBox::question(this, tr("Error in region:"),
                                               theChild->text(0)+tr(" should be followed by binary operator\n")+
                                               tr("Do you want to continue?"),
                                               QMessageBox::Ok|QMessageBox::No, QMessageBox::No); 
            result = false;
            if(button == QMessageBox::No)return false;
          }
        }
      }else if(theChild->text(0)=="op" && theChild->text(1)=="complement"){// unary operator must be followed by region or object
        if(i< (numChild-1)){
          QTreeWidgetItem* nextChild = item->child(i+1);
          if(nextChild->text(0) != "region" && nextChild->text(0) != "object"){
            tree->setCurrentItem(theChild);
            int button = QMessageBox::question(this, tr("Error in region:"),
                                               theChild->text(0)+tr(" should be followed by object or region\n")+
                                               tr("Do you want to continue?"),
                                               QMessageBox::Ok|QMessageBox::No, QMessageBox::No); 
            result = false;
            if(button == QMessageBox::No)return false;
          }
        }else{
          tree->setCurrentItem(theChild);
          int button = QMessageBox::question(this, tr("Error in region:"),
                                             theChild->text(0)+tr(" should not be the last child\n")+
                                             tr("Do you want to continue?"),
                                             QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
          result = false;
          if(button == QMessageBox::No)return false;
        }
      }else{ //binary operator
        if(i==0){//binary operator can not be the first child
           tree->setCurrentItem(theChild);
           int button = QMessageBox::question(this, tr("Error in region:"),
                                              theChild->text(0)+tr(" should not be the first child\n")+
                                             tr("Do you want to continue?"),
                                              QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
           result = false;
          if(button == QMessageBox::No)return false;
        }else if(i == (numChild -1)){//binary operator can not be the last child
          tree->setCurrentItem(theChild);
          int button = QMessageBox::question(this, tr("Error in region:"),
                                             theChild->text(0)+tr(" should not be the last child\n")+
                                             tr("Do you want to continue?"),
                                             QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
          result = false;
          if(button == QMessageBox::No)return false;
          
        }else{//binary should between two children
          QTreeWidgetItem* previousChild = item->child(i-1); //previoud should be region or object
          if(previousChild->text(0) != "region" && previousChild->text(0) != "object"){
            tree->setCurrentItem(theChild);
            int button = QMessageBox::question(this, tr("Error in region:"),
                                               theChild->text(0)+tr(" should follow child object or region\n")+
                                               tr("Do you want to continue?"),
                                               QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
            result = false;
            if(button == QMessageBox::No)return false; 
          }
          
          QTreeWidgetItem* nextChild = item->child(i+1); //next should be region, object or unary operator
          if(nextChild->text(0) != "region"
             && nextChild->text(0) != "object"
             && !(nextChild->text(0) == "op" && nextChild->text(1)=="complement")){
            tree->setCurrentItem(theChild);
            int button = QMessageBox::question(this, tr("Error in region:"),
                                               theChild->text(0)+tr(" should be followed by object or region or op complement\n")+
                                               tr("Do you want to continue?"),
                                               QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
            result = false;
            if(button == QMessageBox::No)return false;
          }
        }
      }
      
    }//finish between children

    //validate each child
  for(int i = 0;i < numChild; i++){
    QTreeWidgetItem* theChild = item->child(i); 
    if(theChild->text(0)=="region")result = (validateRegion(theChild)) && result;
    else if(theChild->text(0)=="object")result = (validateObject(theChild)) && result;
    //else if(theChild->text(0)=="op")validateOp(theChild);
  }
  
  return result;
  //QMessageBox::information(this, "validation", tr("This region passed the validation"));
  
}

bool VCutWindow::validateObject(QTreeWidgetItem* item){
  bool result = true;
  tree->setCurrentItem(item);
  int numChild = item->childCount();
  //the children of object can only be 'transform' or 'shape'
  for(int i = 0;i < numChild; i++){
    QTreeWidgetItem* theChild = item->child(i); 
    if(theChild->text(0)!="transform" &&
       theChild->text(0)!="shape"){
      tree->setCurrentItem(theChild);
      int button = QMessageBox::question(this, tr("Error:"),
                                         tr("children of object can only be transform or shape.\n")+
                                         tr("Illegal child of object: ")+ theChild->text(0)+
                                         tr("\nDo you want to continue?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
      result = false;
      if(button == QMessageBox::No)return false;
    }
  }
  //between children
  int numShape = 0;
  for(int i = 0;i < numChild; i++){
    QTreeWidgetItem* theChild = item->child(i); 
    //transform should not be the last child
    if(theChild->text(0)=="transform" && i == (numChild-1)){
      tree->setCurrentItem(theChild);
      int button = QMessageBox::question(this, tr("Error:"),
                                         tr("transform should not be the last child\n")+
                                         tr("\nDo you want to continue?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
      result = false;
      if(button == QMessageBox::No)return false;
    }
    if(theChild->text(0)=="shape"){
      numShape++;
    }
  }
  if(numShape !=1){
    tree->setCurrentItem(item);
    int button = QMessageBox::question(this, tr("Error:"),
                                       tr("each object should have one and only one shape\n")+
                                       tr("\nDo you want to continue?"),
                                       QMessageBox::Ok|QMessageBox::No, QMessageBox::No);
    result = false;
    if(button == QMessageBox::No)return false;
  }
  // validate each child
  // for(int i = 0;i < numChild; i++){
  //     QTreeWidgetItem* theChild = item->child(i); 
  //     if(theChild->text(0)=="transform"){
  //validateTransform(theChild);
  //     }else if(theChild->text(0)=="shape")validateShape(theChild);
  //   }
  
  return result;
}

    
        
void  VCutWindow::validateTree(){
  if( validateRegion(tree->topLevelItem(0))){
    QMessageBox::information(this, "validation", tr("This tree passed the validation"));
  }
  
}
 
void VCutWindow::helpClicked(){
  HelpWindow* helpwindow = new HelpWindow("page_fvmadapt.html");
  helpwindow->show();
}
QTreeWidgetItem* VCutWindow::getRoot(){
  return root;
}

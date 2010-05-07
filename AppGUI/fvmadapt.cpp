/////////////////////////////////////////////////
//  Filename: fvmadapt.cpp
//
//  Contains: Implementation of Fvmadapt class
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
#include "fvmadapt.h"
#include "grid.h"
#include "pages.h"


void Shape::reset(){
  switch(tp){
  case  SPHERE:
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 1;
    
    break;
  
  case CONE:

   
   
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 1;
    para[4] = 0;
    para[5] = 1;
    break;
   
  case CYLINDER:
    para[0] = 0;
    para[1] = 0;
    para[2] = 1;
    para[3] = 0;
    para[4] = 1;
    break;
    
  case BOX:
   
     para[0] = 0;
     para[1] = 0;
     para[2] = 0;
     para[3] = 1;
     para[4] = 1;
     para[5] = 1;
     break;
  default:
    para[0] = 0;
  }
}
 
ParaPage::ParaPage(Shape* s, QWidget* parent):QGroupBox(tr(""),parent),shape(s){
  if(shape==0)return;
  
  signalMapper = new QSignalMapper(this);
  QBoxLayout *mainLayout;
  
  switch(shape->tp){
  case SPHERE:
    {

      setTitle("parameters of sphere:");
      objs.resize(4);

      
      objs[0] = new FloatSlider(tr("x0"));
      objs[1] = new FloatSlider(tr("y0"));
      objs[2] = new FloatSlider(tr("z0"));
      objs[3] = new FloatSlider(tr("r"));
      for(int i = 0; i < 4; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));


      QGroupBox *centerGroup = new QGroupBox("center", this);
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      hBoxLayout->addWidget(objs[1]);
      hBoxLayout->addWidget(objs[2]);
      centerGroup->setLayout(hBoxLayout);

      QGroupBox *radiusGroup = new QGroupBox("radius");
      QVBoxLayout *hBoxLayout2 = new QVBoxLayout;
      hBoxLayout2->addWidget(objs[3]);
      radiusGroup->setLayout(hBoxLayout2);

      
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(centerGroup);
      mainLayout->addWidget(radiusGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
    }
    break;
  case CONE:
     {

      setTitle("parameters of cone:");
      objs.resize(6);

      
      objs[0] = new FloatSlider(tr("x0"));
      objs[1] = new FloatSlider(tr("y0"));
      objs[2] = new FloatSlider(tr("z0"));
      objs[3] = new FloatSlider(tr("r"));
      objs[4] = new FloatSlider(tr("z1: (z1>z0)"));
      objs[5] = new FloatSlider(tr("z2: (z2>z1)"));
      
      
      for(int i = 0; i < 6; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *centerGroup = new QGroupBox("position of vertex:");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      hBoxLayout->addWidget(objs[1]);
      hBoxLayout->addWidget(objs[2]);
      centerGroup->setLayout(hBoxLayout);

      QGroupBox *radiusGroup = new QGroupBox("radius of base/ distance from center of base to vertex)");
      QVBoxLayout *hBoxLayout2 = new QVBoxLayout;
      hBoxLayout2->addWidget(objs[3]);
      radiusGroup->setLayout(hBoxLayout2);
      
      QGroupBox *z1Group = new QGroupBox("z value of base");
      QVBoxLayout *hBoxLayout3 = new QVBoxLayout;
      hBoxLayout3->addWidget(objs[4]);
      z1Group->setLayout(hBoxLayout3);
      
      QGroupBox *z2Group = new QGroupBox("z value of top");
      QVBoxLayout *hBoxLayout4 = new QVBoxLayout;
      hBoxLayout4->addWidget(objs[5]);
      z2Group->setLayout(hBoxLayout4);
      
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(centerGroup);
      mainLayout->addWidget(radiusGroup);
      mainLayout->addWidget(z1Group);
      mainLayout->addWidget(z2Group);
      mainLayout->addStretch(2);
      setLayout(mainLayout);
     }
     break;
  case CYLINDER:
    {
      setTitle("parameters of cylinder:");
      objs.resize(5);
      
      
      objs[0] = new FloatSlider(tr("x0"));
      objs[1] = new FloatSlider(tr("y0"));
      objs[2] = new FloatSlider(tr("r"));
      objs[3] = new FloatSlider(tr("z1"));
      objs[4] = new FloatSlider(tr("z2"));
    
      
      for(int i = 0; i < 5; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *centerGroup = new QGroupBox("center of top");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      hBoxLayout->addWidget(objs[1]);
      centerGroup->setLayout(hBoxLayout);

      QGroupBox *radiusGroup = new QGroupBox("radius");
      QVBoxLayout *hBoxLayout2 = new QVBoxLayout;
      hBoxLayout2->addWidget(objs[2]);
      radiusGroup->setLayout(hBoxLayout2);

      QGroupBox *z1Group = new QGroupBox("z value of the base");
      QVBoxLayout *hBoxLayout3 = new QVBoxLayout;
      hBoxLayout3->addWidget(objs[3]);
      z1Group->setLayout(hBoxLayout3);
      
       QGroupBox *z2Group = new QGroupBox("z value of the top");
       QVBoxLayout *hBoxLayout4 = new QVBoxLayout;
       hBoxLayout4->addWidget(objs[4]);
       z2Group->setLayout(hBoxLayout4);
      
       mainLayout = new QVBoxLayout;
       mainLayout->addWidget(centerGroup);
       mainLayout->addWidget(radiusGroup);
       mainLayout->addWidget(z1Group);
       mainLayout->addWidget(z2Group);
         mainLayout->addStretch(2);
       setLayout(mainLayout);
    }
    break;
  case BOX:

    {
      setTitle("parameters of box:");
      objs.resize(6);
      
      
      objs[0] = new FloatSlider(tr("x1"));
      objs[1] = new FloatSlider(tr("y1"));
      objs[2] = new FloatSlider(tr("z1"));
      objs[3] = new FloatSlider(tr("x2"));
      objs[4] = new FloatSlider(tr("y2"));
      objs[5] = new FloatSlider(tr("z2"));
      
      for(int i = 0; i < 6; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("vertex with minimum x,y, z values");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      hBoxLayout->addWidget(objs[1]);
      hBoxLayout->addWidget(objs[2]);
      minGroup->setLayout(hBoxLayout);

       QGroupBox *maxGroup = new QGroupBox("vertex with maximum x,y, z values");
       QVBoxLayout *hBoxLayout2 = new QVBoxLayout;
      hBoxLayout2->addWidget(objs[3]);
      hBoxLayout2->addWidget(objs[4]);
      hBoxLayout2->addWidget(objs[5]);
      maxGroup->setLayout(hBoxLayout2);
      
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
      mainLayout->addWidget(maxGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
    }
    
    break;
  case PXPLANE:
 {
   setTitle("parameters of x+plane:");
   objs.resize(1);
      
      
      objs[0] = new FloatSlider(tr("x1"));
     
      
      for(int i = 0; i < 1; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("position of the plane");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      minGroup->setLayout(hBoxLayout);

            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
 }
    
 break;
 

 case NXPLANE:
 {
   setTitle("parameters of x-plane:");
   objs.resize(1);
      
      
      objs[0] = new FloatSlider(tr("x1"));
     
      
      for(int i = 0; i < 1; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("position of the plane");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      minGroup->setLayout(hBoxLayout);

            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
 }
    
 break;
case PYPLANE:
 {
   setTitle("parameters of y+plane:");
   objs.resize(1);
      
      
      objs[0] = new FloatSlider(tr("y1"));
     
      
      for(int i = 0; i < 1; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("position of the plane");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      minGroup->setLayout(hBoxLayout);

            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
 }
    
 break;
    case NYPLANE:
 {
   setTitle("parameters of y-plane:");
   objs.resize(1);
      
      
      objs[0] = new FloatSlider(tr("y1"));
     
      
      for(int i = 0; i < 1; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("position of the plane");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      minGroup->setLayout(hBoxLayout);

            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
      mainLayout->addStretch(2);
      setLayout(mainLayout);
 }
    
 break;
 case PZPLANE:
 {
   setTitle("parameters of z+plane:");
   objs.resize(1);
      
      
      objs[0] = new FloatSlider(tr("z1"));
     
      
      for(int i = 0; i < 1; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));

      QGroupBox *minGroup = new QGroupBox("position of the plane");
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      minGroup->setLayout(hBoxLayout);

            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(minGroup);
      mainLayout->addStretch(2);
      setLayout(mainLayout);
 }
    
 break;
 case NZPLANE:
 {
   setTitle("parameters of z-plane:");
   objs.resize(1);
   
   
   objs[0] = new FloatSlider(tr("z1"));
   
   
   for(int i = 0; i < 1; i++){
     objs[i]->setValue(shape->para[i]);
     signalMapper->setMapping(objs[i], i);
     connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
   }
   connect(signalMapper, SIGNAL(mapped(int)),
           this, SLOT(setValue(int)));

   QGroupBox *minGroup = new QGroupBox("position of the plane");
   QVBoxLayout *hBoxLayout = new QVBoxLayout;
   hBoxLayout->addWidget(objs[0]);
   minGroup->setLayout(hBoxLayout);

   
   mainLayout = new QVBoxLayout;
   mainLayout->addWidget(minGroup);
   mainLayout->addStretch(2);
   setLayout(mainLayout);
 }
 
 break;
  default:
    break;
  }
  
}
void ParaPage::setValue(int i){
  if(i < 0 || i >= (int)objs.size() || i >= (int)(shape->para.size()))return;
  shape->para[i] = objs[i]->value();
  emit valueChanged();
}
void ParaPage::showValue(){
  for(unsigned int i = 0; i < objs.size(); i++){
    objs[i]->setValue(shape->para[i]);
  }
}
void ParaPage::reset(){
  shape->reset();
  showValue();
}

//////////////////////////////////////////////////////////
//  
//////////////////////////////////////////////////////////

Transform::Transform( QWidget *parent)
  :  QGroupBox(tr("transform:"),parent) 
{
  
 
  translate = new VectSlider(tr("translate"));
  rotateCenter  = new VectSlider(tr("rotateCenter"));
  rotateAngle = new VectSlider(tr("rotateAngle"));
  scale = new  VectSlider(tr("scale"));
  scale->setValue(positions3d(1.0, 1.0, 1.0));
  translate->setRange(-1e5, 1e5);
  rotateCenter->setRange(-1e5, 1e5);
  scale->setRange(-1e5, 1e5);
  rotateAngle->setRange(-180, 180);
  
   
  connect(translate, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(rotateCenter, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(rotateAngle, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(scale, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));

  
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addWidget(translate);
  mainLayout->addWidget(rotateCenter);
  mainLayout->addWidget(rotateAngle);
  mainLayout->addWidget(scale);

  setLayout(mainLayout);
}

void Transform::clear(){
 
 
  translate->setValue(positions3d(0.0, 0.0, 0.0));
  rotateCenter->setValue(positions3d(0.0, 0.0, 0.0));
  rotateAngle->setValue(positions3d(0.0, 0.0, 0.0));
  scale->setValue(positions3d(1.0, 1.0, 1.0));
}





TranCoef Transform::value(){
  TranCoef f;
  f.translate = translate->value();
  f.rotateCenter = rotateCenter->value();
  f.rotateAngle = rotateAngle->value();
  f.scale = scale->value();
  return f;
}

void Transform::setValue(positions3d* p){
  translate->setValue(p[0]);
  rotateCenter->setValue(p[1]);
  rotateAngle->setValue(p[2]);
  scale->setValue(p[3]);
}



void Transform::setInfo(){
 
 //  tc.rotateAngle = rotateAngle->value();
//   tc.rotateCenter = rotateCenter->value();
//   tc.translate = translate->value();
//   tc.scale = scale->value();

  //  IDMatrix2 idM = IDMatrix2(gridId, currentM()); 
  emit tcChanged();
  
}



void FVMAdapt::createToolBar(){
  int spacing =2;  
 
  toolbar = new QGroupBox("");
  toolbar->setFlat(true);
  QHBoxLayout* barLayout = new QHBoxLayout;
   barLayout->addStretch(10);
  QPushButton* shapeButton = new QPushButton(tr("Show Shapes"), this);
  barLayout->addWidget(shapeButton);
  connect(shapeButton, SIGNAL(clicked()), this, SIGNAL(showShapesClicked()));
  
  barLayout->addSpacing(spacing);
  QPushButton* nodesButton = new QPushButton(tr("Show Marked Nodes"), this);
  barLayout->addWidget(nodesButton);
  connect(nodesButton, SIGNAL(clicked()), this, SIGNAL(showNodesClicked()));
  
  //barLayout->addSpacing(spacing);
 
  
  toolbar->setLayout(barLayout);
}

void FVMAdapt::createTreeBar(){
  int spacing =2;  
 
  treebar = new QGroupBox("build a  tree");
  // flowbar->setFlat(true);
  QHBoxLayout* barLayout = new QHBoxLayout;

  QPushButton* addShapeButton = new QPushButton(tr("Add Shape"), this);
  barLayout->addWidget(addShapeButton);
  connect(addShapeButton, SIGNAL(clicked()), this, SLOT(addShape()));

  barLayout->addSpacing(spacing);
  QPushButton* addTransButton = new QPushButton(tr("Add Transform"), this);
  barLayout->addWidget(addTransButton);
  connect(addTransButton, SIGNAL(clicked()), this, SLOT(addTransform()));
  
  barLayout->addSpacing(spacing);
  QPushButton* addOpButton = new QPushButton(tr("Add Operator"), this);
  barLayout->addWidget(addOpButton);
  connect(addOpButton, SIGNAL(clicked()), this, SLOT(addOp()));
  
  

  barLayout->addSpacing(spacing);
  QPushButton* removeNodeButton = new QPushButton(tr("Remove Node"), this);
  barLayout->addWidget(removeNodeButton);
  connect(removeNodeButton, SIGNAL(clicked()), this, SLOT(removeNode()));
  
 


  barLayout->addSpacing(spacing);
  QPushButton* addRegionButton = new QPushButton(tr("Add Region "), this);
  barLayout->addWidget(addRegionButton);
  connect(addRegionButton, SIGNAL(clicked()), this, SLOT(addRegion()));
  
  barLayout->addSpacing(spacing);
  QPushButton* helpButton = new QPushButton(tr("Help"), this);
  barLayout->addWidget(helpButton);
  connect(helpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  barLayout->addStretch(10);
  treebar->setLayout(barLayout);
   
}
  


  
void FVMAdapt::createFlowBar(){
  int spacing =2;  
 
  flowbar = new QGroupBox("save xml file and refine grids");
   QHBoxLayout* barLayout = new QHBoxLayout;

  barLayout->addSpacing(spacing);
  QPushButton* saveButton = new QPushButton(tr("Save Xml File"), this);
  barLayout->addWidget(saveButton);
  connect(saveButton, SIGNAL(clicked()), this, SLOT(saveXml()));
  
  barLayout->addSpacing(spacing);
  QPushButton* refineButton = new QPushButton(tr("Refine Grid"), this);
  barLayout->addWidget(refineButton);
  connect(refineButton, SIGNAL(clicked()), this, SIGNAL(refineGrids()));


  
  barLayout->addSpacing(spacing);
  QPushButton* doneButton = new QPushButton(tr("Done"), this);
  barLayout->addWidget(doneButton);
  connect(doneButton, SIGNAL(clicked()), this, SLOT(done()));

  barLayout->addStretch(10);
  flowbar->setLayout(barLayout);
   
}
  
void FVMAdapt::done(){

  close();
}


FVMAdapt::FVMAdapt(QString fileName, QWidget *parent):QWidget(parent),filename(fileName){
  QWidget::setAttribute(Qt::WA_DeleteOnClose, true);

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
    vector<double> vdefault(6);
    vdefault[0] = 0;
    vdefault[1] = 0;
    vdefault[2] = 0;
    vdefault[3] = 1;
    vdefault[4] = 0;
    vdefault[5] = 1;
    Shape* shape = new Shape(CONE, vdefault);
    defaultShapes.push_back(shape);
   }
  //CYLINDER
 {
    vector<double> vdefault(5);
    vdefault[0] = 0;
    vdefault[1] = 0;
    vdefault[2] = 1;
    vdefault[3] = 0;
    vdefault[4] = 1;
   
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
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(PXPLANE, vdefault);
     defaultShapes.push_back(shape);
   }
   {
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(NXPLANE, vdefault);
     defaultShapes.push_back(shape);
   }
   {
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(PYPLANE, vdefault);
     defaultShapes.push_back(shape);
   }
   {
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(NYPLANE, vdefault);
     defaultShapes.push_back(shape);
   }
   {
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(PZPLANE, vdefault);
     defaultShapes.push_back(shape);
   }
   {
     vector<double> vdefault(1);
     vdefault[0] = 0;
     Shape* shape = new Shape(NZPLANE, vdefault);
     defaultShapes.push_back(shape);
   }

   paraPages = new QStackedWidget;

   for(unsigned int i =0; i< defaultShapes.size(); i++){
     ParaPage *paraPage = new ParaPage(defaultShapes[i]);
     paraPages->addWidget(paraPage);
     connect(paraPage, SIGNAL(valueChanged()), this, SIGNAL(valueChanged()));
     connect(paraPage, SIGNAL(valueChanged()), this, SLOT(updateShape()));     
   }

  

   trans = new Transform( this);
   connect(trans, SIGNAL(tcChanged()), this, SIGNAL(valueChanged()));
   connect(trans, SIGNAL(tcChanged()), this, SLOT(updateTransform()));

   connect(tree, SIGNAL(itemClicked(QTreeWidgetItem*, int)),
           this, SLOT(showData(QTreeWidgetItem*)));
   
   QHBoxLayout* objLayout = new QHBoxLayout; 
   objLayout->addWidget(paraPages);
   objLayout->addWidget(trans);
   objLayout->addWidget(tree,2);


   QHBoxLayout* barLayout = new QHBoxLayout;
   //createToolBar();
   createFlowBar();
   createTreeBar();
   
   barLayout->addWidget(treebar);
   barLayout->addWidget(flowbar);
   
     QVBoxLayout* mainLayout = new QVBoxLayout;
     
   mainLayout->addLayout(barLayout);
   mainLayout->addLayout(objLayout);
   setLayout(mainLayout);
  
}

void FVMAdapt::showData(QTreeWidgetItem* item ){

  if(item->text(0)=="shape"){
    QStringList items;
    items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
      tr("x_plus_plane") << tr("x_minus_plane")<<tr("y_plus_plane")<<tr("y_minus_plane")<<tr("z_plus_plane")<<tr("z_minus_plane");
  
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

  
      

void FVMAdapt::updateShape(){
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
  


}



  
void FVMAdapt::changePage(int i){
  
  paraPages->setCurrentIndex(i);
  emit valueChanged();
}


void FVMAdapt::addTransform(){

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
}
void FVMAdapt::updateTransform(){

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
}



void FVMAdapt::addOp(){
  if(tree->currentItem()==0
     ||(tree->currentItem()->text(0) !="object"
        && tree->currentItem()->text(0) !="region")){
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
        tree->currentItem()->text(0) == "region")){
      
      
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

void FVMAdapt::addShape(){
  QStringList items;
  items << tr("sphere") << tr("cone") << tr("cylinder") << tr("box")<<
    tr("x_plus_plane") << tr("x_minus_plane")<<tr("y_plus_plane")<<tr("y_minus_plane")<<tr("z_plus_plane")<<tr("z_minus_plane");

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
           
           QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
           z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           z1Node->setText(0, "z1");
           z1Node->setText(1, QString("%1").arg(ashape->para[4]));
      
           QTreeWidgetItem* z2Node = new QTreeWidgetItem(valueNode);
           z2Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           z2Node->setText(0, "z2");
           z2Node->setText(1, QString("%1").arg(ashape->para[5]));
           

      
         }





      
         break;
         
       case CYLINDER:
         {
           valueNode->setText(0, "cylinder");
           QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
           x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           x0Node->setText(0, "x0");
           x0Node->setText(1, QString("%1").arg(ashape->para[0]));
           
           QTreeWidgetItem* y0Node = new QTreeWidgetItem(valueNode);
           y0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           y0Node->setText(0, "y0");
           y0Node->setText(1, QString("%1").arg(ashape->para[1]));
           
           
           QTreeWidgetItem* rNode = new QTreeWidgetItem(valueNode);
           rNode->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           rNode->setText(0, "r");
           rNode->setText(1, QString("%1").arg(ashape->para[2]));
           
           QTreeWidgetItem* z1Node = new QTreeWidgetItem(valueNode);
           z1Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           z1Node->setText(0, "z1");
           z1Node->setText(1, QString("%1").arg(ashape->para[3]));
           
           QTreeWidgetItem* z2Node = new QTreeWidgetItem(valueNode);
           z2Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
           z2Node->setText(0, "z2");
           z2Node->setText(1, QString("%1").arg(ashape->para[4]));
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
  case PXPLANE:
    {
      valueNode->setText(0, "x_plus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "x1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;

 case NXPLANE:
    {
      valueNode->setText(0, "x_minus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "x1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;
  case PYPLANE:
    {
      valueNode->setText(0, "y_plus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "y1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;
  case NYPLANE:
    {
      valueNode->setText(0, "y_minus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "y1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;
  case PZPLANE:
    {
      valueNode->setText(0, "z_plus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "z1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;
  case NZPLANE:
    {
      valueNode->setText(0, "z_minus_plane");
      QTreeWidgetItem* x0Node = new QTreeWidgetItem(valueNode);
      x0Node->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled );
      x0Node->setText(0, "z1");
      x0Node->setText(1, QString("%1").arg(ashape->para[0]));
    }
    break;
  default:
    break;
  }











       
     }
     tree->expandAll();
     tree->resizeColumnToContents(0);
}



void FVMAdapt::addRegion(){

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


void FVMAdapt::removeNode(){
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
}
   

  
bool FVMAdapt::saveXml(){
  QDomDocument doc = toDom();
  if(doc.isNull())return false;
  QString fileName = filename.left(filename.lastIndexOf('.'))+".xml";

    fileName = QFileDialog::getSaveFileName(this, tr("Save .xml File"),
                                                 fileName,
                                                  tr("xml Files (*.xml)"));
 
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

  return true;
}


FVMAdapt::~FVMAdapt(){
   for(unsigned int i = 0 ; i< defaultShapes.size(); i++){
    if(defaultShapes[i]){
      delete defaultShapes[i];
      defaultShapes[i] = 0;
    }
  }
  defaultShapes.clear();
 
}


QDomNode makeElement(QDomDocument& doc, QTreeWidgetItem* item){

  if(item->text(0)=="translate"){
    double x0 = item->child(0)->text(1).toDouble();
    double y0 = item->child(1)->text(1).toDouble();
    double z0 = item->child(2)->text(1).toDouble();
    if(x0 == 0 && y0 == 0 && z0==0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      if(x0 != 0){
        QDomElement childelem = doc.createElement("x0");
         QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(y0 != 0){
        QDomElement childelem = doc.createElement("y0");
         QDomNode txtNode = doc.createTextNode(item->child(1)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(z0 != 0){
        QDomElement childelem = doc.createElement("z0");
        QDomNode txtNode = doc.createTextNode(item->child(2)->text(1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      return elem;
    }
  }else if(item->text(0)=="scale"){
    double x0 = item->child(0)->text(1).toDouble();
    double y0 = item->child(1)->text(1).toDouble();
    double z0 = item->child(2)->text(1).toDouble();
    if(x0 == 1 && y0 == 1 && z0==1)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      if(x0 != 1){
        QDomElement childelem = doc.createElement("x0");
         QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(y0 != 1){
        QDomElement childelem = doc.createElement("y0");
         QDomNode txtNode = doc.createTextNode(item->child(1)->text(1));
         childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      if(z0 != 1){
        QDomElement childelem = doc.createElement("z0");
        QDomNode txtNode = doc.createTextNode(item->child(2)->text(1));
        childelem.appendChild(txtNode);
        if(!childelem.isNull()) elem.appendChild(childelem);
      }
      return elem;
    }
  }else if(item->text(0)=="rotateX"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }else if(item->text(0)=="rotateY"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }else if(item->text(0)=="rotateZ"){
    double x0 = item->child(0)->text(1).toDouble();
    
    if(x0 == 0)return QDomNode();
    else{
      QDomElement elem = doc.createElement(item->text(0));
      QDomElement childelem = doc.createElement("theta");
      QDomNode txtNode = doc.createTextNode(item->child(0)->text(1));
      childelem.appendChild(txtNode);
      if(!childelem.isNull()) elem.appendChild(childelem);
      return elem;
    }
  }
      
  
  QDomElement elem = doc.createElement(item->text(0));
  if(item->text(1)!=""){
    QDomNode elt = doc.createTextNode(item->text(1));
    elem.appendChild(elt);
    return elem;
  }else if(item->childCount()!=0){
    for(int i = 0; i < item->childCount(); i++){
      QDomNode childElem = makeElement(doc, item->child(i));
      if(!childElem.isNull())elem.appendChild(childElem);
    }
    if(item->text(0)=="transform" && elem.firstChildElement().isNull()) return QDomNode();
    return elem;
  }
  return QDomNode();
}


QDomDocument FVMAdapt::toDom(){
  QDomDocument doc;
  QDomNode rootElem = makeElement(doc, root);
  doc.appendChild(rootElem);
  return doc;
  //QString xml = doc.toString();
  //qDebug() << xml;
}

void FVMAdapt::helpClicked(){
}

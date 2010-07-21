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
#include <QtDebug>
#include <QGridLayout>
#include <QCheckBox>
#include <QMessageBox>
#include <QButtonGroup>

#include <QSignalMapper>
#include <QInputDialog>
#include "pbwindow.h"
#include "grid.h"
#include "pages.h"
#include "progressdialog.h"

#include <QMainWindow>
#include <QToolBar>
#include <QRadioButton>
#include <QButtonGroup>

#define PI 3.14159265358979323846264338327950




//////////////////////////////////////////////////////////
//  



void PbWindow::createToolBar(){
  int spacing =2;  
  toolbar = addToolBar(tr("tree&vis"));
 
  QGroupBox* visbar = new QGroupBox();
  visbar->setFlat(true);

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

  toolbar->setLayoutDirection(Qt::RightToLeft);
  toolbar->addWidget(visbar);

}
  


  
void PbWindow::createFlowBar(){
  int spacing =20;  
  //create flowbar
 
  
  QGroupBox* flowbar = new QGroupBox("flow bar");
  QVBoxLayout* barLayout = new QVBoxLayout;
  
  QPushButton* loadButton = new QPushButton(tr("Load Grid"), this);
  barLayout->addWidget(loadButton);
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadGrid()));
  
  barLayout->addSpacing(spacing);
  

  QPushButton* setParButton = new QPushButton(tr("Set Parameters"), this);
  barLayout->addWidget(setParButton);
  connect(setParButton, SIGNAL(clicked()), this, SLOT(setParameter()));
  
  barLayout->addSpacing(spacing);

  QPushButton* splitButton = new QPushButton(tr("Split Faces"), this);
  barLayout->addWidget(splitButton);
  connect(splitButton, SIGNAL(clicked()), this, SLOT(split()));

  barLayout->addSpacing(spacing);
  QPushButton* doneButton = new QPushButton(tr("Done\n"), this);
  barLayout->addWidget(doneButton);
  connect(doneButton, SIGNAL(clicked()), this, SLOT(done()));

  barLayout->addStretch(10);
  flowbar->setLayout(barLayout);
  
  
  QToolBar* flowToolBar = new QToolBar;
  addToolBar(Qt::LeftToolBarArea,flowToolBar );
  flowToolBar->addWidget(flowbar);
}


void PbWindow::loadGrid(){
  filename =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                 QDir::currentPath(),
                                 tr("vog Files (*.vog)"));
 
  if(!filename.isEmpty()){
    if(viewer->load_boundary(filename, bnames, bids)) updateBoundaryView(bnames, bids);
  }
}
void PbWindow::done(){

  close();
}
void PbWindow::split(){
  QString initialPath = filename.section('/', 0, -2)+"/";
  QString outFile = QFileDialog::getSaveFileName(this,  tr("output file name"),
                                                 initialPath,
                                                 tr("Volume Grid File(*.vog)")
                                                 );
  if(outFile.isEmpty())return;
  QString inFile =  filename.section('.', 0, -2);
  QString command = "pbv2v -o " + outFile + QString(" -b %1 %2").arg(bids[b1]).arg(bids[b2]);
  
  
  int optionIndex = buttonGroup->checkedId();
  
  if(optionIndex ==0){
    positions3d t  = translate->value();
    command += QString(" -t %1,%2,%3").arg(t.x).arg(t.y).arg(t.z);
    
  }else if(optionIndex==1){
    positions3d center = rotateCenter->value();
    positions3d axis = rotateAxis->value();
    double angle = rotateAngle->value();
    command += QString(" -r %1,%2,%3 %4,%5,%6 %7").arg(center.x).arg(center.y).arg(center.z)
      .arg(axis.x).arg(axis.y).arg(axis.z).arg(angle);
  }

  command += " " +inFile;
  
  ProgressDialog* progress = new ProgressDialog(command, initialPath, false);
  progress->show();
  
}
void PbWindow::selectCurrent(int row){
  
  QColor value= modBoundaries->item(row, 0)->background().color();
  QModelIndex index =qobject_cast<const QAbstractItemModel*>(modBoundaries)->index(row, 1);
  qobject_cast<QAbstractItemView*>( boundaryView)->setCurrentIndex(index);
    
  viewer->setCurrentObj(row, value);
  if(boundaryGroup->checkedId()==0){
    boundaryGroup->button(0)->setText("Boundary 1: "+bnames[row]);
    b1 = row;
  }else{
    boundaryGroup->button(1)->setText("Boundary 2: "+bnames[row]);
    b2 = row;
  }
 
}
void PbWindow::showBoundary(QModelIndex top, QModelIndex ){
 
  if(top.column() ==3){//visibility item
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
  if(boundaryGroup->checkedId()==0) {
    boundaryGroup->button(0)->setText("Boundary 1: "+bnames[top.row()]);
    b1 = top.row();
  
  }else{
    boundaryGroup->button(1)->setText("Boundary 2: "+bnames[top.row()]);
    b2 = top.row();
   
  }

}

void PbWindow::setCurrent(QModelIndex top){
  QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), 0)->background().color();
  viewer->setCurrentObj(top.row(), value);
  if(boundaryGroup->checkedId()==0){
    boundaryGroup->button(0)->setText("Boundary 1: "+bnames[top.row()]);
    b1 =top.row();
   
  }else{
    boundaryGroup->button(1)->setText("Boundary 2: "+bnames[top.row()]);
    b2 = top.row();
   
  }
}



void PbWindow::updateBoundaryView(const QStringList& bdnames, const QList<int>& bdids){
  //clean up the data
  if(modBoundaries){
    delete modBoundaries;
    modBoundaries = 0;
  }
  if(boundaryView){
    objLayout->removeWidget(boundaryView);
    delete boundaryView;
    boundaryView = 0;
  }
      
   
  
  selectBoundary(bdnames, bdids);
  
}

bool PbWindow::selectBoundary(const QStringList& bdnames, const QList<int> &bdids){
      
  // Get boundary names from topo file
  if(bdnames.empty()){
    QMessageBox::warning(this, tr("select Boundary"),
                         tr("please  load boundaries first"));
    return false;
  }
  
 
  modBoundaries = new QStandardItemModel(bdnames.size(), 4, this);
  
  modBoundaries->setHeaderData(0, Qt::Horizontal, QObject::tr("color"));
  modBoundaries->setHeaderData(1, Qt::Horizontal, QObject::tr("boundary name"));
  modBoundaries->setHeaderData(2, Qt::Horizontal, QObject::tr("boundary id"));
  modBoundaries->setHeaderData(3, Qt::Horizontal, QObject::tr("show/hide"));
  

    
    
  
  for (int i = 0; i < bdnames.size(); ++i) {
    QColor newColor = default_color[i%12];
    QStandardItem* colorItem = new QStandardItem("");
    QStandardItem* nameItem = new QStandardItem(bdnames[i]);
    QStandardItem* idItem = new QStandardItem(QString("%1").arg(bdids[i]));
    QStandardItem* showItem = new QStandardItem("show");
    colorItem->setBackground(QBrush(newColor));
      
    nameItem->setFlags(Qt::ItemIsSelectable | 
                       Qt::ItemIsUserCheckable | 
                       Qt::ItemIsEnabled);
      
    modBoundaries->setItem(i, 0, colorItem);
    modBoundaries->setItem(i, 1, nameItem);
    modBoundaries->setItem(i, 2, idItem);
    modBoundaries->setItem(i, 3, showItem);
  }
    
    

  QItemSelectionModel *selections = new QItemSelectionModel(modBoundaries, this);

 

  // Construct and show dock widget
  showDelegate* delBoundaries = new showDelegate(this);
  colorDelegate* delColor = new colorDelegate(this);
 

  boundaryView = new QTableView(this);
  objLayout->insertWidget(0, boundaryView);

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
  
  connect(viewer, SIGNAL(pickCurrent(int)), this, SLOT(selectCurrent(int)));
 

  // selectCurrent(currentIndex);


  return true;
}


PbWindow::PbWindow( QWidget *parent):QMainWindow(parent){
  QWidget::setAttribute(Qt::WA_DeleteOnClose, true);

  QGroupBox* central = new QGroupBox;
  central->setFlat(true);
 
  QGroupBox* boundaries = new QGroupBox;
  QHBoxLayout* hLayout = new QHBoxLayout;
  QPushButton* boundary1 = new QPushButton(tr("Boundary 1: " ));
  QPushButton* boundary2 = new QPushButton(tr("Boundary 2: "));
  boundary1->setCheckable(true);
  boundary2->setCheckable(true);
  boundary1->setChecked(true);
  boundaryGroup = new QButtonGroup(this);
  
  boundaryGroup->addButton(boundary1, 0);
  boundaryGroup->addButton(boundary2, 1);
  
  
  
  hLayout->addWidget(boundary1);
  hLayout->addWidget(boundary2);
  boundaries->setLayout(hLayout);
  



  
  QRadioButton* tButton = new QRadioButton(tr("Extrusion"));
  QRadioButton* rButton = new QRadioButton(tr("Rotation"));
 
  
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(tButton);
  buttonLayout->addWidget(rButton);
 
  buttonGroup = new QButtonGroup(this);
  buttonGroup->addButton(tButton, 0);
  buttonGroup->addButton(rButton, 1);
 
  tButton->setChecked(true);
  
  

  QGroupBox* translateGroup = new QGroupBox("Extrusion");
  {
    QHBoxLayout* translateLayout = new QHBoxLayout;
    translate = new VectSpBox("Translate");
    translateLayout->addWidget(translate);
    translateGroup->setLayout(translateLayout);
  }
    


  QGroupBox* rotationGroup = new QGroupBox("Rotation");
  {
    rotateCenter = new VectSpBox("Rotate Center");
    rotateAxis= new VectSpBox("Rotate Axis");
    rotateGroup = new QGroupBox("Rotate Angle");
    QHBoxLayout* rLayout = new QHBoxLayout;
    rotateAngle = new DoubleEdit;
    rLayout->addWidget(rotateAngle);
    rotateGroup->setLayout(rLayout);
  
    QVBoxLayout* rotationLayout = new QVBoxLayout;
    rotationLayout->addWidget(rotateCenter);
    rotationLayout->addWidget(rotateAxis);
    rotationLayout->addWidget(rotateGroup);
    rotationGroup->setLayout(rotationLayout);
  }
    
  QGroupBox* transform = new QGroupBox(tr("Transform"));
  QVBoxLayout* transLayout = new QVBoxLayout;
  transLayout->addLayout(buttonLayout);
  transLayout->addWidget(translateGroup);
  transLayout->addWidget(rotationGroup);
  transform->setLayout(transLayout);
  
  viewer = new GLViewer();
  boundaryView = new QTableView(this);
  
  objLayout = new QVBoxLayout;
  objLayout->addWidget(boundaryView);
  objLayout->addWidget(boundaries);
  objLayout->addWidget(transform);

  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->addLayout(objLayout);
  mainLayout->addWidget(viewer);
  central->setLayout(mainLayout);
  

   
  createFlowBar();
  createToolBar();
  setCentralWidget(central);
  setWindowTitle(tr("PbWindow"));
  setMinimumSize(1000, 700);
}


      

void PbWindow::changePage(int index){
  if(index == 0){
    translate->setDisabled(true);
    rotateAxis->setDisabled(true);
    rotateGroup->setDisabled(true);
  }else if(index==1){
    translate->setDisabled(true);
    rotateAxis->setDisabled(false);
    rotateGroup->setDisabled(true);
  }else{
    translate->setDisabled(false);
    rotateAxis->setDisabled(false);
    rotateGroup->setDisabled(false); 
  }
}

void PbWindow::setParameter(){
  int optionIndex = buttonGroup->checkedId();
  if(optionIndex ==0){
    translate->setValue(viewer->getTranslate(b1, b2));
  }else if(optionIndex==1){
    positions3d axis = positions3d(0, 0, 1);
    positions3d center = positions3d(0, 0, 0);
    double angle = 0;
    viewer->getRotate(b1, b2, angle, axis, center);
    if(fabs(angle)<1e-3){
       QMessageBox::warning(this, tr("rotation parameter"),
                            tr("rotation angle is zero, extrusion is used"));

       buttonGroup->button(0)->setChecked(true);
       setParameter();
       return;
    }
    rotateAngle->setValue(angle);
    rotateAxis->setValue(axis);
    rotateCenter->setValue(center);
  }
}
  

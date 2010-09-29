#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <QToolBar>
#include <QDockWidget>
#include <QScrollArea>
#include <stdlib.h>
#include <unistd.h>
#include "vmergewindow.h"
#include "pages.h"
#include "helpwindow.h"
using namespace std;


VMOption::VMOption(int id, const QString &gridname, const QStringList & bcnames, QWidget *parent ) : QGroupBox(gridname, parent), gridId(id){
 
  tc.push_back(TranCoef2());
  currentCoef = 0;
  gridName = gridname;
  bdnames.clear();
  needGlue.clear();
  //currentBd = -1;
  for(int i = 0; i < bcnames.size(); i++){
    bdnames.push_back(pair<QString, QString>(bcnames[i], ""));
    needGlue.push_back(false);
  }

  QGroupBox* transferBox = new QGroupBox(tr("Transform")); 
  QVBoxLayout* transferLayout = new QVBoxLayout;


  translateVec = new VectSpBox("translation");

  QGroupBox* rotateBox = new QGroupBox("rotation");
  angleVec  =  new VectSpBox("angle");
  angleVec->setRange(-360.0, 360.0);
  centerVec  =  new VectSpBox("center");
  
  QVBoxLayout* rotateLayout = new QVBoxLayout;
  centerVec  =  new VectSpBox("center");
  rotateLayout->addWidget(centerVec);
  rotateLayout->addWidget(angleVec);
  rotateBox->setLayout(rotateLayout);
  
  scaleVec =   new VectSpBox("scale");
  scaleVec->setValue(positions3d(1.0, 1.0, 1.0));


  
  QGroupBox* mirrorBox = new QGroupBox(tr("mirror"));
  QVBoxLayout* mirrorLayout = new QVBoxLayout;
  mButtonGroup = new QButtonGroup;
  
  QHBoxLayout* mButtonLayout = new QHBoxLayout;
  QRadioButton* radion = new QRadioButton("no mirror");
  QRadioButton* radiox = new QRadioButton("x=0 plane");
  QRadioButton* radioy = new QRadioButton("y=0 plane");
  QRadioButton* radioz = new QRadioButton("z=0 plane");
  radiog = new QRadioButton("general plane");
  mButtonLayout->addWidget(radion);
  mButtonLayout->addWidget(radiox);
  mButtonLayout->addWidget(radioy);
  mButtonLayout->addWidget(radioz);
  mButtonLayout->addWidget(radiog);
  mButtonGroup->addButton(radion, 0);
  mButtonGroup->addButton(radiox, 1);
  mButtonGroup->addButton(radioy, 2);
  mButtonGroup->addButton(radioz, 3);
  mButtonGroup->addButton(radiog, 4);
  mirrorLayout->addLayout(mButtonLayout);
  radion->setChecked(true);

  
  connect(radiog, SIGNAL(clicked()), this, SLOT(gRadioClicked()));
  planeBox = new QGroupBox(tr("3d plane")); 
  
  QVBoxLayout* planeLayout = new QVBoxLayout;
  planeOriginVec = new VectSpBox("origin");
  planeNormalVec = new VectSpBox("normal");
  planeNormalVec->setValue(positions3d(1.0, 0.0, 0.0));
  planeLayout->addWidget(planeOriginVec);
  planeLayout->addWidget(planeNormalVec);
  planeBox->setLayout(planeLayout);
  
  mirrorLayout->addWidget(planeBox);
  mirrorBox->setLayout(mirrorLayout);                                  
  planeBox->hide();
 
 
  

  QHBoxLayout* buttonLayout = new QHBoxLayout;
  QPushButton* acceptButton = new QPushButton(tr("Add"));
  QPushButton* cancelButton = new QPushButton(tr("Remove"));
  QPushButton* clearButton = new QPushButton(tr("Clear All"));
  QPushButton* previousButton = new QPushButton(tr("Previous"));
  QPushButton* nextButton = new QPushButton(tr("Next"));
  buttonLayout->addWidget(acceptButton);
  buttonLayout->addWidget(cancelButton);
  buttonLayout->addWidget(previousButton);
  buttonLayout->addWidget(nextButton);
  buttonLayout->addWidget(clearButton);
  connect(acceptButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(cancel()));
  connect(previousButton, SIGNAL(clicked()), this, SLOT(previous()));
  connect(nextButton, SIGNAL(clicked()), this, SLOT(next()));
  connect(clearButton, SIGNAL(clicked()), this, SLOT(clear()));
  
  transferLayout->addWidget(translateVec,1);
  transferLayout->addWidget(rotateBox,2);
  transferLayout->addWidget(scaleVec,1);
  transferLayout->addWidget(mirrorBox,3);
  transferLayout->addLayout(buttonLayout,0);
  
  transferBox->setLayout(transferLayout);
  
  
  QGroupBox* tagBox = new QGroupBox(tr("tag")); 
  QHBoxLayout* tagLayout = new QHBoxLayout;
 
  tagEditor = new QLineEdit;
  tagLayout->addWidget(new QLabel(tr("   ")));
  tagLayout->addWidget(tagEditor);
  tagBox->setLayout(tagLayout);

  tolBox = new QGroupBox(tr("tolerance for gluing")); 
  QHBoxLayout* tolLayout = new QHBoxLayout;
  
  tolEdit = new DoubleEdit;
  tolLayout->addWidget(new QLabel(tr("   ")));
  tolLayout->addWidget(tolEdit);
  tolBox->setLayout(tolLayout);
  tolEdit->setBottom(0.0);
  tolBox->setCheckable(true);
  tolBox->setChecked(false); 
  
  
  
  connect(translateVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(angleVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(centerVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(scaleVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(planeOriginVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(planeNormalVec, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(mButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(setInfo()));
  // connect(tagEditor, SIGNAL(textChanged(const QString&)), this, SLOT(setInfo()));

  

  modBoundaries = new QStandardItemModel(bdnames.size(), 5, this);
  modBoundaries->setHeaderData(0, Qt::Horizontal, QObject::tr("color"));
  modBoundaries->setHeaderData(1, Qt::Horizontal, QObject::tr("name"));
  modBoundaries->setHeaderData(2, Qt::Horizontal, QObject::tr("glue"));
  modBoundaries->setHeaderData(3, Qt::Horizontal, QObject::tr("new name"));
  modBoundaries->setHeaderData(4, Qt::Horizontal, QObject::tr("show/hide"));
 
 
  for (int i = 0; i < bdnames.size(); ++i) {
    QColor newColor = default_color[i%12];
    QStandardItem* colorItem = new QStandardItem("");
    QStandardItem* nameItem = new QStandardItem(bdnames[i].first);
    QStandardItem* showItem = new QStandardItem("show");
    QStandardItem* glueItem = new QStandardItem("no glue"); 
    QStandardItem* newNameItem = new QStandardItem("");
    colorItem->setBackground(QBrush(newColor));
      
    nameItem->setFlags(Qt::ItemIsSelectable | 
                       Qt::ItemIsUserCheckable | 
                       Qt::ItemIsEnabled);
      
    modBoundaries->setItem(i, 0, colorItem);
    modBoundaries->setItem(i, 1, nameItem);
    modBoundaries->setItem(i, 2, glueItem);
    modBoundaries->setItem(i, 3, newNameItem);
    modBoundaries->setItem(i, 4, showItem);
  }
    
  QItemSelectionModel *selections = new QItemSelectionModel(modBoundaries, this);

 

  // Construct and show dock widget
  showDelegate* delBoundaries = new showDelegate(this);
  colorDelegate* delColor = new colorDelegate(this);
  showDelegate* delGlue = new showDelegate(this);
 

  boundaryView = new QTableView(this);

  boundaryView->setModel(modBoundaries);
  boundaryView->setSelectionModel(selections);
  boundaryView->setSelectionMode(QAbstractItemView::SingleSelection);
  boundaryView->setItemDelegateForColumn(4,delBoundaries);
  boundaryView->setItemDelegateForColumn(0,delColor);
  boundaryView->setItemDelegateForColumn(2,delGlue);
  
  boundaryView->setColumnWidth(0, 20);
  boundaryView->setWordWrap(false);
  boundaryView->resizeRowsToContents();
  boundaryView->resizeColumnsToContents();
 
  
  connect(modBoundaries, SIGNAL(dataChanged( const QModelIndex&, const QModelIndex&)),
          this, SLOT(showBoundary(QModelIndex, QModelIndex)));
  
  connect(boundaryView, SIGNAL(clicked( const QModelIndex&)),
          this, SLOT(setCurrentObj(QModelIndex)));

  QVBoxLayout *leftLayout = new QVBoxLayout;
  leftLayout->addWidget(tagBox);
   leftLayout->addWidget(tolBox);
  leftLayout->addWidget(boundaryView);
  
  QHBoxLayout *mainLayout =  new QHBoxLayout;
 
  mainLayout->addWidget(transferBox);
  mainLayout->addLayout(leftLayout);

  setLayout(mainLayout);
 
}
void VMOption::accept(){
  tc.push_back(TranCoef2());
  currentCoef = tc.size()-1;
  cancel();
 
}

void VMOption::cancel(){
  translateVec->setValue(positions3d(0.0, 0.0, 0.0));
  angleVec->setValue(positions3d(0.0, 0.0, 0.0));
  centerVec->setValue(positions3d(0.0, 0.0, 0.0));
  scaleVec->setValue(positions3d(1.0, 1.0, 1.0));
  planeOriginVec->setValue(positions3d(0.0, 0.0, 0.0));
  planeNormalVec->setValue(positions3d(1.0, 0.0, 0.0));
  mButtonGroup->button(0)->setChecked(true);

}

void VMOption::clear(){
  tc.clear();
  tc.push_back(TranCoef2());
  currentCoef = 0;
  cancel();
}
  
  
void VMOption::previous(){
  if(currentCoef == 0) return;
 
  currentCoef--;

  setValue(tc[currentCoef]);
 
}

void VMOption::setValue(const TranCoef2& t){
  translateVec->setValue(t.translate);
  angleVec->setValue(t.rotateAngle);
  centerVec->setValue(t.rotateCenter);
  scaleVec->setValue(t.scale);
  planeOriginVec->setValue(t.mirrorOrigin);
  planeNormalVec->setValue(t.mirrorNormal);
  mButtonGroup->button(t.checkedId)->setChecked(true);
}



  
void VMOption::next(){
  if(currentCoef >=(int)(tc.size()-1))return;
  currentCoef++;
  setValue(tc[currentCoef]);
}
  

void VMOption::showBoundary(QModelIndex top, QModelIndex ){

  if(top.column() ==4){//visibility item
    QString value = top.data(Qt::EditRole).toString();
    
    if(value == "show"){
      IDVisibility idvis = IDVisibility(gridId, top.row(), true);  
      emit setCurrentVisibility(idvis);
    }else{
      IDVisibility idvis = IDVisibility(gridId, top.row(), false);  
      emit setCurrentVisibility(idvis);
      
    }
  }else if(top.column() ==0) {//color item
   
    QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), top.column())->background().color();
     IDColor ic = IDColor(gridId, top.row(), value);
     emit setCurrentColor(ic);
  }else if(top.column() ==3){
      QString value = top.data(Qt::EditRole).toString();
      bdnames[top.row()].second = value;
  }else if(top.column() ==2){
    QString value = top.data(Qt::EditRole).toString();
    needGlue[top.row()] = (value=="glue");
  }

}

void VMOption::setCurrentObj(QModelIndex top){
     QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), 0)->background().color();
     IDColor ic = IDColor(gridId, top.row(), value);
     emit setCurrentColor(ic);
}

void VMOption::setCurrentBid(int bid){
    QColor value= modBoundaries->item(bid, 0)->background().color();
    QModelIndex index =qobject_cast<const QAbstractItemModel*>(modBoundaries)->index(bid, 1);
    qobject_cast<QAbstractItemView*>( boundaryView)->setCurrentIndex(index);

    IDColor ic = IDColor(gridId, bid, value);
    emit setCurrentColor(ic);
   
}

void VMOption::setInfo(){
 
  tc[currentCoef].rotateAngle = angleVec->value();
  tc[currentCoef].rotateCenter = centerVec->value();
  tc[currentCoef].translate = translateVec->value();
  tc[currentCoef].scale = scaleVec->value();
  tc[currentCoef].checkedId = mButtonGroup->checkedId();
  tc[currentCoef].mirrorOrigin = planeOriginVec->value();
  tc[currentCoef].mirrorNormal = planeNormalVec->value();
  
  IDMatrix idM = IDMatrix(gridId, currentM()); 
  emit tcChanged(idM);
  
}

affineMapping2 VMOption::currentM(){
  affineMapping2 gridXform;
 
  for(unsigned int i = 0; i < tc.size(); i++){
    if(norm(tc[i].translate)!=0) gridXform.translate(tc[i].translate) ;


    if(norm(tc[i].rotateAngle)!=0){
      if(norm(tc[i].rotateCenter)!=0) gridXform.translate(tc[i].rotateCenter) ;
      if(tc[i].rotateAngle.x !=0)gridXform.rotateX(tc[i].rotateAngle.x);
      if(tc[i].rotateAngle.y !=0) gridXform.rotateY(tc[i].rotateAngle.y);
      if(tc[i].rotateAngle.z !=0)gridXform.rotateZ(tc[i].rotateAngle.z);
      if(norm(tc[i].rotateCenter)!=0)
        gridXform.translate(vector3d<double>(-tc[i].rotateCenter.x,
                                             -tc[i].rotateCenter.y,
                                             -tc[i].rotateCenter.z)) ;
    }
    
    if(tc[i].scale.x != 1 ||tc[i].scale.y != 1||tc[i].scale.z != 1 ) gridXform.scale(tc[i].scale);
    if(tc[i].checkedId !=0){
      switch(tc[i].checkedId){
      case 1:
        gridXform.mirror(positions3d(0.0, 0.0, 0.0), positions3d(1.0, 0.0, 0.0));
        break;
      case 2:
        gridXform.mirror(positions3d(0.0, 0.0, 0.0), positions3d(0.0, 1.0, 0.0));
        break;
      case 3:
         gridXform.mirror(positions3d(0.0, 0.0, 0.0), positions3d(0.0, 0.0, 1.0));
         break;
      case 4:
        gridXform.mirror(tc[i].mirrorOrigin, tc[i].mirrorNormal);
        break;
      default:
        break;
      }
    }
  }
  
  return gridXform;
}

void VMOption::gRadioClicked(){
  if(radiog->isChecked())planeBox->show();
  else planeBox->hide();
}

QString VMOption::currentText(){
  if(bdnames.size() == 0)return QString();
  QString text;
  text += " -g " + gridName;
  for(unsigned int i = 0; i < tc.size(); i++){
    if(tc[i].translate.x != 0) text += QString(" -xshift %1").arg(tc[i].translate.x);
    if(tc[i].translate.y != 0) text += QString(" -yshift %1").arg(tc[i].translate.y);
    if(tc[i].translate.z != 0) text += QString(" -zshift %1").arg(tc[i].translate.z);
    
    if(norm(tc[i].rotateAngle)!=0){
      if(tc[i].rotateCenter.x != 0) text += QString(" -xshift %1").arg(tc[i].rotateCenter.x);
      if(tc[i].rotateCenter.y != 0) text += QString(" -yshift %1").arg(tc[i].rotateCenter.y);
      if(tc[i].rotateCenter.z != 0) text += QString(" -zshift %1").arg(tc[i].rotateCenter.z);
      
      if(tc[i].rotateAngle.x != 0) text += QString(" -xrotate %1").arg(tc[i].rotateAngle.x);
      if(tc[i].rotateAngle.y != 0) text += QString(" -yrotate %1").arg(tc[i].rotateAngle.y);
      if(tc[i].rotateAngle.z != 0) text += QString(" -zrotate %1").arg(tc[i].rotateAngle.z);
      
      if(tc[i].rotateCenter.x != 0) text += QString(" -xshift %1").arg(-tc[i].rotateCenter.x);
      if(tc[i].rotateCenter.y != 0) text += QString(" -yshift %1").arg(-tc[i].rotateCenter.y);
      if(tc[i].rotateCenter.z != 0) text += QString(" -zshift %1").arg(-tc[i].rotateCenter.z);
    }
    if(tc[i].scale.x != 1) text += QString(" -xscale %1").arg(tc[i].scale.x);
    if(tc[i].scale.y != 1) text += QString(" -yscale %1").arg(tc[i].scale.y);
    if(tc[i].scale.z != 1) text += QString(" -zscale %1").arg(tc[i].scale.z);
    if(tc[i].checkedId !=0){
      switch(tc[i].checkedId){
      case 1:
        text += " -mirrorx ";
        break;
      case 2:
         text += " -mirrory ";
        break;
      case 3:
          text += " -mirrorz ";
         break;
      case 4:
        text +=  QString(" -mirror %1,%2,%3,%4,%5,%6")
          .arg(tc[i].mirrorOrigin.x).arg(tc[i].mirrorOrigin.y).arg(tc[i].mirrorOrigin.z)
          .arg(tc[i].mirrorNormal.x).arg(tc[i].mirrorNormal.y).arg(tc[i].mirrorNormal.z);
       
        break;
      default:
        break;
      }
    }

      
  }
  if(!tagEditor->text().isEmpty()) text += " -tag " + tagEditor->text();
  if(tolBox->isChecked()) text += " -tol " + QString("%1").arg(tolEdit->value());


  for(int i =0; i < bdnames.size(); i++){
    if(!bdnames[i].second.isEmpty()){
      text += " -bc " +bdnames[i].first + "," + bdnames[i].second;
    }
    if(needGlue[i]) text += " -glue " + bdnames[i].first;
  }     
  return text;
}
      


void VMergeWindow::createFlowBar(){
  int spacing =2;  
  //create flowbar
 
  
  QGroupBox* flowbar = new QGroupBox("flow bar");
  QVBoxLayout* barLayout = new QVBoxLayout;
  
  QPushButton *loadButton = new QPushButton(tr("&load grid"));
  QPushButton *mergeButton = new QPushButton(tr("&vogmerge"));
  QPushButton *clearButton = new QPushButton(tr("&restart" ));
  
  
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadGridClicked()));
  connect(mergeButton, SIGNAL(clicked()), this, SLOT(vmClicked()));
  connect(clearButton, SIGNAL(clicked()), this, SLOT(clearAll()));
  
  barLayout->addSpacing(spacing);
  barLayout->addWidget(loadButton);
  
  
  barLayout->addSpacing(spacing);
  barLayout->addWidget(mergeButton);
  
  barLayout->addSpacing(spacing);
  barLayout->addWidget(clearButton);
  
  barLayout->addStretch(10);
  flowbar->setLayout(barLayout);
  
  
  QToolBar* flowToolBar = new QToolBar;
  addToolBar(Qt::LeftToolBarArea,flowToolBar );
  flowToolBar->addWidget(flowbar);
}



VMergeWindow::VMergeWindow( QWidget* parent)
  :QMainWindow(parent)
{
  setAttribute(Qt::WA_DeleteOnClose, true);
  
     
  QGroupBox* central = new QGroupBox;
  central->setFlat(true);
  
 
  typesWidget = new QComboBox;
  pagesWidget = new QStackedWidget;
  QVBoxLayout* objLayout = new QVBoxLayout;
  objLayout->addWidget(typesWidget);
  objLayout->addWidget(pagesWidget);
  
  mgviewer = new MGViewer();
  mgviewer->setMinimumSize(700, 700);
  
  
  connect(typesWidget,
          SIGNAL(currentIndexChanged(int)),
          this, SLOT(changePage(int)));
 
  
  central->setLayout(objLayout);

  
  
  

 
  QDockWidget*  viewerDock  = new QDockWidget("Graphics Viewer", this); 
  viewerDock->setAllowedAreas(Qt::RightDockWidgetArea );
  viewerDock->setWidget(mgviewer);
  viewerDock->setFloating(true);
  addDockWidget(Qt::RightDockWidgetArea, viewerDock);

  
  
 
  connect(mgviewer, SIGNAL(gridLoaded(const QStringList&)), this, SLOT(gridLoaded(const QStringList&)));
  connect(mgviewer, SIGNAL(pickCurrent(const IDOnly&)), this, SLOT(selectCurrent(const IDOnly&)));
  connect(this, SIGNAL(tcChanged(const IDMatrix&)), mgviewer, SLOT(transGrid(const IDMatrix&)));
 
  connect(this, SIGNAL(setCurrentColor(const IDColor&)), mgviewer, SLOT(setCurrentColor(const IDColor&)));
  connect(this, SIGNAL(setCurrentVisibility(const IDVisibility&)),
          mgviewer, SLOT(setCurrentVisibility(const IDVisibility&)));
 
  
  
  createFlowBar();
  createToolBar();
  
  setCentralWidget(central);
  setWindowTitle(tr("vogmerge window"));
  setMinimumSize(1000, 700);
 
}
void VMergeWindow::createToolBar(){
 
  QToolBar* toolbar = addToolBar(tr("tree&vis"));
 
  QGroupBox* visbar = new QGroupBox();
  visbar->setFlat(true);

  QHBoxLayout* visLayout = new QHBoxLayout;
  
  QPushButton* helpButton = new QPushButton(tr("help"));
  connect(helpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  visLayout->addWidget(helpButton); 
  
  QPushButton *clearBoundaryAct = new QPushButton(tr("Clear"), this);
  visLayout->addWidget(clearBoundaryAct);
  connect(clearBoundaryAct, SIGNAL(clicked()),
          mgviewer, SLOT(clearCurrent())); 
   
  QPushButton* resetAct = new QPushButton(tr("Reset"), this);
  visLayout->addWidget(resetAct);
  connect(resetAct, SIGNAL(clicked()),
          mgviewer, SLOT(reset()));
 
  QPushButton* fitAct = new QPushButton(tr("Fit"), this);
  visLayout->addWidget(fitAct);
  connect(fitAct, SIGNAL(clicked()),
          mgviewer, SLOT(fit()));
  
  visbar->setLayout(visLayout);

  toolbar->setLayoutDirection(Qt::RightToLeft);
  toolbar->addWidget(visbar);

}
  

 void VMergeWindow::helpClicked(){
   HelpWindow* helpwindow = new HelpWindow("page_vogmerge.html");
   helpwindow->show();
 }     

void VMergeWindow::changePage(int current)
{
  
  pagesWidget->setCurrentIndex(current);
}










void VMergeWindow::loadGridClicked(){
  
  QString format = "volume grid file(*.vog)";
  QString initialPath =QDir::currentPath();
  
  QString fileName  = initialPath;
 
  
  fileName = QFileDialog::getOpenFileName(this, tr("Load Grid"),
                                          fileName,
                                          format);
   if(fileName==""){
    //no error message in  case of 'cancel' is pressed 
      return;
   }

   //   emit loadGrid(fileName);
    mgviewer->load_boundary(fileName);
}

void VMergeWindow::gridLoaded(const QStringList & names){
  if(names.isEmpty()) return;
  QString filename = names[0];
 
  QStringList bdnames;
  for(int i = 1; i < names.size(); i++)bdnames<<names.at(i);
  typesWidget->addItem(filename);
  int id = typesWidget->count()-1;
  VMOption* agrid = new VMOption(id, filename, bdnames);
  connect(agrid, SIGNAL(setCurrentColor(const IDColor&)), this, SIGNAL(setCurrentColor(const IDColor&)));
  connect(agrid, SIGNAL(setCurrentVisibility(const IDVisibility&)), this, SIGNAL(setCurrentVisibility(const IDVisibility&)));
  connect(agrid, SIGNAL(tcChanged(const IDMatrix&)), this, SIGNAL(tcChanged(const IDMatrix&)));
  pagesWidget->addWidget(agrid);
  typesWidget->setCurrentIndex(id);
  
}
  
void VMergeWindow::selectCurrent(const IDOnly& id){
  
  typesWidget->setCurrentIndex(id.gridId);
  qobject_cast<VMOption*>(pagesWidget->widget(id.gridId))->setCurrentBid(id.boundId); 
 
}


void VMergeWindow::vmClicked(){
  QString initialPath =QDir::currentPath();
  
  QString fileName  = initialPath;
 
  
   fileName = QFileDialog::getSaveFileName(this, tr("Merged Vog File"),
                                                  fileName,
                                                  tr("Volume Grid files (*.vog)"));
  
  if(fileName.section('.', -1, -1)!="vog")fileName+=".vog";
  
  if(!(fileName.section('/', -1, -1).section('.',0,0).isEmpty())){
    QString command = QString("vogmerge");
    for(int i = 0; i < typesWidget->count(); i++){
      command +=qobject_cast<VMOption*>( pagesWidget->widget(i))->currentText();
    }
    command += " -o " + fileName;
   
 
    ProgressDialog* progress = new ProgressDialog(command,false);
    progress->show();
    connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus, QString)), this, SLOT(afterMerge(QString, QProcess::ExitStatus, QString)));
  }
    
 }

void VMergeWindow::clearAll(){
  clear();
  //emit getGrid(QString());
  mgviewer->get_boundary(QString());
}

void VMergeWindow::clear(){
  
  transcoef.clear();
  typesWidget->clear();
  QWidget* tmpWidget = pagesWidget->widget(0);
  while(tmpWidget){
    pagesWidget->removeWidget(tmpWidget);
    delete tmpWidget;
    tmpWidget = pagesWidget->widget(0);
  }

 
}

void VMergeWindow::afterMerge(QString command, QProcess::ExitStatus status, QString directory){
  if(status==QProcess::NormalExit){
    clear();
    QString filename =directory+ command.section(' ', -1, -1);
    // emit getGrid(filename);
     mgviewer->get_boundary(filename);
  }
}

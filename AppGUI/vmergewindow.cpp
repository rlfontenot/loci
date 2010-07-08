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
using namespace std;


VMOption::VMOption(int id, const QString &gridname, const QStringList & bcnames, QWidget *parent ) : QGroupBox(gridname, parent), gridId(id){
 
  tc.push_back(TranCoef());
  currentCoef = 0;
  gridName = gridname;
  bdnames.clear();
  //currentBd = -1;
  for(int i = 0; i < bcnames.size(); i++){
    bdnames.push_back(pair<QString, QString>(bcnames[i], ""));
  }

  QGroupBox* transferBox = new QGroupBox(tr("Transfer")); 
  QVBoxLayout* transferLayout = new QVBoxLayout;
  
  QGroupBox* translateBox = new QGroupBox(tr("translation")); 
  QGridLayout* translate = new QGridLayout;
  
  translate->addWidget(new QLabel(tr("x:")), 0, 1);
  translate->addWidget(new QLabel(tr("y:")), 0, 2);
  translate->addWidget(new QLabel(tr("z:")), 0, 3);
  
 
  
  xEditor1 = new DoubleEdit(0.0); 
  yEditor1 = new DoubleEdit(0.0);
  zEditor1 = new DoubleEdit(0.0); 

 
  translate->addWidget(xEditor1, 1, 1);
  translate->addWidget(yEditor1, 1, 2);
  translate->addWidget(zEditor1, 1, 3);
  translate->addWidget(new QLabel(tr("    ")), 1, 0);
  translate->addWidget(new QLabel(tr("    ")), 2, 0);
  
  translateBox->setLayout(translate);
  
  
  QGroupBox* rotateBox = new QGroupBox(tr("rotation")); 
  QGridLayout* rotate = new QGridLayout;
  
 
  rotate->addWidget(new QLabel(tr("x:")), 0, 1);
  rotate->addWidget(new QLabel(tr("y:")), 0, 2);
  rotate->addWidget(new QLabel(tr("z:")), 0, 3);
  
  //angles
  xEditor2 = new DoubleEdit(0.0); 
  yEditor2 = new DoubleEdit(0.0);
  zEditor2 = new DoubleEdit(0.0);
  xEditor2->setRange(-360.0, 360.0);
  yEditor2->setRange(-360.0, 360.0);
  zEditor2->setRange(-360.0, 360.0);
  xEditor4 = new DoubleEdit(0.0); 
  yEditor4 = new DoubleEdit(0.0);
  zEditor4 = new DoubleEdit(0.0);
  
  rotate->addWidget(xEditor2, 1, 1);
  rotate->addWidget(yEditor2, 1, 2);
  rotate->addWidget(zEditor2, 1, 3);

  rotate->addWidget(xEditor4, 2, 1);
  rotate->addWidget(yEditor4, 2, 2);
  rotate->addWidget(zEditor4, 2, 3);
  rotate->addWidget(new QLabel(tr("angle")), 1, 0);
  rotate->addWidget(new QLabel(tr("center")), 2, 0);
                                  
  
  rotateBox->setLayout(rotate);

  QGroupBox* scaleBox = new QGroupBox(tr("scale")); 
  QGridLayout* scale = new QGridLayout;
  scale->addWidget(new QLabel(tr("x:")), 0, 1);
  scale->addWidget(new QLabel(tr("y:")), 0, 2);
  scale->addWidget(new QLabel(tr("z:")), 0, 3);
  
 
  xEditor3 = new DoubleEdit(1.0); 
  yEditor3 = new DoubleEdit(1.0);
  zEditor3 = new DoubleEdit(1.0); 
  
 
  scale->addWidget(xEditor3, 1, 1);
  scale->addWidget(yEditor3, 1, 2);
  scale->addWidget(zEditor3, 1, 3);

  scale->addWidget(new QLabel(tr("    ")), 1, 0);
  scale->addWidget(new QLabel(tr("     ")), 2, 0);
  scaleBox->setLayout(scale);

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
  
  transferLayout->addWidget(translateBox);
  transferLayout->addWidget(rotateBox);
  transferLayout->addWidget(scaleBox);
  transferLayout->addLayout(buttonLayout);
  transferBox->setLayout(transferLayout);
  
  
  QGroupBox* tagBox = new QGroupBox(tr("tag")); 
  QHBoxLayout* tagLayout = new QHBoxLayout;
 
  tagEditor = new QLineEdit;
  tagLayout->addWidget(new QLabel(tr("   ")));
  tagLayout->addWidget(tagEditor);
  tagBox->setLayout(tagLayout);
  
  connect(xEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(xEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo())); 
  connect(xEditor3, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor3, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zEditor3, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(xEditor4, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor4, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zEditor4, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(tagEditor, SIGNAL(textChanged(const QString&)), this, SLOT(setInfo()));

  modBoundaries = new QStandardItemModel(bdnames.size(), 4, this);
  modBoundaries->setHeaderData(0, Qt::Horizontal, QObject::tr("color"));
  modBoundaries->setHeaderData(1, Qt::Horizontal, QObject::tr("boundary name"));
  modBoundaries->setHeaderData(2, Qt::Horizontal, QObject::tr("show/hide"));
  modBoundaries->setHeaderData(3, Qt::Horizontal, QObject::tr("new boundary name"));
 
  for (int i = 0; i < bdnames.size(); ++i) {
    QColor newColor = default_color[i%12];
    QStandardItem* colorItem = new QStandardItem("");
    QStandardItem* nameItem = new QStandardItem(bdnames[i].first);
    QStandardItem* showItem = new QStandardItem("show");
    QStandardItem* newNameItem = new QStandardItem("");
    colorItem->setBackground(QBrush(newColor));
      
    nameItem->setFlags(Qt::ItemIsSelectable | 
                       Qt::ItemIsUserCheckable | 
                       Qt::ItemIsEnabled);
      
    modBoundaries->setItem(i, 0, colorItem);
    modBoundaries->setItem(i, 1, nameItem);
    modBoundaries->setItem(i, 2, showItem);
    modBoundaries->setItem(i, 3, newNameItem);
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
          this, SLOT(setCurrentObj(QModelIndex)));

  QVBoxLayout *leftLayout = new QVBoxLayout;
  leftLayout->addWidget(tagBox);
  leftLayout->addWidget(boundaryView);
  
  QHBoxLayout *mainLayout =  new QHBoxLayout;
 
  mainLayout->addWidget(transferBox);
  mainLayout->addLayout(leftLayout);

  setLayout(mainLayout);
  //  setInfo();
}
void VMOption::accept(){
  tc.push_back(TranCoef());

  currentCoef = tc.size()-1;
  xEditor1->setText("0.0");
  yEditor1->setText("0.0");
  zEditor1->setText("0.0");

  xEditor2->setText("0.0");
  yEditor2->setText("0.0");
  zEditor2->setText("0.0");
  
  xEditor3->setText("1.0");
  yEditor3->setText("1.0");
  zEditor3->setText("1.0");
  
  xEditor4->setText("0.0");
  yEditor4->setText("0.0");
  zEditor4->setText("0.0");

 
}

void VMOption::cancel(){
 
 
  
  xEditor1->setText("0.0");
  yEditor1->setText("0.0");
  zEditor1->setText("0.0");

  xEditor2->setText("0.0");
  yEditor2->setText("0.0");
  zEditor2->setText("0.0");
  
  xEditor3->setText("1.0");
  yEditor3->setText("1.0");
  zEditor3->setText("1.0");
  
  xEditor4->setText("0.0");
  yEditor4->setText("0.0");
  zEditor4->setText("0.0");


}

void VMOption::clear(){
  tc.clear();
  tc.push_back(TranCoef());
  currentCoef = 0;
  xEditor1->setText("0.0");
  yEditor1->setText("0.0");
  zEditor1->setText("0.0");

  xEditor2->setText("0.0");
  yEditor2->setText("0.0");
  zEditor2->setText("0.0");
  
  xEditor3->setText("1.0");
  yEditor3->setText("1.0");
  zEditor3->setText("1.0");
  
  xEditor4->setText("0.0");
  yEditor4->setText("0.0");
  zEditor4->setText("0.0");

}
  
  
void VMOption::previous(){
  if(currentCoef == 0) return;
 
  currentCoef--;

  vector3d<double> translate = tc[currentCoef].translate;
  vector3d<double> rotateAngle = tc[currentCoef].rotateAngle;
  vector3d<double> scale = tc[currentCoef].scale;
  vector3d<double> rotateCenter = tc[currentCoef].rotateCenter;

  xEditor1->setValue(translate.x);
  yEditor1->setValue(translate.y);
  zEditor1->setValue(translate.z);
  
  xEditor2->setValue(-rotateAngle.x);
  yEditor2->setValue(-rotateAngle.y);
  zEditor2->setValue(-rotateAngle.z);
  
  xEditor3->setValue(scale.x);
  yEditor3->setValue(scale.y);
  zEditor3->setValue(scale.z);
  
  xEditor4->setValue(-rotateCenter.x);
  yEditor4->setValue(-rotateCenter.y);
  zEditor4->setValue(-rotateCenter.z);
  
}


void VMOption::next(){
 if(currentCoef >=(int)(tc.size()-1))return;
 currentCoef++;
  
 vector3d<double> translate = tc[currentCoef].translate;
 vector3d<double> rotateAngle = tc[currentCoef].rotateAngle;
 vector3d<double> scale = tc[currentCoef].scale;
 vector3d<double> rotateCenter = tc[currentCoef].rotateCenter;
 
  
  
 xEditor1->setValue(translate.x);
 yEditor1->setValue(translate.y);
 zEditor1->setValue(translate.z);
 
 xEditor2->setValue(-rotateAngle.x);
 yEditor2->setValue(-rotateAngle.y);
 zEditor2->setValue(-rotateAngle.z);
 
 xEditor3->setValue(scale.x);
 yEditor3->setValue(scale.y);
 zEditor3->setValue(scale.z);
 
 xEditor4->setValue(-rotateCenter.x);
 yEditor4->setValue(-rotateCenter.y);
 zEditor4->setValue(-rotateCenter.z);
  

}
  

void VMOption::showBoundary(QModelIndex top, QModelIndex ){

  if(top.column() ==2){//visibility item
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
 
  vector3d<double> translate =  vector3d<double>(xEditor1->value(),
                                                 yEditor1->value(), zEditor1->value());
 
  vector3d<double> rotateAngle = vector3d<double>(-xEditor2->value(),
                                                  -yEditor2->value(), -zEditor2->value());

 
  vector3d<double> rotateCenter = vector3d<double>(-xEditor4->value(),
                                                   -yEditor4->value(), -zEditor4->value());

 
  vector3d<double> scale = vector3d<double>(xEditor3->value(),
                                            yEditor3->value(), zEditor3->value());

  
 
  tc[currentCoef].rotateAngle = rotateAngle;
  tc[currentCoef].rotateCenter = rotateCenter;
  tc[currentCoef].translate = translate;
  tc[currentCoef].scale = scale;


  
  
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
    if(norm(tc[i].scale)!=1) gridXform.scale(tc[i].scale);
  }
 
  
  return gridXform;
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
  }
  if(!tagEditor->text().isEmpty()) text += " -tag " + tagEditor->text();
  for(int i =0; i < bdnames.size(); i++){
    if(!bdnames[i].second.isEmpty()){
      text += " -bc " +bdnames[i].first + "," + bdnames[i].second;
    }
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
  QPushButton *clearButton = new QPushButton(tr("clear all" ));
  
  
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
  
  //  QScrollArea* centralScrollArea = new QScrollArea;
//    centralScrollArea->setBackgroundRole(QPalette::Dark); 
   
   QGroupBox* central = new QGroupBox;
  central->setFlat(true);
  //QHBoxLayout *mainLayout = new QHBoxLayout;
 
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
 
  //mainLayout->addLayout(objLayout);
  // mainLayout->addWidget(mgviewer);
   
  
  // mainLayout->addStretch(1);
  central->setLayout(objLayout);

  
  // centralScrollArea->setWidget(central);
  

 
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
  // setCentralWidget(centralScrollArea);
    setCentralWidget(central);
  setWindowTitle(tr("vogmerge window"));
  setMinimumSize(1000, 700);
 
}

void VMergeWindow::createToolBar()
{

  // int spacing =2;  
  QToolBar* toolbar = addToolBar(tr("tree&vis"));
  
   addToolBarBreak();
  QGroupBox* visbar = new QGroupBox;
  visbar->setFlat(true);
  QHBoxLayout* visLayout = new QHBoxLayout;
  visLayout->addSpacing(800);
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

   toolbar->addWidget(visbar);

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
  //QListWidgetItem *bdCndiButton = new QListWidgetItem(typesWidget);
  //bdCndiButton->setText(filename);
  //bdCndiButton->setTextAlignment(Qt::AlignHCenter);
  //bdCndiButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
  
  QStringList bdnames;
  for(int i = 1; i < names.size(); i++)bdnames<<names.at(i);
  typesWidget->addItem(filename);
  int id = typesWidget->count()-1;
  VMOption* agrid = new VMOption(id, filename, bdnames);
  connect(agrid, SIGNAL(setCurrentColor(const IDColor&)), this, SIGNAL(setCurrentColor(const IDColor&)));
  connect(agrid, SIGNAL(setCurrentVisibility(const IDVisibility&)), this, SIGNAL(setCurrentVisibility(const IDVisibility&)));
  connect(agrid, SIGNAL(tcChanged(const IDMatrix&)), this, SIGNAL(tcChanged(const IDMatrix&)));
  pagesWidget->addWidget(agrid);
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
    connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus)), this, SLOT(afterMerge(QString, QProcess::ExitStatus)));
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

void VMergeWindow::afterMerge(QString command, QProcess::ExitStatus status){
  if(status==QProcess::NormalExit){
    clear();
    QString filename = command.section(' ', -1, -1);
    // emit getGrid(filename);
     mgviewer->get_boundary(filename);
  }
}

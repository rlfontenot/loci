#include <QtGui>
#include <QApplication>
#include <QToolBar>
#include <QAction>
#include <QFileDialog>
#include <QSlider>
#include <QStandardItemModel>
#include <QTableView>
#include <QDockWidget>
#include <QTextDocument>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMenuBar>
#include <QStackedWidget>
#include <QDomElement>
#include <QApplication>
#include "mainwindow.h"
#include "glviewer.h"
#include "cutdialog.h"
#include "bdcndwindow.h"
#include "solverwindow.h"
#include "initcndwindow.h"
#include "importwindow.h"
#include "physicswindow.h"
#include "mgviewer.h"
#include "refdialog.h"
#include "qualitydialog.h"
#include "progressdialog.h"
#include <cstdlib>
#include <QString>
#include <string>
#include <utility>
#include <algorithm>

using std::string;
using std::vector;
using std::sort;
using std::pair;



void MainWindow::createMenu(){
  tb = new QToolBar(tr("tool bar"),this);
  tb->setWindowTitle(tr("File Actions"));
  addToolBar(tb);
  
  
  QAction *openScaAct = new QAction(tr("&scalar file"), this);
  connect(openScaAct, SIGNAL(triggered()), this, SLOT(openSca()));
  


  QAction *openVogAct = new QAction(tr("vo&g file"), this);
  connect(openVogAct, SIGNAL(triggered()), this, SLOT(openVog()));
  
  QAction *newCaseAct = new QAction(QIcon( pngpath+ "filenew.png"), tr("&New case"), this);
  newCaseAct->setShortcut(QKeySequence::New);
  connect(newCaseAct, SIGNAL(triggered()), this, SLOT(newCase()));
 
  QAction *openCaseAct = new QAction(QIcon( pngpath+ "fileopen.png"), tr("&Open case"), this);
  newCaseAct->setShortcut(QKeySequence::Open);
  connect(openCaseAct, SIGNAL(triggered()), this, SLOT(openCase()));
  
      
 

  QAction *saveXmlAct = new QAction( QIcon( pngpath+"save.png"), tr("save &Case"), this);
  connect(saveXmlAct, SIGNAL(triggered()), this, SLOT(saveXml()));
  
  QAction *saveImageAct = new QAction(  QIcon(  pngpath+"save.png"), tr("&image file"), this);
  connect(saveImageAct, SIGNAL(triggered()), this, SLOT(snapshot()));
  
  
  QAction *exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  connect(exitAct, SIGNAL(triggered()), qApp, SLOT(quit()));
  
  QAction *aboutPreprocessAct = new QAction(tr("P&reprocess"), this);
  connect(aboutPreprocessAct, SIGNAL(triggered()), this, SLOT(aboutPreprocess()));

  QAction *aboutPostprocessAct = new QAction(tr("P&ostprocess"), this);
  connect(aboutPostprocessAct, SIGNAL(triggered()), this, SLOT(aboutPostprocess()));
  
  QMenu* loadMenu = new QMenu(tr("Load grid"), this);
  // loadMenu->addAction(openTopoAct);
  //loadMenu->addAction(openScaAct);
  loadMenu->addAction(openVogAct);
  
  QMenu* saveAsMenu = new QMenu(tr("&Save"), this);
 
  saveAsMenu->addAction(saveXmlAct);
  saveAsMenu->addAction(saveImageAct);
  
  
  QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newCaseAct);
  fileMenu->addAction(openCaseAct);
  fileMenu->addMenu(loadMenu);
  fileMenu->addMenu(saveAsMenu);
  fileMenu->addAction(exitAct);
 
  menuBar()->addSeparator();
  QMenu* helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutPreprocessAct);
  helpMenu->addAction(aboutPostprocessAct);

  menuBar()->addSeparator();
  viewMenu = menuBar()->addMenu(tr("&View"));
 
  tb->addAction(newCaseAct);
  tb->addAction(openCaseAct);
  tb->addAction(saveXmlAct);

  QAction* showStatusAct = new QAction(tr("Show Status"),this);
  connect(showStatusAct, SIGNAL(triggered()), this, SLOT(toggleShowStatus()));
  tb->addSeparator();
  tb->addAction(showStatusAct);   

  menuBar()->addSeparator();
  QAction *quitAct = new QAction(tr("&Quit"), this);
  menuBar()->addAction(quitAct);

 connect(quitAct, SIGNAL(triggered()),
         qApp, SLOT(quit()));  
 connect(openScaAct, SIGNAL(triggered()),
         this, SLOT(resetSlider()));

 viewMenu->addAction(tb->toggleViewAction()); 
}




void MainWindow::createVisBar(){
  
  visbar = new QToolBar("Visualization");

  addToolBar(Qt::TopToolBarArea,visbar );
  insertToolBarBreak(visbar);

  QAction *selectBoundaryAct = new QAction(tr("Select Boundary"), this);
  visbar->addAction(selectBoundaryAct);
  connect(selectBoundaryAct, SIGNAL(triggered()),
          this, SLOT(selectBoundaryPressed())); 
  visbar->addSeparator();

  QAction *clearBoundaryAct = new QAction(tr("Clear"), this);
  visbar->addAction(clearBoundaryAct);
  connect(clearBoundaryAct, SIGNAL(triggered()),
          this, SLOT(clearCurrent())); 
  visbar->addSeparator();
  
  QAction *showBoundariesAct = new QAction(tr("show Boundaries"), this);
  visbar->addAction(showBoundariesAct);
  connect(showBoundariesAct, SIGNAL(triggered()),
          viewer, SLOT(showBoundaries())); 
  visbar->addSeparator();
  
  QAction* resetAct = new QAction(tr("Reset"), this);
  visbar->addAction(resetAct);
  connect(resetAct, SIGNAL(triggered()),
          this, SLOT(reset()));
  visbar->addSeparator();
  QAction* fitAct = new QAction(tr("Fit"), this);
  visbar->addAction(fitAct);
  connect(fitAct, SIGNAL(triggered()),
          this, SLOT(fit()));


  visbar->addSeparator();
  QAction* snapshotAct = new QAction(tr("Snapshot"), this);
  visbar->addAction(snapshotAct);
  connect(snapshotAct, SIGNAL(triggered()),
          this, SLOT(snapshot()));
  
  viewMenu->addAction(visbar->toggleViewAction());
  tb->addSeparator();
  tb->addAction(visbar->toggleViewAction());
 
  viewerAct = new QAction(tr("Viewer"), this);
  connect(viewerAct, SIGNAL(triggered()),
          this, SLOT(toggleViewer()));
  tb->addSeparator();
   tb->addAction(viewerAct);
 
}
void MainWindow::reset(){
  if(central->currentWidget()==viewer)viewer->reset();
  else  if(central->currentWidget()==mgviewer) mgviewer->reset();
}
void MainWindow::clearCurrent(){
  if(central->currentWidget()==viewer)viewer->clearCurrent();
  else  if(central->currentWidget()==mgviewer) mgviewer->clearCurrent();
}

void MainWindow::fit(){
  if(central->currentWidget()==viewer)viewer->fit();
  else  if(central->currentWidget()==mgviewer) mgviewer->fit();
}

  
void MainWindow::toggleViewer(){
  if(central->currentWidget()==viewer && previousWidget!=0){
    central->setCurrentWidget(previousWidget);
    bdock->hide();
  }else{
    previousWidget = central->currentWidget();
    central->setCurrentWidget(viewer);
    bdock->show();
  }
}

void MainWindow::toggleShowStatus(){
  displayStatus = !displayStatus;
  emit showStatus(displayStatus);
}

void MainWindow::selectBoundaryPressed(){

  if(boundaryView == 0){
    bool success = selectBoundary();
    if(!success) return;
  }
  //boundaryView->show();
 
  if(central->currentWidget()==viewer && previousWidget!=0){
  }else{
    previousWidget = central->currentWidget();
    central->setCurrentWidget(viewer);
  }
  
 
}
    
void MainWindow::updateConditionView(){
  QModelIndex index = boundaryView->currentIndex();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  QDomElement myroot = doc.documentElement();
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
}  
void MainWindow::createFlowBar(){
  
  if(flowbar){
    delete flowbar;
    flowbar = 0;
    
  }
  if(flowbarButtons){
    delete flowbarButtons;
    flowbarButtons = 0;
  }
  flowbar = new QToolBar;
  flowbarButtons = new QButtonGroup(this);
 
  addToolBar(Qt::RightToolBarArea,flowbar );
  QDomElement theroot = doc.documentElement();
  QDomElement theelem = theroot.firstChildElement("mainWindow");

  QDomElement elem = theelem.firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                           elem.tagName()+ tr(" has no child")
                         );
    return;
  }

  //add hard-coded buttons
 

  flowbar->addSeparator();
  QPushButton* adaptButton = new QPushButton(tr("FVMAdapt"), this);
  flowbar->addWidget(adaptButton);
  connect(adaptButton, SIGNAL(clicked()), this, SLOT(adaptClicked()));


 flowbar->addSeparator();
  QPushButton* vmButton = new QPushButton(tr("Vogmerge"), this);
  flowbar->addWidget(vmButton);
  connect(vmButton, SIGNAL(clicked()), this, SLOT(vmClicked()));

  flowbar->addSeparator();
  QPushButton* vcheckButton = new QPushButton(tr("Vogcheck"), this);
  flowbar->addWidget(vcheckButton);
  connect(vcheckButton, SIGNAL(clicked()), this, SLOT(vcheck()));
  
  
  int count=0;
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) { 
    flowbar->addSeparator();
    
    QPushButton* newButton = new QPushButton(elem.attribute("buttonTitle"), this);
    flowbar->addWidget(newButton);
    flowbarButtons->addButton(newButton, count);
    elem.setAttribute("buttonIndex", count);
    newButton->setToolTip(elem.attribute("toolTip"));
    newButton->setWhatsThis(elem.attribute("whatsThis"));
    newButton->setStatusTip(elem.attribute("status"));
 
    QWidget* newWindow=0;
    
    if(elem.attribute("element")=="importWindow"){
      newWindow=new ImportWindow(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
    }else  if(elem.attribute("element")=="solverWindow"){
      newWindow=new SolverWindow(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
      
    }else  if(elem.attribute("element")=="physicsWindow"){
      newWindow=new PhysicsWindow(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(newWindow, SIGNAL(stateChanged()), this, SIGNAL(stateChanged()));
      connect(newWindow, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
       connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    }else  if(elem.attribute("element")=="initialWindow"){
      
      initWindow=new InitCndWindow(elem, theroot);
      connect(initWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(initWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), initWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), initWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), initWindow, SLOT(updateShowStatus(const bool&)));
      connect(initWindow, SIGNAL( valueChanged(const QTreeWidgetItem*)),
              viewer, SLOT(updateDoc(const QTreeWidgetItem*)));
    }else  if(elem.attribute("element")=="panel"){
      newWindow=new OptionPage(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
       connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    } else  if(elem.attribute("element")=="page"){
      newWindow=new Page(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    } 
      
    if(newWindow && elem.attribute("inDock")!="true"){
      int i = central->addWidget(newWindow);
      elem.setAttribute("widgetIndex", i);
    }
  }
  connect(flowbarButtons, SIGNAL(buttonClicked(int)), this, SLOT(changePage(int)));

  flowbar->addSeparator();
  QPushButton* ppsButton = new QPushButton(tr("Post\nProcessing"), this);
  flowbar->addWidget(ppsButton);
  connect(ppsButton, SIGNAL(clicked()), this, SLOT(cut()));
}
void MainWindow::changePage(int index){

  
  // dock->setWidget(statusEdit);
 
  if(mgviewer){
    central->removeWidget(mgviewer);
    delete mgviewer;
    mgviewer = 0;
  }
  
  if(vmwindow){
    delete vmwindow;
    vmwindow = 0;
  }

  
  QDomElement theroot = doc.documentElement();
  QDomElement elem = theroot.firstChildElement("mainWindow");
  elem=elem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), "main xml file",
                         tr(" mainWindow has no child" ));
    return;
  }
  for(int i=0; i< index; i++)elem=elem.nextSiblingElement();

  if(elem.tagName()=="saveVar"){
    saveVar();
    return;
  }
  
  if(elem.tagName()=="gridSetup"){
    setGrid(elem);
    return;
  }
  if(elem.attribute("element")=="boundaryWindow"){
    setBoundary(elem);
    
  }
  
  if(elem.hasAttribute("widgetIndex")) 
    central->setCurrentIndex(elem.attribute("widgetIndex").toInt());
  else  central->setCurrentIndex(0);
  if(elem.attribute("inDock")=="true"&&elem.attribute("element")=="boundaryWindow"){
    dock->setWidget(bdWindow);
    dock->show();
    dock->setFloating(true);
  }else if(elem.attribute("inDock")=="true"&&elem.attribute("element")=="initialWindow"){
    dock->setWidget(initWindow);
    dock->show();
    dock->setFloating(true);
  } else{
    dock->setWidget(statusWindow);
    dock->setFloating(false);
  }
  if(central->currentWidget() == viewer){
    bdock->show();
  }else{
    bdock->hide();
  }
   
}


  

void MainWindow::createDockWindow(){
  //status window
  statusEdit = new QTextEdit(this);
  statusWindow = new QGroupBox(this);
  statusWindow->setFlat(true);
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  QPushButton* clearAllButton = new QPushButton(tr("Clear All"));
  QPushButton* clearLastButton = new QPushButton(tr("Clear Last"));
  buttonLayout->addWidget(clearAllButton);
  buttonLayout->addWidget(clearLastButton);
  QVBoxLayout* mainLayout=new QVBoxLayout;
  
  
  mainLayout->addWidget(statusEdit);
  mainLayout->addLayout(buttonLayout);
  statusWindow->setLayout(mainLayout);
  connect(clearAllButton, SIGNAL(clicked()), this, SLOT(clearAllStatus()));
  connect(clearLastButton, SIGNAL(clicked()), this, SLOT(clearLastStatus()));

  //boundary dock
  bdock  = new QDockWidget("boundary dock window", this);
  bdock->setAllowedAreas(Qt::TopDockWidgetArea );

  if(boundaryView)bdock->setWidget(boundaryView);
  addDockWidget(Qt::TopDockWidgetArea, bdock);
  viewMenu->addAction(bdock->toggleViewAction());
  tb->addSeparator();
  tb->addAction(bdock->toggleViewAction());   
  //dock  
  dock = new QDockWidget("status/boundary-condition dock window", this);
  dock->setAllowedAreas(Qt::LeftDockWidgetArea );
  dock->setWidget(statusWindow);
  addDockWidget(Qt::LeftDockWidgetArea, dock);
  viewMenu->addAction(dock->toggleViewAction());
  tb->addSeparator();
  tb->addAction(dock->toggleViewAction());   
}

void MainWindow::clearAllStatus(){
  statusEdit->clear();
}
void MainWindow::clearLastStatus(){
  statusEdit->undo();
}
void MainWindow::updateStatus(const QString& s){
  if (s.isEmpty())
    return;
  QTextCursor cursor(statusEdit->textCursor());
  if (cursor.isNull())
    return;
  cursor.beginEditBlock();
  cursor.insertBlock();
  cursor.insertText(s);
  cursor.insertBlock();
  cursor.endEditBlock();
}

void MainWindow::updateStatusTip(int button){
  QDomElement theroot = doc.documentElement();
  theroot.setAttribute("updated", "true");
  QDomElement theelem = theroot.firstChildElement("mainWindow");

  QDomElement elem = theelem.firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         elem.tagName()+ tr(" has no child")
                         );
    return;
  }
  
  
  for (int i=0; i < button; i++)elem = elem.nextSiblingElement();
  if(elem.attribute("buttonIndex").toInt() != button){
    
   QMessageBox::warning(window(), ".xml",
                        QString("can not find button: %L1 ").arg(button)
                         );
   return;
  }
  if(elem.hasAttribute("status"))qobject_cast<QWidget*>(flowbarButtons->button(button))->setStatusTip(elem.attribute("status"));
}
MainWindow::MainWindow()
{

  char* resourcepath = getenv("CHEMDEMOPATH");
  xmlpath = "./xml/";
  pngpath = "./png/";
  
  if(resourcepath){
    xmlpath = QString(resourcepath).append("xml/");
    pngpath = QString(resourcepath).append( "png/");
  }
  
  //  setAttribute(Qt::WA_DeleteOnClose); //will cause SF
  QString filename= xmlpath+"main.xml";
  QFile file(filename);
  if(!file.exists()){
    QMessageBox::information(window(), filename,
                             filename+ tr(" doesn't exist"));
    return;
  }
  if (!file.open(QIODevice::ReadOnly)){

    QMessageBox::information(window(), "main.xml",
                             tr("cannot open ")+ filename + tr("for reading"));
    
    return;
  }
  
  QString errorStr;
  int errorLine;
  int errorColumn;
  if (!doc.setContent(&file, true, &errorStr, &errorLine,
                      &errorColumn)) {
    QMessageBox::information(window(), filename,
                             tr("Parse error at line %1, column %2:\n%3")
                             .arg(errorLine)
                             .arg(errorColumn)
                             .arg(errorStr));
    file.close();
    return;
  }
  
  
  file.close();
   
   
  isNewCase= true;
 

  viewer = new GLViewer(this);
  central = new QStackedWidget(this);
  central->addWidget(viewer);
  
  modBoundaries = 0;
  boundaryView = 0;
  bdock = 0;
  bdWindow = 0;
  initWindow = 0;
  cutdialog = 0;

  dock = 0;
  adaptwindow = 0;
  refdialog = 0;
  

  flowbar =0;
  flowbarButtons = 0;
 
  visbar =0;
  bdButtonDown = false;
  setCentralWidget(central);
  previousWidget = 0;
  statusWindow = 0;
  displayStatus = false;
  mgviewer = 0;
  vmwindow = 0;
  waitForQualityFile = false;
  createMenu();
  createVisBar();
  //createDisplayBar();
  
 
  createDockWindow();
  createFlowBar();
  hideVisBar();
 
    
  setWindowTitle(tr("chem demo"));
  updateStatus(tr("Please use 'Grid Setup' to load  grid information, then start a new case or use file menu to open a case"));
  
  statusBar()->showMessage(tr("Ready"));
}
////////////////////////////////////////
//  public:
//    QSize sizeHint() const;
//
//  Requests a default size of 800x600.
////////////////////////////////////////

QSize MainWindow::sizeHint() const
{
  return QSize(800, 800);
}



void MainWindow::snapshot(){
  
  QImage pic ;
  if(central->currentWidget()==viewer)pic = viewer->grabFrameBuffer(true);
  else  if(central->currentWidget()==mgviewer) pic = mgviewer->grabFrameBuffer(true);

 

  QString format = "png";
  QString initialPath = QDir::currentPath() + tr("/untitled.") + format;

  QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"),
                                                  initialPath,
                                                  tr("%1 Files (*.%2);;All Files (*)")
                                                  .arg(format.toUpper())
                                                  .arg(format));

  if (!fileName.isEmpty()){
    bool filesaved =  pic.save(fileName, format.toAscii());
    if(filesaved) updateStatus(fileName+tr(" saved"));
    else {
      QMessageBox::information(window(), "mainwindow",
                               tr("Saving ")+ fileName + tr(" failed"));
      
      return;
    }
  }else{
    QMessageBox::information(window(), "mainwindow",
                             tr("Please specify filename for saving" ));
    return;
  }
}


  
void MainWindow::setGrid(QDomElement& theelem)
{
  
  if(cutdialog){
    delete cutdialog;
    cutdialog = 0;
  }

    
  //set up the default filename

  QString format = "volume grid file(*.vog)";
  if(theelem.hasAttribute("format")) format = theelem.attribute("format");
  QString initialPath =QDir::currentPath();
  if(theelem.hasAttribute("initialPath"))initialPath = theelem.attribute("initialPath");
  QString fileName  = initialPath;
  if(theelem.hasAttribute("casename")) fileName =  theelem.attribute("directory")+"/"+theelem.attribute("casename")+".vog";

  //select a file
  fileName = QFileDialog::getOpenFileName(this, tr("Load Grid"),
                                          fileName,
                                          format);
  
 
  if(fileName==""){
    //no error message in  case of 'cancel' is pressed 
    
    return;
  }
  
  
  //load in grid
 
     
    QString surfFileName = fileName.section('.', 0, -2)+".surface";
  
 

   
   
     QString caseName =  (fileName.section('.', 0, -2)).section('/', -1);
    
    QFileInfo surfInfo(surfFileName);
    QFileInfo vogInfo(fileName);

     QString command2 = "vog2surf -surface " + surfFileName+" " + fileName.section('.', 0, -2);
    if(!(surfInfo.exists()) || surfInfo.created() < vogInfo.created()){
      //int first= fileName.lastIndexOf('/');
      //int last = fileName.lastIndexOf('.');
      //QString casename = fileName.mid(first+1, last-first-1);
      //QString directory = fileName.left(first);
      
     //  QString script_filename = directory + "/output/vog2surf_"+caseName;
//       QString out_filename= directory + "/output/vog2surf_"+caseName+".out";
      
//       QFile outfile(script_filename);
//       if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
//         QMessageBox::information(window(), "mainwindow",
//                                  tr("Cannot open ") + script_filename + tr(" for writing!"));
//         return;
//       }
      
//        QTextStream out(&outfile);
    
       //replace this command later, no directory, and vog2surf -surface will be used 
      // QString command2 = "./vog2surf -surface " + surfFileName+" " + fileName.section('.', 0, -2);
      
     //   out <<"#!/bin/bash"<<endl;
//        out <<"exec 6>&1"<<endl;
//        out <<"exec 7>&2"<<endl;
//        out<< "exec &> "<< out_filename <<endl;
//        out << command2 << endl;
//        out<<"exec 1>&6 6>&- " << endl;
//        out<<"exec 2>&7 7>&- " << endl;

//       outfile.close();
//       QString command3 = "chmod 777 " + script_filename;
      
//       int ret =  system(command3.toStdString().c_str());



//       if(!WIFEXITED(ret))
//         {
//           if(WIFSIGNALED(ret))
//             {
//               QMessageBox::information(window(), "mainwindow",
//                                        command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
//                theelem.removeAttribute("casename");
//               return;
//             }
//         }

//        emit updateStatus("COMMAND:    " + command2+ "     STARTED");
//        ret = system(script_filename.toStdString().c_str());
//        if(!WIFEXITED(ret)){
//          if(WIFSIGNALED(ret)){
//            QMessageBox::information(window(), "mainwindow",
//                                     script_filename + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
//            theelem.removeAttribute("casename");
//            return;
//          }
//        }
//        emit updateStatus("COMMAND:    " + command2+ "     FINISHED");


      
   //    QFile file(out_filename);
//       if (!file.open(QFile::ReadOnly | QFile::Text)) {

//         QMessageBox::information(window(), "mainwindow",
//                                  tr("Cannot open ") + out_filename + tr(" for reading!"));

//         theelem.removeAttribute("casename");
        
//         return;
//       }
    
//       QTextStream in(&file);
//       QApplication::setOverrideCursor(Qt::WaitCursor);
     
     //  emit updateStatus(in.readAll());
//       QApplication::restoreOverrideCursor();
//       file.close();
       
      ProgressDialog* progress = new ProgressDialog(command2, true);
       progress->show();
       connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus)), this, SLOT(loadGrid(QString, QProcess::ExitStatus)));
       
    }else{
      loadGrid(command2, QProcess::NormalExit);
    }
  
}
void MainWindow::loadGrid(QString command, QProcess::ExitStatus ){
  QDomElement theroot = doc.documentElement();
  QDomElement theelem = theroot.firstChildElement("mainWindow");
  theelem=theelem.firstChildElement("gridSetup");
  if(theelem.isNull()){
    QMessageBox::warning(window(), "main xml file",
                         tr(" mainWindow has no child gridSetup" ));
    return;
  }
  


  
  QString surfFileName = command.section(' ', 2, 2);
  QString fileName =   command.section(' ', -1, -1)+".vog";
  
   QString caseName =  (fileName.section('.', 0, -2)).section('/', -1);
  // must setCurrentWidget(viewer)first, then load_boundary
  // must leave these two lines at the end of function
  
   
  QStringList boundary_names;   
  central->setCurrentWidget(viewer);   
  bool loaded =   viewer->load_boundary(surfFileName, boundary_names); // and setup the GLWidget.
  viewer->show();
  viewer->reset();
  

  
  
  if(loaded){
    //if different case, remind the user to save the case
    if(theelem.hasAttribute("casename") && theelem.attribute("casename")!=caseName) {
      int button = QMessageBox::question(this, tr("a new grid is loaded"),
                                         tr("Do you want to save the old case ?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveXml();
      
      button = QMessageBox::question(this, tr("a new grid is loaded"),
                                         tr("Do you want to save the .vars file of the old case ?"),
                                     QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveVar();

      //re-create the boundary_condition node
      QDomElement root= doc.documentElement();
      QDomNode cndNode = root.firstChildElement("boundary_conditions");
      if(cndNode.isNull()){
        QMessageBox::warning(window(), ".xml",
                           tr("can not find element 'boundary_conditions'")
                             );
        return;
      }
      
      QDomNode tmp = root.removeChild(cndNode);
      if(tmp.isNull()){
        QMessageBox::warning(window(), "main window",
                             tr("remove boundary conditions node failed")
                             );
        return;
      }
      QDomNode newNode = doc.createElement("boundary_conditions");
      tmp = root.appendChild(newNode);
      if(tmp.isNull()){
        QMessageBox::warning(window(), "main window",
                           tr("append boundary conditions node failed")
                             );
        return;
      }

    
      
    }

    //set up the new case
    
    int first= fileName.lastIndexOf('/');
    int last = fileName.lastIndexOf('.');
    QString casename = fileName.mid(first+1, last-first-1);
    QString directory = fileName.left(first);
    theelem.setAttribute("casename", casename);
    theelem.setAttribute("directory", directory);
    theelem.setAttribute("boundary_names", boundary_names.join(","));
      
    updateStatus(fileName + tr(" loaded"));
    theelem.setAttribute("status", "done");

    //clean up the data
       if(modBoundaries){
      delete modBoundaries;
      modBoundaries = 0;
    }
    
    if(boundaryView){
      delete boundaryView;
      boundaryView = 0;
    }
    
    
    
    if(bdWindow){
      delete bdWindow;
      bdWindow=0;
    }
  // update status
    updateStatusTip(theelem.attribute("buttonIndex").toInt());
  

    if(boundaryView == 0){
      selectBoundary();
      updateConditionView();
    }
    bdock->show();
   
    
    
  }else{
    // theelem.removeAttribute("casename");
    updateStatus(fileName + tr(" not loaded"));
  }
  
  viewer->clearCurrent();
  
  
}

void MainWindow::openVog(){
  
}

void MainWindow::setBoundary(QDomElement& elem){
  bdButtonDown = true;
  dock->show();
  bdock->show();
  
  if(bdWindow ==0){
    QDomElement theroot = doc.documentElement();
    QDomElement root = theroot.firstChildElement("mainWindow");
    QDomElement elem_bdname = root.firstChildElement("gridSetup");
    if(elem_bdname.isNull()){
      QMessageBox::warning(window(), tr("main, boundary.xml"),
                           tr("can not find gridSetup element")
                           );
      return;
    }
   
    if(!elem_bdname.hasAttribute("boundary_names")){
      
      QMessageBox::warning(window(), tr("main, boundary.xml"),
                           tr("no bounday names, please use Grid Setup to load grid first")
                           );
      return;
    }
    QStringList bdnames = elem_bdname.attribute("boundary_names").split(",");
    if(boundaryView==0)selectBoundary();
    else boundaryView->show();
   
    bdWindow = new BdCndWindow(elem, theroot, bdnames, boundaryView);
    connect(this, SIGNAL(setCurrent(QModelIndex)), bdWindow, SLOT(setCurrent(QModelIndex)));
    connect(bdWindow, SIGNAL(closed()), this, SLOT(bdWindowClosed())); 
    connect(bdWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
    connect(bdWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
    connect(this, SIGNAL(stateChanged()), bdWindow, SLOT(changeState()));
    connect(this, SIGNAL(componentsChanged()), bdWindow, SIGNAL(componentsChanged()));
    connect(this, SIGNAL(showStatus(const bool &)), bdWindow, SIGNAL(showStatus(const bool &)));
    connect( bdWindow, SIGNAL(updateConditionView()), this, SLOT(updateConditionView()));
    
    updateStatus(tr("start bounadry conditions setup ..."));
  }

  dock->setWidget(bdWindow);

 }

 void MainWindow::bdWindowClosed(){
  bdButtonDown = false;
  viewer->clearCurrent();
  dock->setWidget(statusWindow);
  bdock->hide();
  }








bool MainWindow::selectBoundary(){
  
    QDomElement theroot = doc.documentElement();
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

    
    
    theroot = doc.documentElement();
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
          this, SIGNAL(setCurrent(QModelIndex)));
  connect(boundaryView, SIGNAL(clicked( const QModelIndex&)),
          this, SLOT(setCurrentObj(QModelIndex)));

  connect(viewer, SIGNAL(pickCurrent(int)), this, SLOT(selectCurrent(int)));
  if(!bdButtonDown)boundaryView->show();

  theroot = doc.documentElement();
  theroot = theroot.firstChildElement("boundary_conditions");
  int currentIndex = theroot.attribute("currentIndex").toInt();
  selectCurrent(currentIndex);

  bdock->show();
  bdock->setWidget(boundaryView);
  return true;
}

// void MainWindow::setCondition(const QTextDocument& c){
//   if(boundaryView){
//     QModelIndex index;
//     if(boundaryView->currentIndex().column()==3){
//       boundaryView->model()->setData(boundaryView->currentIndex(), c);
//     }else{
//       index = boundaryView->currentIndex().sibling( boundaryView->currentIndex().row(), 3);
//       boundaryView->model()->setData(index, c);
//     }
//   }
// }



    
void MainWindow::showBoundary(QModelIndex top, QModelIndex ){
 
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


void MainWindow::showVisBar(){
 if(visbar) visbar->show();
}


void MainWindow::hideVisBar(){
 if(visbar) visbar->hide();
}


void MainWindow::setCurrentObj(QModelIndex top){
     QColor value= qobject_cast<const QStandardItemModel*>(top.model())->item(top.row(), 0)->background().color();
     viewer->setCurrentObj(top.row(), value);
 }

void MainWindow::selectCurrent(int row){
  
    QColor value= modBoundaries->item(row, 0)->background().color();
    QModelIndex index =qobject_cast<const QAbstractItemModel*>(modBoundaries)->index(row, 1);
    qobject_cast<QAbstractItemView*>( boundaryView)->setCurrentIndex(index);
    emit setCurrent(index);
    
    viewer->setCurrentObj(row, value);
  
 }

void MainWindow::cut(){
  

   
 
  if(cutdialog) {
    delete cutdialog;
    cutdialog = 0;
  }
  
  if(viewer ==0) return;
  QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
  elem = elem.firstChildElement("gridSetup");
  LoadInfo ldinfo;
  ldinfo.casename = elem.attribute("casename");
  ldinfo.directory = elem.attribute("directory");
  if(ldinfo.casename.isEmpty()){
    QMessageBox::warning(window(), "post-processing",
                         tr("No casename, please load in grid with 'Grid Setup'")
                         );
    return;
  }

  QDir dir(ldinfo.directory+"/output/");
  QStringList filters;
  filters << "grid_pos.*_" + ldinfo.casename;
  QStringList gridposFiles = dir.entryList(filters);
  if(gridposFiles.size()==0){
   
    int ret = QMessageBox::question(this, "post-processing",
                                    tr("No scalar value, do you want to run vogcheck? "),
                                    QMessageBox::Ok | QMessageBox::Cancel);
    switch(ret){
    case QMessageBox::Ok:
      waitForQualityFile=true;
      check(ldinfo.directory+'/'+ldinfo.casename);
      break;
    default:
      return;
    }
  }else{
    cutdialog = new CutDialog(ldinfo, viewer->boundaryBoxSize(), viewer);
    cutdialog->show();
  }
}

void MainWindow::vcheck(){
  QString fileName =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                  QDir::currentPath(),
                                 tr("vog Files (*.vog)"));
  
  fileName = fileName.section('.', 0, -2);
  check(fileName);
}






  
void MainWindow::check(const QString& fn){

  QString importFileName = fn;
  QString casename = importFileName.section('/', -1, -1);
  
  QFile exist_test(importFileName+tr(".vog"));
  if(!(exist_test.exists())){
    QMessageBox::warning(window(), tr("vogcheck"),
                         tr("Please convert the file to volume grid format first")
                         );
    return;
  }
  
  //  QString script_filename = "./output/check_"+casename;
  QString out_filename="./output/check_"+casename+".out";
  



//   QFile outfile(script_filename);
//   if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
//      QMessageBox::warning(this, tr("file io "),
//                             tr("Cannot write file %1:\n%2.")
//                             .arg(script_filename)
//                           .arg(outfile.errorString()));
//     return;
//   }
  
//   QTextStream out(&outfile);
  
  
   QString command2 = "vogcheck " +casename;
 
 
  
//   out <<"#!/bin/bash"<<endl;
//   out <<"exec 6>&1"<<endl;
//   out <<"exec 7>&2"<<endl;
//   out<< "exec &> "<< out_filename <<endl;
//   out<<"cd " << importFileName.section('/',0,  -2)<<endl;
//   out << command2 << endl;
//   out<<"exec 1>&6 6>&- " << endl;
//   out<<"exec 2>&7 7>&- " << endl;

//   outfile.close();
//   QString command3 = "chmod 777 " + script_filename;

 
 
//   int ret = system(command3.toStdString().c_str());

//   if(!WIFEXITED(ret))
//     {
//       if(WIFSIGNALED(ret))
//         {
//           QMessageBox::information(window(), "mainwindow",
//                                    command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
//           return;
//         }
//       exit(0);
//     }

//    emit updateStatus("COMMAND:    " + command2 + "    STARTED");
//   ret = system(script_filename.toStdString().c_str());

//   if(!WIFEXITED(ret))
//     {
//       if(WIFSIGNALED(ret))
//         {
//           QMessageBox::information(window(), "mainwindow",
//                                    command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
//           return;
//         }
//       exit(0);
//     }
//   emit updateStatus("COMMAND:    " + command2 + "    FINISHED");

   ProgressDialog* progress = new ProgressDialog(command2);
   connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus)), this, SLOT(showQuality(QString, QProcess::ExitStatus)));
   progress->show();
  
   
   
 //   QFile file(out_filename);
//    if (!file.open(QFile::ReadOnly | QFile::Text)) {
//      QMessageBox::warning(this, tr(" file io "),
//                           tr("Cannot read file %1:\n%2.")
//                           .arg(out_filename)
//                           .arg(file.errorString()));
//     return;
//   }
  
//   QTextStream in(&file);
//   QApplication::setOverrideCursor(Qt::WaitCursor);
//   updateStatus(in.readAll());
//   QApplication::restoreOverrideCursor();
//   file.close();


  
 //  QString qualityFileName = importFileName+tr(".quality");
//   QualityDialog qualityDialog(qualityFileName, this);
//   qualityDialog.exec();
  
    

}
void MainWindow::showQuality(QString command, QProcess::ExitStatus status){
  if(status==QProcess::CrashExit)return;
  
  QString filename = command.section(' ',-1, -1)+".quality";
  QualityDialog qualityDialog(filename, this);
  qualityDialog.exec();

  if(waitForQualityFile){
    waitForQualityFile = false;
    QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
    elem = elem.firstChildElement("gridSetup");
    LoadInfo ldinfo;
    ldinfo.casename = elem.attribute("casename");
    ldinfo.directory = elem.attribute("directory");
    cutdialog = new CutDialog(ldinfo, viewer->boundaryBoxSize(), viewer);
    cutdialog->show();
    
  }
}
  
  
// void MainWindow::markVolumeNodes(){

//   QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
//   elem = elem.firstChildElement("gridSetup");
//   if(elem.attribute("directory")=="" || elem.attribute("casename")==""){
//     QMessageBox::information(window(), "mainwindow",
//                                 tr("Please use 'GridSetup' reading in grid first"));
//     return;
//   }
//   QString fileName = elem.attribute("directory")+"/"+elem.attribute("casename")+".vog";
//   viewer->markVolumeNodes(fileName);
// }

void MainWindow::refineGrids(){
  if(refdialog){
    delete refdialog;
    refdialog = 0;
  }
  
  QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
  elem = elem.firstChildElement("gridSetup");
  if(elem.attribute("directory")=="" || elem.attribute("casename")==""){
    QMessageBox::information(window(), "mainwindow",
                             tr("Please use 'GridSetup' reading in grid first"));
    return;
  }
  QString fileName = elem.attribute("directory")+"/"+elem.attribute("casename")+".vog";
 
  refdialog = new RefDialog(fileName);
  refdialog->show();
}


void MainWindow::resetSlider(){
}
void MainWindow::openSca()
{
     QString fileName =
             QFileDialog::getOpenFileName(this, tr("Open scalar File"),
                                          QDir::currentPath(),
                                          tr("scalar Files (*_sca.*)"));
     if (fileName.isEmpty()){
       QMessageBox::information(window(), "mainwindow",
                                tr("Please specify file name for reading"));
       return;
     }

     QFile file(fileName);
     if (!file.open(QFile::ReadOnly | QFile::Text)) {
       QMessageBox::warning(this, tr("Open scalar File"),
                            tr("Cannot read file %1:\n%2.")
                            .arg(fileName)
                            .arg(file.errorString()));
         return;
     }
     

       updateStatus(tr("Scalar File loaded"));
}


  
void MainWindow::openCase(){
  QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
  elem = elem.firstChildElement("gridSetup");
  if(!elem.attribute("casename").isEmpty()){
  
    if(doc.documentElement().attribute("updated")=="true"){
      
      int button = QMessageBox::question(this, tr("Open a New Case"),
                                         tr("Do you want to save the current case?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveXml();
      
      button = QMessageBox::question(this, tr("Open a New Case"),
                                     tr("Do you want to save .vars case?"),
                                     QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveVar();
  
      
    }
  }
  
  
   QString filename =
     QFileDialog::getOpenFileName(this, tr("Open a case"),
                                  QDir::currentPath(),
                                  tr("xml files (*.xml)"));
   if (filename.isEmpty()){
     emit updateStatus("open case cancelled");
     return;
   }

   QFile file(filename);
   if (!file.open(QIODevice::ReadOnly)){
     QMessageBox::information(window(), "mainwindow",
                              tr("Cannot open ") + filename + tr(" for reading!"));
     return;
   }
  QString errorStr;
  int errorLine;
  int errorColumn;
  if (!doc.setContent(&file, true, &errorStr, &errorLine,
                      &errorColumn)) {
    QMessageBox::information(window(), filename,
                             tr("Parse error at line %1, column %2:\n%3")
                             .arg(errorLine)
                             .arg(errorColumn)
                             .arg(errorStr));
    file.close();
    return;
  }
  file.close();

  updateStatus(filename+tr(" loaded"));
 
  if(viewer){
    delete viewer;
    viewer = new GLViewer(this);
  }
  if(central){
    delete central;
    central = new QStackedWidget(this);
    central->addWidget(viewer);
  }
  if(modBoundaries){
    delete modBoundaries;
    modBoundaries = 0;
  }
  if(boundaryView){
    delete boundaryView;
    boundaryView = 0;
  }
 
  
  if(bdWindow){
    delete bdWindow;
    bdWindow = 0;
  }
  setCentralWidget(central);
  previousWidget = 0;
  createFlowBar();
}
void MainWindow::newCase(){
  //if a case already exists, remind user to save the case
  QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
  elem = elem.firstChildElement("gridSetup");
  if(!elem.attribute("casename").isEmpty()){
  
    if(doc.documentElement().attribute("updated")=="true"){
      
      int button = QMessageBox::question(this, tr("Open a New Case"),
                                         tr("Do you want to save the current case?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveXml();
      
      button = QMessageBox::question(this, tr("Open a New Case"),
                                     tr("Do you want to save .vars case?"),
                                     QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveVar();
      
      
    }
  }
  

  //read in main.xml
  char* resourcepath = getenv("CHEMDEMOPATH");
    xmlpath = "./xml/";
   
  
    if(resourcepath){
      xmlpath = QString(resourcepath).append("xml/");
      pngpath = QString(resourcepath).append( "png/");
    }
    
    //  setAttribute(Qt::WA_DeleteOnClose); //will cause SF
    QString filename= xmlpath+"main.xml";
    QFile file(filename);
    if(!file.exists()){
      QMessageBox::information(window(), filename,
                               filename+ tr(" doesn't exist"));
      return;
    }
    if (!file.open(QIODevice::ReadOnly)){
      
      QMessageBox::information(window(), "main.xml",
                               tr("cannot open ")+ filename + tr("for reading"));
      
      return;
    }
  
    QString errorStr;
    int errorLine;
    int errorColumn;
    if (!doc.setContent(&file, true, &errorStr, &errorLine,
                        &errorColumn)) {
      QMessageBox::information(window(), filename,
                               tr("Parse error at line %1, column %2:\n%3")
                               .arg(errorLine)
                               .arg(errorColumn)
                               .arg(errorStr));
      file.close();
      return;
    }
  
  
    file.close();
   
   
    isNewCase= true;

    //clean up everything and recreate flowBar
 if(viewer){
    delete viewer;
    viewer = new GLViewer(this);
  }
  if(central){
    delete central;
    central = new QStackedWidget(this);
    central->addWidget(viewer);
  }
  if(modBoundaries){
    delete modBoundaries;
    modBoundaries = 0;
  }
  if(boundaryView){
    delete boundaryView;
    boundaryView = 0;
  }


  
  if(bdWindow){
    delete bdWindow;
    bdWindow = 0;
  }
  
  setCentralWidget(central);
  previousWidget = 0;
 
  createFlowBar();
   
    
    
}



bool MainWindow::saveVar()
{
  QDomElement root = doc.documentElement();
  root = root.firstChildElement("mainWindow");

 
  QString warningText;
  
  for( QDomElement elem = root.firstChildElement(); !elem.isNull(); elem=elem.nextSiblingElement()){
    if( elem.attribute("status")!="done"){
      warningText += elem.attribute("buttonTitle") + " ";
    }
  }
  
  if(!warningText.isEmpty()){


  
    warningText.replace("\n", "_");
    
    QMessageBox msgBox;
    msgBox.setText(warningText);
    msgBox.setInformativeText("Do you want to save ?");
    msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Save);
    
    QPushButton *gotoButton = msgBox.addButton(tr("go to the unfinished page"), QMessageBox::ActionRole);
    int ret = msgBox.exec();
    if(  msgBox.clickedButton() == gotoButton) {
      emit showStatus(true);
      for( QDomElement elem = root.firstChildElement(); !elem.isNull(); elem=elem.nextSiblingElement()){
        if( elem.attribute("status")!="done" && elem.hasAttribute("buttonIndex")){
          
          changePage(elem.attribute("buttonIndex").toInt());
          
          if(elem.attribute("element")=="boundaryWindow"){
            QDomElement theroot = doc.documentElement();
           
            QDomElement cndNode = theroot.firstChildElement("boundary_conditions");
            if(cndNode.isNull()){
              QMessageBox::information(window(), "mainwindow",
                                       "can not reach boundary_conditions node");
              return false;
            }
            int count = 0;
            QDomElement  elt = cndNode.firstChildElement();
            for(;!elt.isNull(); elt = elt.nextSiblingElement()){
              if(elt.firstChildElement().isNull() ||
                 elt.firstChildElement().attribute("status")!= "done")break;
              count++;
              
            }
            if(!elt.isNull()){
              // cndNode.setAttribute("currentIndex", count);
              selectCurrent(count);
             emit showStatus(true);
             }
            

          return false;
          // saveVar
        }
        }
      }

      
      
    }else{

     switch (ret) {
     case QMessageBox::Save:
       {
         QDomElement elem = root.firstChildElement("gridSetup");
         
         if(elem.isNull()|| !(elem.hasAttribute("casename"))){
           QMessageBox::warning(this, tr("saveVar"),
                              tr("no case to save")
                                );
           return false;
         }
         QString fileName = QDir::currentPath()+"/" +elem.attribute("casename")+tr(".vars");
         if(elem.hasAttribute("directory")) fileName = elem.attribute("directory")+"/" +elem.attribute("casename")+tr(".vars");
         
       
       
         fileName = QFileDialog::getSaveFileName(this, tr("Save .vars File"),
                                                 fileName,
                                                 tr("variable Files (*.vars)"));
         if(fileName.section('.', -1, -1)!="vars")fileName+=".vars";

         QFileInfo vogInfo(fileName);
         if(vogInfo.exists()){
           QString command = "mv " + fileName+ " " + fileName+".bak";
           int ret =  system(command.toStdString().c_str());
           if(!WIFEXITED(ret))
             {
             if(WIFSIGNALED(ret))
               {
                 QMessageBox::information(window(), "mainwindow",
                                          command + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
                 return false;
               }
             }
         
         }
         if(fileName.isNull()){
           emit updateStatus(" no file name specified for saving .vars"); 
           return false;
         }
         
         
  
         
         QFile file(fileName);
         if (!file.open(QFile::WriteOnly | QFile::Text)) {
           QMessageBox::warning(this, tr("save .vars file "),
                                tr("Cannot write file %1:\n%2.")
                                .arg(fileName)
                                .arg(file.errorString()));
           return false;
         }
       
       
         QTextStream out(&file);
         out<<"{"<<endl;
         
         elem = root.firstChildElement();
         vector<pair<int, QString> > nameMap;
         int count  = 100;
         for(; !elem.isNull(); elem=elem.nextSiblingElement(), count++){
           if(elem.hasAttribute("currentText")) {
             if(elem.hasAttribute("saveIndex"))
               nameMap.push_back(std::make_pair(elem.attribute("saveIndex").toInt(), elem.tagName()));
             else  nameMap.push_back(std::make_pair(count, elem.tagName()));
           }
         }
         sort(nameMap.begin(), nameMap.end());
         for(vector<pair<int, QString> > ::const_iterator p = nameMap.begin(); p!=nameMap.end(); p++){
           QDomElement elem = root.firstChildElement(p->second);
           out << elem.attribute("currentText") << endl;   
         }
         
         out<<"}"<<endl;
         updateStatus(fileName + tr(" saved"));
       
       }
         return true;
              
       // Save was clicked
       break;
       
     case QMessageBox::Cancel:
       // Cancel was clicked
       return false;
       break;
     default:
       return false;
       // should never be reached
       break;
     }

   }
  }else{
  
   QDomElement elem = root.firstChildElement("gridSetup");
  
  if(elem.isNull()|| !(elem.hasAttribute("casename"))){
    QMessageBox::warning(this, tr("saveVar"),
                         tr("no case to save")
                         );
    return false;
  }
  QString fileName = QDir::currentPath()+"/" +elem.attribute("casename")+tr(".vars");
  if(elem.hasAttribute("directory")) fileName = elem.attribute("directory")+"/" +elem.attribute("casename")+tr(".vars");

  
  
  fileName = QFileDialog::getSaveFileName(this, tr("Save .vars File"),
                                          fileName,
                                          tr("variable Files (*.vars)"));
  if(fileName.section('.', -1, -1)!="vars")fileName+=".vars";

  QFileInfo vogInfo(fileName);
  if(vogInfo.exists()){
    QString command = "mv " + fileName+ " " + fileName+".bak";
    int ret =  system(command.toStdString().c_str());
    if(!WIFEXITED(ret))
      {
        if(WIFSIGNALED(ret))
          {
            QMessageBox::information(window(), "mainwindow",
                                     command + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
            return false;
          }
      }

  }
  if(fileName.isNull()){
    emit updateStatus(" no file name specified for saving .vars"); 
    return false;
  }


  
 
  QFile file(fileName);
  if (!file.open(QFile::WriteOnly | QFile::Text)) {
    QMessageBox::warning(this, tr("save .vars file "),
                         tr("Cannot write file %1:\n%2.")
                         .arg(fileName)
                         .arg(file.errorString()));
    return false;
  }

  
  QTextStream out(&file);
  out<<"{"<<endl;
  
  elem = root.firstChildElement();
  vector<pair<int, QString> > nameMap;
  int count  = 100;
  for(; !elem.isNull(); elem=elem.nextSiblingElement(), count++){
    if(elem.hasAttribute("currentText")) {
      if(elem.hasAttribute("saveIndex"))
        nameMap.push_back(std::make_pair(elem.attribute("saveIndex").toInt(), elem.tagName()));
      else  nameMap.push_back(std::make_pair(count, elem.tagName()));
    }
  }
  sort(nameMap.begin(), nameMap.end());
  for(vector<pair<int, QString> > ::const_iterator p = nameMap.begin(); p!=nameMap.end(); p++){
    QDomElement elem = root.firstChildElement(p->second);
    out << elem.attribute("currentText") << endl;   
  }
    
  out<<"}"<<endl;
  updateStatus(fileName + tr(" saved"));

  file.close();
  return true;
  }
  return true;
}

 bool MainWindow::saveXml()
 {

   QDomElement root = doc.documentElement();
   root = root.firstChildElement("mainWindow");
   QDomElement elem = root.firstChildElement("gridSetup");
   if(elem.isNull()||!elem.hasAttribute("casename")){
     QMessageBox::warning(this, tr("saveXml"),
                          tr("no case to save")
                          );
     return false;
   }

   QString fileName = QDir::currentPath()+"/" +elem.attribute("casename")+tr(".xml");
    if(elem.hasAttribute("directory"))fileName = elem.attribute("directory")+"/" +elem.attribute("casename")+tr(".xml");
   
   
   fileName = QFileDialog::getSaveFileName(this, tr("Save .xml File"),
                                           fileName,
                                           tr("xml Files (*.xml)"));
   if(fileName.section('.', -1, -1)!="xml")fileName+=".xml";

   
   QFileInfo vogInfo(fileName);
   if(vogInfo.exists()){
     QString command = "mv " + fileName+ " " + fileName+".bak";
     int ret =  system(command.toStdString().c_str());
     if(!WIFEXITED(ret))
       {
         if(WIFSIGNALED(ret))
           {
             QMessageBox::information(window(), "mainwindow",
                                     command + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
             return false;
           }
       }
   }

   if(fileName.isNull()){
     emit updateStatus(" no file name specified for saving .xml"); 
     return false;
   }
   QFile file(fileName);
   if (!file.open(QFile::WriteOnly | QFile::Text)) {
     QMessageBox::warning(this, tr("Application"),
                          tr("Cannot write file %1:\n%2.")
                          .arg(fileName)
                          .arg(file.errorString()));
     return false;
   }
    doc.documentElement().setAttribute("updated","false");
   QTextStream out(&file);
   QApplication::setOverrideCursor(Qt::WaitCursor);
   doc.save(out, 2, QDomNode::EncodingFromDocument);
   QApplication::restoreOverrideCursor();
   updateStatus(fileName + tr(" saved"));
   file.close();
   return true;
 }



 bool MainWindow::saveImage()
 {
   return true;
   //   statusBar()->showMessage(tr("image printed"), 2000);
}




 void MainWindow::aboutPreprocess()
 {
   QMessageBox::about(this, tr("Preprocessing"),
                                   tr("the preprocessing set up boundary conditions, "
                                      "initial conditions, and other variables needed  "
                                      "for simulation"));
}


void MainWindow::aboutPostprocess()
{
  QMessageBox::about(this, tr("Postprocessing"),
                     tr("the postprocessing show boundaries, "
                        "cutting planes, or the whole grid  "
                        "after simulation"));
}



void MainWindow::adaptClicked()
{

  QDomElement root = doc.documentElement();
  root = root.firstChildElement("mainWindow");
  QDomElement elem = root.firstChildElement("gridSetup");
  if(elem.isNull()||!elem.hasAttribute("casename")){
    QMessageBox::warning(this, tr("FVMadapt"),
                          tr("Please first load in a grid from 'Grid Setup'")
                          );
     return;
   }

  if(adaptwindow){
    delete adaptwindow;
    adaptwindow = 0;
  }


    QString fileName = elem.attribute("directory")+"/"+elem.attribute("casename")+".vog";
 
    adaptwindow = new FVMAdapt(fileName);
    adaptwindow->show();
  
    // viewer->setAdaptWindow(adaptwindow);
  central->setCurrentWidget(viewer);
 
  connect(adaptwindow, SIGNAL(destroyed()), viewer, SLOT(cleanDoc()));
  connect(adaptwindow, SIGNAL(refineGrids()), this, SLOT(refineGrids()));
  connect(adaptwindow, SIGNAL(valueChanged(const QTreeWidgetItem*)), viewer, SLOT(updateDoc(const QTreeWidgetItem*)));
}


void MainWindow::vmClicked()
{
  
  if(mgviewer){
    central->removeWidget(mgviewer);
    delete mgviewer;
    mgviewer = 0;
  }
  if(vmwindow){
    delete vmwindow;
   vmwindow = 0;
  }
 
  
  mgviewer = new MGViewer(this);
  central->addWidget(mgviewer);
  central->setCurrentWidget(mgviewer);
  vmwindow = new VMergeWindow();
  vmwindow->show();
  
 connect(vmwindow, SIGNAL(loadGrid(QString)), mgviewer, SLOT(load_boundary(QString)));
 connect(vmwindow, SIGNAL(getGrid(QString)), mgviewer, SLOT(get_boundary(QString)));
 
 connect(mgviewer, SIGNAL(gridLoaded(const QStringList&)), vmwindow, SLOT(gridLoaded(const QStringList&)));
 connect(mgviewer, SIGNAL(pickCurrent(const IDOnly&)), vmwindow, SLOT(selectCurrent(const IDOnly&)));
 connect(vmwindow, SIGNAL(tcChanged(const IDMatrix&)), mgviewer, SLOT(transGrid(const IDMatrix&)));
 
 connect(vmwindow, SIGNAL(setCurrentColor(const IDColor&)), mgviewer, SLOT(setCurrentColor(const IDColor&)));
 connect(vmwindow, SIGNAL(setCurrentVisibility(const IDVisibility&)),
         mgviewer, SLOT(setCurrentVisibility(const IDVisibility&)));

 connect(vmwindow, SIGNAL(destroyed()), this, SLOT(vmClosed()));
 
}
void MainWindow::vmClosed(){
   if(mgviewer){
    central->removeWidget(mgviewer);
    delete mgviewer;
    mgviewer = 0;
   }
   central->setCurrentWidget(viewer); 
}


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
#include "mainwindow.h"
#include "glviewer.h"
#include "cutdialog.h"
#include "bdcndwindow.h"
#include "solverwindow.h"
#include "initcndwindow.h"
#include "importwindow.h"
#include "physicswindow.h"

#include <cstdlib>
#include <QString>
#include <string>
#include <utility>
#include <algorithm>

using std::string;
using std::vector;
using std::sort;
using std::pair;


///////////////////////////////////////////////////////////////////////////
//  public:
//    QWidget* createEditor(QWidget* parent,
//                          const QStyleOptionViewItem &option,
//                          const QModelIndex &index) const;
//
//  Simply inverts the visibility state of the boundary double-clicked on
//  in the boundary visibility model.
///////////////////////////////////////////////////////////////////////////


QWidget* showDelegate::createEditor(QWidget*, const QStyleOptionViewItem&,
				    const QModelIndex &index) const
{
  QString value = index.data(Qt::EditRole).toString();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  if (value == "hide") {
    model->setData(index, "show");
    
  } else {
    model->setData(index, "hide");
  }

  return NULL;
}






  
                                                              



bool colorDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
                                const QStyleOptionViewItem &,
                                    const QModelIndex &index){
  
                                                              
 
  if (event->type() == QEvent::MouseButtonPress) {
    QColor oldColor =   qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->background().color();
    QColor color = QColorDialog::getColor(oldColor);
    
    if(color.isValid())qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->setBackground( QBrush(color));
    else qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->setBackground( QBrush(oldColor));
    return false; //so that the selection can change
  }
  
  return true;
}



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
 
  //  tb->addAction(toolbar->toggleViewAction());   

  menuBar()->addSeparator();
  QAction *quitAct = new QAction(tr("&Quit"), this);
  menuBar()->addAction(quitAct);

 connect(quitAct, SIGNAL(triggered()),
         qApp, SLOT(quit()));  
 connect(openScaAct, SIGNAL(triggered()),
         this, SLOT(resetSlider()));

 viewMenu->addAction(tb->toggleViewAction()); 
}

void MainWindow::createDisplayBar(){

 
  QAction *showContoursAct = new QAction(tr("Show Contours"), this);
  QAction *showGridAct = new QAction(tr("Show  Coordinate Grid"), this);
  QAction *showShadingAct = new QAction(tr("Show Shading"), this);
  
  QAction *shadeType1 = new QAction(tr("Blue to Red"), this);
  QAction *shadeType2 = new QAction(tr("Blackbody"), this);
  QAction *shadeType3 = new QAction(tr("Pressure"), this);
  

  
  toolbar = addToolBar(tr("Cutplane Display"));
  insertToolBarBreak(toolbar);
 
 
  toolbar->addSeparator();
  toolbar->addAction(showContoursAct);
  toolbar->addSeparator();
  toolbar->addAction(showGridAct);
  toolbar->addSeparator();
  toolbar->addAction(showShadingAct);
    
  toolbar->addSeparator();
  toolbar->addAction(shadeType1);
  toolbar->addAction(shadeType2);
  toolbar->addAction(shadeType3);
  toolbar->addSeparator();
  
  //  slider = new QSlider(Qt::Horizontal, toolbar);
  //slider->setRange(5, 50);
  //slider->setValue(10);
  // toolbar->addWidget(slider);
  
 
 
  connect(showGridAct, SIGNAL(triggered()),
          viewer, SLOT(toggleGrid()));
  connect(showContoursAct, SIGNAL(triggered()),
          viewer, SLOT(toggleContours()));
  connect(showShadingAct, SIGNAL(triggered()),
          viewer, SLOT(toggleShading()));
  connect(shadeType1, SIGNAL(triggered()),
          viewer, SLOT(setShadeType1()));
  connect(shadeType2, SIGNAL(triggered()),
          viewer, SLOT(setShadeType2()));
  connect(shadeType3, SIGNAL(triggered()),
          viewer, SLOT(setShadeType3()));
 
  
  // connect(slider, SIGNAL(valueChanged(int)),
  //      viewer, SLOT(changeContours(int)));
 
  viewMenu->addAction(toolbar->toggleViewAction());
  tb->addSeparator();
  tb->addAction(toolbar->toggleViewAction());
  tb->setStyleSheet("* { color: rgb(120,60, 0) }");
  toolbar->setStyleSheet("* { color: rgb(120,60, 0) }");
}


void MainWindow::createVisBar(){
  
  visbar = new QToolBar("Visualization");

  addToolBar(Qt::TopToolBarArea,visbar );
  insertToolBarBreak(visbar);

  QAction *selectBoundaryAct = new QAction(tr("Select Boundary"), this);
  visbar->addAction(selectBoundaryAct);
  connect(selectBoundaryAct, SIGNAL(triggered()),
          this, SLOT(selectBoundary())); 
  visbar->addSeparator();  
  QAction *showBoundariesAct = new QAction(tr("show Boundaries"), this);
  showBoundariesAct->setCheckable(true);
  visbar->addAction(showBoundariesAct);
  connect(showBoundariesAct, SIGNAL(triggered(bool)),
          viewer, SLOT(showBoundaries(bool))); 
  visbar->addSeparator();
  
  QAction* resetAct = new QAction(tr("Reset"), this);
  visbar->addAction(resetAct);
  connect(resetAct, SIGNAL(triggered()),
          viewer, SLOT(reset()));
  visbar->addSeparator();
  QAction* fitAct = new QAction(tr("Fit"), this);
  visbar->addAction(fitAct);
  connect(fitAct, SIGNAL(triggered()),
          viewer, SLOT(fit()));

  cutAct = new QAction(tr("Cut"), this);  
  visbar->addSeparator();
  visbar->addAction(cutAct);
  connect(cutAct, SIGNAL(triggered()),
          this, SLOT(cut()));
  //  connect(cutAct, SIGNAL(clicked()),
  //      this, SLOT(resetSlider())); 
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
  visbar->setStyleSheet("* { color: rgb(120, 60, 0) }");  
 }
  
void MainWindow::toggleViewer(){
  if(central->currentWidget()==viewer && previousWidget!=0){
    central->setCurrentWidget(previousWidget);
  }else{
    previousWidget = central->currentWidget();
    central->setCurrentWidget(viewer);
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
    }else  if(elem.attribute("element")=="physicsWindow"){
      newWindow=new PhysicsWindow(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(newWindow, SIGNAL(stateChanged()), this, SIGNAL(stateChanged()));
      connect(newWindow, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
    }else  if(elem.attribute("element")=="initialWindow"){
      newWindow=new InitCndWindow(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
    }else  if(elem.attribute("element")=="panel"){
      newWindow=new OptionPage(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
    } else  if(elem.attribute("element")=="page"){
      newWindow=new Page(elem, theroot);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
    }   
    if(newWindow && elem.attribute("inDock")!="true"){
      int i = central->addWidget(newWindow);
      elem.setAttribute("widgetIndex", i);
    }
  }
  connect(flowbarButtons, SIGNAL(buttonClicked(int)), this, SLOT(changePage(int)));
  flowbar->setStyleSheet("QPushButton { color: darkGreen }");

}
void MainWindow::changePage(int index){

  QDomElement theroot = doc.documentElement();
  QDomElement elem = theroot.firstChildElement("mainWindow");
  elem=elem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), "main xml file",
                         tr(" mainWindow has no child" ));
    return;
  }
  for(int i=0; i< index; i++)elem=elem.nextSiblingElement();
  if(elem.hasAttribute("widgetIndex")) 
    central->setCurrentIndex(elem.attribute("widgetIndex").toInt());
  else  central->setCurrentIndex(0);
  if(elem.tagName()=="saveVar")saveVar();
  if(elem.tagName()=="gridSetup")setGrid(elem);
  if(elem.attribute("element")=="boundaryWindow")setBoundary(elem);
  if(elem.attribute("inDock")=="true") dock->setWidget(bdWindow);
  else  dock->setWidget(statusWindow);
}


  

void MainWindow::createDockWindow(){
  dock = new QDockWidget("dock window", this);
  dock->setAllowedAreas(Qt::LeftDockWidgetArea );
  statusEdit = new QTextEdit(this);
  // progressEdit = new QTextEdit;
  
  statusWindow = new QGroupBox(this);
  statusWindow->setFlat(true);
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  QPushButton* clearAllButton = new QPushButton(tr("Clear All"));
  QPushButton* clearLastButton = new QPushButton(tr("Clear Last"));
  buttonLayout->addWidget(clearAllButton);
  buttonLayout->addWidget(clearLastButton);
  QVBoxLayout* mainLayout=new QVBoxLayout;
  // mainLayout->addWidget(progressEdit);
  
  mainLayout->addWidget(statusEdit);
  
  mainLayout->addLayout(buttonLayout);
  statusWindow->setLayout(mainLayout);
  
   
  dock->setWidget(statusWindow);
  addDockWidget(Qt::LeftDockWidgetArea, dock);
  viewMenu->addAction(dock->toggleViewAction());
  connect(clearAllButton, SIGNAL(clicked()), this, SLOT(clearAllStatus()));
  connect(clearLastButton, SIGNAL(clicked()), this, SLOT(clearLastStatus()));
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
  bdWindow = 0;
  cutdialog = 0;
  slider = 0;
  dock = 0;
  
  cutAct = 0;
  toolbar=0;
  flowbar =0;
  flowbarButtons = 0;
  bdWindow = 0;
  visbar =0;
  bdButtonDown = false;
  setCentralWidget(central);
  previousWidget = 0;
  statusWindow = 0;
  createMenu();
  createDisplayBar();
  
  createVisBar();
  createDockWindow();
  createFlowBar();
  hideVisBar();
  hideDisplayBar();
  //  statusBar()->showMessage(tr("please load a file with file menu and select a model"));
   
   setWindowTitle(tr("chem demo"));
   resize(800, 600);
   
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

  QPixmap originalPixmap = QPixmap(); // clear image for low memory situations
  // on embedded devices.
  originalPixmap = QPixmap::grabWidget(viewer);
  QString format = "png";
  QString initialPath = QDir::currentPath() + tr("/untitled.") + format;

  QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"),
                                                  initialPath,
                                                  tr("%1 Files (*.%2);;All Files (*)")
                                                  .arg(format.toUpper())
                                                  .arg(format));

  if (!fileName.isEmpty()){
    bool filesaved =  originalPixmap.save(fileName, format.toAscii());
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

  

  QStringList boundary_names;
  QString format = "volume grid file(*.vog)";
  if(theelem.hasAttribute("format")) format = theelem.attribute("format");
  QString initialPath =QDir::currentPath();
  if(theelem.hasAttribute("initialPath"))initialPath = theelem.attribute("initialPath");
  QString fileName  = initialPath;
  if(theelem.hasAttribute("casename")) fileName =  theelem.attribute("directory")+"/"+theelem.attribute("casename")+".vog";

 
  fileName = QFileDialog::getOpenFileName(this, tr("Load Grid"),
                                          fileName,
                                          format);
  
  QString caseName;
  if(fileName==""){
    //no error message in  case of 'cancel' is pressed 
    
    return;
  }
  
  QStringList formatList = fileName.split('.');
  if(formatList.size()!=2){
    QMessageBox::warning(this, tr("Application"),
                         fileName + tr(" has no postfix"));
    return;
  }
  
  bool loaded = false;
  if(formatList[1]=="vog"){
    
    QString surfFileName = fileName.section('.', 0, 0)+".surface";
     caseName =  (fileName.section('.', 0, 0)).section('/', -1);
  
     QFileInfo surfInfo(surfFileName);
     QFileInfo vogInfo(fileName);
     if(!(surfInfo.exists()) || surfInfo.created() < vogInfo.created()){
    
      QString script_filename = "./output/vog2surface_"+caseName;
      QString out_filename="./output/vog2surface_"+caseName+".out";
      
      QFile outfile(script_filename);
      if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::information(window(), "mainwindow",
                                 tr("Cannot open ") + script_filename + tr(" for writing!"));
        return;
      }
      
      QTextStream out(&outfile);
    
    
      QString command2 = "vog2surface " + fileName.section('.', 0, 0);
      emit updateStatus(command2);
      out <<"#!/bin/bash"<<endl;
      out <<"exec 6>&1"<<endl;
      out <<"exec 7>&2"<<endl;
      out<< "exec &> "<< out_filename <<endl;
      out << command2 << endl;
      out<<"exec 1>&6 6>&- " << endl;
      out<<"exec 2>&7 7>&- " << endl;

      outfile.close();
      QString command3 = "chmod 777 " + script_filename;

      int ret =  system(command3.toStdString().c_str());



      if(!WIFEXITED(ret))
        {
          if(WIFSIGNALED(ret))
            {
              QMessageBox::information(window(), "mainwindow",
                                       command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
               theelem.removeAttribute("casename");
              return;
            }
        }
        
      ret = system(script_filename.toStdString().c_str());
      if(!WIFEXITED(ret)){
        if(WIFSIGNALED(ret)){
          QMessageBox::information(window(), "mainwindow",
                                   script_filename + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
           theelem.removeAttribute("casename");
          return;
        }
      }



      
      QFile file(out_filename);
      if (!file.open(QFile::ReadOnly | QFile::Text)) {

        QMessageBox::information(window(), "mainwindow",
                                 tr("Cannot open ") + out_filename + tr(" for reading!"));

        theelem.removeAttribute("casename");
        
        return;
      }
    
      QTextStream in(&file);
      QApplication::setOverrideCursor(Qt::WaitCursor);
      //      msgBox.setDetailedText(in.readAll());
      emit updateStatus(in.readAll());
      QApplication::restoreOverrideCursor();
      file.close();
    }
    


  
   
    // must setCurrentWidget(viewer)first, then load_boundary
    // must leave these two lines at the end of function
    //  central->setCurrentWidget(viewer);   
    loaded =   viewer->load_boundary(surfFileName, boundary_names); // and setup the GLWidget.
    viewer->show();
  }else{
    
    // central->setCurrentWidget(viewer);   
    loaded =  viewer->load_boundary(fileName, boundary_names); // and setup the GLWidget.
    viewer->show(); 
  }
  viewer->reset();
  if(loaded){
    if(theelem.hasAttribute("casename") && theelem.attribute("casename")!=caseName) {
      int button = QMessageBox::question(this, tr("a new grid is loaded"),
                                         tr("Do you want to save the old case ?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveXml();
      
      button = QMessageBox::question(this, tr("a new grid is loaded"),
                                         tr("Do you want to save the .vars file of the old case ?"),
                                     QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveVar();

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


    
    int first= fileName.lastIndexOf('/');
    int last = fileName.lastIndexOf('.');
    QString casename = fileName.mid(first+1, last-first-1);
    QString directory = fileName.left(first);
    theelem.setAttribute("casename", casename);
    theelem.setAttribute("directory", directory);
    theelem.setAttribute("boundary_names", boundary_names.join(","));
      
    updateStatus(fileName + tr(" loaded"));
    theelem.setAttribute("status", "done");
   
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
  

  
   
    
    
  }else{
    // theelem.removeAttribute("casename");
    updateStatus(fileName + tr(" not loaded"));
  }
  
  
}

void MainWindow::openVog(){
  
}

void MainWindow::setBoundary(QDomElement& elem){
  bdButtonDown = true;
  dock->show();
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
   
    bdWindow = new BdCndWindow(elem, theroot, bdnames, boundaryView);
    connect(this, SIGNAL(setCurrent(QModelIndex)), bdWindow, SLOT(setCurrent(QModelIndex)));
    connect(bdWindow, SIGNAL(closed()), this, SLOT(bdWindowClosed())); 
    connect(bdWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
    connect(bdWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
    connect(this, SIGNAL(stateChanged()), bdWindow, SLOT(changeState()));
    connect(this, SIGNAL(componentsChanged()), bdWindow, SIGNAL(componentsChanged()));
    
    updateStatus(tr("start bounadry conditions setup ..."));
  }

  dock->setWidget(bdWindow);

 }

 void MainWindow::bdWindowClosed(){
  bdButtonDown = false;
  viewer->clearCurrent();
  dock->setWidget(statusWindow);
  // hideVisBar();
  //  hideDisplayBar();
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

void MainWindow::showDisplayBar(){
 if(toolbar) toolbar->show();
}


void MainWindow::hideDisplayBar(){
 if(toolbar) toolbar->hide();
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
  

   
  showDisplayBar();
  if(cutdialog) {
    delete cutdialog;
    cutdialog = 0;
  }
  
  if(viewer ==0) return;
  cutdialog = new CutDialog(viewer->boundaryBoxSize());
  cutdialog->show();
  connect(cutdialog, SIGNAL(cutInfoChanged(cutplane_info&)), viewer, SLOT(previewCut(cutplane_info&)));
  connect(cutdialog, SIGNAL(cutPressed()), viewer, SLOT(cut()));
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

  QDomElement elem = root.firstChildElement();
  QString warningText;
  for(; !elem.isNull(); elem=elem.nextSiblingElement()){
    if( elem.attribute("status")!="done"){
      warningText += elem.attribute("buttonTitle") + " ";
    }
  }
  if(!warningText.isEmpty()){
    warningText.replace("\n", "_");
    QMessageBox::warning(this, tr("saveVar"),
                         warningText+ tr(" not done yet")
                         );
    
  }
  
  elem = root.firstChildElement("gridSetup");
  
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




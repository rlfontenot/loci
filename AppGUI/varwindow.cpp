#include <QtGui>
#include <QApplication>
#include <QRegExp>
#include "varwindow.h"
#include "glviewer.h"
#include "bdcndwindow.h"
#include "solverwindow.h"
#include "initcndwindow.h"
#include "physicswindow.h"


/*!
  \class VarWindow
  
  \brief VarWindow is the main window for var file generation.
  
  
  
  VarWindow is a main window, it includes a file menu,
  a tool bar and a visualization bar on the top,  a flow bar in  the left,
  and a dock window in the right. The central area is the user input interface,
  the dock window is for graphics viewer.

  On the bottom is a status bar, it contains a temporary status tip and permanent status label.
  When the cursor is on one of the buttons in the flow bar, the status tip will show if this step
  is done or not. The permanent label shows the current status.

    
*/




/*!
    Create the menu bar. It includes file menu,  
*/
void VarWindow::createMenu(){
  char* resourcepath = getenv("CHEMDEMOPATH");
  QString xmlpath = "./xml/";
  QString pngpath = "./png/";
  
  if(resourcepath){
    xmlpath = QString(resourcepath).append("xml/");
    pngpath = QString(resourcepath).append( "png/");
  }
  
  tb = new QToolBar(tr("tool bar"),this);
  tb->setWindowTitle(tr("File Actions"));
  addToolBar(tb);
 
  
  QAction *newCaseAct = new QAction(QIcon( pngpath+ "filenew.png"), tr("&New case"), this);
  newCaseAct->setShortcut(tr("Ctrl+N"));
  connect(newCaseAct, SIGNAL(triggered()), this, SLOT(newCase()));
 
  QAction *openCaseAct = new QAction(QIcon( pngpath+ "fileopen.png"), tr("&Open case"), this);
  openCaseAct->setShortcut(tr("Ctrl+O"));
  connect(openCaseAct, SIGNAL(triggered()), this, SLOT(openCase()));

  QAction *saveXmlAct = new QAction( QIcon( pngpath+"save.png"), tr("&Save Case"), this);
  saveXmlAct->setShortcut(tr("Ctrl+S"));
  connect(saveXmlAct, SIGNAL(triggered()), this, SLOT(saveXml()));
  
  QAction *saveImageAct = new QAction(  QIcon(pngpath+"snapshot.png"), tr("snapshot"), this);
  connect(saveImageAct, SIGNAL(triggered()), this, SLOT(snapshot()));
  
  
  QAction *exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  connect(exitAct, SIGNAL(triggered()), qApp, SLOT(quit()));
  
   
  QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newCaseAct);
  fileMenu->addAction(openCaseAct);
  fileMenu->addAction(saveXmlAct);
  fileMenu->addAction(exitAct);
  

 
  
 
  menuBar()->addSeparator();
  viewMenu = menuBar()->addMenu(tr("&View"));
 
  tb->addAction(newCaseAct);
  tb->addAction(openCaseAct);
  tb->addAction(saveXmlAct);
  tb->addAction(saveImageAct);

  QAction* showStatusAct = new QAction(tr("Show Status"),this);
  connect(showStatusAct, SIGNAL(triggered()), this, SLOT(toggleShowStatus()));
  tb->addSeparator();
  tb->addAction(showStatusAct);

 //  QAction* updateModelAct = new QAction(tr("Update Models"),this);
//   connect(showStatusAct, SIGNAL(triggered()), this, SLOT(updateModels()));
//   tb->addSeparator();
//   tb->addAction(updateModelAct);
  
  
 
  menuBar()->addSeparator();
  QAction *quitAct = new QAction(tr("&Quit"), this);
  menuBar()->addAction(quitAct);

  connect(quitAct, SIGNAL(triggered()),
          qApp, SLOT(quit()));  
  

  viewMenu->addAction(tb->toggleViewAction());
    
    
  menuBar()->addSeparator();
  QMenu* helpMenu = menuBar()->addMenu(tr("&Help"));
  {
    QSignalMapper *helpMapper = new QSignalMapper(this); 
    QDomElement theroot = doc.documentElement();
    QDomElement elem = theroot.firstChildElement("help");
    for(QDomElement elt = elem.firstChildElement();
        !elt.isNull(); elt= elt.nextSiblingElement()){
      QAction* action= new QAction(elt.tagName(), this);
      connect(action, SIGNAL(triggered()), helpMapper, SLOT(map()));
      helpMapper->setMapping(action, elt.tagName());
      helpMenu->addAction(action);
    }
    connect(helpMapper, SIGNAL(mapped(const QString&)), this, SLOT(help(const QString&)));
  }
   

}

QStringList getComp(QString filename){
 
  QStringList comp;
  QFile file(filename);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    qDebug()<<"Cannot open "<< filename << " for reading!";
    
    return comp;
  }
  
  QTextStream in(&file);
  QString contents = in.readAll();
  contents = contents.trimmed();
  contents = contents.section("species", 1, -1);
  contents = contents.trimmed();
  contents = contents.section("=", 1, -1);
  contents = contents.trimmed();
  contents = contents.section("{", 1, -1);
  contents = contents.trimmed();
  contents = contents.section("}", 0, 0);
  contents = contents.trimmed();
 
  if(contents.size()==0)
    {
      qDebug()<< "can not get 'species' in " << filename;
      return comp;
    }
  
  QStringList tmpList = contents.split(";", QString::SkipEmptyParts);
  for(int i = 0; i < tmpList.count(); i++){
    QString name;
    QRegExp rx("(=|:)");
    int firstIndex =rx.indexIn(tmpList[i]);
    if(firstIndex >=0) name = tmpList[i].left(firstIndex);
    else name = tmpList[i];
    name = name.trimmed();
    comp << name;
  }
 
  return comp;
}
void VarWindow::updateModels(){
  char* resourcepath = getenv("CHEMISTRY_DATABASE");
  QDir directory = QDir::currentPath();
  QStringList nameFilters;
  nameFilters<<QString("*.mdl");
  QStringList files = directory.entryList(nameFilters, QDir::Files|QDir::NoSymLinks);
  if(files.count()==0){
    directory = QString(resourcepath).append("/data_base/models/");
    files= directory.entryList(nameFilters, QDir::Files|QDir::NoSymLinks);
  }
  if(files.count()==0){
    QMessageBox::warning(window(), "update models",
                         tr("can not find any model files")
                         );
    return;
  }
  
  QDomElement modelElem = doc.documentElement().firstChildElement("models");
  QDomElement singleElem = modelElem.firstChildElement("single_component_model");
  QDomElement multiElem = modelElem.firstChildElement("multiple_components_model");
 
 
  
  if(!singleElem.isNull()){
    for(QDomElement celem =singleElem.firstChildElement();
        !celem.isNull(); celem = celem.nextSiblingElement())singleElem.removeChild(celem);
  }else{
    QDomElement newSingleNode = doc.createElement("single_component_model");
    modelElem.appendChild(newSingleNode);
    singleElem = newSingleNode;
  }

  if(!multiElem.isNull()){
     for(QDomElement celem =multiElem.firstChildElement();
        !celem.isNull(); celem = celem.nextSiblingElement())multiElem.removeChild(celem);
  }else{
    QDomElement newMultiNode = doc.createElement("multiple_components_model");
    modelElem.appendChild(newMultiNode);
    multiElem = newMultiNode;
  }
  
  
  for(int i = 0; i < files.count(); i++){
    QStringList components = getComp(directory.absolutePath()+"/"+files[i]);
    if(components.count()==0)continue;
     
    int first= files[i].lastIndexOf('/');
    int last = files[i].lastIndexOf('.');
    QString casename = files[i].mid(first+1, last-first-1);
    QDomElement mdlnode = doc.createElement(casename);
    
    
    if(components.count()==1){
      singleElem.appendChild(mdlnode);
    }else{
      
      mdlnode.setAttribute("components", components.join(tr(","))); 
      multiElem.appendChild(mdlnode);
    }
  }
}




  
void VarWindow::createVisBar(){
  
  visbar = new QToolBar("Visualization");
  addToolBar(Qt::TopToolBarArea,visbar );
  insertToolBarBreak(visbar);
  QAction *clearBoundaryAct = new QAction(tr("Clear"), this);
  visbar->addAction(clearBoundaryAct);
  connect(clearBoundaryAct, SIGNAL(triggered()),
          viewer, SLOT(clearCurrent())); 
  visbar->addSeparator();
  
  QAction *showBoundariesAct = new QAction(tr("show Boundaries"), this);
  visbar->addAction(showBoundariesAct);
  connect(showBoundariesAct, SIGNAL(triggered()),
          viewer, SLOT(showBoundaries())); 
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


 
  viewMenu->addAction(visbar->toggleViewAction());
  tb->addSeparator();
  tb->addAction(visbar->toggleViewAction());


  statusLabel = new QLabel(this);
  statusBar()->addPermanentWidget(statusLabel);


}







void VarWindow::toggleShowStatus(){
  displayStatus = !displayStatus;
  emit showStatus(displayStatus);
}



void VarWindow::createFlowBar(){

  if(central){
    delete central;
  }
  central = new QStackedWidget;
  setCentralWidget(central);
    
  //create flowbar
  if(flowbar){
    delete flowbar;
    flowbar = 0;
    
  }
  flowbar = new QToolBar;
  addToolBar(Qt::LeftToolBarArea,flowbar );
  
  //go to element 'mainWindow'
  QDomElement theroot = doc.documentElement();
  QDomElement theelem = theroot.firstChildElement("mainWindow");

  QDomElement elem = theelem.firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         elem.tagName()+ tr(" has no child")
                         );
    return;
  }
  
  //create flowbarButtons
  if(flowbarButtons){
    delete flowbarButtons;
    flowbarButtons = 0;
  }
  flowbarButtons = new QButtonGroup(this);
  flowbarButtons->setExclusive(true);


  QGroupBox *flowGroup = new QGroupBox(tr("Flowbar"));
  QVBoxLayout* flowLayout = new QVBoxLayout;
  
  int count=0;
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) { 
 
    QPushButton* newButton = new QPushButton(elem.attribute("buttonTitle"), this);
    newButton->setCheckable(true);
    if(count==0) newButton->setChecked(true);
    flowLayout->addWidget(newButton);
    flowbarButtons->addButton(newButton, count);
    elem.setAttribute("buttonIndex", count);
    newButton->setToolTip(elem.attribute("toolTip"));
    newButton->setWhatsThis(elem.attribute("whatsThis"));
    newButton->setStatusTip(elem.attribute("status"));
 
    QWidget* newWindow=0;
    
    if(elem.attribute("element")=="solverWindow"){
      newWindow=new SolverWindow(elem);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
      
    }else  if(elem.attribute("element")=="physicsWindow"){
      newWindow=new PhysicsWindow(elem);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(newWindow, SIGNAL(stateChanged()), this, SIGNAL(stateChanged()));
      connect(newWindow, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    }else  if(elem.attribute("element")=="initialWindow"){
       
      newWindow=new InitCndWindow(elem);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
      connect(newWindow, SIGNAL( valueChanged(const QTreeWidgetItem*)),
              viewer, SLOT(updateDoc(const QTreeWidgetItem*)));
    }else  if(elem.attribute("element")=="boundaryWindow"){
      newWindow = new QWidget();//set it up later
      bdWindow = newWindow;
    }else if(!elem.hasAttribute("element")){
      newWindow = new QWidget();//dummy widget  
      
    }else  if(elem.attribute("element")=="panel"){
      newWindow=new VarPanel(elem);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    } else  if(elem.attribute("element")=="page"){
      newWindow=new VarPage(elem);
      connect(newWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
      connect(newWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
      connect(this, SIGNAL(stateChanged()), newWindow, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), newWindow, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), newWindow, SLOT(updateShowStatus(const bool&)));
    } 
      
    if(newWindow){
      central->addWidget(newWindow);
    }
  }
  
  connect(flowbarButtons, SIGNAL(buttonClicked(int)), this, SLOT(changePage(int)));
  flowGroup->setLayout(flowLayout);
  flowbar->addWidget(flowGroup);
}


void VarWindow::changePage(int index){
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
    setGrid();
  }
 
  if(elem.attribute("inDock")=="true"){
    viewerDock->setFloating(true);
    viewerDock->show();
    viewerDock->raise();
  }else{
    viewerDock->setFloating(false);
    viewerDock->hide();
  }
  
  central->setCurrentIndex(index);
}
  


  



void VarWindow::updateStatus(const QString& s){
  statusLabel->setText(s);
}

void VarWindow::updateStatusTip(int button){
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


VarWindow::VarWindow()
{
  QWidget::setAttribute(Qt::WA_DeleteOnClose, true);
  
  //first use main.xml set up doc
  char* resourcepath = getenv("CHEMDEMOPATH");
  QString xmlpath = "./xml/";
  QString pngpath = "./png/";
  
  if(resourcepath){
    xmlpath = QString(resourcepath).append("xml/");
    pngpath = QString(resourcepath).append( "png/");
  }
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
  updateModels();
  
  //initialize member data 
  bdWindow = 0;
  flowbar =0;
  flowbarButtons = 0;
  displayStatus = false;
  visbar =0;
  central = 0;

  //create viewer and viewerDock
  viewer = new GLViewer(this);
  
  viewerDock  = new QDockWidget("Graphics Viewer", this);
  viewerDock->setAllowedAreas(Qt::RightDockWidgetArea );
  viewerDock->setWidget(viewer);
  addDockWidget(Qt::RightDockWidgetArea, viewerDock);
  
  

  //create menus and tool bars
  
  createMenu();
  createVisBar();
  createFlowBar();
  hideVisBar();
  
 
  tb->addSeparator();
  tb->addAction(viewerDock->toggleViewAction());   

 
     
  setWindowTitle(tr("Generate a .var file"));
  updateStatus(tr("Please use 'Grid Setup' to load  grid information, or use file menu to open a case"));
  
  statusBar()->showMessage(tr("Ready"));
 
 
   setGrid();
}


////////////////////////////////////////
//  public:
//    QSize sizeHint() const;
//
//  Requests a default size of 800x600.
////////////////////////////////////////

QSize VarWindow::sizeHint() const
{
  return QSize(1024, 768);
}



void VarWindow::snapshot(){
  
  QImage pic ;
  pic = viewer->grabFrameBuffer(true);
 

 

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


  
  
void VarWindow::setGrid()
{

  QDomElement theelem = doc.documentElement().firstChildElement("mainWindow").firstChildElement("gridSetup");
  if(theelem.isNull()){
    QMessageBox::information(window(), tr("main.xml"),
                             tr("can not find element 'gridSetup'"));  
    return;
  }
  
  //set up the default filename

  QString format = "volume grid file(*.vog)";
  QString initialPath =QDir::currentPath();
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

  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
 
  //viewer load in grid 
  viewerDock->show();
  viewerDock->raise();
  viewerDock->setFloating(true);
 
  QStringList boundary_names;
  bool loaded =   viewer->load_boundary(fileName, boundary_names); // and setup the GLWidget.
  
  
  if(loaded){
    
  
    
    //if different case, remind the user to save the case
    if(theelem.hasAttribute("casename") ) {
      int button = QMessageBox::question(this, tr("a new grid is loaded"),
                                         tr("Do you want to save the old case ?"),
                                         QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveXml();
      
      button = QMessageBox::question(this, tr("a new grid is loaded"),
                                     tr("Do you want to save the .vars file of the old case ?"),
                                     QMessageBox::Ok|QMessageBox::No, QMessageBox::Ok); 
      if(button == QMessageBox::Ok) saveVar();
    }

    //set up the new case
    theelem.setAttribute("casename", casename);
    theelem.setAttribute("directory", directory);
    theelem.setAttribute("boundary_names", boundary_names.join(","));
    
    updateStatus(fileName + tr(" loaded"));
    theelem.setAttribute("status", "done");

    //set up the cndNode
    QDomElement oldCndNode = doc.documentElement().firstChildElement("boundary_conditions");
    if(oldCndNode.isNull()){
      QMessageBox::warning(window(), ".xml",
                           tr("can not find element 'boundary_conditions'")
                           );
      return;
    }

    //create new cndNode
    QDomElement cndNode =  doc.createElement("boundary_conditions");
    
    for(int i =0; i < boundary_names.size(); i++){
      QDomElement  aNode = doc.createElement(boundary_names[i]);
      cndNode.appendChild(aNode);
    }
    cndNode.setAttribute("currentIndex", "0");
    //replace old cndNode
    doc.documentElement().replaceChild(cndNode, oldCndNode);
    
    //create new boundary window
    setBoundary();
    
  }else{
    // theelem.removeAttribute("casename");
    updateStatus(fileName + tr(" not loaded"));
  }
  
  viewer->clearCurrent();//do not select any boundary 
}


/*!
  Create a new boundary window, connect its signal and slots with this,
  and replace the old boundary window in central.
*/
void VarWindow::setBoundary(){
  int index = central->indexOf(bdWindow);
  
  if(index < 0){
    QMessageBox::warning(window(), "var file generation",
                         tr("boundary condition window doesn't exist in central" ));
    return;
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
  if(elem.isNull()||elem.attribute("element")!="boundaryWindow"){
    QMessageBox::warning(window(), "main xml file",
                         tr(" boundary condition element doesn't exist " ));
    return;
  }

  bdWindow = new BdCndWindow(elem, viewer);
 
  connect(bdWindow, SIGNAL(updateStatus(const QString&)), this, SLOT(updateStatus(const QString&)));
  connect(bdWindow, SIGNAL(updateStatusTip(int)), this, SLOT(updateStatusTip(int)));
  connect(this, SIGNAL(stateChanged()), bdWindow, SLOT(changeState()));
  connect(this, SIGNAL(componentsChanged()), bdWindow, SIGNAL(componentsChanged()));
  connect(this, SIGNAL(showStatus(const bool &)), bdWindow, SIGNAL(showStatus(const bool &)));
  

  QWidget* dummyWidget= central->widget(index);
  central->removeWidget(dummyWidget);
  central->insertWidget(index,bdWindow);
  delete dummyWidget;
}












    


void VarWindow::showVisBar(){
  if(visbar) visbar->show();
}


void VarWindow::hideVisBar(){
  if(visbar) visbar->hide();
}







  




/*!
  Open a case, which is an xml file. The flow bar will be recreated.
  if casename in new xml file is empty, then set up grid.
  otherwise, just let viewer load in grid.
*/
  
void VarWindow::openCase(){
  //remind the user to save the old case
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
  

  //read in new xml file
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
 
  //recreate flow bar 
  createFlowBar();

  elem = doc.documentElement().firstChildElement("mainWindow");
  elem = elem.firstChildElement("gridSetup");
  if(elem.attribute("casename").isEmpty()){
    setGrid();
  }else{
    //viewer load in grid 
    viewerDock->show();
    viewerDock->raise();
    viewerDock->setFloating(true);
    QString fileName =  elem.attribute("directory")+"/"+elem.attribute("casename")+".vog";
    QStringList boundary_names;
    bool loaded = viewer->load_boundary(fileName, boundary_names); // and setup the GLWidget. 
    if(!loaded){
      QMessageBox::information(window(), "Open Case",
                               tr("Can not load in ")+fileName);
    }
    //create new boundary window
    setBoundary();
    //go to current page
     elem = doc.documentElement().firstChildElement("mainWindow");
     elem = elem.firstChildElement("gridSetup");
     for(; !elem.isNull(); elem = elem.nextSiblingElement()){
       if(elem.attribute("status")!="done")break;
     }
     if((!elem.isNull())){
       if(elem.hasAttribute("buttonIndex")){
         flowbarButtons->button(elem.attribute("buttonIndex").toInt())->click();
       }
     }
                               
  }
}


/*!
  Start a new case using file 'main.xml'. The flow bar will be recreated. This command is
  necessary only if the original case is start with a file other that 'main.xml'.
*/

void VarWindow::newCase(){
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
  QString xmlpath = "./xml/";
  QString pngpath = "./png/";
  
  if(resourcepath){
    xmlpath = QString(resourcepath).append("xml/");
    pngpath = QString(resourcepath).append( "png/");
  }
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

  //recreate flow bar
  createFlowBar();
   
  //set grid    
  setGrid();
}



bool VarWindow::saveVar()
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
          
          flowbarButtons->button(elem.attribute("buttonIndex").toInt())->click();
                    
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
              qobject_cast<BdCndWindow*>(bdWindow)->selectCurrent(count);
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

bool VarWindow::saveXml()
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



bool VarWindow::saveImage()
{
  return true;
  //   statusBar()->showMessage(tr("image printed"), 2000);
}



void VarWindow::help(const QString& name){
  QDomElement root = doc.documentElement();
  root = root.firstChildElement("help").firstChildElement(name);
  if(root.isNull())return;
  QMessageBox::about(this, name,
                     root.text());
}





#include <QDomDocument>
#include "mainwindow.h"
#include "grid.h"
#include "cutdialog.h"
#include "vmergewindow.h"
#include "fvmadapt.h"
#include "varwindow.h"
#include "importwindow.h"
#include "qualitydialog.h"
#include "progressdialog.h"

MyPushButton::MyPushButton(const QString& title, QWidget* parent):QPushButton(title, parent){
  QWidget::setMouseTracking(true);
}
void MyPushButton::mouseMoveEvent(QMouseEvent */*event*/)
{
  emit showText(whatsThis());
}
MainWindow::MainWindow(QWidget* parent):QWidget(parent)
{
  QWidget::setAttribute(Qt::WA_DeleteOnClose, true);
  QWidget::setMouseTracking(true);
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
  QDomDocument doc;
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

  QDomElement helpNode = doc.documentElement().firstChildElement("help");
  
   buttonGroup = new QButtonGroup(this);
  
  QHBoxLayout* buttonLayout = new QHBoxLayout;

  int count = 0; 
  for(QDomElement elt = helpNode.firstChildElement();
      !elt.isNull(); elt = elt.nextSiblingElement(), count++){
    
    MyPushButton* button = new MyPushButton(elt.hasAttribute("title")?elt.attribute("title"):elt.tagName());
    //helpInfo << elt.text();
    button->setCheckable(true);
    button->setWhatsThis(elt.text());
    buttonGroup->addButton(button, count);
    buttonLayout->addWidget(button);
    connect(button, SIGNAL(showText(const QString&)), this, SLOT(showText(const QString&)));
    
  }
  connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(processItem(int)));
  display = new QTextEdit;
  display->setText("Please select a module"); 
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addLayout(buttonLayout);
  mainLayout->addWidget(display);
  
  setLayout(mainLayout);
  
}

void MainWindow::showText(const QString& text){
  if(text.isEmpty()) display->setText("Please select a module");
  else display->setText(text);
}

void MainWindow::processItem(int index){
  QPushButton* button = qobject_cast<QPushButton*>(buttonGroup->button(index));
  if(!button->isChecked()) return;
  // display->setText(helpInfo[index]);
  
  QString item = button->text();
  if ( !item.isEmpty()){
    if(item=="Import"){
      import();
    }else if(item=="VogCheck"){
      vcheck();
   }else if(item=="VogMerge"){
      vmerge();
    }else if(item=="FVMadapt"){
      fvmAdapt();
    }else if(item=="Post-Processing"){
      postProcess();
    }else{
      generateVar();
    }
  }
}



void MainWindow::generateVar(){
  VarWindow* varWindow = new VarWindow();
  varWindow->show();
}


void MainWindow::postProcess(){
  

   
 
  QString fileName =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                 QDir::currentPath(),
                                 tr("vog Files (*.vog)"));

  if(fileName=="")return;
  LoadInfo ldinfo;
  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
  ldinfo.casename = casename;
  ldinfo.directory = directory;
  
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
   CutDialog* cutdialog = new CutDialog(ldinfo);
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
  QString directory = importFileName.section('/', 0, -2)+"/";
  
  QFile exist_test(importFileName+tr(".vog"));
  if(!(exist_test.exists())){
    QMessageBox::warning(window(), tr("vogcheck"),
                         tr("Please convert the file to volume grid format first")
                         );
    return;
  }
  QString out_filename="./output/check_"+casename+".out";
  QString command2 = "vogcheck " +casename;
  ProgressDialog* progress = new ProgressDialog(command2, directory);
  connect(progress, SIGNAL(progressFinished(QString, QProcess::ExitStatus, QString)), this, SLOT(showQuality(QString, QProcess::ExitStatus, QString)));
  progress->show();
}
void MainWindow::showQuality(QString command, QProcess::ExitStatus status, QString directory){
  if(status==QProcess::CrashExit)return;
  
  QString filename = directory+command.section(' ',-1, -1)+".quality";
  QualityDialog qualityDialog(filename, this);
  qualityDialog.exec();

  //the following block is for post-processing
 //  if(waitForQualityFile){
//     waitForQualityFile = false;
//     QDomElement elem = doc.documentElement().firstChildElement("mainWindow");
//     elem = elem.firstChildElement("gridSetup");
//     LoadInfo ldinfo;
//     ldinfo.casename = elem.attribute("casename");
//     ldinfo.directory = elem.attribute("directory");
//     CutDialog* cutdialog = new CutDialog(ldinfo, viewer->boundaryBoxSize(), viewer);
//     cutdialog->show();
    
//   }
}
  


  

void MainWindow::fvmAdapt()
{
  QString fileName =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                 QDir::currentPath(),
                                 tr("vog Files (*.vog)"));
  FVMAdapt* adaptwindow = new FVMAdapt(fileName);
  adaptwindow->show();

}

void MainWindow::import()
{
  char* resourcepath = getenv("CHEMDEMOPATH");
  QString xmlpath = "./xml/";
  QString pngpath = "./png/";
  
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
  QDomDocument doc;
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
   
    
  QDomElement theroot = doc.documentElement();
  QDomElement elem = theroot.firstChildElement("importPage");
  if(elem.attribute("element")=="importWindow"){
    ImportWindow* newWindow=new ImportWindow(elem);
    newWindow->show();
   
  }
}


void MainWindow::vmerge()
{

  VMergeWindow*  vmwindow = new VMergeWindow();
  vmwindow->show();
 
}


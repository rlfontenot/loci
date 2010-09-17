/////////////////////////////////////////////////
//  Filename: cutdialog.cpp
//
//  Contains: Implementation of CutDialog class
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
#include <QToolBar>
#include <QDockWidget>
#include "qualitydialog.h"
#include "progressdialog.h"
#include "cutdialog.h"
#include "grid.h"
#include "pages.h"
#include "helpwindow.h"
//////////////////////////////////////////////////////////
//  public:
//    CutDialog(dialog_info *info, QWidget *parent = 0);
//
//  Assembles the dialog shown by the "load grid " in 'file' menu.
//////////////////////////////////////////////////////////

CutDialog::CutDialog(QWidget *parent):QMainWindow(parent)
{

  setWindowTitle("postprocessing");
  
  setAttribute(Qt::WA_DeleteOnClose, true);
  QGroupBox* central= new QGroupBox;
  central->setFlat(true);
  viewer = new GLViewer;
  
  
  size = 1.0; 
  
   buttonGroup = new QButtonGroup(this);
  buttonGroup->setExclusive(true);
  QPushButton* xy_button = new QPushButton(tr("xy_plane"));
  QPushButton* yz_button = new QPushButton(tr("yz_plane"));
  QPushButton* xz_button = new QPushButton(tr("xz_plane"));
  buttonGroup->addButton(xy_button);
  buttonGroup->addButton(yz_button);
  buttonGroup->addButton(xz_button);
  xy_button->setCheckable(true);
  yz_button->setCheckable(true);
  xz_button->setCheckable(true);
  xy_button->setChecked(true);
  buttonGroup->setId(xy_button, 1);
  buttonGroup->setId(yz_button, 2);
  buttonGroup->setId(xz_button, 3);
  QHBoxLayout* hlayout= new QHBoxLayout;
  hlayout->addWidget(xy_button);
  hlayout->addWidget(yz_button);
  hlayout->addWidget(xz_button);
  connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(planeSelected(int)));
  
  

  
  QGroupBox* translateBox = new QGroupBox(tr("translation along normal")); 
  QGridLayout* translate = new QGridLayout;
  
  translate->addWidget(new QLabel(tr("  ")), 0, 0);
  
  
  
  zslider1 = new QSlider;
  zslider1->setMinimum(-1000);
  zslider1->setMaximum(1000);
  zslider1->setValue(0);
  zslider1->setOrientation(Qt::Horizontal);
  
  zEditor1 = new DoubleEdit(0.0); 
  zEditor1->setRange(-size/2.0, size/2.0);
  translate->addWidget(zslider1, 0, 1);
  translate->addWidget(zEditor1, 0, 2);
  connect(zslider1, SIGNAL(valueChanged(int)), zEditor1, SLOT(mapValue(int))); 
 
  translateBox->setLayout(translate);
  
  
  rotateBox = new QGroupBox(tr("rotation")); 
  QGridLayout* rotate = new QGridLayout;
  rotateBox->setCheckable(true);
  rotateBox->setChecked(true);
  rotate->addWidget(new QLabel(tr("x:")), 0, 0);
  rotate->addWidget(new QLabel(tr("y:")), 1, 0);
  rotate->addWidget(new QLabel(tr("z:")), 2, 0);
  
  xslider2 = new QSlider;
  yslider2 = new QSlider;
  zslider2 = new QSlider;
  xslider2->setMinimum(-1000);
  xslider2->setMaximum(1000);
  yslider2->setMinimum(-1000);
  yslider2->setMaximum(1000);
  zslider2->setMinimum(-1000);
  zslider2->setMaximum(1000);
  xslider2->setValue(0);
  yslider2->setValue(0);
  zslider2->setValue(0);
  xslider2->setOrientation(Qt::Horizontal);
  yslider2->setOrientation(Qt::Horizontal);
  zslider2->setOrientation(Qt::Horizontal);
  xEditor2 = new DoubleEdit(0.0); 
  yEditor2 = new DoubleEdit(0.0);
  zEditor2 = new DoubleEdit(0.0); 
  xEditor2->setRange(-90, 90);
  yEditor2->setRange(-90, 90);
  zEditor2->setRange(-90, 90);
  rotate->addWidget(xslider2, 0, 1);
  rotate->addWidget(xEditor2, 0, 2);
  rotate->addWidget(yslider2, 1, 1);
  rotate->addWidget(yEditor2, 1, 2);
  rotate->addWidget(zslider2, 2, 1);
  rotate->addWidget(zEditor2, 2, 2);
  connect(xslider2, SIGNAL(valueChanged(int)), xEditor2, SLOT(mapValue(int)));
  connect(yslider2, SIGNAL(valueChanged(int)), yEditor2, SLOT(mapValue(int)));
  connect(zslider2, SIGNAL(valueChanged(int)), zEditor2, SLOT(mapValue(int))); 
  rotateBox->setLayout(rotate);
  
  connect(zEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(xEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zEditor2, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  

   // Casename group
  QGroupBox *caseGroup = new QGroupBox(tr("Casename"));
  QHBoxLayout *caseLayout = new QHBoxLayout;
  caseLabel = new QLabel(ld_info.casename);
  caseLayout->addWidget(caseLabel);
  caseGroup->setLayout(caseLayout);
  
  
  // Iteration and variable group
  QGroupBox *iterVarGroup = new QGroupBox(tr("Select scalar file: "));
  QFormLayout *iterVarLayout = new QFormLayout;
  comboIter = new QComboBox;
  comboVar = new QComboBox;
  iterVarLayout->addRow(new QLabel(tr("Iteration:")), comboIter);
  iterVarLayout->addRow(new QLabel(tr("Variable:")), comboVar);
  iterVarGroup->setLayout(iterVarLayout);


  connect(comboIter, SIGNAL(currentIndexChanged(QString)),
	  this, SLOT(updateVars(QString)));


  
  
  //updateCase();
 





  // 'Cancel' & 'Okay' buttons
  QHBoxLayout *cutButtons = new QHBoxLayout;
  QPushButton *cancel = new QPushButton(tr("reset cutplane"));
 
 
  
 
  cancel->setDefault(false);


  
  


  cutButtons->addWidget(cancel);
  connect(cancel, SIGNAL(clicked()),
	  this, SLOT(reset()));

 
  
 
  
  QGroupBox* cutGroup = new QGroupBox(tr("Define cut plane:"));
  QVBoxLayout* cutLayout = new QVBoxLayout;
  cutLayout->addLayout(hlayout);
  cutLayout->addWidget(rotateBox);
  cutLayout->addWidget(translateBox);
  
  cutLayout->addLayout(cutButtons);
  cutGroup->setLayout(cutLayout);


  QSlider* extrSlider = new QSlider;
  extrSlider->setMinimum(-1000);
  extrSlider->setMaximum(1000);
  extrSlider->setOrientation(Qt::Horizontal);
  extrEdit  = new DoubleEdit;
  extrEdit->setRange(0.0, 1.0);
  extrEdit->setValue(0.0);
  connect(extrSlider, SIGNAL(valueChanged(int)), extrEdit, SLOT(mapValue(int)));
  QPushButton* clearExtrButton = new QPushButton("Clear");
  connect(clearExtrButton, SIGNAL(clicked()), viewer, SLOT(clearExtrema()));

  QHBoxLayout *buttons = new QHBoxLayout;
  QGroupBox* buttonsGroup = new QGroupBox(tr("Extreme nodes"));
  
  
  buttons->addWidget(extrSlider);
  buttons->addWidget(extrEdit);
  buttons->addWidget(clearExtrButton); 
  buttonsGroup->setLayout(buttons);
 
  connect(extrEdit, SIGNAL(valueChanged(double)), viewer, SLOT(setExtrema(double)));
  
  toolbar = new QToolBar(tr("display"), this);
  addToolBar(toolbar);
  createToolBar();
  createFlowBar();
  createVisBar();

  QVBoxLayout *mainLayout =  new QVBoxLayout;
  mainLayout->addWidget(caseGroup);
  mainLayout->addWidget(iterVarGroup);
  mainLayout->addWidget(buttonsGroup);
  mainLayout->addWidget(cutGroup);

  QDockWidget*  viewerDock  = new QDockWidget("Graphics Viewer", this); 
  viewerDock->setAllowedAreas(Qt::RightDockWidgetArea );
  viewerDock->setWidget(viewer);
  addDockWidget(Qt::RightDockWidgetArea, viewerDock);

  central->setLayout(mainLayout);
  setCentralWidget(central);

}
void CutDialog::createVisBar()
{

  
    
  addToolBarBreak();
  QGroupBox* visbar = new QGroupBox(tr("visualization"));
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

   QPushButton* helpButton = new QPushButton(tr("help"));
  connect(helpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  visLayout->addWidget(helpButton);
  
  visbar->setLayout(visLayout);
  
  //  toolbar->addStretch(20);

  toolbar->addWidget(visbar);
  
} 

void CutDialog::updateCase(){

 // If a case name has previously been selected, make default
  if (!ld_info.casename.isEmpty() && !ld_info.directory.isEmpty()) {
    caseLabel->setText("    " + ld_info.casename);

    QDir dir(ld_info.directory+"/output/");
    QStringList filters;
    filters << "grid_pos.*_" + ld_info.casename;
    QStringList gridposFiles = dir.entryList(filters);
   
    
    // Load in valid iteration numbers (variables loaded automatically)
    for (int i = 0; i < gridposFiles.size(); ++i) {
      QString gridpos = gridposFiles.at(i);
      gridpos.remove(0, 9);
      gridpos.remove(gridpos.size() - ld_info.casename.size() - 1, 
		     ld_info.casename.size() + 1);
      comboIter->addItem(gridpos);
    }

    // If iteration previously selected, make default
    if (!ld_info.iteration.isEmpty()) {
      comboIter->setCurrentIndex(comboIter->findText(ld_info.iteration));

      // If variable previously selected, make default
      if (!ld_info.variable.isEmpty())
	comboVar->setCurrentIndex(comboVar->findText(ld_info.variable));
    }
  }

}
void CutDialog::updateSize(){
     zEditor1->setRange(-size/2.0, size/2.0);
}

void CutDialog::createFlowBar(){
  int spacing =2;  
  //create flowbar
 
  
  QGroupBox* flowbar = new QGroupBox("flow bar");
  
  QVBoxLayout* barLayout = new QVBoxLayout;

  barLayout->addSpacing(spacing);
  QPushButton* loadButton = new QPushButton(tr("Load Grid"), this);
  barLayout->addWidget(loadButton);
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadGrid()));
  
  barLayout->addSpacing(spacing);
  QPushButton* loadScaButton = new QPushButton(tr("Load Scalar Value"), this);
  barLayout->addWidget(loadScaButton);
  connect(loadScaButton, SIGNAL(clicked()), this, SLOT(loadSca()));


  QPushButton *cutButton = new QPushButton(tr("Cut"));
  barLayout->addWidget(cutButton);
  connect(cutButton, SIGNAL(clicked()), this, SLOT(cut()));
  
 
  
  barLayout->addSpacing(spacing);
  QPushButton* doneButton = new QPushButton(tr("Done"), this);
  barLayout->addWidget(doneButton);
  connect(doneButton, SIGNAL(clicked()), this, SLOT(close()));

  barLayout->addStretch(10);
  flowbar->setLayout(barLayout);
  
  
  QToolBar* flowToolBar = new QToolBar;
  addToolBar(Qt::LeftToolBarArea,flowToolBar );
  flowToolBar->addWidget(flowbar);
}


void CutDialog::loadGrid(){
  //get file name
  QString fileName =
    QFileDialog::getOpenFileName(this, tr("Get File"),
                                 QDir::currentPath(),
                                 tr("vog Files (*.vog)"));

  if(fileName=="")return;
  //LoadInfo ldinfo;
  int first= fileName.lastIndexOf('/');
  int last = fileName.lastIndexOf('.');
  QString casename = fileName.mid(first+1, last-first-1);
  QString directory = fileName.left(first);
  
  QDir dir(directory+"/output/");
  QStringList filters;
  filters << "grid_pos.*_" + casename;
  QStringList gridposFiles = dir.entryList(filters);
  if(gridposFiles.size()==0){
   
    int ret = QMessageBox::question(this, "post-processing",
                                    tr("No scalar value, do you want to run vogcheck? "),
                                    QMessageBox::Ok | QMessageBox::Cancel);
    switch(ret){
    case QMessageBox::Ok:

      check(directory+'/'+casename);
      break;
    default:
      ld_info.casename = casename;
      ld_info.directory = directory;
      //QString fileName = ld_info.directory+"/"+ld_info.casename+".vog";
      QStringList bnames;
      if(!(viewer->load_boundary(fileName, bnames))){
        QMessageBox::information(window(), tr("post-processing"),
                                 tr("can not load grid ")+fileName);
        return;
      }
      size = viewer->boundaryBoxSize();
      updateSize();
      updateCase();
      
      return;
    }
  }else{
    ld_info.casename = casename;
    ld_info.directory = directory;
    //QString fileName = ld_info.directory+"/"+ld_info.casename+".vog";
    QStringList bnames;
    if(!(viewer->load_boundary(fileName, bnames))){
      QMessageBox::information(window(), tr("post-processing"),
                               tr("can not load grid ")+fileName);
      return;
    }
    size = viewer->boundaryBoxSize();
    updateSize();
    updateCase();

  }
}



  void CutDialog::check(const QString& fn){
    
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
  void CutDialog::showQuality(QString command, QProcess::ExitStatus status, QString directory){
    if(status==QProcess::CrashExit)return;
    
    QString filename = directory+command.section(' ',-1, -1)+".quality";
    QualityDialog qualityDialog(filename, this);
    qualityDialog.exec();

    //the following block is for post-processing

      //LoadInfo ldinfo;
    ld_info.casename =command.section(' ',-1, -1);
    ld_info.directory = directory;

    QString fileName = ld_info.directory+"/"+ld_info.casename+".vog";
    QStringList bnames;
    if(!(viewer->load_boundary(fileName, bnames))){
      QMessageBox::information(window(), tr("post-processing"),
                                 tr("can not load grid ")+fileName);
        return;
    }
    size = viewer->boundaryBoxSize();
    updateSize();

    updateCase();

  }
  

  
void CutDialog::createToolBar(){
  QPushButton *showBoundaryAct = new QPushButton(tr("Show\nBoundaries"), this);
  QPushButton *showBShadingAct = new QPushButton(tr("Boundary\nShading"), this);
  QPushButton *showContoursAct = new QPushButton(tr("Contours"), this);
  QPushButton *showGridAct = new QPushButton(tr("Show\nGrid"), this);
  QPushButton *showShadingAct = new QPushButton(tr("Cutplane\nShading"), this);
  
  QPushButton *shadeType1 = new QPushButton(tr("Blue\nto Red"), this);
  QPushButton *shadeType2 = new QPushButton(tr("Black\nBody"), this);
  QPushButton *shadeType3 = new QPushButton(tr("Pressure\n"), this);
  
  QGroupBox* bdryBox = new QGroupBox(tr("boundary display"));
  QHBoxLayout* bdryLayout = new QHBoxLayout;
  bdryLayout->addWidget(showBoundaryAct);
  bdryLayout->addWidget(showBShadingAct);
  bdryBox->setLayout(bdryLayout);
  
  QGroupBox* ctplBox = new QGroupBox(tr("cutplane display"));
  QHBoxLayout* ctplLayout = new QHBoxLayout;
 
  ctplLayout->addWidget(showShadingAct);
  ctplLayout->addWidget(showGridAct);
   QVBoxLayout* vlayout = new QVBoxLayout;
  vlayout->addWidget(showContoursAct);
  QSlider* slider = new QSlider(Qt::Horizontal);
  slider->setRange(5, 50);
  slider->setValue(10);
  vlayout->addWidget(slider);

  ctplLayout->addLayout(vlayout);
  ctplBox->setLayout(ctplLayout);
  
  QGroupBox* typeBox = new QGroupBox("color mapping");
  QHBoxLayout* typeLayout = new QHBoxLayout;
  typeLayout->addWidget(shadeType1);
  typeLayout->addWidget(shadeType2);
  typeLayout->addWidget(shadeType3);
  typeBox->setLayout(typeLayout);
    
 
  toolbar->addWidget(bdryBox);
  toolbar->addSeparator();
  toolbar->addWidget(ctplBox);
  toolbar->addSeparator();
  toolbar->addWidget(typeBox);

 //  QGroupBox* emptyBox = new QGroupBox();
//   emptyBox->setFlat(true);
//   QHBoxLayout* emptyLayout = new QHBoxLayout;
//   emptyLayout->addStretch(2);
//   emptyBox->setLayout(emptyLayout);
//   toolbar->addWidget(emptyBox);
 
  
 
 

  
  connect(showBoundaryAct,SIGNAL(clicked()),
          viewer, SLOT(showBoundaries())); 
  connect(showBShadingAct, SIGNAL(clicked()),
          viewer, SLOT(toggleBoundaryShading()));
  connect(showGridAct, SIGNAL(clicked()),
          viewer, SLOT(toggleGrid()));
  connect(showContoursAct, SIGNAL(clicked()),
          viewer, SLOT(toggleContours()));
  connect(showShadingAct, SIGNAL(clicked()),
          viewer, SLOT(toggleShading()));
  connect(shadeType1, SIGNAL(clicked()),
          viewer, SLOT(setShadeType1()));
  connect(shadeType2, SIGNAL(clicked()),
          viewer, SLOT(setShadeType2()));
  connect(shadeType3, SIGNAL(clicked()),
          viewer, SLOT(setShadeType3()));
  connect(slider, SIGNAL(valueChanged(int)),
       viewer, SLOT(changeContours(int)));
 
  connect(this, SIGNAL(destroyed()), viewer, SLOT(uncut()));
}


 void CutDialog::helpClicked(){
   HelpWindow* helpwindow = new HelpWindow("page_postprocess.html");
   helpwindow->show();
 }     
 

// void CutDialog::showExtremeNodes(double value){
//   viewer->setExtrema();
  
  
//   }

//////////////////////////////////////////////////////////////////////////////
//  private slots:
//    void setInfo();
//
//  Called when the 'Okay' button is pressed. Updates the data structure passed
//  from the main window with all the information in the widgets before
//  destroying the dialog.
//////////////////////////////////////////////////////////////////////////////


void CutDialog::setInfo()
{
  positions3d translate =  positions3d(0, 0, -zEditor1->value());
  positions3d rotate = positions3d(-xEditor2->value(), -yEditor2->value(), -zEditor2->value());
  info.rotate = rotate;
  info.translate = translate;
  
  
  ld_info.variable = comboVar->currentText();
  ld_info.iteration = comboIter->currentText();
  
  
  viewer->previewCut(info);  
  viewer->setLoadInfo(ld_info);  
 
}

void CutDialog::cut(){
  setInfo();
  viewer->cut();
}

void CutDialog::loadSca(){
  setInfo();
 
  viewer->loadSca();
  extrEdit->setBottom(viewer->get_min_val());
  extrEdit->setTop(viewer->get_max_val());
  // extraSlider->setValue();
  extrEdit->setValue(viewer->get_mid_val());
}

void CutDialog::reset(){
  buttonGroup->button(1)->setChecked(true);
  zslider1->setValue(0);
  
  xslider2->setValue(0);
  yslider2->setValue(0);
  zslider2->setValue(0);
  

  zEditor1->setValue(0.0); 
  xEditor2->setValue(0.0);
  yEditor2->setValue(0.0);
  zEditor2->setValue(0.0);
  rotateBox->setChecked(true);
  setInfo();
}


void CutDialog::planeSelected(int id){
  switch(id){
  case 1:
    xslider2->setValue(0);
    yslider2->setValue(0);
    zslider2->setValue(0);
    
    break;
  case 2:
    
    rotateBox->setChecked(false);
    xslider2->setValue(0);
    yslider2->setValue(1000);
    zslider2->setValue(0);
    break;
  case 3:
    
    xslider2->setValue(-1000);
    yslider2->setValue(0);
    zslider2->setValue(0); 
    rotateBox->setChecked(false);
    break;
  default:
    qDebug()<<"invalid id in planeSelected()";
  }

 
  
  setInfo();
}
    


/////////////////////////////////////////////////////////////////////////////
//  private slots:
//    void updateVars(QString iter);
//
//  Called when the iteration number is changed. Finds all possible variable
//  names under chosen iteration.
/////////////////////////////////////////////////////////////////////////////

void CutDialog::updateVars(QString iter)
{
  comboVar->clear();

  // Find all possible variable names under iteration
  ld_info.iteration = iter;
  QDir dir(ld_info.directory+"/output/");
  QStringList filters;
  filters << "*_sca." + iter + "_" + ld_info.casename;
  QStringList varFiles = dir.entryList(filters);

  // Load possible variable names into comboVar
  int end = 6 + iter.size() + ld_info.casename.size();
  for (int i = 0; i < varFiles.size(); ++i) {
    QString var = varFiles.at(i);
    var.remove(var.size() - end, end);
    comboVar->addItem(var);
  }

}

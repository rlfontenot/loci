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
#include "cutdialog.h"
#include "grid.h"
#include "pages.h"
//////////////////////////////////////////////////////////
//  public:
//    CutDialog(dialog_info *info, QWidget *parent = 0);
//
//  Assembles the dialog shown by the "load grid " in 'file' menu.
//////////////////////////////////////////////////////////

CutDialog::CutDialog(LoadInfo ldinfo,  QWidget *parent):QMainWindow(parent),ld_info(ldinfo)
{


  setAttribute(Qt::WA_DeleteOnClose, true);
  

  QGroupBox* central= new QGroupBox;
  central->setFlat(true);
  

  viewer = new GLViewer;
  QString fileName = ld_info.directory+"/"+ld_info.casename+".vog";
  QStringList bnames;
  if(!(viewer->load_boundary(fileName, bnames))){
      QMessageBox::information(window(), tr("post-processing"),
                               tr("can not load grid ")+fileName);
      return;
    }
  size = viewer->boundaryBoxSize();
  
  QButtonGroup* buttonGroup = new QButtonGroup(this);
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
  
  

  
  QGroupBox* translateBox = new QGroupBox(tr("translation")); 
  QGridLayout* translate = new QGridLayout;
  
  translate->addWidget(new QLabel(tr("x:")), 0, 0);
  translate->addWidget(new QLabel(tr("y:")), 1, 0);
  translate->addWidget(new QLabel(tr("z:")), 2, 0);
  
  xslider1 = new QSlider;
  yslider1 = new QSlider;
  zslider1 = new QSlider;
  xslider1->setMinimum(-1000);
  xslider1->setMaximum(1000);
  yslider1->setMinimum(-1000);
  yslider1->setMaximum(1000);
  zslider1->setMinimum(-1000);
  zslider1->setMaximum(1000);
  xslider1->setValue(0);
  yslider1->setValue(0);
  zslider1->setValue(0);
  xslider1->setOrientation(Qt::Horizontal);
  yslider1->setOrientation(Qt::Horizontal);
  zslider1->setOrientation(Qt::Horizontal);
  
  xEditor1 = new DoubleEdit(0.0); 
  yEditor1 = new DoubleEdit(0.0);
  zEditor1 = new DoubleEdit(0.0); 
  xEditor1->setRange(-size/2.0, size/2.0);
  yEditor1->setRange(-size/2.0, size/2.0);
  zEditor1->setRange(-size/2.0, size/2.0);
  translate->addWidget(xslider1, 0, 1);
  translate->addWidget(xEditor1, 0, 2);
  translate->addWidget(yslider1, 1, 1);
  translate->addWidget(yEditor1, 1, 2);
  translate->addWidget(zslider1, 2, 1);
  translate->addWidget(zEditor1, 2, 2);
  connect(xslider1, SIGNAL(valueChanged(int)), xEditor1, SLOT(mapValue(int)));
  connect(yslider1, SIGNAL(valueChanged(int)), yEditor1, SLOT(mapValue(int)));
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
  
  connect(xEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yEditor1, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
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
  QGroupBox *iterVarGroup = new QGroupBox(tr("Iteration and variables"));
  QFormLayout *iterVarLayout = new QFormLayout;
  comboIter = new QComboBox;
  comboVar = new QComboBox;
  iterVarLayout->addRow(new QLabel(tr("Iteration:")), comboIter);
  iterVarLayout->addRow(new QLabel(tr("Variable:")), comboVar);
  iterVarGroup->setLayout(iterVarLayout);


  connect(comboIter, SIGNAL(currentIndexChanged(QString)),
	  this, SLOT(updateVars(QString)));


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







  // 'Cancel' & 'Okay' buttons
  QHBoxLayout *cutButtons = new QHBoxLayout;
  QPushButton *cancel = new QPushButton(tr("Cancel"));
  QPushButton *cut = new QPushButton(tr("Cut"));
  QPushButton *close = new QPushButton(tr("Close"));
  QPushButton *load = new QPushButton(tr("Load Scalar file"));
  
  extrSpinBox = new QSpinBox;
  extrSpinBox->setRange(0, 50);
  extrSpinBox->setValue(0);
  
  
  cancel->setDefault(false);
  close->setDefault(false);
  cut->setDefault(false);
  
  
  QHBoxLayout *buttons = new QHBoxLayout;
  
  QGroupBox* buttonsGroup = new QGroupBox;
  
  buttons->addWidget(load);
  
  buttons->addWidget(new QLabel(tr("Show Extreme Nodes:  percentage")));
  buttons->addWidget(extrSpinBox);
  
  buttonsGroup->setLayout(buttons);
  
  cutButtons->addWidget(cut);
  cutButtons->addWidget(cancel);
  buttons->addWidget(close);
 
  
  

  connect(cancel, SIGNAL(clicked()),
	  this, SLOT(reset()));
  connect(cut, SIGNAL(clicked()), this, SLOT(cut()));
  connect(close, SIGNAL(clicked()),
	  this, SLOT(close()));
  connect(load, SIGNAL(clicked()), this, SLOT(load()));
  connect(extrSpinBox, SIGNAL(valueChanged(int)), viewer, SLOT(setPercentage(int)));
  connect(extrSpinBox, SIGNAL(valueChanged(int)), this, SLOT(showExtremeNodes(int)));
  
  QGroupBox* cutGroup = new QGroupBox(tr("define cut plane:"));
  QVBoxLayout* cutLayout = new QVBoxLayout;
  cutLayout->addLayout(hlayout);
  cutLayout->addWidget(translateBox);
  cutLayout->addWidget(rotateBox);
  cutLayout->addLayout(cutButtons);
  cutGroup->setLayout(cutLayout);
 
  createToolBar();  
  QVBoxLayout *mainLayout =  new QVBoxLayout;
 
  mainLayout->addWidget(buttonsGroup);
  
  mainLayout->addWidget(caseGroup);
  mainLayout->addWidget(iterVarGroup);
 
  mainLayout->addWidget(cutGroup);
  //mainLayout->addWidget(viewer);
  QDockWidget*  viewerDock  = new QDockWidget("Graphics Viewer", this); 
  viewerDock->setAllowedAreas(Qt::RightDockWidgetArea );
  viewerDock->setWidget(viewer);
  addDockWidget(Qt::RightDockWidgetArea, viewerDock);

  central->setLayout(mainLayout);
  setCentralWidget(central);
  setInfo();
  
  
}
void CutDialog::createToolBar(){
  QPushButton *showBoundaryAct = new QPushButton(tr("Show Boundaries"), this);
  QPushButton *showBShadingAct = new QPushButton(tr("Boundary Shading"), this);
  QPushButton *showContoursAct = new QPushButton(tr("Show Contours"), this);
  QPushButton *showGridAct = new QPushButton(tr("Show Grid"), this);
  QPushButton *showShadingAct = new QPushButton(tr("Cutplane Shading"), this);
  
  QPushButton *shadeType1 = new QPushButton(tr("Blue to Red"), this);
  QPushButton *shadeType2 = new QPushButton(tr("Blackbody"), this);
  QPushButton *shadeType3 = new QPushButton(tr("Pressure"), this);

  // QPushButton *showBorderAct = new QPushButton(tr("Show Outline"), this);
 
  
  
  
  
  toolbox = new QGroupBox(tr("Display Menu"));
  QHBoxLayout* barlayout = new QHBoxLayout;
  barlayout->addWidget(showBoundaryAct);
  barlayout->addSpacing(5);;
  barlayout->addWidget(showBShadingAct);
 
 
  barlayout->addSpacing(5);;
  barlayout->addWidget(showGridAct);
  barlayout->addSpacing(5);;
  barlayout->addWidget(showShadingAct);
    
  barlayout->addSpacing(5);;
  barlayout->addWidget(shadeType1);
  barlayout->addWidget(shadeType2);
  barlayout->addWidget(shadeType3);
  

    
  barlayout->addSpacing(5);;
 
  
  

    
  barlayout->addSpacing(5);;
  barlayout->addWidget(showContoursAct);
  QSlider* slider = new QSlider(Qt::Horizontal);
  slider->setRange(5, 50);
  slider->setValue(10);
  barlayout->addWidget(slider);


  toolbox->setLayout(barlayout);

  
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
  
  //connect(showBorderAct, SIGNAL(clicked()),
  //      viewer, SLOT(toggleBorder()));
  
  connect(slider, SIGNAL(valueChanged(int)),
       viewer, SLOT(changeContours(int)));
 
  connect(this, SIGNAL(destroyed()), viewer, SLOT(uncut()));
  
  QToolBar* toolbar = addToolBar(tr("display"));
  toolbar->addWidget(toolbox);
 
}



 

void CutDialog::showExtremeNodes(int i){
  viewer->setShading(i==0);
  }

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
  positions3d translate =  positions3d(-xEditor1->value(), -yEditor1->value(), -zEditor1->value());
  positions3d rotate = positions3d(-xEditor2->value(), -yEditor2->value(), -zEditor2->value());
  info.rotate = rotate;
  info.translate = translate;
  
  
  ld_info.variable = comboVar->currentText();
  ld_info.iteration = comboIter->currentText();
  
  if(ld_info.variable=="cellVol")extrSpinBox->setRange(-50, 0);
  else extrSpinBox->setRange(0, 50);

  viewer->previewCut(info);  
  viewer->setLoadInfo(ld_info);  
 
}

void CutDialog::cut(){
  setInfo();
  viewer->cut();
}

void CutDialog::load(){
  setInfo();
  extrSpinBox->setValue(0);
  viewer->loadSca();
}

void CutDialog::reset(){
  xslider1->setValue(0);
  yslider1->setValue(0);
  zslider1->setValue(0);
  
  xslider2->setValue(0);
  yslider2->setValue(0);
  zslider2->setValue(0);
  
  xEditor1->setValue(0.0);
  yEditor1->setValue(0.0);
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
    xEditor2->setValue(0.0);
    yEditor2->setValue(0.0);
    zEditor2->setValue(0.0);
    rotateBox->setChecked(false);
    break;
  case 2:
    
    rotateBox->setChecked(false);
    xEditor2->setValue(0.0);
    yEditor2->setValue(-90.0);
    zEditor2->setValue(-90.0);
    break;
  case 3:
    xEditor2->setValue(-90.0);
    yEditor2->setValue(0.0);
    zEditor2->setValue(0.0); 
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

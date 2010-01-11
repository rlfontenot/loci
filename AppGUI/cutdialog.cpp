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
#include "cutdialog.h"
#include "grid.h"
#include "pages.h"
//////////////////////////////////////////////////////////
//  public:
//    CutDialog(dialog_info *info, QWidget *parent = 0);
//
//  Assembles the dialog shown by the "load grid " in 'file' menu.
//////////////////////////////////////////////////////////

CutDialog::CutDialog( LoadInfo ldInfo, float size, QWidget *parent)
  :  QWidget(parent),ld_info(ldInfo)
{
 
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
  
  xEditor1 = new FloatEdit(0.0); 
  yEditor1 = new FloatEdit(0.0);
  zEditor1 = new FloatEdit(0.0); 
  xEditor1->setRange(-size/2.0, size/2.0);
  yEditor1->setRange(-size/2.0, size/2.0);
  zEditor1->setRange(-size/2.0, size/2.0);
  translate->addWidget(xslider1, 0, 1);
  translate->addWidget(xEditor1, 0, 2);
  translate->addWidget(yslider1, 1, 1);
  translate->addWidget(yEditor1, 1, 2);
  translate->addWidget(zslider1, 2, 1);
  translate->addWidget(zEditor1, 2, 2);
  connect(xslider1, SIGNAL(valueChanged(int)), xEditor1, SLOT(setValue(int)));
  connect(yslider1, SIGNAL(valueChanged(int)), yEditor1, SLOT(setValue(int)));
  connect(zslider1, SIGNAL(valueChanged(int)), zEditor1, SLOT(setValue(int))); 
 
  translateBox->setLayout(translate);
  
  
  rotateBox = new QGroupBox(tr("ratation")); 
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
  xEditor2 = new FloatEdit(0.0); 
  yEditor2 = new FloatEdit(0.0);
  zEditor2 = new FloatEdit(0.0); 
  xEditor2->setRange(-90, 90);
  yEditor2->setRange(-90, 90);
  zEditor2->setRange(-90, 90);
  rotate->addWidget(xslider2, 0, 1);
  rotate->addWidget(xEditor2, 0, 2);
  rotate->addWidget(yslider2, 1, 1);
  rotate->addWidget(yEditor2, 1, 2);
  rotate->addWidget(zslider2, 2, 1);
  rotate->addWidget(zEditor2, 2, 2);
  connect(xslider2, SIGNAL(valueChanged(int)), xEditor2, SLOT(setValue(int)));
  connect(yslider2, SIGNAL(valueChanged(int)), yEditor2, SLOT(setValue(int)));
  connect(zslider2, SIGNAL(valueChanged(int)), zEditor2, SLOT(setValue(int))); 
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
  QHBoxLayout *buttons = new QHBoxLayout;
  QPushButton *cancel = new QPushButton(tr("Cancel"));
  QPushButton *cut = new QPushButton(tr("Cut"));
  QPushButton *close = new QPushButton(tr("Close"));
  cancel->setDefault(false);
  close->setDefault(false);
  cut->setDefault(false);
  buttons->addWidget(cut);
  buttons->addWidget(cancel);
  buttons->addWidget(close);

  connect(cancel, SIGNAL(clicked()),
	  this, SLOT(reset()));
  connect(cut, SIGNAL(clicked()), this, SLOT(cut()));
  connect(close, SIGNAL(clicked()),
	  this, SLOT(close()));

  setInfo();
  QVBoxLayout *mainLayout =  new QVBoxLayout;
  mainLayout->addLayout(hlayout);
  mainLayout->addWidget(translateBox);
  mainLayout->addWidget(rotateBox);
  mainLayout->addWidget(caseGroup);
  mainLayout->addWidget(iterVarGroup);
  mainLayout->addLayout(buttons);

  setLayout(mainLayout);
 

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
  emit cutInfoChanged(info);
  emit loadInfoChanged(ld_info);
}

void CutDialog::cut(){
  setInfo();
  emit cutPressed();
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

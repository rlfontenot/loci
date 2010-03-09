/////////////////////////////////////////////////
//  Filename: cutdialog.cpp
//
//  Contains: Implementation of RefDialog class
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
#include <QRadioButton>
#include <QTextEdit>
#include "refdialog.h"
#include "pages.h"
//////////////////////////////////////////////////////////
//  public:
//    RefDialog(dialog_info *info, QWidget *parent = 0);
//
//  Assembles the dialog shown by the "load grid " in 'file' menu.
//////////////////////////////////////////////////////////

RefDialog::RefDialog(QString bfile,  QWidget *parent)
  :  QWidget(parent),basefile(bfile)
{
  //  setAttribute(Qt::WA_DeleteOnClose, true);
  if(bfile=="")return;
  fileButtonGroup = new QButtonGroup(this);
  fileButtonGroup->setExclusive(false);
  
  QButtonGroup* optButtonGroup = new QButtonGroup(this);
  QRadioButton* xmlRadio = new QRadioButton(tr("use -xml option"));
  QRadioButton* tagRadio = new QRadioButton(tr("use -tag option"));
  optButtonGroup->addButton(xmlRadio, 0);
  optButtonGroup->addButton(tagRadio, 1);
  xmlRadio->setChecked(true);
  xmlOpt = true;
  connect(optButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(optChanged(int)));
  

  
  QLabel* baseLabel = new QLabel(tr("baseline grid File:"));
  QPushButton* baseButton = new QPushButton(basefile);
  
  QLabel* xmlLabel = new QLabel(tr("xml file:"));
  QPushButton* xmlButton = new QPushButton(tr("please specify filename"));
  
  QLabel* tagLabel = new QLabel(tr("tag file:"));
  QPushButton* tagButton = new QPushButton(tr("please specify filename"));
  
  QLabel* rPlanLabel = new QLabel(tr("plan file from last cycle:"));
  QPushButton* rPlanButton = new QPushButton(tr("please specify filename"));


  QLabel* tolLabel = new QLabel(tr("tolerance:"));
  FloatEdit* tolEdit = new FloatEdit;
  tol = 1e-5;
  tolEdit->setValue(tol);
  tolEdit->setBottom(0);
  connect(tolEdit, SIGNAL(valueChanged(double)), this, SLOT(setTol(double)));
  
  
  QLabel* levelLabel = new QLabel(tr("levels:"));
  IntEdit* levelEdit = new IntEdit;
  levels = 1;
  levelEdit->setValue(levels);
  levelEdit->setBottom(1);
  connect(levelEdit, SIGNAL(valueChanged(int)), this, SLOT(setLevels(int)));
  
  QLabel* foldLabel = new QLabel(tr("maximum face folding allowed:"));
  fold = 1.5708;
  FloatEdit* foldEdit = new FloatEdit(fold);
  connect(foldEdit, SIGNAL(valueChanged(double)), this, SLOT(setFold(double)));
  
  QGroupBox* modeGroup = new QGroupBox(tr("split mode:"));
  QVBoxLayout* modeLayout = new QVBoxLayout;
  QButtonGroup* modeButtonGroup = new QButtonGroup(this);
  modeButtonGroup->setExclusive(true);
  QRadioButton* modeRadio0 = new QRadioButton(tr("anisotropic split according to edge length"));
  
  modeButtonGroup->addButton(modeRadio0, 0);
  modeLayout->addWidget(modeRadio0);
  
  QRadioButton* modeRadio1 = new QRadioButton(tr("don't refine in z direction"));
  modeButtonGroup->addButton(modeRadio1, 1);
  modeLayout->addWidget(modeRadio1);

  QRadioButton* modeRadio2 = new QRadioButton(tr("fully isotropic split"));
  modeButtonGroup->addButton(modeRadio2, 2);
  modeLayout->addWidget(modeRadio2);
  modeRadio2->setChecked(true);
  modeGroup->setLayout(modeLayout);
  connect(modeButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(modeChanged(int)));
  mode = 2;

  QGroupBox* balanceGroup = new QGroupBox(tr("balance option:"));
  QVBoxLayout* balanceLayout = new QVBoxLayout;
  QButtonGroup* balanceButtonGroup = new QButtonGroup(this);
  balanceButtonGroup->setExclusive(true);

  QRadioButton* balanceRadio0 = new QRadioButton(tr("no edge's depth is greater than 1"));
  balanceButtonGroup->addButton(balanceRadio0, 0);
  balanceLayout->addWidget(balanceRadio0);
  
  QRadioButton* balanceRadio1 = new QRadioButton(tr("no edge's depth is greater than 1 \n and no cell has more than half of its faces split"));
  balanceButtonGroup->addButton(balanceRadio1, 1);
  balanceLayout->addWidget(balanceRadio1);
  
  QRadioButton* balanceRadio2 = new QRadioButton(tr("no edge's depth is greater than 1 \n and no cell has both the opposite faces split"));
  balanceButtonGroup->addButton(balanceRadio2, 2);
  balanceLayout->addWidget(balanceRadio2);
  balanceRadio0->setChecked(true);
  balanceGroup->setLayout(balanceLayout);
  connect(balanceButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(balanceChanged(int)));
  balanceOpt = 0;
  
  QLabel* outLabel = new QLabel(tr("output File:"));
  QPushButton* outButton = new QPushButton(tr("please specify filename"));
 
  QLabel* planLabel = new QLabel(tr("plan File:"));
  QPushButton* planButton = new QPushButton(tr("please specify filename"));
  
  fileButtonGroup->addButton(xmlButton, 0);
  fileButtonGroup->addButton(tagButton, 1);
  fileButtonGroup->addButton(rPlanButton, 2);
  fileButtonGroup->addButton(planButton, 3);
  
  fileButtonGroup->addButton(outButton, 4);
  connect(fileButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(getfile(int)));

  display = new QTextEdit("marker             \nrefmesh             ");
 
  QPushButton* saveButton = new QPushButton(tr("save script"));
  QPushButton* runButton = new QPushButton(tr("run script"));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  connect(runButton, SIGNAL(clicked()), this, SLOT(runScript()));

  
  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->addWidget(xmlRadio, 0, 0, 1, 1);
  mainLayout->addWidget(tagRadio, 0, 1, 1, 1);
  
  mainLayout->addWidget(baseLabel, 1, 0, 1, 1);
  mainLayout->addWidget(baseButton, 1, 1, 1, 2);
  mainLayout->addWidget(xmlLabel, 2, 0, 1, 1);
  mainLayout->addWidget(xmlButton, 2, 1, 1, 2);
  mainLayout->addWidget(tagLabel, 3, 0, 1, 1);
  mainLayout->addWidget(tagButton, 3, 1, 1, 2);

 mainLayout->addWidget(rPlanLabel, 4, 0, 1, 1);
  mainLayout->addWidget(rPlanButton, 4, 1, 1, 2);

  
  mainLayout->addWidget(planLabel, 5, 0, 1, 1);
  mainLayout->addWidget(planButton, 5, 1, 1, 2);
  mainLayout->addWidget(outLabel, 6, 0, 1, 1);
  mainLayout->addWidget(outButton, 6, 1, 1, 2);
  mainLayout->addWidget(tolLabel, 7, 0, 1, 1);
  mainLayout->addWidget(tolEdit, 7, 1, 1, 2);
  mainLayout->addWidget(levelLabel, 8, 0, 1, 1);
  mainLayout->addWidget(levelEdit, 8, 1, 1, 2);
  mainLayout->addWidget(foldLabel, 9, 0, 1, 1);
  mainLayout->addWidget(foldEdit, 9, 1, 1, 2);

  
  mainLayout->addWidget(modeGroup, 10, 0, 3, 3);
  mainLayout->addWidget(balanceGroup, 13, 0, 3, 3);
  mainLayout->addWidget(display, 16, 0, 4, 3);
  mainLayout->addWidget(saveButton, 20, 0, 1, 1);
  mainLayout->addWidget(runButton, 20, 1, 1, 1);
  
  setLayout(mainLayout);
 
  setText();
}

void RefDialog::runScript(){

  if(planfile ==""){
     QMessageBox::warning(this, tr("Save Script"),
                          tr("please specify plan file  first"));
     return;
  }
  
 if(outfile ==""){
     QMessageBox::warning(this, tr("Save Script"),
                          tr("please specify output grid file first"));
     return;
  }
  
 if(xmlOpt && xmlfile ==""){
    QMessageBox::warning(this, tr("Save Script"),
                          tr("xml option: please specify xml file first"));
    return;
 }

 if(!xmlOpt && tagfile ==""){
   QMessageBox::warning(this, tr("Save Script"),
                          tr("tag option: please specify tag file first"));
   return;
 }

 
  QString command3 = display->toPlainText();
  
  int ret =  system(command3.toStdString().c_str());
  

  
  if(!WIFEXITED(ret))
    {
      if(WIFSIGNALED(ret))
        {
          QMessageBox::information(window(), "mainwindow",
                                   command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
          return;
        }
    }
        

}
void RefDialog::getfile(int i){

  switch(i){
  case 0:
    {
      QString filename =basefile.section('.', 0, 0)+".xml";
      QString format = "xml file(*.xml)";
      xmlfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                              filename,
                                              format);
      QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(xmlfile==""){
        button->setText("Please specify filename");
      }else{
        button->setText(xmlfile);
      }
      
      break;
    }
 case 1:
    {
      QString filename =basefile.section('.', 0, 0)+".tag";
      QString format = "tag file(*.tag)";
      tagfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                             filename,
                                             format);
      QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(tagfile==""){
        button->setText("Please specify filename");
      }else{
        button->setText(tagfile);
      }
      
      break;
    }
  case 3:
    {
      QString filename =basefile.section('.', 0, 0)+".plan";
      QString format = "plan file(*.plan)";
      planfile = QFileDialog::getSaveFileName(this, tr("Save File"),
                                             filename,
                                              format);
     
      QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(planfile==""){
        button->setText("Please specify filename");
      }else{
        button->setText(planfile);
      }
      
      break;
    }

 case 2:
    {
      QString filename =basefile.section('.', 0, 0)+".plan";
      QString format = "plan file(*.plan)";
      rplanfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                               filename,
                                               format);
      
      QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(rplanfile==""){
        button->setText("Please specify filename");
      }else{
        button->setText(rplanfile);
      }
      
      break;
    }
    

    
  case 4:
    {
      int first= basefile.lastIndexOf('/');
      QString directory = basefile.left(first);
      
      QString format = "volume grid file(*.vog)";
      outfile = QFileDialog::getSaveFileName(this, tr("Save File"),
                                             directory,
                                             format);
      QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(outfile==""){
        button->setText("Please specify filename");
      }else{
        button->setText(outfile);
      }
      
      break;
    }
  default:
    break;
  }
  setText();
}


void RefDialog::balanceChanged(int i){
  balanceOpt = i;
  setText();
  
  
}
void RefDialog::optChanged(int i){
  xmlOpt=(i==0);
  setText();
  
}
void RefDialog::modeChanged(int i){
  mode = i;
  setText();
  
}
void  RefDialog::setTol(double d){
  tol = d;
  setText();
  
}
void RefDialog::setFold(double d){
  fold = d;
  setText();
}
void RefDialog::setLevels(int d){
  levels = d;
  setText();
}

void RefDialog::setText(){
  QString marker ="marker -g "+basefile; 
  if(xmlOpt)marker += " -xml " + xmlfile;
  else marker += " -tag " + tagfile;
  if(!xmlOpt && rplanfile !="")marker += " -r " + rplanfile;

  
  marker += QString(" -tol %1").arg(tol);

  if(!xmlOpt)marker += QString(" -levels %1").arg(levels);
  marker += QString(" -fold %1").arg(fold);
  marker += QString(" -mode %1").arg(mode);
  marker += QString(" -balance %1").arg(balanceOpt);
  marker += " -o " + planfile +"\n";


  
  QString refmesh = "refmesh -g " + basefile + " -r " + planfile + " -o " + outfile;
  
  display->setPlainText(marker + refmesh);
  
}

void RefDialog::save(){
  if(planfile ==""){
     QMessageBox::warning(this, tr("Save Script"),
                          tr("please specify plan file  first"));
     return;
  }
  
 if(outfile ==""){
     QMessageBox::warning(this, tr("Save Script"),
                          tr("please specify output grid file first"));
     return;
  }
  
 if(xmlOpt && xmlfile ==""){
    QMessageBox::warning(this, tr("Save Script"),
                          tr("xml option: please specify xml file first"));
    return;
 }

 if(!xmlOpt && tagfile ==""){
   QMessageBox::warning(this, tr("Save Script"),
                          tr("tag option: please specify tag file first"));
   return;
 }
 
  
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save script File"),
                                                  "",
                                                  tr("all Files (*)"));
 
  QFile file(fileName);
  if (!file.open(QFile::WriteOnly | QFile::Text)) {
     QMessageBox::warning(this, tr("Application"),
                          tr("Cannot write file %1:\n%2.")
                          .arg(fileName)
                          .arg(file.errorString()));
     return ;
  }
 
  QTextStream out(&file);

  out << display->toPlainText();
  
  return ;
}

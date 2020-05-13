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
#include <QStackedWidget>
#include "refdialog.h"
#include "progressdialog.h"
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
  setWindowTitle("Refine Grid");
  setAttribute(Qt::WA_DeleteOnClose, true);
  QString tmpname = basefile.left(basefile.lastIndexOf('.'));
  cycle = tmpname.section('_', -1, -1).toInt();
  cycle++;
  if(cycle > 1) basefile = tmpname.section('_', 0, -2)+".vog";
  
  if(bfile=="")return;
  fileButtonGroup = new QButtonGroup(this);
  fileButtonGroup->setExclusive(false);
  
  QButtonGroup* optButtonGroup = new QButtonGroup(this);
  QRadioButton* xmlRadio = new QRadioButton(tr("use -xml option"));
  QRadioButton* tagRadio = new QRadioButton(tr("use -tag option"));

  QRadioButton* tolRadio = new QRadioButton(tr("use -par option"));
  

  optButtonGroup->addButton(xmlRadio, 0);
  optButtonGroup->addButton(tagRadio, 1);
  optButtonGroup->addButton(tolRadio, 2);
  xmlRadio->setChecked(true);
  inputOpt = 0;
  connect(optButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(optChanged(int)));
  

  
  QLabel* baseLabel = new QLabel(tr("baseline grid File:"));
  baseButton = new QPushButton(basefile);
  if(inputOpt==0 && cycle >1)baseButton->setText( basefile.left(basefile.lastIndexOf('.'))+QString("_%1.vog").arg(cycle-1));
  
  pagesWidget = new QStackedWidget;
  {
    QGroupBox *xmlGroup = new QGroupBox;
    xmlGroup->setFlat(true);
    QLabel* xmlLabel = new QLabel(tr("xml file:"));
    xmlButton = new QPushButton(tr("please specify filename"));
    {
      if(cycle > 1) xmlfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.xml").arg(cycle-1);
      else xmlfile = basefile.left(basefile.lastIndexOf('.'))+".xml";
      
      QFile file(xmlfile);
      if(file.exists())xmlButton->setText(xmlfile);
      else xmlfile="";
    }
    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(xmlLabel);
    hlayout->addWidget(xmlButton);
    xmlGroup->setLayout(hlayout);
    pagesWidget->addWidget(xmlGroup);
  }
  
  {
    QGroupBox *tagGroup = new QGroupBox;
    tagGroup->setFlat(true);
    QLabel* tagLabel = new QLabel(tr("tag file:"));
    tagButton = new QPushButton(tr("please specify filename"));
    {
      if(cycle > 1) tagfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.tag").arg(cycle-1);
      else tagfile = basefile.left(basefile.lastIndexOf('.'))+QString(".tag");
      QFile file(tagfile);
      if(file.exists())tagButton->setText(tagfile);
      else tagfile="";
    }
    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(tagLabel);
    hlayout->addWidget(tagButton);
    tagGroup->setLayout(hlayout);
    pagesWidget->addWidget(tagGroup);
  }
  {
    QGroupBox *parGroup = new QGroupBox;
    parGroup->setFlat(true);
    QLabel* parLabel = new QLabel(tr("parameter file:"));
    parButton = new QPushButton(tr("please specify filename"));
    {
      if(cycle > 1) parfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.par").arg(cycle-1);
      else parfile = basefile.left(basefile.lastIndexOf('.'))+QString(".par");
      QFile file(parfile);
      if(file.exists())parButton->setText(parfile);
      else parfile="";
    }
    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(parLabel);
    hlayout->addWidget(parButton);
    parGroup->setLayout(hlayout);
    pagesWidget->addWidget(parGroup);
  }
  
  
  QLabel* rPlanLabel = new QLabel(tr("plan file from last cycle:"));
  rPlanButton = new QPushButton(tr("please specify filename"));
  {
    if(cycle >1){
      rplanfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.plan").arg(cycle-1);
    QFile file(rplanfile);
    if(file.exists())rPlanButton->setText(rplanfile);
    else rplanfile="";
    }else{
      rplanfile="";
      rPlanButton->setText(rplanfile);
      rPlanButton->setEnabled(false);
    }
  }

  QLabel* tolLabel = new QLabel(tr("tolerance:"));
  DoubleEdit* tolEdit = new DoubleEdit;
  tol = 1e-5;
  tolEdit->setValue(tol);
  tolEdit->setBottom(0);
  connect(tolEdit, SIGNAL(valueChanged(double)), this, SLOT(setTol(double)));
  
  
  QLabel* levelLabel = new QLabel(tr("levels:"));
   levelEdit = new IntEdit;
  levels = 1;
  levelEdit->setValue(levels);
  levelEdit->setBottom(1);
  connect(levelEdit, SIGNAL(valueChanged(int)), this, SLOT(setLevels(int)));
  
  QLabel* foldLabel = new QLabel(tr("maximum face folding allowed:"));
  fold = 1.5708;
  DoubleEdit* foldEdit = new DoubleEdit(fold);
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
  outfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.vog").arg(cycle);
  outButton = new QPushButton(outfile); 
  
  QLabel* planLabel = new QLabel(tr("plan File:"));
  planfile = basefile.left(basefile.lastIndexOf('.'))+QString("_%1.plan").arg(cycle);
  planButton = new QPushButton(planfile);
  
  tagButton->setEnabled(false);
  rPlanButton->setEnabled(false);
  levelEdit->setEnabled(false);   
  
  fileButtonGroup->addButton(xmlButton, 0);
  fileButtonGroup->addButton(tagButton, 1);
  fileButtonGroup->addButton(rPlanButton, 2);
  fileButtonGroup->addButton(planButton, 3);
  fileButtonGroup->addButton(outButton, 4);
  fileButtonGroup->addButton(parButton, 5);
  connect(fileButtonGroup, SIGNAL(buttonClicked(int)), this, SLOT(getfile(int)));

  display = new QTextEdit("./marker             \n./refmesh             ");
 
  QPushButton* saveButton = new QPushButton(tr("save script"));
  QPushButton* runButton = new QPushButton(tr("run script"));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  connect(runButton, SIGNAL(clicked()), this, SLOT(runScript()));

  
  QGridLayout *mainLayout = new QGridLayout;
  mainLayout->addWidget(xmlRadio, 0, 0, 1, 1);
  mainLayout->addWidget(tagRadio, 0, 1, 1, 1);
  mainLayout->addWidget(tolRadio, 0, 2, 1, 1);
  mainLayout->addWidget(baseLabel, 1, 0, 1, 1);
  mainLayout->addWidget(baseButton, 1, 1, 1, 2);
  // mainLayout->addWidget(xmlLabel, 2, 0, 1, 1);
  //mainLayout->addWidget(xmlButton, 2, 1, 1, 2);
  // mainLayout->addWidget(tagLabel, 3, 0, 1, 1);
  //mainLayout->addWidget(tagButton, 3, 1, 1, 2);
  
  mainLayout->addWidget(pagesWidget, 2, 0, 1, 3);
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
     QMessageBox::warning(this, tr("Run Script"),
                          tr("please specify plan file  first"));
     return;
  }
  
 if(outfile ==""){
     QMessageBox::warning(this, tr("Run Script"),
                          tr("please specify output grid file first"));
     return;
  }
  
 if(inputOpt==0 && xmlfile ==""){
    QMessageBox::warning(this, tr("Run Script"),
                          tr("xml option: please specify xml file first"));
    return;
 }

 if(inputOpt != 0 && tagfile ==""){
   QMessageBox::warning(this, tr("Run Script"),
                        tr("tag option: please specify tag file first"));
   return;
 }

 
 QString command3 = display->toPlainText();
 
 



 
 ProgressDialog* dialog = new ProgressDialog(command3,false);
 dialog->show();

 
  
  
  
//   int ret =  system(command3.toStdString().c_str());
  

  
//   if(!WIFEXITED(ret))
//     {
//       if(WIFSIGNALED(ret))
//         {
//           QMessageBox::information(window(), "mainwindow",
//                                    command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) );
//           return;
//         }
//     }
  
  
  

}
void RefDialog::getfile(int i){

  switch(i){
  case 0:
    {
      QString filename =basefile.section('.', 0, -2)+".xml";
      QString format = "xml file(*.xml)";
      QString tmpXmlfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        filename,
                                                        format);
      if(tmpXmlfile.isEmpty())return;
      xmlfile = tmpXmlfile;
      // QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(xmlfile==""){
        xmlButton->setText("Please specify filename");
      }else{
        xmlButton->setText(xmlfile);
      }
      
      break;
    }
 case 1:
    {
      QString filename =basefile.section('.', 0, -2)+".tag";
      QString format = "tag file(*.tag)";
      QString tmpTagfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        filename,
                                                        format);
      if(tmpTagfile.isEmpty())return;
      tagfile = tmpTagfile;
      //QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(tagfile==""){
        tagButton->setText("Please specify filename");
      }else{
        tagButton->setText(tagfile);
      }
      
      break;
    }
  case 3:
    {
      QString filename =basefile.section('.', 0, -2)+".plan";
      QString format = "plan file(*.plan)";
      QString tmpPlanfile = QFileDialog::getSaveFileName(this, tr("Save File"),
                                              filename,
                                              format);
      
      if(tmpPlanfile.isEmpty())return;
      planfile = tmpPlanfile;
      if(planfile.section('.', -1, -1)!="plan")planfile+=".plan";
      
      //QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(planfile==""){
        planButton->setText("Please specify filename");
      }else{
        planButton->setText(planfile);
      }
      
      break;
    }

 case 2:
    {
      QString filename =basefile.section('.', 0, -2)+".plan";
      QString format = "plan file(*.plan)";
      QString tmpRplanfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                               filename,
                                               format);

      if(tmpRplanfile.isEmpty())return;
      rplanfile = tmpRplanfile;
      //QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(rplanfile==""){
        rPlanButton->setText("Please specify filename");
      }else{
        rPlanButton->setText(rplanfile);
      }
      
      break;
    }
    

    
  case 4:
    {
      int first= basefile.lastIndexOf('/');
      QString directory = basefile.left(first);
      
      QString format = "volume grid file(*.vog)";
      QString tmpOutfile = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                     directory,
                                                     format);
      if(tmpOutfile.isEmpty())return;
      outfile = tmpOutfile;
      if(outfile.section('.', -1, -1)!="vog")outfile+=".vog";
      
      //QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(outfile==""){
        outButton->setText("Please specify filename");
      }else{
        outButton->setText(outfile);
      }
      
      break;
    }
  case 5:
    {
      QString filename =basefile.section('.', 0, -2)+".par";
      QString format = "parameter file(*.par)";
      QString tmpParfile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        filename,
                                                        format);
      if(tmpParfile.isEmpty())return;
      parfile = tmpParfile;
      //QPushButton* button = qobject_cast<QPushButton*>(fileButtonGroup->button(i));
      
      if(parfile==""){
        parButton->setText("Please specify filename");
      }else{
        parButton->setText(parfile);
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
  inputOpt = i;
  pagesWidget->setCurrentIndex(i);
  // QPushButton* tagButton = qobject_cast<QPushButton*>(fileButtonGroup->button(1));
  //QPushButton* xmlButton = qobject_cast<QPushButton*>(fileButtonGroup->button(0));
  //QPushButton* rplanButton = qobject_cast<QPushButton*>(fileButtonGroup->button(2));

  if(i==0){
    if(cycle > 1)baseButton->setText(basefile.left(basefile.lastIndexOf('.'))+QString("_%1.vog").arg(cycle-1));
    // tagButton->setEnabled(false);
    rPlanButton->setEnabled(false);
    //xmlButton->setEnabled(true);
    levelEdit->setEnabled(false);
  }else if(i==1){
    baseButton->setText(basefile);
    //tagButton->setEnabled(true);
    if(cycle>1)rPlanButton->setEnabled(true);
    //xmlButton->setEnabled(false);
    levelEdit->setEnabled(true);
  }else if(i==2){
    baseButton->setText(basefile);
    rPlanButton->setEnabled(false); 
    //parButton->setEnabled(true);
    //if(cycle>1)rplanButton->setEnabled(true);
    //xmlButton->setEnabled(false);
    levelEdit->setEnabled(false);
  }
    
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
  QString marker ="./marker -g "+baseButton->text(); 
  if(inputOpt==0)marker += " -xml " + xmlfile;
  else if(inputOpt==1)marker += " -tag " + tagfile;
  else if(inputOpt==2)marker += " -par " + parfile;
  if(inputOpt==0 && rplanfile !="")marker += " -r " + rplanfile;

  
  if(inputOpt != 2)marker += QString(" -tol %1").arg(tol);

  if(inputOpt == 1)marker += QString(" -levels %1").arg(levels);
  marker += QString(" -fold %1").arg(fold);
  marker += QString(" -mode %1").arg(mode);
  marker += QString(" -balance %1").arg(balanceOpt);
  marker += " -o " + planfile +"\n";


  
  QString refmesh = "./refmesh -g " + baseButton->text() + " -r " + planfile + " -o " + outfile;
  
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
  
 if(inputOpt==0 && xmlfile ==""){
    QMessageBox::warning(this, tr("Save Script"),
                          tr("xml option: please specify xml file first"));
    return;
 }

 if(inputOpt!=0 && tagfile ==""){
   QMessageBox::warning(this, tr("Save Script"),
                          tr("tag option: please specify tag file first"));
   return;
 }
 
 
 QString casename = basefile.left(basefile.lastIndexOf('.')).section('/', -1, -1);
 QString fileName = basefile.left(basefile.lastIndexOf('/')+1)+"ref_"+casename+QString("_%1").arg(cycle);

 
 
 QString tmpFileName= QFileDialog::getSaveFileName(this, tr("Save script File"),
                                        fileName,
                                        tr("all Files (*)"));
 if(tmpFileName.isEmpty())return;
 fileName = tmpFileName;
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

#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <stdlib.h>
#include <unistd.h>
#include "importwindow.h"
#include "qualitydialog.h"

XdrOption::XdrOption() : QGroupBox(tr("other options")){}

XdrOption::XdrOption( QDomElement& myelem, QWidget *parent ) : QGroupBox(tr("convert to xdr options") , parent){
  QDomElement elem = myelem;
  signalMapper = new QSignalMapper(this);
  QGridLayout* mylayout = new QGridLayout;    
  
  int count = 0;
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) {
    
    QCheckBox* checkBox = new QCheckBox(elem.text());
    tag<<"-"+elem.tagName();
    objs<<checkBox;
    optionList <<"";
    
    mylayout->addWidget(checkBox, count, 0);
    
    if(elem.attribute("value")!="bool"){
      OpGroup* optionGroup = new OpGroup(elem, myelem);
      connect(optionGroup, SIGNAL(textChanged(const QString&)),signalMapper, SLOT(map()));
      signalMapper->setMapping(optionGroup, count);
      mylayout->addWidget(optionGroup, count, 1);
      options << optionGroup;
    }else{
      options << 0;
    }
    
  }
  connect(signalMapper, SIGNAL(mapped(int)),
          this, SLOT(update(int)));  
  
  setLayout(mylayout);      

}

void XdrOption::update(int i){
  optionList[i] = options[i]->currentText();
}
QString XdrOption::currentText(){
  if(objs.size()==0) return "";
  QString text;
  
  for(int i =0; i < objs.size(); i++){
    if(objs[i]->isChecked()){
      text += " " +tag[i] + " " +optionList[i] ;
      //  QStringList tmpList = optionList[i].split("=", QString::SkipEmptyParts);
      // if(tmpList.size()>1) text += " " +tmpList[1];
     
    }
  }
  return text;
}
      

VogOption::VogOption(  QDomElement& myelem, QWidget *parent ) : QGroupBox(tr("convert to vog options") , parent){
  QDomElement elem = myelem.firstChildElement();
  QVBoxLayout* mylayout = new QVBoxLayout;    
  radioGroup = 0;
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr("import.xml"),
                         tr("can not find any element in the children of vogOptions")
                         );
    return;
  }

  for (; !elem.isNull(); elem = elem.nextSiblingElement()) {
    if(elem.attribute("value")=="bool"){
      QCheckBox* checkBox = new QCheckBox(elem.text());
      tag<<elem.tagName();
      objs<<checkBox;
      mylayout->addWidget(checkBox);

    }else if(elem.tagName()=="unit"){
      QGroupBox* unitGroup = new QGroupBox(tr("unit option"));
      alist = elem.text().split(",");
      if(radioGroup){
        delete radioGroup;
        radioGroup = 0;
      }
      radioGroup = new QButtonGroup(this);
      radioGroup->setExclusive(true);
      QHBoxLayout* unitlayout = new QHBoxLayout;
      for(int i=0; i<alist.size(); i++){
        
        QCheckBox* checkBox  = new QCheckBox(alist[i]);
        if(i==0)checkBox->setChecked(true);
        unitlayout->addWidget(checkBox);
        radioGroup->addButton(checkBox);
        radioGroup->setId(checkBox, i);
      }
      connect(radioGroup, SIGNAL(buttonClicked(int)), this, SLOT(unitButtonClicked(int)));
      unitGroup->setLayout(unitlayout);
      mylayout->addWidget(unitGroup);
    }
  }
  QDomElement Lref_elem = myelem.firstChildElement("options");
  Lref_elem = Lref_elem.firstChildElement("Lref");
  LrefBox = new OpGroup(Lref_elem, myelem, this);
  LrefBox->hide();
  mylayout->addWidget(LrefBox); 
  setLayout(mylayout);      
}



void VogOption::unitButtonClicked(int i){
  if(i==(alist.size()-1))LrefBox->show();
  else LrefBox->hide();
}


QString VogOption::currentText(){
  QString text;
  for(int i = 0; i < objs.size(); i++){
    if(objs[i]->isChecked()) text += " -"+tag[i];
  }
  text += " -"+alist[radioGroup->checkedId()];
  if(radioGroup->checkedId()==(alist.size()-1)){
    text += " " +LrefBox->currentText();
  }
  return text;
}
  













ImportWindow::ImportWindow(QDomElement& theelem, QDomElement& theroot, QWidget* parent)
  :GeneralWindow(theelem, theroot, parent)
{
  currentRow = 0;
  importFileName="";

  typesWidget = new QListWidget;
  pagesWidget = new QStackedWidget;
  QGroupBox* typeGroup = new QGroupBox(myelem.attribute("label"));

  QPushButton* usuageButton = new QPushButton(tr("Usuage"));
  QHBoxLayout* aLayout = new QHBoxLayout;
  aLayout->addWidget(typesWidget);
  aLayout->addWidget(usuageButton);
  connect(usuageButton, SIGNAL(clicked()), this, SLOT(usuageButtonClicked()));
  typeGroup->setLayout(aLayout);
  
 
  
  QDomElement elt = myelem.firstChildElement("gridtypes");
  
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'gridtypes' in the children of 'import'")
                         );
    return;
  }
  
  QDomElement elem = elt.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         myelem.tagName()+ tr(" has no child")
                         );
    return;
  }

  QDomElement opt_elem = myelem.firstChildElement("vogOptions");
   if(opt_elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr(" no child vogOptions")
                         );
    return;
   }
   
  option = new VogOption(opt_elem);
 
  for (; !elem.isNull(); elem = elem.nextSiblingElement()) {   
    gridTypes << elem.tagName();
    toXdr << elem.attribute("toXdr");
    toVog << elem.attribute("toVog");
    QListWidgetItem *bdCndiButton = new QListWidgetItem(typesWidget);
    bdCndiButton->setText(elem.attribute("label"));
    bdCndiButton->setTextAlignment(Qt::AlignHCenter);
    bdCndiButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    bdCndiButton->setToolTip(elem.attribute("whatsThis"));
    bdCndiButton->setStatusTip(elem.attribute("whatsThis"));
    bdCndiButton->setWhatsThis(elem.attribute("whatsThis"));
    if(elem.hasAttribute("nameFilter")){
      QGroupBox *aBox = new QGroupBox(tr("Please specify the file name"));
      QVBoxLayout* aLayout = new QVBoxLayout;
      getFileWindow = new GetFileWindow(elem.attribute("nameFilter"), importFileName);
      QDomElement dir_elt = myelem.firstChildElement("directory");
      if(!dir_elt.isNull() && dir_elt.hasAttribute("dir")){
        QString dir = dir_elt.attribute("dir");
        getFileWindow->addDirectory(dir);
      }
      connect(getFileWindow, SIGNAL(fileNameSelected(QString)), this, SLOT(updateFileName(QString)));
      aLayout->addWidget(getFileWindow);
     
      if(!elem.firstChildElement().isNull())
        {
          QDomElement tmpNode = elem.firstChildElement();
        XdrOption*  otherOption = new XdrOption(tmpNode);
        otherOptions << otherOption;
        otherOption->show();
         aLayout->addWidget(otherOption);
        }
      else{
        XdrOption* otherOption  = new XdrOption();
        otherOptions << otherOption;
        otherOption->hide();
        aLayout->addWidget(otherOption);
      }
      
     
      aBox->setLayout(aLayout);
      pagesWidget->addWidget(aBox);
     
    }
    
  }

  connect(typesWidget,
          SIGNAL(currentItemChanged(QListWidgetItem *, QListWidgetItem *)),
          this, SLOT(changePage(QListWidgetItem *, QListWidgetItem*)));
     

  convertButton = new QPushButton(tr("&Convert to Vog"));
  QPushButton *checkButton = new QPushButton(tr("&Check grid quality"));
 

  //  QPushButton *mergeButton = new QPushButton(tr("Merge")); 
  //  QPushButton *extractButton = new QPushButton(tr("Extract"));
  typesWidget->setCurrentRow(0);
  

  
  //connect(mergeButton, SIGNAL(clicked(bool)), this, SLOT(merge(bool)));
  //connect(extractButton, SIGNAL(clicked()), this, SLOT(vog2surface()));
 
 
  connect(convertButton, SIGNAL(clicked()), this, SLOT(convert()));
  connect(checkButton, SIGNAL(clicked()), this, SLOT(check()));

  
  QVBoxLayout *horizontalLayout = new QVBoxLayout;
   horizontalLayout->addWidget(typeGroup);
  horizontalLayout->addWidget(pagesWidget);
     
  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  buttonsLayout->addStretch(1);
 
  buttonsLayout->addWidget(convertButton);
  buttonsLayout->addWidget(checkButton);
  // buttonsLayout->addWidget(mergeButton);
  //  buttonsLayout->addWidget(extractButton);
 
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(horizontalLayout);
  mainLayout->addWidget(option);
  mainLayout->addStretch(1);
  mainLayout->addLayout(buttonsLayout);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  
}
void ImportWindow::updateFileName(QString s){
 
  importFileName = s.section('.', 0, 0);
      
}
void ImportWindow::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
  if (!current)
    current = previous;
  currentRow = typesWidget->row(current);
  if(currentRow == 0){
    convertButton->hide();
    option->hide();
  }else{
    convertButton->show();
    option->show();
  }
  pagesWidget->setCurrentIndex(currentRow);
}

void ImportWindow::convert(){
  
  if(currentRow==0)return;
  QString script_filename = "./output/convert_"+importFileName.section('/', -1);
  QString out_filename="./output/convert_"+importFileName.section('/', -1)+".out";




  QFile outfile(script_filename);
  if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
    QMessageBox::warning(this, tr("savee var file "),
                         tr("Cannot write file %1:\n%2.")
                         .arg(script_filename)
                         .arg(outfile.errorString()));
     return;
  }
  
  QTextStream out(&outfile);
  
  
  QString command2 = toVog[currentRow] +option->currentText() + " " +importFileName ;;
 
  // QMessageBox msgBox;
  //msgBox.setSizeGripEnabled(true);
  //msgBox.setText(tr("converting ")+ importFileName+"                                                       ");

  out <<"#!/bin/bash"<<endl;
  out <<"exec 6>&1"<<endl;
  out <<"exec 7>&2"<<endl;
  out<< "exec &> "<< out_filename <<endl;

  if(toXdr[currentRow]!=""){
    QString command1 = toXdr[currentRow]+otherOptions[currentRow]->currentText()+" " + importFileName;
    

    out << command1 << endl;
    out << command2 << endl;
    // msgBox.setInformativeText(command1 + "\n" + command2);
     emit updateStatus(command1);
     emit updateStatus(command2);
    
  }else{
    out << command2 << endl;
    // msgBox.setInformativeText(command2);
    emit updateStatus(command2);
  }
  out<<"exec 1>&6 6>&- " << endl;
  out<<"exec 2>&7 7>&- " << endl;

  outfile.close();
  QString command3 = "chmod 777 " + script_filename;

  int ret = system(command3.toStdString().c_str());

  if(!WIFEXITED(ret))
        {
          if(WIFSIGNALED(ret))
            {
              QMessageBox::information(window(), "system",
                                       command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
              return;
            }
        }
  
  ret = system(script_filename.toStdString().c_str());
  if(!WIFEXITED(ret)){
    if(WIFSIGNALED(ret)){
      QMessageBox::information(window(), "system",
                               script_filename + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
      return;
    }
 }
 
  QFile file(out_filename);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
     QMessageBox::information(window(), "file io",
                                 tr("Cannot open ") + out_filename + tr(" for reading!"));
     return;
  }
  
  QTextStream in(&file);
  QApplication::setOverrideCursor(Qt::WaitCursor);
  // msgBox.setDetailedText(in.readAll());
  emit updateStatus(in.readAll());
  QApplication::restoreOverrideCursor();
  file.close();

   //update DOM tree
  QDomElement buttons_elem = myroot.firstChildElement(tr("buttons"));
  QDomElement convert_elem = buttons_elem.firstChildElement(tr("convert"));
  convert_elem.setAttribute(tr("status"), tr("completed"));
  convert_elem.setAttribute(tr("dir"),  importFileName.section('/', 0, -2));
  convert_elem.setAttribute(tr("filename"), importFileName.section('/', -1));
  convert_elem.setAttribute(tr("cmd"), command2);
  //  msgBox.exec();
  
}
// void ImportWindow::vog2surface(){
//   if(importFileName==""){
//     QMessageBox::warning(window(), tr(".xml"),
//                          tr("Please specify a filename first")
//                          );
//     return;
//   }
//   QFile exist_test(importFileName+tr(".vog"));
//   if(!(exist_test.exists())){
//     QMessageBox::warning(window(), tr(".xml"),
//                          tr("Please convert the file to volume grid format first")
//                          );
//     return;
//   }





//   QString script_filename = "./output/vog2surface_"+importFileName.section('/', -1);
//   QString out_filename="./output/vog2surface_"+importFileName.section('/', -1)+".out";
  
//   QFile outfile(script_filename);
//   if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
//     return;
//   }
  
//   QTextStream out(&outfile);
  
  
//   QString command2 = "vog2surface " +importFileName;
//   emit updateStatus(command2);
//   out <<"#!/bin/bash"<<endl;
//   out <<"exec 6>&1"<<endl;
//   out <<"exec 7>&2"<<endl;
//   out<< "exec &> "<< out_filename <<endl;
//   out << command2 << endl;
//   out<<"exec 1>&6 6>&- " << endl;
//   out<<"exec 2>&7 7>&- " << endl;

//   outfile.close();
//   QString command3 = "chmod 777 " + script_filename;

//   system(command3.toStdString().c_str());
//   system(script_filename.toStdString().c_str());
//   QFile file(out_filename);
//   if (!file.open(QFile::ReadOnly | QFile::Text)) {
//     return;
//   }
  
//   QTextStream in(&file);
//   QApplication::setOverrideCursor(Qt::WaitCursor);
//   // msgBox.setDetailedText(in.readAll());
//   emit updateStatus(in.readAll());
//   QApplication::restoreOverrideCursor();
//   file.close();

//   //update DOM tree
//   QDomElement buttons_elem = myroot.firstChildElement(tr("buttons"));
//   QDomElement extract_elem = buttons_elem.firstChildElement(tr("extract"));
//   extract_elem.setAttribute(tr("status"), tr("completed"));
//   extract_elem.setAttribute(tr("dir"),  importFileName.section('/', 0, -2));
//   extract_elem.setAttribute(tr("filename"), importFileName.section('/', -1));
//   extract_elem.setAttribute(tr("cmd"), command2);
//   //if extract is completed, import is completed
//   myroot.setAttribute(tr("status"), tr("completed"));

//   updateStatus(tr("surface file generated"));
// }






  
void ImportWindow::check(){
  if(importFileName==""){
      QMessageBox::warning(window(), tr("import"),
                           tr("Please specify a filename first")
                           );
      return;
  }
  QFile exist_test(importFileName+tr(".vog"));
  if(!(exist_test.exists())){
    QMessageBox::warning(window(), tr("import"),
                         tr("Please convert the file to volume grid format first")
                         );
    return;
  }
  
  QString script_filename = "./output/check_"+importFileName.section('/', -1);
  QString out_filename="./output/check_"+importFileName.section('/', -1)+".out";




  QFile outfile(script_filename);
  if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
     QMessageBox::warning(this, tr("file io "),
                            tr("Cannot write file %1:\n%2.")
                            .arg(script_filename)
                          .arg(outfile.errorString()));
    return;
  }
  
  QTextStream out(&outfile);
  
  
  QString command2 = "vogcheck " +importFileName.section('/', -1);
  emit updateStatus(command2);
 
  
  out <<"#!/bin/bash"<<endl;
  out <<"exec 6>&1"<<endl;
  out <<"exec 7>&2"<<endl;
  out<< "exec &> "<< out_filename <<endl;
  out<<"cd " << importFileName.section('/',0,  -2)<<endl;
  out << command2 << endl;
  out<<"exec 1>&6 6>&- " << endl;
  out<<"exec 2>&7 7>&- " << endl;

  outfile.close();
  QString command3 = "chmod 777 " + script_filename;

  int ret = system(command3.toStdString().c_str());

  if(!WIFEXITED(ret))
    {
      if(WIFSIGNALED(ret))
        {
          QMessageBox::information(window(), "mainwindow",
                                   command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
          return;
        }
      exit(0);
    }
  ret = system(script_filename.toStdString().c_str());

  if(!WIFEXITED(ret))
    {
      if(WIFSIGNALED(ret))
        {
          QMessageBox::information(window(), "mainwindow",
                                   command3 + tr(" was terminated with the signal %d") + WTERMSIG(ret) ); 
          return;
        }
      exit(0);
    }
  
  QFile file(out_filename);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    QMessageBox::warning(this, tr(" file io "),
                         tr("Cannot read file %1:\n%2.")
                         .arg(out_filename)
                         .arg(file.errorString()));
    return;
  }
  
  QTextStream in(&file);
  QApplication::setOverrideCursor(Qt::WaitCursor);
  // msgBox.setDetailedText(in.readAll());
  updateStatus(in.readAll());
  QApplication::restoreOverrideCursor();
  file.close();

  //update DOM tree
  QDomElement buttons_elem = myroot.firstChildElement(tr("buttons"));
  QDomElement check_elem = buttons_elem.firstChildElement(tr("check"));
  check_elem.setAttribute(tr("status"), tr("completed"));
  check_elem.setAttribute(tr("dir"),  importFileName.section('/', 0, -2));
  check_elem.setAttribute(tr("filename"), importFileName.section('/', -1));

  
  QString qualityFileName = importFileName+tr(".quality");
  //QString  qualityFileName="bad.quality";
  QualityDialog qualityDialog(qualityFileName, this);
  qualityDialog.exec();
  
    

}

// void ImportWindow::merge(bool checked){
//   if(!checked) return;
//   return ;
// }

// void ImportWindow::extract(){
//   QString path = importFileName.section('/', 0, -2);
//   QString fname = importFileName.section('/', -1);
  
//   QString script_filename = "./output/extract_"+importFileName.section('/', -1);
//   QString out_filename="./output/extract_"+importFileName.section('/', -1)+".out";
  



//   QFile outfile(script_filename);
//   if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
//     return;
//   }
  
//   QTextStream out(&outfile);
  
  
//   QString command2 = "extract -en " +fname;
 

//   emit updateStatus(command2);
  
//   out <<"#!/bin/bash"<<endl;
//   out <<"exec 6>&1"<<endl;
//   out <<"exec 7>&2"<<endl;
//   out<< "exec &> "<< out_filename <<endl;
//   if(!path.isEmpty()) out<<"cd "<<path<<endl;
//   out << command2 << endl;
  
//   out<<"exec 1>&6 6>&- " << endl;
//   out<<"exec 2>&7 7>&- " << endl;

//   outfile.close();
//   QString command3 = "chmod 777 " + script_filename;

//   system(command3.toStdString().c_str());
//   system(script_filename.toStdString().c_str());
  

//   QFile file(out_filename);
//   if (!file.open(QFile::ReadOnly | QFile::Text)) {
//     return;
//   }
  
//   QTextStream in(&file);
//   QApplication::setOverrideCursor(Qt::WaitCursor);
//   emit updateStatus(in.readAll());
//   QApplication::restoreOverrideCursor();
//   file.close();

//   //update DOM tree
//   QDomElement buttons_elem = myroot.firstChildElement(tr("buttons"));
//   QDomElement convert_elem = buttons_elem.firstChildElement(tr("extract"));
//   extract_elem.setAttribute(tr("status"), tr("completed"));
//   extract_elem.setAttribute(tr("dir"), path);
//   extract_elem.setAttribute(tr("filename"), fname);
  
// }



void ImportWindow::usuageButtonClicked(){
  
  QDomElement elt = myroot.firstChildElement("usuage");
  
  if(elt.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'usuage' in the children of 'import'")
                          );
    return;
  }
  
  QDomElement elem = elt.firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("'usuage' has no child")
                         );
    return;
  }
   
  for (int i = 0;  i<currentRow; i++) elem = elem.nextSiblingElement();
  QMessageBox msgBox;
  msgBox.setText(elem.text());

  msgBox.exec();
     
     
    
  
}






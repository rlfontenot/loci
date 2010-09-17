#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <stdlib.h>
#include <unistd.h>
#include "importwindow.h"
#include "qualitydialog.h"
#include "progressdialog.h"
#include "helpwindow.h"
XdrOption::XdrOption() : QGroupBox(tr("other options")){}

XdrOption::XdrOption( QDomElement& myelem, QWidget *parent ) : QGroupBox(tr("other options") , parent){
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
      VarGBox* optionGroup = new VarGBox(elem);
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
  QString text = QString(" -") + options[i]->currentText();
  text.replace(QString("="), QString(" "));
  optionList[i] = text;
}
QString XdrOption::currentText(){
  if(objs.size()==0) return "";
  QString text;
  
  for(int i =0; i < objs.size(); i++){
    if(objs[i]->isChecked()){
      text += optionList[i] ;
    }
  }
  return text;
}
      

VogOption::VogOption(  QDomElement& myelem, QWidget *parent ) : QGroupBox(tr("convert to vog options") , parent){
  QDomElement elem = myelem.firstChildElement();
  
  QVBoxLayout* mylayout = new QVBoxLayout;    
  radioGroup = 0;
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr("main.xml"),
                         tr("can not find any element in the children of vogOptions")
                         );
    return;
  }

  for (; !elem.isNull(); elem = elem.nextSiblingElement()) {
   
    if(elem.attribute("element")=="bool"){
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
  LrefBox = new VarGBox(Lref_elem, this);
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
  




ImportWindow::ImportWindow(QDomElement& theelem,  QWidget* parent)
  :GeneralGroup(theelem, parent)
{

  setAttribute(Qt::WA_DeleteOnClose, true);


  QHBoxLayout* helpLayout = new QHBoxLayout;
  QPushButton* helpButton = new QPushButton(tr("help"));
  connect(helpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  helpLayout->addWidget(helpButton);
  helpLayout->setAlignment(helpButton, Qt::AlignRight);
  
  
  currentRow = 0;
  importFileName="";
  
  typesWidget = new QListWidget;
  pagesWidget = new QStackedWidget;
  QGroupBox* typeGroup = new QGroupBox(myelem.attribute(tr("label")));

  QPushButton* usageButton = new QPushButton(tr("Usage"));
  QHBoxLayout* aLayout = new QHBoxLayout;
  aLayout->addWidget(typesWidget);
  aLayout->addWidget(usageButton);
  connect(usageButton, SIGNAL(clicked()), this, SLOT(usageButtonClicked()));
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
    if(elem.hasAttribute("toXdr")) toXdr << elem.attribute("toXdr");
    else   toXdr <<QString();
    
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
      getFileWindow = new FindFileWindow(elem, importFileName);
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
          SIGNAL(currentRowChanged(int)),
          this, SLOT(changePage(int)));
     

  convertButton = new QPushButton(tr("&Convert to Vog"));
  typesWidget->setCurrentRow(0);
  

  
  
 
 
  connect(convertButton, SIGNAL(clicked()), this, SLOT(convert()));
  
  
  QVBoxLayout *horizontalLayout = new QVBoxLayout;
  horizontalLayout->addWidget(typeGroup);
  horizontalLayout->addWidget(pagesWidget);
     
  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  buttonsLayout->addStretch(1);
 
  buttonsLayout->addWidget(convertButton);
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(helpLayout);
  mainLayout->addLayout(horizontalLayout);
  mainLayout->addWidget(option);
  mainLayout->addStretch(1);
  mainLayout->addLayout(buttonsLayout);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  if(elt.firstChildElement().attribute("vogOptions")=="false")option->hide();
}
void ImportWindow::updateFileName(QString s){
  importFileName = s.section('.', 0, -2);
}
void ImportWindow::changePage(int currentR)
{
  currentRow = currentR;
  pagesWidget->setCurrentIndex(currentR);
  QDomElement elt = myelem.firstChildElement("gridtypes");
  
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'gridtypes' in the children of 'import'")
                         );
    return;
  }
  
  QDomElement elem = elt.firstChildElement();
  for (int i = 0; i < currentR; i++) elem = elem.nextSiblingElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         myelem.tagName()+ tr(" has no child")
                         );
    return;
  }
  QString s = elem.attribute("current");
  
    
  importFileName = s.section('.', 0, -2);
 
  if(elem.attribute("vogOptions")=="false")option->hide();
  else option->show();
   
}
void ImportWindow::helpClicked(){
  HelpWindow* helpwindow = new HelpWindow("page_import.html");
  helpwindow->show();
}
void ImportWindow::convert(){
  if(importFileName.isEmpty()&&gridTypes[currentRow]!="cfd"){
    QMessageBox::warning(window(), tr("convert to vog"),
                         tr("Please specify filename first")
                         );
    return;
  } 


  bool vop = true;

  QDomElement elem = myelem.firstChildElement("gridtypes");
  elem = elem.firstChildElement();
  for(int i = 0; i < currentRow; i++)elem = elem.nextSiblingElement();
  if(elem.attribute("vogOptions")=="false") vop = false;
  

  
  QString command2;
  if(toXdr[currentRow]=="")command2= toVog[currentRow]+" "+ otherOptions[currentRow]->currentText();
  else command2= toVog[currentRow];
  

  if(vop) command2 += " "+option->currentText()+ " " +importFileName ;
  else command2 += " " +importFileName ;


  QString command1;
  if(toXdr[currentRow]!=""){
    command1 = toXdr[currentRow]+otherOptions[currentRow]->currentText()+" " + importFileName;
  }
  
  QString command;
  if(toXdr[currentRow]!=""){
    command = command1+'\n'+command2;
  }else{
    command = command2;
  }

  if(gridTypes[currentRow]=="cfd"){
    command1 = toXdr[currentRow]+otherOptions[currentRow]->currentText();
    command2 = toVog[currentRow]+ " "+option->currentText() + " grid.xdr";
    command = command1+'\n'+command2;
    qDebug() << command;
    ProgressDialog* dialog = new ProgressDialog(command, importFileName.section('/', 0, -2));
    dialog->show(); 
  }else{
    qDebug() << command;
    ProgressDialog* dialog = new ProgressDialog(command, QString());
    dialog->show();
  } 







  
  //update DOM tree

  QDomElement myroot = myelem.ownerDocument().documentElement();
  QDomElement buttons_elem = myroot.firstChildElement(tr("buttons"));
  QDomElement convert_elem = buttons_elem.firstChildElement(tr("convert"));
  convert_elem.setAttribute(tr("status"), tr("completed"));
  convert_elem.setAttribute(tr("dir"),  importFileName.section('/', 0, -2));
  convert_elem.setAttribute(tr("filename"), importFileName.section('/', -1));
  convert_elem.setAttribute(tr("cmd"), command2);

  
}






void ImportWindow::usageButtonClicked(){
  
 
  
  QDomElement elem = myelem.firstChildElement("gridtypes").firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("'usage' has no child")
                         );
    return;
  }
  
  for (int i = 0;  i<currentRow; i++) elem = elem.nextSiblingElement();

  QDomElement elt = myelem.firstChildElement("usage").firstChildElement(elem.tagName());
  
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'usage' in the children of 'import'")
                         );
    return;
  }
  
  QMessageBox msgBox;
  msgBox.setMinimumSize(1000, 700);
  msgBox.setText(elt.text());
  

  msgBox.exec();
     
     
    
  
}






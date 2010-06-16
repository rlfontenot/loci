#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <stdlib.h>
#include <unistd.h>
#include "importwindow.h"
#include "qualitydialog.h"
#include "progressdialog.h"
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
  optionList[i] = options[i]->currentText();
}
QString XdrOption::currentText(){
  if(objs.size()==0) return "";
  QString text;
  
  for(int i =0; i < objs.size(); i++){
    if(objs[i]->isChecked()){
      text += " " +tag[i] + " " +optionList[i] ;
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
  
  currentRow = 0;
  importFileName="";

  typesWidget = new QListWidget;
  pagesWidget = new QStackedWidget;
  QGroupBox* typeGroup = new QGroupBox(myelem.attribute("label"));

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
    typesWidget->setCurrentRow(0);
  

  
  
 
 
  connect(convertButton, SIGNAL(clicked()), this, SLOT(convert()));
  
  
  QVBoxLayout *horizontalLayout = new QVBoxLayout;
   horizontalLayout->addWidget(typeGroup);
  horizontalLayout->addWidget(pagesWidget);
     
  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  buttonsLayout->addStretch(1);
 
  buttonsLayout->addWidget(convertButton);
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(horizontalLayout);
  mainLayout->addWidget(option);
  mainLayout->addStretch(1);
  mainLayout->addLayout(buttonsLayout);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  
}
void ImportWindow::updateFileName(QString s){
   importFileName = s.section('.', 0, -2);
}
void ImportWindow::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
  if (!current)
    current = previous;
  currentRow = typesWidget->row(current);
  pagesWidget->setCurrentIndex(currentRow);
}

void ImportWindow::convert(){
  
   QString command2 = toVog[currentRow] +option->currentText() + " " +importFileName ;;
 
  QString command1;
  if(toXdr[currentRow]!=""){
    command1 = toXdr[currentRow]+otherOptions[currentRow]->currentText()+" " + importFileName;
     
      
  }else{
  }
  
  QString command;
  if(toXdr[currentRow]!=""){
    command = command1+'\n'+command2;
  }else{
    command = command2;
  }
  ProgressDialog* dialog = new ProgressDialog(command, QString());
  dialog->show();
 







  
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
  
  QDomElement elt = myelem.firstChildElement("usage");
  
  if(elt.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'usage' in the children of 'import'")
                          );
    return;
  }
  
  QDomElement elem = elt.firstChildElement();
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("'usage' has no child")
                         );
    return;
  }
   
  for (int i = 0;  i<currentRow; i++) elem = elem.nextSiblingElement();
  QMessageBox msgBox;
  msgBox.setText(elem.text());

  msgBox.exec();
     
     
    
  
}






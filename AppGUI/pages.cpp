#include <QtGui>
//#include <QtCore>
#include <QHBoxLayout>
#include <QComboBox>
#include <QGridLayout>
#include <QLabel>
#include <QMessageBox>
#include <QtDebug>
#include <QLineEdit>
#include <QListWidget>
#include <QStackedWidget>
#include <QListWidgetItem>
#include <QString>
#include "pages.h"
bool conditionIsSatisfied(const QDomElement& theroot, const QString& condition){
  
  QDomElement elem= theroot.firstChildElement("mainWindow");
  elem = elem.firstChildElement("physicsPage");
  if(!elem.hasAttribute("state")) return false;
  if(condition.isEmpty()) return false;
  QStringList state = elem.attribute("state").split(";");
  QStringList list1=condition.split("==");
  if(list1.size()!=2)return false;
  for(int i = 0; i < state.size(); i++){
    if(state[i].contains(list1[0]) && state[i].contains(list1[1])) return true;
  }
  return false;
}

FloatEdit::FloatEdit(QWidget *parent) : QLineEdit(parent){
  
   validator = new QDoubleValidator(this);
   setValidator(validator);
   connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
   
 }

FloatEdit::FloatEdit(double d, QWidget*parent):QLineEdit(parent){
  validator = new QDoubleValidator(this);
  setValidator(validator);
  QString str;
  str.setNum(d);
  setText(str);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
  emit valueChanged(d);
}
void FloatEdit::setValue(double d){
  QString str;
  str.setNum(d);
  setText(str);
  emit valueChanged(d);
}
void FloatEdit::setValue(int value){
  double d = validator->bottom() +(value+1000)/2000.0*(validator->top()-validator->bottom());
  QString str;
  str.setNum(d);
  setText(str);
  emit valueChanged(d);
}
void FloatEdit::changeValue(const QString& str){
  double d = str.toDouble();
  emit valueChanged(d);
}
  
double FloatEdit::value(){
  QString str;
  str = text();
  return str.toDouble();
}

void FloatEdit::setBottom(double d){
  validator->setBottom(d);
}

void FloatEdit::setTop(double d){
  validator->setTop(d);
}

void FloatEdit::setRange(double d1, double d2){
  validator->setRange(d1, d2);
}

















IntEdit::IntEdit(QWidget *parent) : QLineEdit(parent){
   validator = new QIntValidator(this);
   setValidator(validator);
   connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
   
 }

IntEdit::IntEdit(int d, QWidget*parent):QLineEdit(parent){
  validator = new QIntValidator(this);
  setValidator(validator);
  QString str;
  str.setNum(d);
  setText(str);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
  emit valueChanged(d);
}
void IntEdit::setValue(int d){
  QString str;
  str.setNum(d);
  setText(str);
  emit valueChanged(d);
}

void IntEdit::changeValue(const QString& str){
  int d = str.toInt();
  emit valueChanged(d);
}
  
int IntEdit::value(){
  QString str;
  str = text();
  return str.toInt();
}

void IntEdit::setBottom(int d){
  validator->setBottom(d);
}

void IntEdit::setTop(int d){
  validator->setTop(d);
}

void IntEdit::setRange(int d1, int d2){
  validator->setRange(d1, d2);
}

GeneralGroup::GeneralGroup( QDomElement& my_elem, QDomElement& theroot, QWidget *parent ) : QGroupBox(my_elem.hasAttribute("title")?my_elem.attribute("title"):my_elem.tagName(), parent),myelem(my_elem),myroot(theroot){
  
  if(myelem.hasAttribute("whatsThis"))setWhatsThis(myelem.attribute("whatsThis"));
  if(myelem.hasAttribute("toolTip"))setToolTip(myelem.attribute("toolTip"));
  if(myelem.hasAttribute("statusTip"))setStatusTip(myelem.attribute("statusTip"));
  
  connect(this, SIGNAL(toggled(bool)), this, SLOT(updateChecked()));
  if(myelem.hasAttribute("checked")){
    bool isChecked = myelem.attribute("checked").toInt()==1;
    setCheckable(true);
    setChecked(isChecked);
  }
  
}
QString GeneralGroup::currentText(){
  return QString();
}

void GeneralGroup::updateChecked(){
  if(isChecked())myelem.setAttribute("checked", 1);
  else myelem.setAttribute("checked", 0);
  emit textChanged(currentText());
  
}


OpGroup::OpGroup(  QDomElement& my_elem, QDomElement& theroot, QWidget *parent ) : GeneralGroup(my_elem, theroot, parent){

  

  if(myelem.attribute("element")=="vector"){

    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->setSpacing(0);
    QButtonGroup *radioGroup = new QButtonGroup(this);
    radioGroup->setExclusive(true);
    QRadioButton* radio0 = new QRadioButton(" Cartesian  ");
    QRadioButton* radio1 = new QRadioButton("   Polar    ");
    QRadioButton* radio2 = new QRadioButton("   Scale    ");
    if(myelem.attribute("format")=="polar"){
      radio1->setChecked(true);
      
    }else if(myelem.attribute("format")=="scale"){
      radio2->setChecked(true);
      
    }else{
      radio0->setChecked(true);
      myelem.setAttribute("format","cartesian"); 
    }
    radioGroup->addButton(radio0);
    radioGroup->addButton(radio1);
    radioGroup->addButton(radio2);
    
    radioGroup->setId(radio0, 0);
    radioGroup->setId(radio1, 1);
    radioGroup->setId(radio2, 2);
 
    connect(radioGroup, SIGNAL(buttonClicked(int)),this, SLOT(updateLabels(int)));   
    buttonLayout->addWidget(radio0);
    buttonLayout->addWidget(radio1);
    buttonLayout->addWidget(radio2);

    
    
    QVBoxLayout* myLayout = new QVBoxLayout;
    myLayout->setSpacing(1);
    
    QHBoxLayout* xLayout = new QHBoxLayout;


   
   QPointer<QLabel> xLabel = new QLabel("x");
   QPointer<QLabel> yLabel = new QLabel("y");
   QPointer<QLabel> zLabel = new QLabel("z");

   if(myelem.attribute("format")=="polar"){
     xLabel->setText("mag");
     yLabel->setText("theta");
     zLabel->setText("phi");
   }
   QPointer<FloatEdit> xEditor = new FloatEdit;
   QPointer<FloatEdit> yEditor = new FloatEdit;
   QPointer<FloatEdit> zEditor = new FloatEdit;
    
   connect(xEditor, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentX(const QString&)));
   connect(yEditor, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentY(const QString&)));
   connect(zEditor, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentZ(const QString&)));



    xLabel->setAlignment(Qt::AlignHCenter);
    yLabel->setAlignment(Qt::AlignHCenter);
    zLabel->setAlignment(Qt::AlignHCenter);
    labels<<xLabel;
    labels<<yLabel;
    labels<<zLabel;
    mfs<<xEditor;
    mfs<<yEditor;
    mfs<<zEditor;
    xLayout->addWidget(xLabel);
    xLayout->addWidget(xEditor);
   
    xLayout->addWidget(yLabel);
    xLayout->addWidget(yEditor);

     xLayout->addWidget(zLabel);
     xLayout->addWidget(zEditor);

     
    myLayout->addLayout(buttonLayout);
    myLayout->addLayout(xLayout);
       
    if(myelem.hasAttribute("currentX"))xEditor->setText(myelem.attribute("currentX"));
    if(myelem.hasAttribute("currentY"))yEditor->setText(myelem.attribute("currentY"));
    if(myelem.hasAttribute("currentZ"))zEditor->setText(myelem.attribute("currentZ"));    
    if(myelem.attribute("format")=="scale"){
    
     yEditor->setValue(0.0);
     zEditor->setValue(0.0);
     yEditor->setEnabled(false);
     zEditor->setEnabled(false);
   }


    if(!myelem.attributeNode("unit").isNull()){
      QStringList alist = myelem.attribute("unit").split(",");
      QComboBox *unitCombo = new QComboBox;
      unitCombo->setEditable(true);
      unitCombo->setMinimumContentsLength(5);
      unitCombo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLength);
      unitCombo->setInsertPolicy(QComboBox::InsertAtBottom);
      unitCombo->addItems(alist);
      if(myelem.hasAttribute("currentUnit")) {
        int index = unitCombo->findText(myelem.attribute("currentUnit"));
        if(index >= 0)unitCombo->setCurrentIndex(index);
        else unitCombo->setEditText(myelem.attribute("currentUnit"));
      }else{
        myelem.setAttribute("currentUnit",unitCombo->currentText());
      }
      connect(unitCombo, SIGNAL(editTextChanged(const QString&)), this, SLOT(updateUnitList(const QString&)));
      connect(unitCombo, SIGNAL(currentIndexChanged(const QString &)),this, SLOT(updateCurrentUnit(const QString&)));
      myLayout->addWidget(unitCombo);
          
    }
    setLayout(myLayout);      
            
  }else if(myelem.attribute("element")=="dvector"){
    // setTitle(myelem.text());
    
    if(!myroot.firstChildElement("components").attribute("components").isEmpty()){
      QStringList components = myroot.firstChildElement("components").attribute("components").split(",");
    
     QVBoxLayout* myLayout = new QVBoxLayout;
     
      for(int i =0; i < components.size(); i++){ 
        QPointer<QLabel> aLabel = new QLabel(components[i]);
        QPointer<FloatEdit> aEditor = new FloatEdit;
        QHBoxLayout* hLayout = new QHBoxLayout;
        hLayout->addWidget(aLabel);
        hLayout->addWidget(aEditor);
        connect(aEditor, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrent(const QString&)));
        aLabel->setAlignment(Qt::AlignHCenter);
        labels<<aLabel;
        mfs<<aEditor;
        qobject_cast<QVBoxLayout*>(myLayout)->addLayout(hLayout);
      }
      setLayout(myLayout);
       
      if(myelem.hasAttribute("current")){
        QStringList valueList = myelem.attribute("current").split(",");
        for(int i = 0; i < valueList.size(); i++){
          if(i < components.size()) qobject_cast<FloatEdit*>(mfs[i])->setValue(valueList[i].toDouble());
        }
      }  
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
        this->hide();
    }
  }else if(myelem.attribute("element")=="selection"){
    setFlat(false);
    QDomElement opt = myelem.firstChildElement();
    if(opt.isNull())return;
    QButtonGroup* buttonGroup = new QButtonGroup(this);
    buttonGroup->setExclusive(false);
    QHBoxLayout* myLayout = new QHBoxLayout;
    // myLayout->setSpacing(1);
    int count = 0;  
    for(; !opt.isNull(); opt=opt.nextSiblingElement(), count++){
      QPointer<QCheckBox> checkBox  = new QCheckBox(opt.tagName());
      buttonGroup->addButton(checkBox);
      buttonGroup->setId(checkBox, count);
       myLayout->addWidget(checkBox);
      if(opt.hasAttribute("checked"))checkBox->setChecked(opt.attribute("checked").toInt()==1);
      mfs<<checkBox;
    }
    connect(buttonGroup, SIGNAL(buttonClicked(int)),this, SLOT(updateSelection(int))); 
    setLayout(myLayout);
  }else{
    QHBoxLayout* myLayout = new QHBoxLayout(this);  
    if(myelem.attribute("element")=="int"){
      QPointer<QSpinBox> spinBox = new QSpinBox;
      if(!myelem.attributeNode("bottom").isNull()){
        spinBox->setMinimum(myelem.attribute("bottom").toInt());
      }
      if(!myelem.attributeNode("top").isNull()){
        spinBox->setMaximum(myelem.attribute("top").toInt());
      }
      if(!myelem.attributeNode("range").isNull()){
        QString range = myelem.attribute("range");
        QStringList rangeList = range.split(",");
        if(rangeList.size()!=2){
          QMessageBox::warning(window(), tr(".xml"),
                               tr("incorrect format for range, should be range='bottom, top'")
                               );
          return;
        }
        int bottom = rangeList[0].toInt();
        int top = rangeList[1].toInt();
        if(bottom>top){
          int tmp = bottom;
          bottom = top;
          top = tmp;
        }
            
        spinBox->setRange(bottom, top);
      }
      if(myelem.hasAttribute("current"))spinBox->setValue(myelem.attribute("current").toInt());
      connect(spinBox, SIGNAL(valueChanged(const QString&)), this, SLOT(updateCurrent(const QString&)));
      myLayout->addWidget(spinBox);
      mfs << spinBox;
    }else if(myelem.attribute("element")=="float"){
      QPointer<FloatEdit> valueEdit = new FloatEdit;
      if(myelem.hasAttribute("bottom"))valueEdit->setBottom(myelem.attribute("bottom").toDouble());
      if(myelem.hasAttribute("top"))valueEdit->setTop(myelem.attribute("top").toDouble());
      if(myelem.hasAttribute("range")){
        QString range = myelem.attribute("range");
        QStringList rangeList = range.split(",");
        if(rangeList.size()!=2){
          QMessageBox::warning(window(), tr(".xml"),
                               tr("incorrect format for range, should be range='bottom, top'")
                               );
          return;
        }
        double bottom = rangeList[0].toDouble();
        double top = rangeList[1].toDouble();
        if(bottom>top){
          double tmp = bottom;
          bottom = top;
          top = tmp;
        }
            
        valueEdit->setRange(bottom, top);
      }
      if(myelem.hasAttribute("current"))valueEdit->setValue(myelem.attribute("current").toDouble());
      myLayout->addWidget(valueEdit);
      connect(valueEdit, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrent(const QString&)));
      mfs<<valueEdit;
    }else if(myelem.attribute("element")=="string"){
      QPointer<QLineEdit> valueEdit = new QLineEdit;
      if(myelem.hasAttribute("current"))valueEdit->setText(myelem.attribute("current"));
      myLayout->addWidget(valueEdit);
      connect(valueEdit, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrent(const QString&)));
      mfs<<valueEdit;
    }else if(myelem.attribute("element")=="novalue"){


      //setTitle("");
       setFlat(true);
      
    }else{
       QMessageBox::warning(window(), tr(".xml"),
                            tr("options: don't know how to handle it yet: ") + myelem.attribute("element")
                            );
       return;
    }
    
    if(!myelem.attributeNode("unit").isNull()){
      QStringList alist = myelem.attribute("unit").split(",");
      QComboBox *unitCombo = new QComboBox;
      unitCombo->setEditable(true);
      unitCombo->setMinimumContentsLength(5);
      unitCombo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLength);
      unitCombo->setInsertPolicy(QComboBox::InsertAtBottom);
      unitCombo->addItems(alist);
      if(myelem.hasAttribute("currentUnit")) {
        int index = unitCombo->findText(myelem.attribute("currentUnit"));
        if(index >= 0)unitCombo->setCurrentIndex(index);
        else unitCombo->setEditText(myelem.attribute("currentUnit"));
      }else{
        myelem.setAttribute("currentUnit",unitCombo->currentText());
      }
      
      connect(unitCombo, SIGNAL(editTextChanged(const QString&)), this, SLOT(updateUnitList(const QString&)));
      connect(unitCombo, SIGNAL(currentIndexChanged(const QString &)),this, SLOT(updateCurrentUnit(const QString&)));
      myLayout->addWidget(unitCombo);
      
    }      
    setLayout(myLayout);

  }//else not vectors
  emit textChanged(currentText());
}

void OpGroup::updateSelection(int id){
  QDomElement opt = myelem.firstChildElement();
  for(int i=0; i < id; i++)opt = opt.nextSiblingElement();
  if(opt.attribute("checked").toInt()== 1)opt.setAttribute("checked", 0);
  else opt.setAttribute("checked", 1);
  emit textChanged(currentText());
}
void OpGroup::updateCurrent(const QString& c){
  if(myelem.attribute("element")=="string"
     ||myelem.attribute("element")=="float"
     ||myelem.attribute("element")=="int")myelem.setAttribute("current", c);
  else if(myelem.attribute("element")=="dvector"){
    QString tmp;
    for(int i = 0; i < mfs.size(); i++){
      tmp += qobject_cast<FloatEdit*>(mfs[i])->text();
      if(i != mfs.size()-1)tmp +=",";
    }
    myelem.setAttribute("current", tmp);
  }

    
  emit textChanged(currentText());
}


void OpGroup::updateCurrentX(const QString& c){
  myelem.setAttribute("currentX", c);
  emit textChanged(currentText());
}
void OpGroup::updateCurrentY(const QString& c){
  myelem.setAttribute("currentY", c);
  emit textChanged(currentText());
}

void OpGroup::updateCurrentZ(const QString& c){
  myelem.setAttribute("currentZ", c);
  emit textChanged(currentText());
}
void OpGroup::updateCurrentUnit(const QString& c){
  myelem.setAttribute("currentUnit", c);
   emit textChanged(currentText());
}
void OpGroup::updateUnitList(const QString& c){
  myelem.setAttribute("unit", myelem.attribute("unit")+","+c);
  myelem.setAttribute("currentUnit", c);
  emit textChanged(currentText());
}
QString  OpGroup::currentText(){
  
  QString text;
  if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";

  QString prefix = "=";
  if(myelem.hasAttribute("prefix")) prefix = myelem.attribute("prefix");
  QString name = myelem.tagName();
  if(myelem.hasAttribute("name")) name = myelem.attribute("name");

  if(name!="" && prefix!="" && myelem.attribute("element")!="selection" && myelem.attribute("element")!="novalue") text += name + prefix;

  QString unit="";
  if(myelem.hasAttribute("currentUnit"))unit = myelem.attribute("currentUnit") ;
  
  if(myelem.attribute("element")=="float" ||myelem.attribute("element")=="int"){
    if(!myelem.attribute("current").isEmpty()) myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status"); 
    text += myelem.attribute("current") + unit;
  }else if(myelem.attribute("element")=="vector"){
    if(myelem.attribute("format")=="cartesian" ){
      text +=  tr("[") + myelem.attribute("currentX")+unit+tr(", ") +
        myelem.attribute("currentY")+unit+tr(", ")+ myelem.attribute("currentZ")+unit+tr("]");    
    }else if(myelem.attribute("format")=="scale" ){
      text +=  tr("[") + myelem.attribute("currentX")+unit+tr(", 0, 0]"); 
    }else{
      text +=  tr("polar(") + myelem.attribute("currentX")+unit+tr(", ") +
        myelem.attribute("currentY")+unit+tr(", ")+ myelem.attribute("currentZ")+unit+tr(")"); 
    }
    if(myelem.attribute("format")=="scale"){
      if(!myelem.attribute("currentX").isEmpty())myelem.setAttribute("status", "done");
      else myelem.removeAttribute("status");
    }else if(!myelem.attribute("currentX").isEmpty() &&
             !myelem.attribute("currentY").isEmpty() &&
             !myelem.attribute("currentZ").isEmpty()) myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status"); 

  }else if(myelem.attribute("element")=="dvector"){
    double total = 0;
    text += "[";
    for(int i = 0; i < labels.size(); i++){
      text += labels[i]->text() + tr("= %1").arg(qobject_cast<FloatEdit*>(mfs[i])->value());
      if(i != labels.size() -1) text+=", ";
      else text += "]";
      total += qobject_cast<FloatEdit*>(mfs[i])->value();
    }
    if(total > 1.00000000000001){
      QMessageBox::warning(this, "OptionGroup", tr("the sum is greater than 1!"));
    }
    if((myelem.attribute("conditionSatisfied")=="false")||
      fabs(total-1)<0.00001) myelem.setAttribute("status", "done");
    else  myelem.removeAttribute("status"); 
  }else if(myelem.attribute("element")=="selection"){
    
    for(QDomElement opt=myelem.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){
      if(opt.attribute("checked").toInt()== 1)text+=opt.tagName()+ ", ";
    }
    
    if(text.size() > 2) text.remove(text.size()-2, 1);  
     myelem.setAttribute("status", "done"); 
  }else if(myelem.attribute("element")=="novalue"){
    text +=  myelem.tagName(); 
    myelem.setAttribute("status", "done");   
  }else{
    text +=  myelem.attribute("current");
    if(!myelem.attribute("current").isEmpty()) myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status"); 

  }
  
 
  
  myelem.setAttribute("currentText", text);
  return text;
}

void OpGroup::changeState(){
 
  if(myelem.hasAttribute("condition")){
    
    QString condition = myelem.attribute("condition");
    bool isSatisfied = conditionIsSatisfied(myroot, condition);
    
    if(isSatisfied){
      myelem.setAttribute("conditionSatisfied", "true");
       if(!myelem.hasAttribute("index")) this->show();
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
       if(!myelem.hasAttribute("index")) this->hide();
    }
  }else if(myelem.hasAttribute("defaultCondition")){
    QString condition = myelem.attribute("defaultCondition");
    bool isSatisfied = conditionIsSatisfied(myroot, condition);
    if(isSatisfied){
      myelem.setAttribute("conditionSatisfied", "true");
     
      this->setDefault(true);
       myelem.setAttribute("current", myelem.attribute("default"));
      //  this->setEnabled(false);
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
      //this->setEnabled(true);
      this->setDefault(false);
    }
  }
}

void OpGroup::setDefault(bool satisfied){
  if(myelem.attribute("element")=="float"){
    if(mfs.size() >= 1){
      if(satisfied)qobject_cast<FloatEdit*>(mfs[0])->setValue(myelem.attribute("default").toDouble());
     
      qobject_cast<FloatEdit*>(mfs[0])->setEnabled(!satisfied);
    }
  }else if(myelem.attribute("element")=="int"){
    if(mfs.size() >= 1){
      if(satisfied) qobject_cast<QSpinBox*>(mfs[0])->setValue(myelem.attribute("default").toInt());
      qobject_cast<QSpinBox*>(mfs[0])->setEnabled(!satisfied);
    }
  }else if(myelem.attribute("element")=="string"){
    if(mfs.size() >= 1){
      if(satisfied)qobject_cast<QLineEdit*>(mfs[0])->setText(myelem.attribute("default"));
      qobject_cast<QLineEdit*>(mfs[0])->setEnabled(!satisfied);
    }
  }else if(myelem.attribute("element")=="selection"){
    QStringList defaultList =myelem.attribute("default").split(",");
    if(satisfied){
      for(int j = 0; j < mfs.size(); j++){
        qobject_cast<QCheckBox*>(mfs[j])->setChecked(false);
      }
      for(int i =0; i < defaultList.size(); i++){
        bool found = false;
        for(int j = 0; j < mfs.size(); j++){
          if( qobject_cast<QCheckBox*>(mfs[j])->text()==defaultList[i]){
            if(satisfied) qobject_cast<QCheckBox*>(mfs[j])->setChecked(true);
            found = true;
            break;
          }
        }
        if(!found){
          QMessageBox::warning(window(), "OptionGroup", tr("illegal default value ")+defaultList[i] + " in " +myelem.tagName());
        }
      }
    }
    for(int j = 0; j < mfs.size(); j++){
      qobject_cast<QCheckBox*>(mfs[j])->setEnabled(!satisfied);
    }
    
  }else if(myelem.attribute("element")=="vector"){
    QStringList defaultList =myelem.attribute("default").split(",");
    if(defaultList.size()<=3 && mfs.size()==3){
      for(int i =0; i < defaultList.size(); i++){
        if(satisfied) qobject_cast<FloatEdit*>(mfs[i])->setText(defaultList[i]);
        qobject_cast<FloatEdit*>(mfs[i])->setEnabled(!satisfied); 
      }
    }
  }
  emit textChanged(currentText());
}



void OpGroup::updateComponents(){

  if(myelem.attribute("element")!="dvector")return;
    QStringList components = myroot.firstChildElement("components").attribute("components").split(",");
   
    QLayout* theLayout=layout();
    if(theLayout) delete theLayout;
    for(int i = 0; i < mfs.size(); i++){
      if(mfs[i])delete mfs[i];
    }
    for(int i = 0; i < labels.size(); i++){
      if(labels[i])delete labels[i];
    }
    
    mfs.clear();
    labels.clear();
    
    QVBoxLayout* myLayout = new QVBoxLayout;
    setLayout(myLayout);
   
    if(components.size()!=0){
       for(int i =0; i < components.size(); i++){
    
         QPointer<QLabel> aLabel = new QLabel(components[i]);
         QPointer<FloatEdit> aEditor = new FloatEdit;
         QHBoxLayout* hLayout=new QHBoxLayout;
           hLayout->addWidget(aLabel);
          hLayout->addWidget(aEditor);
         connect(aEditor, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrent(const QString&)));
         aLabel->setAlignment(Qt::AlignHCenter);
         labels<<aLabel;
         mfs<<aEditor;
         myLayout->addLayout(hLayout);
        
      }
   
    this->show();
  }else{  
    
    myelem.setAttribute("conditionSatisfied", "false");
    this->hide();
  }
    
  emit textChanged(currentText());
}




  

void OpGroup::updateLabels(int i){
  switch(i){
  case 0:
    labels[0]->setText("x");
    labels[1]->setText("y");
    labels[2]->setText("z");

    myelem.setAttribute("format", "cartesian");
      for(int i =1; i < 3; i++){
        qobject_cast<FloatEdit*>(mfs[i])->setEnabled(true); 
      }
    break;
  case 1:
    labels[0]->setText("mag");
    labels[1]->setText("theta");
    labels[2]->setText("phi");

    
    myelem.setAttribute("format", "polar");
     for(int i =1; i < 3; i++){
        qobject_cast<FloatEdit*>(mfs[i])->setEnabled(true); 
      }
    break;
  case 2:
    labels[0]->setText("x");
    labels[1]->setText("y");
    labels[2]->setText("z");

    

    
    myelem.setAttribute("format", "scale");
    for(int i =1; i < 3; i++){
      qobject_cast<FloatEdit*>(mfs[i])->setValue(0.0); 
      qobject_cast<FloatEdit*>(mfs[i])->setEnabled(false); 
    }
    break;
    
 default:
   QMessageBox::warning(window(), "OptionGroup", "i out of range");
  }
   emit textChanged(currentText());
}





          



AllVBGroup::AllVBGroup(  QDomElement& elem, QDomElement& theroot, bool isVertical, QWidget *parent )
  : GeneralGroup(elem, theroot, parent){
  
 

    signalMapper = new QSignalMapper(this);
  QBoxLayout *mainLayout;
  if(isVertical)mainLayout = new QVBoxLayout;
  else mainLayout = new QHBoxLayout;
  
  QDomElement elem_opt = myelem.firstChildElement();
  if(elem_opt.isNull()){
    QMessageBox::warning(window(), tr(" .xml"),
                         tr("can not find any child element ")+myelem.tagName()+tr(" groups")
                         );
    return;
  }
  
  int count=0;
  for (; !elem_opt.isNull(); elem_opt = elem_opt.nextSiblingElement()) {
    
    //if element is specified, take it, otherwise go to options to find and replace it
    if(!elem_opt.hasAttribute("element")){
      QString title = elem_opt.attribute("title");
      QString text = elem_opt.attribute("name");
      
      
      QDomElement option_elem = myroot.firstChildElement("options");
      if(option_elem.isNull()){
        QMessageBox::warning(window(), tr(" .xml"),
                             tr("can not find element 'options'")
                             );
        return;
      }
      QDomElement newChild = option_elem.firstChildElement(elem_opt.tagName()).cloneNode().toElement();
      if(newChild.isNull()){
        QMessageBox::warning(window(), tr(" .xml"),
                             tr("'options' has no child ")+elem_opt.tagName()
                             );
        return;
      }
     
      QDomNode tmp= myelem.replaceChild(newChild,elem_opt);
      if(tmp.isNull()){
        QMessageBox::warning(window(), tr(" .xml"),
                             tr("'replace ")+elem_opt.tagName()+tr("failed")
                             );
        return;
      }
      elem_opt=newChild;
      //put the textNode of oldChild into the newChild's attribute
      
      if(text!=""){
        elem_opt.setAttribute("name", text);
        if(elem_opt.attribute("element")=="vector") elem_opt.setAttribute("title", text);
      }
      if(title!="")elem_opt.setAttribute("title", title);
      
      
    }
    
    
    if(elem_opt.hasAttribute("defaultCondition") ){
      QString tmp =  elem_opt.attribute("defaultCondition");
      bool  conditionSatisfied = conditionIsSatisfied(theroot, tmp);
      if(conditionSatisfied){
        QString tmp = elem_opt.attribute("default");
        elem_opt.setAttribute("current", tmp);
      }
    }
    
    
    OpGroup *myGroup = new OpGroup(elem_opt, myroot, this);
    myGroup->setFlat(true);
    if(myelem.attribute("element")=="all"){
      myGroup->setCheckable(false);
      if(elem_opt.hasAttribute("checked")){
        QMessageBox::warning(window(), tr(" .xml"),
                             elem_opt.tagName()+tr("should not have attribute 'checked'")
                             );
        elem_opt.removeAttribute(tr("checked"));
      }
    }else if(myelem.attribute("element")=="selection"){
      bool isChecked = elem_opt.hasAttribute("checked") && elem_opt.attribute("checked").toInt()== 1; 
      myGroup->setCheckable(true);
      if(isChecked){
        myGroup->setChecked(true);
      }else{
        myGroup->setChecked(false);
        
      }
    }else if(myelem.attribute("element")== "2of3"){
      //tricky here, because for QGroupBox, setCheckable will autoamtically setChecked
      bool isChecked = elem_opt.hasAttribute("checked") && elem_opt.attribute("checked").toInt()== 1;
      myGroup->setCheckable(true);
      if(isChecked ){
        myGroup->setChecked(true);
      }else {
        myGroup->setChecked(false);
      }
    }else{
      QMessageBox::warning(window(), tr(" .xml"),
                           tr("Don't know how to handle ")+elem_opt.tagName() + " " +elem_opt.attribute("element"));
    }   
    
    connect(myGroup, SIGNAL(clicked()),signalMapper, SLOT(map()));
    
    objs<<myGroup;  
    signalMapper->setMapping(myGroup, count);        
    connect(myGroup, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
    connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
    connect(this, SIGNAL(componentsChanged()), myGroup, SLOT(updateComponents()));
    myGroup->setFlat(true);
    mainLayout->addWidget(myGroup);
    count++;
    
    if(elem_opt.hasAttribute("condition")){
      QString tmp = myelem.attribute("condition");
      bool conditionSatisfied = conditionIsSatisfied(theroot,tmp);
      if(conditionSatisfied){
        myGroup->show();
        elem_opt.setAttribute("checked", 1);
      }else{
        myGroup->hide();
        elem_opt.setAttribute("checked", 0);
      }
    }
    if(elem_opt.hasAttribute("defaultCondition") ){
      bool conditionSatisfied = conditionIsSatisfied(theroot, elem_opt.attribute("defaultCondition"));
      if(conditionSatisfied){
        myGroup->setDisabled(true);
      }else{
        myGroup->setDisabled(false);
      }
      
    }
    
    
  }//for(; elt..)
  
  connect(signalMapper, SIGNAL(mapped(int)),
          this, SLOT(childrenClicked(int)));  
  
  emit textChanged(currentText()); 
  setLayout(mainLayout);
}


void AllVBGroup::childrenClicked(int i){
  
  if(myelem.attribute("element")=="2of3"){
    bool isChecked =  qobject_cast<OpGroup*>(objs[i])->isChecked();
    if(isChecked){
      if( ( qobject_cast<OpGroup*>(objs[(i+1)%3])->isChecked())
          && ( qobject_cast<OpGroup*>(objs[(i+2)%3])->isChecked())){
        qobject_cast<OpGroup*>(objs[i])->setChecked(false);
      }
    }
  }
  
  
  emit textChanged(currentText());
}

QString AllVBGroup::currentText(){
  QString infix=", ";
  if(myelem.hasAttribute("infix")) infix = myelem.attribute("infix");
  
  QString text ;
if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";
  
  for(QDomElement opt = myelem.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){ 
    if(!opt.hasAttribute("checked")||opt.attribute("checked").toInt()== 1){
      if(opt.hasAttribute("condition")){
        if(opt.attribute("conditionSatisfied")=="true") text += opt.attribute("currentText")+ infix;
      }else{
        text += opt.attribute("currentText")+ infix;
      }
    }
  }
  text.remove(text.size()-infix.size(), infix.size());
  myelem.setAttribute("currentText", text);
  if(myelem.attribute("element")=="2of3"){
    int count = 0;
    for(QDomElement opt = myelem.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){  
      if(!opt.hasAttribute("checked")||opt.attribute("checked").toInt()== 1){
        if(opt.attribute("status")=="done") count++;
      }
    }
    if(count == 2) myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status");
  }else{
    bool done = true;
    for(QDomElement opt = myelem.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){ 
      if(!opt.hasAttribute("checked")||opt.attribute("checked").toInt()== 1){
        if(opt.attribute("status")!="done") done=false;
      }
    }
    if(done) myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status");
  }
  return text;
  
}

void AllVBGroup::updateCurrentText(){
  
  emit textChanged(currentText());
}


  

StackGroup::StackGroup( QDomElement& elem, QDomElement& root, QWidget *parent )
  : GeneralGroup(elem, root, parent){

  
  typesWidget = new QComboBox;
  pagesWidget = new QStackedWidget;
 
  QDomElement elem_grp = myroot.firstChildElement("groups");
  
  if(elem_grp.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'groups' in the children of root")
                         );
    return;
  }
 
 
 
 QRegExp rx("*Group");
 rx.setPatternSyntax(QRegExp::Wildcard);

 QDomElement elem_opt = myelem.firstChildElement();
 if(elem_opt.isNull()){
   QMessageBox::warning(window(), tr(".xml"),
                        myelem.tagName()+ tr(" has no child ")
                        );
   return;
 }
 
 for (; !elem_opt.isNull(); elem_opt = elem_opt.nextSiblingElement()) {
   QString title = elem_opt.attribute("title");
   QString whatsThis=elem_opt.attribute("whatsThis");
   QString toolTip = elem_opt.attribute("toolTip");
   QString name =  elem_opt.attribute("name");  
  
   
   
   QGroupBox *myGroup;
   //if elem_opt is a group
   if(rx.exactMatch(elem_opt.tagName())){
     //if the elem_opt is  a group and it needs copy
     if(elem_opt.firstChildElement().isNull()){
       QDomElement elem_opt2 = elem_grp.firstChildElement(elem_opt.tagName()).cloneNode().toElement();
       if(elem_opt2.isNull()){
         QMessageBox::warning(window(), tr(".xml"),
                              tr("can not find element ") + elem_opt.tagName() + tr(" in the children of groups")
                              );
         return;
       }
       QDomNode tmp= myelem.replaceChild(elem_opt2, elem_opt);
       if(tmp.isNull()){
         QMessageBox::warning(window(), tr(" .xml"),
                              tr("'replace ")+elem_opt2.tagName()+tr("failed")
                              );
         return;
       }
       elem_opt=elem_opt2;
       if(title!="")elem_opt.setAttribute("title", title);
       if(name!="")elem_opt.setAttribute("name", name);
       if(whatsThis!="")elem_opt.setAttribute("whatsThis", whatsThis);
       if(toolTip!="")elem_opt.setAttribute("toolTip", toolTip);
       
     }//finish copy group
     
     if(elem_opt.attribute("element")=="all"||elem_opt.attribute("element")=="selection"||elem_opt.attribute("element")=="2of3"){
       myGroup = new AllVBGroup(elem_opt, myroot, false,this);
       connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
     }else{
       qDebug()<<"don't know how to handle it yet";
     }
   }else{//if it's not a group
     //if the options needs copy
     if(!elem_opt.hasAttribute("element")){
       QDomElement elem_options = myroot.firstChildElement("options");
       if(elem_options.isNull()){
         QMessageBox::warning(window(), tr(".xml"),
                              tr("can not find element 'options'")
                              );
         return;
       }
       QDomElement elem_opt2 = elem_options.firstChildElement(elem_opt.tagName()).cloneNode().toElement();
       if(elem_opt2.isNull()){
         QMessageBox::warning(window(), tr(".xml"),
                              tr("can not find element ") + elem_opt.tagName() + tr(" in the children of 'options'")
                              );
         return;
       }
       
       // QString name = elem_opt.text();
       QDomNode tmp= myelem.replaceChild(elem_opt2, elem_opt);
       if(tmp.isNull()){
         QMessageBox::warning(window(), tr(" .xml"),
                              tr("'replace ")+elem_opt2.tagName()+tr("failed")
                              );
         return;
       }
       elem_opt=elem_opt2;
       if(name!="") elem_opt.setAttribute("name", name);
       if(title!="")elem_opt.setAttribute("title", title);
       if(whatsThis!="")elem_opt.setAttribute("whatsThis", whatsThis);
       if(toolTip!="")elem_opt.setAttribute("toolTip", toolTip);
     
     }//finish copy option
     myGroup = new OpGroup(elem_opt, myroot, this);
     if(elem_opt.attribute("element")=="novalue"){
       myGroup->setTitle("");
       myGroup->setFlat(true);
     }
     connect(this, SIGNAL(componentsChanged()), myGroup, SLOT(updateComponents()));
   }
   
   typesWidget->addItem(elem_opt.attribute("title"));
   whatsThisList << elem_opt.attribute("whatsThis");
   toolTipList << elem_opt.attribute("toolTip");
   connect(myGroup, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
   connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
   
   myGroup->setCheckable(false);
   pagesWidget->addWidget(myGroup);
  
 }//for(; elem_opt..)

 connect(typesWidget,
         SIGNAL(currentIndexChanged(int)),
         this, SLOT(changePage(int)));
 
 if(myelem.hasAttribute("currentIndex"))typesWidget->setCurrentIndex(myelem.attribute("currentIndex").toInt());
 else typesWidget->setCurrentIndex(0);

 if(!(whatsThisList[typesWidget->currentIndex()].isNull()))typesWidget->setWhatsThis(whatsThisList[typesWidget->currentIndex()]);
 if(!(toolTipList[typesWidget->currentIndex()].isNull()))typesWidget->setToolTip(toolTipList[typesWidget->currentIndex()]);
    

 QHBoxLayout *mainLayout = new QHBoxLayout;
 mainLayout->addWidget(typesWidget);
 mainLayout->addWidget(pagesWidget);

 QVBoxLayout *vboxLayout = new QVBoxLayout;
 vboxLayout->addLayout(mainLayout);
 setLayout(vboxLayout);
 emit(textChanged(currentText()));
}



QString StackGroup::currentText(){
  
  int current = myelem.attribute("currentIndex").toInt();
  QString text ;
  if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";
  QDomElement opt = myelem.firstChildElement();
  for(int i=0; i < current; i++)if(!opt.isNull()) opt = opt.nextSiblingElement();
  if(opt.attribute("status")=="done")myelem.setAttribute("status", "done");
  else myelem.removeAttribute("status");
  text =  opt.attribute("currentText");
  myelem.setAttribute("currentText", text);
  return text;
}


void  StackGroup::clearText(){

}
     
void StackGroup::changePage(int index)
{
  
  pagesWidget->setCurrentIndex(index);
  
  myelem.setAttribute("currentIndex", index);
  QString tmp = currentText();
  emit textChanged(tmp);
  if(!(whatsThisList[index].isNull()))typesWidget->setWhatsThis(whatsThisList[index]);
  else typesWidget->setWhatsThis("");
  if(!(toolTipList[index].isNull()))typesWidget->setToolTip(toolTipList[index]);
  else typesWidget->setToolTip("");
}

void StackGroup::updateCurrentText(){
  
  emit textChanged(currentText());
}


OptionPage::OptionPage(   QDomElement& my_elem, QDomElement &my_root, QWidget *parent )
  :GeneralGroup(my_elem, my_root, parent)
{
  if(!toolTip().isEmpty())setToolTip("");
  if(!statusTip().isEmpty())setStatusTip("");
  if(!whatsThis().isEmpty())setWhatsThis("");
 

  buttonGroup = 0;  
  QDomElement elt = myelem.firstChildElement();
  if(elt.isNull()){
    QVBoxLayout *mainLayout = new QVBoxLayout;
    QLabel* aLabel = new QLabel("No Options");
    aLabel->setAlignment(Qt::AlignCenter);
    mainLayout->addWidget(aLabel);
    setLayout(mainLayout);
    updateCurrentText();
    return;
  }

  int numColumn= 1;
  if(myelem.hasAttribute("numColumn")) numColumn = myelem.attribute("numColumn").toInt();
  QGridLayout *mainLayout = new QGridLayout;
  QHBoxLayout* advancedLayout = new QHBoxLayout;
  for(int i = 0; i < numColumn; i++){
    mainLayout->setColumnMinimumWidth(i, 30);
    mainLayout->setColumnStretch(i, 1);
  }
  QRegExp rx("*Group");
  rx.setPatternSyntax(QRegExp::Wildcard);
  
  QDomElement elem_grp = myroot.firstChildElement("groups");

  buttonGroup = new QButtonGroup(this);
  buttonGroup->setExclusive(false);
       
  
  int elt_count=0;
  int button_count = 0;
  for (; !elt.isNull(); elt = elt.nextSiblingElement(), elt_count++) {
   QPointer<QWidget> myGroup = 0;
   QPointer<QPushButton> myButton = 0;
    if(rx.exactMatch(elt.tagName())){//groups

       if(elem_grp.isNull()){
         QMessageBox::warning(window(), tr(".xml"),
                              tr("can not find element 'groups' in the children of root")
                              );
         return;
       }
      if(elt.firstChildElement().isNull()){//need copy
        QString title = elt.attribute("title");
        
        QDomElement elem_opt2 = elem_grp.firstChildElement(elt.tagName()).cloneNode().toElement();
        if(elem_opt2.isNull()){
          QMessageBox::warning(window(), tr(".xml"),
                               tr("can not find element ") + elt.tagName() +
                               tr(" in the children of here ") + elem_grp.tagName()
                               );
          return;
        }
     
        QDomNode tmp= myelem.replaceChild(elem_opt2, elt);
        if(tmp.isNull()){
          QMessageBox::warning(window(), tr(" .xml"),
                               tr("'replace ")+elem_opt2.tagName()+tr("failed")
                               );
          return;
        }
        elt=elem_opt2;
         if(title!="")elt.setAttribute("title", title);
      }
      if(elt.attribute("element")=="all" ||elt.attribute("element")=="2of3"  ){
        myGroup = new AllVBGroup(elt, myroot, false, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="stack"){
        myGroup = new StackGroup(elt, myroot, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="choice"){
        myGroup = new ChoiceGroup(elt, myroot, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));

      }else {
        QMessageBox::warning(window(), tr(".xml"),
                             tr("don't know how to handle it yet ") +elt.tagName() + " " + elt.attribute("element")) ;
      }
      if(elt.hasAttribute("buttonTitle")){
        myButton = new QPushButton(elt.attribute("buttonTitle"), this);
        buttonGroup->addButton(myButton);
        buttonGroup->setId(myButton, button_count);
        myButton->setMaximumWidth(200);
        button_count++;
        elt_count--;
        advancedLayout->addWidget(myButton);
        myAdvancedGroup<<myGroup;
        myGroup->hide();
      }
    }else {//options, panels
      
      if(!elt.hasAttribute("element")){//need copy
        QString title = elt.attribute("title");
        QString text = elt.attribute("name");
       
        QDomElement option_elem = myroot.firstChildElement("options");
        if(option_elem.isNull()){
          QMessageBox::warning(window(), tr(" .xml"),
                               tr("can not find element 'options'")
                               );
          return;
        }
        QDomElement newChild = option_elem.firstChildElement(elt.tagName()).cloneNode().toElement();
        if(newChild.isNull()){
          QMessageBox::warning(window(), tr(" .xml"),
                               tr("'options' has no child ")+elt.tagName()
                               );
          return;
        }
        
        QDomNode tmp= myelem.replaceChild(newChild,elt);
        if(tmp.isNull()){
          QMessageBox::warning(window(), tr(" .xml"),
                               tr("'replace ")+elt.tagName()+tr("failed")
                               );
          return;
        }
        elt=newChild;
        //put the textNode of oldChild into the newChild's attribute
        
        if(text!=""){
          elt.setAttribute("name", text);
          if(elt.attribute("element")=="vector") elt.setAttribute("title", text);
        }
        if(title!="")elt.setAttribute("title", title);
 
      
      }
      if(elt.attribute("element")=="choice"){
        myGroup = new ChoiceGroup(elt, myroot);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="all"){
        myGroup = new AllVBGroup(elt, myroot);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
         connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="stack"){
        myGroup = new StackGroup(elt, myroot);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
         connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="panel"){
        myGroup = new OptionPage(elt, myroot);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
         connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else{
        myGroup = new OpGroup(elt, myroot);
        connect(this, SIGNAL(componentsChanged()), myGroup, SLOT(updateComponents()));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }
     
      if(elt.hasAttribute("buttonTitle")){
        myButton = new QPushButton(elt.attribute("buttonTitle"), this);
        buttonGroup->addButton(myButton);
        buttonGroup->setId(myButton, button_count);
        myButton->setMaximumWidth(200);
        button_count++;
        elt_count--;
        advancedLayout->addWidget(myButton);
        myAdvancedGroup<<myGroup;
        myGroup->hide();
      }
     
    }
    if(!elt.hasAttribute("buttonTitle")) mainLayout->addWidget(myGroup, elt_count/numColumn, elt_count%numColumn, 1, 1);   
    connect(myGroup, SIGNAL(textChanged(const QString&)),this, SLOT(updateCurrentText())); 
  }
 if(buttonGroup) connect(buttonGroup, SIGNAL(buttonClicked(int)),this, SLOT(advancedButtonClicked(int))); 
  
  mainLayout->addLayout(advancedLayout, (elt_count/numColumn +1), 0, 1, numColumn);
  
    
  if(myelem.attribute("feedback")!="off"){
    currentLabel = new QLabel;
    currentLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    currentLabel->setWordWrap(true);
    if(myelem.hasAttribute("currentText"))currentLabel->setText(myelem.attribute("currentText"));
    connect(this, SIGNAL(textChanged(const QString&)), currentLabel, SLOT(setText(const QString&)));
    mainLayout->addWidget(currentLabel, (elt_count/numColumn +2), 0, 1, numColumn);
  }
  setLayout(mainLayout);
  
 
 
  updateCurrentText();

   
}
  


void OptionPage::clearText(){
}


void OptionPage::updateCurrentText(){
  QString tmp = currentText();
  myelem.setAttribute("currentText", tmp);
  emit textChanged(tmp);
}

  
QString OptionPage::currentText(){
  if(myelem.firstChildElement().isNull() && myelem.tagName()!="notset"){
    myelem.setAttribute("status","done");
    
    if(myelem.hasAttribute("buttonIndex"))emit(updateStatusTip(myelem.attribute("buttonIndex").toInt()));
    return myelem.tagName();
  }else if(myelem.tagName()=="notset"){
    
    myelem.setAttribute("status","new");
    return QString();
  }
  
  QString text ;
  if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";
  QString prefix ="";
  if(myelem.hasAttribute("prefix")) prefix = myelem.attribute("prefix");
  QString postfix = "";
  if(myelem.hasAttribute("postfix")) postfix = myelem.attribute("postfix");
  QString infix = ", ";
  if(myelem.hasAttribute("infix")) infix = myelem.attribute("infix");
  else if(prefix=="" && postfix=="") infix = "\n";
  
  
 
  for(QDomElement elt = myelem.firstChildElement(); !elt.isNull(); elt=elt.nextSiblingElement()){
    if(prefix=="" &&postfix =="" && elt.attribute("element")!="panel")elt.setAttribute("prefix", ": ");
    if(elt.attribute("element")!="advanced"){
      if((!elt.hasAttribute("checked"))||elt.attribute("checked")== "1"){
        if(elt.attribute("currentText")!="" ){
          if(elt.hasAttribute("condition")){
            if(elt.attribute("conditionSatisfied")=="true")text += elt.attribute("currentText")+infix;
          }else{
            text += elt.attribute("currentText")+infix;
          }
           
        }      
      }
    }else{
      for(QDomElement opt = elt.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){
        if(opt.attribute("checked").toInt()==1){
          if(opt.hasAttribute("condition")){
            if(opt.attribute("conditionSatisfied")=="false")
              text += opt.attribute("currentText")+infix;
          }else{
            text += opt.attribute("currentText")+infix;
          }
        }
      }
    }
  }
  text.remove(text.size()-infix.size(), infix.size());

  QString name = myelem.tagName();
  if(myelem.hasAttribute("name"))name=myelem.attribute("name");
  if(prefix !="") text = name + prefix + text + postfix;

  int count = 0;
  int count_done = 0;
  for(QDomElement elt = myelem.firstChildElement(); !elt.isNull(); elt=elt.nextSiblingElement()){
    if(elt.attribute("element")!="advanced"){
      if((!elt.hasAttribute("checked"))||elt.attribute("checked")== "1"){
        count++;
        if(elt.attribute("status")=="done") count_done++;
      }
    }else{
      for(QDomElement opt = elt.firstChildElement(); !opt.isNull(); opt=opt.nextSiblingElement()){
        if(opt.attribute("checked").toInt()==1){  
          count++;
          if(opt.attribute("status")=="done") count_done++;
        }
      }
    }
  } 
 
 
  QString tmp =  QString("%1 out of ").arg(count_done)+ QString("%1 finished").arg(count);
 
  if(count==count_done)myelem.setAttribute("status", "done");
  else myelem.setAttribute("status",tmp); 
  if(myelem.hasAttribute("buttonIndex"))emit(updateStatusTip(myelem.attribute("buttonIndex").toInt()));
  
  return text;



}
  
void OptionPage::advancedButtonClicked(int index){
  if(index<myAdvancedGroup.size())  myAdvancedGroup[index]->show();
}

OptionPage::~OptionPage(){
  for(int i= 0; i < myAdvancedGroup.size(); i++){
    if(myAdvancedGroup[i])delete myAdvancedGroup[i];
    myAdvancedGroup[i] = 0;
  }
  if(buttonGroup){
    delete buttonGroup;
    buttonGroup = 0;
  }
}









ChoiceGroup::ChoiceGroup(   QDomElement& my_elem,  QDomElement& theroot, QWidget *parent )
  : GeneralGroup(my_elem, theroot, parent)
{

  editGroup = 0;
   int currentChoice = 0;
  if(myelem.hasAttribute("current")) currentChoice=myelem.attribute("current").toInt();
  else myelem.setAttribute("current", 0);
 
  comboBox = new QComboBox;
 
 
  editButton = new QPushButton;
  editButton->hide();
  QHBoxLayout *mainLayout = new QHBoxLayout;
  mainLayout->addWidget(comboBox);
  mainLayout->addWidget(editButton);
  setLayout(mainLayout);
  QDomElement elt = myelem.firstChildElement();
  if(elt.isNull()){
    return;
  }
  
    
  
  QStringList choices;
  int buttonCount = 0;
  
  for (; !elt.isNull(); elt = elt.nextSiblingElement(), buttonCount++) {
  
    choices<<elt.tagName();
    whatsThisList << elt.attribute("whatsThis");
    toolTipList << elt.attribute("toolTip");
    
    if(!elt.attributeNode("editable").isNull()){
      
      if(buttonCount==currentChoice){
        QDomElement edit_elem = myroot.firstChildElement("editable");
        edit_elem = edit_elem.firstChildElement(elt.attribute("editable"));
        
        if(edit_elem.hasAttribute("title"))editButton->setText(edit_elem.attribute("title"));
        else editButton->setText(elt.attribute("editable"));
        editButton->show();
      }
      editItems <<elt.attribute("editable");
    }else{
      editItems <<"";;
    }
  }
  comboBox->addItems(choices);
  comboBox->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLength);
  comboBox->setInsertPolicy(QComboBox::InsertAtBottom);
  comboBox->setCurrentIndex(currentChoice);
  if(!(whatsThisList[currentChoice].isNull()))comboBox->setWhatsThis(whatsThisList[currentChoice]);
  if(!(toolTipList[currentChoice].isNull()))comboBox->setToolTip(toolTipList[currentChoice]);
  
  
  connect(comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(update(int)));
  connect(editButton, SIGNAL(clicked()),
          this, SLOT(editButtonPressed()));  
  updateCurrentText();
}
ChoiceGroup::~ChoiceGroup(){
  if(editGroup){
    delete editGroup;
    editGroup = 0;
  }
}


QString ChoiceGroup::currentText(){
  int i = comboBox->currentIndex();
 
  QDomElement elt = myelem.firstChildElement();
  for(int j = 0; j < i; j++)elt=elt.nextSiblingElement();
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         myelem.tagName()+tr(" has no ") + i + tr("th child")
                         );
    return QString();
  }
  
  
  QString name;
  if(myelem.hasAttribute("name"))name=myelem.attribute("name");
  else name=myelem.tagName();

  QString prefix =": ";
  if(myelem.hasAttribute("prefix"))prefix=myelem.attribute("prefix");

  QString tmp;

  if(myelem.hasAttribute("comments")) tmp += myelem.attribute("comments")+"\n";
   tmp += name + prefix + elt.tagName();
  
  
  if(elt.hasAttribute("editable")){
    QDomElement elt_editable= myroot.firstChildElement("editable");
    if(elt_editable.isNull()){
      QMessageBox::warning(window(), tr(".xml"),
                           tr("can not find element 'editable' in the children of root")
                           );
      return QString();
    }
    
    QDomElement elt_item = elt_editable.firstChildElement(elt.attribute("editable"));
    if(elt_item.isNull()){
      QMessageBox::warning(window(), tr(".xml"),
                           tr("can not find element" )+editItems[comboBox->currentIndex()]+tr(" in the children of 'editable'")
                           );
      return QString();
    }
    
    if(!elt_item.hasAttribute("currentText")){
      if(editGroup){
        delete editGroup;
        editGroup =0;
      }
      if(elt_item.attribute("element")=="panel"){
        editGroup = new  OptionPage(elt_item, myroot);
      }else{
        editGroup = new OpGroup(elt_item, myroot);
      }
    }
      
    if((!(elt_item.hasAttribute("checked")) || elt_item.attribute("checked").toInt()==1)
       && elt_item.attribute("currentText")!=""){
      

      if(elt_item.hasAttribute("condition")){
        if(elt_item.attribute("conditionSatisfied")=="true"){
          tmp += "\n";
          tmp += elt_item.attribute("currentText");
        }
      }else{
        tmp += "\n";
        tmp += elt_item.attribute("currentText");
      }
    }

    
    if(elt_item.attribute("status")=="done")myelem.setAttribute("status", "done");
    else myelem.removeAttribute("status");
  }else{
    myelem.setAttribute("status", "done");
  }
  return tmp;
}
void ChoiceGroup::updateCurrentText(){
  QString tmp = currentText();
      
  myelem.setAttribute("currentText", tmp);
  emit textChanged(tmp);
}







void ChoiceGroup::update(int i){
  myelem.setAttribute("current", i);
  if(editItems[i]!=""){
    editButton->setText(editItems[i]);
    editButton->show();
  }else{
    editButton->hide();
  }
 
  if(!(whatsThisList[i].isNull()))comboBox->setWhatsThis(whatsThisList[i]);
  else comboBox->setWhatsThis("");
  if(!(toolTipList[i].isNull()))comboBox->setToolTip(toolTipList[i]);
  else comboBox->setToolTip("");

  myelem.setAttribute("currentText", currentText()); 
  emit textChanged(currentText());
}

void ChoiceGroup::editButtonPressed(){
  //no edit item, return
  if(editItems[comboBox->currentIndex()]=="") return;
  
  //go to editable to find the item
  QDomElement elt= myroot.firstChildElement("editable");
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element 'editable' in the children of root")
                         );
    return;
  }
  elt = elt.firstChildElement(editItems[comboBox->currentIndex()]);
  if(elt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element" )+editItems[comboBox->currentIndex()]+tr(" in the children of 'editable'")
                         );
    return;
  }
  
  if(editGroup){
    delete editGroup;
    editGroup =0;
  }
  if(elt.attribute("element")=="panel"){
    editGroup = new  OptionPage(elt, myroot);
    connect(this, SIGNAL(componentsChanged()), editGroup, SIGNAL(componentsChanged()));
    connect(this, SIGNAL(stateChanged()), editGroup, SLOT(changeState()));
  }else{
    editGroup = new OpGroup(elt, myroot);
    connect(this, SIGNAL(componentsChanged()), editGroup, SLOT(updateComponents()));
    connect(this, SIGNAL(stateChanged()), editGroup, SLOT(changeState()));
  }
  connect(editGroup, SIGNAL(textChanged(QString)), this, SLOT(updateCurrentText()));
  editGroup->show();
  updateCurrentText();
 
}
AllVBWindow::AllVBWindow(  QDomElement& elem,  QDomElement& root,  QWidget *parent):GeneralWindow(elem, root, parent){
  if(myelem.hasAttribute("whatsThis"))setWhatsThis(myelem.attribute("whatsThis"));
  if(myelem.hasAttribute("toolTip"))setToolTip(myelem.attribute("toolTip"));
  if(myelem.hasAttribute("statusTip"))setStatusTip(myelem.attribute("statusTip"));
  


  
  AllVBGroup* thegroup = new AllVBGroup(elem, root, this);
  connect(this, SIGNAL(componentsChanged()), thegroup, SIGNAL(componentsChanged()));
  connect(this, SIGNAL(stateChanged()), thegroup, SLOT(changeState()));
  
  QPushButton* saveButton = new QPushButton(tr("          save the module        "));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  
  
  QVBoxLayout* mainLayout = new QVBoxLayout;

  mainLayout->addWidget(thegroup);
  mainLayout->addWidget(saveButton);
  setLayout(mainLayout);
}


bool AllVBWindow::save(){
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                            "/home/qxue/chemdemo/untitled.mdl",
                            tr("module file (*.mdl)"));
  if (fileName.isEmpty())
    {
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
  
  QTextStream out(&file);
  QApplication::setOverrideCursor(Qt::WaitCursor);
  out<<"species={" << endl;
  if(myelem.tagName()=="specified_ideal_gas_model")out<<"_gas=<";
  else if(myelem.tagName()=="gamma_model")out<<"_Gamma"<<(int)(myelem.firstChildElement("gamma").attribute("current").toDouble()*10)<<"=< ";
  QDomElement elt = myelem.firstChildElement();
  for(; !elt.isNull(); elt=elt.nextSiblingElement()){
    out<<elt.tagName()<<"="<<elt.attribute("current")<<",";
  }
  if(myelem.tagName()=="specified_ideal_gas_model"||myelem.tagName()=="gamma_model"){
    out<<"mf=1," << endl;
    out<<"href=0, sref=0, Tref=298, Pref=1e5>;"<<endl;
    out<<"};"<<endl;
    out<<"reactions={"<<endl<<"};"<<endl;
  }
  
  QApplication::restoreOverrideCursor();
 file.close();
 myelem.setAttribute("current", fileName);
 return true;  
     
}
 
StackGroup2::StackGroup2(  QDomElement& elem, QDomElement& theroot, const QStringList& pState, QWidget *parent )
  : GeneralGroup(elem, theroot, parent), parentState(pState){
 

  name = elem.tagName();
  if(myelem.hasAttribute("name"))name = myelem.attribute("name");
  
  pagesWidget = new QStackedWidget;
  QDomElement elem_opt= elem.firstChildElement();
  QDomElement elem_models= myroot.firstChildElement("models");
  
  if(elem_opt.isNull()){
         QMessageBox::warning(window(), tr(".xml"),
                         myelem.tagName() + tr("stack element has no child")
                              );
         return;
  }
  if(elem_models.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element models")
                         );
    return;
   }
  
  QHBoxLayout* mainLayout= new QHBoxLayout;
 
  QButtonGroup* buttonGroup =  new QButtonGroup(this);
  QVBoxLayout* buttonLayout= new QVBoxLayout;
  mainLayout->addLayout(buttonLayout);
  int count=0;
  current = 0;
  if(elem.hasAttribute("current"))current=elem.attribute("current").toInt();
  else elem.setAttribute("current", 0);
 
  
 //each child is a button and a page 
 for (; !elem_opt.isNull(); elem_opt = elem_opt.nextSiblingElement(), count++) {
   
   //set up the page
   QGroupBox* aPage = new QGroupBox;
   aPage->setFlat(true);
   QHBoxLayout* pageLayout= new QHBoxLayout;
   aPage->setLayout(pageLayout);
   pagesWidget->addWidget(aPage);
   
   //add the button
   QRadioButton* aButton = new QRadioButton(elem_opt.hasAttribute("title")?elem_opt.attribute("title"):elem_opt.tagName(), this);
   if(elem_opt.hasAttribute("whatsThis"))aButton->setWhatsThis(elem_opt.attribute("whatsThis"));
   if(elem_opt.hasAttribute("toolTip"))aButton->setToolTip(elem_opt.attribute("toolTip"));
   if(elem_opt.hasAttribute("statusTip"))aButton->setStatusTip(elem_opt.attribute("statusTip"));
     

   
   buttonLayout->addWidget(aButton);
   buttonGroup->addButton(aButton);
   buttonGroup->setId(aButton, count);
   
   if(count==current)aButton->setChecked(true);
   values <<elem_opt.tagName();
   
   if(elem_opt.attribute("element")=="stateStack"){
     StackGroup2* stackGroup = new StackGroup2(elem_opt, myroot, pState, parent);
     stacks<<stackGroup;
     connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(updateState(QString)));
     connect(stackGroup, SIGNAL(componentsChanged()), this, SLOT(updateComponents()));
     pageLayout->addWidget(stackGroup);
   }else if(elem_opt.attribute("element")=="models"){
       QDomElement elem_mod= elem_opt.firstChildElement();
       QGroupBox* modelGroup = new QGroupBox;
       
       if(!elem_mod.isNull()){
     
       
      
       QHBoxLayout* modelsLayout = new QHBoxLayout;
       modelGroup->setLayout(modelsLayout);
       modelGroup->setFlat(true);
       for (; !elem_mod.isNull(); elem_mod = elem_mod.nextSiblingElement()){
         QDomElement the_model = elem_models.firstChildElement(elem_mod.tagName());
         ChoiceGroup* modelChoice = new ChoiceGroup(the_model, myroot);
         if(elem_mod.hasAttribute("condition")){
          
           conditionalModels << modelChoice;
           conditionalElems << elem_mod;
           if( !conditionIsSatisfied(myroot, elem_mod.attribute("condition"))){
             modelChoice->hide();
             elem_mod.setAttribute("conditionSatisfied", "false");
           }else{
             modelChoice->show();
             elem_mod.setAttribute("conditionSatisfied", "true");
           }
         }else{
           unconditionalModels << modelChoice;
         }
         
         modelsLayout->addWidget(modelChoice);
         connect(modelChoice, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
       }
       }
      
       pageLayout->addWidget(modelGroup);
   }
   
   if(elem_opt.hasAttribute("define")){
     
     QPushButton* addButton=new QPushButton("add a module");
     pageLayout->addWidget(addButton);
     connect(addButton, SIGNAL(clicked()), this, SLOT(add()));
   }
 }


 mainLayout->addWidget(pagesWidget);
 pagesWidget->setCurrentIndex(current);
 connect(buttonGroup,
         SIGNAL(buttonClicked(int)),
         this, SLOT(changePage(int)));
 
 setLayout(mainLayout);

 if(values.size()!=0)parentState<<name+"="+values[0];
 else parentState<<name;

 emit textChanged(currentText());
}
void StackGroup2::updateComponents(){
  
  int currentIndex= myelem.attribute("current").toInt();
  QDomElement elt = myelem.firstChildElement();
  for(int i = 0; i < currentIndex; i++)elt=elt.nextSiblingElement();
  if(elt.hasAttribute("components")){
    if(myelem.attribute("components")!=elt.attribute("components")){
      myelem.setAttribute("components", elt.attribute("components"));
      myroot.firstChildElement("components").setAttribute("components", elt.attribute("components"));
      emit componentsChanged();
    }
  }
 
}
QString StackGroup2::myState(){
  
  if(stacks.size()==0){
    if(values.size()==0) return name;
    return name+"="+values[current];
  }else{
    QString childState = stacks[current]->myState();
    return name+"="+childState.replace('=',',');
  }
}

void StackGroup2::changeState(){
}
void StackGroup2::updateCurrentText(){
  
  emit textChanged(currentText());
}

QString StackGroup2::currentText(){
  int currentIndex= myelem.attribute("current").toInt();
  QDomElement elt = myelem.firstChildElement();
  for(int i = 0; i < currentIndex; i++)elt=elt.nextSiblingElement(); 
  QString tmp="";
  if(myelem.hasAttribute("comments")) tmp += myelem.attribute("comments")+"\n";

  bool done = true;
  if(elt.attribute("element")=="models"){
    QDomElement elem_models= myroot.firstChildElement("models");
    QDomElement elem_mod= elt.firstChildElement();
    for (; !elem_mod.isNull(); elem_mod = elem_mod.nextSiblingElement()){
      QDomElement the_model = elem_models.firstChildElement(elem_mod.tagName());
      if(elem_mod.hasAttribute("condition")){
        if(elem_mod.attribute("conditionSatisfied")=="true"){
          tmp += the_model.attribute("currentText") + "\n";
          if(the_model.attribute("status")!="done")done = false;
          
          QDomElement cmpnt_elem = the_model.firstChildElement();
          for(int i = 0; i < the_model.attribute("current").toInt(); i++)cmpnt_elem=cmpnt_elem.nextSiblingElement();
          if(cmpnt_elem.hasAttribute("components")){
            if(cmpnt_elem.attribute("components") != myelem.attribute("components")){
              myelem.setAttribute("components", cmpnt_elem.attribute("components"));
              myroot.firstChildElement("components").setAttribute("components", elt.attribute("components"));
              emit componentsChanged();
            }
          }
        }   
      }else{
        tmp += the_model.attribute("currentText") + "\n";
        if(the_model.attribute("status")!="done")done = false;
        QDomElement cmpnt_elem = the_model.firstChildElement();
        for(int i = 0; i < the_model.attribute("current").toInt(); i++)cmpnt_elem=cmpnt_elem.nextSiblingElement();
        if(cmpnt_elem.hasAttribute("components")){
          if(cmpnt_elem.attribute("components") != myelem.attribute("components")){
            myelem.setAttribute("components", cmpnt_elem.attribute("components"));
             myroot.firstChildElement("components").setAttribute("components", elt.attribute("components"));
            emit componentsChanged();
          }
        }
      }
    }
  
  }else if(elt.attribute("element")=="stateStack"){
    if(elt.hasAttribute("components")){
      if(myelem.attribute("components")!=elt.attribute("components")){
        myelem.setAttribute("components", elt.attribute("components"));
        myroot.firstChildElement("components").setAttribute("components", elt.attribute("components"));
        emit componentsChanged();
      }
    }
    tmp += elt.attribute("currentText") +"\n";
    if(elt.attribute("status")!="done")done = false; 
  }
  myelem.setAttribute("currentText", tmp);
  if(done) myelem.setAttribute("status", "done");
  else myelem.removeAttribute("status");
  return tmp;
}
  



void StackGroup2:: parentStateChanged(QString stat){
  QStringList list1=stat.split("=");
   //update parentstate
  bool Found = false;
  for(int i = 0; i < parentState.size(); i++){
    if(parentState[i].contains(list1[0])){
      parentState[i] = stat;
      Found = true;
    }
  }
  if(!Found){
    QMessageBox::warning(window(), tr("physicswindow.cpp"),
                         tr("illegal state")+stat
                         );
    return;
  }
   //update conditional models  
  for(int i=0; i<conditionalElems.size(); i++){
    if(conditionalElems[i].attribute("condition").contains(list1[0])){
      if(conditionIsSatisfied(myroot,conditionalElems[i].attribute("condition"))){
        conditionalModels[i]->show();
        conditionalElems[i].setAttribute("conditionSatisfied","true");
      } else {
        conditionalModels[i]->hide();
        conditionalElems[i].setAttribute("conditionSatisfied", "false");
        
      }
    }
  }
  emit textChanged(currentText());
 }


     
void StackGroup2::changePage(int i)
{
  if(current==i) return;
  current = i;
  myelem.setAttribute("current", current);
  pagesWidget->setCurrentIndex(i);
  QString tmp = myState();
  emit stateChanged(tmp);
  emit textChanged(currentText());
}

void StackGroup2::updateState(QString stat)
{
  
  QString tmp = name+"="+stat.replace('=',',');
  emit stateChanged(tmp);
  emit textChanged(currentText());
}



void StackGroup2::setParentState(const QStringList& ps){
  parentState = ps;
}


GeneralWindow:: GeneralWindow(QDomElement &elem,
                              QDomElement &root,
                              QWidget *parent):QWidget(parent),
                                               myelem(elem),myroot(root){}


Page::Page(QDomElement& theelem, QDomElement& theroot, QWidget* parent): GeneralWindow(theelem, theroot, parent)
{
  int numColumn= 1;
  if(myelem.hasAttribute("numColumn")) numColumn = myelem.attribute("numColumn").toInt();
  QGridLayout *mainLayout = new QGridLayout;

  for(int i = 0; i < numColumn; i++){
    mainLayout->setColumnMinimumWidth(i, 30);
    mainLayout->setColumnStretch(i, 1);
  }  
  
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+tr(" has no child")
                         );
    return;
  }
  int elt_count = 0;
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), elt_count++) {   
    OptionPage* bdCndPage;
    if(elem.attribute("element")=="panel"){
      bdCndPage = new OptionPage(elem, myroot);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
    }else{
      QMessageBox::warning(window(), elem.tagName(),
                           tr(" don't know how to handle it yet: ")
                           +elem.attribute("element")
                           );
          
    }
    mainLayout->addWidget(bdCndPage, elt_count/numColumn, elt_count%numColumn, 1, 1);
  }//for(; elt..)
        

  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  updateCurrentText();
}

void GeneralWindow::changeState(){

  if(myelem.hasAttribute("condition")){
    bool isSatisfied = conditionIsSatisfied(myroot, myelem.attribute("condition"));
    if(isSatisfied){
      myelem.setAttribute("conditionSatisfied", "true");
      if(!myelem.hasAttribute("index")) this->show();
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
      if(!myelem.hasAttribute("index"))  this->hide();
    }
  }
  emit stateChanged();
}

void GeneralGroup::changeState(){
  if(myelem.hasAttribute("condition")){
    bool isSatisfied = conditionIsSatisfied(myroot, myelem.attribute("condition"));
    if(isSatisfied){
      myelem.setAttribute("conditionSatisfied", "true");
       if(!myelem.hasAttribute("index")) this->show();
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
      if(!myelem.hasAttribute("index"))  this->hide();
    }
  }
  emit stateChanged();
}


void Page::updateCurrentText(){
  QString text;
  if(myelem.hasAttribute("comments")) text += myelem.attribute("comments")+"\n";
  QString infix = "\n";
  int count = 0;
  int count_done = 0;
  for(QDomElement elt = myelem.firstChildElement(); !elt.isNull(); elt = elt.nextSiblingElement()){
    text +=  elt.attribute("currentText")+infix;
    if(elt.attribute("status")=="done")count_done++;
    count++;
  }
  myelem.setAttribute("currentText", text);
  emit textChanged(text);
  
  if(count==count_done) myelem.setAttribute("status", "done");
  else myelem.setAttribute("status", QString("%1  out of ").arg(count_done)+QString(" %1 finished").arg(count));
  if(myelem.hasAttribute("buttonIndex")) emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
}


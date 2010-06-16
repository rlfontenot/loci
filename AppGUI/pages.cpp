#include <QtGui>
#include <QDomNamedNodeMap>
#include "pages.h"
#include <vector>
#include "mdlwindow.h"
using std::vector;


#define PI 3.14159265358979323846264338327950

/*!
  copy the attributes and child elements of \a fromElement to \a toElement.
  For each attribute \a att of \a fromElement, if \a toElement has no attribute \a att, att will be copied to \a toElement.
  If \a toElement has no child element, all child elements of \a fromElement will be copied to \a toElement.
 */
bool copy_element(const QDomElement& fromElement, QDomElement& toElement){
  if(fromElement.isNull() || toElement.isNull()) return false;
  QDomNamedNodeMap fromatt = fromElement.attributes();
  for(int i = 0; i < fromatt.length(); i++){
    if(!toElement.hasAttribute(fromatt.item(i).nodeName())){
      toElement.setAttribute(fromatt.item(i).nodeName(), fromatt.item(i).nodeValue());
    }
  }
  if((!fromElement.firstChildElement().isNull())&&toElement.firstChildElement().isNull()){
    for(QDomElement elt = fromElement.firstChildElement();
        !elt.isNull(); elt = elt.nextSiblingElement()){
      QDomNode tmp= toElement.appendChild(elt.cloneNode().toElement());
      if(tmp.isNull())return false;
    }
  }
  return true; 
}
      
  





//order of rotation y first(heading), then z(atitude), then x(bank)
bool angleBetween(positions3d v1,positions3d v2, double& heading, double& attitude, double& bank) {
  double norm1 = norm(v1);
  double norm2 = norm(v2);
  double angle = 0;
  positions3d axis = positions3d(0, 0, 1);
  
  if(norm1 > 1e-34 && norm2 > 1e-34){
        
    // turn vectors into unit vectors 
    positions3d n1 = v1/norm1;
    positions3d n2 = v2/norm2; 
    angle = acos(dot(n1,n2)); 
	
    if( fabs(angle) > 1e-5 && fabs(angle-PI)>1e-5){
      axis = cross(n1, n2);
    }else if(fabs(angle) <= 1e-5){
      // if no noticable rotation is available return zero rotation
      // this way we avoid Cross product artifacts 
      axis = positions3d(0, 0, 1);
      angle = 0;
    }else if(fabs(angle-PI) < 1e-5){
      // in this case there are 2 lines on the same axis 
      // there are an infinite number of normals 
      // in this case. Anyone of these normals will be 
      // a valid rotation (180 degrees). so I pick one axis that is perpendicular to n1
      angle = PI;
      if(fabs(fabs(n1.z)-1.0)<1e-5)axis=positions3d(1, 0, 0);
      else axis = positions3d(-1*n1.y, n1.x, 0);
    }
  }else{
    return false;
  }
  
  //from axis angle to Euler
  double s=sin(angle);
  double c=cos(angle);
  double t=1-c;
  //  if axis is not already normalized then uncomment this
  double magnitude = norm(axis);
  if (magnitude==0){
    cerr << " maginitude is zero" << endl;
    return false;
  }
  axis = axis/magnitude;
  double x = axis.x;
  double y = axis.y;
  double z = axis.z;
  
  if ((x*y*t + z*s) > 0.998) { // north pole singularity detected
    heading = 2*atan2(x*sin(angle/2.0),cos(angle/2.0));
    attitude = PI/2;
    bank = 0;
    return true;
  }
  if ((x*y*t + z*s) < -0.998) { // south pole singularity detected
    heading = -2*atan2(x*sin(angle/2),cos(angle/2));
    attitude = -PI/2;
    bank = 0;
    return true;
  }
  heading = atan2(y * s- x * z * t , 1 - (y*y+ z*z ) * t);
  attitude = asin(x * y * t + z * s) ;
  bank = atan2(x * s - y * z * t , 1 - (x*x + z*z) * t);
  
  return true;
}

bool conditionIsSatisfied(const QDomElement& theroot, const QString& condition){
  
  QDomElement elem= theroot.firstChildElement("mainWindow");
  elem = elem.firstChildElement("physicsPage");
  if(!elem.hasAttribute("state")) {
    return false;
  }
  if(condition.isEmpty()) {
    return false;
  }
  
  QStringList state = elem.attribute("state").split(";");
  QStringList list1=condition.split("==");
  if(list1.size()!=2){
    return false;
  }
  
 
  for(int i = 0; i < state.size(); i++){
    if(state[i].contains(list1[0]) && state[i].contains(list1[1])){
      return true;
    }
  }
  return false;
}


/*!
  \class DoubleEdit
  
  \brief The DoubleEdit a line editor that takes doubles .
  
   A DoubleEdit allows the user to enter a double value. The range of the
  value can be set up using setBottom(), setTop() or setRange() methods.
  The value of the number is given by value().

  The value of the double can be set using setValue().

  The valueChanged() signal is emitted by slot changeValue() when the text in the line editor changes.

    
*/




/*!
    Constructs a new DoubleEdit with the given \a parent.
*/
DoubleEdit::DoubleEdit(QWidget *parent) : QLineEdit(parent){

  validator = new QDoubleValidator(this);
  setValidator(validator);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
  setMinimumWidth(100);
}

/*!
    Constructs a new DoubleEdit with the given \a parent and a double \a d.
*/

DoubleEdit::DoubleEdit(double d, QWidget*parent):QLineEdit(parent){
  double d0 = d;
  validator = new QDoubleValidator(this);
  setValidator(validator);
  QString str;
  str.setNum(d0);
  setText(str);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
  emit valueChanged(d0);
}

void DoubleEdit::setValue(double d){
  double d0 = d;
  QString str;
  str.setNum(d0);
  setText(str);
  emit valueChanged(d0);
}

/*!
  This function is used when the DoubleEdit is connected with an interger editor, such as slider or spin box. Whenever the
  integer value changes, the double value of the DoubleEdit is set by mapping the range of the interger editor  to the range of the DoubleEdit. 

  The range of the interger editor must be set to [-1000, 1000].
  
  For example, if \a value is  500, the range of double is[0, 1.0], then the value shown on DoubleEdit  is (500+1000.)/2000.*(1.0-0), which is 0.75
  
*/
void DoubleEdit::mapValue(int value){
  double d = validator->bottom() +(value+1000)/2000.0*(validator->top()-validator->bottom());
  QString str;
  str.setNum(d);
  setText(str);
  emit valueChanged(d);
}

/*!
  Whenever the textChanged() signal is emitted, also emit the valueChanged() signal
*/
void DoubleEdit::changeValue(const QString& str){
  double d = str.toDouble();
  emit valueChanged(d);
}
  
double DoubleEdit::value(){
  QString str;
  str = text();
  return str.toDouble();
}

void DoubleEdit::setBottom(double d){
  double d0 = d;
  validator->setBottom(d0);
}

void DoubleEdit::setTop(double d){
  double d0 = d;
  validator->setTop(d0);
}

void DoubleEdit::setRange(double d1, double d2){
  double dd1 = d1;
  double dd2 = d2;
  validator->setRange(dd1, dd2);
}

/*!
  \class DoubleSpinBox
  
  \brief class DoubleSpinBox is a QDoubleSpinBox that allow the user to change the value
  of step and decimal using keyPressEvent.
  
   A DoubleSpinBox reimplement keyPressEvent() so that step and decimal values can be adjusted from the interface.

  key "S": increase step value  10 times
  
  key "s": decrease step value 0.1 times

  key "D": increase decimal value by 1

  key "d": decrease decimal value by 1
    
*/

DoubleSpinBox::DoubleSpinBox(QWidget* parent):QDoubleSpinBox(parent){};
void DoubleSpinBox::keyPressEvent(QKeyEvent *event){
  QString key = event->text();
  
  if(key=="S"){
    double step = singleStep()*10;
    setSingleStep(step);
    emit paraChanged();
  }else if(key== "s"){
    double step = singleStep()*0.1;
    setSingleStep(step);
    emit paraChanged();
      
  }else if(key =="D"){
    int decimal = decimals()+1;
    setDecimals(decimal);
    emit paraChanged();
           
  }else if(key=="d"){
     
    int decimal = decimals()-1;
    if(decimal>=0){
      setDecimals(decimal);
      emit paraChanged();
    }
    
  }else{
    QAbstractSpinBox::keyPressEvent(event);
  }
}

/*!
  \class LabeledDoubleSpBox
  
  \brief class LabeledDoubleSpBox is a DoubleSpinBox with a label in front of it.
  Usually the label displays the name of the variable.
  
   Except for the name label to the right of spin box, LabbeledDoubleSpBox has another hidden label to the left of it. when the key "S", "s", "D" or "d" is pressed, this label will show the step value and decimal value. When the vale of the spin box changed, this label will hide.  
  
    
*/

LabeledDoubleSpBox::LabeledDoubleSpBox(const QString& title,  QWidget *parent):QWidget(parent)
{
 
  edit = new DoubleSpinBox;
  edit->setSingleStep(0.01);
  edit->setRange(-1e5, 1e5);
  edit->setDecimals(8);
  connect(edit, SIGNAL(valueChanged(double)), this, SIGNAL(valueChanged(double)));
  connect(edit, SIGNAL(paraChanged()), this, SLOT(display()));
  connect(edit,  SIGNAL(valueChanged(double)), this, SLOT(undisplay()));
  
  QLabel* label = new QLabel(title);
  valueLabel = new QLabel;
  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->addWidget(label);
  mainLayout->addWidget(edit);
  mainLayout->addWidget(valueLabel);
  valueLabel->hide();
  setLayout(mainLayout);
}

void LabeledDoubleSpBox::display(){
  valueLabel->setText(QString("Decimals: %1   Single Step: %2")
                      .arg(edit->decimals())
                      .arg(edit->singleStep()));
  emit paraChanged();
  valueLabel->show();
}
void LabeledDoubleSpBox::undisplay(){
  valueLabel->hide();
}

void LabeledDoubleSpBox::setValue(double d){
  edit->setValue(d);
}
double LabeledDoubleSpBox::value(){
  return edit->value();
}
void LabeledDoubleSpBox::setRange(double d1, double d2){
  edit->setRange(d1, d2);
}
/*!
  \class VectSpBox
  
  \brief class VectSpBox is a group box that allows the user to
  input the x, y, z values of a vector using spin boxes..
  
   The range of the spin box can be set universely using setRange()
  or separately using setXRange(), setYRange() or setZRange().
  
  Since x, y, z components usually share the same precision, whenever the single
  step value or decimals value of x component changed, the values of y and z
  components will be set as the same.   
  
*/

VectSpBox::VectSpBox( const QString& title, QWidget *parent):QGroupBox(title, parent){
  
  xedit = new LabeledDoubleSpBox("x");
  yedit = new LabeledDoubleSpBox("y");
  zedit = new LabeledDoubleSpBox("z");
  connect(xedit, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(yedit, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(zedit, SIGNAL(valueChanged(double)), this, SLOT(setInfo()));
  connect(xedit, SIGNAL(paraChanged()), this, SLOT(setPara()));
  
  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->addWidget(xedit);
  mainLayout->addSpacing(10);
  mainLayout->addWidget(yedit);
  mainLayout->addSpacing(10);
  mainLayout->addWidget(zedit);
  setLayout(mainLayout);
}

void VectSpBox::setPara(){
  yedit->setSingleStep(xedit->singleStep());
  zedit->setSingleStep(xedit->singleStep());
  yedit->setDecimals(xedit->decimals());
  zedit->setDecimals(xedit->decimals());
  
}
  
void VectSpBox::setValue(const positions3d& p){
  xedit->setValue(p.x);
  yedit->setValue(p.y);
  zedit->setValue(p.z);
}
void VectSpBox::setInfo(){
  emit valueChanged(value());
}
positions3d VectSpBox::value(){
  return positions3d(xedit->value(), yedit->value(), zedit->value());
}

void VectSpBox::setRange(double d1, double d2){
  xedit->setRange(d1, d2);
  yedit->setRange(d1, d2);
  zedit->setRange(d1, d2);
}
void VectSpBox::setXRange(double d1, double d2){
  xedit->setRange(d1, d2);
}

void VectSpBox::setYRange(double d1, double d2){
  yedit->setRange(d1, d2);
}
void VectSpBox::setZRange(double d1, double d2){
  zedit->setRange(d1, d2);
}


/*!
  \class IntEdit
  
  \brief The IntEdit a line editor that takes doubles .
  
   A IntEdit allows the user to enter a integer value. The range of the
  value can be set up using setBottom(), setTop() or setRange() methods.
  The value of the number is given by value().

  The value of the integer can be set using setValue().
 

  The valueChanged() signal is emitted by slot changeValue() when the text in line editor changes.

    
*/


IntEdit::IntEdit(QWidget *parent) : QLineEdit(parent){
  validator = new QIntValidator(this);
  setValidator(validator);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
   
}

IntEdit::IntEdit(int d, QWidget*parent):QLineEdit(parent){
  int d0 = d;
  validator = new QIntValidator(this);
  setValidator(validator);
  QString str;
  str.setNum(d0);
  setText(str);
  connect(this, SIGNAL(textChanged(const QString&)), this, SLOT(changeValue(const QString&)));
  emit valueChanged(d0);
}
void IntEdit::setValue(int d){
  int d0 = d;
  QString str;
  str.setNum(d0);
  setText(str);
  emit valueChanged(d0);
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
  int d0 = d;
  validator->setBottom(d0);
}

void IntEdit::setTop(int d){
  int d0 = d;
  validator->setTop(d0);
}

void IntEdit::setRange(int d1, int d2){
  int dd1 = d1;
  int dd2 = d2;
  validator->setRange(dd1, dd2);
}


/*!
  \class GeneralGroup
  
  \brief The GeneralGroup is a base class of most chemdemo objects. It's a QGroupBox.
  
   GeneralGroup is a group box that use a QDomElement \a myelem to define its attributes
  and child widgets, and all the user input from
  GeneralGroup will update the attributes and child elements of \a myelem. 
    
*/

/*!
    Constructs a GeneralGroup with the given \a my_elem and \a parent.
    The attributes "title", "checked", "whatsThis", "toolTip" and "statusTip" of
    \a my_elem are used to define the attributes of the group box. 
    The \a parent argument is passed on to the QGroupBox constructor.
    The \a my_elem argument is assigned to member data \a myelem.

    The signal toggled() of groupbox is connected to updateChecked() slot of this.
*/
GeneralGroup::GeneralGroup( QDomElement& my_elem, QWidget *parent ) : QGroupBox(my_elem.hasAttribute("title")?my_elem.attribute("title"):my_elem.tagName(), parent),myelem(my_elem){
  
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

/*!
  When stateChanged() signal from other GeneralGroup(usually from parant or mainwindow) is received, if \a myelem
  has attribute "condition", check if the condition is satisfied, and update the
  attribute "conditionSatisfied" of \a myelem.
  
  If \a myelem has no attribute 'index', means it's not in a stacked widget or tabbed widget, then
  this  will hide or show according to if the conditions is satisfied.

  If \a myelem does have attribute 'index', then its parent should decide what to do,for example, it can be disabled/enabled in tabbed widget.

 Signal stateChanged() is emitted to notify all child widgets. 
*/
void GeneralGroup::changeState(){
  if(myelem.hasAttribute("condition")){
    QDomElement myroot = myelem.ownerDocument().documentElement();
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

/*!
  This method will be reimplemented by subclasses to return the current text.
*/
QString GeneralGroup::currentText(){
  return QString();
}

/*!
  When the groupbox is toggled, this slot will update the attribute "checked"
  of \a myelem and emit signal textChanged().
*/
void GeneralGroup::updateChecked(){
  if(isChecked())myelem.setAttribute("checked", 1);
  else myelem.setAttribute("checked", 0);
  emit textChanged(currentText());
  
}

/*!
  When GeneralGroup receives the signal showStatus()(usually from \a parent or mainwindow),
  this slot will  decide whether the child widgets need to show their status or not, and then
  emit signal showStatus() to notify the child widgets.

  This slot will be reimplemented by some subclasses to actually show the status.
  
*/
void GeneralGroup::updateShowStatus(const bool& show){

  if(show){//parent says show
    if(((!myelem.hasAttribute("checked")) || myelem.attribute("checked").toInt()==1) &&
       myelem.attribute("status")!="done" )emit showStatus(show); //tell my children to show
    else  emit showStatus(!show); //tell my children don't show
  }else  emit showStatus(show); //parent no show, my children won't show
  
}

/*!
  \class VarGBox
  
  \brief The VarGBox is a basic unit of chemdemo GUI. It allows the user to input a variable.
 
  
   The element of variable can be: vector, int, float, dvector, selection, string.

  The variable has a value, and optionally has a unit.

  VarGBox allows the user to enter the value of the variable, and select the unit from a combo box,
  it also allows the user to edit the unit list from the combo box.

  The name of the variable that appears in .var file and label appears on
  interface can be different, and both can be specified in xml file.

  For a vector, VarGBox allows the user to select Cartesian, Polar or Scale format.
  Cartesian format will take x, y, z components of the vector, Polar format allows the user to enter
  magnitude, theta and phi components, Scale format  needs x component only.
  
  dvector is used for variable \a components in multiple components model. Each component has a name and a double value.
  
  selection is for flags, it's a group of checkable widgets, the item checked will appear in the output.

  
  
  \sa \ref xml_variable 

*/




VarGBox::VarGBox(  QDomElement& my_elem, QWidget *parent ) : GeneralGroup(my_elem, parent){
  QDomElement myroot = myelem.ownerDocument().documentElement();
  

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
    QPointer<DoubleEdit> xEditor = new DoubleEdit;
    QPointer<DoubleEdit> yEditor = new DoubleEdit;
    QPointer<DoubleEdit> zEditor = new DoubleEdit;
    
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
        QPointer<DoubleEdit> aEditor = new DoubleEdit;
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
          if(i < components.size()) qobject_cast<DoubleEdit*>(mfs[i])->setValue(valueList[i].toDouble());
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
      QPointer<DoubleEdit> valueEdit = new DoubleEdit;
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
                           tr("variables: don't know how to handle it yet: ") + myelem.attribute("element")
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
  connect(this, SIGNAL(showStatus(const bool &)), this, SLOT(updateShowStatus(const bool &)));
}

/*!
  In case  the element is 'selection', when one checkable widget is toggled, update attribute
  "checked" in the xml file and emit signal textChanged()
*/
void VarGBox::updateSelection(int id){
  QDomElement opt = myelem.firstChildElement();
  for(int i=0; i < id; i++)opt = opt.nextSiblingElement();
  if(opt.attribute("checked").toInt()== 1)opt.setAttribute("checked", 0);
  else opt.setAttribute("checked", 1);
  emit textChanged(currentText());
}

/*!
  In case  the element is 'string', 'float', 'int' or 'dvector', when the value of the variable is entered,
  update the attribute "current" in  xml file and emit signal textChanged()
*/
void VarGBox::updateCurrent(const QString& c){
  if(myelem.attribute("element")=="string"
     ||myelem.attribute("element")=="float"
     ||myelem.attribute("element")=="int")myelem.setAttribute("current", c);
  else if(myelem.attribute("element")=="dvector"){
    QString tmp;
    for(int i = 0; i < mfs.size(); i++){
      tmp += qobject_cast<DoubleEdit*>(mfs[i])->text();
      if(i != mfs.size()-1)tmp +=",";
    }
    myelem.setAttribute("current", tmp);
  }
 
  emit textChanged(currentText());
}

/*!
  In signal showStatus() is received, this slot will change the background
  color of the unfinished field. 
*/
void VarGBox::updateShowStatus(const bool& show){
  if(show) myelem.setAttribute("showStatus", "1");
  else myelem.setAttribute("showStatus", "0");

    
  if(myelem.attribute("element")=="string"
     ||myelem.attribute("element")=="float"
     ||myelem.attribute("element")=="vector"
     ||myelem.attribute("element")=="dvector"
     ||myelem.attribute("element")=="int"
     ){

    if(!show){//not show
      for(int i = 0; i < mfs.size(); i++){
        
        QPalette palette;
        palette.setColor(mfs[i]->backgroundRole(), QColor(255, 255, 255));
        mfs[i]->setPalette(palette);
      }
    }else{//show
      if(myelem.attribute("status")=="done"){
        for(int i = 0; i < mfs.size(); i++){
          QPalette palette;
          palette.setColor(mfs[i]->backgroundRole(), QColor(255, 255, 255));
          mfs[i]->setPalette(palette);
        }
      }else  if(!myelem.hasAttribute("checked")||myelem.attribute("checked").toInt()== 1){
        if(myelem.attribute("element")=="string"
           ||myelem.attribute("element")=="float"
           ||myelem.attribute("element")=="vector"
           ||myelem.attribute("element")=="dvector"){
          
          for(int i = 0; i < mfs.size(); i++){
            if(qobject_cast<QLineEdit*>(mfs[i])->text()==""){
              QPalette palette;
              palette.setColor(mfs[i]->backgroundRole(), QColor(255, 0, 0));
              mfs[i]->setPalette(palette);
             
            }
          }
        }else if(myelem.attribute("element")=="int"){
          for(int i = 0; i < mfs.size(); i++){
            if(qobject_cast<QSpinBox*>(mfs[i])->cleanText()==""){
              QPalette palette;
              palette.setColor(mfs[i]->backgroundRole(), QColor(255, 0, 0));
              mfs[i]->setPalette(palette);
            }
          }
        }
      }else if(myelem.parentNode().toElement().attribute("element")=="2of3" && myelem.parentNode().toElement().attribute("status")!="done"){
        for(int i = 0; i < mfs.size(); i++){
          if(qobject_cast<QLineEdit*>(mfs[i])->text()==""){
            QPalette palette;
            palette.setColor(mfs[i]->backgroundRole(), QColor(255, 0, 0));
            mfs[i]->setPalette(palette);
             
          }
        } 
      }
    }
  }

}



/*!
  In case the element is 'vector', when the first component value is entered,
  update the attribute 'currentX' in xml file and emit textChanged()
*/ 
void VarGBox::updateCurrentX(const QString& c){
  myelem.setAttribute("currentX", c);
  emit textChanged(currentText());
}


/*!
  In case the element is 'vector', when the second component value is entered,
  update the attribute 'currentY' in xml file and emit textChanged()
*/ 
void VarGBox::updateCurrentY(const QString& c){
  myelem.setAttribute("currentY", c);
  emit textChanged(currentText());
}

/*!
  In case the element is 'vector', when the first component value is entered,
  update the attribute 'currentX' in xml file and emit textChanged()
*/ 
void VarGBox::updateCurrentZ(const QString& c){
  myelem.setAttribute("currentZ", c);
  emit textChanged(currentText());
}


/*!
  When the current index of  unit comboBox changed , 
  update the attribute 'currentUnit' in xml file and emit textChanged()
*/ 
void VarGBox::updateCurrentUnit(const QString& c){
  myelem.setAttribute("currentUnit", c);
  emit textChanged(currentText());
}

/*!
  When the edit text  of  unit comboBox changed , 
  update the attribute 'currentUnit' in xml file to \arg c and add \arg
  c to attribute 'unit'.
  
  emit signal textChanged().
*/ 
void VarGBox::updateUnitList(const QString& c){
  myelem.setAttribute("unit", myelem.attribute("unit")+","+c);
  myelem.setAttribute("currentUnit", c);
  emit textChanged(currentText());
}

/*!
  return the current text as appears in .var file
*/ 
QString  VarGBox::currentText(){
  
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
      text += labels[i]->text() + tr("= %1").arg(qobject_cast<DoubleEdit*>(mfs[i])->value());
      if(i != labels.size() -1) text+=", ";
      else text += "]";
      total += qobject_cast<DoubleEdit*>(mfs[i])->value();
    }
    if(total > 1.00000000000001){
      QMessageBox::warning(this, "VarGBox", tr("the sum is greater than 1!"));
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
  
 
  if(myelem.attribute("showStatus").toInt()==1)updateShowStatus(true);
    
  myelem.setAttribute("currentText", text);
  return text;
}

/*!
  If this variable is conditional, check if condition or defaultCondition is satisfied, update the attribute
  "conditionSatisfied" in xml file.
  
  if condition is satisfied, show this, otherwise hide this.
  
  if it's defaultCondition, call setDefault(). 
*/ 
void VarGBox::changeState(){
  QDomElement myroot = myelem.ownerDocument().documentElement();
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
    }else{
      myelem.setAttribute("conditionSatisfied", "false");
      this->setDefault(false);
    }
  }
}


/*!
 
  When physics state changed, if defaultCondition is satisfied,
  set values to default as described in xml file, and deactivate the user input field.
  Otherwise, activate the user input field.
  
  emit signal textChanged()
*/ 
void VarGBox::setDefault(bool satisfied){
  if(myelem.attribute("element")=="float"){
    if(mfs.size() >= 1){
      if(satisfied)qobject_cast<DoubleEdit*>(mfs[0])->setValue(myelem.attribute("default").toDouble());
     
      qobject_cast<DoubleEdit*>(mfs[0])->setEnabled(!satisfied);
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
        if(satisfied) qobject_cast<DoubleEdit*>(mfs[i])->setText(defaultList[i]);
        qobject_cast<DoubleEdit*>(mfs[i])->setEnabled(!satisfied); 
      }
    }
  }
  emit textChanged(currentText());
}



/*!
  When the components changed, if the element of the variable is 'dvector',
  remove the layout and all child widgets of VarGBox, and rebuild the interface.
 
*/ 
void VarGBox::updateComponents(){
 
  if(myelem.attribute("element")!="dvector")return;

  QDomElement myroot = myelem.ownerDocument().documentElement();
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
      QPointer<DoubleEdit> aEditor = new DoubleEdit;
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
  
  const QString str;
  updateCurrent(str);
  emit textChanged(currentText());
}




  
/*!
 update the label of components if the user change format for a vector.
 
*/ 
void VarGBox::updateLabels(int i){
  switch(i){
  case 0:
    labels[0]->setText("x");
    labels[1]->setText("y");
    labels[2]->setText("z");

    myelem.setAttribute("format", "cartesian");
    for(int i =1; i < 3; i++){
      qobject_cast<DoubleEdit*>(mfs[i])->setEnabled(true); 
    }
    break;
  case 1:
    labels[0]->setText("mag");
    labels[1]->setText("theta");
    labels[2]->setText("phi");

    
    myelem.setAttribute("format", "polar");
    for(int i =1; i < 3; i++){
      qobject_cast<DoubleEdit*>(mfs[i])->setEnabled(true); 
    }
    break;
  case 2:
    labels[0]->setText("x");
    labels[1]->setText("y");
    labels[2]->setText("z");

    

    
    myelem.setAttribute("format", "scale");
    for(int i =1; i < 3; i++){
      qobject_cast<DoubleEdit*>(mfs[i])->setValue(0.0); 
      qobject_cast<DoubleEdit*>(mfs[i])->setEnabled(false); 
    }
    break;
    
  default:
    QMessageBox::warning(window(), "VarGBox", "i out of range");
  }
  emit textChanged(currentText());
}

/*!
  \class AllGroup
  
  \brief AllGroup  allows the user to input a set of related variables. The set of variables
  is describles in xml file as one Dom element under parent element "groups", and each varaiable
  is defined as one child dom element of the group element. 
  
  
    AllGroup is used for a group whose element is "all" or "2of3". Its child VarGBoxes are layed out either horizontally or vertically.

  If the element is "all", each variable in the group need to be specified.

  For AllGroup whose element is "2of3", two and only two of three variables will be output in .var file.
  The pressure, density and temperature is in such a group. 
  
  
*/          



AllGroup::AllGroup(  QDomElement& elem,  bool isVertical, QWidget *parent )
  : GeneralGroup(elem, parent){
  
  QDomElement myroot = myelem.ownerDocument().documentElement();

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
    
    //if element is specified, take it, otherwise go to variables to find and replace it
    if(!elem_opt.hasAttribute("element")){
      QString title = elem_opt.attribute("title");
      QString text = elem_opt.attribute("name");
      
      
      QDomElement option_elem = myroot.firstChildElement("variables");
      if(option_elem.isNull()){
        QMessageBox::warning(window(), tr(" .xml"),
                             tr("can not find element 'variables'")
                             );
        return;
      }
      QDomElement newChild = option_elem.firstChildElement(elem_opt.tagName()).cloneNode().toElement();
      if(newChild.isNull()){
        QMessageBox::warning(window(), tr(" .xml"),
                             tr("'variables' has no child ")+elem_opt.tagName()
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
      bool  conditionSatisfied = conditionIsSatisfied(myroot, tmp);
      if(conditionSatisfied){
        QString tmp = elem_opt.attribute("default");
        elem_opt.setAttribute("current", tmp);
      }
    }
    
    
    VarGBox *myGroup = new VarGBox(elem_opt,this);
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
    connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
    myGroup->setFlat(true);
    mainLayout->addWidget(myGroup);
    count++;
    
    if(elem_opt.hasAttribute("condition")){
      QString tmp = myelem.attribute("condition");
      bool conditionSatisfied = conditionIsSatisfied(myroot,tmp);
      if(conditionSatisfied){
        myGroup->show();
        elem_opt.setAttribute("checked", 1);
      }else{
        myGroup->hide();
        elem_opt.setAttribute("checked", 0);
      }
    }
    if(elem_opt.hasAttribute("defaultCondition") ){
      bool conditionSatisfied = conditionIsSatisfied(myroot, elem_opt.attribute("defaultCondition"));
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

/*!
  For a group whose element is "2of3", if  child widget \a i is clicked,
  check if the other two child widgets. If both of them are checked, child widget \i
  will be unchecked.
 
*/ 
void AllGroup::childrenClicked(int i){
  
  if(myelem.attribute("element")=="2of3"){
    bool isChecked =  qobject_cast<VarGBox*>(objs[i])->isChecked();
    if(isChecked){
      if( ( qobject_cast<VarGBox*>(objs[(i+1)%3])->isChecked())
          && ( qobject_cast<VarGBox*>(objs[(i+2)%3])->isChecked())){
        qobject_cast<VarGBox*>(objs[i])->setChecked(false);
      }
    }
  }
  
  
  emit textChanged(currentText());
}


/*!
  compose and return the current text, at the same time check the status of the group,
  update the attribute "currentText" and "status" in xml file. 
 
*/ 
QString AllGroup::currentText(){
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
/*!
 emit signal textChanged().
 
*/ 
void AllGroup::updateCurrentText(){
  
  emit textChanged(currentText());
}


/*!
  \class StackGroup
  
  \brief StackGroup  allows the user to input a group of related variables/groups that  
  only one of them need to be specified. On the
  interface, There is an exclusive radio button group, if the user checks one
  of them, the interface of the selected variable/group will be brought up. 

  
 
 
*/          
  

StackGroup::StackGroup( QDomElement& elem, QWidget *parent )
  : GeneralGroup(elem, parent){
  QDomElement myroot = myelem.ownerDocument().documentElement();
  
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
  
    QGroupBox *myGroup;
    //if elem_opt is a group
    if(rx.exactMatch(elem_opt.tagName())){
      //if the elem_opt is  a group and it needs copy
      if(elem_opt.firstChildElement().isNull()){
        copy_element(elem_grp.firstChildElement(elem_opt.tagName()), elem_opt);
      }//finish copy group
     
      if(elem_opt.attribute("element")=="all"||elem_opt.attribute("element")=="selection"||elem_opt.attribute("element")=="2of3"){
        myGroup = new AllGroup(elem_opt, false,this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
      }else{
        qDebug()<<"don't know how to handle it yet";
      }
    }else{//if it's not a group
      //if the variables needs copy
      if(!elem_opt.hasAttribute("element")){
        copy_element(myroot.firstChildElement("variables").firstChildElement(elem_opt.tagName()), elem_opt);
      }//finish copy option
      myGroup = new VarGBox(elem_opt, this);
    
      connect(this, SIGNAL(componentsChanged()), myGroup, SLOT(updateComponents()));
      connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
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


VarPanel::VarPanel(   QDomElement& my_elem,  QWidget *parent )
  :GeneralGroup(my_elem, parent)
{
  QDomElement myroot = myelem.ownerDocument().documentElement();
  if(!toolTip().isEmpty())setToolTip("");
  if(!statusTip().isEmpty())setStatusTip("");
  if(!whatsThis().isEmpty())setWhatsThis("");
  showStatus = false;

  buttonGroup = 0;  
  QDomElement elt = myelem.firstChildElement();
  if(elt.isNull()){
    QVBoxLayout *mainLayout = new QVBoxLayout;
    QLabel* aLabel = new QLabel("No Options");
    aLabel->setAlignment(Qt::AlignCenter);
    mainLayout->addWidget(aLabel);
    mainLayout->addStretch(10);
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
        copy_element(elem_grp.firstChildElement(elt.tagName()),elt);        
      }
      if(elt.attribute("element")=="all" ||elt.attribute("element")=="2of3"  ){
        myGroup = new AllGroup(elt, false, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="stack"){
        myGroup = new StackGroup(elt, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="choice"){
        myGroup = new ChoiceGroup(elt, this);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
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
    }else {//variables, panels
      
      if(!elt.hasAttribute("element")){//need copy
        copy_element(myroot.firstChildElement("variables").firstChildElement(elt.tagName()), elt);
      }
      if(elt.attribute("element")=="choice"){
        myGroup = new ChoiceGroup(elt);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="all"){
        myGroup = new AllGroup(elt);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="stack"){
        myGroup = new StackGroup(elt);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else if(elt.attribute("element")=="panel"){
       
        myGroup = new VarPanel(elt);
        connect(this, SIGNAL(componentsChanged()), myGroup, SIGNAL(componentsChanged()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
        connect(this, SIGNAL(stateChanged()), myGroup, SLOT(changeState()));
      }else{
        myGroup = new VarGBox(elt);
        connect(this, SIGNAL(componentsChanged()), myGroup, SLOT(updateComponents()));
        connect(this, SIGNAL(showStatus(const bool &)), myGroup, SLOT(updateShowStatus(const bool &)));
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
  }else{
    mainLayout->addWidget(new QLabel, (elt_count/numColumn +2), 0, 1,  numColumn);
  }
  mainLayout->setRowStretch( elt_count/numColumn +2, 10);

  setLayout(mainLayout);
  
 
 
  updateCurrentText();

   
}
  
void VarPanel::updateShowStatus(const bool& show){
  showStatus = show;
  if(buttonGroup){
    for( int i = 0; i <  buttonGroup->buttons().size(); i++){
      buttonGroup->button(i)->setPalette(this->palette());
    }
    
    if(show){//parent says show
    
      if(((!myelem.hasAttribute("checked")) || myelem.attribute("checked").toInt()==1) &&
         myelem.attribute("status")!="done" ){
        //emit showStatus(show); //tell my children to show
      
        //if advanced not done, bring up the page
     
        int count = 0;
        for(QDomElement elt = myelem.firstChildElement(); !elt.isNull(); elt = elt.nextSiblingElement()){
          if(elt.hasAttribute("buttonTitle")){
            if(elt.attribute("status")!= "done"){
              if(count < buttonGroup->buttons().size()){
                QPalette palette;
                palette.setColor(buttonGroup->button(count)->backgroundRole(), QColor(255, 0, 0));
                buttonGroup->button(count)->setPalette(palette);
              }
            }
        
            count++;
          
          }
      
        
        }
      }
    }
  }
  GeneralGroup::updateShowStatus(show);   
     
        
}


void VarPanel::clearText(){
}


void VarPanel::updateCurrentText(){
  QString tmp = currentText();
  myelem.setAttribute("currentText", tmp);
  emit textChanged(tmp);
  if(showStatus) updateShowStatus(showStatus);
}

  
QString VarPanel::currentText(){
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
  
void VarPanel::advancedButtonClicked(int index){
  if(index<myAdvancedGroup.size())  {
    myAdvancedGroup[index]->show();
  }
}

VarPanel::~VarPanel(){
  for(int i= 0; i < myAdvancedGroup.size(); i++){
    myAdvancedGroup[i]->close();
    if(myAdvancedGroup[i]){
      delete myAdvancedGroup[i];
      myAdvancedGroup[i] = 0;
       
    }
  }
  if(buttonGroup){
    delete buttonGroup;
    buttonGroup = 0;
  }
  destroy();
}







/*!
  \class ChoiceGroup
  
  \brief ChoiceGroup  actually allows the user to define the value of one variable by selecting one
  from a set of variables. However, once the choice is made, there might be other variables need to be specified with it.
  The dom element of ChoiceGroup is no
  
   

  
 
 
*/          
  

ChoiceGroup::ChoiceGroup(   QDomElement& my_elem, QWidget *parent )
  : GeneralGroup(my_elem, parent)
{
  
  QDomElement myroot = myelem.ownerDocument().documentElement();
  showStatus = false;
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
void ChoiceGroup::updateShowStatus(const bool& show){
  showStatus = show;
  if(editButton){
    if(show){//parent says show
   
      if(((!myelem.hasAttribute("checked")) || myelem.attribute("checked").toInt()==1) &&
         myelem.attribute("status")!="done" ){
      
        QPalette palette;
        palette.setColor(editButton->backgroundRole(), QColor(255, 0, 0));
        editButton->setPalette(palette);
   
        // editButton->click();       
      }else{
        // emit showStatus(!show); //tell my children don't show
      
     
        editButton->setPalette(this->palette());
    
      }
    }else{
      //emit showStatus(show); //parent no show, my children won't show
   
      editButton->setPalette(this->palette());
    
    }
  }
  GeneralGroup::updateShowStatus(show);
  
}

QString ChoiceGroup::currentText(){
  QDomElement myroot = myelem.ownerDocument().documentElement();
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
        editGroup = new  VarPanel(elt_item);
      }else{
        editGroup = new VarGBox(elt_item);
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
  if(showStatus)updateShowStatus(showStatus);
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
  QDomElement myroot = myelem.ownerDocument().documentElement();
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
    editGroup = new  VarPanel(elt);
    connect(this, SIGNAL(componentsChanged()), editGroup, SIGNAL(componentsChanged()));
    connect(this, SIGNAL(showStatus(const bool &)), editGroup, SLOT(updateShowStatus(const bool &)));
    connect(this, SIGNAL(stateChanged()), editGroup, SLOT(changeState()));
  }else{
    editGroup = new VarGBox(elt);
    connect(this, SIGNAL(componentsChanged()), editGroup, SLOT(updateComponents()));
    connect(this, SIGNAL(showStatus(const bool &)), editGroup, SLOT(updateShowStatus(const bool &)));
    connect(this, SIGNAL(stateChanged()), editGroup, SLOT(changeState()));
  }
  connect(editGroup, SIGNAL(textChanged(QString)), this, SLOT(updateCurrentText()));
  editGroup->show();
  updateCurrentText();
 
}
AllVBWindow::AllVBWindow(  QDomElement& elem, QWidget *parent):GeneralGroup(elem, parent){
  if(myelem.hasAttribute("whatsThis"))setWhatsThis(myelem.attribute("whatsThis"));
  if(myelem.hasAttribute("toolTip"))setToolTip(myelem.attribute("toolTip"));
  if(myelem.hasAttribute("statusTip"))setStatusTip(myelem.attribute("statusTip"));
  


  
  AllGroup* thegroup = new AllGroup(elem, this);
  connect(this, SIGNAL(componentsChanged()), thegroup, SIGNAL(componentsChanged()));
  connect(this, SIGNAL(showStatus(const bool &)), thegroup, SLOT(updateShowStatus(const bool &))); 
  connect(this, SIGNAL(stateChanged()), thegroup, SLOT(changeState()));
  
  QPushButton* saveButton = new QPushButton(tr("          save the model        "));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  
  
  QVBoxLayout* mainLayout = new QVBoxLayout;

  mainLayout->addWidget(thegroup);
  mainLayout->addWidget(saveButton);
  setLayout(mainLayout);
}


bool AllVBWindow::save(){
  QString initialPath = QDir::currentPath()+"/untitled.mdl";
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                  initialPath,
                                                  tr("model file (*.mdl)"));
  if (fileName.isEmpty())
    {
      return false;
    }
  if(fileName.section('.', -1, -1)!="mdl")fileName+=".mdl";
  
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
/*!
  \class StateStackGroup
  
  \brief StateStackGroup  allows the user to define a state variable by selecting a state. On the
  interface, There is an exclusive radio button group, if the user checks one
  of them, the interface of the selected state will be brought up.
  The state variable will not be written into .var file, it's used to control the conditional variables. 

  The user is also able to add a new module.

  All models are defined under Dom elements "models"

  If \a myelem has attribute "components", then this StateStackGroup will affact  components too. 
  
  
 
 
*/          
   
StateStackGroup::StateStackGroup(  QDomElement& elem,  QWidget *parent )
  : GeneralGroup(elem,  parent){
  //find the DOM elements
  QDomElement myroot = myelem.ownerDocument().documentElement();
  QDomElement elem_models= myroot.firstChildElement("models");
  if(elem_models.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find element models")
                         );
    return;
  }
 
  QDomElement elem_opt= elem.firstChildElement();
  if(elem_opt.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         myelem.tagName() + tr("stack element has no child")
                         );
    return;
  }
  
  //build the interface
  pagesWidget = new QStackedWidget;
  QHBoxLayout* mainLayout= new QHBoxLayout;
  QButtonGroup* buttonGroup =  new QButtonGroup(this);
  QVBoxLayout* buttonLayout= new QVBoxLayout;
  mainLayout->addLayout(buttonLayout);
  int count=0;
  int current = 0;
  if(elem.hasAttribute("current"))current=elem.attribute("current").toInt();
  else elem.setAttribute("current", 0);
 
  
  //each child is a button and a page 
  for (; !elem_opt.isNull(); elem_opt = elem_opt.nextSiblingElement(), count++) {
   
    //set up the blank page
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
    buttonGroup->addButton(aButton, count);
    if(count==current)aButton->setChecked(true);
    
   
    //if there are models, put models on the page
    if(elem_opt.attribute("element")=="models"){
      QDomElement elem_mod= elem_opt.firstChildElement();
      QGroupBox* modelGroup = new QGroupBox;
      if(!elem_mod.isNull()){
        QHBoxLayout* modelsLayout = new QHBoxLayout;
        modelGroup->setLayout(modelsLayout);
        modelGroup->setFlat(true);
        for (; !elem_mod.isNull(); elem_mod = elem_mod.nextSiblingElement()){
          QDomElement the_model = elem_models.firstChildElement(elem_mod.tagName());
          ChoiceGroup* modelChoice = new ChoiceGroup(the_model);
          if(elem_mod.hasAttribute("condition")){
            if( !conditionIsSatisfied(myroot, elem_mod.attribute("condition"))){
              modelChoice->hide();
              elem_mod.setAttribute("conditionSatisfied", "false");
            }else{
              modelChoice->show();
              elem_mod.setAttribute("conditionSatisfied", "true");
            }
          }else{
          }
          
          modelsLayout->addWidget(modelChoice);
          connect(modelChoice, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
          connect(this, SIGNAL(showStatus(const bool &)), modelChoice, SLOT(updateShowStatus(const bool &)));
          connect(this, SIGNAL(stateChanged()), modelChoice, SLOT(changeState()));
        }
      }
      
      pageLayout->addWidget(modelGroup);
    }
   
    if(elem_opt.hasAttribute("define")){
      QPushButton* addButton=new QPushButton("Add a Model");
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
  emit textChanged(currentText());
}

 /*!
   Add a new model. The user first select which model to define, then the window will show.
   Three kinds of elements for the new model is allowed: all, CpWindow and ChemistryMdl
 */ 
void StateStackGroup::add(){
  QDomElement myroot = myelem.ownerDocument().documentElement();
  QDomElement elm = myelem.firstChildElement();
  if(elm.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("stack elememt has no child")
                         );
     return;
  }
  for(int i = 0; i < myelem.attribute("current").toInt(); i++)elm = elm.nextSiblingElement();
  
  if(elm.hasAttribute("define")){
    QStringList define_models = elm.attribute("define").split(",");
    if(define_models.size() == 0) return;
    bool ok;
    QString item = QInputDialog::getItem(this, tr("select a model to define:"),
                                         tr("models:"), define_models, 0, false, &ok);
    
    if (ok && !item.isEmpty()){
      QDomElement mdl_elem= myroot.firstChildElement("models");
      if(mdl_elem.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find element 'models' ")
                               );
        return;
      }
      mdl_elem = mdl_elem.firstChildElement(item);
      if(mdl_elem.isNull()){
        QMessageBox::warning(window(), tr(".xml"),
                             tr("can not find the model ")+item
                             );
        return;
      }
        
      
      if(mdl_elem.attribute("element")=="all"){
        
        AllVBWindow* window= new AllVBWindow(mdl_elem);
        window->show();
      }else if(mdl_elem.attribute("element")=="CpWindow"){
        CpWindow* window = new CpWindow;
        window->show();
      }else if(mdl_elem.attribute("element")=="ChemistryMdl"){
        ChemistryMdl* window = new ChemistryMdl;
        window->show();
      }
    
    }
  }
}

/*!
  return a state in the form of "variableName=value",
  \a variableName is the tagName of \a myelem, or the attribute "name" of \a myelem if its defined.
  \a value is the tagName of current child element. 
 */
QString StateStackGroup::myState(){
  QString name =myelem.tagName();
  if(myelem.hasAttribute("name"))name = myelem.attribute("name");
  int current = 0;
  if(myelem.hasAttribute("current"))current=myelem.attribute("current").toInt();
  QDomElement elt = myelem.firstChildElement();
  for(int i = 0; i < current; i++)elt = elt.nextSiblingElement();
  return name+"="+elt.tagName();
}

/*!
  update the attributes \a currentText, \a status, and \a components,
  emit signal componentsChanged() and textChanged()
 */
void StateStackGroup::updateCurrentText(){
  emit textChanged(currentText());
}
/*!
  update the attributes \a currentText, \a status, and \a components,
  emit signal componentsChanged()
 */
QString StateStackGroup::currentText(){
  //go to current child element
  QDomElement myroot = myelem.ownerDocument().documentElement();
  int currentIndex= myelem.attribute("current").toInt();
  QDomElement elt = myelem.firstChildElement();
  for(int i = 0; i < currentIndex; i++)elt=elt.nextSiblingElement(); 
  QString tmp="";
  //add comments
  if(myelem.hasAttribute("comments")) tmp += myelem.attribute("comments")+"\n";
  bool done = true;

  if(elt.attribute("element")=="models"){
    QDomElement elem_models= myroot.firstChildElement("models");

    QDomElement elem_mod= elt.firstChildElement();
    for (; !elem_mod.isNull(); elem_mod = elem_mod.nextSiblingElement()){
      QDomElement the_model = elem_models.firstChildElement(elem_mod.tagName());
      if((the_model.hasAttribute("condition")
          && the_model.attribute("conditionSatisfied")=="true")||
         (!the_model.hasAttribute("condition"))){
          tmp += the_model.attribute("currentText") + "\n";
          if(the_model.attribute("status")!="done")done = false;//if any model not done, done is false
          
          QDomElement cmpnt_elem = the_model.firstChildElement();
          
          for(int i = 0; i < the_model.attribute("current").toInt(); i++)cmpnt_elem=cmpnt_elem.nextSiblingElement();
         
          if(cmpnt_elem.hasAttribute("components")){
            if(cmpnt_elem.attribute("components") != myroot.firstChildElement("components").attribute("components")){
              myroot.firstChildElement("components").setAttribute("components", cmpnt_elem.attribute("components"));
              emit componentsChanged();
            }
          }
        }   
    }
  
  }
  myelem.setAttribute("currentText", tmp);
  if(done) myelem.setAttribute("status", "done");
  else myelem.removeAttribute("status");
  return tmp;
}
  




 /*!
   change page, update the attribute \a current, 
   emit signal stateChanged(), textChanged(), and indirectly componentsChanged()
 */    
void StateStackGroup::changePage(int i)
{
  myelem.setAttribute("current", i);
  pagesWidget->setCurrentIndex(i);
  QString tmp = myState();
  emit stateChanged(tmp);
  emit textChanged(currentText());
}




VarPage::VarPage(QDomElement& theelem,  QWidget* parent): GeneralGroup(theelem,parent)
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
    VarPanel* bdCndPage;
    if(elem.attribute("element")=="panel"){
      bdCndPage = new VarPanel(elem);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(updateCurrentText()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool &)), bdCndPage, SLOT(updateShowStatus(const bool &)));
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



void VarPage::updateCurrentText(){
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

///////////////////////////////////////////////////////////////////////////
//  public:
//    QWidget* createEditor(QWidget* parent,
//                          const QStyleOptionViewItem &option,
//                          const QModelIndex &index) const;
//
//  Simply inverts the visibility state of the boundary double-clicked on
//  in the boundary visibility model.
///////////////////////////////////////////////////////////////////////////


QWidget* showDelegate::createEditor(QWidget*, const QStyleOptionViewItem&,
				    const QModelIndex &index) const
{
  QString value = index.data(Qt::EditRole).toString();
  QAbstractItemModel* model = const_cast<QAbstractItemModel*>(index.model());
  if (value == "hide") {
    model->setData(index, "show");
    
  } else {
    model->setData(index, "hide");
  }

  return NULL;
}






  
                                                              



bool colorDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
                                const QStyleOptionViewItem &,
                                const QModelIndex &index){
  
                                                              
 
  if (event->type() == QEvent::MouseButtonPress) {
    QColor oldColor =   qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->background().color();
    QColor color = QColorDialog::getColor(oldColor);
    
    if(color.isValid())qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->setBackground( QBrush(color));
    else qobject_cast<const QStandardItemModel*>(model)->item(index.row(), index.column())->setBackground( QBrush(oldColor));
    return false; //so that the selection can change
  }
  
  return true;
}



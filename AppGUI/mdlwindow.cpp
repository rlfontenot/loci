#include <QtGui>
#include <QFile>
#include <QString>


#include "mdlwindow.h"
#include "pages.h"
#include "getfile.h"









DoubleEditDelegate::DoubleEditDelegate(QObject *parent)
  : QItemDelegate(parent)
{
}

QWidget *DoubleEditDelegate::createEditor(QWidget *parent,
                                         const QStyleOptionViewItem &/* option */,
                                         const QModelIndex &/* index */) const
{
  DoubleEdit *editor = new DoubleEdit(parent);
  return editor;
}

void DoubleEditDelegate::setEditorData(QWidget *editor,
                                      const QModelIndex &index) const
{
  double value = index.model()->data(index, Qt::EditRole).toDouble();
  
     DoubleEdit *floatEdit = static_cast<DoubleEdit*>(editor);
     floatEdit->setValue(value);
}

void DoubleEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                    const QModelIndex &index) const
{
  DoubleEdit *floatEdit = static_cast<DoubleEdit*>(editor);
    
  double value = floatEdit->value();
  
  model->setData(index, value, Qt::EditRole);
}

void DoubleEditDelegate::updateEditorGeometry(QWidget *editor,
                                             const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
  editor->setGeometry(option.rect);
}

ListGroup::ListGroup(const QString& title, QWidget* parent):QGroupBox(title, parent){

  QLabel *aLabel = new QLabel(tr("add an item"));
  fEdit = new DoubleEdit;
  
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(aLabel);
  buttonLayout->addWidget(fEdit);
  
  QPushButton * removeButton = new QPushButton("remove an item");
  
  display = new QLabel("");
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(buttonLayout);
  mainLayout->addWidget(removeButton);
  mainLayout->addWidget(display);
  connect(fEdit, SIGNAL(returnPressed()), this, SLOT(add()));
  connect(removeButton, SIGNAL(clicked()), this, SLOT(remove()));

  setLayout(mainLayout);
}

void ListGroup::add(){
  QString text = display->text();
  if(text!="")text += ","+fEdit->text();
  else text += fEdit->text();
  fEdit->setText("");
  display->setText(text);
  emit textChanged();
}
void ListGroup::clear(){
  fEdit->clear();
  display->setText("");
}
void ListGroup::remove(){
  QString text = display->text();
  QStringList items = text.split(",");
  items.removeLast();
  display->setText(items.join(","));
  emit textChanged();
}

QString ListGroup::text(){
  QString text = display->text();
  if(text.contains(',')) text = "["+text+"]";
  return text;
}

SpeciesGroup::SpeciesGroup(const QString& title, QWidget* parent):QGroupBox(title, parent){
  
  QGroupBox* nameGroup= new QGroupBox(tr("name of the gas"));
  nameBGroup = new QButtonGroup(this);
  nameBGroup->setExclusive(true);
  QRadioButton* nradio0 = new QRadioButton("named");
  QRadioButton* nradio1 = new QRadioButton("chemical");
  nradio1->setChecked(true);
  nameBGroup->addButton(nradio0);
  nameBGroup->addButton(nradio1);
  nameBGroup->setId(nradio0, 0);
  nameBGroup->setId(nradio1, 1);
  nameEdit = new QLineEdit;

  connect(nameBGroup, SIGNAL(buttonClicked(int)), this, SLOT(nameBGroupClicked(int)));
  connect(nameEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  QGridLayout *nameLayout = new QGridLayout;
  
  nameLayout->addWidget(nradio0, 0, 0, 1, 1);
  nameLayout->addWidget(nradio1, 0, 1, 1, 1);
  nameLayout->addWidget(nameEdit,1, 0, 1, 2);
  nameGroup->setLayout(nameLayout);
  
  QGroupBox  *replaceGroup = new QGroupBox(tr("replacement or augment"));
  QHBoxLayout *replaceLayout = new QHBoxLayout;
  replaceBGroup = new QButtonGroup(this);
  replaceBGroup->setExclusive(true);
  QRadioButton* rradio0 = new QRadioButton("replace");
  QRadioButton* rradio1 = new QRadioButton("augument");
  rradio0->setChecked(true);
  replaceBGroup->addButton(rradio0);
  replaceBGroup->addButton(rradio1);
  replaceBGroup->setId(rradio0, 0);
  replaceBGroup->setId(rradio1, 1);
  
  replaceLayout->addWidget(rradio0);
  replaceLayout->addWidget(rradio1);
  replaceGroup->setLayout(replaceLayout);


  connect(replaceBGroup, SIGNAL(buttonClicked(int)), this, SLOT(replaceBGroupClicked(int)));
  mGroup = new QGroupBox(tr("molecular mass"));
  QHBoxLayout *mLayout = new QHBoxLayout;
  mEdit = new DoubleEdit;
  mEdit->setBottom(0.0);
  mLayout->addWidget(mEdit);
  mGroup->setLayout(mLayout);
  mGroup->setCheckable(true);
  mGroup->setChecked(false);
  connect(mGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(mEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  sGroup = new QGroupBox(tr("reference entropy"));
  QHBoxLayout *sLayout = new QHBoxLayout;
  srefEdit = new DoubleEdit;
  sLayout->addWidget(srefEdit);
  sGroup->setLayout(sLayout);
  sGroup->setCheckable(true);
  sGroup->setChecked(false);
  connect(sGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(srefEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  

  tGroup = new QGroupBox(tr("reference temperature"));
  tGroup->setWhatsThis("reference temperature used for Kc");
  QHBoxLayout *tLayout = new QHBoxLayout;
  trefEdit = new DoubleEdit;
  tLayout->addWidget(trefEdit);
  tGroup->setLayout(tLayout);
  tGroup->setCheckable(true);
  tGroup->setChecked(false);
  connect(tGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(trefEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  pGroup = new QGroupBox(tr("reference pressure"));
  pGroup->setWhatsThis("reference pressure used for Kc");
  QHBoxLayout *pLayout = new QHBoxLayout;
  prefEdit = new DoubleEdit;
  pLayout->addWidget(prefEdit);
  pGroup->setLayout(pLayout);
  pGroup->setCheckable(true);
  pGroup->setChecked(false);
  connect(pGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(prefEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  hGroup = new QGroupBox(tr("reference enthalpy"));
  hGroup->setWhatsThis("reference enthalpy used for Kc");  
  QHBoxLayout *hLayout = new QHBoxLayout;
  hrefEdit = new DoubleEdit;
  hLayout->addWidget(hrefEdit);
  hGroup->setLayout(hLayout);
  hGroup->setCheckable(true);
  hGroup->setChecked(false);
  connect(hGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(hrefEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  mfGroup = new QGroupBox(tr("mass fraction"));
  mfGroup->setWhatsThis("mass fraction used to specify default mixture");
  QHBoxLayout *mfLayout = new QHBoxLayout;
  mfEdit = new DoubleEdit;
  mfLayout->addWidget(mfEdit);
  mfGroup->setLayout(mfLayout);
  mfGroup->setCheckable(true);
  mfGroup->setChecked(false);
  connect(mfGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(mfEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  thetaCpGroup = new QGroupBox(tr("vibrational temperatures or curve fit for cp "));
  thetaCpGroup->setCheckable(true);
  thetaCpGroup->setChecked(false);
  QHBoxLayout *thetaCpLayout = new QHBoxLayout;
  thetaCpBGroup = new QButtonGroup(this);
  thetaCpBGroup->setExclusive(true);
  QRadioButton* tradio0 = new QRadioButton("theta_v");
  QRadioButton* tradio1 = new QRadioButton("cp");
  tradio0->setChecked(true);
  thetaCpBGroup->addButton(tradio0);
  thetaCpBGroup->addButton(tradio1);
  thetaCpBGroup->setId(tradio0, 0);
  thetaCpBGroup->setId(tradio1, 1);
  
  thetaCpLayout->addWidget(tradio0);
  thetaCpLayout->addWidget(tradio1);
  
  
  thetavGroup = new ListGroup;
  cpGroup = new CpGroup;
  
  
  pagesWidget = new QStackedWidget;
  pagesWidget->addWidget(thetavGroup);
  pagesWidget->addWidget(cpGroup);
  
  QVBoxLayout *tcLayout = new QVBoxLayout;
  tcLayout->addLayout(thetaCpLayout);
  tcLayout->addWidget(pagesWidget);
  thetaCpGroup->setLayout(tcLayout);
  connect(thetaCpBGroup, SIGNAL(buttonClicked(int)), this, SLOT(changePage(int)));
  connect(thetavGroup, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(cpGroup, SIGNAL(textChanged()), this, SLOT(setText()));
  
  QVBoxLayout* mainLayout = new QVBoxLayout;
  QGridLayout* vBoxLayout = new QGridLayout;
  
  vBoxLayout->addWidget(nameGroup, 0, 0, 1, 1);
  vBoxLayout->addWidget(replaceGroup, 0,1,1,1);
  vBoxLayout->addWidget(mGroup,0, 2, 1, 1);
  vBoxLayout->addWidget(sGroup, 1, 0, 1, 1);
  vBoxLayout->addWidget(tGroup, 1, 1, 1, 1);
  vBoxLayout->addWidget(pGroup, 1, 2, 1, 1); 
  vBoxLayout->addWidget(hGroup, 1, 3, 1, 1);
  vBoxLayout->addWidget(mfGroup, 0, 3, 1, 1);
  mainLayout->addLayout(vBoxLayout);
  mainLayout->addWidget(thetaCpGroup);
  display = new QLabel;
  mainLayout->addWidget(display);
  QPushButton* acceptButton = new QPushButton(tr("&Accept"));
  connect(acceptButton, SIGNAL(clicked()), this, SLOT(clear()));
  mainLayout->addWidget(acceptButton);
  setLayout(mainLayout);
}
void SpeciesGroup::clear(){
  if(nameEdit->text()=="") return;
  emit accept();
  mEdit->clear();
  hrefEdit->clear();
  srefEdit->clear();
  trefEdit->clear();
  prefEdit->clear();
  mfEdit->clear();
  thetavGroup->clear();
  cpGroup->clear();
  nameEdit->clear();
  mGroup->setChecked(false);
  mfGroup->setChecked(false);
  sGroup->setChecked(false);
  tGroup ->setChecked(false);
  pGroup->setChecked(false);
  hGroup->setChecked(false) ;
  thetaCpGroup->setChecked(false);
  display->setText("");

}

  

void SpeciesGroup::changePage(int i){
  pagesWidget->setCurrentIndex(i);
  setText();
}
void SpeciesGroup::setText(){
  display->setText(text());
}
QString SpeciesGroup::text(){
  QString text_name = nameEdit->text();
  QString text_infix;
  if(replaceBGroup->checkedId() == 0) text_infix  = " = ";
  else text_infix = " : ";
  QString clean_text ;
  if(mGroup->isChecked()) clean_text  += "m=" + mEdit->text()+", ";
  if(hGroup->isChecked()) clean_text += "href=" + hrefEdit->text() + ", ";
  if(sGroup->isChecked()) clean_text += "sref=" + srefEdit->text() + ", ";
  if(tGroup->isChecked()) clean_text += "Tref=" + trefEdit->text() + ", ";
  if(pGroup->isChecked()) clean_text += "Pref=" + prefEdit->text() + ", ";
  if(mfGroup->isChecked()) clean_text += "mf=" +mfEdit->text() + ", ";
  if(thetaCpGroup->isChecked()){
    if(thetaCpBGroup->checkedId() == 0)clean_text += "\n     theta_v=" + thetavGroup->text();
    else clean_text += "\n" + cpGroup->text();
  }
  if(clean_text.right(2)==QString(", "))clean_text = clean_text.section(',', 0, -2);
  if(clean_text != "") return text_name + text_infix + "<"+clean_text+">"+";\n";
  else return text_name + ";\n";
}

void CpWindow::setText(){
  display->setText(text());
}


CpGroup ::CpGroup(QWidget* parent):QGroupBox(tr("curve fit for Cp"),parent){
 
  QGroupBox *equationGroup = new QGroupBox(tr("please specify equation type "));
  QHBoxLayout *equationLayout = new QHBoxLayout;
  shomateButton = new QRadioButton(tr("shomate"));
  polyButton = new QRadioButton(tr("polynomial"));
  shomateButton->setChecked(true);                               
  equationLayout->addWidget(shomateButton);
  equationLayout->addWidget(polyButton);
  equationGroup->setLayout(equationLayout);

  connect(shomateButton, SIGNAL(toggled(bool)), this, SIGNAL(textChanged()));
  QGroupBox *nGroup = new QGroupBox(tr("please specify number of  temperature intervals"));
  QHBoxLayout *nLayout = new QHBoxLayout;
  nEdit = new IntEdit;
  nEdit->setRange(1, 100);
  nEdit->setValue(3);
  
  nLayout->addWidget(nEdit);
  nGroup->setLayout(nLayout);
  connect(nEdit, SIGNAL(valueChanged(int)), this, SLOT(setNumberOfIntervals(int))); 

  
  

  
  
  model = new  QStandardItemModel();
  
  //   }
  model->setColumnCount(6);
  model->setHeaderData(0, Qt::Horizontal, QObject::tr("Temperature"));
  model->setHeaderData(1, Qt::Horizontal, QObject::tr("      A      "));
  model->setHeaderData(2, Qt::Horizontal, QObject::tr("      B      "));
  model->setHeaderData(3, Qt::Horizontal, QObject::tr("      C      "));
  model->setHeaderData(4, Qt::Horizontal, QObject::tr("      D      "));
  model->setHeaderData(5, Qt::Horizontal, QObject::tr("      E      "));
  
  QItemSelectionModel *selections = new QItemSelectionModel(model);
  DoubleEditDelegate* delegate = new DoubleEditDelegate(this);
  
  tableView = new QTableView(this);
  tableView->setModel(model);
  tableView->setSelectionModel(selections);
  tableView->setSelectionMode(QAbstractItemView::SingleSelection);
  tableView->setItemDelegate(delegate);
  tableView->resizeColumnsToContents();
  
  numberOfIntervals =3;
  model->setRowCount(4);

  
  connect(model, SIGNAL(itemChanged(QStandardItem*)), this, SIGNAL(textChanged()));
  QVBoxLayout* mainLayout = new QVBoxLayout;
  
 
  mainLayout->addWidget(equationGroup);
  mainLayout->addWidget(nGroup);
  mainLayout->addWidget(tableView);
  setLayout(mainLayout);

}
void CpGroup::setNumberOfIntervals(int n){
  numberOfIntervals = n;
  model->setRowCount(n+1);
  tableView->setModel(model);
  tableView->show();
}
void CpGroup::clear(){
  model->setRowCount(0);
  
  nEdit->setValue(3);
}

QString CpGroup::text(){
  QString tab = "     ";
  QString text = tab + "cp=[";
  QString prefix;
  tab += "    "; 
     if(shomateButton->isChecked())prefix =tab+"shomate";
     else prefix=tab+"poly";
     
     for(int i =0; i < (nEdit->value()); i++){
       if(i==0) text += (model->index(i, 0)).data(Qt::EditRole).toString()+",\n";
       else text += tab+(model->index(i, 0)).data(Qt::EditRole).toString()+",\n";
       text += prefix+"(";
       for(int j=1; j<5; j++){
         text += (model->index(i, j)).data(Qt::EditRole).toString()+",";
       }
       text += (model->index(i, 5)).data(Qt::EditRole).toString()+"),\n";
     }
     text += tab+(model->index(nEdit->value(), 0)).data(Qt::EditRole).toString()+"]";
     return text;
}


void SpeciesGroup::replaceBGroupClicked(int bid){
  if(bid==0)mGroup->show();
  else{
    mGroup->setChecked(false);
    mGroup->hide();
  }
  setText();
}

void SpeciesGroup::nameBGroupClicked(int bId){
  if(bId==0){
    mGroup->show();
    mGroup->setChecked(true);
    nameEdit->setInputMask("_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
  }else{
    mGroup->show();
    mGroup->setChecked(false);
    nameEdit->setInputMask(QString());
  }
  setText();
}



QString SpeciesGroup::getSpecies(){
  return nameEdit->text();
}

ExpNu::ExpNu(const QStringList& sp, QWidget* parent):
  QGroupBox(tr("Exp_nu"), parent),species(sp){
  setCheckable(true);
  setChecked(false);   
  bool floatAllowed = true;
  QGridLayout* mainLayout = new QGridLayout;
  
  for(int i = 0; i < sp.count(); i++){
    NumberString *aSp = new NumberString(sp.at(i),floatAllowed);
    mainLayout->addWidget(aSp, i/5, i%5, 1, 1);
    exp_nu<<aSp;
    connect(aSp, SIGNAL(textChanged()), this, SIGNAL(textChanged()));
    connect(aSp, SIGNAL(toggled(bool)), this, SIGNAL(textChanged()));
  }
  setLayout(mainLayout);
}
QString ExpNu::text(){
   QString text="[";
   for(int i = 0; i < species.count(); i++){
     if(exp_nu.at(i)->isChecked()){
       text += exp_nu.at(i)->cleanText()+", ";
     }
   }
   if(text.right(2)==QString(", "))text = text.section(',', 0, -2); 
   text += "]";
   return text;
}
     
void ExpNu::clear(){
  for(int i = 0; i < species.count(); i++){
    exp_nu.at(i)->clear();
  }
  setChecked(false);
}

RateModifier::RateModifier(const QString& vname, int n, QWidget* parent):
  QGroupBox(tr("Rate Modifier:")+vname, parent), name(vname), numValues(n){
  setCheckable(true);
  setChecked(false);   
 
  QGridLayout* mainLayout = new QGridLayout;
  
  for(int i = 0; i < n; i++){
    DoubleEdit *editor = new DoubleEdit;
    mainLayout->addWidget(editor, i/2, i%2, 1, 1);
    edits<<editor;
    connect(editor, SIGNAL(valueChanged(double)), this, SIGNAL(textChanged()));
   
  }
  setLayout(mainLayout);
}

QString RateModifier::text(){
   QString text=name+"(";
   for(int i = 0; i < numValues; i++){
     text += edits[i]->text()+", ";
     if((i+1)%2==0) text+="\n";
   }
   text = text.section(',', 0, -2); 
   text += ")";
   return text;
}
     
void RateModifier::clear(){
  for(int i = 0; i < numValues; i++){
    edits[i]->clear();
  }
  setChecked(false);
}

KFKC::KFKC(const QString& title, QWidget* parent):
  QGroupBox(title, parent){
  setCheckable(true);
  setChecked(false);
  arr = new QRadioButton(tr("Arrhennius"));
  arr->setChecked(true);
  ther = new QRadioButton(tr("Thermodynamic"));
  f1 = new DoubleEdit;
  f2 = new DoubleEdit;
  f3 = new DoubleEdit;

  connect(arr, SIGNAL(toggled(bool)), this, SIGNAL(textChanged()));
  connect(arr, SIGNAL(toggled(bool)), this, SLOT(arrToggled(bool)));
  connect(ther, SIGNAL(toggled(bool)), this, SIGNAL(textChanged()));
  connect(f1, SIGNAL(textChanged(const QString&)), this, SIGNAL(textChanged()));
  connect(f2, SIGNAL(textChanged(const QString&)), this, SIGNAL(textChanged()));
  connect(f3, SIGNAL(textChanged(const QString&)), this, SIGNAL(textChanged()));
  
  QHBoxLayout* hBoxLayout = new QHBoxLayout;
  hBoxLayout->addWidget(arr);
  hBoxLayout->addWidget(f1);
  hBoxLayout->addWidget(f2);
  hBoxLayout->addWidget(f3);
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(hBoxLayout);
  mainLayout->addWidget(ther);
  setLayout(mainLayout);
}



    
  
void KFKC::arrToggled(bool checked){
  if(checked){
    f1->setDisabled(false);
    f2->setDisabled(false);
    f3->setDisabled(false);
  }else{
    f1->setDisabled(true);
    f2->setDisabled(true);
    f3->setDisabled(true); 
  }
}
QString KFKC::text(){
 
  if(ther->isChecked()) return "Thermodynamic";
  else return "Arrhenius(" + f1->text() + "," + f2->text() + "," + f3->text() + ")";
}
    
void KFKC::checkTher(){
  ther->setChecked(true);
}

NumberString::NumberString(const QString& title, bool fa, QWidget* parent):
  QGroupBox(title, parent),floatAllowed(fa){
  setCheckable(true);
  setChecked(false);
  if(floatAllowed){
    edit=new DoubleEdit;
    qobject_cast<DoubleEdit*>(edit)->setBottom(0.0);
    connect(edit, SIGNAL(textChanged(const QString&)), this, SIGNAL(textChanged()));
  }else{
    edit = new QSpinBox;
    qobject_cast<QSpinBox*>(edit)->setMinimum(1);
    connect(edit, SIGNAL(valueChanged(const QString&)), this, SIGNAL(textChanged()));
  }

  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->addWidget(edit);
  setLayout(mainLayout);

}
void NumberString::clear(){
 if(isCheckable()) setChecked(false);
 if(floatAllowed) qobject_cast<DoubleEdit*>(edit)->clear();
 else qobject_cast<QSpinBox*>(edit)->setValue(1);
}

QString NumberString::text(){
  QString text;
  if(floatAllowed){
    
    if(qobject_cast<DoubleEdit*>(edit)->value() != 1)  text += qobject_cast<DoubleEdit*>(edit)->text()+" "+title();
    else if(qobject_cast<DoubleEdit*>(edit)->value() != 0) text +=title();
                
    
  }else{
    
    if(qobject_cast<QSpinBox*>(edit)->value() != 1) text += qobject_cast<QSpinBox*>(edit)->cleanText()+" "+title();
    else text += title();
  }
  return text;
}
 
 QString NumberString::cleanText(){
   QString text;
   if(floatAllowed)
     text += title() + "=" +qobject_cast<DoubleEdit*>(edit)->text();
   else text += title() + "="+ qobject_cast<QSpinBox*>(edit)->cleanText();
   return text;
 }

      
Reactants::Reactants(const QStringList& sp, const QString& title, QWidget* parent):
  QGroupBox(title, parent){
  bool floatAllowed = false;
  QGridLayout* mainLayout = new QGridLayout;
  
  for(int i = 0; i < sp.count(); i++){
    if(sp.at(i).at(0)=='_'){
      floatAllowed = true;
      break;
    }
  }

  for(int i = 0; i < sp.count(); i++){
    NumberString *aSp = new NumberString(sp.at(i), floatAllowed);
    mainLayout->addWidget(aSp, i/5, i%5, 1, 1);
    species<<aSp;
    connect(aSp, SIGNAL(textChanged()), this, SIGNAL(textChanged()));
    connect(aSp, SIGNAL(toggled(bool)), this, SIGNAL(textChanged()));
  }
  setLayout(mainLayout);
}


void Reactants::clear(){

  for(int i = 0; i < species.count(); i++){
    species.at(i)->clear();
  }
}

QString Reactants::text(){
  QString text;
  for(int i =0; i < species.count(); i++){
    if(species.at(i)->isChecked()) {
      if(text != "")text += " + " + species.at(i)->text();
      else text +=  species.at(i)->text();
    }
  }
  return text;
}



ReactionGroup::ReactionGroup(const QStringList& species, const QString& title, QWidget* parent):QGroupBox(title, parent){
  
  if(species.count()==0){
    QMessageBox::warning(this, tr("Reaction Group"), tr("Please specify species first")); 
    return;
  }
  
  QHBoxLayout* hBoxLayout = new QHBoxLayout;
  reactants = new Reactants(species, tr("please specify reactants:"), this);
  products = new Reactants(species, tr("please specify products:"), this);

  direction = new QComboBox;
  direction->addItem(" -> ");
  direction->addItem(" <-> ");

  connect(reactants, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(products, SIGNAL(textChanged()), this, SLOT(setText()));

  hBoxLayout->addWidget(reactants);
 
  hBoxLayout->addWidget(products);

  mbody = new QCheckBox(tr("is M-body reactions?"));
  mbody->setChecked(false);
  
  QHBoxLayout* vBoxLayout = new QHBoxLayout;
  vBoxLayout->addWidget(mbody);
  vBoxLayout->addStretch(10);
  vBoxLayout->addWidget(new QLabel(tr("direction:")));
  vBoxLayout->addWidget(direction);
  vBoxLayout->addStretch(20);
                        

  connect(mbody, SIGNAL(stateChanged(int)), this, SLOT(setText()));
  connect(direction, SIGNAL(currentIndexChanged(int)), this, SLOT(setText()));
  connect(direction, SIGNAL(currentIndexChanged(int)), this, SLOT(directionChanged()));

  kf = new KFKC(tr("please specify Kf:"));
  kc = new KFKC(tr("please specify Kc:"));
  kc->checkTher();
  kc->hide();

  expnu= new ExpNu(species);
  minMf = new NumberString(tr("MinMF"), true);
  mdf = new RateModifier(tr("pressure"),2);
  
  connect(kf, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(kc, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(kf, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(kc, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(expnu, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(expnu, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(minMf, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(minMf, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(mdf, SIGNAL(textChanged()), this, SLOT(setText()));
  connect(mdf, SIGNAL(clicked(bool)), this, SLOT(setText()));

  QVBoxLayout* mainLayout = new QVBoxLayout;
  display = new QLabel;
  QPushButton* acceptButton = new QPushButton(tr("accept"));

  connect(acceptButton, SIGNAL(clicked()), this, SLOT(clear()));

  mainLayout->addLayout(vBoxLayout);
  mainLayout->addLayout(hBoxLayout);
  mainLayout->addWidget(kf);
  mainLayout->addWidget(kc);
  mainLayout->addWidget(expnu);
  mainLayout->addWidget(minMf);
  mainLayout->addWidget(mdf); 
  mainLayout->addWidget(display);
  mainLayout->addWidget(acceptButton);
  setLayout(mainLayout);
}
void ReactionGroup::clear(){
  emit accept();
  reactants->clear();
  products->clear();
  kf->clear();
  kc->clear();
  expnu->clear();
  minMf->clear();
  mdf->clear();
  display->setText("");
}
void KFKC::clear(){
  f1->clear();
  f2->clear();
  f3->clear();
}
QString ReactionGroup::text(){
  QString text = reactants->text();
  if(mbody->isChecked()) text += " +M";
  text += direction->currentText();
  text += products->text();
  if(mbody->isChecked()) text += " +M";
  if(kf->isChecked() || kc->isChecked()
     ||expnu->isChecked()||minMf->isChecked()
     ||mdf->isChecked()) text += " = <";
  if(kf->isChecked()) text += "Kf = " + kf->text();
  
  if(kc->isChecked()){
    if(kf->isChecked()) text += ",\n                     ";
    text += "Kc = " +  kc->text();
  }
  if(expnu->isChecked()){
    if(kf->isChecked()||kc->isChecked())text += ",\n                     ";
    text += "exp_nu = " +  expnu->text();
  }
  if(minMf->isChecked()){
    if(kf->isChecked()||kc->isChecked()||expnu->isChecked())text += ",\n                     ";
    text +=   minMf->cleanText();
  }
  if(mdf->isChecked()){
    if(kf->isChecked()||kc->isChecked()||expnu->isChecked()||minMf->isChecked())text += ",\n                     ";
    text += "rate_modifier = " +  mdf->text();
  }
  
  if(kf->isChecked() || kc->isChecked()||expnu->isChecked()
     ||minMf->isChecked()||mdf->isChecked()) text += ">";
  text += ";\n";
  
  return text;
}
void ReactionGroup::setText(){
  display->setText(text());
}
void ReactionGroup::directionChanged(){
  if(direction->currentIndex() == 0){
    kc->setChecked(false);
    kc->hide();
  }else{
    kc->show();
  }
}


SpeciesWindow::SpeciesWindow( const QString& title, QWidget* parent):QGroupBox(title, parent){
  speciesList  = new QListWidget;
  speciesPage = new SpeciesGroup;
  
  connect(speciesPage, SIGNAL(accept()), this, SLOT(add()));
  QPushButton *clearSpeciesBtn= new QPushButton(tr("clear all species"));
  connect(clearSpeciesBtn, SIGNAL(clicked()), this, SLOT(clear()));
  
 
  
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(clearSpeciesBtn);
  
  QVBoxLayout* listLayout = new QVBoxLayout;
  listLayout->addWidget(speciesList);
  listLayout->addLayout(buttonLayout);
  QHBoxLayout* mainLayout = new QHBoxLayout;
  mainLayout->addWidget(speciesPage);
  mainLayout->addLayout(listLayout);
  setLayout(mainLayout);
}

void SpeciesWindow::add(){
  QString item = speciesPage->text();
  species << speciesPage->getSpecies();
  speciesList->addItem(item);
}

void SpeciesWindow::clear(){
  speciesList->clear();
  species.clear();
}
QStringList SpeciesWindow:: getSpecies(){
  return species;
}





ReactionWindow::ReactionWindow( const QStringList& sps, const QString& title,QWidget* parent):QGroupBox(title, parent){
  species = sps;
  reactionList  = new QListWidget;
  reactionPage = new ReactionGroup(species);;
  
  connect(reactionPage, SIGNAL(accept()), this, SLOT(add()));
  QPushButton *clearReactionBtn= new QPushButton(tr("clear all reaction"));
  connect(clearReactionBtn, SIGNAL(clicked()), this, SLOT(clear()));
  
 
  
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(clearReactionBtn);
  
  QVBoxLayout* listLayout = new QVBoxLayout;
  listLayout->addWidget(reactionList);
  listLayout->addLayout(buttonLayout);
  mainLayout = new QHBoxLayout;
  mainLayout->addWidget(reactionPage, 1);
  mainLayout->addLayout(listLayout, 2);
  setLayout(mainLayout);
}

QString ReactionWindow::text(){
  QString text = "{\n" ;
  for(int i = 0; i < reactionList->count(); i++){
    text += reactionList->item(i)->text();
  }
  text += "};\n";
  return text;
}
QString SpeciesWindow::text(){
  QString text = "{\n" ;
  for(int i = 0; i < speciesList->count(); i++){
    text += speciesList->item(i)->text();
  }
  text += "};\n";
  return text;
} 

void ReactionWindow::add(){
  QString item = reactionPage->text();
  reactionList->addItem(item);
}

void ReactionWindow::clear(){
  reactionList->clear();
}
QStringList ReactionWindow:: getSpecies(){
  return species;
}
void ReactionWindow::resetSpecies(const QStringList& sps){
  mainLayout->removeWidget(reactionPage);
  species = sps; 
  delete reactionPage;
  reactionPage = new ReactionGroup(species);
  connect(reactionPage, SIGNAL(accept()), this, SLOT(add()));
  mainLayout->insertWidget(0,reactionPage);
}
ChemistryMdl::ChemistryMdl(bool noReaction, QWidget* parent):QWidget(parent){
  QButtonGroup* buttonGroup = new QButtonGroup;
  buttonGroup->setExclusive(true);
  
  QRadioButton* speciesButton = new QRadioButton(tr("Add Species"));
   reactionButton = new QRadioButton(tr("Add Reaction"));
  QRadioButton* saveButton = new QRadioButton(tr("save the model"));
 
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(speciesButton);
  buttonLayout->addWidget(reactionButton);
  buttonLayout->addWidget(saveButton);
 if(noReaction) reactionButton->hide();
  buttonGroup->addButton(speciesButton);
  buttonGroup->addButton(reactionButton);
  buttonGroup->addButton(saveButton);
  buttonGroup->setId(speciesButton, 0);
  buttonGroup->setId(reactionButton, 1);
  buttonGroup->setId(saveButton, 2);
  speciesButton->setChecked(true);
  
  connect(buttonGroup, SIGNAL(buttonClicked(int)), this, SLOT(changePage(int)));
  speciesWindow = new SpeciesWindow;
  center  = new QStackedWidget;
  center->addWidget(speciesWindow);
  reactionWindow = 0;
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addLayout(buttonLayout);
  mainLayout->addWidget(center);
  setLayout(mainLayout);
}
void  ChemistryMdl::changePage(int id){
  QStringList species = speciesWindow->getSpecies();
  switch(id){
  case 0:
    center->setCurrentWidget(speciesWindow);
    break;
  case 1:
    if(species.count()==0){
       QMessageBox::warning(this, tr("Reaction Window"), tr("Please specify species first")); 
       return;
    }else{
      if(reactionWindow == 0){
        reactionWindow = new ReactionWindow(species);
        center->addWidget(reactionWindow);
      }else if(species != reactionWindow->getSpecies()){
        reactionWindow->resetSpecies(species);
      }
      center->setCurrentWidget(reactionWindow);
    }
    break;
  case 2:
    save();
    break;
  }
}

void ChemistryMdl::hideReaction(){
  reactionButton->hide();
}
  
bool ChemistryMdl::save(){
    QString initialPath = QDir::currentPath()+"/untitled.mdl";
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                  initialPath,
                                                    tr("model file (*.mdl)"));
    
    if (fileName.isEmpty()) return false;
      

    if(fileName.section('.', -1, -1)!="mdl")fileName += ".mdl";
    
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
    out << "species = " << speciesWindow->text();
    
    if(reactionWindow) out << "reactions = " << reactionWindow->text();
    else out << "reactions = { \n};";
           
    QApplication::restoreOverrideCursor(); 
    file.close();
    return true;
  }
CpWindow::CpWindow(const QString& title, const QString& pf, QWidget* parent):QGroupBox(title, parent),prefix(pf){
  
  QGroupBox* nameGroup= new QGroupBox(tr("name of the gas"));
  nameBGroup = new QButtonGroup(this);
  nameBGroup->setExclusive(true);
  QRadioButton* nradio0 = new QRadioButton("named");
  QRadioButton* nradio1 = new QRadioButton("chemical");
  nradio1->setChecked(true);
  nameBGroup->addButton(nradio0);
  nameBGroup->addButton(nradio1);
  nameBGroup->setId(nradio0, 0);
  nameBGroup->setId(nradio1, 1);
  nameEdit = new QLineEdit;

  connect(nameBGroup, SIGNAL(buttonClicked(int)), this, SLOT(nameBGroupClicked(int)));
  connect(nameEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));
  
  QGridLayout *nameLayout = new QGridLayout;
  
  nameLayout->addWidget(nradio0, 0, 0, 1, 1);
  nameLayout->addWidget(nradio1, 0, 1, 1, 1);
  nameLayout->addWidget(nameEdit,1, 0, 1, 2);
  nameGroup->setLayout(nameLayout);
  
  QGroupBox  *replaceGroup = new QGroupBox(tr("replacement or augment"));
  QHBoxLayout *replaceLayout = new QHBoxLayout;
  replaceBGroup = new QButtonGroup(this);
  replaceBGroup->setExclusive(true);
  QRadioButton* rradio0 = new QRadioButton("replace");
  QRadioButton* rradio1 = new QRadioButton("augument");
  rradio0->setChecked(true);
  replaceBGroup->addButton(rradio0);
  replaceBGroup->addButton(rradio1);
  replaceBGroup->setId(rradio0, 0);
  replaceBGroup->setId(rradio1, 1);
  
  replaceLayout->addWidget(rradio0);
  replaceLayout->addWidget(rradio1);
  replaceGroup->setLayout(replaceLayout);

  connect(replaceBGroup, SIGNAL(buttonClicked(int)), this, SLOT(replaceBGroupClicked(int)));
  
  mGroup = new QGroupBox(tr("molecular mass"));
  QHBoxLayout *mLayout = new QHBoxLayout;
  mEdit = new DoubleEdit;
  mEdit->setBottom(0.0);
  mLayout->addWidget(mEdit);
  mGroup->setLayout(mLayout);
  mGroup->setCheckable(true);
  mGroup->setChecked(false);
  connect(mGroup, SIGNAL(clicked(bool)), this, SLOT(setText()));
  connect(mEdit, SIGNAL(textChanged(const QString&)), this, SLOT(setText()));     

  cpGroup = new CpGroup;
   connect(cpGroup, SIGNAL(textChanged()), this, SLOT(setText()));  
  QVBoxLayout* mainLayout = new QVBoxLayout;
  QGridLayout* vBoxLayout = new QGridLayout;
  
  vBoxLayout->addWidget(nameGroup, 0, 0, 1, 1);
  vBoxLayout->addWidget(replaceGroup, 0,1,1,1);
  vBoxLayout->addWidget(mGroup,0, 2, 1, 1);
  
  mainLayout->addLayout(vBoxLayout);
  mainLayout->addWidget(cpGroup);
  
  display = new QLabel;
  mainLayout->addWidget(display);
   
  QPushButton* acceptButton = new QPushButton(tr("&Save the model"));

  connect(acceptButton, SIGNAL(clicked()), this, SLOT(save()));
  mainLayout->addWidget(acceptButton);
  setLayout(mainLayout);
}

QString CpWindow::text(){
  QString text ="species = {\n";
  QString text_name = nameEdit->text();
  QString text_infix;
  if(replaceBGroup->checkedId() == 0) text_infix  = " = ";
  else text_infix = " : ";
  QString clean_text ;
  if(mGroup->isChecked()) clean_text  += "m=" + mEdit->text() +", ";
  clean_text += prefix+"\n"  + cpGroup->text();
  text += text_name + text_infix + "<"+clean_text+">"+" ;\n} ;\n\nreactions = {\n} ;";
  return text;
}
bool CpWindow::save(){
  QString initialPath = QDir::currentPath()+"/untitled.mdl";
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                  initialPath,
                                                  tr("model file (*.mdl)"));

  if (fileName.isEmpty()) return false;
    
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

  out<<text()<<endl;
  QApplication::restoreOverrideCursor(); 
  file.close();
  return true;
  
}

void CpWindow::replaceBGroupClicked(int bid){
  if(bid==0)mGroup->show();
  else{
    mGroup->setChecked(false);
    mGroup->hide();
  }
  setText();
}

void CpWindow::nameBGroupClicked(int bId){
  if(bId==0){
    mGroup->show();
    mGroup->setChecked(true);
    nameEdit->setInputMask("_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
  }else{
    mGroup->show();
    mGroup->setChecked(false);
    nameEdit->setInputMask(QString());
  }
  setText();
}

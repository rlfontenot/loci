#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>

#include "physicswindow.h"
#include "pages.h"
#include "getfile.h"









FloatEditDelegate::FloatEditDelegate(QObject *parent)
  : QItemDelegate(parent)
{
}

QWidget *FloatEditDelegate::createEditor(QWidget *parent,
                                         const QStyleOptionViewItem &/* option */,
                                         const QModelIndex &/* index */) const
{
  FloatEdit *editor = new FloatEdit(parent);
  return editor;
}

void FloatEditDelegate::setEditorData(QWidget *editor,
                                      const QModelIndex &index) const
{
  double value = index.model()->data(index, Qt::EditRole).toDouble();
  
     FloatEdit *floatEdit = static_cast<FloatEdit*>(editor);
     floatEdit->setValue(value);
}

void FloatEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                    const QModelIndex &index) const
{
  FloatEdit *floatEdit = static_cast<FloatEdit*>(editor);
    
  double value = floatEdit->value();
  
  model->setData(index, value, Qt::EditRole);
}

void FloatEditDelegate::updateEditorGeometry(QWidget *editor,
                                             const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
  editor->setGeometry(option.rect);
}

CpWindow ::CpWindow(QDomElement& elem, QWidget* parent):QWidget(parent), myelem(elem){
 
  QGroupBox *nameGroup= new QGroupBox(tr("please specify name of the gas"));
  QHBoxLayout *nameLayout = new QHBoxLayout;
  nameEdit = new QLineEdit(tr("_gas"));
  nameLayout->addWidget(nameEdit);
  nameGroup->setLayout(nameLayout);
  
  QGroupBox  *replaceGroup = new QGroupBox(tr("please specify replacement or augment"));
  QHBoxLayout *replaceLayout = new QHBoxLayout;
  replaceButton = new QRadioButton(tr("replace"));
  augmentButton = new QRadioButton(tr("augment"));
  replaceButton->setChecked(true);                               
  replaceLayout->addWidget(replaceButton);
  replaceLayout->addWidget(augmentButton);
  replaceGroup->setLayout(replaceLayout);
  connect(replaceButton, SIGNAL(toggled(bool)), this, SLOT(replaceButtonToggled(bool)));
  
  mGroup = new QGroupBox(tr("please specify molecular mass"));
  QHBoxLayout *mLayout = new QHBoxLayout;
  mEdit = new FloatEdit;
  //  mEdit->setMinimum(0.0);
  mLayout->addWidget(mEdit);
  mGroup->setLayout(mLayout);
  
  QGroupBox *equationGroup = new QGroupBox(tr("please specify equation type "));
  QHBoxLayout *equationLayout = new QHBoxLayout;
  shomateButton = new QRadioButton(tr("shomate"));
  polyButton = new QRadioButton(tr("polynomial"));
  shomateButton->setChecked(true);                               
  equationLayout->addWidget(shomateButton);
  equationLayout->addWidget(polyButton);
  equationGroup->setLayout(equationLayout);
  
  QGroupBox *nGroup = new QGroupBox(tr("please specify number of  temperature intervals"));
  QHBoxLayout *nLayout = new QHBoxLayout;
  nEdit = new IntEdit;
  nEdit->setRange(1, 100);
  nEdit->setValue(3);
  
  nLayout->addWidget(nEdit);
  nGroup->setLayout(nLayout);
  connect(nEdit, SIGNAL(valueChanged(int)), this, SLOT(setNumberOfIntervals(int))); 

  
  

  
  
  model = new  QStandardItemModel();
  //  for(int i = 0; i <=n; i++){
  //     QString str;
  //     str.arg(i);
  //     model->setHeaderData(i, Qt::Horizontal, str);
  //   }
  model->setColumnCount(6);
  model->setHeaderData(0, Qt::Horizontal, QObject::tr("Temperature"));
  model->setHeaderData(1, Qt::Horizontal, QObject::tr("      A      "));
  model->setHeaderData(2, Qt::Horizontal, QObject::tr("      B      "));
  model->setHeaderData(3, Qt::Horizontal, QObject::tr("      C      "));
  model->setHeaderData(4, Qt::Horizontal, QObject::tr("      D      "));
  model->setHeaderData(5, Qt::Horizontal, QObject::tr("      E      "));
  
  QItemSelectionModel *selections = new QItemSelectionModel(model);
  FloatEditDelegate* delegate = new FloatEditDelegate(this);
  
  tableView = new QTableView(this);
  tableView->setModel(model);
  tableView->setSelectionModel(selections);
  tableView->setSelectionMode(QAbstractItemView::SingleSelection);
  tableView->setItemDelegate(delegate);
  // tableView->resizeRowsToContents();
  tableView->resizeColumnsToContents();
  
  numberOfIntervals =3;
  model->setRowCount(4);

  QPushButton* saveButton = new QPushButton(tr("          save the module        "));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  
  
  QVBoxLayout* mainLayout = new QVBoxLayout;
  
  mainLayout->addWidget(nameGroup);
  mainLayout->addWidget(replaceGroup);
  mainLayout->addWidget(mGroup);
  mainLayout->addWidget(equationGroup);
  mainLayout->addWidget(nGroup);
  mainLayout->addWidget(tableView);
  mainLayout->addWidget(saveButton);
  setLayout(mainLayout);
  //    tableView->hide();
}

void CpWindow::replaceButtonToggled(bool checked){
  if(checked)mGroup->show();
  else mGroup->hide();
}
  
bool CpWindow::save(){
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                  "/home/qxue/chemdemo/untitled.mdl",
                                                  tr("module file (*.mdl)"));

  if (fileName.isEmpty())
    {
      QMessageBox::warning(this, tr("CPWindow"), tr("Please specify filename")); 
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
  out << nameEdit->text();
     QString tab=QString(nameEdit->text().size()+8, ' ');
     QString firstTab=QString(nameEdit->text().size()+3, ' ');
     if(replaceButton->isChecked())out<<" = <";
     else out <<" : <";
     if(replaceButton->isChecked()){
       out<<"m=";
       out<<mEdit->value();
       out<<", href=55749, sref=130751, Tref=300, Pref=101325.0, mf=1,"<<endl;
       out<<firstTab+ "cp=[ ";
     }else{
       out<< "cp=[ ";
     }
     QString prefix;
     if(shomateButton->isChecked())prefix =tab+"shomate";
     else prefix=tab+"poly";
     
     for(int i =0; i < (nEdit->value()-1); i++){
       if(i==0) out<<(model->index(i, 0)).data(Qt::EditRole).toString()+","<<endl;
       else out<<tab+(model->index(i, 0)).data(Qt::EditRole).toString()+","<<endl;
       out<<prefix+"(";
       for(int j=1; j<5; j++){
         out<<(model->index(i, j)).data(Qt::EditRole).toString()+",";
       }
       out<<(model->index(i, 5)).data(Qt::EditRole).toString()+"),"<<endl;
     }
     out<<tab+(model->index(nEdit->value()-1, 0)).data(Qt::EditRole).toString()+"]>;"<<endl;
     
     QApplication::restoreOverrideCursor();
     file.close();
     myelem.setAttribute("current", fileName);
     return true;
     
}
void CpWindow::setNumberOfIntervals(int n){
  numberOfIntervals = n;
  model->setRowCount(n+1);
  tableView->setModel(model);
  tableView->show();
}





PhysicsWindow::PhysicsWindow(QDomElement& theelem, QDomElement& theroot, QWidget* parent)
  :GeneralWindow(theelem, theroot, parent)
{

  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+ tr(" has no child")
                         );
    return;
  }
  int count=0;   
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) {   
    if(elem.attribute("element")=="stateStack"){
      StackGroup2* stackGroup = new StackGroup2(elem, myroot, state, this);
      stacks << stackGroup;
      state << stackGroup->myState();
      connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(updateState(QString)));
      connect(this, SIGNAL(parentStateChanged(QString)), stackGroup, SLOT(parentStateChanged(QString)));
      connect(this, SIGNAL(stateUpdated(const QStringList&)), stackGroup, SLOT(setParentState(const QStringList&)));
       connect(this, SIGNAL(showStatus(const bool &)), stackGroup, SLOT(updateShowStatus(const bool &)));
      connect(stackGroup, SIGNAL(stateChanged(QString)), this, SLOT(checkStatus()));
      connect(stackGroup, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(stackGroup, SIGNAL(componentsChanged()), this, SIGNAL(componentsChanged()));
      mainLayout->addWidget(stackGroup);
    }else{
      QMessageBox::warning(window(), myelem.tagName(),
                           tr("don't know how to handle it yet ")+ myelem.attribute("element") 
                           );
    }
  }
    emit stateUpdated(state);
    setLayout(mainLayout);
    setWindowTitle(myelem.attribute("title"));
    checkStatus();
}

void PhysicsWindow::updateState(QString stat){
  QStringList list1=stat.split("=");
  bool Found = false;
  for(int i = 0; i < state.size(); i++){
    if(state[i].contains(list1[0])){
      if(state[i]!=stat)emit parentStateChanged(stat);
      state[i] = stat;
      
      Found = true;
    }
  }
  if(Found){
    QString tmp = state.join(";");
    if(tmp != myelem.attribute("state")){
      myelem.setAttribute("state", tmp);
      if(!tmp.isEmpty())emit stateChanged();
    }
  }
  if(!Found){
    QMessageBox::warning(window(), myelem.tagName(),
                          tr("illegal state ")+ stat 
                         );
    return;
  }
  
  
}







 void PhysicsWindow::checkStatus(){
  int count = 0;
  int count_done = 0;
  QString text;
  for( QDomElement elt = myelem.firstChildElement(); !elt.isNull();
       elt = elt.nextSiblingElement()){
    if(elt.attribute("status")=="done")count_done++;
    count++;
    text += elt.attribute("currentText");
  }
  


  
  
  myelem.setAttribute("currentText", text);
  if(count==count_done) myelem.setAttribute("status", "done");
  else myelem.setAttribute("status", QString("%1  out of ").arg(count_done)+QString(" %1 finished").arg(count));


//   QString tmp = state.join(";");
//   myelem.setAttribute("state", tmp);
//    if(!tmp.isEmpty())emit stateChanged();
  emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
 }



void StackGroup2::add(){
  
   QDomElement elm = myelem.firstChildElement();
   if(elm.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                              tr("stack elememt has no child")
                          );
     return;
   }
   for(int i = 0; i < current; i++)elm = elm.nextSiblingElement();
    if(elm.isNull()){
     QMessageBox::warning(window(), tr(".xml"),
                              tr("stack elememt has not enough children")
                          );
     return;
    }
  
    QDomElement elt_page = myroot.firstChildElement("models");
  
    if(elt_page.isNull()){
      QMessageBox::warning(window(), tr(".xml"),
                           tr("can not find element 'models' ")
                           );
      return;
    }
   
  QDomElement elem= elt_page.firstChildElement(elm.attribute("define"));
  
  if(elem.isNull()){
    QMessageBox::warning(window(), tr(".xml"),
                         tr("can not find the model ")+elm.attribute("define")
                         );
    return;
  }

     
  if(elem.attribute("element")=="all"){
    
      AllVBWindow* window= new AllVBWindow(elem, myroot);
      
      window->show();
         
    }else if(elem.attribute("element")=="CpWindow"){
      CpWindow* window=new CpWindow(elem);
      window->show();
    
    }
 
}
        




        


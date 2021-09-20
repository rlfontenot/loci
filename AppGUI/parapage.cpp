
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
#include <QListWidget>
#include <QStackedWidget>
#include <QSignalMapper>
#include <QScrollArea>
#include <QToolBar>


#include "parapage.h"
#include "pages.h"

#define PI 3.14159265358979323846264338327950
void Shape::reset(){
  switch(tp){
  case  SPHERE:
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 1;
    
    break;
  
  case CONE:

   
   
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 0;
    para[4] = 0;
    para[5] = 1;
    para[6] = 1;
    para[7] = 0;
    
    break;
   
  case CYLINDER:
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 0;
    para[4] = 0;
    para[5] = 1;
    para[6] = 1;
    break;
    
  case BOX:
   
     para[0] = 0;
     para[1] = 0;
     para[2] = 0;
     para[3] = 1;
     para[4] = 1;
     para[5] = 1;
     break;
  default:
    para[0] = 0;
    para[1] = 0;
    para[2] = 0;
    para[3] = 0;
    para[4] = 0;
    para[5] = 1;
  }
}
 
ParaPage::ParaPage(Shape* s, QWidget* parent):QGroupBox(tr(""),parent),shape(s){
  if(shape==0)return;
  
  signalMapper = new QSignalMapper(this);
  QBoxLayout *mainLayout;
  
  switch(shape->tp){
  case SPHERE:
    {

      setTitle("parameters of sphere:");
      objs.resize(4);

      
      objs[0] = new LabeledDoubleSpBox(tr("x0"));
      objs[1] = new LabeledDoubleSpBox(tr("y0"));
      objs[2] = new LabeledDoubleSpBox(tr("z0"));
      objs[3] = new LabeledDoubleSpBox(tr("radius"));
      for(int i = 0; i < 4; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));


      QGroupBox *centerGroup = new QGroupBox("center", this);
      QVBoxLayout *hBoxLayout = new QVBoxLayout;
      hBoxLayout->addWidget(objs[0]);
      hBoxLayout->addWidget(objs[1]);
      hBoxLayout->addWidget(objs[2]);
      centerGroup->setLayout(hBoxLayout);

      QGroupBox *radiusGroup = new QGroupBox("radius");
      QVBoxLayout *hBoxLayout2 = new QVBoxLayout;
      hBoxLayout2->addWidget(objs[3]);
      radiusGroup->setLayout(hBoxLayout2);

      
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(centerGroup);
      mainLayout->addWidget(radiusGroup);
        mainLayout->addStretch(2);
      setLayout(mainLayout);
    }
    break;
  case CONE:
     {

      setTitle("parameters of cone:");
      objs.resize(8);

      QGroupBox *p1Group = new QGroupBox("p1", this);
      QVBoxLayout *p1BoxLayout = new QVBoxLayout;
      
      objs[0] = new LabeledDoubleSpBox(tr("x1"));
      objs[1] = new LabeledDoubleSpBox(tr("y1"));
      objs[2] = new LabeledDoubleSpBox(tr("z1"));
      p1BoxLayout->addWidget(objs[0]);
      p1BoxLayout->addWidget(objs[1]);
      p1BoxLayout->addWidget(objs[2]);
      p1Group->setLayout(p1BoxLayout);

       QGroupBox *p2Group = new QGroupBox("p2", this);
      QVBoxLayout *p2BoxLayout = new QVBoxLayout;
      objs[3] = new LabeledDoubleSpBox(tr("x2"));
      objs[4] = new LabeledDoubleSpBox(tr("y2"));
      objs[5] = new LabeledDoubleSpBox(tr("z2"));
      p2BoxLayout->addWidget(objs[3]);
      p2BoxLayout->addWidget(objs[4]);
      p2BoxLayout->addWidget(objs[5]);
      p2Group->setLayout(p2BoxLayout);
      
      objs[6] = new LabeledDoubleSpBox(tr("r1"));
      objs[7] = new LabeledDoubleSpBox(tr("r2"));
      
      
      
      for(int i = 0; i < 8; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));
      
           
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(p1Group);
      mainLayout->addWidget(p2Group);
      mainLayout->addWidget(objs[6]);
      mainLayout->addWidget(objs[7]);
      mainLayout->addStretch(2);
      setLayout(mainLayout);
     }
     break;
  case CYLINDER:
    {
      setTitle("parameters of cylinder:");
      objs.resize(7);
      QGroupBox *p1Group = new QGroupBox("p1", this);
      QVBoxLayout *p1BoxLayout = new QVBoxLayout;
      
      objs[0] = new LabeledDoubleSpBox(tr("x1"));
      objs[1] = new LabeledDoubleSpBox(tr("y1"));
      objs[2] = new LabeledDoubleSpBox(tr("z1"));
      p1BoxLayout->addWidget(objs[0]);
      p1BoxLayout->addWidget(objs[1]);
      p1BoxLayout->addWidget(objs[2]);
      p1Group->setLayout(p1BoxLayout);
      
       QGroupBox *p2Group = new QGroupBox("p2", this);
       QVBoxLayout *p2BoxLayout = new QVBoxLayout;
       objs[3] = new LabeledDoubleSpBox(tr("x2"));
       objs[4] = new LabeledDoubleSpBox(tr("y2"));
       objs[5] = new LabeledDoubleSpBox(tr("z2"));
       p2BoxLayout->addWidget(objs[3]);
       p2BoxLayout->addWidget(objs[4]);
       p2BoxLayout->addWidget(objs[5]);
       p2Group->setLayout(p2BoxLayout);
       
       objs[6] = new LabeledDoubleSpBox(tr("r1"));
      
      
      
      
      for(int i = 0; i < 7; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));
      
           
      

      
       mainLayout = new QVBoxLayout;
       mainLayout->addWidget(p1Group);
       mainLayout->addWidget(p2Group);
       mainLayout->addWidget(objs[6]);
       mainLayout->addStretch(2);
       setLayout(mainLayout);
    }
    break;
  case BOX:

    {
      setTitle("parameters of box:");
      objs.resize(6);
      QGroupBox *p1Group = new QGroupBox("p1", this);
      QVBoxLayout *p1BoxLayout = new QVBoxLayout;
      
      objs[0] = new LabeledDoubleSpBox(tr("x1"));
      objs[1] = new LabeledDoubleSpBox(tr("y1"));
      objs[2] = new LabeledDoubleSpBox(tr("z1"));
      p1BoxLayout->addWidget(objs[0]);
      p1BoxLayout->addWidget(objs[1]);
      p1BoxLayout->addWidget(objs[2]);
      p1Group->setLayout(p1BoxLayout);
      
      QGroupBox *p2Group = new QGroupBox("p2", this);
      QVBoxLayout *p2BoxLayout = new QVBoxLayout;
      objs[3] = new LabeledDoubleSpBox(tr("x2"));
      objs[4] = new LabeledDoubleSpBox(tr("y2"));
      objs[5] = new LabeledDoubleSpBox(tr("z2"));
      p2BoxLayout->addWidget(objs[3]);
       p2BoxLayout->addWidget(objs[4]);
       p2BoxLayout->addWidget(objs[5]);
       p2Group->setLayout(p2BoxLayout);
       
      
       
      for(int i = 0; i < 6; i++){
        objs[i]->setValue(shape->para[i]);
        signalMapper->setMapping(objs[i], i);
        connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
      }
      connect(signalMapper, SIGNAL(mapped(int)),
              this, SLOT(setValue(int)));
      
            
      mainLayout = new QVBoxLayout;
      mainLayout->addWidget(p1Group);
      mainLayout->addWidget(p2Group);
      mainLayout->addStretch(2);
      setLayout(mainLayout);
    }
    
    break;
 

 case LEFTPLANE:
 {
   setTitle("parameters of leftPlane:");
   objs.resize(6);
   QGroupBox *p1Group = new QGroupBox("point", this);
   QVBoxLayout *p1BoxLayout = new QVBoxLayout;
   
   objs[0] = new LabeledDoubleSpBox(tr("x1"));
   objs[1] = new LabeledDoubleSpBox(tr("y1"));
   objs[2] = new LabeledDoubleSpBox(tr("z1"));
   p1BoxLayout->addWidget(objs[0]);
   p1BoxLayout->addWidget(objs[1]);
   p1BoxLayout->addWidget(objs[2]);
   p1Group->setLayout(p1BoxLayout);
      
   QGroupBox *p2Group = new QGroupBox("normal", this);
   QVBoxLayout *p2BoxLayout = new QVBoxLayout;
   objs[3] = new LabeledDoubleSpBox(tr("x2"));
   objs[4] = new LabeledDoubleSpBox(tr("y2"));
   objs[5] = new LabeledDoubleSpBox(tr("z2"));
   p2BoxLayout->addWidget(objs[3]);
   p2BoxLayout->addWidget(objs[4]);
   p2BoxLayout->addWidget(objs[5]);
   p2Group->setLayout(p2BoxLayout);
   
   
   
   for(int i = 0; i < 6; i++){
     objs[i]->setValue(shape->para[i]);
     signalMapper->setMapping(objs[i], i);
     connect(objs[i], SIGNAL(valueChanged(double)),signalMapper, SLOT(map()));
   }
   connect(signalMapper, SIGNAL(mapped(int)),
           this, SLOT(setValue(int)));
   
   
   mainLayout = new QVBoxLayout;
   mainLayout->addWidget(p1Group);
   mainLayout->addWidget(p2Group);
   mainLayout->addStretch(2);
   setLayout(mainLayout);
 }
    
 break;

  default:
    break;
  }
  
}
void ParaPage::setValue(int i){
  if(i < 0 || i >= (int)objs.size() || i >= (int)(shape->para.size()))return;
  shape->para[i] = objs[i]->value();
   emit valueChanged();
}
void ParaPage::showValue(){
  for(unsigned int i = 0; i < objs.size(); i++){
    objs[i]->setValue(shape->para[i]);
  }
}
void ParaPage::reset(){
  shape->reset();
  showValue();
}

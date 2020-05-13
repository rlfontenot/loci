//////////////////////////////////////////////////////////
//  
//////////////////////////////////////////////////////////
#include <QVBoxLayout>
#include "transform.h"

Transform::Transform( QWidget *parent)
  :  QGroupBox(tr("transform:"),parent) 
{
  
 
  translate = new VectSpBox(tr("translate"));
  rotateCenter  = new VectSpBox(tr("rotateCenter"));
  rotateAngle = new VectSpBox(tr("rotateAngle"));
  scale = new  VectSpBox(tr("scale"));
  scale->setValue(positions3d(1.0, 1.0, 1.0));
  translate->setRange(-1e5, 1e5);
  rotateCenter->setRange(-1e5, 1e5);
  scale->setRange(-1e5, 1e5);
  rotateAngle->setRange(-180, 180);
  
   
  connect(translate, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(rotateCenter, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(rotateAngle, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));
  connect(scale, SIGNAL(valueChanged(const positions3d&)), this, SLOT(setInfo()));

  
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addWidget(translate);
  mainLayout->addWidget(rotateCenter);
  mainLayout->addWidget(rotateAngle);
  mainLayout->addWidget(scale);
  mainLayout->addStretch(1);
  setLayout(mainLayout);
}

void Transform::clear(){
 
 
  translate->setValue(positions3d(0.0, 0.0, 0.0));
  rotateCenter->setValue(positions3d(0.0, 0.0, 0.0));
  rotateAngle->setValue(positions3d(0.0, 0.0, 0.0));
  scale->setValue(positions3d(1.0, 1.0, 1.0));
}





TranCoef Transform::value(){
  TranCoef f;
  f.translate = translate->value();
  f.rotateCenter = rotateCenter->value();
  f.rotateAngle = rotateAngle->value();
  f.scale = scale->value();
  return f;
}

void Transform::setValue(positions3d* p){
  translate->setValue(p[0]);
  rotateCenter->setValue(p[1]);
  rotateAngle->setValue(p[2]);
  scale->setValue(p[3]);
}



void Transform::setInfo(){
 
 
   emit tcChanged();
  
}




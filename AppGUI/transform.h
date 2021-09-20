#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "pages.h"
struct IDMatrix2{
  int gridId;
  affineMapping matrix;
  IDMatrix2(int id, const affineMapping& m):gridId(id), matrix(m){}
};
  
class Transform: public QGroupBox
{
  Q_OBJECT
  public:
  Transform(QWidget *parent = 0);
  TranCoef value();
 
public slots:
  void setInfo();
  void clear();
  void setValue(positions3d* p);
signals:
  void tcChanged();
  
private:
  VectSpBox* translate;
  VectSpBox* rotateCenter;
  VectSpBox* rotateAngle;
  VectSpBox* scale;
};
#endif

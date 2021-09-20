#ifndef PARAPAGE_H
#define PARAPAGE_H
#include "pages.h"
class QGroupBox;
class QSignalMapper;

enum BASIC_SHAPES{SPHERE, CONE, CYLINDER, BOX, LEFTPLANE};
class Shape{
public:
  Shape(BASIC_SHAPES t, vector<double>& p):tp(t),para(p){};
  Shape(){};
  void reset();
public:
  BASIC_SHAPES tp;
  vector<double> para;
};


class ParaPage:public QGroupBox{
  Q_OBJECT
  public:
  ParaPage(Shape* s, QWidget* parent = 0);
                                         
public slots:
  void setValue(int);
  void showValue();
  void reset();
signals:
  void valueChanged();
public:
  
  vector<LabeledDoubleSpBox*> objs;
  Shape* shape;
  QSignalMapper* signalMapper;
};
#endif

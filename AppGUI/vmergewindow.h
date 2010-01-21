#ifndef VMERGEWINDOW_H
#define VMERGEWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QModelIndex>
#include <QGroupBox>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QTableView>
#include <QList>
#include <QHBoxLayout>
#include <utility>
#include "pages.h"
#include "grid.h"

class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QCheckBox;
class QButtonGroup;
using namespace std;

//typedef affineMapping affineMapping2;


struct affineMapping2 {
  double M[4][4] ;
  double determinant() { return
      (M[0][0]*M[1][1]*M[2][2]+
       M[1][0]*M[2][1]*M[0][2]+
       M[2][0]*M[0][1]*M[1][2]) -
      (M[0][0]*M[2][1]*M[1][2]+
       M[1][0]*M[0][1]*M[2][2]+
       M[2][0]*M[1][1]*M[0][2]) ;
  }

  bool leftHanded() {return (determinant() < 0) ; }
  affineMapping2() {
    for(int i=0;i<4;++i) {
      for(int j=0;j<4;++j)
        M[i][j] = 0 ;
      M[i][i] = 1 ;
    }
  }
  void Combine(affineMapping2 a) {
    double Mtmp[4][4] ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) {
        double mtmp = 0 ;
        for(int k=0;k<4;++k)
          mtmp += a.M[i][k]*M[k][j] ;
        Mtmp[i][j] = mtmp ;
      }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        M[i][j] = Mtmp[i][j] ;
  }
  void translate(vector3d<double> tv) {
    affineMapping2 tmp ;
    tmp.M[0][3] = tv.x ;
    tmp.M[1][3] = tv.y ;
    tmp.M[2][3] = tv.z ;
    Combine(tmp) ;
  }
  void scale(vector3d<double> tv) {
    affineMapping2 tmp ;
    tmp.M[0][0] = tv.x ;
    tmp.M[1][1] = tv.y ;
    tmp.M[2][2] = tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping2 tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  sth ;
    tmp.M[2][1] = -sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping2 tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][2] = -sth ;
    tmp.M[2][0] =  sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping2 tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][1] =  sth ;
    tmp.M[1][0] = -sth ;
    tmp.M[1][1] =  cth ;
    Combine(tmp) ;
  }
  vector3d<double> Map(vector3d<double> v) {
    double tmp[4] ;
    tmp[0] = v.x ;
    tmp[1] = v.y ;
    tmp[2] = v.z ;
    tmp[3] = 1. ;
    double res[4] ;
    for(int i=0;i<4;++i)
      res[i] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        res[i] += M[i][j]*tmp[j] ;
    vector3d<double> r(res[0],res[1],res[2]) ;
    return r ;
  }
   
} ;

struct TranCoef{
  // int gridId;
  vector3d<double> translate;
  vector3d<double> rotateAngle;
  vector3d<double> rotateCenter;
  vector3d<double> scale;
  TranCoef(const vector3d<double> &v1, const vector3d<double> &v2,
           const vector3d<double> &v4, const vector3d<double> &v3):
    translate(v1), rotateAngle(v2), rotateCenter(v4), scale(v3){}

  TranCoef(){
    
    translate = vector3d<double>(0.0, 0.0, 0.0);
    rotateAngle = vector3d<double>(0.0, 0.0, 0.0);
    rotateCenter = vector3d<double>(0.0, 0.0, 0.0);
    scale = vector3d<double>(1.0, 1.0, 1.0);
  }
  
  
};

struct IDMatrix{
  int gridId;
  affineMapping2 matrix;
  IDMatrix(int id, const affineMapping2& m):gridId(id), matrix(m){}
};

struct IDColor{
  int gridId;
  int boundId;
  QColor color;
  IDColor(int id1, int id2, QColor value):gridId(id1), boundId(id2), color(value){}
  
};

struct IDVisibility{
  int gridId;
  int boundId;
  bool show;
  IDVisibility(int id1, int id2, bool value):gridId(id1), boundId(id2), show(value){}
  
};


class VMOption: public QGroupBox{
  Q_OBJECT
public:
  VMOption(int gridId, const QString &gridname, const QStringList &bcnames, QWidget *parent=0 );

  signals:
  void tcChanged(const IDMatrix&);
  void  setCurrentColor(const IDColor&);
  void setCurrentVisibility(const IDVisibility&);
  void setCurrent(QModelIndex);
  
public slots:
//void addGrid();
  void setCurrentObj(QModelIndex);
 
  void setInfo();
  void showBoundary(QModelIndex, QModelIndex);
  void accept();
  void cancel();
  void previous();
  void next();
  void clear();
public:
   QString currentText();
private:
  affineMapping2 currentM();
private:
  //translate
  FloatEdit* xEditor1;
  FloatEdit* yEditor1;
  FloatEdit* zEditor1;
  //rotate angle
  FloatEdit* xEditor2;
  FloatEdit* yEditor2;
  FloatEdit* zEditor2;
  //scale
  FloatEdit* xEditor3;
  FloatEdit* yEditor3;
  FloatEdit* zEditor3;
  //roate center
  FloatEdit* xEditor4;
  FloatEdit* yEditor4;
  FloatEdit* zEditor4;
  //tag
  QLineEdit* tagEditor;

  QStandardItemModel* modBoundaries;  // Boundary condition model
  QTableView* boundaryView;  // Use for boundary Select
  

  QString gridName;
  int gridId;
  vector<TranCoef> tc;
  int currentCoef;
  QString tag;
  QList<pair<QString, QString> > bdnames;
  // int currentBd;
};


class VMergeWindow : public QWidget
{
  Q_OBJECT
  
public:
  VMergeWindow(QWidget* parent = 0);
 
 public slots:
 void vmClicked();
  void loadGridClicked();
  void gridLoaded(const QStringList &);
  void changePage(QListWidgetItem *, QListWidgetItem *);
 
  void clearAll();
  
  signals:
  void currentGridChanged(int);
  void loadGrid(QString);//add a grid
  void getGrid(QString);//load in merged grid
  void tcChanged(const IDMatrix&);
  void  setCurrentColor(const IDColor&);
  void setCurrentVisibility(const IDVisibility&);
    void done();
private:
  void clear();
private:
  int currentRow;
  // VMOption* vmOptions;
  vector<vector<TranCoef> > transcoef;

  
  QListWidget *typesWidget;
  QStackedWidget *pagesWidget;
 
};

#endif

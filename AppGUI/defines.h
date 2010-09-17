#ifndef DEFINES_H
#define DEFINES_H

#include <map>
#include <list>
#include <vector>
#include <string>
#include <QString>
#include <QColor>
#include <iostream>
#include <math.h>
#include <utility>
#include <QDomDocument>
class QTreeWidgetItem;
using namespace std;


static const QColor default_color[] = {Qt::red, Qt::green, Qt::blue, Qt::cyan, Qt::magenta, Qt::yellow,
                                     Qt::darkRed, Qt::darkGreen, Qt::darkBlue,
                                     Qt::darkCyan, Qt::darkMagenta, Qt::darkYellow};



// Data structure passed to load dialog
struct LoadInfo
{
  //  bool load_boundary;
  //bool load_grid;
  //bool load_scalar;

  QString casename;
  QString iteration;
  QString variable;
  QString directory;
  LoadInfo():casename(""),iteration(""),variable(""),directory(""){};
};

struct edges {
  size_t l,r ;
  edges(size_t li,size_t ri) : l(li),r(ri) {}
  edges() {}
 
};

inline bool operator==(const edges &e1, const edges &e2)
{
 return (e1.l == e2.l && e1.r == e2.r) || (e1.l == e2.r && e1.r == e2.l) ;
}

inline bool operator<( const edges &e1,  const edges &e2)
{
  if (e1.l != e2.l)
    return e1.l < e2.l;
  else
    return e1.r < e2.r;
}

struct triangles {
  int t1,t2,t3 ;
  triangles(int ti1,int ti2, int ti3) :
    t1(ti1),t2(ti2),t3(ti3) {}
  triangles() {}
} ;

struct positions {
  double x,y ;
  positions(double xi,double yi): x(xi),y(yi) {}
  positions() {}
} ;


 template <class T> 
    struct vector3d {
      T x,y,z ;
      vector3d() {} 
      vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
      vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
    {
      s >> v.x >> v.y >> v.z ;
      return s ;
    }

  template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
  }

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
    return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                       v1.z*ra2[0]-v1.x*ra2[2],
                       v1.x*ra2[1]-v1.y*ra2[0]) ;
  }
  template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
    return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                       ra1[2]*v2.x-ra1[0]*v2.z,
                       ra1[0]*v2.y-ra1[1]*v2.x) ;
  }

template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, float val) {
  target.x *= val ;
  target.y *= val ;
  target.z *= val ;
  return target ;
}

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    target.z += val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    target.z -= val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(float r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(long double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

typedef  vector3d<double> positions3d;



inline positions operator+(const positions &p1, const positions &p2)
{ return positions(p1.x+p2.x,p1.y+p2.y) ; }

inline positions operator-(const positions &p1, const positions &p2)
{ return positions(p1.x-p2.x,p1.y-p2.y) ; }

inline positions operator*(const positions &p, double s)
{ return positions(p.x*s,p.y*s) ; }

inline positions operator*(double s,const positions &p)
{ return positions(p.x*s,p.y*s) ; }

inline positions operator/(const positions &p, double s)
{ return positions(p.x/s,p.y/s) ; }

inline positions operator/(double s,const positions &p)
{ return positions(p.x/s,p.y/s) ; }



//typedef affineMapping affineMapping2;
struct affineMapping {
  double M[4][4] ;
  //unit matrix
  affineMapping() {
    for(int i=0;i<4;++i) {
      for(int j=0;j<4;++j)
        M[i][j] = 0 ;
      M[i][i] = 1 ;
    }
  }
  
  void Combine(affineMapping a) {
    double Mtmp[4][4] ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) {
        double mtmp = 0 ;
        for(int k=0;k<4;++k)
          mtmp += M[i][k]*a.M[k][j] ;
        Mtmp[i][j] = mtmp ;
      }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        M[i][j] = Mtmp[i][j] ;
  }
  void translate(positions3d tv) {
    affineMapping tmp ;
    tmp.M[0][1] = -tv.x ;
    tmp.M[0][2] = -tv.y ;
    tmp.M[0][3] = -tv.z ;
    Combine(tmp) ;
  }
  void scale(positions3d tv) {
    affineMapping tmp ;
    tmp.M[1][1] = 1.0/tv.x ;
    tmp.M[2][2] = 1.0/tv.y ;
    tmp.M[3][3] = 1.0/tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[2][2] =  cth ;
    tmp.M[2][3] =  -sth ;
    tmp.M[3][2] = sth ;
    tmp.M[3][3] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][3] =  sth ;
    tmp.M[3][1] =  -sth ;
    tmp.M[3][3] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  -sth ;
    tmp.M[2][1] = sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  positions3d MapNode(positions3d v) {
    double tmp[4] ;
    tmp[0] = 1 ;
    tmp[1] = v.x ;
    tmp[2] = v.y ;
    tmp[3] = v.z ;
    double res[4] ;
    for(int i=0;i<4;++i)
      res[i] = 0 ;
    for(int j=0;j<4;++j)
      for(int i=0;i<4;++i)
        res[j] += M[i][j]*tmp[i] ;
    positions3d r(res[1],res[2],res[3]) ;
    return r ;
  }
   
} ;

template<class T, int S>
class Array
{
public:
  Array() {}

  T &operator[](size_t indx) { return array[indx]; }
  const T &operator[](size_t indx) const { return array[indx] ; }
  //private:
  T array[S];
};

template <class T, int S>
bool operator>( const Array<T,S> &arr1,  const Array<T,S> &arr2) {
  for (int i = 0; i < S; ++i) {
    if (arr1[i] == arr2[i]) continue;
    if (arr1[i] > arr2[i]) return true;
    break;
  }
  return false;
}


template <class T, int S>
bool operator<( const Array<T,S> &arr1, const Array<T,S> &arr2) {
  for (int i = 0; i < S; ++i) {
    if (arr1[i] == arr2[i]) continue;
    if (arr1[i] < arr2[i]) return true;
    break;
  }
  return false;
}

template <class T, int S>
bool operator==(const Array<T,S> &arr1, const Array<T,S> &arr2) {
  for (int i = 0; i < S; ++i) {
    if (arr1[i] == arr2[i]) continue;
    if (arr1[i] != arr2[i]) return false;
    break;
  }
  return true;
}




//the definition of affineMapping2 is copied from vogmerge.cc, should always
//be consistent with vogmerge.cc
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
  void Combine(affineMapping a) {
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
  void translate(positions3d tv) {
    affineMapping tmp ;
    tmp.M[0][3] = tv.x ;
    tmp.M[1][3] = tv.y ;
    tmp.M[2][3] = tv.z ;
    Combine(tmp) ;
  }
  void scale(positions3d tv) {
    affineMapping tmp ;
    tmp.M[0][0] = tv.x ;
    tmp.M[1][1] = tv.y ;
    tmp.M[2][2] = tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  sth ;
    tmp.M[2][1] = -sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][2] = -sth ;
    tmp.M[2][0] =  sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][1] =  sth ;
    tmp.M[1][0] = -sth ;
    tmp.M[1][1] =  cth ;
    Combine(tmp) ;
  }
  //mirror with a 3d plane with origin p and  normal n  
  void mirror(positions3d p, positions3d n){
    if(norm(n) > 1e-37) n = n/norm(n);   
    double k = dot(p, n);
    affineMapping tmp ;
    tmp.M[0][0] = 1-2.*n.x*n.x;
    tmp.M[1][1] = 1-2.*n.y*n.y;
    tmp.M[2][2] = 1-2.*n.z*n.z;
    
    tmp.M[0][1]  = tmp.M[1][0] = -2.*n.x*n.y;
    tmp.M[0][2] = tmp.M[2][0] = -2.*n.x*n.z;
    tmp.M[0][3] = 2.*n.x*k;
    tmp.M[3][0] = 0;

    tmp.M[1][2] = tmp.M[2][1] = -2.*n.y*n.z;
    tmp.M[1][3] = 2.*n.y*k;
    tmp.M[2][3] = 2.*n.z*k;
    Combine(tmp) ;
  }
    
    
  positions3d Map(positions3d v) {
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
    positions3d r(res[0],res[1],res[2]) ;
    return r ;
  }
   
} ;


struct TranCoef{  //in fvmadapt and other places
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

struct TranCoef2{  //used in vogmerge
  vector3d<double> translate;
  vector3d<double> rotateAngle;
  vector3d<double> rotateCenter;
  vector3d<double> scale;
  int checkedId;
  vector3d<double> mirrorOrigin;
  vector3d<double> mirrorNormal;
  
  TranCoef2(const vector3d<double> &v1, const vector3d<double> &v2,
            const vector3d<double> &v4, const vector3d<double> &v3,
            int id, const vector3d<double> &v5, const vector3d<double> &v6):
    translate(v1), rotateAngle(v2), rotateCenter(v4), scale(v3),
    checkedId(id), mirrorOrigin(v5), mirrorNormal(v6){}
  
  TranCoef2(){
    
    translate = vector3d<double>(0.0, 0.0, 0.0);
    rotateAngle = vector3d<double>(0.0, 0.0, 0.0);
    rotateCenter = vector3d<double>(0.0, 0.0, 0.0);
    scale = vector3d<double>(1.0, 1.0, 1.0);
    checkedId = 0;
    mirrorOrigin = vector3d<double>(0.0, 0.0, 0.0);
    mirrorNormal = vector3d<double>(1.0, 0.0, 0.0);
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

struct IDOnly{
  int gridId;
  int boundId;
  IDOnly(const pair<int, int>& id):gridId(id.first), boundId(id.second){}
};


struct IDVisibility{
  int gridId;
  int boundId;
  bool show;
  IDVisibility(int id1, int id2, bool value):gridId(id1), boundId(id2), show(value){}
  
};


QDomDocument tree2dom(const QTreeWidgetItem* root);



#endif

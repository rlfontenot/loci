#ifndef GRID_H
#define GRID_H

#include <map>
#include <list>
#include <vector>
#include <string>
#include <QString>
#include <QColor>
#include <iostream>
#include <math.h>

using std::map;
using std::list;
using std::vector;
using std::string;




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

struct segments {
  positions p1,p2 ;
  segments(positions p1i,positions p2i) : p1(p1i),p2(p2i) {}
} ;

struct contour_info {
  positions mingradpos ;
  double mingrad ;
  double value ;
  contour_info() {mingrad = 1e30; mingradpos = positions(0,0) ;}
  contour_info(double val)
    {mingrad = 1e30; mingradpos = positions(0,0) ; value = val ;}
} ;

struct rule {
  positions start ;
  double length ;
  bool xdir ;
  rule(const positions &p1,const positions &p2) {
    if(fabs(p1.x-p2.x) > fabs(p1.y-p2.y)) {
      xdir = true ;
      double y = 0.5*(p1.y+p2.y) ;
      start = positions(std::min(p1.x,p2.x),y) ;
      length = fabs(p1.x-p2.x) ;
    } else {
      xdir = false ;
      double x = 0.5*(p1.x+p2.x) ;
      start = positions(x,std::min(p1.y,p2.y)) ;
      length = fabs(p1.y-p2.y) ;
    }
  }
} ;


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




struct cutplane_info
{
  positions3d rotate;
  positions3d translate;

  // QString casename;
  // QString iteration;
  // QString variable;
  //QString directory;
};


  


struct grid {
  
  std::vector<positions> pos ;
  positions minpos,maxpos ;
  positions minview,maxview ;
  double size;

  bool has_values, valid ;
  double min_val,max_val,med_val ;
  std::vector<double> val ;
  
  std::vector<edges> edge_list ;
  unsigned int interior ;

  std::vector<triangles> triangle_list ;

  std::string filename ;

  double contour_spacing ;
  std::vector<segments> contour_curves ;
  std::vector<contour_info> contour_values ;

  std::vector<rule> rulers ;

  grid() { has_values = false ; filename = "none"; valid = false ; }
  void optimize_edge_list() ;
  void generate_contour_curves(int number = 0) ;

  void set_value_range() ;

  void input_2dgv(std::istream &in,bool read_values) ;
  void input_generalized(std::istream &in,bool read_values) ;
  void input_cobalt(std::string filename, bool read_values) ;
  void input(std::istream &in,bool read_values) ;
  void input(const char *filename, bool read_values = true) ;
  void cut(cutplane_info &info, const LoadInfo& load_info, const positions3d& center);
};

#endif

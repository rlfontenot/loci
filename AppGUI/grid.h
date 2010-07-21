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
#include "defines.h"
using std::map;
using std::list;
using std::vector;
using std::string;
using std::cerr;
using std::endl;



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



struct cutplane_info
{

  positions3d rotate;
  positions3d translate;
};


  


struct VolGrid {
  
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

  VolGrid() { has_values = false ; filename = "none"; valid = false ; }
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

#include "draw.h"
#include "debug.h"
#include "tools.h"
#include "pens.h"

#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

using std::min ;
using std::max ;
//using std::abs ;
using std::cerr ;
using std::cout ;
using std::endl ;

inline double sgn(double x) { return x>0?1.0:-1.0; }

const double EPSILON = 1e-20 ;

void drawfigure::scalefigure()
{
  const double minx = figure.minview.x ;
  const double miny = figure.minview.y ;
  const double maxx = figure.maxview.x ;
  const double maxy = figure.maxview.y ;
  double ticksize = currdev->getTickSize() ;
  double sgnx = sgn(currdev->getDrawBox().xright-currdev->getDrawBox().xleft) ;
  double sgny = sgn(currdev->getDrawBox().ytop-currdev->getDrawBox().ybottom) ;
  const double xleft   = currdev->getDrawBox().xleft   + sgnx*ticksize ;
  const double xright  = currdev->getDrawBox().xright  - sgnx*ticksize ;
  const double ytop    = currdev->getDrawBox().ytop    - sgny*ticksize ;
  const double adj = (currdev->isLabelDev())?8.0:1.0 ;
  const double ybottom = currdev->getDrawBox().ybottom + sgny*adj*ticksize ;
  titlepos = positions(xleft,ybottom-sgny*7.0*ticksize) ;
  legendpos = positions(xright-sgnx*30.0*ticksize,(ytop+ybottom)/2.0+sgny*double(MAXPENS/2)*ticksize*2.0) ;
  
  const double tscalex = (xright-xleft)/(maxx-minx+EPSILON) ;
  const double tscaley = (ytop-ybottom)/(maxy-miny+EPSILON) ;
  const double scale = min(fabs(tscalex),fabs(tscaley)) ;
  const double scalex = scale *sgn(tscalex) ;
  const double scaley = scale *sgn(tscaley) ;

  // Set scaleing factor
  figtransform.identity() ;
  figtransform.scale_transform(scalex,scaley) ;

  double basex = minx ;
  double basey = miny ;
  double offsetx = xleft ;
  double offsety = ybottom ;
  if(fabs(tscalex) < fabs(tscaley)) 
    offsety += (ytop-ybottom-fabs(maxy-miny)*scale*sgny)/2.0 ;
  else 
    offsetx += (xright-xleft-fabs(maxx-minx)*scale*sgnx)/2.0 ;

  // Pan transform to center of drawing space 
  figtransform.pan_transform(positions(offsetx-basex*scalex,
                                       offsety-basey*scaley)) ;
}

void drawfigure::drawgrid() {
  int nedges = figure.edge_list.size() ;

  SelectPen(0) ;
  for(int e=0;e<nedges;++e) {
    const edges &ed = figure.edge_list[e] ;
    FATAL(ed.l >= npos || ed.r >= npos) ;
    DrawLine(figure.pos[ed.l],figure.pos[ed.r]) ;
  }
}

void drawfigure::drawborder() {
  int nedges = figure.edge_list.size() ;

  SelectPen(0) ;
  for(int e=figure.interior;e<nedges;++e) {
    const edges &ed = figure.edge_list[e] ;
    FATAL(ed.l >= npos || ed.r >= npos) ;
    DrawLine(figure.pos[ed.l],figure.pos[ed.r]) ;
  }
}

void drawfigure::drawtriangles() {
  int ntri = figure.triangle_list.size() ;

  SelectPen(0) ;
  for(int t=0;t<ntri;++t) {
    const triangles &td = figure.triangle_list[t] ;
    FATAL(td.t1 >= npos || td.t2 >= npos  || td.t3>=npos) ;
    DrawLine(figure.pos[td.t1],figure.pos[td.t2]) ;
    DrawLine(figure.pos[td.t2],figure.pos[td.t3]) ;
    DrawLine(figure.pos[td.t3],figure.pos[td.t1]) ;
  }
}

void drawfigure::drawcontours()
{
  SelectPen(0) ;
  for(size_t i=0;i<figure.contour_curves.size();++i) {
    const segments &si=figure.contour_curves[i] ;
    DrawLine(si.p1,si.p2) ;
  }
  
  for(size_t i=0;i<figure.contour_values.size();++i) {
    if(figure.contour_values[i].mingrad < 100.0) {
      char buf[256] ;
      bzero(buf,256) ;
      snprintf(buf,255,"%s",fourSigDigs(figure.contour_values[i].value)) ;
      DrawLabel(buf,figure.contour_values[i].mingradpos,0,10,0) ;
    }
  }
}

void drawfigure::drawshading()
{
  for(int i=0;i<MAXPENS;++i) {
    if(figure.pnts[i].size() > 0) {
      SelectPen(i) ;
      vector<positions>::const_iterator pi = figure.pnts[i].begin();
      for(vector<int>::const_iterator vi=figure.pntindex[i].begin()+1;
          vi!=figure.pntindex[i].end();++vi) {
        int num_pnts = *vi - *(vi-1)  ;
        FillPoly(pi,num_pnts) ;
        pi = pi + num_pnts ;
      }
    }
  }    
}

void drawfigure::drawtitle() {
  char buf[512] ;
  char val[256] ;
  char date[21] ;
  time_t tloc ;

  time(&tloc) ;
  strncpy(buf,ctime(&tloc),4) ;
  strncpy(date,buf+4,12) ;
  strncpy(date+12,buf+19,5) ;
  date[17]='\0' ;

  if(!figure.valid) {
    snprintf(buf,512,"No grid available") ;
  } else if(figure.has_values) {
    snprintf(val,256,"%s",fourSigDigs(figure.contour_spacing)) ;
    snprintf(buf,512,"%s contour spacing = %s [Function: max=%4g min=%4g] Date: %s",
            figure.filename.c_str(),val,figure.max_val,figure.min_val,date) ;
  } else {
    snprintf(buf,512,"Grid filename: %s        Date: %s",figure.filename.c_str(),date) ;
  }
  currdev->DrawLabel(buf,titlepos,3,10,0) ;
}

void drawfigure::drawlegend()
{
  double ticksize = currdev->getTickSize()*2.0 ;
  positions rect[4] ;
  for(int i=1;i<MAXPENS;++i) {
    SelectPen(i) ;
    for(int j=0;j<4;++j) 
      rect[j] = positions(legendpos.x,legendpos.y+ticksize*double(i-1)) ;
    
    rect[1].x += ticksize*2.0 ;
    rect[2].x += ticksize*2.0 ;
    rect[2].y += ticksize ;
    rect[3].y += ticksize ;
    currdev->FillPoly(rect,4) ;
  }
  double maxv = figure.max_val ;
  double minv = figure.min_val ;
  for(size_t i=0;i<figure.contour_values.size();++i) {
    double v = figure.contour_values[i].value ;
    double s = (v-minv)/(maxv-minv) ;
    if(s<0.0 || s>1.0)
      continue ;
    positions p1(legendpos.x,legendpos.y+s*ticksize*double(MAXPENS-1.01)) ;
    positions p2(p1.x+ticksize*2.0,p1.y) ;
    positions p3(p2.x+ticksize,p1.y) ;
    SelectPen(int(1+s*double(MAXPENS-1.01))) ;
    currdev->DrawLine(p1,p2) ;
    char buf[256] ;
    snprintf(buf,256,"%s",fourSigDigs(figure.contour_values[i].value)) ;
    currdev->DrawLabel(buf,p3,0,10,0) ;
  }
  
}

void drawfigure::drawrulers()
{
  double ticksize = currdev->getTickSize() ;
  positions t1(0.,0.),t2(1.,0.),p1,p2 ;
  p1 = figtransform.apply_transform(t1) ;
  p2 = figtransform.apply_transform(t2) ;
  t2 = p1-p2 ;
  const double xscale = sqrt(t2.x*t2.x+t2.y*t2.y) ;
  t2 = positions(0.,1.) ;
  p2 = figtransform.apply_transform(t2) ;
  t2 = p1-p2 ;
  const double yscale = sqrt(t2.x*t2.x+t2.y*t2.y) ;
  const double num_graduations = 10. ;
  const double rangex = (figure.maxview.x-figure.minview.x)/num_graduations ;
  const double basex = pow(10.0,floor(log10(rangex))) ;
  const double incx  = floor(rangex/basex)*basex ;
  const double rangey = (figure.maxview.y-figure.minview.y)/num_graduations ;
  const double basey = pow(10.0,floor(log10(rangey))) ;
  const double incy  = floor(rangey/basey)*basey ;
  char buf[512] ;
  std::vector<rule>::const_iterator i ;
  for(i=figure.rulers.begin();i!=figure.rulers.end();++i) {
    p1 = i->start ;
    double len = i->length ;
    if(i->xdir) {
      p2 = p1 + len*positions(1.,0.) ;
      DrawLine(p1,p2) ;
      p1.x = floor(p1.x/incx)*incx ;
      if(p1.x < i->start.x)
        p1.x += incx ;
      while(p1.x <= p2.x) {
        if(fabs(p1.x*1000.0) < incx)
          p1.x = 0.;
        positions p3 = positions(p1.x,p1.y-ticksize/yscale) ;
        DrawLine(p1,p3) ;
        snprintf(buf,256,"%s",fourSigDigs(p1.x)) ;
        DrawLabel(buf,positions(p1.x,p1.y-4.*ticksize/yscale),5,10,0) ;
        p1.x += incx ;
      }
    } else {
      p2 = p1 + len*positions(0.,1.) ;
      DrawLine(p1,p2) ;
      p1.y = floor(p1.y/incy)*incy ;
      if(p1.y < i->start.y)
        p1.y += incy ;
      while(p1.y <= p2.y) {
        if(fabs(p1.y*1000.0) < incy)
          p1.y = 0. ;
        positions p3 = positions(p1.x+ticksize/xscale,p1.y) ;
        DrawLine(p1,p3) ;
        snprintf(buf,256,"%s",fourSigDigs(p1.y)) ;
        DrawLabel(buf,positions(p1.x+2.*ticksize/xscale,p1.y),5,10,0) ;
        p1.y += incy ;
      }
    }
  }
}

void drawfigure::draw()
{
  if(figure.show_shading)
    drawshading() ;
  drawborder() ;
  if(figure.show_grid)
    drawgrid() ;
  if(figure.show_contours)
    drawcontours() ;
  // Label the plot
  if(currdev->isLabelDev())
    drawtitle() ;
  if(figure.rulers.size() > 0)
    drawrulers() ;

  //  drawlegend() ;
}

void draw(const grid &drawgrid, grDevice &drawdev)
{
  drawfigure fig(drawgrid) ;
  drawdev.OpenDevice() ;
  fig.selectdev(drawdev) ;
  if(drawgrid.valid) {
    fig.scalefigure() ;
    fig.draw() ;
  } else {
    fig.drawtitle() ;
  }
  drawdev.CloseDevice() ;
}

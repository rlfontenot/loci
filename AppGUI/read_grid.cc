#include "cutplane.h"
#include "grid.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <strings.h>
#include <algorithm>


using std::ifstream ;
using std::ios ;
using std::cerr ;
using std::cin ;
using std::endl ;

using std::cout ;
using std::cerr ;
using std::endl ;
using std::istream ;
using std::vector ;
using std::list ;

using std::max ;
using std::min ;
using std::string ;




const double epsilon = 1e-30 ;

void VolGrid::cut(cutplane_info &info, const LoadInfo &load_info, const positions3d& center)
{
  CutPlane cutter;
  cutter.cut(info, load_info, center, this);
}

void VolGrid::optimize_edge_list() {
  const size_t num_nodes = pos.size() ;
  const size_t num_edges = edge_list.size() ;

  vector<int> node2edgec(num_nodes) ;
  for(size_t i=0;i<num_nodes;++i)
    node2edgec[i] = 0 ;

  for(size_t i=0;i<num_edges;++i) {
    node2edgec[edge_list[i].l]++ ;
    node2edgec[edge_list[i].r]++ ;
  }

  vector<int> node2edgeo(num_nodes+1) ;

    int count = 0 ;
  for(size_t i=0;i<num_nodes;++i) {
    node2edgeo[i] = count ;
    count += node2edgec[i] ;
  }
  node2edgeo[num_nodes] = count ;
  vector<size_t> node2edge(count) ;

  //  for(size_t i=0;i<num_edges;++i)
  //cout << node2edgec[i] << '\t' << node2edgeo[i] << endl;

  for(size_t i=0;i<num_edges;++i) {
    size_t l = edge_list[i].l ;
    size_t r = edge_list[i].r ;
    node2edge[node2edgeo[l]+(--node2edgec[l])] = i ;
    node2edge[node2edgeo[r]+(--node2edgec[r])] = i ;
  }

  vector<bool> edge_visited(num_edges) ;
  vector<edges> interior_edges,exterior_edges ;

  for(size_t i=0;i<num_nodes;++i) {
    for(int j=node2edgeo[i];j<node2edgeo[i+1];++j) {
      size_t je = node2edge[j] ;
      //      FATAL(je>=edge_visited.size()) ;
      if(!edge_visited[je]) {
        edge_visited[je] = true ;
        bool unique = true ;
        for(int k=j+1;k<node2edgeo[i+1];++k) 
          if(edge_list[je] == edge_list[node2edge[k]]) {
            interior_edges.push_back(edge_list[je]) ;
            edge_visited[node2edge[k]] = true ;
            unique = false ;
            break ;
          }
        if(unique)
          exterior_edges.push_back(edge_list[je]) ;
      }
    }
  }

  interior = interior_edges.size() ;
  if(!pos.empty()) {
    minpos = pos[0] ;
    maxpos = minpos ;
  } else {
    minpos = positions(0,0) ;
    maxpos = positions(0,0) ;
  }
  for(size_t i=0;i<exterior_edges.size();++i) {
    interior_edges.push_back(exterior_edges[i]) ;
    maxpos.x = max(maxpos.x,pos[exterior_edges[i].l].x) ;
    maxpos.x = max(maxpos.x,pos[exterior_edges[i].r].x) ;
    minpos.x = min(minpos.x,pos[exterior_edges[i].l].x) ;
    minpos.x = min(minpos.x,pos[exterior_edges[i].r].x) ;
    maxpos.y = max(maxpos.y,pos[exterior_edges[i].l].y) ;
    maxpos.y = max(maxpos.y,pos[exterior_edges[i].r].y) ;
    minpos.y = min(minpos.y,pos[exterior_edges[i].l].y) ;
    minpos.y = min(minpos.y,pos[exterior_edges[i].r].y) ; 
  }
  minview = minpos ;
  maxview = maxpos ;
  size = max(maxview.x - minview.x, maxview.y - minview.y);

  edge_list = interior_edges ;
}

// void VolGrid::input_2dgv(istream &in,bool read_values) {
//   bool multi = false ;
//   valid = true ;
//   has_values = read_values ;
//   if(parse::is_name(in)) {
//     string key = parse::get_name(in) ;
//     if(key != "multi") {
//       cerr << "improper key in 2dgv file '" << key << "'" << endl ;
//       valid = false ;
//       return ;
//     }
//     multi = true ;
//   }
//   do {

//     int ni=-1,nj=-1 ;
//     in >> ni ;
//     if(in.eof()) break ;
//     if(in.peek() == ',')
//       in.get() ;
//     in >> nj ;
//     if(ni<0 || nj<0 || ni>10000 || nj> 10000) {
//       cerr << "unreasonable block size (" << ni <<","<<nj<<")"<< endl;
//       ni=0;
//       nj=0;
//       valid = false ;
//       return ;
//     }
//     int pbase = pos.size() ;
    
//     for(int i=0;i<ni*nj;++i) {
//       double x=0,y=0;
//       in >> x ;
//       if(in.peek() == ',')
//         in.get() ;
//       in >> y ;
//       pos.push_back(positions(x,y)) ;
//     }

//     if(read_values) {
//       for(int i=0;i<ni*nj;++i) {
//         double v = 0 ;
//         in >> v ;
//         val.push_back(v) ;
//       }
      
//       WARN(vbase != pbase) ;
//     }
    

//     for(int i=0;i<ni-1;++i)
//       for(int j=0;j<nj-1;++j) {
//         int n1 = pbase + i + j*ni ;
//         int n2 = pbase + (i+1) + j*ni ;
//         int n3 = pbase + (i+1) + (j+1)*ni ;
//         int n4 = pbase + i + (j+1)*ni ;
//         positions center = (pos[n1]+pos[n2]+pos[n3]+pos[n4])*0.25 ;
//         double vc = (val[n1]+val[n2]+val[n3]+val[n4])*0.25 ;
//         int nc = pos.size() ;
//         pos.push_back(center) ;
//         if(read_values) {
//           val.push_back(vc) ;
//           triangle_list.push_back(triangles(n1,n2,nc)) ;
//           triangle_list.push_back(triangles(n2,n3,nc)) ;
//           triangle_list.push_back(triangles(n3,n4,nc)) ;
//           triangle_list.push_back(triangles(n4,n1,nc)) ;
//         }
//         edge_list.push_back(edges(n1,n2)) ;
//         edge_list.push_back(edges(n2,n3)) ;
//         edge_list.push_back(edges(n3,n4)) ;
//         edge_list.push_back(edges(n4,n1)) ;
//       }
//   } while(multi) ;
// }

// void VolGrid::input_generalized(std::istream &in,bool read_values) {
//   valid = true ;
//   has_values = read_values ;
//   if(parse::is_name(in)) {
//     string key = parse::get_name(in) ;
//     if(key != "general") {
//       cerr << "improper key in 2dgv file '" << key << "'" << endl ;
//       valid = false ;
//       return ;
//     }
//   } else {
//     cerr << "confused in input_generalized method" << endl ;
//     return ;
//   }
//   size_t nnodes, nodes_start ;
//   in >> nnodes >> nodes_start ;
//   pos.reserve(nnodes) ;
//   for(size_t i=0;i<nnodes;++i) {
//     double x = 0, y = 0 ;
//     in >> x >> y ;
//     pos.push_back(positions(x,y)) ;
//   }
//   int nfaces, faces_start, ncells, cells_start ;
//   in >> nfaces >> faces_start >> ncells >> cells_start ;
//   pos.reserve(nnodes+ncells) ;

//   edge_list.reserve(nfaces) ;
//   vector<edges> bedges ;
//   vector<edges> bf2c ;
//   vector<edges> f2c ;
//   f2c.reserve(nfaces) ;
//   bf2c.reserve(nfaces) ;
//   bedges.reserve(nfaces) ;
//   for(int i=0;i<nfaces;++i) {
//     int n1,n2,c1,c2 ;
//     in >> n1 >> n2 >> c1 >> c2 ;
//     if(c1<0)
//       std::swap(c1,c2) ;
//     c1 = c1-cells_start+nnodes ;
//     c2 = c2<0?c2:c2-cells_start+nnodes ;
//     n1 -= nodes_start ;
//     n2 -= nodes_start ;
//     if(c2 > 0) {
//       f2c.push_back(edges(c1,c2)) ;
//       edge_list.push_back(edges(n1,n2)) ;
//     } else {
//       bf2c.push_back(edges(c1,c2)) ;
//       bedges.push_back(edges(n1,n2)) ;
//     }
//   }
//   interior = edge_list.size() ;

//   if(!pos.empty()) {
//     minpos = pos[0] ;
//     maxpos = minpos ;
//   } else {
//     minpos = positions(0,0) ;
//     maxpos = positions(0,0) ;
//   }

//   for(size_t i=0;i<bedges.size();++i) {
//     maxpos.x = max(maxpos.x,pos[bedges[i].l].x) ;
//     maxpos.x = max(maxpos.x,pos[bedges[i].r].x) ;
//     minpos.x = min(minpos.x,pos[bedges[i].l].x) ;
//     minpos.x = min(minpos.x,pos[bedges[i].r].x) ;
//     maxpos.y = max(maxpos.y,pos[bedges[i].l].y) ;
//     maxpos.y = max(maxpos.y,pos[bedges[i].r].y) ;
//     minpos.y = min(minpos.y,pos[bedges[i].l].y) ;
//     minpos.y = min(minpos.y,pos[bedges[i].r].y) ; 

//     edge_list.push_back(bedges[i]) ;
//     f2c.push_back(bf2c[i]) ;
//   }
//   minview = minpos ;
//   maxview = maxpos ;
//   size = max(maxview.x - minview.x, maxview.y - minview.y);

//   val.reserve(nnodes+ncells) ;
//   if(parse::is_name(in)) {
//     string key = parse::get_name(in) ;
//     if(key != "data") {
//       cerr << "improper data key in 2dgv file '" << key << "'" << endl ;
//       has_values = false ;
//       return ;
//     }
//     int nvec=1, vecnum =1;
//     if(parse::is_int(in))
//       nvec = parse::get_int(in) ;
//     if(parse::is_int(in))
//       vecnum = parse::get_int(in) ;
//     if(!parse::is_string(in)) {
//       cerr << "file name not specified on data statement" << endl ;
//       has_values = false ;
//       return ;
//     }
//     string dname = parse::get_string(in) ;
//     ifstream data(dname.c_str(),ios::in) ;
//     if(data.fail()) {
//       cerr << "unable to open data file `" << dname << "'" << endl ;
//       has_values = false;
//       return ;
//     }
//     while(data.peek() == '#') {
//       int ch ;
//       while((ch = data.get()) != EOF && ch != '\n')
//         /*NULL STATEMENT */ ;
//     }
//     vector<double> v(nvec) ;
//     for(size_t i=0;i<nnodes;++i) {
//       if(data.fail() || data.eof()) {
//         has_values = false ;
//         cerr << "ran out of values" << endl ;
//         break ;
//       }
//       for(int j=0;j<nvec;++j)
//         data >> v[j] ;
//       val.push_back(v[vecnum]) ;
//     }
//   } else {
//     for(size_t i=0;i<nnodes;++i) {
//       if(in.fail() || in.eof()) {
//         has_values = false ;
//         break ;
//       }
//       double v ;
//       in >> v ;
//       val.push_back(v) ;
//     }
//   }

//   for(int i=0;i<ncells;++i) {
//     pos.push_back(positions(0,0)) ;
//     if(has_values)
//       val.push_back(double(0)) ;
//   }
//   vector<double> weight(pos.size()) ;
//   triangle_list.reserve(nfaces*2) ;
//   for(int i=0;i<nfaces;++i) {
//     int n1 = edge_list[i].l, n2 = edge_list[i].r ;
    
//     positions fc = .5*(pos[n1]+pos[n2]) ;
//     double vc = .5*(val[n1]+val[n2]) ;
//     double len = sqrt(pos[n1].x*pos[n1].x + pos[n2].y*pos[n2].y) ;
//     if(f2c[i].l >= nnodes && f2c[i].l< pos.size()) {
//       pos[f2c[i].l] = pos[f2c[i].l] + len*fc ;
//       if(has_values)
//         val[f2c[i].l] += len*vc ;
//       weight[f2c[i].l] += len ;
//       triangle_list.push_back(triangles(n1,n2,f2c[i].l)) ;
//     }
//     if(f2c[i].r >= nnodes && f2c[i].r< pos.size()) {
//       pos[f2c[i].r] = pos[f2c[i].r] + len*fc ;
//       weight[f2c[i].r] += len ;
//       if(has_values)
//         val[f2c[i].r] += len*vc ;
//       triangle_list.push_back(triangles(n1,n2,f2c[i].r)) ;
//     }
//   }
//   for(size_t i=nnodes;i<pos.size();++i) {
//     double rweight = 1./weight[i] ;
//     pos[i] = pos[i]*rweight ;
//     val[i] *= rweight ;
//   }
// }

// void VolGrid::input_cobalt(std::string filename,bool read_values) {
//   valid = true ;
//   has_values = read_values ;

//   ifstream in(filename.c_str(),ios::in) ;
//   if(in.fail()) {
//     cerr << "Open failed for " << filename << endl ;
//     has_values = false ;
//     return ;
//   }
//   int ndm, nzones, npatch ;
//   int npnts, nfaces, ncells, maxppf, maxfpc ;
//   in >> ndm >> nzones >> npatch;
//   in >> npnts >> nfaces >> ncells >> maxppf >> maxfpc ;

//   if(ndm != 2) {
//     cerr << "can only read two dimensional cobalt grids!" ;
//     has_values = false ;
//     return ;
//   }
//   if(nzones != 1) {
//     cerr << "can only read one zone cobalt grids!" ;
//     has_values = false ;
//     return ;
//   }
  
//   size_t nnodes, nodes_start ;
//   nnodes = npnts ;
//   nodes_start = 1 ;
//   pos.reserve(nnodes+ncells) ;

//   int faces_start, cells_start ;
//   faces_start = 1 ;
//   cells_start = 1 ;
  
//   for(size_t i=0;i<nnodes;++i) {
//     double x = 0, y = 0 ;
//     in >> x >> y ;
//     pos.push_back(positions(x,y)) ;
//   }

//   edge_list.reserve(nfaces) ;
//   vector<edges> bedges ;
//   vector<edges> bf2c ;
//   vector<edges> f2c ;
//   f2c.reserve(nfaces) ;
//   bf2c.reserve(nfaces) ;
//   bedges.reserve(nfaces) ;
//   for(int i=0;i<nfaces;++i) {
//     int fsz, n1,n2,c1,c2 ;
//     in >> fsz ;
//     if(fsz != 2) {
//       cerr << "two dimensional grids can only have faces with 2 nodes!"
//            << endl ;
//       return ;
//     }
//     in >> n1 >> n2 >> c1 >> c2 ;
//     if(c1<0)
//       std::swap(c1,c2) ;
//     c1 = c1-cells_start+nnodes ;
//     c2 = c2<0?c2:c2-cells_start+nnodes ;
//     n1 -= nodes_start ;
//     n2 -= nodes_start ;
//     if(c2 > 0) {
//       f2c.push_back(edges(c1,c2)) ;
//       edge_list.push_back(edges(n1,n2)) ;
//     } else {
//       bf2c.push_back(edges(c1,c2)) ;
//       bedges.push_back(edges(n1,n2)) ;
//     }
//   }
//   interior = edge_list.size() ;

//   if(!pos.empty()) {
//     minpos = pos[0] ;
//     maxpos = minpos ;
//   } else {
//     minpos = positions(0,0) ;
//     maxpos = positions(0,0) ;
//   }
  
//   for(size_t i=0;i<bedges.size();++i) {
//     maxpos.x = max(maxpos.x,pos[bedges[i].l].x) ;
//     maxpos.x = max(maxpos.x,pos[bedges[i].r].x) ;
//     minpos.x = min(minpos.x,pos[bedges[i].l].x) ;
//     minpos.x = min(minpos.x,pos[bedges[i].r].x) ;
//     maxpos.y = max(maxpos.y,pos[bedges[i].l].y) ;
//     maxpos.y = max(maxpos.y,pos[bedges[i].r].y) ;
//     minpos.y = min(minpos.y,pos[bedges[i].l].y) ;
//     minpos.y = min(minpos.y,pos[bedges[i].r].y) ; 

//     edge_list.push_back(bedges[i]) ;
//     f2c.push_back(bf2c[i]) ;
//   }
//   minview = minpos ;
//   maxview = maxpos ;
//   size = max(maxview.x - minview.x, maxview.y - minview.y);

//   val.reserve(nnodes+ncells) ;

//   string pixfilename ;
//   for(int i=filename.size()-1;i>=0;--i) {
//     if(filename[i] == '.') {
//       for(int j=0;j<i;++j)
//         pixfilename += filename[j] ;
//       i=0 ;
//     }
//   }
//   if(pixfilename == "") {
//     cerr << "unable to get pixfile name!" << endl ;
//     has_values = false ;
//     return ;
//   }
//   pixfilename += ".pix" ;
  
//   ifstream data(pixfilename.c_str(),ios::in) ;
//   if(data.fail()) {
//     cout << "unable to open pix file `" << pixfilename << "'" << endl ;
//     has_values = false;
//     return ;
//   }
//   while(data.peek() == '#') {
//     int ch ;
//     while((ch = data.get()) != EOF && ch != '\n')
//       /*NULL STATEMENT */ ;
//   }
//   int pix_nodes,v1,v2,v3,v4 ;
//   int nvec = 1 ;
//   data >> pix_nodes >> nvec >> v1 ;
//   data >> v2 >> v3 >> v4 ;

//   for(size_t i=0;i<nnodes;++i) {
//     if(data.fail() || data.eof()) {
//       has_values = false ;
//       cerr << "ran out of values" << endl ;
//       break ;
//     }
//     double vin ;
//     data >> vin ;
//     val.push_back(vin) ;
//   }

  
//   for(int i=0;i<ncells;++i) {
//     pos.push_back(positions(0,0)) ;
//     if(has_values)
//       val.push_back(double(0)) ;
//   }

//   vector<double> weight(pos.size()) ;
//   triangle_list.reserve(nfaces*2) ;
//   for(int i=0;i<nfaces;++i) {
//     int n1 = edge_list[i].l, n2 = edge_list[i].r ;
    
//     positions fc = (pos[n1]+pos[n2]) ;
//     double vc = val[n1]+val[n2] ;
//     if(f2c[i].l >= nnodes && f2c[i].l< pos.size()) {
//       pos[f2c[i].l] = pos[f2c[i].l] + fc ;
//       if(has_values)
//         val[f2c[i].l] += vc ;
//       weight[f2c[i].l] += 2.0 ;
//       triangle_list.push_back(triangles(n1,n2,f2c[i].l)) ;
//     }
//     if(f2c[i].r >= nnodes && f2c[i].r< pos.size()) {
//       pos[f2c[i].r] = pos[f2c[i].r] + fc ;
//       weight[f2c[i].r] += 2.0 ;
//       if(has_values)
//         val[f2c[i].r] += vc ;
//       triangle_list.push_back(triangles(n1,n2,f2c[i].r)) ;
//     }
//   }
//   for(size_t i=nnodes;i<pos.size();++i) {
//     double rweight = 1./weight[i] ;
//     pos[i] = pos[i]*rweight ;
//     val[i] *= rweight ;
//   }
// }

void VolGrid::set_value_range()
{
  vector<float> tmpVal;

  min_val = (val.size() == 0)?0:val[0] ;
  max_val = min_val ;
  for(size_t i=0;i<val.size();++i) {
    min_val = min(min_val,val[i]) ;
    max_val = max(max_val,val[i]) ;

    tmpVal.push_back(val[i]);
  }

  nth_element(tmpVal.begin(), tmpVal.begin()+tmpVal.size()/2, tmpVal.end());
  med_val = tmpVal[tmpVal.size()/2];
}

// void VolGrid::input(std::istream &in, bool read_values ) {
//   valid = true ;
//   int ichar = in.get() ;
//   if(in.eof()) {
//     cerr << "premature EOF while reading grid" << endl ;
//     valid = false ;
//   }
//   in.putback(ichar) ;
//   if(ichar == 'g') {
//     input_generalized(in,read_values) ;
//     set_value_range() ;
//   } else {
//     input_2dgv(in,read_values) ;
//     set_value_range() ;
//     optimize_edge_list() ;
//   }
//   if(has_values) {
//     generate_contour_curves(10);
//   }
// }

// void VolGrid::input(const char *file,bool read_values) {
//   filename = file ;
//   char *newfilename = tilde(file) ;
//   filename = newfilename ;
//   char *rp = rindex(newfilename,'.') ;
//   if(rp!=0 && !strcmp(rp,".cfg")) {
//     input_cobalt(string(filename),read_values) ;
//     set_value_range() ;
//     if(has_values) {
//       generate_contour_curves(10);
//     }
//   } else {
//     ifstream fio(newfilename,ios::in) ;
//     if(fio.fail()) {
//       cerr << "Open failed for " << newfilename << endl ;
//     } else {
//       free(newfilename) ;
//       input(fio,read_values) ;
//     }
//   }
// }

struct vertigo {
  positions p ;
  double v ;
  vertigo(positions pi, double vi) : p(pi),v(vi) {}
} ;

double sqr(double v) { return v*v ; }

void VolGrid::generate_contour_curves(int num_contours) {
  if (num_contours == 0)
    num_contours = 10 ;

  double range = (max_val-min_val)/ static_cast<double>(num_contours) ;
  double base = range==0?1:pow(10.0,double(floor(log10(range)))) ;
  contour_spacing = floor(range/base)*base;

  //  if(!has_values)
  //  return ;
  //  int num_contours = int(ceil((max_val-min_val)/contour_spacing)) ;
  int contour_base = int(ceil(min_val/contour_spacing)) ;
  if(num_contours <=0)
    num_contours = 1 ;
  num_contours++ ;

  if(num_contours > 500) {
    contour_spacing =(max_val-min_val)/500.0 ;
    num_contours = int(ceil((max_val-min_val)/contour_spacing)) ;
    contour_base = int(ceil(min_val/contour_spacing)) ;
    num_contours++ ;
  }
    
  contour_curves.clear() ;
  contour_values.clear() ;
  for(int i=0;i<num_contours;++i) {
    double currval = double(i+contour_base)*contour_spacing ;
    contour_values.push_back(contour_info(currval)) ;
  }
    
  double bignormal = sqrt(sqr(maxview.x-minview.x)+sqr(maxview.y-minview.y)) ;
  
  for(size_t i=0;i<triangle_list.size();++i) {
    const triangles &tp = triangle_list[i] ;
    vertigo v1(pos[tp.t1],val[tp.t1]),v2(pos[tp.t2],val[tp.t2]),
      v3(pos[tp.t3],val[tp.t3]) ;
    if(v1.v > v2.v)
      std::swap(v1,v2) ;
    if(v2.v > v3.v)
      std::swap(v2,v3) ;
    if(v1.v > v2.v)
      std::swap(v2,v1) ;
    if(v1.v == v3.v)
      continue ; // Skip degenerate case of constant value function
    int contour_index = int(ceil(v1.v/contour_spacing)) - contour_base ;
    double currval = double(contour_index+contour_base)*contour_spacing ;
    const double rdv31 = 1./(v3.v-v1.v+epsilon) ;
    const double rdv32 = 1./(v3.v-v2.v+epsilon) ;
    const double rdv12 = 1./(v1.v-v2.v+epsilon) ;
    
    while(currval <= v3.v) {
      positions p1 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;

      positions p2 ;
      if(currval > v2.v) 
        p2 = v2.p + (v3.p-v2.p)*(currval-v2.v)*rdv32 ;
      else 
        p2 = v2.p + (v1.p-v2.p)*(currval-v2.v)*rdv12 ;

      contour_curves.push_back(segments(p1,p2)) ;

      double gradient = (v3.v-v1.v)*bignormal/
        (contour_spacing*sqrt(sqr(v1.p.x-v3.p.x)+sqr(v1.p.y-v3.p.y))) ;
      if(contour_index < int(contour_values.size())) 
        if(gradient < contour_values[contour_index].mingrad) {
          contour_values[contour_index].mingrad = gradient ;
          contour_values[contour_index].mingradpos = 0.5*(p1+p2) ;
        }
        
      contour_index++ ;
      currval += contour_spacing ;
    }
  }
}

#include <Loci.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std ;
using namespace Loci ;
namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]) ;
  }
}

void kill_comment(istream &s) {
  while(s.peek() == ' ' || s.peek() == '\t' ||
        s.peek() == '\n' || s.peek() == '\r')
    s.get() ;
  
  while(s.peek() == '#') {
    char buf[1024] ;
    s.getline(buf,1023) ;
    while(s.peek() == ' ' || s.peek() == '\t' ||
          s.peek() == '\n' || s.peek() == '\r')
      s.get() ;
  }
}

struct bc_info {
  string name ;
  vector<pair<int,int> > edge_list ;
} ;

int main(int ac, char *av[]) {
  if(ac != 2) {
    cerr << "takes surface file as input" << endl ;
    return -1 ;
  }

  ifstream isurf(av[1],ios::in) ;

  if(isurf.fail()) {
    cerr << "unable to open file " << av[1] << endl ;
    return -1 ;
  }

  kill_comment(isurf) ;
  int npnts = 0 ;
  isurf >> npnts ;
  kill_comment(isurf) ;
  vector<double> x(npnts),y(npnts),z(npnts) ;
  for(int i=0;i<npnts;++i)
    isurf >> x[i] >> y[i] >> z[i] ;
  kill_comment(isurf) ;
  int nfaces = 0 ;
  isurf >> nfaces ;
  vector<int> f2size(nfaces) ;
  vector<int> f2node ;
  kill_comment(isurf) ;
  for(int i=0;i<nfaces;++i) {
    isurf >> f2size[i] ;
    for(int j=0;j<f2size[i];++j) {
      int nd = 0 ;
      isurf >> nd ;
      f2node.push_back(nd) ;
    }
  }

  kill_comment(isurf) ;
  int nbcs = 0;
  isurf >> nbcs ;
  vector<bc_info> bcs(nbcs) ;


  // Get grid edges
  vector<Array<int,3> > edge_map ;
  
  for(int i=0;i<nbcs;++i) {
    kill_comment(isurf) ;
    isurf >> bcs[i].name ;
    cout << "BC: " << bcs[i].name << endl ;
    kill_comment(isurf) ;
    int nedges =0 ;
    isurf >> nedges ;
    kill_comment(isurf) ;
    for(int j=0;j<nedges;++j) {
      int n1, n2 ;
      isurf >> n1 >> n2 ;
      Array<int,3> edge_info ;
      edge_info[0] = min(n1,n2) ;
      edge_info[1] = max(n1,n2) ;
      edge_info[2] = -(i+1) ;
      edge_map.push_back(edge_info) ;
      bcs[i].edge_list.push_back(pair<int,int>(n1,n2)) ;
    }
  }

  int off = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    Array<int,3> edge_info ;
    int n1 = f2node[off] ;
    int n2 = f2node[off+fsz-1] ;
    edge_info[0] = min(n1,n2) ;
    edge_info[1] = max(n1,n2) ;
    edge_info[2] = i+1 ; // Cell number
    edge_map.push_back(edge_info) ;
    for(int j=1;j<fsz;++j) {
      n1 = f2node[off+j-1] ;
      n2 = f2node[off+j] ;
      edge_info[0] = min(n1,n2) ;
      edge_info[1] = max(n1,n2) ;
      edge_info[2] = i+1 ;
      edge_map.push_back(edge_info) ;
    }
    off += fsz ;
  }

  if((edge_map.size() & 1) != 0) {
    cerr << "unmatched edges in surface?  Are all boundaries set?" << endl ;
    return -1 ;
  }
  sort(edge_map.begin(),edge_map.end()) ;

  int emsz = edge_map.size() ;
  for(int i=0;i<emsz;i+=2) {
    if(edge_map[i][0] != edge_map[i+1][0] ||
       edge_map[i][1] != edge_map[i+1][1]) {
      cerr << "unmatched edge in surface?  Are all boundaries set?" << endl ;
    }
  }

  ofstream ofile("grid.cog",ios::out) ;
  if(ofile.fail()) {
    cerr << "can't open grid.cog" << endl ;
    return  -1 ;
  }

  ofile.precision(16) ;
  int npatch = nbcs+2 ;
  ofile << "3 1 " << npatch << endl ;
  int vol_nfaces = nfaces*2 + emsz/2 ;
  int vol_ncells = nfaces ;
  int vol_npnts = npnts*2 ;
  int mxppf = 4 ;
  int mxfpc = 6 ;
  ofile << vol_npnts<< " " << vol_nfaces << " " << vol_ncells << " " << mxppf << " " << mxfpc << endl ;

  // output positions
  for(int i=0;i<npnts;++i) {
    double r = sqrt(y[i]*y[i]+z[i]*z[i]) ;
    y[i] = r ;
    z[i] = -0.05 ;
    ofile << x[i] << ' ' << y[i] << ' ' << z[i] << endl ;
  }
  for(int i=0;i<npnts;++i) {
    z[i] = 0.05 ;
    ofile << x[i] << ' ' << y[i] << ' ' << z[i] << endl ;
  }

  // now output mesh faces, starting with edge extrusion

  for(int i=0;i<emsz;i+=2) {
    // all extruded faces are quads
    ofile << 4 << ' ' ;
    ofile << edge_map[i][0] << ' ' << edge_map[i][1] << ' '
          << edge_map[i][1]+npnts << ' ' << edge_map[i][0] + npnts << ' '
          << edge_map[i][2] << ' ' << edge_map[i+1][2] << endl ;
  }
  // Now build symmetry planes
  int bc = -nbcs-1 ;
  off = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    ofile << fsz  ;
    for(int j=0;j<fsz;++j)
      ofile << ' ' << f2node[off+j] ;
    ofile << ' ' << i+1 << ' ' << bc << endl ;
    off += fsz ;
  }
  bc -= 1 ;
  off = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    ofile << fsz  ;
    for(int j=0;j<fsz;++j)
      ofile << ' ' << f2node[off+j]+npnts ;
    ofile << ' ' << i+1 << ' ' << bc << endl ;
    off += fsz ;
  }
  return 0;
}

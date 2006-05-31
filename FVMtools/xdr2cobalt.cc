#include <Loci.h>
//#include "read_grid.h"
//#include "sciTypes.h"
#include <iostream>
#include <fstream>

#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <string>

#include <rpc/rpc.h>
#include <rpc/xdr.h>

using std::string ;

//using namespace std ;
using std::ifstream ;
using std::cout ;
using std::cin ;
using std::endl ;
using std::cerr ;
using std::vector ;
using std::list ;
using std::ios ;

using std::ofstream ;

//using chem::vect3d ;

typedef vector3d<double> vect3d ;

void write_cb(fact_db &facts) {
  cout << "Writing out a cobalt file" << endl ; 
  Map cl,cr ;
  cl = facts.get_variable("cl") ;
  cr = facts.get_variable("cr") ;
  multiMap face2node ;
  face2node =facts.get_variable("face2node") ;
  constraint geom_cells ;
  geom_cells = facts.get_variable("geom_cells") ;

  
  param<int> maxppf,maxfpc ;
  maxppf = facts.get_variable("maxppf") ;
  maxfpc = facts.get_variable("maxfpc") ;

  int npatch = 0 ;
  entitySet faces = face2node.domain() ;
  for(entitySet::const_iterator ei=faces.begin();ei!=faces.end();++ei)
    npatch = min(cr[*ei],min(cl[*ei],npatch)) ;

  npatch = -npatch ;
  int mxppf = *maxppf;
  int mxfpc = *maxfpc;

  const_store<vect3d> pos ;
  pos = facts.get_variable("pos") ;
  entitySet nodes = pos.domain() ;

  int npnts = nodes.size() ;
  
  ofstream ofile("grid.cog",ios::out) ;
  if(ofile.fail()) {
    cerr << "can't open grid.cog" << endl ;
    return ;
  }

  ofile.precision(16) ;
  ofile << "3 1 " << npatch << endl ;
  int nfaces = faces.size() ;
  int ncells = (*geom_cells).size() ;
  ofile << npnts << " " << nfaces << " " << ncells << " " << mxppf << " " << mxfpc << endl ;

  for(entitySet::const_iterator ei=nodes.begin();ei!=nodes.end();++ei)
    ofile << pos[*ei] << endl ;
  
  for(entitySet::const_iterator i=faces.begin();i!=faces.end();++i) {
    ofile <<face2node.end(*i)-face2node.begin(*i) << " " ;
    for(const int* fp = face2node.begin(*i); fp != face2node.end(*i); ++fp)
      ofile << *fp << " " ;
    ofile << cl[*i] << " "  << cr[*i] <<endl ;
  }

  ofile.close();
}


bool readXDR(fact_db &facts,string filename) {
  entitySet nodes,faces,cells ;
    // First read in header information
    // Note:  Only processor 0 reads the file, it then sends the results
    // to other processors.
    int ndm ;
    int npnts, nfaces, ncells ;
    int dummy1,dummy2 ;
    FILE* FP = 0 ;
    XDR xdr_handle ;

    cerr << "opening file " << filename << endl ;
    FP = fopen(filename.c_str(), "r") ;
    if(FP == NULL) 
      return false ;

    cerr << "reading header info " << endl ;
    xdrstdio_create(&xdr_handle, FP, XDR_DECODE) ;
    if(!xdr_int(&xdr_handle, &ndm))
      return false ;
    if(!xdr_int(&xdr_handle, &dummy1))
      return false ;
    if(!xdr_int(&xdr_handle, &dummy2))
      return false ;
    if(!xdr_int(&xdr_handle, &npnts))
      return false ;
    if(!xdr_int(&xdr_handle, &nfaces))
      return false ;
    if(!xdr_int(&xdr_handle, &ncells))
      return false ;
    if(!xdr_int(&xdr_handle, &dummy1))
      return false ;
    if(!xdr_int(&xdr_handle, &dummy2))
      return false ;

    cerr << "reading positions" << endl ;
    // First read in node positions
    nodes = interval(1,npnts) ;
    store<vect3d> pos ;
    pos.allocate(nodes) ;
    for(entitySet::const_iterator ei=nodes.begin();ei!=nodes.end();++ei) {
      double tmp_pos[3] ;
      for(int i=0;i<3;++i)
        if(!xdr_double(&xdr_handle,&tmp_pos[i]))
          return false ;
      pos[*ei] = vect3d(tmp_pos[0],tmp_pos[1],tmp_pos[2]) ;
    }

    cerr << "reading cl,cr,count" << endl ;
    // Read in face left and right cells and count
    faces = interval(1,nfaces) ;
    Map cl,cr ;
    store<int> offset,count ;
    
    cl.allocate(faces) ;
    cr.allocate(faces) ;
    count.allocate(faces) ;
    Entity fp1 = faces.Max()+1 ;
    entitySet facesp1 = interval(1,fp1) ;
    offset.allocate(facesp1) ;
    for(entitySet::const_iterator ei=faces.begin();ei!=faces.end();++ei) {
      if(!xdr_int(&xdr_handle, &offset[*ei]))
        return false ;
      if(!xdr_int(&xdr_handle, &cl[*ei]))
        return false ;
      if(!xdr_int(&xdr_handle, &cr[*ei]))
        return false ;

      if(cl[*ei] < 0) 
        if(cr[*ei] < 0) {
          cerr << " boundary condition on both sides of a face?" << endl ;
          exit(1) ;
        } else {
          int tmp_swap = cr[*ei] ;
          cr[*ei] = cl[*ei] ;
          cl[*ei] = tmp_swap ;
        }
    }
    if(!xdr_int(&xdr_handle, &offset[fp1]))
      return false ;
    for(entitySet::const_iterator ei=faces.begin();ei!=faces.end();++ei) {
      count[*ei] = offset[*ei+1]-offset[*ei] ;
    }
    cout << "faces = " << faces << endl ;
    
    multiMap face2node ;
    face2node.allocate(count) ;
    Loci::MapRepP lp = Loci::MapRepP(cl.Rep()) ;
    Loci::MapRepP rp = Loci::MapRepP(cr.Rep()) ;
    cells = lp->image(faces) + rp->image(faces) ;
    store<int> ccount ;
    ccount.allocate(cells) ;
    for(entitySet::const_iterator ei=cells.begin();ei!=cells.end();++ei) {
      ccount[*ei] = 0 ;
    }
  
    int maxppf_val = 0 ;
    int maxfpc_val = 0 ;

    int face_stats[100] ;
    for(int i=0;i<100;++i)
      face_stats[i] = 0 ;
    
    
    for(entitySet::const_iterator ei=faces.begin();ei!=faces.end();++ei) {
      maxppf_val = max(maxppf_val,count[*ei]) ;
      face_stats[count[*ei]]++ ;
      ccount[cl[*ei]]++ ;
      ccount[cr[*ei]]++ ;
      for(int i=0;i<count[*ei];++i) {
        if(!xdr_int(&xdr_handle,&face2node[*ei][i]))
          return false ;
        face2node[*ei][i] += 1 ;
      }
    }

    for(int i=0;i<=maxppf_val;++i) {
      if(face_stats[i] != 0) {
        cout << "faces with " << i <<" edges = " << face_stats[i] << endl ;
      }
    }
    
    cells = cells & interval(0,Loci::UNIVERSE_MAX) ;
    for(entitySet::const_iterator ei=cells.begin();ei!=cells.end();++ei) {
      maxfpc_val = max(maxfpc_val,ccount[*ei]) ;
    }

    facts.create_fact("pos",pos) ;
    facts.create_fact("cl",cl) ;
    facts.create_fact("cr",cr) ;
    facts.create_fact("face2node",face2node) ;

    param<int> maxppf,maxfpc ;
    *maxppf = maxppf_val ;
    *maxfpc = maxfpc_val ;

    cout << "maxppf = " << maxppf_val << ", maxfpc = " << maxfpc_val << endl ;
    facts.create_fact("maxppf",maxppf) ;
    facts.create_fact("maxfpc",maxfpc) ;

    constraint geom_cells;
    *geom_cells = cells ;
    facts.create_fact("geom_cells",geom_cells) ;
    
    return true ;
}

int main(int ac, char *av[]) {
  if(ac != 2) {
    cerr << "provide grid file as argument." << endl ;
  }
  char *problem_name = av[1] ;
  ac-- ;
  av++ ;
  fact_db facts ;
  string gridfile = string(problem_name) + string(".xdr") ;
  if(!readXDR(facts,gridfile)) {
    cerr << "Unable to read file '" << gridfile << "'" << endl
         << "unable to continue." << endl ;
    return -1 ;
  }

  write_cb(facts) ;
  return 0 ;
}

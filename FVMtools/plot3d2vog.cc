//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci.h>
#include "vogtools.h"
//#include "Tools/xdr.h"
#include <Tools/stream.h>
//#include <Tools/ftrn_reader.h>
//#include <Tools/tools.h>
#include <algorithm>
using std::sort ;
#include <map>
using std::map ;
#include <vector>
using std::vector ;
using std::string ;
typedef double real ;
typedef Loci::vector3d<real> vect3d ;
#include <Tools/digraph.h>
using Loci::digraph ;
using Loci::component_sort ;

static char bcnames[6][4] = {"IJ1","IJN","IK1","IKN","JK1","JKN"} ;


class block_topo {
public:
  int ni,nj,nk ;
  int nodes_base ;
  int faces_base ;
  int ewfaces_base ;
  int nsfaces_base ;
  int tbfaces_base ;
  int cells_base ;
  block_topo() { ni=0,nj=0,nk=0,nodes_base=-1 ;}
  block_topo(int ni,int nj,int nk) ;
  ~block_topo() {};
  int num_nodes() { return ni*nj*nk ; }
  int num_faces() { return ni*(nj-1)*(nk-1)+(ni-1)*nj*(nk-1)+
		      (ni-1)*(nj-1)*(nk) ; }
  int num_cells() { return (ni-1)*(nj-1)*(nk-1) ; }
  
  int naddr(int i,int j,int k) { return i+j*ni+k*ni*nj + nodes_base; }
  int caddr(int i,int j,int k) 
  { return (i-1)+(j-1)*(ni-1)+(k-1)*(ni-1)*(nj-1)+cells_base ; }
  int ewaddr(int i,int j,int k) 
  { return i+(j-1)*(ni)+(k-1)*(ni)*(nj-1)+ewfaces_base ; }
  int nsaddr(int i,int j,int k) 
  { return (i-1)+j*(ni-1)+(k-1)*(ni-1)*(nj)+nsfaces_base ; }
  int tbaddr(int i,int j,int k) 
  { return (i-1)+(j-1)*(ni-1)+k*(ni-1)*(nj-1)+tbfaces_base; }
} ;

block_topo::block_topo(int NI, int NJ, int NK) {
  block_topo::ni = NI ;
  block_topo::nj = NJ ;
  block_topo::nk = NK ;
  
  cout << "block {" << ni << "," << nj << "," << nk << "}" << endl ;
  nodes_base = -1 ;
}

struct sort_order {
  const_store<vect3d> pos ;
  sort_order(const sort_order &ref) : pos(ref.pos.Rep()) {}
  explicit sort_order(const store<vect3d> &pin) :pos(pin.Rep()) { }
  inline bool operator()(int p1, int p2) {
    vect3d vp1 = pos[p1] ;
    vect3d vp2 = pos[p2] ;
    real norm1 = fabs(vp1.x)+fabs(vp1.y)+fabs(vp1.z) ;
    real norm2 = fabs(vp2.x)+fabs(vp2.y)+fabs(vp2.z) ;
    return norm1 < norm2 ;
  }
} ;

struct face_ident {
  int n[4] ;
  face_ident() {} ;
  face_ident(int n1,int n2, int n3, int n4) {
    n[0] = n1 ;
    n[1] = n2 ;
    n[2] = n3 ;
    n[3] = n4 ;
    sort(&n[0],&n[0]+4) ;
  }
  bool operator==(const face_ident &fi) const {
    return n[0] == fi.n[0] && n[1] == fi.n[1] &&
      fi.n[2] == fi.n[2] && n[3] == fi.n[3] ;
  }
} ;

struct boundary_segment_info {
  int block, face, Is,Ie,Js,Je ;
} ;

struct boundary_face_info {
  vector<boundary_segment_info> boundary_segments ;
  int tag ;
  entitySet boundary_faces ;
} ;

void usage() {
  cout << "Plot3d to VOG file converter usage" << endl
       << "This converter assumes that the grid file has a .grd postfix" << endl
       << " If only a single argument is given, then this will be " << endl
       << " grid file sans the .grd postfix, this will convert the" << endl
       << " grid to an unstructured VOG format grid file.  If the" << endl
       << " input file has more than one block, then the point matching" <<endl
       << " faces will be glued, and the remaining faces will be tagged" <<endl
       << " with a unique tag for each block and face. " << endl << endl ;

  cout << "If you want to control how tags are assigned to boundary" <<endl
       << " faces, then you can do so by providing a supplemental" << endl
       << " boundary specification file with the flag -bc" << endl << endl
       << "Example:" << endl
       << " plot3d2vog -bc grid.bc grid" << endl << endl ;

  cout << "The boundary condition specification file is contains the" <<endl
       << " following information:" << endl << endl
       << " <Number of boundary groups> " << endl
       << " For each group: " << endl
       << " <boundary tag> <number of segments> <boundary name> " << endl
       << " For each segment:" << endl
       << " <block number> <faceid> <index1 start> <index1 end> <index2 start> <index2 end>" << endl
       << "  Where:" << endl
       << " <faceid> is one of six strings: [IJ1,IJN,JK1,JKN,IK1,IKN]" << endl
       << " indices are given in the order indicated by the faceid string." 
       << endl <<endl;
  
  cout << "Use options -in, -ft, -cm, -m, or -Lref to set grid units."<<endl<<endl<<endl;
    
}

int main(int ac, char* av[]) {
  using namespace Loci ;
  using namespace VOG ;

  bool optimize = true ;
  
  Loci::Init(&ac,&av) ;

  if(ac == 1) {
    usage() ;
    exit(1) ;
  }
  constraint nodes ;
  store<vect3d> pos ;
  entitySet faces ;
  constraint quad_faces ;
  
  Map cl, cr ;
  Map pmap ;
  Map pci,pco ;
  MapVec<4> qnds ;
  Map face_flags ;
  
  constraint cells, geom_cells ;
  constraint geom_nodes ;
  Map north, south, east, west, top, bottom ;
  Map node_copy_input, node_copy_output ;
  
  char buf[512] ;
  entitySet degenerate,fo_reflecting ;
  entitySet faces_removed, faces_joined ;

  bool boundary_file = false ;
  string boundary_filename ;
  vector<boundary_face_info> boundaries_desc ;
  
  vector<string> combine_bc ;

  string Lref = "NOSCALE" ;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
        Loci::Finalize() ;
        exit(0) ;
      }
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-in")) {
      Lref = "1 inch" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-ft")) {
      Lref = "1 foot" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-cm")) {
      Lref = "1 centimeter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-m")) {
      Lref = "1 meter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mm")) {
      Lref = "1 millimeter" ;
      ac-- ;
      av++ ;
    }else if(ac >= 2 && !strcmp(av[1],"-c")) {
      combine_bc.push_back(string(av[2])) ;
      ac -= 2 ;
      av += 2 ;
    } else if (ac >= 2 && !strcmp(av[1],"-bc")) {
      boundary_file = true ;
      boundary_filename = string(av[2]) ;
      ac -= 2 ;
      av += 2 ;
    }else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }

  if(Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units." << endl ;
    exit(-1) ;
  }

if(Lref == "")
    Lref = "1 meter" ;
  
  if(!isdigit(Lref[0])) {
    Lref = string("1") + Lref ;
  }

   Loci::UNIT_type tp ;
  istringstream iss(Lref) ;
  iss >> tp ;
  double posScale = tp.get_value_in("meter") ;
  
  

  map<int, string> bcnamelist ;
  if(boundary_file) {
    ifstream bfile(boundary_filename.c_str(),ios::in) ;
    if(bfile.fail() || bfile.eof()) {
      cerr << "error opening boundary file " << boundary_filename << endl ;
      exit(-1) ;
    }
    int numbcs = 0;
    bfile >> numbcs ;
    Loci::parse::kill_white_space(bfile) ;
    
    for(int i=0;i<numbcs;++i) {
      if(bfile.fail() || bfile.eof()) {
        cerr << "error reading boundary file " << boundary_filename << endl ;
        exit(-1) ;
      }
      int boundary_flag = 0 ;
      bfile >> boundary_flag ;
      int boundary_segments = 0 ;
      bfile >> boundary_segments ;
      
      string boundary_name = Loci::parse::get_name(bfile) ;
      bcnamelist[boundary_flag] = boundary_name ;
      
      boundary_face_info fi ;
      fi.tag = boundary_flag ;
      for(int j=0;j<boundary_segments;++j) {
        int block, face, Is,Ie,Js,Je ;
        bfile >> block ;
        Loci::parse::kill_white_space(bfile) ;
        string fname = Loci::parse::get_name(bfile) ;
        face = -1 ;
        for(int f=0;f<6;++f)
          if(fname == bcnames[f])
            face = f ;
        if(face == -1) {
          cerr << "unknown face type " << fname
               << " in file " << boundary_filename
               << endl ;
          exit(-1) ;
        }
        Loci::parse::kill_white_space(bfile) ;
        bfile >> Is >> Ie >> Js >> Je ;
        boundary_segment_info si ;
        si.block = block ;
        si.face = face ;
        si.Is = Is ;
        si.Ie = Ie ;
        si.Js = Js ;
        si.Je = Je ;
        fi.boundary_segments.push_back(si) ;
      }
      boundaries_desc.push_back(fi) ;
    }
  }
  int base = 0 ;
  int num_blocks ;
  block_topo *blk ;
  char* filename = av[1] ;
  sprintf(buf,"%s.grd",filename) ;
  ifstream infile(buf,ios::in) ;
  if(infile.fail()) {
    cerr << "can't open '" << buf << "'" << endl ;
    exit(-1) ;
  }
  infile >> num_blocks ;
  if(num_blocks < 1 || num_blocks > 10000) {
    cerr << "ridiculous number of blocks, check file '" << buf <<"'" << endl ;
    exit(-1) ;
  }
  blk = new block_topo[num_blocks] ;
  for(int b=0;b<num_blocks;++b) {
    int ni,nj,nk ;
    infile >> ni >> nj >> nk ;
    blk[b] = block_topo(ni,nj,nk) ;
  }
  int sz = 0 ;
  int totsz = 0 ;
  for(int b=0;b<num_blocks;++b) {
    blk[b].nodes_base = totsz ;
    sz = max(sz,blk[b].num_nodes()) ;
    totsz += blk[b].num_nodes() ;
  }
  nodes = interval(base,base+totsz-1) ;
  base += totsz ;
  pos.allocate(*nodes) ;

  *geom_nodes = *nodes ;
  
  for(int b=0;b<num_blocks;++b) {
    sz = blk[b].num_nodes() ;
    int o = blk[b].nodes_base ;
    for(int i=0;i<sz;++i)
      infile >> pos[i+o].x ;
    for(int i=0;i<sz;++i)
      infile >> pos[i+o].y ;
    for(int i=0;i<sz;++i)
      infile >> pos[i+o].z ;
  }

  if(infile.fail()) {
    cerr << "had problems reading '" << buf << "', check file." << endl ;
    exit(-1) ;
  }
  for(int b=0;b<num_blocks;++b) {
    const int ni = blk[b].ni ;
    const int nj = blk[b].nj ;
    const int nk = blk[b].nk ;
    blk[b].faces_base = base ;
    blk[b].ewfaces_base = blk[b].faces_base ;
    blk[b].nsfaces_base = blk[b].faces_base + ni*(nj-1)*(nk-1) ;
    blk[b].tbfaces_base = blk[b].nsfaces_base + (ni-1)*(nj)*(nk-1) ;
    base += blk[b].num_faces() ;
  }
  
  quad_faces = interval(blk[0].faces_base,base-1) ;
  faces = *quad_faces ;
  for(int b=0;b<num_blocks;++b) {
    blk[b].cells_base = base ;
    base += blk[b].num_cells() ;
  }
  
  geom_cells = interval(blk[0].cells_base,base-1) ;
  
  qnds.allocate(faces) ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;
  for(entitySet::const_iterator fi=faces.begin();fi!=faces.end();++fi) {
    cl[*fi] = -1 ;
    cr[*fi] = -1 ;
    qnds[*fi][0] = -1 ;
    qnds[*fi][1] = -1 ;
    qnds[*fi][2] = -1 ;
    qnds[*fi][3] = -1 ;
  }
  
  north.allocate(*geom_cells) ;
  south.allocate(*geom_cells) ;
  east.allocate(*geom_cells) ;
  west.allocate(*geom_cells) ;
  top.allocate(*geom_cells) ;
  bottom.allocate(*geom_cells) ;
  
  entitySet fri,fli ;


  entitySet dontglue ;
  
  map<int,entitySet> boundary_tag_map ;
  if(boundary_file)
    for(size_t I=0;I<boundaries_desc.size();++I) {
      entitySet bcfcs ;
      int tag = boundaries_desc[I].tag ;
      for(size_t J=0;J<boundaries_desc[I].boundary_segments.size();++J) {
        const boundary_segment_info &bi =
          boundaries_desc[I].boundary_segments[J];
        int b = bi.block-1 ;
        int fc = bi.face ;
        const int ni = blk[b].ni ;
        const int nj = blk[b].nj ;
        const int nk = blk[b].nk ;
        int is,ie,js,je,i,j,k ;
        is = bi.Is ;
        ie = bi.Ie ;
        js = bi.Js ;
        je = bi.Je ;
        entitySet sfaces ;
        switch(fc) {
        case 0: // IJ1
          if(is < 1) {
            is = 1 ;
            ie = ni ;
            js = 1 ;
            je = nj ;
          }
          for(i=is;i<ie;++i)
            for(j=js;j<je;++j) {
              sfaces += blk[b].tbaddr(i,j,0) ;
            }

          
          break ;
        case 1: // IJN
          if(is < 1) {
            is = 1 ;
            ie = ni ;
            js = 1 ;
            je = nj ;
          }
            
          for(i=is;i<ie;++i)
            for(j=js;j<je;++j) {
              sfaces += blk[b].tbaddr(i,j,nk-1) ;
            }

          break ;
        case 2: // IK1
          if(is < 1) {
            is = 1 ;
            ie = ni ;
            js = 1 ;
            je = nk ;
          }
          for(i=is;i<ie;++i)
            for(k=js;k<je;++k) {
              sfaces += blk[b].nsaddr(i,0,k) ;
            }
          break ;
        case 3: // IKN
          if(is < 1) {
            is = 1 ;
            ie = ni ;
            js = 1 ;
            je = nk ;
          }
          for(i=is;i<ie;++i)
            for(k=js;k<je;++k) {
              sfaces += blk[b].nsaddr(i,nj-1,k) ;
            }
          break ;
        case 4: // JK1
          if(is < 1) {
            is = 1 ;
            ie = nj ;
            js = 1 ;
            je = nk ;
          }
          for(j=is;j<ie;++j)
            for(k=js;k<je;++k) {
              sfaces += blk[b].ewaddr(0,j,k) ;
            }
          break ;
        case 5: // JKN
          if(is < 1) {
            is = 1 ;
            ie = nj ;
            js = 1 ;
            je = nk ;
          }
          for(j=is;j<ie;++j)
            for(k=js;k<je;++k) {
              sfaces += blk[b].ewaddr(ni-1,j,k) ;
            }
          break ;
        }
        bcfcs += sfaces ;
        if(bi.Is != -2)
          dontglue += sfaces ;
      }
      boundary_tag_map[tag] = bcfcs ;
    }   

  
   
  typedef entitySet BCS[6] ;
  BCS *bcs ;
  bcs = new BCS[num_blocks] ;

  vector<entitySet> bset_faces(num_blocks) ;
  for(int b=0;b<num_blocks;++b) {
    const int ni = blk[b].ni ;
    const int nj = blk[b].nj ;
    const int nk = blk[b].nk ;
    int i,j,k ;
    /* compute topology */
    /* East-West Faces */
    for(i=0;i<ni;++i)
      for(j=1;j<nj;++j)
	for(k=1;k<nk;++k) {
	  int fc = blk[b].ewaddr(i,j,k) ;
	  qnds[fc][0] = blk[b].naddr(i,j-1,k-1) ;
	  qnds[fc][1] = blk[b].naddr(i,j,k-1) ;
	  qnds[fc][2] = blk[b].naddr(i,j,k ) ;
	  qnds[fc][3] = blk[b].naddr(i,j-1, k  ) ;
	}
    /* North-South Faces */
    for(i=1;i<ni;++i)
      for(j=0;j<nj;++j)
	for(k=1;k<nk;++k) {
	  int fc = blk[b].nsaddr(i,j,k) ;
	  qnds[fc][0] = blk[b].naddr(i-1,j,k-1) ;
	  qnds[fc][1] = blk[b].naddr(i-1,j,k  ) ;
	  qnds[fc][2] = blk[b].naddr(i  ,j,k) ;
	  qnds[fc][3] = blk[b].naddr(i,  j,k-1  ) ;
	}
    /* Top-Bottom Faces */
    for(i=1;i<ni;++i)
      for(j=1;j<nj;++j)
	for(k=0;k<nk;++k) {
	  int fc = blk[b].tbaddr(i,j,k) ;
	  qnds[fc][0] = blk[b].naddr(i-1,j-1,k) ;
	  qnds[fc][1] = blk[b].naddr(i,  j-1,k) ;
	  qnds[fc][2] = blk[b].naddr(i,j,  k) ;
	  qnds[fc][3] = blk[b].naddr(i-1,  j,  k) ;
	}
    
    /* First compute East-West Face cursors */
    for(i=1;i<ni;++i)
      for(j=1;j<nj;++j)
	for(k=1;k<nk;++k) {
	  int cc = blk[b].caddr(i,j,k) ;
	  int wside = blk[b].ewaddr(i-1,j,k) ;
	  int eside = blk[b].ewaddr(i,  j,k) ;
	  west[cc] = wside ;
	  east[cc] = eside ;
	  fli += eside ;
	  fri += wside ;
	  cl[eside] = cc ;
	  cr[wside] = cc ;
          bset_faces[b] += eside ;
          bset_faces[b] += wside ;
          
	}
    
    /* North-South Face Cursors */
    for(i=1;i<ni;++i)
      for(j=1;j<nj;++j)
	for(k=1;k<nk;++k) {
	  int cc = blk[b].caddr(i,j,k) ;
	  int nside = blk[b].nsaddr(i,j  ,k) ;
	  int sside = blk[b].nsaddr(i,j-1,k) ;
	  north[cc] = nside ;
	  south[cc] = sside ;
	  fli += nside ;
	  fri += sside ;
	  cl[nside] = cc ;
	  cr[sside] = cc ;
          bset_faces[b] += nside ;
          bset_faces[b] += sside ;
          
	}
    
    /* Top-Bottom Face Cursors */
    for(i=1;i<ni;++i)
      for(j=1;j<nj;++j)
	for(k=1;k<nk;++k) {
	  int cc = blk[b].caddr(i,j,k) ;
	  int tside = blk[b].tbaddr(i,j,k  ) ;
	  int bside = blk[b].tbaddr(i,j,k-1) ;
	  top[cc] = tside ;
	  bottom[cc] = bside ;
	  fli += tside ;
	  fri += bside ;
	  cl[tside] = cc ;
	  cr[bside] = cc ;
          bset_faces[b] += tside ;
          bset_faces[b] += bside ;
	}
    
    // Set up boundary faces information
    for(i=1;i<ni;++i)
      for(j=1;j<nj;++j) {
	bcs[b][0] += blk[b].tbaddr(i,j,0) ;
	bcs[b][1] += blk[b].tbaddr(i,j,nk-1) ;
      }
    for(i=1;i<ni;++i)
      for(k=1;k<nk;++k) {
	bcs[b][2] += blk[b].nsaddr(i,0,k) ;
	bcs[b][3] += blk[b].nsaddr(i,nj-1,k) ;
      }
    
    for(j=1;j<nj;++j)
      for(k=1;k<nk;++k) {
	bcs[b][4] += blk[b].ewaddr(0,j,k) ;
	bcs[b][5] += blk[b].ewaddr(ni-1,j,k) ;
      }
  }
  
  if(num_blocks>=1) {
    // If it is a multiblock grid, we need to search for
    // shared boundaries between blocks.
    entitySet bfaces ;
    for(int b=0;b<num_blocks;++b)
      for(int i=0;i<6;++i)
	bfaces += bcs[b][i] ;
    entitySet bnodes = Loci::MapRepP(qnds.Rep())->image(bfaces) ;
    entitySet allocnodes = bnodes ;
    if(dontglue != EMPTY) {
      entitySet ndg = Loci::MapRepP(qnds.Rep())->image(dontglue) ;
      entitySet glu = Loci::MapRepP(qnds.Rep())->image(bfaces-dontglue) ;
      bnodes -= (ndg-glu) ;
    }
    
    store<real> dr ;
    dr.allocate(allocnodes) ;
    vector<int> node_sort ;
    for(entitySet::const_iterator nd = bnodes.begin();nd!=bnodes.end();++nd) {
      dr[*nd] = 1e10 ;
      node_sort.push_back(*nd) ;
    }
    
    const real delta = 1e-2 ; // If a point is 100 times closer than the
    // smallest distance then assume that it is the same cell
    
    for(entitySet::const_iterator fc = bfaces.begin();fc!=bfaces.end();++fc) {
      const int cc = cl[*fc]<0?cr[*fc]:cl[*fc] ;
      
      int fcmap[3][2] ;
      fcmap[0][0] = north[cc] ;
      fcmap[0][1] = south[cc] ;
      fcmap[1][0] = east[cc] ;
      fcmap[1][1] = west[cc] ;
      fcmap[2][0] = top[cc] ;
      fcmap[2][1] = bottom[cc] ;
      
      real mindr = 1e50 ;
      
      for(int jj=0;jj<3;++jj) {
	real maxdr = 0 ;
	for(int kk=0;kk<2;++kk) {
	  const int fcc = fcmap[jj][kk] ;
	  const int qnds0 = qnds[fcc][0] ;
	  const int qnds1 = qnds[fcc][1] ;
	  const int qnds2 = qnds[fcc][3] ;
	  const int qnds3 = qnds[fcc][2] ;
	  
	  const vect3d d1 = (pos[qnds0]+pos[qnds2])-(pos[qnds1]+pos[qnds3]) ;
	  const real d1m = sqrt(0.25*dot(d1,d1)) ;
	  const vect3d d2 = (pos[qnds0]+pos[qnds1])-(pos[qnds2]+pos[qnds3]) ;
	  const real d2m = sqrt(0.25*dot(d2,d2)) ;
	  maxdr = max(maxdr,min(d1m,d2m)) ;
	}
	mindr = min(mindr,maxdr) ;
      }
      
      if(mindr < 1e-30) { 
	cerr << "mindr = " << mindr << endl ;
	cerr << "cc= " << cc << endl ;
	for(int jj=0;jj<3;++jj) {
	  for(int kk=0;kk<2;++kk) {
	    const int fcc = fcmap[jj][kk] ;
	    const int qnds0 = qnds[fcc][0] ;
	    const int qnds1 = qnds[fcc][1] ;
	    const int qnds2 = qnds[fcc][3] ;
	    const int qnds3 = qnds[fcc][2] ;
	    
	    const vect3d d1 = (pos[qnds0]+pos[qnds2])-
	      (pos[qnds1]+pos[qnds3]) ;
	    const real d1m = sqrt(0.25*dot(d1,d1)) ;
	    const vect3d d2 = (pos[qnds0]+pos[qnds1])-
	      (pos[qnds2]+pos[qnds3]) ;
	    const real d2m = sqrt(0.25*dot(d2,d2)) ;
	    cerr << "jj= " << jj << " kk= " << kk << endl ;
	    cerr.precision(16) ;
	    cerr << pos[qnds0] << endl ;
	    cerr << pos[qnds1] << endl ;
	    cerr << pos[qnds2] << endl ;
	    cerr << pos[qnds3] << endl ;
	    cerr << "d1m = " << d1m << " d2m= " << d2m << endl ;
	  }
	}
	exit(-1) ;
      }   
      mindr *= delta ;
      
      const int qnds0 = qnds[*fc][0] ;
      const int qnds1 = qnds[*fc][1] ;
      const int qnds2 = qnds[*fc][3] ;
      const int qnds3 = qnds[*fc][2] ;
      dr[qnds0] = min(dr[qnds0],mindr) ;
      dr[qnds1] = min(dr[qnds1],mindr) ;
      dr[qnds2] = min(dr[qnds2],mindr) ;
      dr[qnds3] = min(dr[qnds3],mindr) ;
      }
    sort(node_sort.begin(),node_sort.end(),sort_order(pos)) ;
    digraph equal_graph ;
    for(size_t i=0;i<node_sort.size();++i) {
      const int nsi = node_sort[i] ;
      const vect3d p1 = pos[nsi] ;
      const real dr1 = dr[nsi] ;
      for(size_t j=i+1;j<node_sort.size();++j) {
	const int nsj = node_sort[j] ;
	const vect3d p2 = pos[nsj] ;
	const vect3d dv = p1 - p2 ;
	const real dr2 = dr[nsj] ;
	const real dvni = abs(dv.x)+abs(dv.y)+abs(dv.z);
	if(abs((abs(p1.x)+abs(p1.y)+abs(p1.z))-
	       (abs(p2.x)+abs(p2.y)+abs(p2.z))) > dr1)
	  break ;
	if(dvni <= dr1 && dvni <= dr2) {
	  equal_graph.add_edge(nsi,nsj) ;
	  equal_graph.add_edge(nsj,nsi) ;
	}
      }
    }
    vector<digraph::nodeSet>
      components = component_sort(equal_graph).get_components() ;
    
    store<int> noderemap ;
    noderemap.allocate(*nodes) ;
    entitySet newnodes,interface_nodes ;
    for(entitySet::const_iterator ni = (*nodes).begin();ni!=(*nodes).end();++ni)
      noderemap[*ni] = *ni ;
    for(size_t i=0;i<components.size();++i) {
      entitySet eq = components[i] ;
      WARN(eq.size()<=1) ;
      int base_node = *(eq.begin()) ;
      newnodes += base_node ;
      interface_nodes += eq ;
      for(entitySet::const_iterator ni=eq.begin();ni!=eq.end();++ni) 
	noderemap[*ni] = base_node ;
    }
    entitySet deleted_nodes = interface_nodes - newnodes ;
    int num_deleted_nodes = deleted_nodes.size() ;
    
    entitySet copy_node_context ;
    if(num_deleted_nodes > 0)
      copy_node_context = interval(-num_deleted_nodes,-1) ;
    
    node_copy_input.allocate(copy_node_context) ;
    node_copy_output.allocate(copy_node_context) ;
    for(entitySet::const_iterator ei =deleted_nodes.begin();
	ei!=deleted_nodes.end();++ei) {
      node_copy_input[-num_deleted_nodes] = noderemap[*ei] ;
      node_copy_output[-num_deleted_nodes] = *ei ;
      num_deleted_nodes-- ;
    }
    
    *geom_nodes = *nodes - deleted_nodes ;
    for(entitySet::const_iterator fc=faces.begin();fc!=faces.end();++fc)
      for(int i=0;i<4;++i)
	qnds[*fc][i] = noderemap[qnds[*fc][i]] ;
    
    entitySet editfaces = Loci::MapRepP(qnds.Rep())->preimage(newnodes).first ;
    editfaces = editfaces & bfaces ;
    multiMap nd2face ;
    
    inverseMap(nd2face,qnds,newnodes,editfaces) ;
    
    FORALL(fli&editfaces,fc) {
      cr[fc] = -1 ;
    } ENDFORALL ;
    FORALL(fri&editfaces,fc) {
      cl[fc] = -1 ;
    } ENDFORALL ;
    
    
    for(entitySet::const_iterator
	  ni=newnodes.begin();ni!=newnodes.end();++ni) {
      for(const int *fi = nd2face.begin(*ni);fi!=nd2face.end(*ni);++fi) {
	if(cl[*fi] == -1 && cr[*fi] == -1)
	  continue ;
	
	if((qnds[*fi][0] == qnds[*fi][1] && qnds[*fi][2] == qnds[*fi][3]) ||
	   (qnds[*fi][0] == qnds[*fi][3] && qnds[*fi][1] == qnds[*fi][2])) {
	  cl[*fi] = -1 ; // Remove degenerate faces
	  cr[*fi] = -1 ;
	  faces_removed += *fi ;
	  continue ;
	}
	face_ident f1(qnds[*fi][0],qnds[*fi][1],qnds[*fi][2],qnds[*fi][3]) ;
	for(const int *fc = fi+1;fc!=nd2face.end(*ni);++fc) {
	  if((cl[*fc] == -1 && cr[*fc] == -1) || *fc == *fi)
	    continue ;
	  if(f1 == face_ident(qnds[*fc][0],qnds[*fc][1],
			      qnds[*fc][2],qnds[*fc][3])) {
	    const int cc = max(cl[*fc],cr[*fc]) ;
	    WARN(cl[*fi] != -1 && cr[*fi] != -1) ;
	    WARN(cl[*fc] != -1 && cr[*fc] != -1) ;
	    if(cl[*fi] == -1) 
	      cl[*fi] = cc ;
	    else
	      cr[*fi] = cc ;
	    
	    cl[*fc] = -1 ;
	    cr[*fc] = -1 ;
	    for(int i=0;i<4;++i)
	      qnds[*fc][i] = -1 ;
	    if(north[cc] == *fc)
	      north[cc] = *fi ;
	    else if(south[cc] == *fc)
	      south[cc] = *fi ;
	    else if(top[cc] == *fc)
	      top[cc] = *fi ;
	    else if(bottom[cc] == *fc)
	      bottom[cc] = *fi ;
	    else if(east[cc] == *fc)
	      east[cc] = *fi ;
	    else if(west[cc] == *fc)
	      west[cc] = *fi ;
	    else {
	      cerr << "something wrong connecting faces in multiblock" << endl ;
	      cerr << "cc = " << cc << endl ;
	      exit(-1) ;
	    }
	    faces_removed += *fc ;
	    faces_joined += *fi ;
	    break ;
	  }
	}
      }
    }
    entitySet notbc = faces_removed + faces_joined ;
    for(int b=0;b<num_blocks;++b) 
      for(int i=0;i<6;++i) 
	bcs[b][i] -= notbc ;
    fri -= faces_removed ;
    fli -= faces_removed ;
    fri += faces_joined ;
    fli += faces_joined ;

    map<int,entitySet>::iterator mi ;
    for(mi=boundary_tag_map.begin();mi!=boundary_tag_map.end();++mi) {
      if((mi->second & notbc) != EMPTY) {
        cerr << "Warning: faces glued for bc " << bcnamelist[mi->first] << endl ;
      }
      mi->second -= notbc ;
    }
  }
  
  entitySet right_boundary = fli&~fri ;
  entitySet left_boundary = ~fli & fri ;
  constraint boundary_faces ;
  boundary_faces = (right_boundary + left_boundary) ;
  
  // move all boundaries with left inside to right inside
  FORALL(left_boundary, fc) {
    cl[fc] = cr[fc] ;
  } ENDFORALL ;
  
  // now all boundaries have the left side inside and the right side outside
  
  
  entitySet boundary_idents = interval(-1,-(num_blocks*6)) ;
  
  entitySet boundaries ;
  if(boundary_file) {
    map<int,entitySet>::const_iterator mi ;
    for(mi=boundary_tag_map.begin();mi!=boundary_tag_map.end();++mi) {
      int tag = mi->first ;
      entitySet bcfaces = mi->second ;
      if((boundaries & bcfaces) != EMPTY) {
        cerr << "faces identified multiple times in BC specification!" << endl ;
        exit(-1) ;
      }
      boundaries |= bcfaces ;

      FORALL(bcfaces, bfc) {
        cr[bfc] = -tag ;
      } ENDFORALL ;
    }
    
  } else {
    
    for(int b=0;b<num_blocks;++b) {
      for(int i=0;i<6;++i) {
        int bid = -(1+b*6+i) ;
        int rbid = bid ;
        for(size_t j=0;j<combine_bc.size();++j) {
          if(bcnames[i] == combine_bc[j]) {
            rbid = -(1+i) ;
          }
        }
        FATAL((boundaries & bcs[b][i]) != EMPTY) ;
        boundaries |= bcs[b][i] ;
        FORALL(bcs[b][i], bfc) {
          cr[bfc] = rbid ;
        } ENDFORALL ;
      }
    }
  }

  if(boundaries != (left_boundary + right_boundary) ) {
    cerr << "not all boundaries were identified!" << endl ;

    entitySet tmpset = (left_boundary+right_boundary)-boundaries ;
    for(int b=0;b<num_blocks;++b) {
      entitySet overlap = bset_faces[b] & tmpset ;
      if(overlap.size() > 0) {
        cerr << "block number " << b+1 << " has " << overlap.size() << " missing boundary faces " << endl ;
      }
    }
    cerr << "missing boundaries = " << (left_boundary+right_boundary)-boundaries << endl ;
  }
  
  // Remove redundant faces from data-stucture
  entitySet real_faces = faces - faces_removed ;
  entitySet new_faces = interval(faces.Min(),faces.Min()+real_faces.size()-1) ;
  
  
  store<int> sizes ;
  sizes.allocate(new_faces) ;
  int rfc = new_faces.Min() ;
  FORALL(real_faces,fc) {
    int npf = 0 ;
    if(qnds[fc][0] != qnds[fc][1])
      npf++ ;
    if(qnds[fc][1] != qnds[fc][2])
      npf++ ;
    if(qnds[fc][2] != qnds[fc][3])
      npf++ ;
    if(qnds[fc][3] != qnds[fc][0])
      npf++ ;
    sizes[rfc++] = npf ;
  } ENDFORALL ;
  
  multiMap face2node ;
  Map ncl, ncr ;
  ncl.allocate(new_faces) ;
  ncr.allocate(new_faces) ;
  
  entitySet new_fobc ;
  face2node.allocate(sizes) ;
  rfc = new_faces.Min() ;
  
  Map color ;
  color.allocate(real_faces) ;
  FORALL(real_faces,fc) {
    int i=0 ;
    if(qnds[fc][0] != qnds[fc][1]) {
      face2node[rfc][i] = qnds[fc][0] ;
      i++ ;
    }
    if(qnds[fc][1] != qnds[fc][2]) {
      face2node[rfc][i] = qnds[fc][1] ;
      i++ ;
    }
    if(qnds[fc][2] != qnds[fc][3]) {
      face2node[rfc][i] = qnds[fc][2] ;
      i++ ;
    }
    if(qnds[fc][3] != qnds[fc][0]) {
      face2node[rfc][i] = qnds[fc][3] ;
      i++ ;
    }
    if(i!=sizes[rfc])
      cerr << "error, sizes didn't match" << endl ;
    
    ncr[rfc] = cr[fc] ;
    ncl[rfc] = cl[fc] ;
    color[fc] = rfc ;
    rfc++ ;
  } ENDFORALL ;
  
















  int npatch = 0 ;
  entitySet::const_iterator i ;
  for(i=new_faces.begin();i!=new_faces.end();++i) {
    npatch = min(min(npatch,ncl[*i]),ncr[*i]) ;
  }
  npatch = -npatch ;
  
  
  
  
  
  entitySet real_nodes ;
  for(i=new_faces.begin();i!=new_faces.end();++i) 
    for(const int* fp = face2node.begin(*i); fp != face2node.end(*i); ++fp)
      real_nodes += *fp ;
  store<int> nmap ;
  nmap.allocate(real_nodes) ;
  entitySet::const_iterator ei ;
  int npnts = 0 ;
  for(ei=real_nodes.begin();ei!=real_nodes.end();++ei) {
    npnts++ ;
    nmap[*ei] = npnts ;
  }

  
 
  //move the domain of data structures
  store<vector3d<double> > tpos;
  entitySet ndom = interval(0, npnts-1);
  tpos.allocate(ndom);
  entitySet::const_iterator ti = ndom.begin();

  
  for(ei=real_nodes.begin();ei!=real_nodes.end();++ei, ++ti) {
    tpos[*ti] = pos[*ei]*posScale;
   
  }

  //change the value of face2node
  entitySet::const_iterator fi ;
  for(fi=new_faces.begin();fi!=new_faces.end();++fi) {
    for(int k = 0; k < face2node[*fi].size(); k++){
      face2node[*fi][k] =  nmap[face2node[*fi][k]]-1 ;
    }
  }
  

  
 
  // establish face left-right orientation
  if(Loci::MPI_rank == 0) 
    cerr << "orienting faces" << endl ;
  VOG::orientFaces(tpos,ncl,ncr,face2node) ;
    
  if(Loci::MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(tpos,ncl,ncr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0) 
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(tpos,ncl,ncr,face2node) ;
  }
  
  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;

  
  //get boundary names
  vector<pair<int,string> > surf_ids ;
  for(int i=1;i<=npatch;++i){
    if(boundary_file) surf_ids.push_back(pair<int,string>(i, bcnamelist[i]));
    else{
      char buf[512] ;
      sprintf(buf,"BC_%d",i) ; 
      surf_ids.push_back(pair<int,string>(i, string(buf))) ;
    }
  }
  


 
  sprintf(buf,"%s.vog",filename) ;
  string outfile = string(buf);
  Loci::writeVOG(outfile, tpos, ncl, ncr, face2node,surf_ids) ;
 
  
  
  
  return 0 ;
}

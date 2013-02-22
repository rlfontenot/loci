#include <Loci>
#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
#include <string>
using std::string ;
#include <sstream>
using std::istringstream ;
#include <fstream>
using std::ifstream ;
using std::fstream ;
#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <utility>
#include <algorithm>
using std::max ;
using std::min ;
#include <math.h>
#include "dmatrix.h"
using std::pair ;

typedef  Loci::vector3d<double> vect3d;

int Usage() {
  cout << "incorrect command line arguments!" << endl
       << "Usage:" << endl
       << "computegeom <options> case" << endl
       << "  where case is the casename and the <options> may be " << endl
       << endl
       << "  -theta_r <angle> : ridge angle threshold" << endl
       << "  -theta_c <angle> : corner angle threshold" << endl
       << "  -geom_output     : output refined surface using computed geometry" << endl
       << "  -from_surf       : obtain input surface meshes from <case>.surf" 
       << endl
       << "                   : instead of <case>.vog" << endl 
       << endl ;
  exit(-1) ;
}

struct file_info {
  unsigned long numNodes ;
  unsigned long numFaces ;
  unsigned long numCells ;
  vector<short> cluster_sizes ;
} ;

unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern int fillClusterFaceSizes(unsigned char *cluster, int *sizes) ;
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Map &cl, Map &cr, int face_base) ;
  vector<unsigned char>
  encode_face_cluster(const multiMap &face2node,
                      const Map &cl, const Map &cr,
                      entitySet fcluster,
                      entitySet nodeSet,
                      entitySet cellSet) ;
}

struct surface_info {
  string name ;
  int id ;
  vector<Loci::Array<int,3> > trias ;
  vector<Loci::Array<int,4> > quads ;
  vector<vector<int> > gen_faces ;
} ;


void readSurfaces(string filename,
		  vector<surface_info> &surf_list,
		  vector<vect3d > &pos,
                  vector<int>& node_map) {

  surf_list.clear() ;
  pos.clear() ;

  map<int,int> surf_lookup ;

  // read in boundary names.
  vector<pair<int,string> > boundary_ids ;
  Loci::readBCfromVOG(filename,boundary_ids) ;
  map<int,string> surf_id ;
  for(size_t i=0;i<boundary_ids.size();++i)
    surf_id[boundary_ids[i].first] = boundary_ids[i].second ;

  hid_t input_fid ; 
  input_fid = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(input_fid <= 0) {
    cerr << "unable to open file '" << filename << "'"<< endl ;
    Usage() ;
  }
  
  hid_t face_g = H5Gopen(input_fid,"face_info") ;
  
  // Read cluster sizes
  hid_t dataset = H5Dopen(face_g,"cluster_sizes") ;
  hid_t dspace = H5Dget_space(dataset) ;
  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  vector<unsigned short> csizes(size) ;
  hsize_t dimension = size ;
  hsize_t stride = 1 ;
#ifdef H5_INTERFACE_1_6_4
  hsize_t start = 0 ;
#else
  hssize_t start = 0 ;
#endif
  hsize_t count = size ;
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
  int rank = 1 ;
  hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
  hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
  if(err < 0) {
    cerr << "unable to read cluster sizes from file '" <<
      filename << "'" << endl ;
    exit(-1) ;
  }
  H5Dclose(dataset) ;
  H5Sclose(memspace) ;

  // Read in clusters and transform
  dataset = H5Dopen(face_g,"cluster_info") ;
  dspace = H5Dget_space(dataset) ;
  start = 0 ;
  for(size_t c=0;c<size;++c) { // Loop over clusters
    count = csizes[c] ;
    vector<unsigned char> cluster(count) ;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
    dimension = count ;
    memspace = H5Screate_simple(rank,&dimension,NULL) ;
    err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster[0]) ;
    if(err < 0) {
      cerr << "unable to read cluster from file '" << filename << "'" << endl ;
    }
    start += count ;
    // Now scan cluster into local buffers

    int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
    Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
    Loci::store<int> fcnts ;
    fcnts.allocate(fclust) ;
    Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
    Loci::multiMap face2node ;
    Loci::Map cl,cr ;
    face2node.allocate(fcnts) ;
    cl.allocate(fclust) ;
    cr.allocate(fclust) ;
    Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0) ;
    // Now loop over faces in cluster to determine if any are boundary 
    // faces
    for(int f=0;f<nfaces;++f) {
      if(cr[f] < 0) { // boundary face 
	// Now check to see if we have encountered it before
	map<int,int>::const_iterator mi = surf_lookup.find(-cr[f]) ;
	int bc ;
	if(mi == surf_lookup.end() ) { // First time to see this boundary tag
	  // Get name of boundary
	  string bc_name ;
	  map<int,string>::const_iterator si = surf_id.find(-cr[f]) ;
	  if(si == surf_id.end()) {
	    char buf[128] ;
	    bzero(buf,128) ;
	    snprintf(buf,127,"BC_%d",-cr[f]) ;
	    bc_name = buf ;
	  } else {
	    bc_name = si->second ;
	  }
	  int bc_id = surf_list.size() ;
	  surf_list.push_back(surface_info()) ;
	  surf_list[bc_id].name = bc_name ;
          surf_list[bc_id].id = -cr[f] ;
	  surf_lookup[-cr[f]] = bc_id ;
	  bc = bc_id ;
	} else {
	  bc = mi->second ;
	}
	int fsz = face2node[f].size() ;
	if(fsz == 3) {
	  Loci::Array<int,3> tri ;
	  tri[0] = face2node[f][0] ;
	  tri[1] = face2node[f][1] ;
	  tri[2] = face2node[f][2] ;
	  surf_list[bc].trias.push_back(tri) ;
	} else if(fsz == 4) {
	  Loci::Array<int,4> qua ;
	  qua[0] = face2node[f][0] ;
	  qua[1] = face2node[f][1] ;
	  qua[2] = face2node[f][2] ;
	  qua[3] = face2node[f][3] ;
	  surf_list[bc].quads.push_back(qua) ;
	} else {
	  vector<int> tmp(fsz) ;
	  for(int i=0;i<fsz;++i)
	    tmp[i] = face2node[f][i] ;
	  surf_list[bc].gen_faces.push_back(tmp) ;
	}
      }
    }
  }

  // read in positions
  hid_t fi = H5Gopen(input_fid,"file_info") ;
  unsigned long numNodes = readAttributeLong(fi,"numNodes") ;
    
  H5Gclose(fi) ;

  count = numNodes ;

#ifdef H5_INTERFACE_1_6_4
    hsize_t lstart = 0 ;
#else
    hssize_t lstart = 0 ;
#endif
      
  // Read in pos data from file i
  vector<vect3d > pos_dat(numNodes) ;
  hid_t node_g = H5Gopen(input_fid,"node_info") ;
  dataset = H5Dopen(node_g,"positions") ;
  dspace = H5Dget_space(dataset) ;
      
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
  rank = 1 ;
  dimension = count ;
  memspace = H5Screate_simple(rank,&dimension,NULL) ;
  typedef Loci::data_schema_traits<vect3d > traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  hid_t datatype = dp->get_hdf5_type() ;
  err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
		      &pos_dat[0]) ;
  if(err < 0) {
    cerr << "unable to read positions from '" << filename << "'" << endl ;
    exit(-1) ;
  }
  H5Sclose(dspace) ;
  H5Dclose(dataset) ;
  H5Gclose(node_g) ;

  // Now mark all nodes that the faces access

  vector<int> used(numNodes) ;
  for(size_t i=0;i<numNodes;++i)
    used[i] = 0 ;
  
  int ssz = surf_list.size() ;

  
  for(int i=0;i<ssz;++i) {
    for(size_t j=0;j<surf_list[i].trias.size();++j) {
      used[surf_list[i].trias[j][0]] = 1 ;
      used[surf_list[i].trias[j][1]] = 1 ;
      used[surf_list[i].trias[j][2]] = 1 ;
    }
    for(size_t j=0;j<surf_list[i].quads.size();++j) {
      used[surf_list[i].quads[j][0]] = 1 ;
      used[surf_list[i].quads[j][1]] = 1 ;
      used[surf_list[i].quads[j][2]] = 1 ;
      used[surf_list[i].quads[j][3]] = 1 ;
    }
    for(size_t j=0;j<surf_list[i].gen_faces.size();++j) {
      for(size_t k=0;k<surf_list[i].gen_faces[j].size();++k)
	used[surf_list[i].gen_faces[j][k]] = 1 ;
    }
  }

  // Count nodes in the surface mesh
  int nnodes = 0 ;
  for(size_t i=0;i<numNodes;++i)
    if(used[i] != 0) {
      used[i] += nnodes++ ;
    }


  vector<vect3d > ptmp(nnodes) ;
   node_map.resize(nnodes);
  
  for(size_t i=0;i<numNodes;++i)
    {
      if(used[i] != 0){
      ptmp[used[i]-1] = pos_dat[i] ;
      node_map[used[i]-1] = i;
    }
    }
  pos.swap(ptmp) ;

  for(int i=0;i<ssz;++i) {
    for(size_t j=0;j<surf_list[i].trias.size();++j) {
      surf_list[i].trias[j][0] = used[surf_list[i].trias[j][0]] ;
      surf_list[i].trias[j][1] = used[surf_list[i].trias[j][1]] ;
      surf_list[i].trias[j][2] = used[surf_list[i].trias[j][2]] ;
    }
    for(size_t j=0;j<surf_list[i].quads.size();++j) {
      surf_list[i].quads[j][0] = used[surf_list[i].quads[j][0]] ;
      surf_list[i].quads[j][1] = used[surf_list[i].quads[j][1]] ;
      surf_list[i].quads[j][2] = used[surf_list[i].quads[j][2]] ;
      surf_list[i].quads[j][3] = used[surf_list[i].quads[j][3]] ;
    }
    for(size_t j=0;j<surf_list[i].gen_faces.size();++j) {
      for(size_t k=0;k<surf_list[i].gen_faces[j].size();++k)
	surf_list[i].gen_faces[j][k] = used[surf_list[i].gen_faces[j][k]] ;
    }
  }

#ifdef DEBUG
  cout << "nnodes = "<< nnodes << endl ;



  for(int i=0;i<ssz;++i) {
    cout << "surf = " << surf_list[i].name << endl ;
    cout << "ntrias=" << surf_list[i].trias.size() 
    << ",nquads=" << surf_list[i].quads.size()
    << ",ngenfs=" << surf_list[i].gen_faces.size()
         << endl;
  }
#endif

}



struct Tri {
  int t[3] ; // triangle nodes
  vect3d normal ;
  double angle[3] ;
  Tri(int p0, int p1, int p2) { t[0]=p0; t[1]=p1; t[2]= p2; } 
  Tri() {}
} ;

struct geomCoeff {
    vect3d b300;
    vect3d b030;
    vect3d b003;
    vect3d b210;
    vect3d b120;
    vect3d b021;
    vect3d b012;
    vect3d b102;
    vect3d b201;
    vect3d b111;
  vect3d loc(double u, double v) const {
    double w = 1.-u-v ;
    return (w*w*w)*b300 + (u*u*u)*b030 + (v*v*v)*b003+ 
      (3.0*w*w*u)*b210 + (3.0*w*u*u)*b120 + (3.0*w*w*v)*b201 +
      (3.0*u*u*v)*b021 + (3.0*w*v*v)*b102 + (3.0*u*v*v)*b012 +
      (6.0*w*u*v)*b111;
  }
} ;
//Overload ostream and istream (Input/Output) operators for struct geomCoeff
inline std::ostream & operator<<(std::ostream &s, const geomCoeff &g)
{
  s << g.b300 << ' ' << g.b030 << ' ' <<g.b003 << ' ' << g.b210 << ' '
    << g.b120 << ' ' << g.b021 << ' ' <<g.b012 << ' ' << g.b102 << ' '
    << g.b201 << ' ' << g.b111 <<endl;
  return s ;
}

inline std::istream &operator>>(std::istream &s, geomCoeff &g)
{
  s >> g.b300  >> g.b030  >>g.b003  >> g.b210 
    >> g.b120  >> g.b021  >>g.b012  >> g.b102 
    >> g.b201  >> g.b111 ;
  return s ;
}


struct Edge {
  int e[2] ;
  int f[2] ; 
} ;

struct nodeInfo {
  int primary_k ;
  vect3d primary_d ;
  vect3d e[3] ;
  double lambda[3] ;
} ;

struct edgeInfo {
  vect3d pmid[2] ;
} ;

bool operator<(const Edge &e1, const Edge &e2) {
  if(e1.e[0] < e2.e[0])
    return true ;
  if(e1.e[0] == e2.e[0] && e1.e[1] < e2.e[1])
    return true ;
  return false ;
}

void getEdges(vector<Tri> &trias, vector<Edge> &edges) {
  using std::swap ;
  int ntri = trias.size() ;
  vector<Edge> etmp(ntri*3) ;
  for(int i=0;i<ntri;++i) {
    Edge e ;
    for(int j=0;j<3;++j) {
      e.e[0] = trias[i].t[j] ;
      e.e[1] = trias[i].t[(j+1)%3] ;
      e.f[0] = i ;
      e.f[1] = -1 ;
      if(e.e[0] > e.e[1]) {
	swap(e.e[0],e.e[1]) ;
	swap(e.f[0],e.f[1]) ;
      }
      etmp[i*3+j] = e ;
    }
  }
  sort(etmp.begin(),etmp.end()) ;
  edges.clear() ;
  for(int i=0;i<ntri*3-1;++i) {
    if(etmp[i].e[0] == etmp[i+1].e[0] &&
       etmp[i].e[1] == etmp[i+1].e[1]) {
      etmp[i+1].f[0] = max(etmp[i].f[0],etmp[i+1].f[0]) ;
      etmp[i+1].f[1] = max(etmp[i].f[1],etmp[i+1].f[1]) ;
      if(etmp[i+1].f[0] < 0 || etmp[i+1].f[1] < 0) {
	cerr << "surface triangulation not topologically consistent" << endl ;
      }
      i++ ;
    }
    edges.push_back(etmp[i]) ;
  }
  if(etmp[ntri*3-2].e[0] != etmp[ntri*3-1].e[0] ||
     etmp[ntri*3-2].e[1] != etmp[ntri*3-1].e[1]) { 
    // last edge not paired, so insert it into edge list
    edges.push_back(etmp[ntri*3-1]) ;
  }
}

int getNode2Face(const vector<Tri> &trias,
		   const vector<vect3d> &pos,
		   vector<int> &node2face_offsets,
		   vector<pair<int,int> > &node2face) {
  vector<int> n2size(pos.size(),0) ;
  int tsz = trias.size() ;
  for(int i=0;i<tsz;++i)
    for(int j=0;j<3;++j)
      n2size[trias[i].t[j]] += 1 ;
  
  node2face_offsets = vector<int>(pos.size()+1,0) ;
  int nsz = pos.size() ;
  int max_f2n = 0 ;
  for(int i=1;i<nsz+1;++i) {
    max_f2n = max(max_f2n,n2size[i-1]) ;
    node2face_offsets[i] = node2face_offsets[i-1]+n2size[i-1] ;
  }
  int osz = node2face_offsets[nsz] ;
  node2face = vector<pair<int,int> >(osz) ;
  for(int i=0;i<tsz;++i)
    for(int j=0;j<3;++j) {
      int nd = trias[i].t[j] ;
      n2size[nd] -= 1 ;
      int loc = node2face_offsets[nd]+n2size[nd] ;
      node2face[loc].first = i ;
      node2face[loc].second = j ;
    }
  return max_f2n ;
}

//End copy from loci
//read in a .surf file 
void read_surf(string filename,
               vector<vect3d> &pos,
	       vector<Tri> &trias) {
  std::ifstream fin(filename.c_str(),std::ios::in) ;
  if(fin.fail()) {
    cerr << "unable to open input file '" << filename << "'" << endl ;
    exit(0);
  }
  
  int ntri =0, nquads=0, nnodes = 0;
  fin >> ntri >> nquads >> nnodes;
  int nface = ntri + nquads;
  if(nface==0 || nnodes ==0){
    cerr << " data file wrong -check it up" << endl;
    exit(0);
  }
  string s;
  getline(fin, s);
  double x, y, z, d1; 
  for(int i = 0; i < nnodes; i++){
    getline(fin, s);
    istringstream ins(s);
    ins >> x >> y >> z >> d1;   
    pos.push_back(vect3d(x, y, z));
  }

  int p1, p2, p3, p4,surf_id, rec_flag, bc_flag;
  for(int i = 0; i < ntri; i++){
    getline(fin, s);
    istringstream ins(s);
    ins >> p1 >> p2 >> p3 >> surf_id >> rec_flag >> bc_flag;
    trias.push_back(Tri(p1-1,p2-1,p3-1)) ;
  }
  for(int i = 0; i < nquads; i++){
    getline(fin, s);
    istringstream ins(s);
    ins >> p1 >> p2 >> p3 >>p4 >> surf_id >> rec_flag >> bc_flag;
    trias.push_back(Tri(p1-1,p2-1,p3-1)) ;
    trias.push_back(Tri(p1-1,p3-1,p4-1)) ;
  }
  fin.close();
}


void computeTriangleProperties(vector<Tri> &trias, const vector<vect3d> &pos) {
  int ntrias = trias.size() ;
  for(int i=0;i<ntrias;++i) {
    vect3d p0 = pos[trias[i].t[0]] ;
    vect3d p1 = pos[trias[i].t[1]] ;
    vect3d p2 = pos[trias[i].t[2]] ;
    vect3d n = cross(p1-p0,p2-p0) ;
    n *= 1./max(norm(n),1e-30) ;
    trias[i].normal = n ;
    vect3d v01 = p1-p0 ;
    vect3d v02 = p2-p0 ;
    vect3d v11 = p0-p1 ;
    vect3d v12 = p2-p1 ;
    vect3d v21 = p0-p2 ;
    vect3d v22 = p1-p2 ;
    trias[i].angle[0] = acos(max(-1.0,min(1.0,dot(v01,v02)/(norm(v01)*norm(v02))))) ;
    trias[i].angle[1] = acos(max(-1.0,min(1.0,dot(v11,v12)/(norm(v11)*norm(v12))))) ;
    trias[i].angle[2] = acos(max(-1.0,min(1.0,dot(v21,v22)/(norm(v21)*norm(v22))))) ;
  }
}

// Do Eigensystem analysis of face normals surrounding node
// Input:
// N[nn] array of normals
// W[nn] array of weights
// nn number of normal inputs
// e[3] output of 3 eigenvectors
// lambda[3] output of 3 eigenvectors
void eigenAnalysis(const vect3d N[],const double W[], int nn,  
		   vect3d e[3], double lambda[3]) 
{ 

  dmatrix sqrtW_N = dmatrix(1,nn,1,3); //result of sqrtW*N

  for(int i = 0; i < nn; i++){
    double sqrtW = sqrt(W[i]) ;
    sqrtW_N[i+1][1] = sqrtW*N[i].x ;
    sqrtW_N[i+1][2] = sqrtW*N[i].y ;
    sqrtW_N[i+1][3] = sqrtW*N[i].z ;
  }
  //single value decompose SQRT(W)*N = U SQRT(lambda) V^T, 
  // where A = N^T W N = V lambda V^T
  dvector w = dvector(1,3); //eigenvalues from dcv of A
  dmatrix v = dmatrix(1, 3, 1, 3); //eigenvectors from dcv of A 
  svdcmp(sqrtW_N, nn, 3, w, v);

  // svd found sqrt of lambda, so square
  for(int i=0;i<3;++i)
    lambda[i] = w[i+1]*w[i+1] ;
  
  // convert eigenvectors to vect3d
  for(int i=0;i<3;++i)
    e[i] = vect3d(v[1][i+1],v[2][i+1],v[3][i+1]) ;
}

// Classify node:
// 1 = smooth
// 2 = ridge
// 3 = cusp (corner)
// Input
// vect3d N[] : normal array
// double A[] : angle array
// int m : size of N and A arrays
// vect3d e[3] : eigenvectors
// double lambda[3] : eigenvalues
// returns classification
static double Xc = 0.2 ;
static double Xr = 0.03 ;
int classifyNode(vect3d N[], double A[], int m, 
		 vect3d e[3] , double lambda[3]) {
  double sumA = 0;
  for(int i=0;i<m;++i)
    sumA += A[i] ;
  double theta_a = fabs(sumA-2.0*M_PI) ;
  //  double Xc = 0.2 ; // Suggested by Jiao's paper
  //  double Xr = 0.03 ;
  if(theta_a > .45*M_PI || lambda[2] >= Xc*lambda[0])
    return 3 ;
  for(int i=0;i<m;++i)
    if(dot(e[0],N[i]) <= 0)
      return 2 ;
  if(lambda[1] >= Xr*lambda[0])
    return 2 ;
  return 1 ;
}

void getNodeInfo(const vector<Tri> &trias,
		 const vector<int> &node2face_offsets,
		 const vector<pair<int,int> > &node2face,
		 int max_f2n,
		 vector<nodeInfo> &ninfo) {
  vector<vect3d> N(max_f2n) ;
  vector<double> A(max_f2n) ;
  int nn = ninfo.size() ;
  for(int i=0;i<nn;++i) {
    int m = node2face_offsets[i+1]-node2face_offsets[i] ;
    int o = node2face_offsets[i] ;
    for(int j=0;j<m;++j) {
      int t = node2face[o+j].first ;
      int nd = node2face[o+j].second ;
      N[j] = trias[t].normal ;
      A[j] = trias[t].angle[nd] ;
    }
    eigenAnalysis(&N[0],&A[0],m,ninfo[i].e,ninfo[i].lambda) ;
    int primary_k = classifyNode(&N[0],&A[0],m,
				 ninfo[i].e,ninfo[i].lambda) ;
    ninfo[i].primary_k = primary_k ;

    // compute b
    vect3d b = vect3d(0,0,0);
    for(int j=0;j<m;++j)
      b +=  N[j]*A[j] ;
    // Compute primary direction
    vect3d d(0,0,0) ;
    for(int j=0;j<primary_k;++j)
      if(ninfo[i].lambda[j] > 0.003*ninfo[i].lambda[0])
	d += dot(ninfo[i].e[j],b)*ninfo[i].e[j]/ninfo[i].lambda[j] ;
    // Normalize d
    d *= 1./max(norm(d),1e-30) ;
    ninfo[i].primary_d = d ;
  }
}
void edgeReconstruct(const vector<Edge> &edges,
		     const vector<vect3d> &pos,
		     const vector<nodeInfo> &ninfo,
		     const vector<Tri> &trias,
		     vector<edgeInfo> &edge_points) {
  int esz = edges.size() ;
  int nsz = pos.size() ;
  vector<bool> ridge(esz,false) ;
  vector<int> pk(nsz) ;
  for(int i=0;i<nsz;++i)
    pk[i] = ninfo[i].primary_k ;

  for(int e=0;e<esz;++e) {
    int n0 = edges[e].e[0] ;
    int n1 = edges[e].e[1] ;
    if((ninfo[n0].primary_k==2 && ninfo[n1].primary_k > 1) ||
       (ninfo[n1].primary_k==2 && ninfo[n0].primary_k > 1)) {
      double ndotn = 1.0 ;
      if(edges[e].f[0] >= 0 && edges[e].f[1] >=0) 
	ndotn = dot(trias[edges[e].f[0]].normal,
		    trias[edges[e].f[1]].normal) ;
      if(ndotn < 0.95) { // could be a ridge, check tangency
	ridge[e] = true ;
	vect3d dp = pos[n1]-pos[n0] ;
	if(ninfo[n0].primary_k==2) {
	  double angle = fabs(dot(dp,ninfo[n0].e[2])/
			      max(norm(dp)*norm(ninfo[n0].e[2]),1e-30)) ;
	  if(angle < .7)
	    ridge[e] = false ;
	}
	if(ninfo[n1].primary_k==2) {
	  double angle = fabs(dot(dp,ninfo[n1].e[2])/
			      max(norm(dp)*norm(ninfo[n1].e[2]),1e-30)) ;
	  if(angle < .7)
	    ridge[e] = false ;
	}
      }
    }
  }
  vector<int> ncnt(nsz,0) ;
  for(int e=0;e<esz;++e) {
    int n0 = edges[e].e[0] ;
    int n1 = edges[e].e[1] ;
    if(ridge[e]) {
      ncnt[n0] += 1 ;
      ncnt[n1] += 1 ;
    }
  }
  // Find spurious corners
  for(int i=0;i<nsz;++i)
    if(ncnt[i] > 2)
      pk[i] = 3 ;
  
  
    

                        
  for(int e=0;e<esz;++e) {
    int n0 = edges[e].e[0] ;
    int n1 = edges[e].e[1] ;
    vect3d p0 = pos[n0] ;
    vect3d p1 = pos[n1] ;
    vect3d norm0 = ninfo[n0].primary_d ;
    vect3d norm1 = ninfo[n1].primary_d ;
    vect3d dp = p1-p0 ;
    // Now check for corner points
    if(pk[n0]==3 && pk[n1]==3) {
      norm0 = vect3d(0.,0.,0.) ;
      norm1 = vect3d(0.,0.,0.) ;
    }



    // If not a ridge and both sides of edge do not provide good
    // normal info, search triangles for good info
    if(!ridge[e] && pk[n0]!=1 && pk[n1]!=1) {
      norm0 = vect3d(0.,0.,0.) ;
      norm1 = vect3d(0.,0.,0.) ;
      
      int t1 = edges[e].f[0] ;
      int t2 = edges[e].f[1] ;
      vect3d nedge = vect3d(0.,0.,0.) ;
      if(t1 >=0) {
	nedge += trias[t1].normal ;
      }
      if(t2 >=0) {
	nedge += trias[t2].normal ;
      }
      nedge *= 1./norm(nedge) ;

      if(pk[n0] == 2) { 
	// Adjust normal to be orthoganal to edge tangent
	vect3d etan = ninfo[n0].e[2] ;
	etan *= 1./max(norm(etan),1e-30) ;
	norm0 = nedge - dot(nedge,etan)*etan ;
      }
      if(pk[n1] == 2) {
	// Adjust normal to be orthoganal to edge tangent
	vect3d etan = ninfo[n1].e[2] ;
	etan *= 1./max(norm(etan),1e-30) ;
	norm1 = nedge - dot(nedge,etan)*etan ;
      }

#ifdef OLDWAY
      if(t1 >0 && t2 > 0) {

        int nt1 = trias[t1].t[0] ;
        if(trias[t1].t[1] != n0 && trias[t1].t[1] != n1)
          nt1 = trias[t1].t[1] ;
        if(trias[t1].t[2] != n0 && trias[t1].t[2] != n1)
          nt1 = trias[t1].t[2] ;
        int nt2 = trias[t2].t[0] ;
        if(trias[t2].t[1] != n0 && trias[t2].t[1] != n1)
          nt2 = trias[t2].t[1] ;
        if(trias[t2].t[2] != n0 && trias[t2].t[2] != n1)
          nt2 = trias[t2].t[2] ;
        // Now average projected normals from good edges
        double w1 = 0 ;
        double w2 = 0 ;
        if(pk[nt1] == 1)
          w1 = 1.0 ;
        if(pk[nt2] == 1)
          w2 = 1.0 ;
        // note, if search was not able to find a good candidate, then
        // we just revert back to zero normals which gives a linear curve
        // for the edge.  In this case we would really like to go to a more
        // advanced algorithm.  For surfaces meshed for CFD applications, the
        // present algorithm is usually more than sufficient
        w1 *= 1./max(w1+w2,1e-30) ;
        w2 *= 1./max(w1+w2,1e-30) ;
        vect3d dp0 = pos[nt1]-p0 ;
        
        vect3d n1 = ninfo[nt1].primary_d ;
        vect3d n2 = ninfo[nt2].primary_d ;
        n1 *= 1./max(norm(n1),1e-30) ;
        n2 *= 1./max(norm(n2),1e-30) ;
        
        norm0 += w1*(n1 - (2.*dot(dp0,n1)/dot(dp0,dp0))*dp0) ;
        dp0 = pos[nt2]-p0 ;
        
        norm0 += w2*(n2 - (2.*dot(dp0,n2)/dot(dp0,dp0))*dp0) ;
        norm0 *= 1./max(norm(norm0),1e-30) ;
        vect3d dp1 = pos[nt1]-p1 ;
        
        norm1 += w1*(n1 - (2.*dot(dp1,n1)/dot(dp1,dp1))*dp1) ;
        dp1 = pos[nt2]-p1 ;
        
        norm1 += w2*(n2 - (2.*dot(dp1,n2)/dot(dp1,dp1))*dp1) ;
        norm1 *= 1./max(norm(norm1),1e-30) ;
      }
#endif
    }
    // If node is connected connected to ridge or corner then extrapolate
    // normal to ridge/corner
    if(pk[n0]!=1 && pk[n1]==1) {
      norm0 = norm1 - (2.*dot(dp,norm1)/dot(dp,dp))*dp ;
    }
    if(pk[n0]==1 && pk[n1]!=1) {
      norm1 = norm0 - (2.*dot(dp,norm0)/dot(dp,dp))*dp ;
    }
    
    edge_points[e].pmid[0] = (2.*p0 + p1 - dot(dp,norm0)*norm0)/3.0 ;
    edge_points[e].pmid[1] = (2.*p1 + p0 + dot(dp,norm1)*norm1)/3.0 ;
    
    if(ridge[e] && ! (pk[n0] == 3  && pk[n1]==3 ) ) {
      // Compute based on tangent lines
      vect3d tan0 = ninfo[n0].e[2] ;
      tan0 *= 1./max(norm(tan0),1e-30) ;
      if(dot(tan0,dp) < 0.0)
	tan0 *= -1.0 ;
      vect3d tan1 = ninfo[n1].e[2] ;
      tan1 *= 1./max(norm(tan1),1e-30) ;
      if(dot(tan1,dp) < 0.0)
	tan1 *= -1.0 ;
      if(ninfo[n0].primary_k == 3) {
	tan0 = tan1 - 2.*(tan1 - (dot(tan1,dp)/dot(dp,dp))*dp) ;
      }
      if(ninfo[n1].primary_k == 3) {
	tan1 = tan0 - 2.*(tan0 - (dot(tan0,dp)/dot(dp,dp))*dp) ;
      }
      edge_points[e].pmid[0] = p0+tan0*norm(dp)/3. ;
      edge_points[e].pmid[1] = p1-tan1*norm(dp)/3. ;
    }
  }
}
void  geoReconstruct(const vector<Tri> &trias, const vector<Edge> &edges, 
		     const vector<edgeInfo> &edge_points, 
		     const vector<vect3d> &pos,
		     vector<geomCoeff> &trigeo) {
  using std::swap ;
  int fsz = trias.size() ;
  // collect corner points
  for(int f=0;f<fsz;++f) {
    trigeo[f].b300 = pos[trias[f].t[0]] ;
    trigeo[f].b030 = pos[trias[f].t[1]] ;
    trigeo[f].b003 = pos[trias[f].t[2]] ;
  }
  // collect edge curve information
  int esz = edges.size() ;
  for(int e=0;e<esz;++e) {
    const int e0 = edges[e].e[0] ;
    const int e1 = edges[e].e[1] ;
    for(int i=0;i<2;++i) {
      int f = edges[e].f[i] ;
      if(f >= 0) {
	int le0 = 0 ;
	int le1 = 0 ;
	for(int lf=0;lf<3;++lf) {
	  if(e0 == trias[f].t[lf])
	    le0 = lf ;
	  if(e1 == trias[f].t[lf])
	    le1 = lf ;
	}
	vect3d pp0 = edge_points[e].pmid[0] ;
	vect3d pp1 = edge_points[e].pmid[1] ;
	if(le0 > le1) {
	  swap(le0,le1) ;
	  swap(pp0,pp1) ;
	}
	if(le0 == 0) {
	  if(le1 == 1) {
	    // Edge between P1 and P2
	    trigeo[f].b210 = pp0 ;
	    trigeo[f].b120 = pp1 ;
	  } else {
	    // Edge between P1 and P3
	    trigeo[f].b201 = pp0 ;
	    trigeo[f].b102 = pp1 ;
	  }
	} else {
	  // Edge between P2 and P3
	  trigeo[f].b021 = pp0 ;
	  trigeo[f].b012 = pp1 ;
	}
      }
    }
  }
  // compute center point
  for(int f=0;f<fsz;++f) {
    vect3d E = (trigeo[f].b210 + trigeo[f].b120 + 
		trigeo[f].b021 + trigeo[f].b012 + 
		trigeo[f].b201 + trigeo[f].b102)/6. ;
    vect3d V = (trigeo[f].b300 + trigeo[f].b030+trigeo[f].b003)/3. ;
    trigeo[f].b111 = E+0.5*(E-V) ;
  }
}
void outputGeom(string geo_file, 
		const vector<vect3d> &pos,
		const vector<Tri> &trias,
		const vector<geomCoeff> &trigeo) {
  //output the geometry into .coeff file 
  std::ofstream ofile(geo_file.c_str()); 

  int num_nodes = pos.size() ;
  ofile.precision(14) ;
  //write out positions  of nodes
  ofile << num_nodes << endl;
  for(int i = 0; i < num_nodes; i++){
    ofile << pos[i].x<< " "  << pos[i].y<< " " <<pos[i].z<< endl;
  }
  
  int num_geom = trias.size() ;
  //write out face2node
  ofile << num_geom << endl;
  for(int i = 0; i < num_geom; i++){
    ofile << trias[i].t[0] << ' '
          << trias[i].t[1] << ' '
          << trias[i].t[2] << endl;
  }

  //write out GeomCoeff
  ofile.precision(10);
  for(int i = 0; i < num_geom; i++){
    geomCoeff g = trigeo[i];
      ofile << g.b300<<' '
            << g.b030<<' '
            << g.b003<<' '
            << g.b210<<' '
            << g.b120<<' '
            << g.b021<<' '
            << g.b012<<' '
            << g.b102<<' '
            << g.b201<<' '
            << g.b111<<endl;
  }
  ofile.close();
}


void outputGeomSurf(string geo_file, 
		const vector<geomCoeff> &trigeo) {
  //output the geometry into .coeff file 
  std::ofstream ofile(geo_file.c_str()); 
  int ngt = trigeo.size() ;
  ofile << ngt*16 << ' ' << 0 << ' ' << ngt*15 << endl ;
  ofile.precision(14) ;
  // Write out points
  static const double ulist[15] = {0.,0.,0.,0.,0.,
				   .25,.25,.25,.25,
				   .5,.5,.5,
				   .75,.75,
				   1.} ;
  static const double vlist[15] = {0.,.25,.5,.75,1.,
				   0.,.25,.5,.75,
				   0.,.25,.5,
				   0.,.25,
				   0.} ;

  for(int i=0;i<ngt;++i) {
    for(int j=0;j<15;++j) {
      vect3d p = trigeo[i].loc(ulist[j],vlist[j]) ;
      ofile << p.x << ' ' << p.y << ' ' << p.z << ' ' << 0 << endl ;
    }
  }
  static const int t1[16] = { 0, 1, 1, 2, 2, 3, 3, 5, 6, 6, 7, 7, 9,10,10,12} ;
  static const int t2[16] = { 1, 6, 2, 7, 3, 8, 4, 6,10, 7,11, 8,10,13,11,13} ;
  static const int t3[16] = { 5, 5, 6, 6, 7, 7, 8, 9, 9,10,10,11,12,12,13,14} ;
  
  int ploc = 1 ;
  for(int i=0;i<ngt;++i) {
    for(int j=0;j<16;++j) {
      ofile << ploc+t1[j] << ' ' << ploc+t2[j] << ' ' << ploc+t3[j] << ' '
	    << 1 << ' ' << 0 << ' ' << 1 << endl ;
    }
    ploc += 15 ;
  }
  
}

int main(int ac, char *av[]) {
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << av[0] << " is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  double theta_r = 15 ; // ridge angle
  double theta_c = 42 ; // corner angle
  bool from_surf = false ;
  bool output_geom = false ;
  while(ac > 2) {
    if(ac > 2 && !strcmp(av[1],"-geom_output")) {
      output_geom = true ;
      ac-- ;
      av++ ;
    } else if(ac > 2 && !strcmp(av[1],"-from_surf")) {
      from_surf = true ;
      ac-- ;
      av++ ;
    } else if(ac > 3 && !strcmp(av[1],"-theta_r")) {
      theta_r = atof(av[2]) ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac > 3 && !strcmp(av[1],"-theta_c")) {
      theta_c = atof(av[2]) ;
      ac -= 2 ;
      av += 2 ;
    } else {
      cerr << "unrecognized argument " << av[1] << endl ;
      Usage() ;
    }
  }

  if(ac != 2) return Usage();
      
  string surf_file = av[1];
  string geo_file = av[1] ;
  string vog_file = av[1] ;
  vog_file += ".vog" ;
  surf_file += ".surf" ;
  geo_file += ".geom" ;
  Xr = tan((theta_r/2.0)*(M_PI/180.)) ;
  Xr = Xr*Xr ;
  Xr = Xr*2.0 ;
  Xc = tan((theta_c/2.0)*(M_PI/180.)) ;
  Xc = Xc*Xc ;
  Xc = Xc*2.0 ;
  cout << "theta_r = " << theta_r << ", theta_c = " << theta_c << endl ;
  //  cout << "Xr = " << Xr << ", Xc = " << Xc << endl ;
	   
  vector<vect3d> pos; //positions of nodes, change each iteration
  vector<Tri> trias;
  // Read in surf file
  if(from_surf) {
    read_surf(surf_file, pos, trias);
  } else {
    vector<surface_info> tmp_surf ;
    vector<int> tmp_map;
    readSurfaces(vog_file,tmp_surf,pos,tmp_map) ;
    int nsurf = tmp_surf.size() ;
    int ntri=0,nqua=0,ngen=0 ;
    for(int i=0;i<nsurf;++i) {
      ntri += tmp_surf[i].trias.size() ;
      nqua += tmp_surf[i].quads.size() ;
      ngen += tmp_surf[i].gen_faces.size() ;
    }
    if(ngen != 0) {
      cerr << "warning, general faces not supported in current geometry extraction tool" << endl ;
    }
    // create triangle surface mesh
    trias = vector<Tri>(ntri+nqua*2) ;
    int t= 0 ;
    for(int i=0;i<nsurf;++i) {
      for(size_t f=0;f<tmp_surf[i].trias.size();++f) {
	trias[t] = Tri(tmp_surf[i].trias[f][0]-1,
		       tmp_surf[i].trias[f][1]-1,
		       tmp_surf[i].trias[f][2]-1) ;
	t++ ;
      }
      for(size_t f=0;f<tmp_surf[i].quads.size();++f) {
	trias[t] = Tri(tmp_surf[i].quads[f][0]-1,
		       tmp_surf[i].quads[f][1]-1,
		       tmp_surf[i].quads[f][2]-1) ;
	t++ ;
	trias[t] = Tri(tmp_surf[i].quads[f][0]-1,
		       tmp_surf[i].quads[f][2]-1,
		       tmp_surf[i].quads[f][3]-1) ;
	t++ ;
      }
    }
  }



  // compute triangle normals and vertex angles
  cout << "read in " << trias.size() << " triangles" << endl ;
  computeTriangleProperties(trias,pos) ;

  // Create edge data structures
  vector<Edge> edges ;
  getEdges(trias, edges) ;
  cout << "number of edges = " << edges.size() << endl ;
  // create node2face structures
  vector<int> node2face_offsets ;
  vector<pair<int,int> > node2face ;
  int max_f2n = getNode2Face(trias,pos,node2face_offsets,node2face) ;

  // Compute eigenstructure at nodes to obtain principle direction
  // and null space
  vector<nodeInfo> ninfo(pos.size()) ;
  getNodeInfo(trias,node2face_offsets,node2face,max_f2n,ninfo) ;
  vector<edgeInfo> edge_points(edges.size()) ;
  // Reconstruct curves for edges
  edgeReconstruct(edges,pos,ninfo,trias,edge_points) ;
  vector<geomCoeff> trigeo(trias.size()) ;
  // reconstruct geometry for triangles using cubic bezier patches
  geoReconstruct(trias, edges, edge_points, pos,trigeo) ;

  // Write out a geometry file for adjustpos
  outputGeom(geo_file,pos,trias,trigeo) ;

  if(output_geom) {
    string geom_out = string(av[1]) + "_ref.surf" ;
    outputGeomSurf(geom_out,trigeo) ;
  }

  Loci::Finalize() ;
}




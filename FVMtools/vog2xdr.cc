#include <store.h>
#include <DStore.h>
#include <Map.h>
#include <DMap.h>
#include <multiMap.h>
#include <DMultiMap.h>
#include <constraint.h>
#include <distribute.h>
#include <distribute_container.h>
#include <distribute_io.h>
#include <parameter.h>
#include <fact_db.h>
#include <Loci_types.h>
#include <LociGridReaders.h>

#include <iostream>
#include <fstream>

#include <Tools/tools.h>
#include <map>

#include <Tools/xdr.h>

#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;
using std::ofstream ;
using std::ios ;

using Loci::fact_db ;
using Loci::debugout ;
namespace VOG {
  using Loci::entitySet ;
  using Loci::MapRepP ;
  using Loci::storeRepP ;
  using Loci::sequence ;
  using Loci::interval ;
  using Loci::create_intervalSet ;
  using Loci::UNIVERSE_MIN ;
  using Loci::UNIVERSE_MAX ;
  using Loci::store ;
  using Loci::dstore ;
  using Loci::param ;
  using Loci::Map ;
  using Loci::dMap ;
  using Loci::multiMap ;
  using Loci::constraint ;
  using Loci::vector3d ;
  using Loci::real_t ;
  using Loci::EMPTY ;
  using Loci::Entity ;
  using Loci::MPI_processes ;
  using Loci::MPI_rank ;
  using Loci::distributed_inverseMap ;
  using Loci::data_schema_traits ;
  
  template<class T> void readVector(hid_t group_id, const char *vector_name,
                                    vector<T> &v) {
    v.clear() ;
    hid_t dataset = H5Dopen(group_id,vector_name) ;
    if(dataset < 0) {
      H5Eclear() ;
      return ;
    }
    hid_t dspace = H5Dget_space(dataset) ;

    hsize_t size = 0 ;
    H5Sget_simple_extent_dims(dspace,&size,NULL) ;

    v.resize(size) ;

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }

  int getClusterNumFaces(unsigned char *cluster) {
    int num_faces = 0 ;
    while(*cluster != 0) {
      int npnts = *cluster ;
      cluster++ ;
      int nfaces = *cluster ;
      cluster++ ;
      num_faces += nfaces ;
      
      cluster += nfaces * (npnts + 2) ;
    }
    return num_faces ;
  }
  
  int fillClusterFaceSizes(unsigned char *cluster, int *sizes) {
    int num_faces = 0 ;
    while(*cluster != 0) {
      int npnts = *cluster ;
      cluster++ ;
      int nfaces = *cluster ;
      cluster++ ;
      for(int i=0;i<nfaces;++i)
        *sizes++ = npnts ;
      num_faces += nfaces ;
      
      cluster += nfaces * (npnts + 2) ;
    }
    return num_faces ;
  }

  unsigned char *readSignedVal(unsigned char *p, long &val) {
    unsigned char byte = *p++ ;

    int shift = 6 ;
    bool sign = (byte & 0x40)==0x40 ;
    val = byte & 0x3f ;
    while((byte & 0x80) == 0x80) {
      byte = *p++ ;
      int chunk = byte & 0x7f ;
      val += (chunk << shift) ;
      shift += 7 ;
    }
    if(sign)
      val = -val ;
    return p ;
  }
  
  unsigned char *readUnsignedVal(unsigned char *p, long &val) {
    unsigned char byte = *p++ ;
    int shift = 7 ;
    val = byte & 0x7f ;

    while((byte & 0x80) == 0x80) {
      byte = *p++ ;
      int chunk = byte & 0x7f ;
      val += (chunk << shift) ;
      shift += 7 ;
    }
    return p ;
  }

  unsigned char *readTable(unsigned char *p, long table[]) {
    int sz = *p++ ;
    // Get table size... note zero indicates maximum size
    if(sz == 0)
      sz = 256 ;

    // Get first table entry
    p = readSignedVal(p,table[0]) ;

    // Remaining entries are offsets
    for(int i=1;i<sz;++i) {
      long off ;
      p = readUnsignedVal(p,off) ;
      table[i] = table[i-1]+off ;
    }
    return p ;
  }    
      
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Loci::Map &cl, Loci::Map &cr, int face_base) {
    int num_faces = 0 ;
    while(*cluster != 0) {
      int npnts = *cluster ;
      cluster++ ;
      int nfaces = *cluster ;
      cluster++ ;
      for(int i=0;i<nfaces;++i) {
        int fc = face_base+num_faces+i ;
        for(int j=0;j<npnts;++j)
          face2node[fc][j] = *cluster++ ;
        cl[fc] = *cluster++ ;
        cr[fc] = *cluster++ ;
      }
      num_faces += nfaces ;
    }
    cluster++ ;
    // read in tables
    long nodeTable[256] ;
    cluster = readTable(cluster,nodeTable) ;
    long cellTable[256] ;
    cluster = readTable(cluster,cellTable) ;
    for(int i=0;i<num_faces;++i) {
      int fc = face_base+i ;
      int fsz = face2node[fc].size() ;
      for(int j=0;j<fsz;++j)
        face2node[fc][j] = nodeTable[face2node[fc][j]] ;
      cl[fc] = cellTable[cl[fc]] ;
      cr[fc] = cellTable[cr[fc]] ;
    }
    return num_faces ;
  }
}




int main(int ac, char *av[]) {
  using namespace Loci ;
  using namespace VOG ;
  Loci::Init(&ac,&av) ;
  string filename = av[1] ;
  filename += ".vog" ;

  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open vog file: " << filename << endl ;
    Loci::Abort() ;
  }

  hid_t face_g = H5Gopen(file_id,"face_info") ;

  vector<unsigned short> cluster_sizes ;
  vector<unsigned char> cluster_info ;
  readVector(face_g,"cluster_sizes",cluster_sizes) ;
  readVector(face_g,"cluster_info",cluster_info) ;
  vector<int> cluster_offset(cluster_sizes.size()+1) ;
  cluster_offset[0] = 0 ;
  for(size_t i=0;i<cluster_sizes.size();++i)
    cluster_offset[i+1] = cluster_offset[i] + cluster_sizes[i] ;

  int tot_faces = 0 ;
  for(size_t i=0;i<cluster_sizes.size();++i) {
    int nfaces = getClusterNumFaces(&cluster_info[cluster_offset[i]]) ;
    tot_faces += nfaces ;
  }

  entitySet faces = interval(0,tot_faces-1) ;
  store<int> counts ;
  counts.allocate(faces) ;
  tot_faces = 0 ;
  for(size_t i=0;i<cluster_sizes.size();++i) {
    int nfaces = fillClusterFaceSizes(&cluster_info[cluster_offset[i]],
                                      &counts[tot_faces]) ;
    tot_faces += nfaces ;
  }
  multiMap face2node ;
  face2node.allocate(counts) ;
  Loci::Map cl,cr ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;
  int face_base = 0 ;
  for(size_t i=0;i<cluster_sizes.size();++i) {
    int nfaces = fillFaceInfo(&cluster_info[cluster_offset[i]],
                              face2node,cl,cr,face_base) ;
    face_base += nfaces ;
  }

  entitySet cells = cl.image(faces) + cr.image(faces) ;
  cells &= interval(0,UNIVERSE_MAX) ;
  int ncells = cells.size() ;
  

  hid_t node_g = H5Gopen(file_id,"node_info") ;
  
  vector<vector3d<double> > pos ;
  readVector(node_g,"positions",pos) ;

  string outfile = av[1] ;
  outfile += ".xdr" ;
  FILE *FP = fopen(outfile.c_str(), "w") ;
  if(FP == NULL) {
    cerr << "can't open " << outfile <<  endl ;
    Loci::Abort() ;
  }
  XDR xdr_handle ;
  xdrstdio_create(&xdr_handle, FP, XDR_ENCODE) ;

  
  int npatch = 6 ;
  int nfaces = faces.size() ;
  int npnts = pos.size() ;
  int mxppf = 4 ;
  int mxfpc = 6 ;
  int ndim = 3 ;
  int nzones = 1 ;
  xdr_int(&xdr_handle, &ndim) ;
  xdr_int(&xdr_handle, &nzones) ;
  xdr_int(&xdr_handle, &npatch) ;
  xdr_int(&xdr_handle, &npnts) ;
  xdr_int(&xdr_handle, &nfaces) ;
  xdr_int(&xdr_handle, &ncells) ;
  xdr_int(&xdr_handle, &mxppf) ;
  xdr_int(&xdr_handle, &mxfpc) ;


  // Position Information
  for(int i=0;i<npnts;++i) {
    xdr_double(&xdr_handle,&pos[i].x) ;
    xdr_double(&xdr_handle,&pos[i].y) ;
    xdr_double(&xdr_handle,&pos[i].z) ;
  }
  int off = 0 ;
  for(entitySet::const_iterator i=faces.begin();i!=faces.end();++i) {
    xdr_int(&xdr_handle,&off) ;
    int clc = cl[*i]+1 ;
    int crc = cr[*i] ;
    if(crc >=0)
      crc+=1 ;
    xdr_int(&xdr_handle,&clc) ;
    xdr_int(&xdr_handle,&crc) ;
    off += face2node[*i].size() ;
  }
  xdr_int(&xdr_handle,&off) ;
    
  for(entitySet::const_iterator i=faces.begin();i!=faces.end();++i) {
    int sz = face2node[*i].size() ;
    for(int j=0;j<sz;++j)
      xdr_int(&xdr_handle,&face2node[*i][j]) ;
  }
  
  fclose(FP);
  xdr_destroy(&xdr_handle) ;

  

  return 0 ;
}

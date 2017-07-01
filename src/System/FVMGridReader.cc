//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#include "Config/conf.h"
#include <store.h>
#include <DStore.h>
#include <Map.h>
#include <DMap.h>
#include <multiMap.h>
#include <DMultiMap.h>
#include <constraint.h>
#include <distribute.h>
#include <distribute_container.h>
#include <parameter.h>
#include <fact_db.h>
#include <Loci_types.h>
#include <LociGridReaders.h>
#include "loci_globs.h"
#include <sstream>

#include <Tools/tools.h>
#include <map>
using std::map ;

#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
#include <algorithm>
using std::sort ;
using std::pair ;
using std::ostringstream ;
using std::istringstream ;

#ifdef HAS_MALLINFO
#include <malloc.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>

#include <parmetis.h>

#if REALTYPEWIDTH == 32
typedef float metisreal_t ;
#else
typedef double metisreal_t ;
#endif

namespace Loci {
  //#define MEMDIAG
  
#ifdef MEMDIAG
  void *memtop =0;

  class beginexec {
  public:
    beginexec() {
      memtop = sbrk(0) ;
    }
  } ;

  beginexec hackit;
#endif
  extern void ORBPartition(const vector<vector3d<float> > &pnts,
                           vector<int> &procid,
                           MPI_Comm comm) ;

  void memSpace(string s) {
#ifdef MEMDIAG

    unsigned long memsize = (char *)sbrk(0)-(char *)memtop ;
    debugout << s << ": malloc = " << double(memsize)/(1024.*1024) << endl ;
    //#define MEMINFO
#ifdef MEMINFO
    struct mallinfo info = mallinfo() ;
    debugout << s << ": minfo, arena=" << info.arena
             << ", ordblks=" << info.ordblks
             << ", hblks="<< info.hblks
             << ", hblkhd="<< info.hblkhd
             << ", uordblks=" << info.uordblks
             << ", fordblks="<<info.fordblks
             << ", keepcost="<<info.keepcost << endl ;
#endif
    debugout.flush() ;
#endif
  }
  extern bool use_simple_partition ;
  extern bool use_orb_partition ;
  extern bool load_cell_weights ;
  extern string cell_weight_file ;
  //Following assumption is needed for the next three functions
  //Assumption: We follow the convention that boundary cells are always on right side
  //of a face.

  template<class T> void readVectorDist(hid_t group_id,
                                        const char *vector_name,
                                        vector<long> sizes,
                                        vector<T> &v) {

    hid_t dataset = 0 ;
    hid_t dspace = 0 ;

    if(MPI_rank == 0) {
#ifdef H5_USE_16_API
      dataset = H5Dopen(group_id,vector_name) ;
#else
      dataset = H5Dopen(group_id,vector_name,H5P_DEFAULT) ;
#endif
      if(dataset < 0) {
        cerr << "unable to open dataset" << endl ;
        Loci::Abort() ;
      }
      dspace = H5Dget_space(dataset) ;
    }
    v.resize(sizes[MPI_rank]) ;

    if(MPI_rank == 0) { // read in vector from processor 0, send to other
      // processors
      long lsz = sizes[0] ;

      hsize_t dimension = lsz ;
      hsize_t count = dimension ;
      hsize_t stride = 1 ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif

      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count, NULL) ;
      start += dimension  ;

      int rank = 1 ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      Loci::DatatypeP dp = traits_type::get_type() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                          &v[0]) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
        FATAL(err < 0) ;
        Loci::Abort() ;
      }
      H5Sclose(memspace) ;

      // now read in remaining processor segments and send to corresponding
      // processor
      for(int i=1;i<MPI_processes;++i) {
        // read in remote processor data
        long sz = sizes[i] ;
        vector<T> tmp(sz) ;
        if(sz > 0) {
          dimension = sz ;
          count = dimension ;
          H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count, NULL) ;
          start += dimension  ;
          
          memspace = H5Screate_simple(rank,&dimension,NULL) ;
          hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                              &tmp[0]) ;
          if(err < 0) {
            cerr << "H5Dread() failed" << endl ;
            FATAL(err < 0) ;
            Loci::Abort() ;
          }
          H5Sclose(memspace) ;
        }
        // send to remote processor
        MPI_Send(&tmp[0],sz*sizeof(T),MPI_BYTE,1,0,MPI_COMM_WORLD) ;

      }
      H5Dclose(dataset) ;
      H5Sclose(dspace) ;

    } else {
      long size = sizes[MPI_rank] ;
      //      if(size > 0) {
      MPI_Status status ;
      MPI_Recv(&v[0],size*sizeof(T),MPI_BYTE,MPI_rank-1,0,MPI_COMM_WORLD,&status) ;
      //      }
      for(int i=MPI_rank+1;i<MPI_processes;++i) {
        long lsz = sizes[i] ;
        vector<T> tmp(lsz) ;
        MPI_Recv(&tmp[0],lsz*sizeof(T),MPI_BYTE,MPI_rank-1,0,MPI_COMM_WORLD,&status) ;
        MPI_Send(&tmp[0],lsz*sizeof(T),MPI_BYTE,MPI_rank+1,0,MPI_COMM_WORLD) ;
      }
    }
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
                   Map &cl, Map &cr, int face_base) {
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

  bool readBCfromVOG(string filename,
                 vector<pair<int,string> > &boundary_ids) {
    hid_t file_id = 0 ;
    int failure = 0 ; // No failure
    /* Save old error handler */
    H5E_auto_t old_func = 0 ;
    void *old_client_data = 0 ;
    if(MPI_rank == 0) {
#ifdef H5_USE_16_API
      H5Eget_auto(&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(NULL, NULL);
#else
      H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif


      file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
      if(file_id <= 0) 
        failure = 1 ;

      // Check to see if the file has surface info
#ifdef H5_USE_16_API
      hid_t bc_g = H5Gopen(file_id,"surface_info") ;
#else
      hid_t bc_g = H5Gopen(file_id,"surface_info",H5P_DEFAULT) ;
#endif

      boundary_ids.clear() ;
      // If surface info, then check surfaces
      if(bc_g > 0) {
        hsize_t num_bcs = 0 ;
        H5Gget_num_objs(bc_g,&num_bcs) ;
        for(hsize_t bc=0;bc<num_bcs;++bc) {
          char buf[1024] ;
          memset(buf, '\0', 1024) ;
          H5Gget_objname_by_idx(bc_g,bc,buf,sizeof(buf)) ;
          buf[1023]='\0' ;
          
          string name = string(buf) ;
#ifdef H5_USE_16_API
          hid_t sf_g = H5Gopen(bc_g,buf) ;
#else
          hid_t sf_g = H5Gopen(bc_g,buf,H5P_DEFAULT) ;
#endif
          hid_t id_a = H5Aopen_name(sf_g,"Ident") ;
          int ident ;
          H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
          H5Aclose(id_a) ;
          H5Gclose(sf_g) ;
          boundary_ids.push_back(pair<int,string>(ident,name)) ;
        }
        H5Gclose(bc_g) ;
        H5Fclose(file_id) ;
        /* Restore previous error handler */
#ifdef H5_USE_16_API
        H5Eset_auto(old_func, old_client_data);
#else
        H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
#endif
      }
    }
    // Share boundary tag data with all other processors
    int bsz = boundary_ids.size() ;
    MPI_Bcast(&bsz,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(bsz > 0) {
      string buf ;
      if(MPI_rank == 0) {
        ostringstream oss ;
        
        for(int i=0;i<bsz;++i)
          oss << boundary_ids[i].first << ' ' << boundary_ids[i].second << ' ';
        buf = oss.str() ;
      }
      int bufsz = buf.size() ;
      MPI_Bcast(&bufsz,1,MPI_INT,0,MPI_COMM_WORLD) ;
      char *data = new char[bufsz+1] ;
      if(MPI_rank == 0)
        strcpy(data,buf.c_str()) ;
      MPI_Bcast(data,bufsz,MPI_CHAR,0,MPI_COMM_WORLD) ;
      buf = string(data) ;
      istringstream iss(buf) ;
      delete[] data ;
      boundary_ids.clear() ;
      for(int i=0;i<bsz;++i) {
        int id ;
        string name ;
        iss >> id >> name ;
        boundary_ids.push_back(pair<int,string>(id,name)) ;
      }
    }
    if(failure)
      return false ;
    return true ;
  }
  
  unsigned long readAttributeLong(hid_t group, const char *name) {
    hid_t id_a = H5Aopen_name(group,name) ;
    unsigned long val = 0;
    H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
    H5Aclose(id_a) ;
    return val ;
  }
  
  bool readVolTags(hid_t input_fid,
                   vector<pair<string,Loci::entitySet> > &volDat) {
    using namespace Loci ;
    /* Save old error handler */
    H5E_auto_t old_func = 0;
    void *old_client_data = 0 ;
#ifdef H5_USE_16_API
    H5Eget_auto(&old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(NULL, NULL);
#else
    H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif
    
    vector<pair<string,entitySet> > volTags ;
#ifdef H5_USE_16_API
    hid_t cell_info = H5Gopen(input_fid,"cell_info") ;
#else
    hid_t cell_info = H5Gopen(input_fid,"cell_info",H5P_DEFAULT) ;
#endif
    if(cell_info > 0) {
      vector<string> vol_tag ;
      vector<entitySet> vol_set ;
      vector<int> vol_id ;
      
      hsize_t num_tags = 0 ;
      H5Gget_num_objs(cell_info,&num_tags) ;
      for(hsize_t tg=0;tg<num_tags;++tg) {
        char buf[1024] ;
        memset(buf, '\0', 1024) ;
        H5Gget_objname_by_idx(cell_info,tg,buf,sizeof(buf)) ;
        buf[1023]='\0' ;
        
        string name = string(buf) ;
#ifdef H5_USE_16_API
        hid_t vt_g = H5Gopen(cell_info,buf) ;
#else
        hid_t vt_g = H5Gopen(cell_info,buf,H5P_DEFAULT) ;
#endif
        hid_t id_a = H5Aopen_name(vt_g,"Ident") ;
        int ident ;
        H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
        H5Aclose(id_a) ;
        entitySet dom ;
        HDF5_ReadDomain(vt_g,dom) ;
        vol_tag.push_back(name) ;
        vol_set.push_back(dom) ;
        vol_id.push_back(ident) ;
        H5Gclose(vt_g) ;
      }
      int maxi =0 ;
      for(size_t i=0;i<vol_id.size();++i)
        maxi = max(vol_id[i],maxi) ;
      vector<pair<string,entitySet> > tmp(maxi+1) ;
      volTags.swap(tmp) ;
      for(size_t i=0;i<vol_id.size();++i) {
        volTags[vol_id[i]].first = vol_tag[i] ;
        volTags[vol_id[i]].second = vol_set[i] ;
      }
    } else {
#ifdef H5_USE_16_API
      hid_t file_info = H5Gopen(input_fid,"file_info") ;
#else
      hid_t file_info = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
#endif
      long numCells = readAttributeLong(file_info,"numCells") ;
      volTags.push_back(pair<string,entitySet>
                        (string("Main"),
                         entitySet(interval(0,numCells-1)))) ;
      H5Gclose(file_info) ;
    }
  
#ifdef H5_USE_16_API
    /* Restore previous error handler */
    H5Eset_auto(old_func, old_client_data);
#else
    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
#endif
    volDat.swap(volTags) ;
    return true ;
  }

  //Description: Reads grid structures from grid file in the .vog format.
  //Input: file name and max_alloc (starting of entity assignment - node base)
  //Output:
  // local_cells, local_nodes, local_faces: partition of nodes, faces, cells
  // pos: position of nodes
  // cl: Mapping from face to cell on the left side
  // cr: Mapping from face to cell on the right side
  // face2node: MultiMapping from a face to nodes
  bool readGridVOG(vector<entitySet> &local_nodes,
		   vector<entitySet> &local_faces,
		   vector<entitySet> &local_cells,
		   store<vector3d<double> > &pos, Map &cl, Map &cr,
		   multiMap &face2node, 
		   store<string> &boundary_names,
		   store<string> &boundary_tags,
                   vector<pair<string,entitySet> > &volTags,
		   int max_alloc, string filename) {

    local_nodes.resize(Loci::MPI_processes);
    local_faces.resize(Loci::MPI_processes);
    local_cells.resize(Loci::MPI_processes);

    hid_t file_id = 0 ;
    hid_t face_g = 0 ;
    hid_t node_g = 0 ;
    hid_t dataset = 0 ;
    hid_t dspace = 0 ;
    long nnodes = 0 ;
    int failure = 0 ; // No failure
    /* Save old error handler */
    H5E_auto_t old_func = 0 ;
    void *old_client_data = 0 ;
    vector<pair<int,string> > boundary_ids ;
    if(MPI_rank == 0) {
#ifdef H5_USE_16_API
      H5Eget_auto(&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(NULL, NULL);
#else
      H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif

      file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
      if(file_id <= 0) 
        failure = 1 ;
      
#ifdef H5_USE_16_API
      face_g = H5Gopen(file_id,"face_info") ;
      node_g = H5Gopen(file_id,"node_info") ;
      // Check to see if the file has surface info
      hid_t bc_g = H5Gopen(file_id,"surface_info") ;
#else
      face_g = H5Gopen(file_id,"face_info",H5P_DEFAULT) ;
      node_g = H5Gopen(file_id,"node_info",H5P_DEFAULT) ;
      // Check to see if the file has surface info
      hid_t bc_g = H5Gopen(file_id,"surface_info",H5P_DEFAULT) ;
#endif



      // If surface info, then check surfaces
      if(bc_g > 0) {
        hsize_t num_bcs = 0 ;
        H5Gget_num_objs(bc_g,&num_bcs) ;
        for(hsize_t bc=0;bc<num_bcs;++bc) {
          char buf[1024] ;
          memset(buf, '\0', 1024) ;
          H5Gget_objname_by_idx(bc_g,bc,buf,sizeof(buf)) ;
          buf[1023]='\0' ;
          
          string name = string(buf) ;
#ifdef H5_USE_16_API
          hid_t sf_g = H5Gopen(bc_g,buf) ;
#else
          hid_t sf_g = H5Gopen(bc_g,buf,H5P_DEFAULT) ;
#endif
          hid_t id_a = H5Aopen_name(sf_g,"Ident") ;
          int ident ;
          H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
          H5Aclose(id_a) ;
          H5Gclose(sf_g) ;
          boundary_ids.push_back(pair<int,string>(ident,name)) ;
        }
        H5Gclose(bc_g) ;
      }

      // First read in and disribute node positions...
#ifdef H5_USE_16_API
      dataset = H5Dopen(node_g,"positions") ;
#else
      dataset = H5Dopen(node_g,"positions",H5P_DEFAULT) ;
#endif
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = 1 ;
        
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      nnodes = size ;

    }
    int fail_state = 0 ;
    MPI_Allreduce(&failure,&fail_state,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(fail_state != 0)
      return false ;
    // Share boundary tag data with all other processors
    int bsz = boundary_ids.size() ;
    MPI_Bcast(&bsz,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(bsz > 0) {
      string buf ;
      if(MPI_rank == 0) {
        ostringstream oss ;
        
        for(int i=0;i<bsz;++i)
          oss << boundary_ids[i].first << ' ' << boundary_ids[i].second << ' ';
        buf = oss.str() ;
      }
      int bufsz = buf.size() ;
      MPI_Bcast(&bufsz,1,MPI_INT,0,MPI_COMM_WORLD) ;
      char *data = new char[bufsz+1] ;
      if(MPI_rank == 0)
        strcpy(data,buf.c_str()) ;
      MPI_Bcast(data,bufsz+1,MPI_CHAR,0,MPI_COMM_WORLD) ;
      buf = string(data) ;
      istringstream iss(buf) ;
      delete[] data ;
      boundary_ids.clear() ;
      for(int i=0;i<bsz;++i) {
        int id ;
        string name ;
        iss >> id >> name ;
        boundary_ids.push_back(pair<int,string>(id,name)) ;
      }
    }
      
    MPI_Bcast(&nnodes,1,MPI_LONG,0,MPI_COMM_WORLD) ;

    // create node allocation
    long npnts = nnodes ;
    int node_ivl = npnts / Loci::MPI_processes;
    int node_ivl_rem = npnts % Loci::MPI_processes ;
    int node_accum = 0 ;
    int nodes_base = max_alloc ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int node_accum_update = node_accum + node_ivl + ((i<node_ivl_rem)?1:0) ;
      if(i == Loci::MPI_processes-1) {
	local_nodes[i] = interval(nodes_base + node_accum,
                                  nodes_base + npnts - 1) ;
      } else {
	local_nodes[i] = interval(nodes_base + node_accum,
                                  nodes_base + node_accum_update - 1) ;
      }
      node_accum = node_accum_update ;
    }

    pos.allocate(local_nodes[MPI_rank]) ;
    if(MPI_rank == 0) { // read in node positions, send to other processors
      memSpace("Read in Pos") ;
      // read processor zero section first
      int lst = local_nodes[MPI_rank].Min() ;
      int lsz = local_nodes[0].size() ;

      hsize_t dimension = lsz ;
      hsize_t count = dimension ;
      hsize_t stride = 1 ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif

      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count, NULL) ;
      start += dimension  ;

      int rank = 1 ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
      typedef data_schema_traits<vector3d<double> > traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                          &pos[lst]) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
        FATAL(err < 0) ;
        Loci::Abort() ;
      }
      H5Sclose(memspace) ;

      // now read in remaining processor segments and send to corresponding
      // processor
      size_t mxsz = local_nodes[0].size() ;
      for(int i=0;i<MPI_processes;++i)
        mxsz = max(size_t(local_nodes[i].size()),mxsz) ;
      
      for(int i=1;i<MPI_processes;++i) {
        // read in remote processor data
        int sz = local_nodes[i].size() ;
        if(sz == 0) {
          cerr << "sending a zero sized block" << endl ;
        }
        vector<vector3d<double> > tpos(sz) ;
        
        dimension = sz ;
        count = dimension ;
        H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count, NULL) ;
        start += dimension  ;

        memspace = H5Screate_simple(rank,&dimension,NULL) ;
        hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                            &tpos[0]) ;
        if(err < 0) {
          cerr << "H5Dread() failed" << endl ;
          FATAL(err < 0) ;
          Loci::Abort() ;
        }
        H5Sclose(memspace) ;

        // send to remote processor
        MPI_Send(&tpos[0],sz*3,MPI_DOUBLE,1,0,MPI_COMM_WORLD) ;
      }
      
      H5Sclose(dspace) ;
      H5Dclose(dataset) ;
      H5Gclose(node_g) ;
    } else {
      // Receive nodes from root processor
      FATAL(local_nodes[MPI_rank].num_intervals()!=1) ;
      int start = local_nodes[MPI_rank].Min() ;
      int size = local_nodes[MPI_rank].size() ;
      MPI_Status status ;
      MPI_Recv(&pos[start],size*3,MPI_DOUBLE,MPI_rank-1,0,MPI_COMM_WORLD,&status) ;
      int mxsz = 0 ;
      for(int i=MPI_rank+1;i<MPI_processes;++i) {
        int lsz = local_nodes[i].size() ;
        mxsz = max(mxsz,lsz) ;
      }

      vector<vector3d<double> > tpos(mxsz) ;
      // Shift remaining ones into place
      for(int i=MPI_rank+1;i<MPI_processes;++i) {
        int lsz = local_nodes[i].size() ;
        MPI_Recv(&tpos[0],lsz*3,MPI_DOUBLE,MPI_rank-1,0,MPI_COMM_WORLD,&status) ;
        int count = 0 ;
        MPI_Get_count(&status,MPI_DOUBLE,&count) ;
        if(count != lsz*3) {
          cerr << "processor" << MPI_rank << " recieved " << count <<
            " words but was expecting " << lsz*3 << endl ;
        }
        MPI_Send(&tpos[0],lsz*3,MPI_DOUBLE,MPI_rank+1,0,MPI_COMM_WORLD) ;
      }
    }

    vector<unsigned char> cluster_info ;
    vector<unsigned short> cluster_sizes ;
    // Now read in face clusters
    long nclusters = 0 ;
    if(MPI_rank == 0) {
#ifdef H5_USE_16_API
      dataset = H5Dopen(face_g,"cluster_sizes") ;
#else
      dataset = H5Dopen(face_g,"cluster_sizes",H5P_DEFAULT) ;
#endif
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = 1 ;
      
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      nclusters = size ;
    }
    MPI_Bcast(&nclusters,1,MPI_LONG,0,MPI_COMM_WORLD) ;

    vector<long> cluster_dist(MPI_processes,0) ;
    long sum = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      cluster_dist[i] = nclusters/MPI_processes +
        ((nclusters%MPI_processes)>i?1:0);
      sum += cluster_dist[i] ;
    }
    FATAL(sum != nclusters) ;
    memSpace("before face cluster size reading" ) ;
    readVectorDist(face_g,"cluster_sizes",cluster_dist,cluster_sizes) ;

    long cluster_info_size = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i)
      cluster_info_size += cluster_sizes[i] ;


    MPI_Allgather(&cluster_info_size,sizeof(long),MPI_BYTE,
                  &cluster_dist[0],sizeof(long),MPI_BYTE,
                  MPI_COMM_WORLD) ;

    memSpace("before face cluster reading" ) ;
    readVectorDist(face_g,"cluster_info",cluster_dist,cluster_info) ;

    memSpace("after face cluster reading" ) ;

    // Read in volume tag information
    vector<pair<string,Loci::entitySet> > volDat ;
    if(MPI_rank == 0) {
      readVolTags(file_id,volDat) ;
    }
    int nvtags = volDat.size() ;
    MPI_Bcast(&nvtags,1,MPI_INT,0,MPI_COMM_WORLD) ;
    for(int i=0;i<nvtags;++i) {
      int sz = 0 ;
      if(MPI_rank == 0) sz = volDat[i].first.size() ;
      MPI_Bcast(&sz,1,MPI_INT,0,MPI_COMM_WORLD) ;
      char *buf = new char[sz+1] ;
      buf[sz] = '\0' ;
      if(MPI_rank == 0) strcpy(buf,volDat[i].first.c_str()) ;
      MPI_Bcast(buf,sz,MPI_CHAR,0,MPI_COMM_WORLD) ;
      string name = string(buf) ;
      delete[] buf ;
      int nivals = 0 ;
      if(MPI_rank == 0) nivals = volDat[i].second.num_intervals() ;
      MPI_Bcast(&nivals,1,MPI_INT,0,MPI_COMM_WORLD) ;

      int *ibuf = new int[nivals*2] ;
      if(MPI_rank == 0) 
        for(int j=0;j<nivals;++j) {
          ibuf[j*2]= volDat[i].second[j].first ;
          ibuf[j*2+1]= volDat[i].second[j].second ;
        }
      MPI_Bcast(ibuf,nivals*2,MPI_INT,0,MPI_COMM_WORLD) ;
      entitySet set ;
      for(int j=0;j<nivals;++j) 
        set += interval(ibuf[j*2],ibuf[j*2+1]) ;
      if(MPI_rank != 0) 
        volDat.push_back(pair<string,entitySet>(name,set)) ;
      delete[] ibuf ;
    }

    volTags.swap(volDat) ;
    if(MPI_rank == 0) {
      H5Gclose(face_g) ;
      H5Fclose(file_id) ;
      /* Restore previous error handler */
#ifdef H5_USE_16_API
      H5Eset_auto(old_func, old_client_data);
#else
      H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
#endif
    }
    MPI_Allreduce(&failure,&fail_state,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(fail_state != 0)
      return false ;
    
    memSpace("before unpacking clusters") ;
    vector<long> cluster_offset(cluster_sizes.size()+1) ;
    cluster_offset[0] = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i)
      cluster_offset[i+1] = cluster_offset[i] + cluster_sizes[i] ;

    int tot_faces = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      int nfaces = getClusterNumFaces(&cluster_info[cluster_offset[i]]) ;
      tot_faces += nfaces ;
    }

    // Now get a face allocation for each processor
    vector<int> faces_pp(MPI_processes,0) ;
    MPI_Allgather(&tot_faces,1,MPI_INT,&faces_pp[0],1,MPI_INT,
                  MPI_COMM_WORLD) ;

    int face_accum = 0 ;
    int faces_base = max_alloc + nnodes ;
    for(int i = 0; i < MPI_processes; ++i) {
      local_faces[i] = EMPTY ;
      if(faces_pp[i] > 0)
	local_faces[i] = interval(faces_base + face_accum,
				  faces_base + face_accum+faces_pp[i]- 1) ;
      face_accum += faces_pp[i] ;
    }
    int cells_base = faces_base + face_accum ;


    store<int> counts ;
    counts.allocate(local_faces[MPI_rank]) ;
    tot_faces = local_faces[MPI_rank].Min() ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      int nfaces = fillClusterFaceSizes(&cluster_info[cluster_offset[i]],
                                        &counts[tot_faces]) ;
      tot_faces += nfaces ;
    }
    face2node.allocate(counts) ;
    cl.allocate(local_faces[MPI_rank]) ;
    cr.allocate(local_faces[MPI_rank]) ;
    int face_base = local_faces[MPI_rank].Min() ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      int nfaces = fillFaceInfo(&cluster_info[cluster_offset[i]],
                                face2node,cl,cr,face_base) ;
      face_base += nfaces ;
    }


    for(size_t i=0;i<volTags.size();++i)
      volTags[i].second = (volTags[i].second >> cells_base) ;
    entitySet boundary_faces ;
    entitySet boundary_taglist ;
    int max_cell = 0 ;
    FORALL(local_faces[MPI_rank],fc) {
      int fsz = face2node[fc].size() ;
      for(int i=0;i<fsz;++i)
        face2node[fc][i] += nodes_base ;
      max_cell = max(max_cell,cl[fc]) ;
      max_cell = max(max_cell,cr[fc]) ;
      cl[fc] += cells_base ;
      if(cr[fc] >= 0)
        cr[fc] += cells_base ;
      else {
	boundary_faces += fc ;
	boundary_taglist += cr[fc] ;
      }
    } ENDFORALL ;

    int global_max_cell = max_cell ;
    MPI_Allreduce(&max_cell,&global_max_cell,1,MPI_INT,MPI_MAX,
                  MPI_COMM_WORLD) ;

    int ncells = global_max_cell+1 ;
    int cell_ivl = ncells / MPI_processes;
    int cell_ivl_rem = ncells % MPI_processes ;
    int cell_accum = 0 ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int cell_accum_update = cell_accum + cell_ivl + ((i<cell_ivl_rem)?1:0) ;

      if(i == MPI_processes-1) {
	local_cells[i] = interval(cells_base + cell_accum,
				  cells_base + ncells-1) ;
      } else {
	local_cells[i] = interval(cells_base + cell_accum,
                                  cells_base + cell_accum_update - 1) ;
      }
      cell_accum = cell_accum_update ;
    }

    // Now fixup boundary cells
    boundary_taglist = all_collect_entitySet(boundary_taglist) ;
    int ref_base = cells_base+ncells ;
    entitySet refSet = interval(ref_base,ref_base+(boundary_taglist.size()-1)) ;
    
    store<int> boundary_info ;
    boundary_info.allocate(refSet) ;
    map<int,int> boundary_remap ;
    int cnt = 0 ;
    FORALL(boundary_taglist,ii) {
      boundary_remap[ii]=cnt+ref_base ;
      boundary_info[cnt+ref_base]=ii ;
      cnt++ ;
    } ENDFORALL ;
    
    FORALL(boundary_faces,fc) {
      cr[fc] = boundary_remap[cr[fc]] ;
    } ENDFORALL ;

    boundary_names.allocate(refSet) ;
    boundary_tags.allocate(refSet) ;

    map<int,int> bcmap ;
    FORALL(refSet, ii) {
      int bc = boundary_info[ii] ;
      char buf[128] ;
      bzero(buf,128) ;
      snprintf(buf,127,"BC_%d",-bc) ;
      bcmap[-bc] = ii ;
      boundary_tags[ii] = string(buf) ;
      boundary_names[ii] = string(buf) ;
    } ENDFORALL ;

    for(size_t i=0;i<boundary_ids.size();++i) {
      int id = boundary_ids[i].first ;
      map<int,int>::const_iterator mi = bcmap.find(id) ;
      if(mi != bcmap.end())
	boundary_names[mi->second] = boundary_ids[i].second ;
      else 
	debugout << "id " << id << " for boundary surface " << boundary_ids[i].second
	     << " not found!" << endl ;
    }

      
    if(Loci::MPI_rank == 0) {
      Loci::debugout << " boundaries identified as:" ;
      FORALL(refSet, bc) {
	debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      Loci::debugout << endl ;
    }

    memSpace("after unpacking clusters") ;
    return true ;
  }

  extern void distributed_inverseMap(multiMap &result,
                                     vector<pair<Entity,Entity> > &input,
                                     entitySet input_image,
                                     entitySet input_preimage,
                                     const std::vector<entitySet> &init_ptn) ;

  vector<entitySet> newMetisPartitionOfCells(const vector<entitySet> &local_cells,
                                             const Map &cl, const Map &cr,
					     const store<string> &boundary_tags) {


    entitySet dom = cl.domain() & cr.domain() ;
    entitySet refset = boundary_tags.domain() ;
    entitySet bcfaces = cr.preimage(refset).first ;
    dom -= bcfaces ;
    int cnt = dom.size() ;
    
    vector<pair<int,int> > rawMap(cnt*2) ;
    int j = 0 ;
    entitySet::const_iterator ei ;
    for(ei=dom.begin();ei!=dom.end();++ei) {
      rawMap[j++] = pair<int,int>(cl[*ei],cr[*ei]) ;
      rawMap[j++] = pair<int,int>(cr[*ei],cl[*ei]) ;
    }

    sort(rawMap.begin(),rawMap.end()) ;

    multiMap cell2cell ;
    entitySet all_cells ;
    for(int i=0;i<MPI_processes;++i)
      all_cells += local_cells[i] ;

    distributed_inverseMap(cell2cell,rawMap, all_cells,all_cells,local_cells) ;

    vector<pair<int,int> >().swap(rawMap) ; // Free up memory from rawMap
    int count = 0 ;
    idx_t size_map = local_cells[Loci::MPI_rank].size() ;
    vector<idx_t> size_adj(size_map) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_adj[count] = cell2cell[*ei].size() ;
      ++count ;
    }

    vector<idx_t> part(size_map) ;
    vector<idx_t> xadj(size_map+1) ;
    idx_t edgecut ;
    vector<idx_t> vdist(Loci::MPI_processes + 1) ;
    int cmin = local_cells[0].Min();
    for(int i = 0; i < Loci::MPI_processes; i++) {
      cmin = min(local_cells[i].Min(), cmin);
    }
    
    // check local_cells to be consistent with parmetis partitioning
    for(int i=0;i<Loci::MPI_processes;++i) {
      if(local_cells[i].size() == 0 ||
         (int)local_cells[i].size() != local_cells[i].Max()-local_cells[i].Min()+1) {
        cerr << "invalid local cell set, p=" << i
             << ", local_cells=" << local_cells[i] << endl ;
      }
      if(i>0 && local_cells[i-1].Max()+1 != local_cells[i].Min()) {
        cerr << "gap between processors in cell numbering" << endl ;
      }
    }
    
      
    
    edgecut = 0 ;
    xadj[0] = 0 ;
    for(idx_t i = 0; i < size_map; ++i)
      xadj[i+1] = xadj[i] + size_adj[i] ;

    idx_t min_size = size_adj[0] ;
    for(idx_t i = 0; i < size_map; ++i)
      min_size = min(min_size,size_adj[i]) ;
    if(min_size == 0) 
      cerr << "cell with no adjacency!" << endl ;
    
    idx_t tot = xadj[size_map] ;
    vector<idx_t> adjncy(tot) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_t sz = cell2cell[*ei].size() ;
      for(size_t i = 0; i != sz; ++i)        {
	adjncy[count] = cell2cell[*ei][i] - cmin ;
	count ++ ;
      }
    }
    cell2cell.setRep(multiMap().Rep()) ;// Free up memory from multiMap

    vdist[0] = 0 ;
    for(int i = 1; i <= Loci::MPI_processes; ++i)
      vdist[i] = vdist[i-1] + local_cells[i-1].size() ;
    idx_t top = vdist[Loci::MPI_processes] ;
    
    bool trouble = false ;
    for(idx_t i=0;i<tot;++i)
      if(adjncy[i] >= top)
        trouble = true ;
    if(trouble)
      cerr << "adjacency list contains out of bounds reference" << endl ;
    
    MPI_Comm mc = MPI_COMM_WORLD ;
    idx_t nparts = Loci::MPI_processes ; // number of partitions
    idx_t wgtflag = 0 ;
    idx_t numflag = 0 ;
    idx_t options = 0 ;

    // read in additional vertex weights if any
    if(load_cell_weights) {
      // check if the file exists
      int file_exists = 1 ;
      if(Loci::MPI_rank == 0) {
        struct stat buf ;
        if(stat(cell_weight_file.c_str(),&buf) == -1 ||
           !S_ISREG(buf.st_mode))
          file_exists = 0 ;
      }
      MPI_Bcast(&file_exists,1,MPI_INT,0,MPI_COMM_WORLD) ;
      
      if(file_exists == 1) {
        if(Loci::MPI_rank == 0) {
          std::cout << "ParMETIS reading additional cell weights from: "
                    << cell_weight_file << std::endl ;
        }
        
        // create a hdf5 handle
        hid_t file_id = Loci::hdf5OpenFile(cell_weight_file.c_str(),
                                           H5F_ACC_RDONLY, H5P_DEFAULT) ;
        if(file_id < 0) {
          std::cerr << "...file reading failed..., Aborting" << std::endl ;
          Loci::Abort() ;
        }
        
        // read
        entitySet dom = local_cells[Loci::MPI_rank] ;
        store<int> cell_weights ;

        readContainerRAW(file_id,"cell weight", cell_weights.Rep(),
                         MPI_COMM_WORLD) ;

        if(cell_weights.domain() != local_cells[Loci::MPI_rank]) {
          cerr << "cell weights partition inconsistent!" << endl ;
          Loci::Abort() ;
        }
        
        Loci::hdf5CloseFile(file_id) ;

        // compute necessary ParMETIS data-structure
        wgtflag = 2 ;           // weights on the vertices only
        idx_t ncon = 2 ;          // number of weights per vertex
        idx_t tpwgts_len = ncon*nparts ;
        vector<metisreal_t> tpwgts(tpwgts_len) ;

        for(idx_t i=0;i<tpwgts_len;++i)
          tpwgts[i] = 1.0 / double(nparts) ;
        
        vector<metisreal_t> ubvec(ncon) ;
        for(idx_t i=0;i<ncon;++i)
          ubvec[i] = 1.05 ;     // as recommended by the ParMETIS manual

        // now construct the vertex weights
        vector<idx_t> vwgt(ncon*size_map) ;
        int cnt = 0 ;
        for(entitySet::const_iterator
              ei=local_cells[Loci::MPI_rank].begin();
            ei!=local_cells[Loci::MPI_rank].end();++ei,cnt+=ncon) {
          // first weight for cell is 1 (the cell computation)
          vwgt[cnt] = 1 ;
          // the second weight is from the store cell_weights[*ei]
          vwgt[cnt+1] = cell_weights[*ei] ;
        }

        // now call the ParMETIS routine (V3)
        ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],&vwgt[0],NULL,
                             &wgtflag,&numflag,&ncon,&nparts,
                             &tpwgts[0],&ubvec[0],&options,&edgecut,&part[0],
                             &mc) ;

          
      } else {
        // if weight file does not exist, then we would
        // fall back to the non weighted partition
        if(Loci::MPI_rank == 0) {
          std::cout << "ParMETIS cell weight file not found, "
                    << "using non-weighted partition..." << std::endl ;
        }
        idx_t ncon = 1 ;
        idx_t tpwgts_len = ncon*nparts ;
        vector<metisreal_t> tpwgts(tpwgts_len) ;
        for(idx_t i=0;i<tpwgts_len;++i)
          tpwgts[i] = 1.0 / double(nparts) ;
        
        
        metisreal_t ubvec = 1.05 ;
        wgtflag = 0 ;
        ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],NULL,NULL,
                             &wgtflag,&numflag,&ncon,&nparts,
                             &tpwgts[0],&ubvec,&options,&edgecut,
                             &part[0],&mc) ;
      }
      
    } else {
      idx_t ncon = 1 ;
      idx_t tpwgts_len = ncon*nparts ;
      vector<metisreal_t> tpwgts(tpwgts_len) ;
      for(idx_t i=0;i<tpwgts_len;++i)
        tpwgts[i] = 1.0 / double(nparts) ;

      metisreal_t ubvec = 1.05 ;
      wgtflag = 0 ;
      ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],NULL,NULL,
                           &wgtflag,&numflag,&ncon,&nparts,
                           &tpwgts[0],&ubvec,&options,&edgecut,
                           &part[0],&mc) ;
    }

    if(Loci::MPI_rank == 0)
      Loci::debugout << " Parmetis Edge cut   " <<  edgecut << endl ;

    //find the partition ptn given by Metis
    vector<entitySet> ptn ;

    for(int i = 0; i < Loci::MPI_processes; ++i)
      ptn.push_back(EMPTY) ;
    cmin = local_cells[Loci::MPI_rank].Min() ;
    for(int i=0;i<size_map;++i) {
      ptn[part[i]] += i + cmin ;
    }
    return ptn;

  }


  void redistribute_container(const vector<entitySet> &ptn,
                              const vector<entitySet> &ptn_t,
                              entitySet new_alloc,
                              storeRepP inRep,storeRepP outRep) {
    vector<sequence> rdom(MPI_processes) ;
    vector<int> send_sizes(MPI_processes) ;
    // Determine how to redistribute current domain to new processors

    entitySet::const_iterator ei = new_alloc.begin() ;
    for(int i=0;i<MPI_processes;++i) {
      send_sizes[i] = inRep->pack_size(ptn[i]) ;
      sequence s ;
      for(entitySet::const_iterator si=ptn_t[i].begin();si!=ptn_t[i].end();++si) {
        s += *ei ;
        ++ei ;
      }
      rdom[i] = s ;
    }
    WARN(ei != new_alloc.end()) ;

    vector<int> recv_sizes(MPI_processes) ;
    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sizes[i] ;
      size_recv += recv_sizes[i] ;
    }
    //    outRep->allocate(new_alloc) ;
    unsigned char *send_store = new unsigned char[size_send] ;
    unsigned char *recv_store = new unsigned char[size_recv] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sizes[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sizes[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      inRep->pack(send_store, loc_pack, size_send, ptn[i]) ;


    MPI_Alltoallv(send_store,&send_sizes[0], send_displacement , MPI_PACKED,
		  recv_store, &recv_sizes[0], recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      outRep->unpack(recv_store, loc_pack, size_recv, rdom[i]) ;
    }

    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;
  }

  inline bool fieldSort(const std::pair<Entity,Entity> &p1,
                        const std::pair<Entity,Entity> &p2) {
    return p1.first < p2.first ;
  }

  void remapGrid(vector<entitySet> &node_ptn,
                 vector<entitySet> &face_ptn,
                 vector<entitySet> &cell_ptn,
                 vector<entitySet> &node_ptn_t,
                 vector<entitySet> &face_ptn_t,
                 vector<entitySet> &cell_ptn_t,
                 store<vector3d<double> > &t_pos, Map &tmp_cl,
                 Map &tmp_cr, multiMap &tmp_face2node,
		 vector<entitySet> &bcsurf_ptn,
		 store<string> &tmp_boundary_names,
		 store<string> &tmp_boundary_tags,
                 entitySet nodes, entitySet faces, entitySet cells,
                 store<vector3d<double> > &pos, Map &cl, Map &cr,
                 multiMap &face2node,
		 store<string> &boundary_names, 
		 store<string> &boundary_tags, 
		 entitySet bcsurfset,
                 fact_db &facts) {

    pos.allocate(nodes) ;
    cl.allocate(faces) ;
    cr.allocate(faces) ;
    entitySet old_nodes = t_pos.domain() ;
    redistribute_container(node_ptn,node_ptn_t,nodes,t_pos.Rep(),pos.Rep()) ;
    t_pos.allocate(EMPTY) ;
    redistribute_container(face_ptn,face_ptn_t,faces,tmp_cr.Rep(),cr.Rep()) ;
    tmp_cr.allocate(EMPTY) ;
    redistribute_container(face_ptn,face_ptn_t,faces,tmp_cl.Rep(),cl.Rep()) ;
    tmp_cl.allocate(EMPTY) ;

    using std::pair ;
    vector<pair<Entity,Entity> > sortlist(faces.size()) ;

    store<int> count ;
    entitySet infaces = tmp_face2node.domain() ;
    count.allocate(infaces) ;
    for(entitySet::const_iterator ii=infaces.begin();ii!=infaces.end();++ii)
      count[*ii] = tmp_face2node.end(*ii)-tmp_face2node.begin(*ii) ;
    store<int> count_reorder ;
    count_reorder.allocate(faces) ;
    redistribute_container(face_ptn,face_ptn_t,faces,count.Rep(),count_reorder.Rep()) ;

    face2node.allocate(count_reorder) ;
    redistribute_container(face_ptn,face_ptn_t,faces,tmp_face2node.Rep(),
                           face2node.Rep()) ;
    tmp_face2node.allocate(EMPTY) ;

    // sort faces
    int i=0 ;
    FORALL(faces,fc) {
      Entity minc = max(cr[fc],cl[fc]) ;
      sortlist[i++] = pair<Entity,Entity>(minc,fc) ;
    } ENDFORALL ;
    sort(sortlist.begin(),sortlist.end(),fieldSort) ;
    i = 0 ;
    Map convert ;
    convert.allocate(faces) ;
    FORALL(faces,fc) {
      convert[fc] = sortlist[i++].second ;
      count_reorder[fc] = (face2node.end(convert[fc])-
                           face2node.begin(convert[fc])) ;
    } ENDFORALL ;
    Map clt,crt ;
    clt.allocate(faces) ;
    crt.allocate(faces) ;
    FORALL(faces,fc) {
      clt[fc] = cl[convert[fc]] ;
      crt[fc] = cr[convert[fc]] ;
    } ENDFORALL ;
    cl.setRep(clt.Rep()) ;
    cr.setRep(crt.Rep()) ;
    multiMap face2nodet ;
    face2nodet.allocate(count_reorder) ;
    FORALL(faces,fc) {
      int sz = count_reorder[fc] ;
      for(int j=0;j<sz;++j)
        face2nodet[fc][j] = face2node[convert[fc]][j] ;
    } ENDFORALL ;
    face2node.setRep(face2nodet.Rep()) ;

    // update remap from global to file numbering for faces after sorting
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    dMap g2f ;
    g2f = df->g2f.Rep() ;
    vector<pair<int, int> > remap_update(faces.size()) ;
    int cnt=0 ;
    FORALL(faces,fc) {
      remap_update[cnt].second = fc ;
      remap_update[cnt].first = g2f[convert[fc]] ;
      cnt++ ;
    } ENDFORALL ;
    facts.update_remap(remap_update) ;
    
    using std::cout ;
    using std::endl ;

    vector<int> saddr(MPI_processes)  ;
    for(int i=0;i<MPI_processes;++i) {
      saddr[i] = node_ptn[i].size() ;
    }
    vector<int> raddr(MPI_processes) ;
    MPI_Alltoall(&saddr[0],1,MPI_INT,
                 &raddr[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int b = *nodes.begin() ;
    int sum = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      int tmp = raddr[i] ;
      raddr[i] = b+sum ;
      sum += tmp ;
    }
    MPI_Alltoall(&raddr[0],1,MPI_INT,
                 &saddr[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    // Renumber maps (targets nodes and cells)
    dMap remap;
    for(int i=0;i<MPI_processes;++i) {
      int k = 0 ;
      FORALL(node_ptn[i], li) {
        remap[li] = saddr[i]+k ;
        k++ ;
      } ENDFORALL ;
    }

    entitySet orig_cells ;
    for(int i=0;i<MPI_processes;++i) {
      saddr[i] = cell_ptn[i].size() ;
      orig_cells += cell_ptn[i] ;
    }
    MPI_Alltoall(&saddr[0],1,MPI_INT,
                 &raddr[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    b = *cells.begin() ;
    sum = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      int tmp = raddr[i] ;
      raddr[i] = b+sum ;
      sum += tmp ;
    }
    MPI_Alltoall(&raddr[0],1,MPI_INT,
                 &saddr[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    for(int i=0;i<MPI_processes;++i) {
      int k = 0 ;
      FORALL(cell_ptn[i], li) {
        remap[li] = saddr[i]+k ;
        k++ ;
      } ENDFORALL ;
    }

    vector<entitySet> bcptn = all_collect_vectors(bcsurf_ptn[MPI_rank],MPI_COMM_WORLD) ;
    vector<entitySet> bcalloc = all_collect_vectors(bcsurfset,MPI_COMM_WORLD) ;
    for(int i=0;i<MPI_processes;++i) {
      if(bcptn[i].size() != bcalloc[i].size()) {
	cerr << "WARNING, boundary faces information inconsistent in remap"
	     << endl ;
      }
      entitySet::const_iterator fei = bcptn[i].begin() ;
      entitySet::const_iterator tei = bcalloc[i].begin() ;
      while(fei!=bcptn[i].end() && tei!=bcalloc[i].end()) {
	remap[*fei] = * tei ;
	fei++ ;
	tei++ ;
      }
    }
    // WARNING currently copying all boundary data  for all processors,
    // subsequent steps rely on this but it isn't really conforming to 
    // the way other data is handled.
    entitySet bcallalloc = bcalloc[0] ;
    for(int i=1;i<MPI_processes;++i) 
      bcallalloc += bcalloc[i] ;
    
    boundary_names.allocate(bcallalloc) ;
    boundary_tags.allocate(bcallalloc) ;
    entitySet bcdom = tmp_boundary_names.domain() ;
    FORALL(bcdom,ii) {
      boundary_names[remap[ii]]= tmp_boundary_names[ii] ;
      boundary_tags[remap[ii]] = tmp_boundary_tags[ii] ;
    } ENDFORALL ;

    orig_cells += bcdom ;
    entitySet out_of_dom ;
    MapRepP f2n = MapRepP(face2node.Rep()) ;
    out_of_dom += cr.image(cr.domain())-orig_cells ;
    out_of_dom += cl.image(cl.domain())-orig_cells ;
    out_of_dom += f2n->image(f2n->domain())-old_nodes ;
    entitySet old_dom = orig_cells+old_nodes ;
    vector<entitySet> old_ptn = all_collect_vectors(old_dom) ;
    {
      storeRepP PRep = remap.Rep() ;
      fill_clone(PRep,out_of_dom,old_ptn) ;
    }

    MapRepP(face2node.Rep())->compose(remap,faces) ;
    MapRepP(cr.Rep())->compose(remap,faces) ;
    MapRepP(cl.Rep())->compose(remap,faces) ;
    
  }

  //Note: This function is designed for serial version.
  //Input:
  // nodes, faces, cells
  // t_pos: position of nodes(dynamic version)
  // tmp_cl, tmp_cr: mapping from face to cell on left and right side(dynamic version)
  // tmp_face2node: mapping from face to nodes (dynamic version)
  //Output:
  // pos, cl, cr, face2node: static version of structures in the Input
  void copyGridStructures( entitySet nodes, entitySet faces, entitySet cells,
			   entitySet bcset,
			   const store<vector3d<double> > &t_pos,
			   const Map &tmp_cl, const Map &tmp_cr,
			   const multiMap &tmp_face2node,
			   const store<string> &tmp_boundary_names,
			   const store<string> &tmp_boundary_tags,
			   store<vector3d<double> > &pos, Map &cl, Map &cr,
			   multiMap &face2node,
			   store<string> &boundary_names,
			   store<string> &boundary_tags) {

    dMap identity_map;
    FORALL(nodes, ei) {
      identity_map[ei] = ei;
    } ENDFORALL ;
    FORALL(faces, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;
    FORALL(cells, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;
    FORALL(bcset, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;

    pos = t_pos.Rep()->remap(identity_map);
    cl = tmp_cl.Rep()->remap(identity_map);
    cr = tmp_cr.Rep()->remap(identity_map);
    face2node = MapRepP(tmp_face2node.Rep())->get_map();
    boundary_names = tmp_boundary_names.Rep()->remap(identity_map) ;
    boundary_tags = tmp_boundary_tags.Rep()->remap(identity_map) ;
  }

  void assignOwner(vector<pair<int,pair<int,int> > > &scratchPad,
                   vector<entitySet> ptn,
                   vector<entitySet> &out_ptn) ;
  
  void fill_clone_proc( map<int,int> &mapdata, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {

    memSpace("fill_clone start") ;
    vector<entitySet> recv_req(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      if(i!=MPI_rank) 
        recv_req[i] = out_of_dom & init_ptn[i] ;
    
    memSpace("reqcomp") ;
    // send the recieve requests
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    for(int i=0;i<MPI_processes;++i)
      send_count[i] = recv_req[i].num_intervals() * 2 ;
    MPI_Alltoall(send_count,1,MPI_INT, recv_count, 1, MPI_INT,MPI_COMM_WORLD) ;

    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }
    int mp = MPI_processes-1 ;
    int send_sizes = send_displacement[mp]+send_count[mp] ;
    int recv_sizes = recv_displacement[mp]+recv_count[mp] ; 
    memSpace("before allocating send_set") ;
    int * send_set_buf = new int[send_sizes] ;
    int * recv_set_buf = new int[recv_sizes] ;

    for(int i=0;i<MPI_processes;++i) {
      for(size_t j=0;j<recv_req[i].num_intervals();++j) {
        send_set_buf[send_displacement[i]+j*2  ] = recv_req[i][j].first ;
        send_set_buf[send_displacement[i]+j*2+1] = recv_req[i][j].second ;
      }
    }

    MPI_Alltoallv(send_set_buf, send_count, send_displacement , MPI_INT,
		  recv_set_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;

    vector<entitySet> send_set(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i) {
      for(int j=0;j<recv_count[i]/2;++j) {
        int i1 = recv_set_buf[recv_displacement[i]+j*2  ] ;
        int i2 = recv_set_buf[recv_displacement[i]+j*2+1] ;
        send_set[i] += interval(i1,i2) ;
      }
    }
    delete[] recv_set_buf ;
    delete[] send_set_buf ;
    
#ifdef DEBUG
    // Sanity check, no send set should be outside of entities we own
    entitySet mydom = init_ptn[MPI_rank] ;
    for(int i=0;i<MPI_processes;++i) {
      if((send_set[i] & mydom) != send_set[i]) {
        cerr << "problem with partitioning in fill_clone!" ;
        debugout << "send_set["<< i << "] = " << send_set[i]
                 << "not owned = " << (send_set[i]-mydom) << endl ;
      }
    }
#endif
    
    // Now that we know what we are sending and receiving (from recv_req and
    // send_set) we can communicate the actual information...

    // Compute sizes of sending buffers
    for(int i=0;i<MPI_processes;++i) 
      send_count[i] =  send_set[i].size() ;

    // Get sizes needed for receiving buffers
    MPI_Alltoall(send_count,1,MPI_INT, recv_count, 1, MPI_INT,MPI_COMM_WORLD) ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }

    send_sizes = send_displacement[mp]+send_count[mp] ;
    recv_sizes = recv_displacement[mp]+recv_count[mp] ;
    
    memSpace("before allocating send_store") ;
    int *send_store = new int[send_sizes] ;
    int *recv_store = new int[recv_sizes] ;

    for(int i=0;i<send_sizes;++i)
      send_store[i] = 0 ;
    

    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;

      FORALL(send_set[i],ii) {
        send_store[send_displacement[i]+loc_pack] = mapdata[ii] ;
        loc_pack++ ;
      } ENDFORALL ;
    }
    
    MPI_Alltoallv(send_store, send_count, send_displacement, MPI_INT,
		  recv_store, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;
    
    memSpace("before mapdata insert") ;
    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;
      FORALL(recv_req[i],ii) {
        mapdata[ii] = recv_store[recv_displacement[i]+loc_pack] ;
        loc_pack++ ;
      } ENDFORALL ;
    }
    delete[] recv_store ;
    delete[] send_store ;

    delete[] recv_count ;
    delete[] send_count ;
    delete[] send_displacement ;
    delete[] recv_displacement ;
  }
  

  vector<entitySet> partitionFaces(vector<entitySet> cell_ptn, const Map &cl,
                                   const Map &cr,
				   const store<string> &boundary_tags) {
    map<int,int> P ;
    entitySet cells ;
    for(int i=0;i<MPI_processes;++i) {
      FORALL(cell_ptn[i],cc) {
        P[cc] = i ;
      } ENDFORALL ;
      cells+= cell_ptn[i] ;
    }
    vector<entitySet> ptn_cells = all_collect_vectors(cells) ;
    entitySet faces = cl.domain() & cr.domain() ;
    entitySet dom = cl.image(faces) | cr.image(faces) ;


    dom -= cells ;
    
    fill_clone_proc(P,dom,ptn_cells) ;

    entitySet refset = boundary_tags.domain() ;
    entitySet bcfaces = cr.preimage(refset).first ;
    vector<entitySet> face_ptn(MPI_processes) ;
    FORALL(bcfaces,fc) {
      face_ptn[P[cl[fc]]] += fc ;
    } ENDFORALL ;
    faces -= bcfaces ;

    entitySet boundary_faces ; // Boundary between processors
    FORALL(faces,fc) {
      if(P[cl[fc]] == P[cr[fc]])
        face_ptn[P[cl[fc]]] += fc ;
      else
        boundary_faces += fc ;
    } ENDFORALL ;
    memSpace("start face_ptn") ;
    vector<int> curr_sizes(MPI_processes),tot_sizes(MPI_processes) ;


    // Number of balancing steps.  In the balancing step, the faces that
    // share proceessors are allocated by selected processors, then all
    // processors share current face counts.
    int STEPS = min(MPI_processes,13);
    for(int s=0;s<STEPS;++s) {
      
      memSpace("STEPS") ;
      for(int i=0;i<MPI_processes;++i)
        curr_sizes[i] = face_ptn[i].size() ;

      MPI_Allreduce(&curr_sizes[0],&tot_sizes[0],MPI_processes,MPI_INT,MPI_SUM,
                    MPI_COMM_WORLD) ;

      if(MPI_rank%STEPS == s) { // My processors turn to assign faces
        FORALL(boundary_faces,fc) {
          int Pl = P[cl[fc]] ;
          int Pr = P[cr[fc]] ;
          if(tot_sizes[Pl] < tot_sizes[Pr]) {
            tot_sizes[Pl]+=1 ;
            face_ptn[Pl] += fc ;
          } else {
            tot_sizes[Pr]+=1 ;
            face_ptn[Pr] += fc ;
          }
        } ENDFORALL ;
      }
    }

    for(int i=0;i<MPI_processes;++i)
      curr_sizes[i] = face_ptn[i].size() ;

    MPI_Allreduce(&curr_sizes[0],&tot_sizes[0],MPI_processes,MPI_INT,MPI_SUM,
                  MPI_COMM_WORLD) ;

    if(MPI_rank ==0) {
      debugout << "balanced face sizes:" ;
      for(int i=0;i<MPI_processes;++i)
        debugout << ' ' << tot_sizes[i] ;
      debugout << endl ;
    }
    return face_ptn ;
  }


  vector<entitySet> partitionNodes(vector<entitySet> face_ptn, MapRepP face2node,entitySet old_node_dom) {

    // find node_ptn that best matches the face partition.  Loop over faces
    entitySet fdom ;
    for(int i=0;i<MPI_processes;++i)
      fdom += face_ptn[i] ;
    store<int> procmap ;
    procmap.allocate(fdom) ;
    for(int i=0;i<MPI_processes;++i) {
      FORALL(face_ptn[i],fc) {
        procmap[fc] = i ;
      } ENDFORALL ;
    }
    
    multiMap f2n ;
    f2n = storeRepP(face2node) ;

    vector<pair<int,pair<int,int> > > scratchPad ;

    FORALL(fdom,fc) {
      pair<int,int> p2i(procmap[fc],1) ;
      int sz = f2n[fc].size() ;
      for(int ii=0;ii<sz;++ii) {
        scratchPad.push_back(pair<int,pair<int,int> >(f2n[fc][ii],p2i)) ;
      }
    } ENDFORALL ;
    
    vector<entitySet> node_ptn_old = all_collect_vectors(old_node_dom) ;
    vector<entitySet> node_ptn(MPI_processes) ;
    assignOwner(scratchPad,node_ptn_old,node_ptn) ;
    return node_ptn ;
  }

  vector<entitySet> transposePtn(const vector<entitySet> &ptn) {
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = ptn[i].num_intervals()*2 ;
    vector<int> recv_sz(MPI_processes) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }
    //    outRep->allocate(new_alloc) ;
    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  MPI_processes; ++i)
      for(size_t j=0;j<ptn[i].num_intervals();++j) {
        send_store[send_displacement[i]+j*2] = ptn[i][j].first ;
        send_store[send_displacement[i]+j*2+1] = ptn[i][j].second ;
      }


    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;

    vector<entitySet> ptn_t(MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        ptn_t[i] += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return ptn_t ;
  }

  inline bool sortFirstORB(const pair<int,pair<int,int> > &p1,
                           const pair<int,pair<int,int> > &p2) {
    return p1.first < p2.first ;
  }

  void assignOwner(vector<pair<int,pair<int,int> > > &scratchPad,
                   vector<entitySet> ptn,
                   vector<entitySet> &out_ptn) {

    std::sort(scratchPad.begin(),scratchPad.end()) ;
    vector<pair<int,pair<int,int> > >::iterator i1,i2 ;

    // Note, we are assuming scratch pad is non-zero size
    if(scratchPad.size() > 0) {
      i1 = scratchPad.begin() ;
      i2 = i1+1 ;
      while(i2!=scratchPad.end()) {
        while(i2!=scratchPad.end() &&
              i1->first==i2->first && i1->second.first==i2->second.first) {
          i1->second.second += i2->second.second ;
          i2++ ;
        }
        i1++ ;
        if(i2!=scratchPad.end()) {
          *i1 = *i2 ;
          i2++ ;
        }
      }
      scratchPad.erase(i1,scratchPad.end()) ;
    }

    vector<pair<int,pair<int,int> > > nsplits(MPI_processes-1) ;
    for(int i=1;i<MPI_processes;++i) {
      nsplits[i-1].first = ptn[i].Min() ;
      nsplits[i-1].second.first = 0 ;
      nsplits[i-1].second.second = 0 ;
    }

    parSplitSort(scratchPad,nsplits,sortFirstORB,MPI_COMM_WORLD) ;

    for(size_t x=0;x!=scratchPad.size();) {
      size_t y = x+1 ;
      while (y < scratchPad.size() &&
             scratchPad[x].first == scratchPad[y].first)
        ++y ;

      int imax = 0 ;
      int pmax = -1 ;
      for(size_t j=x;j<y;++j) {
        // Sum total number of references from this processor
        int tot = scratchPad[j].second.second ;
        for(size_t k=j+1;k<y;++k)
          if(scratchPad[j].second.first == scratchPad[k].second.first)
            tot += scratchPad[k].second.second ;
        if(tot > imax) {
          imax = tot ;
          pmax = scratchPad[j].second.first ;
        }
      }
      int nd = scratchPad[x].first ;
      out_ptn[pmax] += nd ;
      x = y ;
    }
    // Check to see if any sets are left out.
    entitySet assigned ;
    for(int i=0;i<MPI_processes;++i)
      assigned += out_ptn[i] ;
    entitySet unassigned = ptn[MPI_rank] - assigned ;
    out_ptn[MPI_rank] += unassigned ; // allocate unassigned entities to
    // original processor

  }

    
  void ORB_Partition_Mesh(const vector<entitySet> &local_nodes,
                          const vector<entitySet> &local_faces,
                          const vector<entitySet> &local_cells,
                          const store<vector3d<double> > &pos,
                          const Map &cl, const Map &cr,
                          const multiMap &face2node,
			  const store<string> &boundary_tags,
                          vector<entitySet> &cell_ptn,
                          vector<entitySet> &face_ptn,
                          vector<entitySet> &node_ptn) {

    vector<entitySet> tmp(MPI_processes) ; // Initialize partition vectors
    cell_ptn = tmp ;
    face_ptn = tmp ;
    node_ptn = tmp ;

    // Compute face center
    dstore<vector3d<float> > tmp_pos ;
    FORALL(pos.domain(),pi) {
      tmp_pos[pi] = vector3d<float>(realToFloat(pos[pi].x),
				    realToFloat(pos[pi].y),
				    realToFloat(pos[pi].z)) ;
    } ENDFORALL ;
    entitySet fdom = face2node.domain() ;
    entitySet total_dom =
      Loci::MapRepP(face2node.Rep())->image(fdom) + pos.domain() ;
    Loci::storeRepP sp = tmp_pos.Rep() ;
    vector<entitySet> ptn = local_nodes ;
    fill_clone(sp,total_dom,ptn) ;
    vector<vector3d<float> > fcenter(fdom.size()) ;
    int i=0 ;
    FORALL(fdom,fc) {
      int sz = face2node[fc].size() ;
      vector3d<float> pnt(0.,0.,0.) ;
      for(int ii=0;ii<sz;++ii)
        pnt += tmp_pos[face2node[fc][ii]] ;
      pnt *= 1./float(sz) ;
      fcenter[i++] = pnt ;
    } ENDFORALL ;

    // perform ORB partition of faces
    vector<int> fprocmap ;
    ORBPartition(fcenter,fprocmap,MPI_COMM_WORLD) ;
    i=0 ;
    // Create face_ptn ;
    FORALL(fdom,fc) {
      face_ptn[fprocmap[i++]] += fc ;
    } ENDFORALL ;

    // Now find node_ptn and cell_ptn that best matches the face partition
    vector<pair<int,pair<int,int> > > scratchPad ;
    i=0 ;
    FORALL(fdom,fc) {
      pair<int,int> p2i(fprocmap[i++],1) ;
      int sz = face2node[fc].size() ;
      for(int ii=0;ii<sz;++ii)
        scratchPad.push_back(pair<int,pair<int,int> >(face2node[fc][ii],p2i)) ;
    } ENDFORALL ;

    assignOwner(scratchPad,local_nodes,node_ptn) ;


    entitySet refset = boundary_tags.domain() ;
    scratchPad.clear() ;
    i=0 ;
    FORALL(fdom,fc) {
      pair<int,int> p2i(fprocmap[i++],1) ;
      scratchPad.push_back(pair<int,pair<int,int> >(cl[fc],p2i)) ;
      if(!refset.inSet(cr[fc]))
        scratchPad.push_back(pair<int,pair<int,int> >(cr[fc],p2i)) ;
    } ENDFORALL ;

    assignOwner(scratchPad,local_cells,cell_ptn) ;
  }


  //Description: Reads grid structures in the fact database
  //Input: facts and grid file name
  //Output: true if sucess
  bool readFVMGrid(fact_db &facts, string filename) {
    double t1 = MPI_Wtime() ;
    //	timer_token read_file_timer = new timer_token;
    //	if(collect_perf_data)
    //		read_file_timer = perfAnalysis->start_timer("Reading in FVM Grid");
    memSpace("readFVMGrid Start") ;
    vector<entitySet> local_nodes;
    vector<entitySet> local_cells;
    vector<entitySet> local_faces;

    store<vector3d<double> > t_pos;
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;
    store<string> tmp_boundary_names ;
    store<string> tmp_boundary_tags ;

    int max_alloc = facts.get_max_alloc() ;

    bool useVOG = true ;

    string input_file = filename ;
    string::size_type spos = string::npos ;
    if((spos = input_file.rfind('.')) != string::npos) {
      input_file.erase(0,spos) ;
      if(input_file == ".xdr")
        useVOG = false ;
    }

    vector<pair<string,entitySet> > volTags ;

    if(useVOG) {
      if(!readGridVOG(local_nodes, local_faces, local_cells,
                      t_pos, tmp_cl, tmp_cr, tmp_face2node,
		      tmp_boundary_names,tmp_boundary_tags, volTags,
                      max_alloc, filename))
        return false;
    } else {
      cerr << "XDR Grid File Reading Support has been removed!!!, convert file to VOG format." << endl ;
      Loci::Abort() ;
      return false;
    }
    
    memSpace("after reading grid") ;

    // Identify boundary tags

    entitySet global_boundary_cells = tmp_boundary_tags.domain() ;
    memSpace("after all_collect boundary cells") ;
    
    if(MPI_processes == 1) {

      int npnts = local_nodes[0].size();
      int nfaces = local_faces[0].size();
      int ncells = local_cells[0].size();
      int nbcs = global_boundary_cells.size() ;
      entitySet nodes = facts.get_distributed_alloc(npnts).first ;
      entitySet faces = facts.get_distributed_alloc(nfaces).first ;
      entitySet cells = facts.get_distributed_alloc(ncells).first;
      entitySet bcset = facts.get_distributed_alloc(nbcs).first ;


      store<vector3d<double> > pos ;
      Map cl ;
      Map cr ;
      multiMap face2node ;
      store<string> boundary_names ;
      store<string> boundary_tags ;

      copyGridStructures(nodes, faces, cells, bcset,
                         t_pos, tmp_cl, tmp_cr, tmp_face2node,
			 tmp_boundary_names,tmp_boundary_tags,
                         pos, cl, cr, face2node,
			 boundary_names,boundary_tags);

      facts.create_fact("cl", cl) ;
      facts.create_fact("cr", cr) ;
      facts.create_fact("pos", pos) ;
      facts.create_fact("face2node",face2node) ;
      facts.create_fact("boundary_names", boundary_names) ;
      facts.create_fact("boundary_tags", boundary_tags) ;

      for(size_t i=0;i<volTags.size();++i) {
        param<string> Tag ;
        *Tag = volTags[i].first ;
        Tag.set_entitySet(volTags[i].second) ;
        ostringstream oss ;
        oss << "volumeTag(" << volTags[i].first << ")" ;
        facts.create_fact(oss.str(),Tag) ;
      }
        
      return true ;
    }

    memSpace("before partitioning") ;

    vector<entitySet> cell_ptn,face_ptn,node_ptn ;

    if(use_orb_partition) {
      ORB_Partition_Mesh(local_nodes, local_faces, local_cells,
                         t_pos, tmp_cl, tmp_cr, tmp_face2node,
			 tmp_boundary_tags,
                         cell_ptn,face_ptn,node_ptn) ;
    } else {
      // Partition Cells
      bool useMetis = !use_simple_partition ;
      
      // If the number of cells per processor is too thin, metis
      // no longer is an effective partitioner, so switch to the 
      // simple partitioner instead
      if(local_cells[0].size() < 450)  {
	if(useMetis && MPI_rank==0) 
	  debugout << "switching to simple partition due to small number of cells per proc!"<< endl ;
	useMetis = false ;
      }
      if(useMetis) {
        cell_ptn = newMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr,
					    tmp_boundary_tags) ;
      } else {
        cell_ptn = vector<entitySet>(MPI_processes) ;
        cell_ptn[MPI_rank] = local_cells[MPI_rank] ;
	// read in additional vertex weights if any
	if(load_cell_weights) {
	  // check if the file exists
	  int file_exists = 1 ;
	  if(Loci::MPI_rank == 0) {
	    struct stat buf ;
	    if(stat(cell_weight_file.c_str(),&buf) == -1 ||
	       !S_ISREG(buf.st_mode))
	      file_exists = 0 ;
	  }
	  MPI_Bcast(&file_exists,1,MPI_INT,0,MPI_COMM_WORLD) ;
	  
	  if(file_exists == 1) {
	    if(Loci::MPI_rank == 0) {
	      std::cout << "simple partition reading additional cell weights from: "
			<< cell_weight_file << std::endl ;
	    }
	    
	    // create a hdf5 handle
	    hid_t file_id = Loci::hdf5OpenFile(cell_weight_file.c_str(),
					       H5F_ACC_RDONLY, H5P_DEFAULT) ;
	    if(file_id < 0) {
	      std::cerr << "...file reading failed..., Aborting" << std::endl ;
	      Loci::Abort() ;
	    }
        
	    // read
	    store<int> cell_weights ;
	    
	    readContainerRAW(file_id,"cell weight", cell_weights.Rep(),
			     MPI_COMM_WORLD) ;
	    
	    if(cell_weights.domain() != local_cells[Loci::MPI_rank]) {
	      cerr << "cell weights partition inconsistent!" << endl ;
	      Loci::Abort() ;
	    }
        
	    Loci::hdf5CloseFile(file_id) ;
	    int tot_weight = 0 ;

	    cell_ptn[MPI_rank] = EMPTY ;
	    entitySet dom = cell_weights.domain() ;
	    FORALL(dom,i) {
	      tot_weight += cell_weights[i] ;
	    } ENDFORALL ;
	    vector<int> pweights(MPI_processes,0) ;
	    MPI_Allgather(&tot_weight,1,MPI_INT,&pweights[0],1,MPI_INT,
			  MPI_COMM_WORLD) ;
	    vector<int> woffsets(MPI_processes+1,0) ;
	    for(int i=0;i<MPI_processes;++i)
	      woffsets[i+1] = pweights[i]+woffsets[i] ;

	    // compute the weight per processor
	    int wpp = ((woffsets[MPI_processes]+MPI_processes-1)/MPI_processes) ;
	    // Now compute local weighted sums
	    int ncel = cell_weights.domain().size() ;
	    vector<int> wts(ncel+1,woffsets[MPI_rank]) ;
	    int cnt = 0 ;
	    FORALL(dom,i) {
	      wts[cnt+1] = wts[cnt] + cell_weights[i] ;
	      cnt++ ;
	    } ENDFORALL ;
	    
	    entitySet pdom = local_cells[MPI_rank] ;
	    cnt = 0 ;
	    FORALL(pdom,i) {
	      cell_ptn[wts[cnt]/wpp] += i ;
	      cnt++ ;
	    } ENDFORALL ;
	  }
	  
	}
      }

      memSpace("mid partitioning") ;
      face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr,tmp_boundary_tags) ;
      memSpace("after partitionFaces") ;

      node_ptn = partitionNodes(face_ptn,
                                MapRepP(tmp_face2node.Rep()),
                                t_pos.domain()) ;
    }
    
    vector<entitySet> bcsurf_ptn(MPI_processes) ;
    entitySet refset = tmp_boundary_tags.domain() ;
    
    // round robin allocate boundary faces
    // NOTE!!!! This could be improved.
    int cnt = 0 ;
    //    int refsetsz = refset.size() ;
    FORALL(refset,ii) {
      if(cnt == MPI_rank)
	bcsurf_ptn[cnt] += ii ;
      cnt++ ;
      if(cnt == MPI_processes)
	cnt = 0 ;
    } ENDFORALL ;

    memSpace("after partitioning") ;
      
    vector<entitySet> cell_ptn_t = transposePtn(cell_ptn) ;
    vector<entitySet> face_ptn_t = transposePtn(face_ptn) ;
    vector<entitySet> node_ptn_t = transposePtn(node_ptn) ;

    int newnodes = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newnodes += node_ptn_t[p].size() ;

    vector<int> node_alloc(newnodes) ;
    int i=0;
    for(int p=0;p<MPI_processes;++p)
      FORALL(node_ptn_t[p], ni) {
        node_alloc[i++] = ni ;
      } ENDFORALL;

    entitySet nodes = facts.get_distributed_alloc(node_alloc).first ;
    node_alloc.resize(0) ;

    int newfaces = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newfaces += face_ptn_t[p].size() ;

    vector<int> face_alloc(newfaces) ;
    i = 0 ;
    for(int p=0;p<MPI_processes;++p)
      FORALL(face_ptn_t[p], ni) {
        face_alloc[i++] = ni ;
      }ENDFORALL;

    entitySet faces = facts.get_distributed_alloc(face_alloc).first ;
    face_alloc.resize(0) ;

    int newcells = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newcells += cell_ptn_t[p].size() ;

    vector<int> cell_alloc(newcells) ;
    i = 0 ;
    for(int p=0;p<MPI_processes;++p)
      FORALL(cell_ptn_t[p], ni) {
        cell_alloc[i++] = ni ;
      }ENDFORALL;

    entitySet cells = facts.get_distributed_alloc(cell_alloc).first ;

    Loci::debugout << "nodes = " << nodes << ", size= "
                   << nodes.size() << endl;
    Loci::debugout << "faces = " << faces << ", size = "
                   << faces.size() << endl ;
    Loci::debugout << "cells = " << cells << ", size = "
                   << cells.size() << endl ;

    vector<int> bcsurf_alloc(bcsurf_ptn[MPI_rank].size()) ;
    i=0 ;
    FORALL(bcsurf_ptn[MPI_rank],ii) {
      bcsurf_alloc[i++] = ii ;
    } ENDFORALL ;
    
    entitySet bcsurfset = facts.get_distributed_alloc(bcsurf_alloc).first ;
    
    memSpace("before remapGridStructures") ;
    Map cl, cr ;
    multiMap face2node ;
    store<vector3d<double> > pos ;
    store<string> boundary_names,boundary_tags ;
    remapGrid(node_ptn, face_ptn, cell_ptn,
              node_ptn_t, face_ptn_t, cell_ptn_t,
              t_pos, tmp_cl, tmp_cr, tmp_face2node,
	      bcsurf_ptn,tmp_boundary_names,tmp_boundary_tags,
              nodes, faces, cells,
              pos, cl, cr, face2node,
	      boundary_names,boundary_tags, bcsurfset, facts);
    memSpace("after remapGridStructures") ;

    facts.create_fact("cl", cl) ;
    facts.create_fact("cr", cr) ;
    facts.create_fact("pos", pos) ;
    facts.create_fact("face2node",face2node) ;
    facts.create_fact("boundary_names", boundary_names) ;
    facts.create_fact("boundary_tags", boundary_tags) ;

    // update remap from global to file numbering for faces after sorting
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    dMap g2f ;
    g2f = df->g2f.Rep() ;

    for(size_t i=0;i<volTags.size();++i) {
      param<string> Tag ;
      *Tag = volTags[i].first ;

      // Map entitySet to new ordering
      entitySet tagset ;
      FORALL(cells,cc) {
        if(volTags[i].second.inSet(g2f[cc])) 
          tagset += cc ;
      } ENDFORALL ;
      Tag.set_entitySet(tagset) ;
      ostringstream oss ;
      oss << "volumeTag(" << volTags[i].first << ")" ;
      facts.create_fact(oss.str(),Tag) ;
    }
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(read_file_timer);
    double t2 = MPI_Wtime() ;
    debugout << "Time to read in file '" << filename << ", is " << t2-t1
             << endl ;
    memSpace("returning from FVM grid reader") ;
    return true ;
  }

  enum matrix_coloring_type {COLOR_DEFAULT, COLOR_DFS} ;

  void create_ref(fact_db &facts) {
    store<string> boundary_names ;
    store<string> boundary_tags ;
    boundary_names = facts.get_fact("boundary_names") ;
    boundary_tags = facts.get_fact("boundary_tags") ;
    Map cr ;
    cr = facts.get_fact("cr") ;
    entitySet bdom = boundary_names.domain() ;
    
    constraint boundary_faces ;
    boundary_faces = facts.get_fact("boundary_faces") ;
    entitySet refdom = *boundary_faces ;
    
    Map ref ;
    ref.allocate(refdom) ;
    FORALL(refdom,i1) {
      ref[i1] = cr[i1] ;
    } ENDFORALL ;
    facts.create_fact("ref",ref) ;
  }

  void create_ghost_cells(fact_db &facts) {
    constraint interior_faces,boundary_faces ;
    constraint geom_cells, ghost_cells, cells ;
    Map cl,cr ;
    std::vector<int> vec ;
    std::vector<int>::const_iterator vi ;
    interior_faces = facts.get_variable("interior_faces") ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    FORALL(*interior_faces,fc) {
      vec.push_back(cl[fc]) ;
      vec.push_back(cr[fc]) ;
    } ENDFORALL ;
    std::sort(vec.begin(), vec.end()) ;
    for(vi = vec.begin(); vi != vec.end(); ++vi)
      *geom_cells += *vi ;

    FORALL(*boundary_faces,fc) {
      *geom_cells += cl[fc] ;
    } ENDFORALL ;

    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    entitySet global_geom = all_collect_entitySet(*geom_cells,facts) ;
    *geom_cells = global_geom & init_ptn[ MPI_rank] ;
    *boundary_faces &= init_ptn[ MPI_rank] ;
    std::pair<entitySet, entitySet> ghost_pair = facts.get_distributed_alloc((*boundary_faces).size()) ;
    entitySet tmp_ghost = ghost_pair.first ;
    entitySet::const_iterator ei = tmp_ghost.begin() ;
    FORALL(*boundary_faces,fc) {
      cr[fc] = *ei++;
    } ENDFORALL ;
    *ghost_cells = ghost_pair.first ;

    facts.update_fact("cr",cr) ;

    *cells = *geom_cells + *ghost_cells ;

    facts.create_fact("geom_cells",geom_cells) ;
    facts.create_fact("ghost_cells",ghost_cells) ;
    facts.create_fact("cells",cells) ;

    Loci::debugout << "geom_cells = " << *geom_cells << endl ;
    Loci::debugout << "ghost_cells = " << *ghost_cells << endl ;
    Loci::debugout << "cells = " << *cells << endl ;
  }

  void create_face_info(fact_db &facts) {
    Map cl, cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    constraint faces ;
    faces = (cl.domain() & cr.domain()) ;
    
    //    faces = all_collect_entitySet(*faces) ;

    facts.create_fact("faces",faces) ;
    store<string> boundary_names ;
    boundary_names = facts.get_variable("boundary_names") ;
    entitySet bcset = boundary_names.domain() ;
    entitySet bcfaces = cr.preimage(bcset).first ;
    constraint boundary_faces ;
    boundary_faces = bcfaces ;
    //    boundary_faces = all_collect_entitySet(bcfaces) ;
    constraint interior_faces ;
    interior_faces = (*faces-*boundary_faces) ;
    facts.create_fact("boundary_faces",boundary_faces) ;
    facts.create_fact("interior_faces",interior_faces) ;
  }

  void color_matrix(fact_db &facts, matrix_coloring_type mct) {
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    multiMap c2c ;
    store<int> sizes ;
    constraint geom_cells ;
    geom_cells = facts.get_variable("geom_cells") ;
    sizes.allocate(*geom_cells) ;
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    constraint faces ;
    constraint interior_faces ;
    faces = facts.get_variable("faces") ;
    interior_faces = facts.get_variable("interior_faces") ;

    FORALL(*geom_cells, cc) {
      sizes[cc] = 0 ;
    } ENDFORALL ;
    FORALL(*interior_faces,fc) {
      if((init_ptn[Loci::MPI_rank].inSet(cl[fc]))&& (init_ptn[Loci::MPI_rank].inSet(cr[fc])))
	sizes[cl[fc]]++ ;
      if((init_ptn[Loci::MPI_rank].inSet(cr[fc])) && (init_ptn[Loci::MPI_rank].inSet(cl[fc])))
	sizes[cr[fc]]++ ;
    } ENDFORALL ;
    c2c.allocate(sizes) ;
    FORALL(*interior_faces,fc) {
      //FATAL(sizes[cl[fc]] < 1) ;
      //FATAL(sizes[cr[fc]] < 1) ;
      if((init_ptn[Loci::MPI_rank].inSet(cl[fc]))  && (init_ptn[Loci::MPI_rank].inSet(cr[fc])))
	c2c[cl[fc]][--sizes[cl[fc]]] = cr[fc] ;
      if((init_ptn[Loci::MPI_rank].inSet(cr[fc])) && (init_ptn[Loci::MPI_rank].inSet(cl[fc])))
	c2c[cr[fc]][--sizes[cr[fc]]] = cl[fc] ;
    } ENDFORALL ;
    dstore<int> color ;
    entitySet out_of_dom ;

    entitySet glob_geom = Loci::all_collect_entitySet(*geom_cells) ;
    FORALL(*interior_faces,fc) {
      if(!init_ptn[Loci::MPI_rank].inSet(cl[fc]))
	out_of_dom += cl[fc] ;
      if(!init_ptn[Loci::MPI_rank].inSet(cr[fc]))
	out_of_dom += cr[fc] ;
    } ENDFORALL ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      init_ptn[i] &= glob_geom ;
      entitySet tmp_set =  (init_ptn[i] & out_of_dom) ;
      for(entitySet::const_iterator ei = tmp_set.begin(); ei != tmp_set.end(); ++ei)
	color[*ei] = init_ptn[i].Min() ;
    }
    FORALL(*geom_cells,cc) {
      color[cc] = -1 ;
    } ENDFORALL ;
    int col = (*geom_cells).Min() ;
    vector<int> visited ;
    entitySet left_out = *geom_cells ;
    while(left_out != EMPTY) {
      vector<int> work ;
      work.push_back(left_out.Min()) ;
      while(work.size() != 0) {
	vector<int> working ;
	for(size_t i=0;i<work.size();++i) {
          int cc = work[i] ;
	  if(color[cc] == -1) {
	    color[cc] = col++ ;
	    visited.push_back(cc) ;
	    for(const int *pi = c2c.begin(cc);pi!=c2c.end(cc);++pi)
	      working.push_back(*pi) ;
	  }
	}
        work.swap(working) ;
      }
      entitySet visitSet = create_intervalSet(visited.begin(),visited.end()) ;
      left_out -= visitSet ;
    }

    FORALL(*interior_faces,fc) {
      if(color[cl[fc]] > color[cr[fc]])
        std::swap(cl[fc],cr[fc]) ;
    } ENDFORALL ;

  }
  void make_faces_consistent(fact_db &facts) {
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;
    store<vector3d<double> > fpos ;
    store<vector3d<double> > area ;

    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    constraint faces ;
    constraint interior_faces ;
    faces = facts.get_variable("faces") ;
    interior_faces = facts.get_variable("interior_faces") ;
    entitySet total_dom = Loci::MapRepP(face2node.Rep())->image(*faces) + pos.domain() ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    Loci::storeRepP pos_sp = pos.Rep() ;
    dstore<vector3d<double> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    Loci::storeRepP sp = tmp_pos.Rep() ;
    if(MPI_processes > 1)
      fill_clone(sp, total_dom, init_ptn) ;
    entitySet face_dom = face2node.domain() ;
    fpos.allocate(face_dom) ;
    area.allocate(face_dom) ;

    FORALL(face_dom,fc) {
      int nnodes = face2node.end(fc) - face2node.begin(fc) ;
      vector3d<double> fp(0,0,0) ;
      double w = 0 ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;

        double len = norm(p1-p2) ;

        fp += len*(p1+p2) ;
        w += len ;
      }
      fpos[fc] = fp/(2.*w) ;
      vector3d<double> a(0,0,0) ;
      for(int i=0;i<nnodes;++i) {
        vector3d<double> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<double> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;
        a += cross(p1-fpos[fc],p2-fpos[fc]) ;
      }
      area[fc] = .5*a ;
    } ENDFORALL ;
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    dstore<vector3d<double> > cpos ;
    dstore<double> cnum ;
    constraint geom_cells ;
    geom_cells = facts.get_variable("geom_cells") ;
    entitySet tmp_cells =  cl.image(*faces) | cr.image(*interior_faces) ;
    // Add cells owned by this processor!
    tmp_cells += (*geom_cells)& init_ptn[MPI_rank] ;
    cpos.allocate(tmp_cells) ;
    cnum.allocate(tmp_cells) ;
    FORALL(tmp_cells,cc) {
      cpos[cc] = vector3d<double>(0,0,0) ;
      cnum[cc] = 0 ;
    } ENDFORALL ;
    FORALL(*faces,fc) {
      double A = norm(area[fc]) ;
      cpos[cl[fc]] += A*fpos[fc] ;
      cnum[cl[fc]] += A ;
    } ENDFORALL ;
    FORALL(*interior_faces,fc) {
      double A = norm(area[fc]) ;
      cpos[cr[fc]] += A*fpos[fc] ;
      cnum[cr[fc]] += A ;
    } ENDFORALL ;
    Loci::storeRepP cp_sp = cpos.Rep() ;
    Loci::storeRepP cn_sp = cnum.Rep() ;
    entitySet clone_cells = tmp_cells - *geom_cells ;
    std::vector<Loci::storeRepP> v_cpos = send_global_clone_non(cp_sp, clone_cells, init_ptn) ;
    std::vector<Loci::storeRepP> v_cnum = send_global_clone_non(cn_sp, clone_cells, init_ptn) ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      entitySet dom = v_cpos[i]->domain() & cpos.domain() ;
      dstore<vector3d<double> > tmp_cpos(v_cpos[i]) ;
      dstore<double> tmp_cnum(v_cnum[i]) ;
      FORALL(dom, di) {
	cpos[di] += tmp_cpos[di] ;
	cnum[di] += tmp_cnum[di] ;
      } ENDFORALL ;
    }
    fill_clone(cp_sp, clone_cells, init_ptn) ;
    fill_clone(cn_sp, clone_cells, init_ptn) ;
    FORALL(tmp_cells,cc) {
      cpos[cc] = cpos[cc]/cnum[cc] ;
    } ENDFORALL ;

    vector<int> broken_faces ;

    FORALL(*interior_faces,fc) {
      vector3d<double> dv = cpos[cr[fc]]-cpos[cl[fc]] ;
      vector3d<double> dv2 = fpos[fc]-cpos[cl[fc]] ;
      vector3d<double> dv3 = cpos[cr[fc]]-fpos[fc] ;

      int t1 = (dot(area[fc],dv) <0.0)?1:0 ;
      int t2 = (dot(area[fc],dv2) <0.0)?1:0 ;
      int t3 = (dot(area[fc],dv3) <0.0)?1:0 ;
      int test = t1+t2+t3 ;
      if(test != 3 && test != 0) {
        debugout << "problem with face located at " << fpos[fc]
                 << endl ;
        broken_faces.push_back(fc) ;
      }


      else if(t1 == 1) { // Face oriented incorrectly
	int i = 0 ;
	int j = face2node.end(fc) - face2node.begin(fc) -1 ;
	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
	  i++ ;
	  j-- ;
	}
      }
    } ENDFORALL ;

    entitySet boundary_faces = *faces - *interior_faces ;
    FORALL(boundary_faces,fc) {
      const vector3d<double> center = fpos[fc] ;

      vector3d<double> ccenter ;
      ccenter = cpos[cl[fc]] ;
      vector3d<double> dv = center-ccenter ;
      if(dot(area[fc],dv) < 0.0) {
	int i = 0 ;
	int j = face2node.end(fc) - face2node.begin(fc) -1 ;
	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
	  i++ ;
	  j-- ;
	}
      }

    } ENDFORALL ;

    int rsize = 0 ;
    int size = broken_faces.size() ;

    MPI_Allreduce(&size,&rsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(rsize!=0) {
      if(MPI_rank == 0) {
        cerr << "Bad Grid: Non-Convex Cell (centroid outside cell bounds)" << endl ;
      }
      Loci::Abort() ;
    }
  }



  bool setupFVMGrid(fact_db &facts, string filename) {
    if(!readFVMGrid(facts,filename))
      return false ;

    memSpace("before create_face_info") ;
    create_face_info(facts) ;

    create_ref(facts) ;
    create_ghost_cells(facts) ;

    return true ;
  }


}

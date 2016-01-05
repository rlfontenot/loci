
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
#include <Loci_types.h>
#include "loci_globs.h"
#include <sstream>

#include <gconstraint.h>
#include <partition.h>
#include <gstore.h>
#include <gmap.h>
#include <gmultimap.h>
#include <gfact_db.h>
#include <gparameter.h>

#include <Tools/intervalSet.h>
#include <Tools/tools.h>
#include <Tools/simple_partition_long.h>
#include <hdf5_readwrite_long.h>
#include <distribute_long.h>
#include <map>
#include <fact_db.h>
#include <field_sort.h>

#include <LociGridReaders.h> //real_t is defined here

//for test, in copy_fact()
#include <constraint.h>
#include <store.h>
#include <multiMap.h>
#include <Map.h>
#include <parameter.h>

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
using Loci::gEntity;
using Loci::gEntitySet;

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

//typedef double real_t ;

namespace Loci {
#define MEMDIAG
  typedef int int_type; 
  typedef Loci::int_type gEntity ;
 
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
  //these variables can be put into the string to parse for partition
  //or defined as gParam
  // extern bool use_simple_partition ;
  //extern bool use_orb_partition ;
  extern bool load_cell_weights ;
  extern string cell_weight_file ;


  //Following assumption is needed for the next three functions
  //Assumption: We follow the convention that boundary cells are always on right side
  //of a face.
  //no modification now
  template<class T> void readVectorDist(hid_t group_id,
                                        const char *vector_name,
                                        vector<long> sizes,
                                        vector<T> &v) {

    hid_t dataset = 0 ;
    hid_t dspace = 0 ;
    memSpace("Start readVectorDist, sizes") ;
   
    if(MPI_rank == 0) {
      dataset = H5Dopen(group_id,vector_name) ;
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
      memSpace("close dataset, and start read for other processes") ;
      
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
  
  //inside a cluster, npnts and nfaces are still int
  //but total num_faces is gEntity
  gEntity getClusterNumFaces(unsigned char *cluster) {
    gEntity num_faces = 0 ;
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
  
 

  //no modification
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
  //no modification
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
  //no modification
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

  //no modification, only need by vogmerge.cc
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

  //no modification, only needed by vogmerge.cc
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



  
  //num_faces is gEntity
  //gStore replaced the traditional stores
  //to avoid remapping, read in nodeTable and cellTable first
  gEntity g_fillFaceInfo(unsigned char *cluster, gMultiMap &face2node,
                         gMap &cl, gMap &cr, gEntity face_base,
                         gEntity& max_cell, gEntitySet& boundary_faces, gEntitySet& boundary_tags ) {
    //copy the pointer to the front of the  cluster
    unsigned char *cluster_front = cluster;
    
    gEntity num_faces = 0 ;
    while(*cluster != 0) {
      int npnts = *cluster ;
      cluster++ ;
      int nfaces = *cluster ;
      cluster++ ;
      num_faces += nfaces ;
      cluster += nfaces * (npnts + 2) ;
    }
    cluster++ ;
   
    // read in tables
    long nodeTable[256] ;
    cluster = readTable(cluster,nodeTable) ;
    long cellTable[256] ;
    cluster = readTable(cluster,cellTable) ;

    unsigned char *cluster_back = cluster; 

    num_faces = 0 ;
    cluster = cluster_front;
    while(*cluster != 0) {
      int npnts = *cluster ;
      cluster++ ;
      int nfaces = *cluster ;
      cluster++ ;
      for(int i=0;i<nfaces;++i) {
        int fc = face_base+num_faces+i ;
        for(int j=0;j<npnts;++j){
          gEntity val = nodeTable[*cluster++];
          face2node.insert(gEntity(fc), val) ;
        }
        
        gEntity cl_val = cellTable[*cluster++];
        gEntity cr_val = cellTable[*cluster++];
        cl.insert(fc, cl_val) ;
        cr.insert(fc, cr_val) ;
        max_cell = max(max_cell, cl_val);
        max_cell = max(max_cell, cr_val);
        if(cr_val < 0){
          boundary_faces += fc;
          boundary_tags += cr_val;
        }
      }
      num_faces += nfaces ;
    }
    cluster = cluster_back;
    
    return num_faces ;
  }
  //no modification, only needed by vogmerge.cc
  bool readBCfromVOG(string filename,
                     vector<pair<int,string> > &boundary_ids) {
    hid_t file_id = 0 ;
    int failure = 0 ; // No failure
    /* Save old error handler */
    herr_t (*old_func)(void*) = 0;
    void *old_client_data = 0 ;
    if(MPI_rank == 0) {
      H5Eget_auto(&old_func, &old_client_data);

      /* Turn off error handling */
      H5Eset_auto(NULL, NULL);

      file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
      if(file_id <= 0) 
        failure = 1 ;

      // Check to see if the file has surface info
      hid_t bc_g = H5Gopen(file_id,"surface_info") ;

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
          hid_t sf_g = H5Gopen(bc_g,buf) ;
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
        H5Eset_auto(old_func, old_client_data);
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
  
  
  //no modification,
  //modify later when entites are stored as 64-bit in file 
  bool readVolTags(hid_t input_fid,
                   vector<pair<string,Loci::entitySet> > &volDat) {
    using namespace Loci ;
    /* Save old error handler */
    herr_t (*old_func)(void*) = 0;
    void *old_client_data = 0 ;
    H5Eget_auto(&old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(NULL, NULL);
    
    vector<pair<string,entitySet> > volTags ;
    hid_t cell_info = H5Gopen(input_fid,"cell_info") ;
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
        hid_t vt_g = H5Gopen(cell_info,buf) ;
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
      hid_t file_info = H5Gopen(input_fid,"file_info") ;
      long numCells = readAttributeLong(file_info,"numCells") ;
      volTags.push_back(pair<string,entitySet>
                        (string("Main"),
                         entitySet(interval(0,numCells-1)))) ;
      H5Gclose(file_info) ;
    }
  
    /* Restore previous error handler */
    H5Eset_auto(old_func, old_client_data);
    volDat.swap(volTags) ;
    return true ;
  }
  
  // turn entitySet into gEntitySet
  gEntitySet setUp(const entitySet &s){
    gEntitySet result;
    FORALL(s, e){
      result += gEntity(e);
    }ENDGFORALL;
    return result;
  }

  //get casename from filename,
  //at current stage, return empty sring
  string get_casename(const string& filename){
    string casename;
    return casename;
  }

  //mpi communication of gEntity and gEntity using type MPI_BYTE
  bool readGridVOG_raw(gfact_db& facts,
                       string filename) {

    string casename = get_casename(filename);

    // create basic key spaces;
    gKeySpaceP universe_space = new gKeySpace();
    universe_space->register_space("UniverseSpace", casename);
    gKeySpaceP node_space = new gKeySpace();
    node_space->register_space("NodeSpace", casename);
    gKeySpaceP face_space = new gKeySpace();
    face_space->register_space("FaceSpace", casename);
    gKeySpaceP cell_space = new gKeySpace();
    cell_space->register_space("CellSpace", casename);
    gKeySpaceP bc_space = new gKeySpace();
    bc_space->register_space("BcSpace", casename);
  
    
    memSpace("Start readGridVOG_raw") ;
   
    //declare containers
    gStore<vector3d<real_t> > pos;
    gMap cl;
    gMap cr;
    gMultiMap face2node;
    vector<pair<int,string> > boundary_ids;
    gEntitySet boundary_faces;
    gEntitySet boundary_taglist;
    vector<pair<string,gEntitySet> > volTags;

    memSpace("get ready to open file") ;                 
    //get ready to read from file
    hid_t file_id = 0 ;
    hid_t face_g = 0 ;
    hid_t node_g = 0 ;
    hid_t dataset = 0 ;
    hid_t dspace = 0 ;
    int failure = 0 ; // No failure


    /* Save old error handler */
    herr_t (*old_func)(void*) = 0;
    void *old_client_data = 0 ;


    gEntity nnodes = 0 ;//nnodes type changed
    if(MPI_rank == 0) {
      
      H5Eget_auto( &old_func, &old_client_data);
      
      /* Turn off error handling */
      H5Eset_auto( NULL, NULL);
      memSpace("open file "+ filename) ;
      file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
      if(file_id <= 0) 
        failure = 1 ;
      
      face_g = H5Gopen(file_id,"face_info") ;
      node_g = H5Gopen(file_id,"node_info") ;
      

      // Check to see if the file has surface info
      hid_t bc_g = H5Gopen(file_id,"surface_info") ;
      // If surface info exists, then read in boundary_ids
      if(bc_g > 0) {
        memSpace("Start read surface_info") ;
        hsize_t num_bcs = 0 ;
        H5Gget_num_objs(bc_g,&num_bcs) ;
        for(hsize_t bc=0;bc<num_bcs;++bc) {
          char buf[1024] ;
          memset(buf, '\0', 1024) ;
                   
          H5Gget_objname_by_idx(bc_g,bc,buf,sizeof(buf)) ;
          buf[1023]='\0' ;
          string name = string(buf) ;
          hid_t sf_g = H5Gopen(bc_g,buf) ;
          hid_t id_a = H5Aopen_name(sf_g,"Ident") ;
          int ident = 0 ;
          H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
          H5Aclose(id_a) ;
          H5Gclose(sf_g) ;
          boundary_ids.push_back(pair<int,string>(ident,name)) ;
        }
        H5Gclose(bc_g) ;
      }
      memSpace("finish reading surface_info") ;
      // First read in and disribute node positions...
      //read in a buf tpos first, and then insert the values into gStore
      memSpace("Start read pos") ;
      
      dataset = H5Dopen(node_g,"positions") ;
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = 1 ;
      //read in nnodes    
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      nnodes = size ;
      memSpace("finsih  read pos") ;
    }
   
    
    int fail_state = 0 ;
    MPI_Allreduce(&failure,&fail_state,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(fail_state != 0)
      return false ;
    
   
    MPI_Bcast(&nnodes,sizeof(gEntity),MPI_BYTE,0,MPI_COMM_WORLD) ;//mpi type changed


    
    vector<gEntitySet> local_nodes = g_simple_partition(nnodes,MPI_COMM_WORLD);
    
    int mxsz = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      int lsz = local_nodes[i].size() ;
      mxsz = max(mxsz,lsz) ;
    }
    
    vector<vector3d<double> > tpos(mxsz) ;
        
    if(MPI_rank == 0) { // read in node positions, send to other processors
      memSpace("Read in Pos") ;
      // read processor zero section first
      
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
                          &tpos[0]) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
        FATAL(err < 0) ;
        Loci::Abort() ;
      }
      H5Sclose(memspace) ;
      
      //put data into container pos
      pos.insert(local_nodes[MPI_rank], &tpos[0]);
      
      
      // now read in remaining processor segments and send to corresponding
      // processor
      for(int i=1;i<MPI_processes;++i) {
        // read in remote processor data
        int sz = local_nodes[i].size() ;
        if(sz == 0) {
          cerr << "sending a zero sized block" << endl ;
        }
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
      // Receive nodes from previous processor
      FATAL(local_nodes[MPI_rank].num_intervals()!=1) ;
      int size = local_nodes[MPI_rank].size() ;
      MPI_Status status ;
      MPI_Recv(&tpos[0],size*3,MPI_DOUBLE,MPI_rank-1,0,MPI_COMM_WORLD,&status) ;
      pos.insert(local_nodes[MPI_rank], &tpos[0]);
      // Shift remaining ones to next processor
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
    gEntity nclusters = 0 ;
    if(MPI_rank == 0) {
      dataset = H5Dopen(face_g,"cluster_sizes") ;
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = 1 ;
      
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      nclusters = size ;
    }
    MPI_Bcast(&nclusters,sizeof(gEntity),MPI_BYTE,0,MPI_COMM_WORLD) ;//mpi type changed

      
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

    // Read in volume tag information and broadcast it
    vector<pair<string,Loci::gEntitySet> > volDat ;
    if(MPI_rank == 0) {
      vector<pair<string,Loci::entitySet> > temp_volDat ;
      readVolTags(file_id,volDat) ;
      for(size_t i = 0; i < temp_volDat.size(); i++){
        volDat.push_back(std::pair<string, Loci::gEntitySet>(temp_volDat[i].first, setUp(temp_volDat[i].second)));
      }
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

      gEntity *ibuf = new gEntity[nivals*2] ;
      if(MPI_rank == 0) 
        for(int j=0;j<nivals;++j) {
          ibuf[j*2]= volDat[i].second[j].first ;
          ibuf[j*2+1]= volDat[i].second[j].second ;
        }
      MPI_Bcast(ibuf,nivals*2*sizeof(gEntity),MPI_BYTE,0,MPI_COMM_WORLD) ;
      gEntitySet set ;
      for(int j=0;j<nivals;++j) 
        set += std::pair<gEntity, gEntity>(ibuf[j*2],ibuf[j*2+1]) ;
      if(MPI_rank != 0) 
        volDat.push_back(pair<string,gEntitySet>(name,set)) ;
      delete[] ibuf ;
    }
   
    volTags.swap(volDat) ;
    if(MPI_rank == 0) {
      H5Gclose(face_g) ;
      H5Fclose(file_id) ;
      /* Restore previous error handler */
      H5Eset_auto(old_func, old_client_data);
    }
    MPI_Allreduce(&failure,&fail_state,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(fail_state != 0)
      return false ;
    
    memSpace("before unpacking clusters") ;
    vector<long> cluster_offset(cluster_sizes.size()+1) ;
    cluster_offset[0] = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i)
      cluster_offset[i+1] = cluster_offset[i] + cluster_sizes[i] ;

    gEntity tot_faces = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      size_t nfaces = getClusterNumFaces(&cluster_info[cluster_offset[i]]) ;
      tot_faces += nfaces ;
    }

    // Now get a face allocation for each processor
    vector<gEntity> faces_pp(MPI_processes,0) ;
    MPI_Allgather(&tot_faces, 1,MPI_GENTITY_TYPE,&faces_pp[0],1,MPI_GENTITY_TYPE,
                  MPI_COMM_WORLD) ;

    gEntity face_base = 0 ;
    for(int i = 0; i < MPI_rank; ++i) {
      face_base += faces_pp[i] ;
    }
   

    gEntity max_cell = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      size_t nfaces = g_fillFaceInfo(&cluster_info[cluster_offset[i]],
                                     face2node,cl,cr,face_base,
                                     max_cell, boundary_faces, boundary_taglist) ;
      face_base += nfaces ;
    }
    gEntity global_max_cell = g_GLOBAL_MAX<gEntity>(max_cell);
    gEntity ncells = global_max_cell+1 ;
   
    // Do not need all collect
    boundary_taglist = g_all_collect_entitySet<gEntity>(boundary_taglist) ;
    gEntitySet all_boundary_faces = g_all_collect_entitySet<gEntity>(boundary_faces) ;


    
    memSpace("after unpacking clusters") ;
   
  
    //facts need create
    facts.create_fact("pos", pos, node_space);
  
    facts.create_fact("cl", cl, face_space, cell_space);
    facts.create_fact("cr", cr, face_space, cell_space);
    facts.create_fact("face2node", face2node, face_space, node_space);
    gParam<gEntity> num_cells;
    *num_cells = ncells;
    facts.create_fact("num_cells", num_cells);
    
    gConstraint bfaces;
    bfaces = all_boundary_faces;
    facts.create_fact("boundary_faces", bfaces, face_space);

    gConstraint bcSurf;
    bcSurf = boundary_taglist;
    facts.create_fact("bcSurf", bcSurf, bc_space);
   
    //only process 0 have the whole store
    gStore<string> boundary_names;
    gStore<string> boundary_tags;
    if(MPI_rank==0){
      FORALL(boundary_taglist, ii) {
        int_type bc = ii;
        char buf[128] ;
        bzero(buf,128) ;
        snprintf(buf,127,"BC_%d",-bc) ;
        boundary_tags.insert(ii, string(buf)) ;
        int id = -bc;
        string bname = string(buf); 
        for(size_t i=0;i<boundary_ids.size();++i) {
          if(boundary_ids[i].first == id){
            bname = boundary_ids[i].second;
            break;
          }
        }
        boundary_names.insert(ii, bname);
      }ENDFORALL ;
    }
    facts.create_fact("boundary_names", boundary_names,bc_space);
    facts.create_fact("boundary_tags", boundary_tags, bc_space);

    
    //add allTags parameter to store all the volTag names
    string all_tags;
    gEntitySet all_tagged_cells;
    for(size_t i=0;i<volTags.size();++i) {
      gParam<string> Tag ;
      *Tag = volTags[i].first ;
      all_tags += volTags[i].first + ",";
      all_tagged_cells += volTags[i].second;
      Tag.set_entitySet(volTags[i].second) ;
      ostringstream oss ;
      oss << "volumeTag(" << volTags[i].first << ")" ;
      facts.create_fact(oss.str(),Tag, cell_space) ;
    }
    gParam<string > allTag ;
    *allTag = all_tags;
    allTag.set_entitySet(all_tagged_cells);
    facts.create_fact("allTags",allTag, cell_space) ;
    debugout<< " all tags cells " << all_tagged_cells<< endl;
    memSpace("Finish reading from file");
    return true ;
  }
  namespace{
    vector<string> split(string str, char delim){
      vector<string> result;
      std::istringstream ss(str);
      string token;
      while(std::getline(ss, token, delim)){
        result.push_back(token);
      }
      return result;
    }
  }
 

  //this function parses the user provided string par_str
  //and creates the gParams for partition functions to use 
  void get_partition_strategy(string par_str, parStrategy& s){

    string casename;
    
     
     
    gKeySpaceP cell_space = gKeySpace::get_space("CellSpace", casename);
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    gKeySpaceP node_space  = gKeySpace::get_space("NodeSpace", casename);

    // stragetgy one, test success on NP = 5, 4, 3
    s.space1 = cell_space;
    s.space2 = face_space;
    s.space3 = node_space; 

    s.strategy = "Metis";
    s.map1 = "cl";//use cl and cr to get self-map
    s.map2 = "cl";//use cl and cr to partition space2
    s.from_space1 = false;//space3 is partitioned to best match space2 partition
    s.map3 = "face2node";//use face2node to partition space3

    //The following test is NP = 5
    // //stragetgy two,
    //  //Test /simcenter/data1/qxue/chem/OBJ/quickTest/Inviscid/highSpeed fixedMass FAILED
    //  //Test /simcenter/data1/qxue/chem/OBJ/quickTest/Inviscid/highSpeed wedge FAILED
    //  //Test /simcenter/data1/qxue/chem/OBJ/quickTest/Inviscid/lowSpeed fixedMass FAILED
   

    // s.space1 = face_space;
    // s.space2 = cell_space;
    // s.space3 = node_space;
    

    // s.strategy = "Orb";
    // s.map1 = "facecenter";//use facecenter to partition space1
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = true; //space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3

    //  // stragetgy three, test success
    // s.space1 = cell_space;
    // s.space2 = face_space;
    // s.space3 = node_space;
    

    // s.strategy = "Orb";
    // s.map1 = "cellcenter";//use facecenter to partition space1
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = false; //space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3

    // // stragetgy four, test success
    // s.space1 = cell_space;
    // s.space2 = face_space;
    // s.space3 = node_space;
    

    // s.strategy = "Simple";
    // s.map1 = "cellcenter";//use facecenter to partition space1
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = false; //space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3


    // //// stragetgy five, test failed, old version failed on NP=4
    // s.space1 = face_space;
    // s.space2 = cell_space;
    // s.space3 = node_space;
    

    // s.strategy = "Simple";
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = true; // space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3


    // //stragetgy six, test success
    // s.space1 = cell_space;
    // s.space2 = face_space;
    // s.space3 = node_space;
    

    // s.strategy = "Simple";
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = false; //space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3



    // // stragetgy seven, test success
    // s.space1 = face_space;
    // s.space2 = cell_space;
    // s.space3 = node_space;
    

    // s.strategy = "Metis";
    // s.map1 = "cl"; //use cl and cr to get self-map
    // s.map2 = "cl"; //use cl and cr to partition space2
    // s.from_space1 = true; // space3 is partitioned to best match space1 partition
    // s.map3 = "face2node"; //use face2node to partition space3


    // //  stragetgy eight, test failed
    // s.space1 = node_space;
    // s.space2 = face_space;
    // s.space3 = cell_space;
    

    // s.strategy = "Simple";
    // s.map2 = "face2node"; //use face2node to partition space2
    // s.from_space1 = false; // space3 is partitioned to best match space1 partition
    // s.map3 = "cl"; //use cl and cr to partition space3

    // //  stragetgy nine, test failed
    //  s.space1 = node_space;
    //  s.space2 = face_space;
    //  s.space3 = cell_space;
    

    //  s.strategy = "Orb";
    //  s.map1 = "pos";
    //  s.map2 = "face2node"; //use face2node to partition space2
    //  s.from_space1 = false; // space3 is partitioned to best match space1 partition
    //  s.map3 = "cl"; //use cl and cr to partition space3

     

  }
  
  //after file is read in, this function collects domain of variables,
  //and then sets up the key partition for each space
  void set_file_key_ptn(gfact_db& facts, const std::string& casename ){
     
   
    // node space
    gStoreRepP posRep = facts.get_fact("pos");
    gKeySpaceP nodeSpace  = gKeySpace::get_space("NodeSpace", casename);
    gEntitySet node_dom = posRep->domain();
    // nodeSpace->set_key_ptn(g_simple_partition<gEntity>(g_all_collect_entitySet<gEntity>(node_dom),
    //  nodeSpace->get_mpi_comm()));
    nodeSpace->set_key_ptn(g_all_collect_vectors<gEntity>(node_dom));

    //face space
    gStoreRepP clRep = facts.get_fact("cl");
    gKeySpaceP faceSpace = gKeySpace::get_space("FaceSpace", casename);
    gEntitySet face_dom = clRep->domain();
    gEntitySet all_faces = g_all_collect_entitySet<gEntity>(face_dom);
    // faceSpace->set_key_ptn(g_simple_partition<gEntity>(all_faces,faceSpace->get_mpi_comm()));
    faceSpace->set_key_ptn(g_all_collect_vectors<gEntity>(face_dom));

      
    //cell space
    gParam<gEntity> num_cells;
    num_cells = facts.get_fact("num_cells");
    gKeySpaceP cellSpace = gKeySpace::get_space("CellSpace", casename);
    cellSpace->set_key_ptn(g_simple_partition<gEntity>(*num_cells,cellSpace->get_mpi_comm()));
      
    //bc space
    gStoreRepP bnameRep = facts.get_fact("boundary_names");
    gKeySpaceP bcSurfSpace = gKeySpace::get_space("BcSpace", casename);
    gEntitySet bc_dom = bnameRep->domain();
    //bcSurfSpace->set_key_ptn(g_all_collect_vectors<gEntity>(bc_dom));
    bcSurfSpace->set_key_ptn(g_simple_partition<gEntity>(g_all_collect_entitySet<gEntity>(bc_dom),
                                                         bcSurfSpace->get_mpi_comm()));
      
  }
  

   
    
   

  //remap all out_vars in space using m
  void remap_space(gKeySpaceP space, gfact_db& facts, const gMap& m){
    variableSet vset =  space->get_out_vars(); 
    if(vset != EMPTY){
      for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
        variable v = *vi;
        gStoreRepP vRep = facts.get_fact(v);
        facts.update_fact(v,vRep->remap(m));
      }
    }
  }
  //This funntion remap all out_vars of space using its own f2g map
  void remap_space_f2g(gfact_db& facts, gKeySpaceP space){
    variableSet vset =  space->get_out_vars(); 
    if(vset != EMPTY){
      gMap f2g;
      f2g = space->get_f2g_map(); 
      for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
        variable v = *vi;
        gStoreRepP vRep = facts.get_fact(v);
        facts.update_fact(v,vRep->remap(f2g));
      }
    }
  }
  //each key space redistributes all its out_vars according to its send_ptn,
  //and then maps the result from file numbering to global numbering 
  void redistribute_remap_all(gfact_db& facts){
    vector<gKeySpaceP> spaces = gKeySpace::get_all_spaces();
    for( unsigned int i = 0; i< spaces.size(); i++){
      variableSet vset = spaces[i]->get_out_vars(); 
      if(vset != EMPTY){
        vector<gEntitySet> ptn = spaces[i]->get_send_ptn(); 
        gMap f2g;
        f2g = spaces[i]->get_f2g_map(); 
        for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
          variable v = *vi;
          gStoreRepP vRep = facts.get_fact(v);
          facts.update_fact(v,vRep->redistribute(ptn, f2g));
        }
      }
    }
    return;
  }
  
  //each key space redistributes all its out_vars according to its key_ptn
  //exclude universeSpace
  void split_redistribute_all(gfact_db& facts){
    vector<gKeySpaceP> spaces = gKeySpace::get_all_spaces(); 
    for(unsigned int i = 0; i < spaces.size(); i++){
      variableSet vset = spaces[i]->get_out_vars();
      if(vset != EMPTY){
        vector<gEntitySet> ptn = spaces[i]->get_key_ptn();
        MPI_Comm comm = spaces[i]->get_mpi_comm();
        for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
          variable v = *vi;
          gStoreRepP vRep = facts.get_fact(v);
          facts.update_fact(v,vRep->split_redistribute(ptn, comm));
        }
      }
    }
  }

  //redistribute all out_vars in space using its send_ptn 
  void redistribute_space(gfact_db& facts, gKeySpaceP space){
    variableSet vset = space->get_out_vars(); 
    if(vset != EMPTY){
      vector<gEntitySet> ptn = space->get_send_ptn();
      MPI_Comm comm = space->get_mpi_comm();
      for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
        variable v = *vi;
        gStoreRepP vRep = facts.get_fact(v);
        facts.update_fact(v,vRep->redistribute(ptn,comm));
      }
    }
  }
  
  //each key space recomposes its their in_vars using its f2g map
  void compose_all(gfact_db& facts){
    vector<gKeySpaceP> spaces = gKeySpace::get_all_spaces(); 
    for( unsigned int i = 0; i < spaces.size(); i++){
      variableSet vset = spaces[i]->get_in_vars();
      if(vset != EMPTY){
        gMap f2g;
        f2g = spaces[i]->get_f2g_map();
        MPI_Comm comm = spaces[i]->get_mpi_comm();
        for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
          variable v = *vi;
          gMapRepP vRep = static_cast<gMapRepP>(facts.get_fact(v));
          vRep->inplace_compose(f2g, comm);
        }
      }
    }
   
  }
  /*
    After partition, this method will redistribute all containers
    and then remap the resulted containers from file numbering to global numbering
    faces are sorted according to the lager cell number each face is associated to,
    and then renumbered.
    the image field of maps are then mapped from file numbering to global numbering
  */
  void remapGrid(gfact_db &facts) {
    string casename;
    
    //get key manager and key spaces
    gKeyManagerP key_manager = facts.get_key_manager();
    
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    gKeySpaceP node_space = gKeySpace::get_space("NodeSpace", casename);
    gKeySpaceP cell_space = gKeySpace::get_space("CellSpace", casename);
    gKeySpaceP bc_space = gKeySpace::get_space("BcSpace", casename);
    
    //when recv_ptn is known,  generate globals keys
    node_space->generate_key(key_manager);
    face_space->generate_key(key_manager);
    cell_space->generate_key(key_manager);
    bc_space->generate_key(key_manager);
    
   
  
   
    //all key spaces redistribute all its out_vars according to its send_ptn,
    //and then map the result from file numbering to global numbering 
    redistribute_remap_all(facts);

    //Next, renumber the faces
    //get variables cl, cr
    gMap cl, cr;
    cl = facts.get_fact("cl");
    cr = facts.get_fact("cr");
    gEntitySet faces = face_space->get_my_keys();//global number
    
    // sort faces

    using std::pair ;
    vector<pair<gEntity,gEntity> > interior_sortlist; //pair<larger cell, face>
    vector<pair<gEntity,gEntity> > boundary_sortlist; //pair<larger cell, face>
    //for boundary faces 
    {
      gMap::const_iterator itr1 = cl.begin();
      gMap::const_iterator itr2 = cr.begin();
      for(; itr1 != cl.end(); ){
        if(itr2->second <0){
          gEntity minc = itr2->second ;//associate each face with the cell with larger cell number
          boundary_sortlist.push_back(pair<gEntity,gEntity>(minc,itr1->first)) ;
        }else{
          gEntity minc = max(itr1->second,itr2->second) ;//associate each face with the cell with larger cell number
          interior_sortlist.push_back(pair<gEntity,gEntity>(minc,itr1->first)) ;
        }
        itr1++;
        itr2++;
      }
    }
   
    sort(interior_sortlist.begin(),interior_sortlist.end(),fieldSort1_unique<gEntity>) ;//sort the faces according cell number
    sort(boundary_sortlist.begin(),boundary_sortlist.end(),fieldSort1_unique<gEntity>) ;//sort the faces according cell number
       
    //faces are renumbered according to the order after sorting
    gMap convert; //f2g map allocated on the domain after redistribution

    gEntitySet::const_iterator ei = faces.begin();
    for(unsigned int i = 0; i< interior_sortlist.size(); i++){
      convert.insert(interior_sortlist[i].second, *ei);
      ei++;
    }
    for(unsigned int i = 0; i< boundary_sortlist.size(); i++){
      convert.insert(boundary_sortlist[i].second, *ei);
      ei++;
    }

    convert.local_sort();
    gMap face_g2f;
    face_g2f = face_space->get_g2f_map();
    face_space->set_g2f_map(face_g2f.remap(convert));//update g2f map of face_space
    remap_space(face_space,facts,convert);//remap cl, cr, and face2node

    compose_all(facts);
    //cr is a special case,
    gMap bc_f2g;
    bc_f2g = bc_space->get_f2g_map();
    MPI_Comm comm = bc_space->get_mpi_comm();
    gMapRepP(facts.get_variable("cr"))->inplace_compose(bc_f2g, comm);
    debugout << "finish with remapGrid" << endl; 
  }

 
  //Description: Reads grid structures in the fact database
  //Input: facts and grid file name
  //Output: true if sucess

  //steps:
  //1. read in the raw data from file, create facts. 
  //2. set up the distribution: each keySpace set_key_ptn
  //3. redistribute, each keySpace redistribute its outgoing variables( the variable whose domain is this keySpace)
  // 4. generalize the metis partition, the user specify which space to partition first, which map is used for partition
  // 5. generalize the affinity partition, the user specify which space to partition, which map is used for partition,

  //??  for larger number of process, the information for partition should not be stored in a vector<gEntitySet>,
  //     instead, a process map (i.e. vector<int>) should be used. But in pack and unpack function, gEntitySet is needed. 
  //     so currently, in gKeySpace, all partition info is still in the form  vector<gEntitySet>

  //6. after partition, remapGrid is called to redistribute and remap all containers to global numbering
  
  bool readFVMGrid(gfact_db &facts, string filename) {
    string casename = get_casename(filename);
    double t1 = MPI_Wtime() ;
    //	timer_token read_file_timer = new timer_token;
    //	if(collect_perf_data)
    //		read_file_timer = perfAnalysis->start_timer("Reading in FVM Grid");
    memSpace("readFVMGrid Start") ;
    bool useVOG = true ;
    string input_file = filename ;
    string::size_type spos = string::npos ;
    if((spos = input_file.rfind('.')) != string::npos) {
      input_file.erase(0,spos) ;
      if(input_file == ".xdr")
        useVOG = false ;
    }
    if(useVOG) {
      if(!readGridVOG_raw(facts, filename))
        return false;
    } else {
      cerr << "XDR Grid File Reading Support has been removed!!!, convert file to VOG format." << endl ;
      Loci::Abort() ;
      return false;
    }

    //need redistribute and sort here
   
    set_file_key_ptn(facts, casename);
    if(MPI_processes > 1)split_redistribute_all(facts);
         
    memSpace("after reading grid") ;

     
    if(MPI_processes == 1) {
      //set keys, my keys
      gKeyManagerP key_manager = facts.get_key_manager();
      gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
      gKeySpaceP node_space = gKeySpace::get_space("NodeSpace", casename);
      gKeySpaceP cell_space = gKeySpace::get_space("CellSpace", casename);
      gKeySpaceP bc_space = gKeySpace::get_space("BcSpace", casename);
      
      node_space->set_send_recv_ptn(node_space->get_key_ptn());
      face_space->set_send_recv_ptn(face_space->get_key_ptn());
      cell_space->set_send_recv_ptn(cell_space->get_key_ptn());
      bc_space->set_send_recv_ptn(bc_space->get_key_ptn());
     

     
      node_space->generate_key(key_manager);
      face_space->generate_key(key_manager);
      cell_space->generate_key(key_manager);
      bc_space->generate_key(key_manager);
       
      remap_space_f2g(facts, node_space);
      remap_space_f2g(facts, face_space);
      remap_space_f2g(facts, cell_space);
      remap_space_f2g(facts, bc_space);
      
      compose_all(facts);
      //cr is a special case,
      gMap bc_f2g;
      bc_f2g = bc_space->get_f2g_map();
      MPI_Comm comm = bc_space->get_mpi_comm();
      gMapRepP(facts.get_variable("cr"))->inplace_compose(bc_f2g, comm);       
      return true ;
    }

    memSpace("before partitioning") ;

   

    string par_str;//this should be from a param<string> in facts 
    parStrategy s;
    get_partition_strategy(par_str,s);
    
    vector<int> procmap; //the result of primary partition
    gParam<string> strategy;
    
    if(s.strategy == "Metis"){
      if(MPI_rank==0)std::cout<<"metis partition used " << endl;
      primary_partition_metis(facts, procmap, s);
    }else if(s.strategy == "Orb"){
      if(MPI_rank==0)std::cout<<"orb partition used " << endl;
      primary_partition_orb(facts, procmap, s);
    }else if(s.strategy == "Simple"){
      if(MPI_rank==0)std::cout<<"simple partition used " << endl;
      s.space1->set_simple_partition(procmap);
    }else{
      cerr<<"ERROR: partition strategy "<< s.strategy <<" unrecognized" << endl;
      cerr<<"       acceptable values are: Metis, Orb, Simple" << endl;
      Loci::Abort();
    }
    debugout<<"start affinity partition " << endl;
    affinity_partition(facts, procmap, s);
       
    gKeySpaceP bc_space = gKeySpace::get_space("BcSpace", casename);
    bc_space->set_simple_partition();
    
    memSpace("after partitioning") ;
    remapGrid(facts);
    memSpace("after remapGridStructures") ;
    
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(read_file_timer);
    double t2 = MPI_Wtime() ;
    debugout << "Time to read in file '" << filename << ", is " << t2-t1
             << endl ;
    memSpace("returning from FVM grid reader") ;
    return true ;
  }

  

  enum matrix_coloring_type {COLOR_DEFAULT, COLOR_DFS} ;

  void create_ref(gfact_db &facts) {
    string casename;
    gMap cr ;
    cr = facts.get_fact("cr") ;
    gConstraint boundary_faces ;
    boundary_faces = facts.get_fact("boundary_faces") ;
    gEntitySet refdom = *boundary_faces ;
    gMap ref ;
    gMap::const_iterator itr = cr.begin();
    GFORALL(refdom,i1) {
      while(itr!= cr.end() && itr->first < i1)itr++;
      while(itr!= cr.end() && itr->first == i1){
        ref.insert(i1, itr->second) ;
        itr++;
      }
    } ENDGFORALL ;

    gKeySpaceP bc_space = gKeySpace::get_space("BcSpace", casename);
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    facts.create_fact("ref",ref, face_space, bc_space) ;
  }

  void create_ghost_cells(gfact_db &facts) {
    string casename;
    gKeySpaceP cell_space = gKeySpace::get_space("CellSpace", casename);
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    gKeyManagerP km = facts.get_key_manager();
    
    gConstraint interior_faces,boundary_faces ;
    gConstraint geom_cells, ghost_cells, cells ;

    geom_cells = cell_space->get_my_keys();
    boundary_faces = facts.get_fact("boundary_faces");
    
    gEntitySet tmp_ghost = cell_space->generate_key(km,(*boundary_faces).size());
    
    gEntitySet::const_iterator ei = tmp_ghost.begin() ;
    gMap cr;
    cr = facts.get_fact("cr");
    gMap::iterator itr = cr.begin();
    GFORALL(*boundary_faces,fc) {
      while(itr != cr.end() && itr->first < fc)itr++;
      while(itr != cr.end() && itr->first == fc) {
        itr->second = *ei;
        ei++;
        itr++;
      }
    } ENDGFORALL ;
    
    ghost_cells = tmp_ghost;
    facts.replace_fact("cr",cr, face_space, cell_space) ;
    
    cells = *geom_cells + *ghost_cells ;

    
    facts.create_fact("geom_cells",geom_cells, cell_space) ;
    facts.create_fact("ghost_cells",ghost_cells, cell_space) ;
    facts.create_fact("cells",cells, cell_space) ;

    Loci::debugout << "geom_cells = " << *geom_cells << endl ;
    Loci::debugout << "ghost_cells = " << *ghost_cells << endl ;
    Loci::debugout << "cells = " << *cells << endl ;
  }

  void create_face_info(gfact_db &facts) {
    string casename;
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    gConstraint faces;
    faces = face_space->get_my_keys();
    facts.create_fact("faces",faces) ;
    gConstraint boundary_faces ;
    boundary_faces = facts.get_variable("boundary_faces");
    gConstraint interior_faces ;
    interior_faces = (*faces-*boundary_faces) ;
    facts.create_fact("interior_faces",interior_faces, face_space) ;
  }

 

  void copy_facts(gfact_db& gfacts, fact_db& facts){
   
    string casename;
    gKeySpaceP face_space = gKeySpace::get_space("FaceSpace", casename);
    gKeySpaceP node_space = gKeySpace::get_space("NodeSpace", casename);
    gKeySpaceP cell_space = gKeySpace::get_space("CellSpace", casename);
    gKeySpaceP bc_space = gKeySpace::get_space("BcSpace", casename);
    gKeySpaceP edge_space = gKeySpace::get_space("EdgeSpace", casename);
    //set up the distribute info in facts
   
   
    { gMap g2f;
      g2f =  node_space->get_g2f_map();
      vector<int> alloc_file;
      alloc_file.resize(g2f.size());
      int ind = 0;
      for(gMap::const_iterator mi = g2f.begin(); mi != g2f.end(); mi++){
        alloc_file[ind++] = mi->second;
      }
      entitySet entities_global = facts.get_distributed_alloc(alloc_file).first;
      g2f.clear();
      alloc_file.resize(0);
    }
    { gMap g2f;
      vector<int> alloc_file;
      g2f =  face_space->get_g2f_map();
      alloc_file.resize(g2f.size());
      int ind = 0;
      for(gMap::const_iterator mi = g2f.begin(); mi != g2f.end(); mi++){
        alloc_file[ind++] = mi->second;
      }
      entitySet entities_global = facts.get_distributed_alloc(alloc_file).first;
      g2f.clear();
      alloc_file.resize(0);
    }
    { gMap g2f;
      vector<int> alloc_file;
      g2f =  cell_space->get_g2f_map();
      gConstraint geom_cells;
      geom_cells = gfacts.get_fact("geom_cells");
      
      
      alloc_file.resize((*geom_cells).size());
      
      int ind = 0;
      gMap::const_iterator itr = g2f.begin();
      GFORALL(*geom_cells, c){
        while(itr!=g2f.end() && itr->first<c)itr++;
        while(itr!=g2f.end() && itr->first == c){
          alloc_file[ind++] = itr->second;
          itr++;
        }
      }ENDGFORALL;
      entitySet entities_global = facts.get_distributed_alloc(alloc_file).first;
      
      g2f.clear();
      alloc_file.resize(0);
    }
    { gMap g2f;
      g2f =  bc_space->get_g2f_map();
      vector<int> alloc_file;
      alloc_file.resize(g2f.size());
      int ind = 0;
      for(gMap::const_iterator mi = g2f.begin(); mi != g2f.end(); mi++){
        alloc_file[ind++] = mi->second;
      }
      entitySet entities_global = facts.get_distributed_alloc(alloc_file).first;
      
      g2f.clear();
      alloc_file.resize(0);
    }
    {
      gConstraint  boundary_faces;
      boundary_faces =  gfacts.get_fact("boundary_faces");
      std::pair<entitySet, entitySet> ghost_pair = facts.get_distributed_alloc((*boundary_faces).size()) ;
    }

    if(edge_space != 0){
      gMap g2f;
      g2f =  edge_space->get_g2f_map();
      vector<int> alloc_file;
      alloc_file.resize(g2f.size());
      int ind = 0;
      for(gMap::const_iterator mi = g2f.begin(); mi != g2f.end(); mi++){
        alloc_file[ind++] = mi->second;
      }
      entitySet entities_global = facts.get_distributed_alloc(alloc_file).first;
      g2f.clear();
      alloc_file.resize(0);
    }
      
    gfacts.copy_facts(facts);
  }


  
 
  bool setupFVMGrid(gfact_db &gfacts, string filename) {

    string casename = get_casename(filename);
    if(!readFVMGrid(gfacts,filename))
      return false ;
    debugout << " before create face info " << endl;
    create_face_info(gfacts) ;
    debugout<<"finish create face info " << endl;
    create_ref(gfacts) ;
    debugout<<"finish create ref " << endl;
    create_ghost_cells(gfacts) ;
    debugout<<"finish create ghost cells " << endl;
    debugout<<"finish setupFVMGrid " << endl;
    return true ;
  }
}



#define VERBOSE
#define DEBUG
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

#include <Tools/xdr.h>

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

#include <malloc.h>

extern "C" {
  typedef int idxtype ;
  void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
}

namespace Loci {

  void memSpace(string s) {
    //#define DIAG
#ifdef DIAG
    struct mallinfo info = mallinfo() ;
    debugout << s << ": minfo, arena=" << info.arena
             << ", ordblks=" << info.ordblks
             << ", hblks="<< info.hblks
             << ", hblkhd="<< info.hblkhd
             << ", uordblks=" << info.uordblks
             << ", fordblks="<<info.fordblks
             << ", keepcost="<<info.keepcost << endl ;
    debugout.flush() ;
#endif
  }
  extern bool use_simple_partition ;
  //Following assumption is needed for the next three functions
  //Assumption: We follow the convention that boundary cells are always on right side
  //of a face.

  //Input: Mapping from faces to its right cells.
  //Output: Entities of boundary cells.
  entitySet getBoundaryCells(const MapRepP tmp_cr_rep) {
    entitySet cri = tmp_cr_rep->image(tmp_cr_rep->domain()) ;
    return(cri & interval(Loci::UNIVERSE_MIN,-1)) ;
  }

  //Input: Mapping from faces to its right cells.
  //Output: Entities for the boundary faces in the domain of given map.
  entitySet getBoundaryFaces(const MapRepP tmp_cr_rep) {
    entitySet boundary_cells = getBoundaryCells(tmp_cr_rep);
    return(tmp_cr_rep->preimage(boundary_cells).first);
  }

  //Input: Mapping from faces to its right cells.
  //Output: Entities for the interior faces in the domain of given map.
  entitySet getInteriorFaces(const MapRepP tmp_cr_rep) {
    return(tmp_cr_rep->domain() - getBoundaryFaces(tmp_cr_rep)) ;
  }

  //Description: Reads grid structures from grid file in the .xdr format.
  //Input: file name and max_alloc (starting of entity assignment - node base)
  //Output:
  // local_cells, local_nodes, local_faces: partition of nodes, faces, cells
  // pos: position of nodes
  // cl: Mapping from face to cell on the left side
  // cr: Mapping from face to cell on the right side
  // face2node: MultiMapping from a face to nodes
  bool readGridXDR(vector<entitySet> &local_nodes,
		   vector<entitySet> &local_faces,
		   vector<entitySet> &local_cells,
		   store<vector3d<real_t> > &pos, Map &cl, Map &cr,
		   multiMap &face2node, int max_alloc, string filename) {

    memSpace("begin readGridXDR") ;
    local_nodes.resize(Loci::MPI_processes);
    local_faces.resize(Loci::MPI_processes);
    local_cells.resize(Loci::MPI_processes);

    // First read in header information
    // Note:  Only processor 0 reads the file, it then sends the results
    // to other processors.
    int ndm ;
    int npnts, nfaces, ncells ;
    int dummy1,dummy2 ;
    FILE* FP = 0 ;
    XDR xdr_handle ;
    if(Loci::MPI_rank == 0 ) {
      FP = fopen(filename.c_str(), "r") ;
      if(FP == NULL)
        return false ;

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
      int data[3] ;
      data[0] = npnts ;
      data[1] = nfaces ;
      data[2] = ncells ;

      if(Loci::MPI_processes > 1)
	MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD) ;
    }
    else {
      int data[3] = {0,0,0} ;
      MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD) ;
      npnts = data[0] ;
      nfaces = data[1] ;
      ncells = data[2] ;
    }

    // Create initial allocation of nodes, faces, and cells
    int node_ivl = npnts / Loci::MPI_processes;
    int node_ivl_rem = npnts % Loci::MPI_processes ;
    int node_accum = 0 ;
    int face_ivl = nfaces / Loci::MPI_processes;
    int face_ivl_rem = nfaces % Loci::MPI_processes ;
    int face_accum = 0 ;
    int cell_ivl = ncells / Loci::MPI_processes;
    int cell_ivl_rem = ncells % Loci::MPI_processes ;
    int cell_accum = 0 ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int nodes_base = max_alloc ;
      int faces_base = max_alloc+npnts ;
      int cells_base = max_alloc+npnts+nfaces ;
      int j = Loci::MPI_processes - i - 1 ;
      int node_accum_update = node_accum + node_ivl + ((j<node_ivl_rem)?1:0) ;
      int face_accum_update = face_accum + face_ivl + ((j<face_ivl_rem)?1:0) ;
      int cell_accum_update = cell_accum + cell_ivl + ((j<cell_ivl_rem)?1:0) ;

      if(i == Loci::MPI_processes-1) {
	local_nodes[i] = interval(nodes_base + node_accum,
                                  nodes_base + npnts - 1) ;
	local_faces[i] = interval(faces_base + face_accum,
                                  faces_base + nfaces-1) ;
	local_cells[i] = interval(cells_base + cell_accum,
				  cells_base + ncells-1) ;
      }
      else {
	local_nodes[i] = interval(nodes_base + node_accum,
                                  nodes_base + node_accum_update - 1) ;
	local_faces[i] = interval(faces_base + face_accum,
                                  faces_base + face_accum_update - 1) ;
	local_cells[i] = interval(cells_base + cell_accum,
                                  cells_base + cell_accum_update - 1) ;
      }
      node_accum = node_accum_update ;
      face_accum = face_accum_update ;
      cell_accum = cell_accum_update ;
    }

    memSpace("read file, before node redistribution") ;
    // Distribute positions
    if(Loci::MPI_rank == 0) {
      double *tmp_pos ;
      tmp_pos = new double[local_nodes[Loci::MPI_processes-1].size() * 3] ;
      for(int i = 0; i < 3*local_nodes[Loci::MPI_rank].size(); ++i)
      	if(!xdr_double(&xdr_handle, &tmp_pos[i]))
          return false ;

      int tmp = 0 ;
      pos.allocate(local_nodes[Loci::MPI_rank]) ;
      for(entitySet::const_iterator ei = local_nodes[Loci::MPI_rank].begin(); ei != local_nodes[Loci::MPI_rank].end(); ++ei) {
	vector3d<real_t> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
	tmp += 3 ;
	pos[*ei] = t ;
      }

      for(int i = 1; i < Loci::MPI_processes; ++i) {
	for(int j = 0; j < 3*local_nodes[i].size(); ++j)
	  if(!xdr_double(&xdr_handle, &tmp_pos[j]))
            return false ;
        MPI_Send(tmp_pos,local_nodes[i].size()*3,MPI_DOUBLE, i, 9,
                 MPI_COMM_WORLD) ;
      }
      delete [] tmp_pos ;

    }
    else {
      MPI_Status status ;
      int recv_count = local_nodes[Loci::MPI_rank].size()*3 ;
      double *tmp_pos = new double[recv_count] ;
      MPI_Recv(tmp_pos,recv_count,MPI_DOUBLE,0,9,MPI_COMM_WORLD,&status) ;

      int tmp = 0 ;
      pos.allocate(local_nodes[Loci::MPI_rank]) ;
      for(entitySet::const_iterator ei = local_nodes[Loci::MPI_rank].begin(); ei != local_nodes[Loci::MPI_rank].end(); ++ei) {
	vector3d<real_t> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
	tmp += 3 ;
	pos[*ei] = t ;
      }

      delete [] tmp_pos ;
    }
    cl.allocate(local_faces[Loci::MPI_rank]) ;
    cr.allocate(local_faces[Loci::MPI_rank]) ;

    store<int> local_count ;
    local_count.allocate(local_faces[Loci::MPI_rank]) ;
    std::vector<int> offset ;
    std::vector<int> start_off(Loci::MPI_processes+1) ;
    start_off[0] = 0 ;
    if(Loci::MPI_rank == 0) {
      int *off_cl_cr ;
      off_cl_cr = new int[3*local_faces[Loci::MPI_processes-1].size() + 1] ;

      for(int i = 0; i < (local_faces[0].size() * 3) + 1; ++i)
	if(!xdr_int(&xdr_handle, &off_cl_cr[i]))
          return false ;
      int tmp = 0 ;
      for(entitySet::const_iterator ei = local_faces[0].begin(); ei != local_faces[0].end(); ++ei) {
	offset.push_back(off_cl_cr[tmp++]) ;
	cl[*ei] = off_cl_cr[tmp++] ;
	cr[*ei] = off_cl_cr[tmp++] ;
	if(cl[*ei] < 0)
	  if(cr[*ei] < 0) {
	    cerr << " boundary condition on both sides of a face?" << endl ;
	    exit(1) ;
	  } else {
	    int tmp_swap = cr[*ei] ;
	    cr[*ei] = cl[*ei] ;
	    cl[*ei] = tmp_swap ;
	  }
	cl[*ei] += max_alloc + npnts + nfaces - 1 ;
	if(cr[*ei] > 0)
	  cr[*ei] += max_alloc + npnts + nfaces - 1 ;
      }
      offset.push_back(off_cl_cr[tmp]) ;
      entitySet::const_iterator ii = local_faces[Loci::MPI_rank].begin() ;
      for(size_t i = 1; i < offset.size(); ++i) {
	local_count[*ii] = offset[i] - offset[i-1] ;
	++ii ;
      }

      int init_off = off_cl_cr[tmp] ;
      for(int i = 1; i < Loci::MPI_processes; ++i) {
	off_cl_cr[0] = init_off ;
	start_off[i] = init_off ;
	for(int j = 1; j < (local_faces[i].size() * 3) +1; ++j)
	  if(!xdr_int(&xdr_handle, &off_cl_cr[j]))
            return false ;
	int send_size = local_faces[i].size() * 3 + 1 ;
	MPI_Send(off_cl_cr, send_size, MPI_INT, i, 10, MPI_COMM_WORLD) ;
	init_off = off_cl_cr[local_faces[i].size()*3] ;
	if(i==Loci::MPI_processes-1)
	  start_off[Loci::MPI_processes] = init_off ;
      }
      delete [] off_cl_cr ;
    } else {
      MPI_Status status ;
      int recv_count = local_faces[Loci::MPI_rank].size() * 3 + 1 ;
      int *off_cl_cr = new int[3*local_faces[Loci::MPI_rank].size() + 1] ;
      MPI_Recv(off_cl_cr, recv_count, MPI_INT, 0, 10, MPI_COMM_WORLD, &status) ;
      int tmp = 0 ;
      for(entitySet::const_iterator ei = local_faces[Loci::MPI_rank].begin(); ei != local_faces[Loci::MPI_rank].end(); ++ei) {
	offset.push_back(off_cl_cr[tmp++]) ;
	cl[*ei] = off_cl_cr[tmp++] ;
	cr[*ei] = off_cl_cr[tmp++] ;
	if(cl[*ei] < 0)
	  if(cr[*ei] < 0) {
	    cerr << "2 boundary condition on both sides of a face?" << endl ;
	    exit(1) ;
	  } else {
	    int tmp = cr[*ei] ;
	    cr[*ei] = cl[*ei] ;
	    cl[*ei] = tmp ;
	  }
	cl[*ei] += max_alloc + npnts + nfaces - 1 ;
	if(cr[*ei] > 0)
	  cr[*ei] += max_alloc + npnts + nfaces - 1 ;
      }
      offset.push_back(off_cl_cr[tmp]) ;
      entitySet::const_iterator ii = local_faces[Loci::MPI_rank].begin() ;
      for(size_t i = 1; i < offset.size(); ++i) {
	local_count[*ii] = offset[i] - offset[i-1] ;
	++ii ;
      }
      delete [] off_cl_cr ;
    }

    face2node.allocate(local_count);
    if(Loci::MPI_processes > 1) {
      if(Loci::MPI_rank == 0) {
	int maxf2nsize = 0 ;
	for(int i=0;i<Loci::MPI_processes;++i)
	  maxf2nsize = max(maxf2nsize,start_off[i+1]-start_off[i]) ;
	int* tmp_f2n = new int[maxf2nsize];
	for(int i = 0;
	    i < (start_off[Loci::MPI_rank+1] - start_off[Loci::MPI_rank]);
	    ++i) {
	  if(!xdr_int(&xdr_handle, &tmp_f2n[i]))
	    return false ;
	}
	int tmp = 0 ;
	for(entitySet::const_iterator ei = local_faces[0].begin(); ei != local_faces[0].end(); ++ei)
	  for(int i = 0; i < local_count[*ei]; ++i)
	    face2node[*ei][i] = tmp_f2n[tmp++] + max_alloc ;
	for(int i = 1; i < Loci::MPI_processes; ++i) {
	  int send_size = start_off[i+1] - start_off[i] ;
	  for(int j = 0; j < send_size; ++j)
	    if(!xdr_int(&xdr_handle, &tmp_f2n[j]))
	      return false ;
	  MPI_Send(tmp_f2n, send_size, MPI_INT, i, 11, MPI_COMM_WORLD) ;
	}
	delete [] tmp_f2n ;
      } else {
	MPI_Status status ;
	int recv_count = offset[offset.size()-1] - offset[0] ;
	int *tmp_f2n = new int[recv_count] ;
	MPI_Recv(tmp_f2n, recv_count, MPI_INT, 0, 11, MPI_COMM_WORLD, &status) ;
	int tmp = 0 ;
	for(entitySet::const_iterator ei = local_faces[Loci::MPI_rank].begin(); ei != local_faces[Loci::MPI_rank].end(); ++ei)
	  for(int i = 0; i < local_count[*ei]; ++i)
	    face2node[*ei][i] = tmp_f2n[tmp++] + max_alloc ;
	delete [] tmp_f2n ;
      }
    }
    else {
      for(entitySet::const_iterator ei = local_faces[0].begin(); ei != local_faces[0].end(); ++ei)
	for(int i = 0; i < local_count[*ei]; ++i) {
	  if(!xdr_int(&xdr_handle, &face2node[*ei][i]))
            return false ;
          face2node[*ei][i] += max_alloc ;
        }
    }

    if(Loci::MPI_rank == 0 ) {
      Loci::debugout <<" All the procs finished reading the files " << endl ;
      xdr_destroy(&xdr_handle) ;
      fclose(FP) ;
    }
    memSpace("returning from grid reader") ;
    return true;
  }

  template<class T> void readVectorDist(hid_t group_id,
                                        const char *vector_name,
                                        vector<int> sizes,
                                        vector<T> &v) {

    hid_t dataset = 0 ;
    hid_t dspace = 0 ;

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
      int lsz = sizes[0] ;

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
        int sz = sizes[i] ;
        if(sz > 0) {
          vector<T> tmp(sz) ;
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

          // send to remote processor
          MPI_Send(&tmp[0],sz*sizeof(T),MPI_BYTE,i,0,MPI_COMM_WORLD) ;
        }
      }
      H5Dclose(dataset) ;
      H5Sclose(dspace) ;

    } else {
      int size = sizes[MPI_rank] ;
      if(size > 0) {
        MPI_Status status ;
        MPI_Recv(&v[0],size*sizeof(T),MPI_BYTE,0,0,MPI_COMM_WORLD,&status) ;
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
		   store<vector3d<real_t> > &pos, Map &cl, Map &cr,
		   multiMap &face2node, int max_alloc, string filename,
                   vector<pair<int,string> > &boundary_ids) {
    local_nodes.resize(Loci::MPI_processes);
    local_faces.resize(Loci::MPI_processes);
    local_cells.resize(Loci::MPI_processes);

    hid_t file_id = 0 ;
    hid_t face_g = 0 ;
    hid_t node_g = 0 ;
    hid_t dataset = 0 ;
    hid_t dspace = 0 ;
    int nnodes = 0 ;
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

      face_g = H5Gopen(file_id,"face_info") ;
      node_g = H5Gopen(file_id,"node_info") ;


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
      }

      // First read in and disribute node positions...
      dataset = H5Dopen(node_g,"positions") ;
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = true ;
        
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
      MPI_Bcast(data,bufsz,MPI_CHAR,0,MPI_COMM_WORLD) ;
      buf = string(data) ;
      istringstream iss(buf) ;
      boundary_ids.clear() ;
      for(int i=0;i<bsz;++i) {
        int id ;
        string name ;
        iss >> id >> name ;
        boundary_ids.push_back(pair<int,string>(id,name)) ;
      }
    }
      
    MPI_Bcast(&nnodes,1,MPI_INT,0,MPI_COMM_WORLD) ;

    // create node allocation
    int npnts = nnodes ;
    int node_ivl = npnts / Loci::MPI_processes;
    int node_ivl_rem = npnts % Loci::MPI_processes ;
    int node_accum = 0 ;
    int nodes_base = max_alloc ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = Loci::MPI_processes - i - 1 ;
      int node_accum_update = node_accum + node_ivl + ((j<node_ivl_rem)?1:0) ;
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
      for(int i=1;i<MPI_processes;++i) {
        // read in remote processor data
        int sz = local_nodes[i].size() ;
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
        MPI_Send(&tpos[0],sz*3,MPI_DOUBLE,i,0,MPI_COMM_WORLD) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(dspace) ;
      H5Gclose(node_g) ;
    } else {
      // Receive nodes from root processor
      FATAL(local_nodes[MPI_rank].num_intervals()!=1) ;
      int start = local_nodes[MPI_rank].Min() ;
      int size = local_nodes[MPI_rank].size() ;
      MPI_Status status ;
      MPI_Recv(&pos[start],size*3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status) ;
    }

    vector<unsigned char> cluster_info ;
    vector<unsigned short> cluster_sizes ;
    // Now read in face clusters
    int nclusters = 0 ;
    if(MPI_rank == 0) {
      dataset = H5Dopen(face_g,"cluster_sizes") ;
      dspace = H5Dget_space(dataset) ;
      if(dataset <=0 || dspace <=0)
        failure = true ;
      
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      nclusters = size ;
    }
    MPI_Bcast(&nclusters,1,MPI_INT,0,MPI_COMM_WORLD) ;

    vector<int> cluster_dist(MPI_processes) ;
    int sum = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      cluster_dist[i] = nclusters/MPI_processes +
        ((nclusters%MPI_processes)>i?1:0);
      sum += cluster_dist[i] ;
    }
    FATAL(sum != nclusters) ;
    readVectorDist(face_g,"cluster_sizes",cluster_dist,cluster_sizes) ;

    int cluster_info_size = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i)
      cluster_info_size += cluster_sizes[i] ;

    MPI_Allgather(&cluster_info_size,1,MPI_INT,&cluster_dist[0],1,MPI_INT,
                  MPI_COMM_WORLD) ;
    readVectorDist(face_g,"cluster_info",cluster_dist,cluster_info) ;

    if(MPI_rank == 0) {
      H5Gclose(face_g) ;
      H5Fclose(file_id) ;
      /* Restore previous error handler */
      H5Eset_auto(old_func, old_client_data);
    }
    MPI_Allreduce(&failure,&fail_state,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(fail_state != 0)
      return false ;

    vector<int> cluster_offset(cluster_sizes.size()+1) ;
    cluster_offset[0] = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i)
      cluster_offset[i+1] = cluster_offset[i] + cluster_sizes[i] ;

    int tot_faces = 0 ;
    for(size_t i=0;i<cluster_sizes.size();++i) {
      int nfaces = getClusterNumFaces(&cluster_info[cluster_offset[i]]) ;
      tot_faces += nfaces ;
    }

    // Now get a face allocation for each processor
    vector<int> faces_pp(MPI_processes) ;
    MPI_Allgather(&tot_faces,1,MPI_INT,&faces_pp[0],1,MPI_INT,
                  MPI_COMM_WORLD) ;

    int face_accum = 0 ;
    int faces_base = max_alloc + nnodes ;
    for(int i = 0; i < MPI_processes; ++i) {
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
    } ENDFORALL ;

    int global_max_cell = max_cell ;
    MPI_Allreduce(&max_cell,&global_max_cell,1,MPI_INT,MPI_MAX,
                  MPI_COMM_WORLD) ;

    int ncells = global_max_cell+1 ;
    int cell_ivl = ncells / MPI_processes;
    int cell_ivl_rem = ncells % MPI_processes ;
    int cell_accum = 0 ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = MPI_processes - i - 1 ;
      int cell_accum_update = cell_accum + cell_ivl + ((j<cell_ivl_rem)?1:0) ;

      if(i == MPI_processes-1) {
	local_cells[i] = interval(cells_base + cell_accum,
				  cells_base + ncells-1) ;
      } else {
	local_cells[i] = interval(cells_base + cell_accum,
                                  cells_base + cell_accum_update - 1) ;
      }
      cell_accum = cell_accum_update ;
    }

    return true ;
  }

  extern void distributed_inverseMap(multiMap &result,
                                     vector<pair<Entity,Entity> > &input,
                                     entitySet input_image,
                                     entitySet input_preimage,
                                     const std::vector<entitySet> &init_ptn) ;

  vector<entitySet> newMetisPartitionOfCells(const vector<entitySet> &local_cells,
                                             const Map &cl, const Map &cr) {


    entitySet dom = cl.domain() & cr.domain() ;
    entitySet::const_iterator ei ;
    int cnt = 0 ;
    for(ei=dom.begin();ei!=dom.end();++ei) {
      if(cl[*ei] > 0 && cr[*ei]>0)
        cnt++ ;
    }
    vector<pair<int,int> > rawMap(cnt*2) ;
    int j = 0 ;
    for(ei=dom.begin();ei!=dom.end();++ei) {
      if(cl[*ei] > 0 && cr[*ei]>0) {
        rawMap[j++] = pair<int,int>(cl[*ei],cr[*ei]) ;
        rawMap[j++] = pair<int,int>(cr[*ei],cl[*ei]) ;
      }
    }

    sort(rawMap.begin(),rawMap.end()) ;

    multiMap cell2cell ;
    entitySet all_cells ;
    for(int i=0;i<MPI_processes;++i)
      all_cells += local_cells[i] ;

    distributed_inverseMap(cell2cell,rawMap, all_cells,all_cells,local_cells) ;

    vector<pair<int,int> >().swap(rawMap) ; // Free up memory from rawMap
    int count = 0 ;
    int size_map = local_cells[Loci::MPI_rank].size() ;
    entitySet dom_map = Loci::interval(0, size_map-1) ;
    store<int> size_adj ;
    size_adj.allocate(dom_map) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_adj[count] = cell2cell.end(*ei)-cell2cell.begin(*ei) ;
      ++count ;
    }

    int *part = new int[size_map] ;
    int *xadj = new int[size_map+1] ;
    int edgecut ;
    int *vdist = new int[Loci::MPI_processes + 1] ;
    int cmin = local_cells[0].Min();
    for(int i = 0; i < Loci::MPI_processes; i++)
      cmin = min(local_cells[i].Min(), cmin);

    edgecut = 0 ;
    xadj[0] = 0 ;
    for(int i = 0; i < size_map; ++i)
      xadj[i+1] = xadj[i] + size_adj[i] ;

    int *adjncy = new int[xadj[size_map]] ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_t sz = cell2cell.end(*ei)-cell2cell.begin(*ei) ;
      for(size_t i = 0; i != sz; ++i)        {
	adjncy[count] = cell2cell[*ei][i] - cmin ;
	count ++ ;
      }
    }
    cell2cell.setRep(multiMap().Rep()) ;// Free up memory from multiMap

    vdist[0] = 0 ;
    for(int i = 1; i <= Loci::MPI_processes; ++i)
      vdist[i] = vdist[i-1] + local_cells[i-1].size() ;

#ifndef MPI_STUBB
    MPI_Comm mc = MPI_COMM_WORLD ;
    int num_partitions = Loci::MPI_processes ;
    int wgtflag = 0 ;
    int numflag = 0 ;
    int options = 0 ;
    ParMETIS_PartKway(vdist,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,&num_partitions,&options,&edgecut,part, &mc) ;
#endif
    if(Loci::MPI_rank == 0)
      Loci::debugout << " Parmetis Edge cut   " <<  edgecut << endl ;
    delete [] xadj ;
    delete [] adjncy ;
    delete [] vdist ;

    //find the partition ptn given by Metis
    vector<entitySet> ptn ;

    for(int i = 0; i < Loci::MPI_processes; ++i)
      ptn.push_back(EMPTY) ;
    cmin = local_cells[Loci::MPI_rank].Min() ;
    for(int i=0;i<size_map;++i) {
      ptn[part[i]] += i + cmin ;
    }
    delete [] part ;
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
                 store<vector3d<real_t> > &t_pos, Map &tmp_cl,
                 Map &tmp_cr, multiMap &tmp_face2node,
                 entitySet nodes, entitySet faces, entitySet cells,
                 store<vector3d<real_t> > &pos, Map &cl, Map &cr,
                 multiMap &face2node) {

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
      Entity minc = min(cr[fc],cl[fc]) ;
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
    // Remember to add an update remap!!!

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

    entitySet loc_boundary_cells = getBoundaryCells(MapRepP(cr.Rep()));
    loc_boundary_cells = all_collect_entitySet(loc_boundary_cells) ;

    FORALL(loc_boundary_cells, li) {
      remap[li] = li;
    } ENDFORALL ;


    entitySet out_of_dom ;
    MapRepP f2n = MapRepP(face2node.Rep()) ;
    out_of_dom += cr.image(cr.domain())-(orig_cells+loc_boundary_cells) ;
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
			   const store<vector3d<real_t> > &t_pos,
			   const Map &tmp_cl, const Map &tmp_cr,
			   const multiMap &tmp_face2node,
			   store<vector3d<real_t> > &pos, Map &cl, Map &cr,
			   multiMap &face2node) {

    entitySet boundary_cells = getBoundaryCells(Loci::MapRepP(tmp_cr.Rep()));

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
    FORALL(boundary_cells, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;

    pos = t_pos.Rep()->remap(identity_map);
    cl = tmp_cl.Rep()->remap(identity_map);
    cr = tmp_cr.Rep()->remap(identity_map);
    face2node = MapRepP(tmp_face2node.Rep())->get_map();

  }

  vector<entitySet> partitionFaces(vector<entitySet> cell_ptn, const Map &cl,
                                   const Map &cr) {
    dstore<short> P ;
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
    dom -= interval(UNIVERSE_MIN,-1) ;
    dom -= cells ;
    {
      storeRepP PRep = P.Rep() ;
      fill_clone(PRep,dom,ptn_cells) ;
    }
    vector<entitySet> face_ptn(MPI_processes) ;
    entitySet boundary_faces ; // Boundary between processors
    FORALL(faces,fc) {
      if(cl[fc]<0)
        face_ptn[P[cr[fc]]] += fc ;
      else if(cr[fc]<0)
        face_ptn[P[cl[fc]]] += fc ;
      else if(P[cl[fc]] == P[cr[fc]])
        face_ptn[P[cl[fc]]] += fc ;
      else
        boundary_faces += fc ;
    } ENDFORALL ;


    // Number of balancing steps.  In the balancing step, the faces that
    // share proceessors are allocated by selected processors, then all
    // processors share current face counts.
    int STEPS = min(MPI_processes,7) ;
    for(int s=0;s<STEPS;++s) {

      vector<int> curr_sizes(MPI_processes),tot_sizes(MPI_processes) ;
      for(int i=0;i<MPI_processes;++i)
        curr_sizes[i] = face_ptn[i].size() ;

      MPI_Allreduce(&curr_sizes[0],&tot_sizes[0],MPI_processes,MPI_INT,MPI_SUM,
                    MPI_COMM_WORLD) ;

      if(s == 0 && MPI_rank ==0) {
        debugout << "fixed face sizes:" ;
        for(int i=0;i<MPI_processes;++i)
          debugout << ' ' << tot_sizes[i] ;
        debugout << endl ;
      }

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

    return face_ptn ;
  }

  vector<entitySet> partitionNodes(vector<entitySet> face_ptn, MapRepP face2node,entitySet old_node_dom) {
    dstore<int> np ;
    entitySet nall ;
    for(int i=0;i<MPI_processes;++i) {
      entitySet ntouch = face2node->image(face_ptn[i]) ;
      nall += ntouch ;
      FORALL(ntouch,nn) {
        np[nn] = i ;
      } ENDFORALL ;
    }
    vector<entitySet> node_ptn_old = all_collect_vectors(old_node_dom) ;
    vector<int> send_sz(MPI_processes) ;
    vector<entitySet> sendSets(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i) {
      if(i != MPI_rank)
        sendSets[i] = nall & node_ptn_old[i] ;
      send_sz[i] = sendSets[i].size()*2 ;
    }
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
    for(int i = 0; i <  MPI_processes; ++i) {
      int j = 0 ;
      FORALL(sendSets[i],ss) {
        send_store[send_displacement[i]+j*2] = ss ;
        send_store[send_displacement[i]+j*2+1] = np[ss] ;
        j++ ;
      } ENDFORALL ;
    }
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;

    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        nall += i1 ;
        np[i1] = i2 ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;


    vector<entitySet> node_ptn(MPI_processes) ;

    FATAL(((nall&old_node_dom)-old_node_dom) != EMPTY) ;

    FORALL(old_node_dom,nn) {

      node_ptn[np[nn]] += nn ;
    } ENDFORALL ;

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
      for(int j=0;j<ptn[i].num_intervals();++j) {
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


  //Description: Reads grid structures in the fact database
  //Input: facts and grid file name
  //Output: true if sucess
  bool readFVMGrid(fact_db &facts, string filename) {
    double t1 = MPI_Wtime() ;

    memSpace("readFVMGrid Start") ;
    vector<entitySet> local_nodes;
    vector<entitySet> local_cells;
    vector<entitySet> local_faces;

    store<vector3d<real_t> > t_pos;
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;

    int max_alloc = facts.get_max_alloc() ;

    bool useVOG = false ;

    string input_file = filename ;
    string::size_type spos = string::npos ;
    if((spos = input_file.rfind('.')) != string::npos) {
      input_file.erase(0,spos) ;
      if(input_file == ".vog")
        useVOG = true ;
    }

    vector<pair<int,string> > boundary_ids ;
    if(useVOG) {
      if(!readGridVOG(local_nodes, local_faces, local_cells,
                      t_pos, tmp_cl, tmp_cr, tmp_face2node,
                      max_alloc, filename, boundary_ids))
        return false;
    } else {
      if(!readGridXDR(local_nodes, local_faces, local_cells,
                      t_pos, tmp_cl, tmp_cr, tmp_face2node,
                      max_alloc, filename))
        return false;
    }

    memSpace("after reading grid") ;

    // Identify boundary tags
    entitySet local_boundary_cells = getBoundaryCells(MapRepP(tmp_cr.Rep()));

    entitySet global_boundary_cells = all_collect_entitySet(local_boundary_cells) ;
    if(MPI_processes == 1) {

      int npnts = local_nodes[0].size();
      int nfaces = local_faces[0].size();
      int ncells = local_cells[0].size();

      entitySet nodes = facts.get_distributed_alloc(npnts).first ;
      entitySet faces = facts.get_distributed_alloc(nfaces).first ;
      entitySet cells = facts.get_distributed_alloc(ncells).first;

      store<vector3d<real_t> > pos ;
      Map cl ;
      Map cr ;
      multiMap face2node ;
      copyGridStructures(nodes, faces, cells,
                         t_pos, tmp_cl, tmp_cr, tmp_face2node,
                         pos, cl, cr, face2node);


      store<string> boundary_names ;
      store<string> boundary_tags ;
      boundary_names.allocate(global_boundary_cells) ;
      boundary_tags.allocate(global_boundary_cells) ;
      Loci::debugout << " boundaries identified as:" ;

      
      
      FORALL(global_boundary_cells, bc) {
        char buf[512] ;
        sprintf(buf,"BC_%d",-bc) ;
        boundary_tags[bc] = string(buf) ;
        boundary_names[bc] = string(buf) ;
	debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;

      for(size_t i=0;i<boundary_ids.size();++i) {
        int id = boundary_ids[i].first ;
        if(global_boundary_cells.inSet(-id))
          boundary_names[-id] = boundary_ids[i].second ;
      }

      FORALL(global_boundary_cells, bc) {
	debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      Loci::debugout << endl ;

      facts.create_fact("cl", cl) ;
      facts.create_fact("cr", cr) ;
      facts.create_fact("pos", pos) ;
      facts.create_fact("face2node",face2node) ;
      facts.create_fact("boundary_names", boundary_names) ;
      facts.create_fact("boundary_tags", boundary_tags) ;
      return true ;
    }

    memSpace("before partitioning") ;
    // Partition Cells
    vector<entitySet> cell_ptn ;
    if(!use_simple_partition) {
      cell_ptn = newMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr) ;
    } else {
      cell_ptn = vector<entitySet>(MPI_processes) ;
      cell_ptn[MPI_rank] = local_cells[MPI_rank] ;
    }

    memSpace("mid partitioning") ;
    vector<entitySet> face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr) ;
    memSpace("after partitionFaces") ;
    vector<entitySet> node_ptn = partitionNodes(face_ptn,
                                                MapRepP(tmp_face2node.Rep()),
                                                t_pos.domain()) ;
    //    vector<entitySet> node_ptn(MPI_processes) ;
    //    node_ptn[MPI_rank] = local_nodes[MPI_rank] ;
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


    vector<std::pair<int, int> > boundary_update;
    FORALL(local_boundary_cells, li) {
      boundary_update.push_back(std::make_pair(li, li));
    }ENDFORALL;

    memSpace("before update_remap") ;
    facts.update_remap(boundary_update);
    memSpace("after update_remap") ;

    memSpace("before remapGridStructures") ;
    Map cl, cr ;
    multiMap face2node ;
    store<vector3d<real_t> > pos ;

    remapGrid(node_ptn, face_ptn, cell_ptn,
              node_ptn_t, face_ptn_t, cell_ptn_t,
              t_pos, tmp_cl, tmp_cr, tmp_face2node,
              nodes, faces, cells,
              pos, cl, cr, face2node);
    memSpace("after remapGridStructures") ;

    local_boundary_cells = getBoundaryCells(Loci::MapRepP(cr.Rep()));
    entitySet boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;

    store<string> boundary_names ;
    store<string> boundary_tags ;
    boundary_names.allocate(boundary_cells) ;
    boundary_tags.allocate(boundary_cells) ;
    

    FORALL(boundary_cells, bc) {
      char buf[512] ;
      sprintf(buf,"BC_%d",-bc) ;
      boundary_names[bc] = string(buf) ;
      boundary_tags[bc] = string(buf) ;
    } ENDFORALL ;

    for(size_t i=0;i<boundary_ids.size();++i) {
      int id = boundary_ids[i].first ;
      if(boundary_cells.inSet(-id))
        boundary_names[-id] = boundary_ids[i].second ;
    }

    if(Loci::MPI_rank == 0) {
      Loci::debugout << "boundaries identified as:" ;
      FORALL(boundary_cells, bc) {
        debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      Loci::debugout << endl ;
    }

    facts.create_fact("cl", cl) ;
    facts.create_fact("cr", cr) ;
    facts.create_fact("pos", pos) ;
    facts.create_fact("face2node",face2node) ;
    facts.create_fact("boundary_names", boundary_names) ;
    facts.create_fact("boundary_tags", boundary_names) ;

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

#ifdef DEBUG
    entitySet bdom2 = all_collect_entitySet(bdom,facts) ;
    FATAL(bdom2 != bdom) ;
#endif
    int ndom = bdom.size() ;
    int nloc = ndom/Loci::MPI_processes ;
    if((ndom%MPI_processes) > Loci::MPI_rank)
      nloc++ ;

    pair<entitySet,entitySet> alloc = facts.get_distributed_alloc(nloc) ;
    store<string> bn2 ;
    bn2.allocate(alloc.second) ;
    store<string> bt2 ;
    bt2.allocate(alloc.second) ;
    
    FATAL(bdom.size() != alloc.second.size()) ;
    Map mp ;
    mp.allocate(bdom) ;
    entitySet::const_iterator i1 = bdom.begin() ;
    entitySet::const_iterator i2 = alloc.second.begin() ;
    for(;i1!=bdom.end();++i1,++i2)
      mp[*i1] = *i2 ;

    for(i1=bdom.begin();i1!=bdom.end();++i1) {
      bn2[mp[*i1]] = boundary_names[*i1] ;
      bt2[mp[*i1]] = boundary_tags[*i1] ;
    }

    entitySet refdom = cr.preimage(bdom).first ;
    Map ref ;
    ref.allocate(refdom) ;

    for(i1=refdom.begin();i1!=refdom.end();++i1)
      ref[*i1] = mp[cr[*i1]] ;
    facts.create_fact("ref",ref) ;
    facts.update_fact("boundary_names",bn2) ;
    facts.update_fact("boundary_tags",bt2) ;
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
    entitySet bcset = interval(Loci::UNIVERSE_MIN,-1) ;
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
    store<vector3d<real_t> > pos ;
    pos = facts.get_variable("pos") ;
    store<vector3d<real_t> > fpos ;
    store<vector3d<real_t> > area ;

    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    constraint faces ;
    constraint interior_faces ;
    faces = facts.get_variable("faces") ;
    interior_faces = facts.get_variable("interior_faces") ;
    entitySet total_dom = Loci::MapRepP(face2node.Rep())->image(*faces) + pos.domain() ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    Loci::storeRepP pos_sp = pos.Rep() ;
    dstore<vector3d<real_t> > tmp_pos ;
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
      vector3d<real_t> fp(0,0,0) ;
      real_t w = 0 ;
      for(int i=0;i<nnodes;++i) {
        vector3d<real_t> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<real_t> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;

        real_t len = norm(p1-p2) ;

        fp += len*(p1+p2) ;
        w += len ;
      }
      fpos[fc] = fp/(2.*w) ;
      vector3d<real_t> a(0,0,0) ;
      for(int i=0;i<nnodes;++i) {
        vector3d<real_t> p1 = (tmp_pos[face2node[fc][i]]) ;
        vector3d<real_t> p2 = (tmp_pos[face2node[fc][(i+1)%nnodes]]) ;
        a += cross(p1-fpos[fc],p2-fpos[fc]) ;
      }
      area[fc] = .5*a ;
    } ENDFORALL ;
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    dstore<vector3d<real_t> > cpos ;
    dstore<real_t> cnum ;
    constraint geom_cells ;
    geom_cells = facts.get_variable("geom_cells") ;
    entitySet tmp_cells =  cl.image(*faces) | cr.image(*interior_faces) ;
    // Add cells owned by this processor!
    tmp_cells += (*geom_cells)& init_ptn[MPI_rank] ;
    cpos.allocate(tmp_cells) ;
    cnum.allocate(tmp_cells) ;
    FORALL(tmp_cells,cc) {
      cpos[cc] = vector3d<real_t>(0,0,0) ;
      cnum[cc] = 0 ;
    } ENDFORALL ;
    FORALL(*faces,fc) {
      real_t A = norm(area[fc]) ;
      cpos[cl[fc]] += A*fpos[fc] ;
      cnum[cl[fc]] += A ;
    } ENDFORALL ;
    FORALL(*interior_faces,fc) {
      real_t A = norm(area[fc]) ;
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
      dstore<vector3d<real_t> > tmp_cpos(v_cpos[i]) ;
      dstore<real_t> tmp_cnum(v_cnum[i]) ;
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
      vector3d<real_t> dv = cpos[cr[fc]]-cpos[cl[fc]] ;
      vector3d<real_t> dv2 = fpos[fc]-cpos[cl[fc]] ;
      vector3d<real_t> dv3 = cpos[cr[fc]]-fpos[fc] ;

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
      const vector3d<real_t> center = fpos[fc] ;

      vector3d<real_t> ccenter ;
      ccenter = cpos[cl[fc]] ;
      vector3d<real_t> dv = center-ccenter ;
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

    bool useVOG = false ;

    string input_file = filename ;
    string::size_type spos = string::npos ;
    if((spos = input_file.rfind('.')) != string::npos) {
      input_file.erase(0,spos) ;
      if(input_file == ".vog")
        useVOG = true ;
    }

    if(!useVOG) {
      memSpace("before color matrix") ;
      color_matrix(facts, COLOR_DFS) ;
      memSpace("before make_faces_consistent") ;
      make_faces_consistent(facts) ;
    } else {
      //      make_faces_consistent(facts) ;
    }
    return true ;
  }


}

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


namespace Loci {
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
    return true;
  }

  struct faceCluster {
    int nfaces ;
    int nnodes ;
    int ncells ;
    long nodeMap[256] ;
    long cellMap[256] ;
    int *face_sizes ;
    unsigned char *cluster_info ;
  } ;

  entitySet faceCluster(const multiMap &face2node, const Map &cl, const Map &cr, entitySet faces) {
    entitySet faceSet ;
    entitySet nodeSet ;
    entitySet cellSet ;
    entitySet fcluster ;
    
    int nnodes = 0 ;
    int ncells = 0 ;
    entitySet::const_iterator ei ;
    for(ei = faces.begin();ei!=faces.end();++ei) {
      int ncells_local = 0 ;
      int nnodes_local = 0 ;
      Entity fc = *ei ;
      if(!cellSet.inSet(cl[fc]))
        ncells_local++ ;
      if(!cellSet.inSet(cr[fc]))
        ncells_local++ ;
      int sz = face2node[fc].size() ;
      for(int i=0;i<sz;++i)
        if(!nodeSet.inSet(face2node[fc][i]))
          nnodes_local++ ;
      if(nnodes +nnodes_local >256 ||
         ncells +ncells_local > 256)
        break ;
      cellSet += cl[fc] ;
      cellSet += cr[fc] ;
      for(int i=0;i<sz;++i)
        nodeSet += face2node[fc][i] ;
      nnodes = nodeSet.size() ;
      ncells = cellSet.size() ;
      fcluster += fc ;
    }
    std::cout << "fcluster.size() = " << fcluster.size() << endl ;
    std::cout << "nodeSet = " << nodeSet.Max()-nodeSet.Min()+1 << endl ;
    //    std::cout << "cellSet = " << cellSet << endl ;
    
    return fcluster ;
  }
}
using namespace Loci ;
int main(int ac, char *av[]) {
    string filename = av[1] ;
    vector<entitySet> local_nodes;
    vector<entitySet> local_cells;
    vector<entitySet> local_faces;
    
    store<vector3d<real_t> > t_pos;
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;
    
    fact_db facts ;
    int max_alloc = facts.get_max_alloc() ;
    if(!readGridXDR(local_nodes, local_faces, local_cells,
		    t_pos, tmp_cl, tmp_cr, tmp_face2node,
		    max_alloc, filename))
      cerr << "unable to read grid file " << filename << endl ;

    entitySet faces = tmp_face2node.domain() ;
    while(faces != EMPTY) {
      entitySet fcluster = faceCluster(tmp_face2node,tmp_cl,tmp_cr,faces) ;
      faces -= fcluster ;
    }
    return 0 ;
  }



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
#include "loci_globs.h"

#include <Tools/tools.h>
#include <map>
#include <rpc/rpc.h>
#include <rpc/xdr.h>

#include <string>
using std::string ;
#include <vector>
using std::vector ;

extern "C" {
  typedef int idxtype ;
  void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
}

namespace Loci {

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
		   dstore<vector3d<double> > &pos, dMap &cl, dMap &cr,
		   dmultiMap &face2node, int max_alloc, string filename) {
    
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
      int data[3] ;
      MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD) ;
      npnts = data[0] ;
      nfaces = data[1] ;
      ncells = data[2] ;
    }
    
    // Create initial allocation of nodes, faces, and cells
    int node_ivl = npnts / Loci::MPI_processes;
    int face_ivl = nfaces / Loci::MPI_processes;
    int cell_ivl = ncells / Loci::MPI_processes;
    
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int nodes_base = max_alloc ;
      int faces_base = max_alloc+npnts ;
      int cells_base = max_alloc+npnts+nfaces ;
      if(i == Loci::MPI_processes-1) {
	local_nodes[i] = interval(nodes_base + i*node_ivl,
                                  nodes_base + npnts - 1) ;
	local_faces[i] = interval(faces_base + i*face_ivl,
                                  faces_base + nfaces-1) ;
	local_cells[i] = interval(cells_base + i*cell_ivl,
				  cells_base + ncells-1) ;
      }
      else {
	local_nodes[i] = interval(nodes_base + i*node_ivl,
                                  nodes_base + i*node_ivl + node_ivl - 1) ;
	local_faces[i] = interval(faces_base + i*face_ivl,
                                  faces_base + i*face_ivl + face_ivl - 1) ;
	local_cells[i] = interval(cells_base + i*cell_ivl,
                                  cells_base + i*cell_ivl + cell_ivl - 1) ;
      }
    }

    // Distribute positions
    if(Loci::MPI_rank == 0) {
      double *tmp_pos ;
      tmp_pos = new double[local_nodes[Loci::MPI_processes-1].size() * 3] ;
      for(int i = 0; i < 3*local_nodes[Loci::MPI_rank].size(); ++i)
      	if(!xdr_double(&xdr_handle, &tmp_pos[i]))
          return false ;
      
      int tmp = 0 ;
      for(entitySet::const_iterator ei = local_nodes[Loci::MPI_rank].begin(); ei != local_nodes[Loci::MPI_rank].end(); ++ei) {
	vector3d<double> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
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
      for(entitySet::const_iterator ei = local_nodes[Loci::MPI_rank].begin(); ei != local_nodes[Loci::MPI_rank].end(); ++ei) {
	vector3d<double> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
	tmp += 3 ;
	pos[*ei] = t ; 
      }
      
      delete [] tmp_pos ;
    }

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

  //Note: This function should not be called for serial version
  //Input: 
  // max_alloc(maximum allocated entity for fact_db that is 
  //                 the minimum entity for the returned distribution)
  // number of  nodes, faces and cells
  //Output: Returns naive distribution(init_ptn of fact_db) by allocating number of nodes, faces, cells
  vector<entitySet> getNaiveDistribution(int max_alloc,unsigned int num_nodes,
					 unsigned int num_faces,unsigned int num_cells) {
    fact_db local_facts ;
    local_facts.set_maximum_allocated(max_alloc) ;
    local_facts.get_distributed_alloc(num_nodes) ;
    local_facts.get_distributed_alloc(num_faces) ;
    local_facts.get_distributed_alloc(num_cells) ;
    return(local_facts.get_init_ptn()) ;
  }

  //Note: This function should not be called for serial version
  //Description: It does the preprocessing and creates datastrctures to call metisPartitionOfCells.
  //Input: 
  // tmp_cl, tmp_cr: mapping from faces to their left and right cells respectively
  // local_faces, local_cells: partition of faces and cells
  // init_ptn: partition of all entities
  //Output:
  // left_cells_to_cells, right_cells_to_cells: cell to cell mapping for cells on the left
  //      side and right side of the faces respectively
  // tmp_cl(modified), tmp_cr(modified): expanded for the faces which do not belong to 
  //       a processor but these faces has cells that belong to that processor on left
  //       and right side respectively
  void createCelltoCellMapping(dMap &tmp_cl, dMap &tmp_cr,  dmultiMap &left_cells_to_cells,
			       dmultiMap &right_cells_to_cells, const vector<entitySet> &local_faces,
			       const vector<entitySet> &local_cells, vector<entitySet> init_ptn) {
    entitySet cells, faces;
    for(int i = 0; i < Loci::MPI_processes; i++) {
      cells += local_cells[i];
      faces += local_faces[i];
    }

    entitySet local_interior_faces = getInteriorFaces(MapRepP(tmp_cr.Rep()));
    entitySet global_interior_faces = all_collect_entitySet(local_interior_faces) ;

    Loci::distributed_inverseMap(left_cells_to_cells, tmp_cl, cells, global_interior_faces, init_ptn) ;
    Loci::distributed_inverseMap(right_cells_to_cells, tmp_cr, cells, global_interior_faces, init_ptn) ;

    entitySet cl_out, cr_out ;
    
    //faces of which local cells are on left side
    entitySet cl_inv_ran = Loci::MapRepP(left_cells_to_cells.Rep())->image(local_cells[Loci::MPI_rank]) ;

    //faces of which local cells are on left side but faces don't belong to this processor.
    //Since all faces in left_cells_to_cells are interior faces, there should be another processor
    //which has cells that belong to right side of these faces
    cr_out  = cl_inv_ran - local_faces[Loci::MPI_rank] ;
    
    //Right now we have cells->faces and faces->cells mapping. We need to have cells->cells mapping.  
    //For that we need the information of the all faces for which local cells are either right or left sides
    Loci::storeRepP crsp = Loci::MapRepP(tmp_cr.Rep())->expand(cr_out, init_ptn) ;
    tmp_cr.setRep(crsp);
    
    //same logic as above to expand tmp_cl
    entitySet cr_inv_ran = Loci::MapRepP(right_cells_to_cells.Rep())->image(local_cells[Loci::MPI_rank]) ;
    cl_out  = cr_inv_ran - local_faces[Loci::MPI_rank] ; 
    Loci::storeRepP clsp = Loci::MapRepP(tmp_cl.Rep())->expand(cl_out, init_ptn) ;
    tmp_cl.setRep(clsp);
    
    //Now we can get cells->cells mapping for the metis patitioning of cells
    Loci::MapRepP(left_cells_to_cells.Rep())->compose(tmp_cr, left_cells_to_cells.domain()) ;
    Loci::MapRepP(right_cells_to_cells.Rep())->compose(tmp_cl, right_cells_to_cells.domain()) ;
  }
  
  //Note: This function should not be called for serial version
  //Input:
  // local_cells: current partition of cells
  // cell->cell mappings
  //Output:
  // metis partition of cells
  vector<entitySet> metisPartitionOfCells(const vector<entitySet> &local_cells, 
					  const dmultiMap &left_cells_to_cells, 
					  const dmultiMap &right_cells_to_cells) {

    //create parameters to call ParMetis for partitioning the cells
    dmultiMap dynamic_map ;
    FORALL(local_cells[Loci::MPI_rank], lci) {
      std::vector<int> tmp_int ;
      for(size_t i = 0; i < left_cells_to_cells[lci].size(); ++i)
	tmp_int.push_back(left_cells_to_cells[lci][i]) ;
      for(size_t i = 0; i < right_cells_to_cells[lci].size(); ++i)
	tmp_int.push_back(right_cells_to_cells[lci][i]) ;
      dynamic_map[lci] = tmp_int ;
    } ENDFORALL ;
    int count = 0 ;
    int size_map = local_cells[Loci::MPI_rank].size() ;
    entitySet dom_map = Loci::interval(0, size_map-1) ;
    store<int> size_adj ;
    size_adj.allocate(dom_map) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_adj[count] = dynamic_map[*ei].size() ;  
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
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) 
      for(size_t i = 0; i != dynamic_map[*ei].size(); ++i)        {
	adjncy[count] = dynamic_map[*ei][i] - cmin ;
	count ++ ;
      }
    vdist[0] = 0 ;
    for(int i = 1; i <= Loci::MPI_processes; ++i) 
      vdist[i] = vdist[i-1] + local_cells[i-1].size() ;
      
    MPI_Barrier(MPI_COMM_WORLD) ;
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
  
  //Note: This function should not be called for serial version
  //Input:
  // naive_init_ptn: initial naive partition of cells
  // metis_cell_ptn: partition of cells given by metis
  // global_nodes, global_faces, global_cells: entites on all processors for  nodes, faces, cells
  // t_pos: position vector of nodes
  // tmp_cl, tmp_cr: mapping from faces to cells on left and right side
  // tmp_face2node: mapping from faces to nodes 
  //Output: 
  // partition of all entities based on the cell partition given by metis
  // tmp_cl, tmp_cr, tmp_face2node, t_pos (all modified): Adds information of newly
  // added entities to my ownership
  vector<entitySet> newPartitionUsingCellPartition(vector<entitySet> naive_init_ptn, 
						   vector<entitySet> metis_cell_ptn,
						   entitySet global_nodes,
						   entitySet global_faces,
						   entitySet global_cells,
						   dstore<vector3d<double> > &t_pos, 
						   dMap &tmp_cl, dMap &tmp_cr, 
						   dmultiMap &tmp_face2node,
						   entitySet &naive_extra_comp_ent) {

    std::vector<entitySet> new_init_ptn = naive_init_ptn ;

    //remove old cells and add new cells to the partition
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      metis_cell_ptn[i] = all_collect_entitySet(metis_cell_ptn[i]) ;
      new_init_ptn[i] -= global_cells ;
      new_init_ptn[i] += metis_cell_ptn[i] ; 
    } 
    
    dmultiMap left_cells_to_faces ;
    Loci::distributed_inverseMap(left_cells_to_faces, tmp_cl, global_cells, global_faces, naive_init_ptn) ;

    //cells that are added to my ownership
    entitySet cells_out = metis_cell_ptn[Loci::MPI_rank] - left_cells_to_faces.domain();

    Loci::storeRepP inverse_sp = Loci::MapRepP(left_cells_to_faces.Rep())->expand(cells_out, naive_init_ptn) ;

    entitySet cl_inv_ran = Loci::MapRepP(inverse_sp)->image(metis_cell_ptn[Loci::MPI_rank]) ;
      
    //faces that are added to my ownership
    entitySet faces_out = cl_inv_ran - tmp_cl.domain() ;
      
    Loci::storeRepP clsp = Loci::MapRepP(tmp_cl.Rep())->expand(faces_out, naive_init_ptn) ;
    tmp_cl.setRep(clsp);
    
    std::vector<entitySet> v_req = all_collect_vectors(cl_inv_ran) ;
      
    //remove old faces and add new faces to the partition
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      new_init_ptn[i] -= global_faces ;
      new_init_ptn[i] += v_req[i] ;
    }

    if(duplicate_work) {
      dmultiMap right_cells_to_faces;
      Loci::distributed_inverseMap(right_cells_to_faces, tmp_cr, global_cells, global_faces, naive_init_ptn) ;
      inverse_sp = Loci::MapRepP(right_cells_to_faces.Rep())->expand(cells_out, naive_init_ptn) ;
      entitySet cr_inv_ran = Loci::MapRepP(inverse_sp)->image(metis_cell_ptn[Loci::MPI_rank]) ;
      naive_extra_comp_ent += cr_inv_ran - v_req[Loci::MPI_rank];
    }

    entitySet f2n_out = v_req[Loci::MPI_rank] - tmp_face2node.domain() ;
    Loci::storeRepP f2n_sp = Loci::MapRepP(tmp_face2node.Rep())->expand(f2n_out, naive_init_ptn) ; 
    tmp_face2node.setRep(f2n_sp);
    
    entitySet my_nodes = Loci::MapRepP(tmp_face2node.Rep())->image(v_req[Loci::MPI_rank]) ;
    v_req = all_collect_vectors(my_nodes) ;

    //In case of tie of node ownership, 
    //resolve it by giving node entity to the processor who has lower rank
    //remove old nodes and add new nodes to the partition
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      for(int j = i+1 ; j < Loci::MPI_processes; ++j) {
	entitySet  tmp = v_req[i] & v_req[j] ;
	v_req[j] -= tmp ;
      } 
      new_init_ptn[i] -= global_nodes ;
      new_init_ptn[i] += v_req[i] ;
    }
    
    //add information of newly added faces to tmp_cr and nodes to t_pos
    entitySet cr_out = cl_inv_ran - tmp_cr.domain() ;
    Loci::storeRepP crsp = Loci::MapRepP(tmp_cr.Rep())->expand(cr_out, naive_init_ptn) ;
    tmp_cr.setRep(crsp);
    
    entitySet pos_out = my_nodes - t_pos.domain() ;
    Loci::storeRepP t_pos_sp = t_pos.Rep() ;
    entitySet total_dom = t_pos.domain() + my_nodes ;
    fill_clone(t_pos_sp, pos_out, naive_init_ptn) ;
      
    return new_init_ptn;
  }

  //Note: This function should not be called for serial version.  copyGridStructure is used for that.
  //Input:
  // loc_nodes, loc_faces, loc_cells: old entities of nodes, cells and faces to my ownership
  // tmp_cl, tmp_cr: old mapping from face to cell on left and right side
  // tmp_face2node: old mapping from face to nodes
  // t_pos: position of old nodes

  // nodes, faces, cells: new entities of nodes, faces and cells to my ownership
  //Output:
  // cl, cl: new mapping from face to cell
  // face2node: new mapping from face to nodes
  // pos: position reamapped to new nodes
  void remapGridStructures(entitySet loc_nodes, entitySet loc_faces, entitySet loc_cells, 
			   const dstore<vector3d<double> > &t_pos, const dMap &tmp_cl,
			   const dMap &tmp_cr, const dmultiMap &tmp_face2node, 
			   entitySet nodes, entitySet faces, entitySet cells,
		   	   store<vector3d<double> > &pos, Map &cl, Map &cr, multiMap &face2node) {

    dMap remap;
    entitySet::const_iterator ei = loc_nodes.begin() ;
    FORALL(nodes, li) {
      remap[*ei] = li ;
      ++ei ;
    } ENDFORALL ;
      
    pos = t_pos.Rep()->remap(remap);

    ei = loc_faces.begin() ;
    FORALL(faces, li) {
      remap[*ei] = li ;
      ++ei ;
    } ENDFORALL ;
      
    ei = loc_cells.begin() ;
    FORALL(cells, li) {
      remap[*ei] = li ;
      ++ei ;
    } ENDFORALL ;

    entitySet loc_boundary_cells = getBoundaryCells(Loci::MapRepP(tmp_cr.Rep()));
    FORALL(loc_boundary_cells, li) {
      remap[li] = li;
    } ENDFORALL ;


    std::vector<entitySet> node_ptn = all_collect_vectors(loc_nodes) ;
    std::vector<entitySet> cell_ptn = all_collect_vectors(loc_cells) ;

    //Next few lines adds information of entities of nodes and cells belong to other processors  to the remap
    //That information is needed to call remap() on tmp_face2node, tmp_cl and tmp_cr
    entitySet entities_accessed = Loci::MapRepP(tmp_face2node.Rep())->image(loc_faces) ;
    entitySet remap_out = entities_accessed - loc_nodes ;
    Loci::storeRepP remap_sp = Loci::MapRepP(remap.Rep())->expand(remap_out, node_ptn) ;
    remap.setRep(remap_sp) ; 

    entities_accessed = Loci::MapRepP(tmp_cl.Rep())->image(loc_faces) ;
    remap_out = entities_accessed - loc_cells ;
    entities_accessed = Loci::MapRepP(tmp_cr.Rep())->image(loc_faces) ;
    entities_accessed &= interval(0, Loci::UNIVERSE_MAX) ;
    remap_out += entities_accessed - loc_cells ;
    remap_sp = Loci::MapRepP(remap.Rep())->expand(remap_out, cell_ptn) ;
    remap.setRep(remap_sp) ; 

    cl = tmp_cl.Rep()->remap(remap);
    cr = tmp_cr.Rep()->remap(remap);
    face2node = tmp_face2node.Rep()->remap(remap);
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
			   const dstore<vector3d<double> > &t_pos,
			   const dMap &tmp_cl, const dMap &tmp_cr,
			   const dmultiMap &tmp_face2node,
			   store<vector3d<double> > &pos, Map &cl, Map &cr,
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

  //Description: Reads grid structures in the fact database
  //Input: facts and grid file name
  //Output: true if sucess 
  bool readFVMGrid(fact_db &facts, string filename) {
    
    vector<entitySet> local_nodes;
    vector<entitySet> local_cells;
    vector<entitySet> local_faces;
    
    dstore<vector3d<double> > t_pos;
    dMap tmp_cl, tmp_cr;
    dmultiMap tmp_face2node;
    
    int max_alloc = facts.get_max_alloc();

    if(!readGridXDR(local_nodes, local_faces, local_cells,
		    t_pos, tmp_cl, tmp_cr, tmp_face2node,
		    max_alloc, filename))
      return false;

    int npnts = 0, nfaces = 0, ncells = 0;
    for(int i = 0; i < Loci::MPI_processes; i++) {
      npnts += local_nodes[i].size();
      nfaces += local_faces[i].size();
      ncells += local_cells[i].size();
    }

    entitySet nodes, faces, cells; 
    entitySet naive_loc_nodes, naive_loc_faces, naive_loc_cells ;
    
    if(Loci::MPI_processes > 1) {
      entitySet global_nodes, global_faces, global_cells;
      for(int i = 0; i < Loci::MPI_processes; i++) {
	global_nodes += local_nodes[i];
	global_faces += local_faces[i];
	global_cells += local_cells[i];
      }
      entitySet local_boundary_cells = getBoundaryCells(Loci::MapRepP(tmp_cr.Rep()));
      
      entitySet global_boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;
      store<string> boundary_names ;
      boundary_names.allocate(global_boundary_cells) ;
      if(Loci::MPI_rank == 0) 
	Loci::debugout << "boundaries identified as:" ;
      FORALL(global_boundary_cells, bc) {
	char buf[512] ;
	sprintf(buf,"BC_%d",-bc) ;
	boundary_names[bc] = string(buf) ;
	if(Loci::MPI_rank == 0)
	  Loci::debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      if(Loci::MPI_rank == 0)
	Loci::debugout << endl ;

      Loci::debugout << "original boundaries: " << global_boundary_cells << endl;
      std::vector<entitySet> vset = Loci::all_collect_vectors(local_boundary_cells) ;

      std::vector<entitySet> naive_init_ptn = getNaiveDistribution(max_alloc, local_nodes[Loci::MPI_rank].size(), local_faces[Loci::MPI_rank].size(), local_cells[Loci::MPI_rank].size());
      
      for(int i = 0; i < Loci::MPI_processes; ++i) 
	naive_init_ptn[i] += vset[i] ; 

      dmultiMap left_cells_to_cells, right_cells_to_cells;
      createCelltoCellMapping(tmp_cl, tmp_cr, left_cells_to_cells, right_cells_to_cells, 
			      local_faces, local_cells, naive_init_ptn) ;

      vector<entitySet> metis_cell_ptn = metisPartitionOfCells(local_cells, left_cells_to_cells,
							       right_cells_to_cells);

      entitySet naive_extra_comp_ent;
      vector<entitySet> new_init_ptn = newPartitionUsingCellPartition(naive_init_ptn, metis_cell_ptn, 
								      global_nodes, global_faces, global_cells,
								      t_pos, tmp_cl, tmp_cr, tmp_face2node,
								      naive_extra_comp_ent);

      naive_loc_nodes = global_nodes & new_init_ptn[Loci::MPI_rank] ;
      naive_loc_faces = global_faces & new_init_ptn[Loci::MPI_rank] ;
      naive_loc_cells = global_cells & new_init_ptn[Loci::MPI_rank] ;

      vector<int> alloc_vector;
      FORALL(naive_loc_nodes, ni) {
	alloc_vector.push_back(ni);
      }ENDFORALL;

      nodes = facts.get_distributed_alloc(alloc_vector).first ;

      alloc_vector.resize(0);
      FORALL(naive_loc_faces, ni) {
	alloc_vector.push_back(ni);
      }ENDFORALL;

      faces = facts.get_distributed_alloc(alloc_vector).first ;

      alloc_vector.resize(0);
      FORALL(naive_loc_cells, ni) {
	alloc_vector.push_back(ni);
      }ENDFORALL;

      cells = facts.get_distributed_alloc(alloc_vector).first ;

      vector<std::pair<int, int> > boundary_update;
      FORALL(local_boundary_cells, li) {
	boundary_update.push_back(std::make_pair(li, li));
      }ENDFORALL;

      facts.update_remap(boundary_update);
      if(duplicate_work) {
	fact_db::distribute_infoP df = facts.get_distribute_info();
	dMap remap;
	Loci::storeRepP remap_sp = df->remap.Rep();
	Loci::MapRepP(remap)->copy(remap_sp, remap_sp->domain());
	entitySet comp_out = naive_extra_comp_ent - remap.domain();
	remap_sp = Loci::MapRepP(remap.Rep())->expand(comp_out, new_init_ptn);
	//facts.global_comp_entities += Loci::MapRepP(remap_sp)->image(naive_extra_comp_ent);
      }
    }
    else {
      nodes = facts.get_allocation(npnts) ; 
      faces = facts.get_allocation(nfaces) ;
      cells = facts.get_allocation(ncells);
    }

    Loci::debugout << "nodes = " << nodes << endl;
    Loci::debugout << "faces = " <<faces << endl ;
    Loci::debugout << "cells = " << cells << endl ;
    
    Map cl, cr ;
    multiMap face2node ;
    store<vector3d<double> > pos ;
    if(Loci::MPI_processes > 1)
      remapGridStructures(naive_loc_nodes, naive_loc_faces, naive_loc_cells,
			  t_pos, tmp_cl, tmp_cr, tmp_face2node,
			  nodes, faces, cells,
			  pos, cl, cr, face2node);
    else 
      copyGridStructures(nodes, faces, cells, t_pos, tmp_cl, tmp_cr, tmp_face2node, pos, cl, cr, face2node);
    
    entitySet boundary_cells, local_boundary_cells;

    if(Loci::MPI_processes > 1) {
      local_boundary_cells = getBoundaryCells(Loci::MapRepP(cr.Rep()));
      boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;
    }
    else
      boundary_cells = getBoundaryCells(Loci::MapRepP(cr.Rep()));
    store<string> boundary_names ;
    boundary_names.allocate(boundary_cells) ;
    if(Loci::MPI_rank == 0) {
      if(Loci::MPI_processes > 1) 
	Loci::debugout << " new" ;
      Loci::debugout << " boundaries identified as:" ;
    }      

    FORALL(boundary_cells, bc) {
      char buf[512] ;
      sprintf(buf,"BC_%d",-bc) ;
      boundary_names[bc] = string(buf) ;
      if(Loci::MPI_rank == 0 )
	Loci::debugout << " " << boundary_names[bc] ;
    } ENDFORALL ;
    
    if(Loci::MPI_rank == 0)
      Loci::debugout << endl ;

    if(Loci::MPI_processes > 1) {
      vector<entitySet> vset = all_collect_vectors(local_boundary_cells) ;
      for(int i = 0; i < Loci::MPI_processes; ++i) {
	vector<entitySet> tmp_init_ptn = facts.get_init_ptn() ;
	tmp_init_ptn[i] += vset[i] ;
	Loci::debugout << " init_ptn[" <<i << "] = " << tmp_init_ptn[i] << endl ;
	facts.put_init_ptn(tmp_init_ptn);
      }
    }

    param<int> min_node ;

    if(Loci::MPI_processes > 1)
      *min_node = max_alloc ;
    else 
      *min_node = 0;
    
    facts.create_fact("min_node", min_node) ; 
    facts.create_fact("cl", cl) ;
    facts.create_fact("cr", cr) ;
    facts.create_fact("pos", pos) ;
    facts.create_fact("face2node",face2node) ;
    facts.create_fact("boundary_names", boundary_names) ;
    
    return true ;
  }
}

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
#include <map>

#include <Tools/xdr.h>
#include <sstream>
using std::istringstream ;
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

//using Loci::fact_db ;
using Loci::debugout ;
namespace Loci {
  extern vector<entitySet> newMetisPartitionOfCells(const vector<entitySet> &local_cells,
                                             const Map &cl, const Map &cr) ;
  extern void distributed_inverseMap(multiMap &result,
                                     vector<pair<Entity,Entity> > &input,
                                     entitySet input_image,
                                     entitySet input_preimage,
                                     const std::vector<entitySet> &init_ptn) ;
}
namespace VOG {
  using Loci::MapRepP ;
  using Loci::storeRepP ;
  using Loci::MPI_processes ;
  using Loci::MPI_rank ;
  using Loci::create_intervalSet ;
  using Loci::UNIVERSE_MIN ;
  using Loci::UNIVERSE_MAX ;

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
		   store<vector3d<double> > &pos, Map &cl, Map &cr,
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
      pos.allocate(local_nodes[Loci::MPI_rank]) ;
      for(entitySet::const_iterator ei = local_nodes[Loci::MPI_rank].begin(); ei != local_nodes[Loci::MPI_rank].end(); ++ei) {
	vector3d<double> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
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
	if(cl[*ei] < 0) {
	  if(cr[*ei] < 0) {
	    cerr << " boundary condition on both sides of a face?" << endl ;
	    exit(1) ;
	  } else {
	    int tmp_swap = cr[*ei] ;
	    cr[*ei] = cl[*ei] ;
	    cl[*ei] = tmp_swap ;
	  }
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
	if(cl[*ei] < 0) {
	  if(cr[*ei] < 0) {
	    cerr << "2 boundary condition on both sides of a face?" << endl ;
	    exit(1) ;
	  } else {
	    int tmp = cr[*ei] ;
	    cr[*ei] = cl[*ei] ;
	    cl[*ei] = tmp ;
	  }
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
                 entitySet nodes, entitySet faces, entitySet cells,
                 store<vector3d<double> > &pos, Map &cl, Map &cr,
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
      Entity maxc = max(cr[fc],cl[fc]) ;
      cr[fc] = minc ;
      cl[fc] = maxc ;
      sortlist[i++] = pair<Entity,Entity>(minc,fc) ;
    } ENDFORALL ;
    std::sort(sortlist.begin(),sortlist.end(),fieldSort) ;
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
    loc_boundary_cells = Loci::all_collect_entitySet(loc_boundary_cells) ;

    FORALL(loc_boundary_cells, li) {
      remap[li] = li;
    } ENDFORALL ;

    
    entitySet out_of_dom ;
    MapRepP f2n = MapRepP(face2node.Rep()) ;
    out_of_dom += cr.image(cr.domain())-(orig_cells+loc_boundary_cells) ;
    out_of_dom += cl.image(cl.domain())-orig_cells ;
    out_of_dom += f2n->image(f2n->domain())-old_nodes ;
    entitySet old_dom = orig_cells+old_nodes ;
    vector<entitySet> old_ptn = Loci::all_collect_vectors(old_dom) ;
    {
      storeRepP PRep = remap.Rep() ;
      Loci::fill_clone(PRep,out_of_dom,old_ptn) ;
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
			   const store<vector3d<double> > &t_pos,
			   const Map &tmp_cl, const Map &tmp_cr,
			   const multiMap &tmp_face2node,
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
    vector<entitySet> ptn_cells = Loci::all_collect_vectors(cells) ;
    entitySet faces = cl.domain() & cr.domain() ;
    entitySet dom = cl.image(faces) | cr.image(faces) ;
    dom -= interval(UNIVERSE_MIN,-1) ;
    dom -= cells ;
    {
      storeRepP PRep = P.Rep() ;
      Loci::fill_clone(PRep,dom,ptn_cells) ;
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
    vector<entitySet> node_ptn_old = Loci::all_collect_vectors(old_node_dom) ;
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
  bool readGridFile(fact_db &facts, string filename) {
    double t1 = MPI_Wtime() ;

    vector<entitySet> local_nodes;
    vector<entitySet> local_cells;
    vector<entitySet> local_faces;
    
    store<vector3d<double> > t_pos;
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;
    
    int max_alloc = facts.get_max_alloc() ;

    if(!VOG::readGridXDR(local_nodes, local_faces, local_cells,
		    t_pos, tmp_cl, tmp_cr, tmp_face2node,
		    max_alloc, filename))
      return false;


    // Identify boundary tags
    entitySet local_boundary_cells = getBoundaryCells(MapRepP(tmp_cr.Rep()));
      
    entitySet global_boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;
    if(MPI_processes == 1) {

      int npnts = local_nodes[0].size();
      int nfaces = local_faces[0].size();
      int ncells = local_cells[0].size();

      entitySet nodes = facts.get_distributed_alloc(npnts).first ; 
      entitySet faces = facts.get_distributed_alloc(nfaces).first ;
      entitySet cells = facts.get_distributed_alloc(ncells).first;

      store<vector3d<double> > pos ;
      Map cl ;
      Map cr ;
      multiMap face2node ;
      copyGridStructures(nodes, faces, cells,
                         t_pos, tmp_cl, tmp_cr, tmp_face2node,
                         pos, cl, cr, face2node);
      
      store<string> boundary_names ;
      boundary_names.allocate(global_boundary_cells) ;
      Loci::debugout << " boundaries identified as:" ;
      FORALL(global_boundary_cells, bc) {
        char buf[512] ;
	bzero(buf,512) ;
        snprintf(buf,511,"BC_%d",-bc) ;
        boundary_names[bc] = string(buf) ;
	debugout << " " << boundary_names[bc] ;
      } ENDFORALL ;
    
      Loci::debugout << endl ;

      facts.create_fact("cl", cl) ;
      facts.create_fact("cr", cr) ;
      facts.create_fact("pos", pos) ;
      facts.create_fact("face2node",face2node) ;
      facts.create_fact("boundary_names", boundary_names) ;
      return true ;
        
    
    }


#define METIS_PARTITIONING
#ifdef METIS_PARTITIONING
    vector<entitySet> cell_ptn ;
    cell_ptn = Loci::newMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr) ;
    vector<entitySet> face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr) ;
    vector<entitySet> node_ptn = partitionNodes(face_ptn,
                                                MapRepP(tmp_face2node.Rep()),
                                                t_pos.domain()) ;
#else
    vector<entitySet> cell_ptn = vector<entitySet>(MPI_processes) ;
    cell_ptn[MPI_rank] = local_cells[MPI_rank] ;

    vector<entitySet> face_ptn(MPI_processes) ;
    face_ptn[MPI_rank] = local_faces[MPI_rank] ;
    
    vector<entitySet> node_ptn(MPI_processes) ;
    node_ptn[MPI_rank] = local_nodes[MPI_rank] ;
#endif
    
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

    facts.update_remap(boundary_update);

    Map cl, cr ;
    multiMap face2node ;
    store<vector3d<double> > pos ;

    remapGrid(node_ptn, face_ptn, cell_ptn,
              node_ptn_t, face_ptn_t, cell_ptn_t,
              t_pos, tmp_cl, tmp_cr, tmp_face2node,
              nodes, faces, cells,
              pos, cl, cr, face2node);

    local_boundary_cells = getBoundaryCells(Loci::MapRepP(cr.Rep()));
    entitySet boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;

    store<string> boundary_names ;

    boundary_names.allocate(boundary_cells) ;
    if(Loci::MPI_rank == 0) {
      Loci::debugout << "boundaries identified as:" ;
    }      

    FORALL(boundary_cells, bc) {
      char buf[512] ;
      bzero(buf,512) ;
      snprintf(buf,511,"BC_%d",-bc) ;
      boundary_names[bc] = string(buf) ;
      if(Loci::MPI_rank == 0 )
	Loci::debugout << " " << boundary_names[bc] ;
    } ENDFORALL ;
    
    if(Loci::MPI_rank == 0)
      Loci::debugout << endl ;

    facts.create_fact("cl", cl) ;
    facts.create_fact("cr", cr) ;
    facts.create_fact("pos", pos) ;
    facts.create_fact("face2node",face2node) ;
    facts.create_fact("boundary_names", boundary_names) ;

    double t2 = MPI_Wtime() ;
    debugout << "Time to read in file '" << filename << ", is " << t2-t1
             << endl ;
    return true ;
  }

  void getGridInfo(fact_db &facts, entitySet &faces,
                   entitySet &interior_faces, entitySet &boundary_faces,
                   entitySet &geom_cells) {
    Map cl, cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    faces = (cl.domain() & cr.domain()) ;
    entitySet bcset = interval(Loci::UNIVERSE_MIN,-1) ;
    boundary_faces = cr.preimage(bcset).first ;

    interior_faces = (faces-boundary_faces) ;

    vector<int> vec ;
    FORALL(interior_faces,fc) {
      vec.push_back(cl[fc]) ;
      vec.push_back(cr[fc]) ;
    } ENDFORALL ;
    FORALL(boundary_faces,fc) {
      vec.push_back(cl[fc]) ;
    } ENDFORALL ;

    geom_cells = Loci::create_entitySet(vec.begin(),vec.end()) ;
    entitySet global_geom = Loci::all_collect_entitySet(geom_cells,facts) ;
    std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
    geom_cells = global_geom&init_ptn[MPI_rank] ; 
  }
  


  bool inputXDRGrid(fact_db &facts, string filename) {

    // Read in the xdr grid
    if(MPI_rank == 0)
      cerr << "reading " << filename << endl ;
    if(!readGridFile(facts,filename))
      return false ;

    return true ;
  }
}

int main(int ac, char *av[]) {
  using namespace Loci ;
  using namespace VOG ;

  bool optimize = true ;
  
  Loci::Init(&ac,&av) ;

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
    } else {
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


  if(ac != 2) {
    cerr << "Usage: xdr2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -m : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;
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

  if(Loci::MPI_rank == 0) {
    cout << "input grid file units = " << tp ;
    if(posScale != 1.0) 
      cout << " = " << posScale << " meters " ;
    cout << endl ;
  }

  string filename = av[1] ;
  filename += ".xdr" ;

  string outfile = av[1] ;
  outfile += ".vog" ;

  store<vector3d<double> > pos;
  Map cl, cr;
  multiMap face2node;

  // Read in Mesh
  {
    fact_db facts ;
    if(!inputXDRGrid(facts,filename)) {
      cerr << "unable to read grid file " << filename << endl ;
      Loci::Abort() ;
    }

    pos = facts.get_variable("pos")->getRep() ;
    cl = facts.get_variable("cl")->getRep() ;
    cr = facts.get_variable("cr")->getRep() ;
    face2node = facts.get_variable("face2node")->getRep() ;

  }

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }

  // establish face left-right orientation
  if(MPI_rank == 0) 
    cerr << "orienting faces" << endl ;
  VOG::orientFaces(pos,cl,cr,face2node) ;
    
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0) 
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }
  
  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;

  vector<BC_descriptor> bcs ;

  if(MPI_rank == 0) {
    string tagsfile = string(av[1]) + ".tags" ;
    bcs = readTags(tagsfile) ;
    if(bcs.size() != 0) {
      cout << "boundary faces:"<< endl ;
      for(size_t i=0;i<bcs.size();++i)
        cout << bcs[i].name << ' ' << bcs[i].id << endl ;
    }
  }

  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;
  
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  return 0 ;
}
  

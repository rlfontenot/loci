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
#include <hdf5.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include "sciTypes.h"
#include "defines.h"
#include <LociGridReaders.h>
#include <Tools/tools.h>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;
using std::vector;
using Loci::MPI_rank;
using Loci::MPI_processes;
namespace Loci {

  extern void memSpace(string s) ;

  extern bool use_simple_partition ;
  extern bool use_orb_partition ;
  void ORB_Partition_Mesh(const vector<entitySet> &local_nodes,
                          const vector<entitySet> &local_faces,
                          const vector<entitySet> &local_cells,
                          const store<vector3d<real_t> > &pos,
                          const Map &cl, const Map &cr,
                          const multiMap &face2node,
                          vector<entitySet> &cell_ptn,
                          vector<entitySet> &face_ptn,
                          vector<entitySet> &node_ptn);
  entitySet getBoundaryCells(const Loci::MapRepP tmp_cr_rep);
  void copyGridStructures( entitySet nodes, entitySet faces, entitySet cells,
			   const store<vector3d<real_t> > &t_pos,
			   const Map &tmp_cl, const Map &tmp_cr,
			   const multiMap &tmp_face2node,
			   store<vector3d<real_t> > &pos, Map &cl, Map &cr,
			   multiMap &face2node);
  vector<entitySet> newMetisPartitionOfCells(const vector<entitySet> &local_cells,
                                             const Map &cl, const Map &cr);
  vector<entitySet> partitionFaces(vector<entitySet> cell_ptn, const Map &cl,
                                   const Map &cr);
  vector<entitySet> transposePtn(const vector<entitySet> &ptn);
  vector<entitySet> partitionNodes(vector<entitySet> face_ptn, MapRepP face2node,entitySet old_node_dom);
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
                 multiMap &face2node,
                 fact_db &facts);
 
  void create_face_info(fact_db &facts);
  void create_ref(fact_db &facts) ;
  void create_ghost_cells(fact_db &facts);
  
  std::vector<int> simplePartitionVec(int mn, int mx, int p);
  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm);
  vector<sequence> transposeSeq(const vector<sequence> sv);
}
void colorMatrix(Map &cl, Map &cr, multiMap &face2node) ;

namespace Loci{
  
  // Convert container from global numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator(comm)
  storeRepP Global2FileOrder(storeRepP sp, entitySet dom, int &offset,
                             fact_db::distribute_infoP dist, MPI_Comm comm) {
   
    // Now get global to file numbering
    dMap g2f ;
    g2f = dist->g2f.Rep() ;

    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[i] ;
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    // Find bounds in file numbering from this processor
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;
    
    // Find overall bounds
    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Get partitioning of file numbers across processors
    vector<entitySet> out_ptn = simplePartition(imn,imx,comm) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;

    // Loop over processors and compute sets of entities to send
    // To efficiently compute this mapping, first sort the transpose
    // of the newnum map to quickly find the set of entities to send
    // without searching entire newnum map for each processor
    vector<pair<int,int> > file2num(dom.size()) ;
    size_t cnt = 0 ;
    FORALL(dom,ii) {
      file2num[cnt].first = newnum[ii] ;
      file2num[cnt].second = ii ;
      cnt++ ;
    } ENDFORALL ;
    sort(file2num.begin(),file2num.end()) ;

    // Check each processor, find out which sets to send
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = out_ptn[i].Max() ;
      while(cnt < file2num.size() && file2num[cnt].first <= mxi) {
        send_sets[i] += file2num[cnt].second ;
        cnt++ ;
      }
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;


    // shift by the offset
    offset = out_ptn[prank].Min() ;
    for(int i=0;i<p;++i)
      recv_seqs[i] <<= offset ;

    // Compute allocation domain
    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    // allocate store over shifted domain
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;

    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
		  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    return qcol_rep ;
  } 


  
  std::vector<entitySet> getDist( Loci::entitySet &faces,
                                  Loci::entitySet &cells,
                                  Map &cl, Map &cr, multiMap &face2node) {
  
    // First establish current distribution of entities across processors
    std::vector<Loci::entitySet> ptn(Loci::MPI_processes) ; // entity Partition

    // Get entity distributions
  
    faces = face2node.domain() ;
    entitySet allFaces = Loci::all_collect_entitySet(faces) ;
    std::vector<int> facesizes(MPI_processes) ;
    int  size = faces.size() ;
    MPI_Allgather(&size,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int  cnt = allFaces.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] += interval(cnt,cnt+facesizes[i]-1) ;
      cnt += facesizes[i] ;
    }
    
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
    int mn = geom_cells.Min() ;
    int mx = geom_cells.Max() ;
    std:: vector<int> pl = Loci::simplePartitionVec(mn,mx,MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      ptn[i] += interval(pl[i],pl[i+1]-1) ;
    faces = allFaces ;
    cells = geom_cells ;
    return ptn ;
  }

  bool inputFVMGrid(fact_db &facts,
                    vector<entitySet>& local_nodes,
                    vector<entitySet>& local_faces,
                    vector<entitySet>& local_cells,
                    store<vector3d<double> >& t_pos,
                    Map& tmp_cl,
                    Map& tmp_cr,
                    multiMap& tmp_face2node,
                    vector<pair<int,string> >& boundary_ids,
                    vector<pair<string,entitySet> >& volTags
                    ) {
    double t1 = MPI_Wtime() ;
    // Identify boundary tags
    entitySet local_boundary_cells = getBoundaryCells(MapRepP(tmp_cr.Rep()));
    entitySet global_boundary_cells = all_collect_entitySet(local_boundary_cells) ;
       
    if(Loci::MPI_processes == 1) {

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
      Loci::copyGridStructures(nodes, faces, cells,
                               t_pos, tmp_cl, tmp_cr, tmp_face2node,
                               pos, cl, cr, face2node);
      store<string> boundary_names ;
      store<string> boundary_tags ;
      boundary_names.allocate(global_boundary_cells) ;
      boundary_tags.allocate(global_boundary_cells) ;
      cout << " boundaries identified as:" ;

      
      
      FORALL(global_boundary_cells, bc) {
        char buf[512] ;
	bzero(buf,512) ;
        snprintf(buf,511,"BC_%d",-bc) ;
        boundary_tags[bc] = string(buf) ;
        boundary_names[bc] = string(buf) ;
	cout << " " << boundary_names[bc] ;
      } ENDFORALL ;

      for(size_t i=0;i<boundary_ids.size();++i) {
        int id = boundary_ids[i].first ;
        if(global_boundary_cells.inSet(-id))
          boundary_names[-id] = boundary_ids[i].second ;
      }

      FORALL(global_boundary_cells, bc) {
	cout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      cout << endl ;
      

     

      
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
        std::ostringstream oss ;
        oss << "volumeTag(" << volTags[i].first << ")" ;
        facts.create_fact(oss.str(),Tag) ;
      }
        
      return true ;
    }

    Loci::memSpace("before partitioning") ;

    vector<entitySet> cell_ptn,face_ptn,node_ptn ;

    if(use_orb_partition) {
      ORB_Partition_Mesh(local_nodes, local_faces, local_cells,
                         t_pos, tmp_cl, tmp_cr, tmp_face2node,
                         cell_ptn,face_ptn,node_ptn) ;
    } else {
      // Partition Cells
      if(!use_simple_partition) {
        cell_ptn = newMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr) ;
      } else {
        cell_ptn = vector<entitySet>(MPI_processes) ;
        cell_ptn[MPI_rank] = local_cells[MPI_rank] ;
      }

      Loci::memSpace("mid partitioning") ;
      face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr) ;
      Loci::memSpace("after partitionFaces") ;

      node_ptn = partitionNodes(face_ptn,
                                MapRepP(tmp_face2node.Rep()),
                                t_pos.domain()) ;
    }
    Loci::memSpace("after partitioning") ;
     
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

    debugout << "nodes = " << nodes << ", size= "
             << nodes.size() << endl;
    debugout << "faces = " << faces << ", size = "
             << faces.size() << endl ;
    debugout << "cells = " << cells << ", size = "
             << cells.size() << endl ;


    vector<std::pair<int, int> > boundary_update;
    FORALL(local_boundary_cells, li) {
      boundary_update.push_back(std::make_pair(li, li));
    }ENDFORALL;
    
   
    Loci::memSpace("before update_remap") ;
    facts.update_remap(boundary_update);
    Loci::memSpace("after update_remap") ;
    
   
    Loci::memSpace("before remapGridStructures") ;
    Map cl, cr ;
    multiMap face2node ;
    store<vector3d<real_t> > pos ;

    remapGrid(node_ptn, face_ptn, cell_ptn,
              node_ptn_t, face_ptn_t, cell_ptn_t,
              t_pos, tmp_cl, tmp_cr, tmp_face2node,
              nodes, faces, cells,
              pos, cl, cr, face2node,facts);
    Loci::memSpace("after remapGridStructures") ;
         



    local_boundary_cells = getBoundaryCells(Loci::MapRepP(cr.Rep()));
    entitySet boundary_cells = Loci::all_collect_entitySet(local_boundary_cells) ;
   
    store<string> boundary_names ;
    store<string> boundary_tags ;
    boundary_names.allocate(boundary_cells) ;
    boundary_tags.allocate(boundary_cells) ;
    
 
    FORALL(boundary_cells, bc) {
      char buf[512] ;
      bzero(buf,512) ;
      snprintf(buf,511,"BC_%d",-bc) ;
      boundary_names[bc] = string(buf) ;
      boundary_tags[bc] = string(buf) ;
    } ENDFORALL ;
    
  
    for(size_t i=0;i<boundary_ids.size();++i) {
      int id = boundary_ids[i].first ;
      if(boundary_cells.inSet(-id))
        boundary_names[-id] = boundary_ids[i].second ;
    }
  
    if(Loci::MPI_rank == 0) {
      cout << "boundaries identified as:" ;
      FORALL(boundary_cells, bc) {
        cout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      cout << endl ;
    }
    
    
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
      std::ostringstream oss ;
      oss << "volumeTag(" << volTags[i].first << ")" ;
      facts.create_fact(oss.str(),Tag) ;
    }
  
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(read_file_timer);
    double t2 = MPI_Wtime() ;
    if(MPI_rank == 0)
      cout << "Time to read in repartition grid is " << t2-t1
	   << endl ;
 
    Loci::memSpace("returning from FVM grid reader") ;
    return true ;
  }
 
  bool setupFVMGridFromContainer(fact_db &facts,
                                 vector<entitySet>& local_nodes,
                                 vector<entitySet>& local_faces,
                                 vector<entitySet>& local_cells,
                                 store<vector3d<double> >& t_pos,
                                 Map& tmp_cl,
                                 Map& tmp_cr,
                                 multiMap& tmp_face2node,
                                 vector<pair<int,string> >& boundary_ids,
                                 vector<pair<string,entitySet> >& volTags ) {
       
    if(!inputFVMGrid(facts,
                     local_nodes,
                     local_faces,
                     local_cells,
                     t_pos,
                     tmp_cl,
                     tmp_cr,
                     tmp_face2node,
                     boundary_ids,
                     volTags))
      return false ;
    Loci::memSpace("before create_face_info") ;
    create_face_info(facts) ;
    create_ref(facts) ;
    create_ghost_cells(facts) ;
    return true ;
  }
}

#include "FVMAdapt/defines.h"

namespace Loci {
  //create a new store new_pos from pos and inner_nodes
  //re_number the nodes in pos and inner_nodes,
  //and then redistribute them across the processes
  void createVOGNode(store<vector3d<double> > &new_pos,
                     const store<Loci::FineNodes> &inner_nodes,
                     int& num_nodes,
                     fact_db & facts,//in global numbering
                     vector<entitySet>& nodes_ptn
                     ){
    //get store pos
    store<vector3d<double> > pos;
    pos =  facts.get_variable("pos");
       
    if(MPI_processes == 1){
      //firsr write out numNodes
      long  num_original_nodes  = pos.domain().size();
      long  num_inner_nodes  = 0;
      FORALL(inner_nodes.domain(), cc){
        num_inner_nodes += inner_nodes[cc].size();
      }ENDFORALL;
    
      int node_base = 0;
      long npnts = num_original_nodes + num_inner_nodes;
      entitySet new_domain= interval(node_base, node_base+npnts-1);
      new_pos.allocate(new_domain);
    
      entitySet::const_iterator nei = new_domain.begin();
      entitySet::const_iterator ei = pos.domain().begin();
    
      for(long count = 0; count < num_original_nodes; count++, nei++, ei++){
        new_pos[*nei] = pos[*ei];
      }

      Loci::constraint faces, geom_cells;
      faces = facts.get_variable("faces");
      geom_cells = facts.get_variable("geom_cells");
      entitySet local_edges = facts.get_variable("edge2node")->domain();

      FORALL(local_edges, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL;
      FORALL(*geom_cells, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL; 
      FORALL(*faces, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL;
      num_nodes = new_domain.size();
      nodes_ptn.resize(1);
      nodes_ptn[0] = new_domain;
      return;
    }

 
 
  
    store<vector3d<double> > pos_t; //temp container 
    vector<entitySet> temp_node_ptn;
    entitySet my_temp_nodes;
    {
      //reorder store first, from global to io entities
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      constraint  my_faces, my_geom_cells, my_edges; 
   
      my_faces = facts.get_variable("faces");
      my_geom_cells = facts.get_variable("geom_cells");
      my_edges =  facts.get_variable("edges");
    
      entitySet local_edges = *my_edges;
      entitySet local_faces =  *my_faces;
      entitySet local_cells =   *my_geom_cells;
      entitySet local_nodes =  pos.domain();

      //before write out,create stores for pos and inner_nodes which
      //are ordered across processors in the file numbering, the domain of this container
      //shifted by offset is the actual file numbering. offset will be modified after function call  
     
      int noffset, eoffset, coffset, foffset;
      noffset = 0;
      store<vector3d<double> > pos_io;
      pos_io = Loci::Global2FileOrder(pos.Rep(), local_nodes, noffset, dist, MPI_COMM_WORLD) ;
      entitySet file_nodes = pos_io.domain(); 
    
    
      eoffset = 0;
      store<Loci::FineNodes> edge_inner_nodes;
      edge_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_edges,eoffset,dist,MPI_COMM_WORLD) ;
      entitySet file_edges = edge_inner_nodes.domain();
      
      coffset= 0;
      // Create container 
      store<Loci::FineNodes> cell_inner_nodes;
      cell_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_cells,coffset,dist,MPI_COMM_WORLD) ;
      entitySet file_cells = cell_inner_nodes.domain();
  
      foffset= 0;
      // Create container
      store<Loci::FineNodes> face_inner_nodes;
      face_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_faces,foffset,dist,MPI_COMM_WORLD) ;
      entitySet file_faces = face_inner_nodes.domain();
      
      //Now allocate temp entitySet in file numbering
   
   
      //compute the size of pos
      int num_pos_nodes = file_nodes.size();
      std::vector<int> pos_sizes(Loci::MPI_processes) ;
      pos_sizes = Loci::all_collect_sizes(num_pos_nodes);
      
      //compute the size of inner_nodes
      int num_local_edge_nodes = 0;
      FORALL(file_edges, cc){
        num_local_edge_nodes += edge_inner_nodes[cc].size();
      }ENDFORALL;
  
      int num_local_cell_nodes = 0;
      FORALL(file_cells, cc){
        num_local_cell_nodes += cell_inner_nodes[cc].size();
      }ENDFORALL;
  
      int num_local_face_nodes = 0;
      FORALL(file_faces, cc){
        num_local_face_nodes += face_inner_nodes[cc].size();
      }ENDFORALL;
      std::vector<int> inner_edge_nodes_sizes(Loci::MPI_processes) ;
      inner_edge_nodes_sizes = Loci::all_collect_sizes(num_local_edge_nodes);
    
      std::vector<int> inner_cell_nodes_sizes(Loci::MPI_processes) ;
      inner_cell_nodes_sizes = Loci::all_collect_sizes(num_local_cell_nodes);

      std::vector<int> inner_face_nodes_sizes(Loci::MPI_processes) ;
      inner_face_nodes_sizes = Loci::all_collect_sizes(num_local_face_nodes);
    
      //set up the offset
      int my_id = Loci::MPI_rank;
      int num_procs = Loci::MPI_processes;
      noffset = 0;
      for(int i = 0; i < my_id; i++){
        noffset += pos_sizes[i];
      }
    
      eoffset = 0;
      for(int i = 0; i < num_procs; i++){
        eoffset += pos_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        eoffset += inner_edge_nodes_sizes[i];
      }
    
      coffset = 0;
      for(int i = 0; i < num_procs; i++){
        coffset += pos_sizes[i]+inner_edge_nodes_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        coffset += inner_cell_nodes_sizes[i];
      }
    
    
      foffset = 0;
      for(int i = 0; i < num_procs; i++){
        foffset += pos_sizes[i]+inner_edge_nodes_sizes[i]+inner_cell_nodes_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        foffset += inner_face_nodes_sizes[i];
      }
    

      num_nodes = 0;
      for(int i = 0; i < num_procs; i++){
        num_nodes += pos_sizes[i]+inner_edge_nodes_sizes[i]+
          inner_cell_nodes_sizes[i]+inner_face_nodes_sizes[i];
      }
      
      //create a new store pos_t in file numbering that contains all the nodes  
      entitySet my_original_nodes, my_inner_edge_nodes;
      entitySet my_inner_cell_nodes, my_inner_face_nodes;
      if(pos_sizes[my_id]>0) my_original_nodes = interval(noffset, noffset+pos_sizes[my_id]-1);
      if(inner_edge_nodes_sizes[my_id]>0)my_inner_edge_nodes = interval(eoffset, eoffset+inner_edge_nodes_sizes[my_id]-1);
      if(inner_cell_nodes_sizes[my_id]>0)my_inner_cell_nodes = interval(coffset, coffset+inner_cell_nodes_sizes[my_id]-1);
      if(inner_face_nodes_sizes[my_id]>0)my_inner_face_nodes = interval(foffset, foffset+inner_face_nodes_sizes[my_id]-1);
      my_temp_nodes = my_original_nodes + my_inner_edge_nodes + 
        my_inner_cell_nodes + my_inner_face_nodes;
      pos_t.allocate(my_temp_nodes);
    
      //fill container pos_t will pos_io
      entitySet::const_iterator ti = my_temp_nodes.begin();
      FORALL(file_nodes, ei){
        pos_t[*ti] = pos_io[ei];
        ti++;
      }ENDFORALL;
      //fill container pos_t will edge nodes
      FORALL(file_edges, ei){
        for(unsigned int i = 0; i < edge_inner_nodes[ei].size(); i++){
          pos_t[*ti] = edge_inner_nodes[ei][i];
          ti++;
        }
     
      }ENDFORALL;
      //fill container pos_t will cell nodes
      FORALL(file_cells, ei){
        for(unsigned int i = 0; i < cell_inner_nodes[ei].size(); i++){
          pos_t[*ti] = cell_inner_nodes[ei][i];
          ti++;
        }
      
      }ENDFORALL;
      //fill container pos_t will face nodes
      FORALL(file_faces, ei){
        for(unsigned int i = 0; i < face_inner_nodes[ei].size(); i++){
          pos_t[*ti] = face_inner_nodes[ei][i];
          ti++;
        }
      }ENDFORALL;
      if(ti != my_temp_nodes.end()){
        debugout<<"ERROR in createVOGNOde()" << endl;
      }
      temp_node_ptn = all_collect_vectors(my_temp_nodes);
    }
    //create new partition  
    nodes_ptn.resize(Loci::MPI_processes);
    // create node allocation
    long npnts = num_nodes ;
    int node_ivl = npnts / Loci::MPI_processes;
    int node_ivl_rem = npnts % Loci::MPI_processes ;
    int node_accum = 0 ;
    int nodes_base = 0 ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = Loci::MPI_processes - i - 1 ;
      int node_accum_update = node_accum + node_ivl + ((j<node_ivl_rem)?1:0) ;
      if(i == Loci::MPI_processes-1) {
        nodes_ptn[i] = interval(nodes_base + node_accum,
                                nodes_base + npnts - 1) ;
      } else {
        nodes_ptn[i] = interval(nodes_base + node_accum,
                                nodes_base + node_accum_update - 1) ;
      }
      node_accum = node_accum_update ;
    }
    
    entitySet new_nodes = nodes_ptn[MPI_rank];
    new_pos.allocate(new_nodes) ;
  
    // Now compute where to send data
    int p = Loci::MPI_processes;
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;
    
    // Check each processor, find out which sets to send
    for(int i=0;i<p;++i) {
      send_sets[i] = nodes_ptn[i]&my_temp_nodes;
      sequence s= send_sets[i] ;
      send_seqs[i] = s ;
    }
    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;
  
    storeRepP sp =pos_t.Rep();
    storeRepP qcol_rep =new_pos.Rep();
  
  
    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;
    
    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD ) ;
    
    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  MPI_COMM_WORLD ) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    pos_t.allocate(EMPTY) ;
  }


  //create cl, cr and face2node maps from fine_faces
  //re_number the nodes in pos and inner_nodes,
  //and then redistribute them across the processes
  void createVOGFace(int numNodes,
                     const store<Loci::FineFaces> &fine_faces,
                     fact_db & facts,
                     int& numFaces,
                     int& ncells,
                     Map& cl,
                     Map& cr,
                     multiMap& face2node,
                     vector<entitySet>& local_faces,
                     vector<entitySet>& local_cells
                     ){

    //compute numFaces
    constraint my_faces, my_geom_cells;
    my_faces = facts.get_variable("faces");
    my_geom_cells = facts.get_variable("geom_cells");
    entitySet dom = fine_faces.domain() & (*my_faces+*my_geom_cells);
  
    int local_num_face = 0;
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      local_num_face += fine_faces[*ei].size();
    }
    numFaces = 0;
    std::vector<int> face_sizes= Loci::all_collect_sizes(local_num_face);
    for(int i =0; i < MPI_processes; i++) numFaces += face_sizes[i];
  
    //get face domains and allocate the maps 
    int face_min = numNodes;
    for(int i =0; i < MPI_rank; i++) face_min += face_sizes[i];
    int face_max = face_min + face_sizes[MPI_rank] -1;
    store<int> count;
    entitySet faces = interval(face_min, face_max);
    cl.allocate(faces);
    cr.allocate(faces);
    count.allocate(faces);
    local_faces.resize(MPI_processes);
    local_faces = all_collect_vectors(faces);
    //fill up the maps cl , cr and count 
    int cell_base = numNodes + numFaces;
    int cell_max = std::numeric_limits<int>::min();
    int cell_min = std::numeric_limits<int>::max();
    entitySet::const_iterator fid = faces.begin();
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
        cl[*fid] = fine_faces[*ei][i][0] + cell_base -1;// -1 finefaces cell index start at 1
        if(fine_faces[*ei][i][1]>=0) cr[*fid] = fine_faces[*ei][i][1] + cell_base -1;
        else  cr[*fid] = fine_faces[*ei][i][1];
        cell_max = max(cell_max,cl[*fid]);
        cell_max = max(cell_max,cr[*fid]);
        cell_min = min(cell_min,cl[*fid]);
        if (cr[*fid]>=0) cell_min = min(cell_min,cr[*fid]);
        count[*fid] = fine_faces[*ei][i].size()-2;
        fid++;
      }
    }
  
    //get cells distribution
    int global_max_cell = cell_max ;
    MPI_Allreduce(&cell_max,&global_max_cell,1,MPI_INT,MPI_MAX,
                  MPI_COMM_WORLD) ;
    int global_min_cell = cell_max ;
    MPI_Allreduce(&cell_min,&global_min_cell,1,MPI_INT,MPI_MIN,
                  MPI_COMM_WORLD) ;
  
  
    ncells = global_max_cell - global_min_cell +1;
  

    local_cells.resize(MPI_processes);
    int cell_ivl = ncells / MPI_processes;
    int cell_ivl_rem = ncells % MPI_processes ;
    int cell_accum = 0 ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = MPI_processes - i - 1 ;
      int cell_accum_update = cell_accum + cell_ivl + ((j<cell_ivl_rem)?1:0) ;
    
      if(i == MPI_processes-1) {
        local_cells[i] = interval(cell_base + cell_accum,
                                  cell_base + ncells-1) ;
      } else {
        local_cells[i] = interval(cell_base + cell_accum,
                                  cell_base + cell_accum_update - 1) ;
      }
      cell_accum = cell_accum_update ;
    }
    //fill up face2node  
    face2node.allocate(count);
    fid = faces.begin();
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
        for(int j = 0; j < count[*fid]; j++){
          //vog file node index start with 0
          // face2node[*fid][j] = fine_faces[*ei][i][j+2]-1;
          // if(fine_faces[*ei][i][j+2] < 0) cerr <<"WARNING: negative node index" << endl;
          face2node[*fid][j] = fine_faces[*ei][i][j+2];
        }
        fid++;
      }
    }
    colorMatrix(cl, cr, face2node);
  }
  
}

namespace Loci{
 

  entitySet faceCluster(const multiMap &face2node,
                        const Map &cl, const Map &cr, entitySet faces,
                        vector<unsigned char> &cluster_info,
                        vector<unsigned short> &cluster_sizes) ;
  bool readBCfromVOG(string filename,
                    vector<pair<int,string> > &boundary_ids);
  
  
  hid_t writeVOGOpen(string filename);
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids);
  void writeVOGClose(hid_t file_id) ;


}

void writeVOGNode(hid_t file_id,
                Loci::storeRepP &pos,
                  const_store<Loci::FineNodes> &inner_nodes);


void colorMatrix(Map &cl, Map &cr, multiMap &face2node);
namespace Loci{
  hid_t writeVOGOpen(string filename) {
    hid_t file_id = 0 ;
    if(MPI_rank==0) 
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
    return file_id ;
  }
   void writeVOGClose(hid_t file_id) {
     if(MPI_rank == 0) H5Fclose(file_id) ;
   }
  
       
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids) {
    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      if(surface_ids.size() != 0) {
        group_id = H5Gcreate(file_id,"surface_info",0) ;
        for(size_t i=0;i<surface_ids.size();++i) {
          hid_t bc_id = 0 ;
          bc_id = H5Gcreate(group_id,surface_ids[i].second.c_str(),0) ;
          hsize_t dims = 1 ;
          hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
          
          hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                                   dataspace_id, H5P_DEFAULT) ;
          H5Awrite(att_id,H5T_NATIVE_INT,&surface_ids[i].first) ;
          H5Aclose(att_id) ;
          H5Gclose(bc_id) ;
        }
        H5Gclose(group_id) ;
      }
    }
  }


}
void writeVOGFace(hid_t file_id, Map &cl, Map &cr, multiMap &face2node) {
    // Compute cell set
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
  
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;

  
    long long local_num_faces = face2node.domain().size()  ;

  
    long long num_cells = geom_cells.size() ;
    long long num_faces = 0 ;

    // Reduce these variables
    MPI_Allreduce(&local_num_faces,&num_faces,1,MPI_LONG_LONG_INT,
                  MPI_SUM,MPI_COMM_WORLD) ;

    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      group_id = H5Gopen(file_id,"file_info") ;

      std::cout<< "num_cells = " << num_cells << endl
               << "num_faces = " << num_faces << endl ;

      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
      hid_t att_id = H5Acreate(group_id,"numFaces", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_faces) ;
      H5Aclose(att_id) ;
      att_id = H5Acreate(group_id,"numCells", H5T_STD_I64BE,
                         dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_cells) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;
      group_id = H5Gcreate(file_id,"face_info",0) ;
    }
  
    entitySet faces = face2node.domain() ;
    vector<pair<pair<int,int>, int> > f_ord(faces.size()) ;
    int i = 0 ;
    // For small number of cells, sort to keep bc groupings
    if(num_cells<100000) {
      FORALL(faces,fc) {
        f_ord[i].first.first = cr[fc] ;
        f_ord[i].first.second = cl[fc] ;
        f_ord[i].second = fc ;
        i++ ;
      } ENDFORALL ;
      std::sort(f_ord.begin(),f_ord.end()) ;
    } else {
      FORALL(faces,fc) {
        f_ord[i].first.first = cl[fc] ;
        f_ord[i].first.second = cr[fc] ;
        f_ord[i].second = fc ;
        i++ ;
      } ENDFORALL ;
    }

    i=0 ;
    store<int> count ;
    count.allocate(faces) ;
    FORALL(faces,fc) {
      int nfc = f_ord[i].second ;
      count[fc] = face2node[nfc].size() ;
      i++ ;
    } ENDFORALL ;
    tmp_face2node.allocate(count) ;
    tmp_cl.allocate(faces) ;
    tmp_cr.allocate(faces) ;
    i=0 ;
  
    int mc = (geom_cells).Min() ;
    // Nodes should be adjusted to start from zero also... for the general case
    FORALL(faces,fc) {
      int nfc = f_ord[i].second ;
      tmp_cl[fc] = cl[nfc]-mc ;
      tmp_cr[fc] = cr[nfc] ;
      if(tmp_cr[fc] >= 0)
        tmp_cr[fc] -= mc ;
      for(int j=0;j<count[fc];++j)
        tmp_face2node[fc][j] = face2node[nfc][j] ;
      i++ ;
    } ENDFORALL ;

    vector<unsigned char> cluster_info ;
    vector<unsigned short> cluster_sizes ;
    while(faces != EMPTY) {
      entitySet fcluster = Loci::faceCluster(tmp_face2node,tmp_cl,tmp_cr,faces,
                                             cluster_info,cluster_sizes) ;
      faces -= fcluster ;
    }

    Loci::writeUnorderedVector(group_id,"cluster_sizes",cluster_sizes) ;
    Loci::writeUnorderedVector(group_id,"cluster_info",cluster_info) ;
  
  
    if(MPI_rank == 0) {
      H5Gclose(group_id) ;
    }
  }



// void writeVOG(string filename,store<vector3d<double> > &pos,
//               store<Loci::FineNodes> & inner_nodes,
//               Map &cl, Map &cr, multiMap &face2node,
//               vector<pair<int,string> > surface_ids) {
//     // write grid file
//   hid_t file_id = writeVOGOpen(filename) ;
//   writeVOGSurf(file_id,surface_ids) ;
//   writeVOGNode(file_id,pos) ;
//   writeVOGFace(file_id,cl,cr,face2node) ;
//   writeVOGClose(file_id) ;
// }





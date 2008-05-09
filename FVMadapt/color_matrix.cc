#include <hdf5.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include "sciTypes.h"
#include "defines.h"
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;
using Loci::MPI_processes;
using Loci::MPI_rank;
using Loci::entitySet;
using Loci::UNIVERSE_MAX;
using Loci::UNIVERSE_MIN;

namespace Loci{
  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm);
  std::vector<int> all_collect_sizes(int size);
  std::vector<int> simplePartitionVec(int mn, int mx, int p);
  vector<sequence> transposeSeq(const vector<sequence> sv) ;
  storeRepP collect_reorder_store(storeRepP &sp, fact_db &facts){
    entitySet dom = sp->domain() ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    constraint my_entities ; 
    my_entities = d->my_entities ;
    dom = *my_entities & dom ;


    Map l2g ;
    l2g = d->l2g.Rep() ;
    MapRepP l2gP = MapRepP(l2g.Rep()) ;

    entitySet dom_global = l2gP->image(dom) ;

    FATAL(dom.size() != dom_global.size()) ;

    dMap g2f ;
    g2f = d->g2f.Rep() ;
    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[l2g[i]] ;
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;

    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;
    
    int p = MPI_processes ;

    
    vector<entitySet> out_ptn = simplePartition(imn,imx,MPI_COMM_WORLD) ;

    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;
    for(int i=0;i<p;++i) {
      send_sets[i] = newnum.preimage(out_ptn[i]).first ;
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;

    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;


    int *send_sizes = new int[p] ;
    int *recv_sizes = new int[p] ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int *send_dspl = new int[p] ;
    int *recv_dspl = new int[p] ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    unsigned char *send_store = new unsigned char[send_sz] ;
    unsigned char *recv_store = new unsigned char[recv_sz] ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(send_store, &send_sizes[0], send_dspl, MPI_PACKED,
		  recv_store, &recv_sizes[0], recv_dspl, MPI_PACKED,
		  MPI_COMM_WORLD) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    delete[] recv_store ;
    delete[] send_store ;
    delete[] recv_dspl ;
    delete[] send_dspl ;
    delete[] recv_sizes ;
    delete[] send_sizes ;

    return qcol_rep ;
  } 


}



void writeVOGNode(hid_t file_id,
                Loci::storeRepP &pos,
                const_store<Loci::FineNodes> &inner_nodes){

  hid_t group_id = 0 ;
  
  if(MPI_processes == 1){
    //firsr write out numNodes
    long long num_original_nodes  = pos->domain().size();
    long long num_inner_nodes  = 0;
      FORALL(inner_nodes.domain(), cc){
      num_inner_nodes += inner_nodes[cc].size();
      }ENDFORALL;
      
      long long array_size = num_original_nodes + num_inner_nodes;
      
      group_id = H5Gcreate(file_id,"file_info",0) ;

      cerr << "num_nodes = " << array_size << endl ;
      
      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
      
      hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&array_size) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;

      if(array_size == 0)
        return ;


      //prepare to write positions
      group_id = H5Gcreate(file_id,"node_info",0) ;
      
      //create dataspace and dataset
      int rank = 1 ;
      hsize_t dimension = array_size ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      typedef data_schema_traits<vect3d> traits_type ;
      Loci::DatatypeP dp = traits_type::get_type() ;

#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      
      hid_t dataset = H5Dcreate(group_id,"positions",dp->get_hdf5_type(),
                                dataspace, H5P_DEFAULT) ;

      
      //first write out pos
      hsize_t count = num_original_nodes ;
      if(count != 0) {
        std::vector<vect3d> v_pos(count);
        //put pos_io in a vector
        store<vect3d>pos_io;
        pos_io = pos;
        int index = 0;
        FORALL(pos_io.domain(),cc){
          v_pos[index++] = pos_io[cc];
        }ENDFORALL;
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &v_pos[0]) ;
        H5Sclose(memspace) ;
      }
      
      //put inner_nodes in a vector
      Loci::constraint faces, geom_cells;
      faces = Loci::exec_current_fact_db->get_variable("faces");
      geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
      entitySet local_edges = Loci::exec_current_fact_db->get_variable("edge2node")->domain();
    
   
      
      //next, write out inner_nodes
      start += num_original_nodes  ;
      
      count = num_inner_nodes;
      if(num_inner_nodes != 0) {
        std::vector<vect3d> v_nodes(num_inner_nodes);
        //put inner_nodes in a vector

        long index = 0;
        FORALL(local_edges, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
            v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
        FORALL(*geom_cells, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL; 
        FORALL(*faces, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
                
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &v_nodes[0]) ;
        H5Sclose(memspace) ;
      }

      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Gclose(group_id) ;
      return;
  }
  
  //reorder store first, from local to io entities
  
  store<vect3d> pos_io;
  pos_io =  Loci::collect_reorder_store(pos, *(Loci::exec_current_fact_db));
  entitySet nodes = pos_io.domain();

  //compute the size of pos
  int local_pos_size = nodes.size();
  std::vector<int> pos_sizes(Loci::MPI_processes) ;
  MPI_Gather(&local_pos_size,1,MPI_INT,
             &pos_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;



  //compute the size of inner_nodes

  fact_db::distribute_infoP d =  Loci::exec_current_fact_db->get_distribute_info() ;
  Loci::constraint my_entities, faces, geom_cells, interior_faces, boundary_faces; 
  my_entities = d->my_entities ;
  faces = Loci::exec_current_fact_db->get_variable("faces");
  geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
  interior_faces =  Loci::exec_current_fact_db->get_variable("interior_faces");
  boundary_faces =  Loci::exec_current_fact_db->get_variable("boundary_faces");
  
  entitySet local_edges =*my_entities &
    (Loci::exec_current_fact_db->get_variable("edge2node"))->domain(); 
  entitySet local_faces = *my_entities & *faces;
  entitySet local_cells =  *my_entities & *geom_cells;
  entitySet local_interior_faces = *my_entities & *interior_faces;
  entitySet local_boundary_faces = *my_entities & *boundary_faces;
  

  int num_local_inner_nodes = 0;
  FORALL(local_edges, cc){
    num_local_inner_nodes += inner_nodes[cc].size();
  }ENDFORALL;
  FORALL(local_cells, cc){
    num_local_inner_nodes += inner_nodes[cc].size();
  }ENDFORALL; 
  FORALL(local_faces, cc){
    num_local_inner_nodes += inner_nodes[cc].size();
  }ENDFORALL;
  
  std::vector<int> inner_nodes_sizes(Loci::MPI_processes);
  MPI_Gather(&num_local_inner_nodes,1,MPI_INT,
             &inner_nodes_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;
  
  if(Loci::MPI_rank == 0) {
    //compute array_size
      hsize_t array_size = 0 ;
      for(int i=0;i<MPI_processes;++i)
        array_size += (pos_sizes[i]+ inner_nodes_sizes[i]);
      //first write out numNodes

      group_id = H5Gcreate(file_id,"file_info",0) ;

      cerr << "num_nodes = " << array_size << endl ;
      
      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
      
      hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&array_size) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;


      
      if(array_size == 0)
        return ;



      group_id = H5Gcreate(file_id,"node_info",0) ;

      //create dataspace and dataset
      int rank = 1 ;
      hsize_t dimension = array_size ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      typedef data_schema_traits<vect3d> traits_type ;
      Loci::DatatypeP dp = traits_type::get_type() ;

#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
    
      hid_t dataset = H5Dcreate(group_id,"positions",dp->get_hdf5_type(),
                                dataspace, H5P_DEFAULT) ;

      
      //first write out pos
      hsize_t count = pos_sizes[0] ;
      if(count != 0) {
        std::vector<vect3d> v_pos(count);
        //put pos_io in a vector
        int index = 0;
        FORALL(nodes,cc){
          v_pos[index++] = pos_io[cc];
        }ENDFORALL;
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &v_pos[0]) ;
        H5Sclose(memspace) ;
      }
      for(int i=1;i<MPI_processes;++i) {
        start += pos_sizes[i-1] ;
        if(pos_sizes[i] == 0)
          continue ;
        int flag = 0 ;
        MPI_Send(&flag,1,MPI_INT,i,0,MPI_COMM_WORLD) ;
        std::vector<vect3d> rv(pos_sizes[i]) ;
        MPI_Status mstat ;
        MPI_Recv(&rv[0],sizeof(vect3d)*pos_sizes[i],MPI_BYTE,i,1,MPI_COMM_WORLD,
                 &mstat) ;
        count = pos_sizes[i] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &rv[0]) ;
        H5Sclose(memspace) ;
      }
      //next, write out inner_nodes
      start += pos_sizes[MPI_processes-1] ;
      
      count = num_local_inner_nodes;
      if(num_local_inner_nodes != 0) {
        std::vector<vect3d> v_nodes(num_local_inner_nodes);
        //put inner_nodes in a vector

        long index = 0;
        FORALL(local_edges, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
            v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
        FORALL(local_cells, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL; 
        FORALL(local_faces, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
          v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
                
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &v_nodes[0]) ;
        H5Sclose(memspace) ;
      }

      for(int i=1;i<MPI_processes;++i) {
        start += inner_nodes_sizes[i-1] ;
        if(inner_nodes_sizes[i] == 0)
          continue ;
        int flag = 0 ;
        MPI_Send(&flag,1,MPI_INT,i,2,MPI_COMM_WORLD) ;
        std::vector<vect3d> rv(inner_nodes_sizes[i]) ;
        MPI_Status mstat ;
        MPI_Recv(&rv[0],sizeof(vect3d)*inner_nodes_sizes[i],MPI_BYTE,i,3,MPI_COMM_WORLD,
                 &mstat) ;
        count = inner_nodes_sizes[i] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &rv[0]) ;
        H5Sclose(memspace) ;
      }

        
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      
    } else {
      if(local_pos_size != 0){
          std::vector<vect3d> v_pos(local_pos_size);
          //put pos_io in a vector
          int index = 0;
          FORALL(nodes,cc){
            v_pos[index++] = pos_io[cc];
          }ENDFORALL;
          
          int flag = 0;
          MPI_Status mstat ;
          MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&mstat) ;
          MPI_Send(&v_pos[0],sizeof(vect3d)*local_pos_size,MPI_BYTE,0,1,MPI_COMM_WORLD) ;
      }
      if(num_local_inner_nodes != 0){


        std::vector<vect3d> v_nodes(num_local_inner_nodes);
        //put inner_nodes in a vector
        
        int index = 0;
        FORALL(local_edges, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
            v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
        FORALL(local_cells, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
            v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL; 
        FORALL(local_faces, cc){
          for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
            v_nodes[index++] = inner_nodes[cc][i];
          }
        }ENDFORALL;
        
        
        int flag = 0;
        MPI_Status mstat ;
        MPI_Recv(&flag,1,MPI_INT,0,2,MPI_COMM_WORLD,&mstat) ;
        MPI_Send(&v_nodes[0],sizeof(vect3d)*num_local_inner_nodes,MPI_BYTE,0,3,MPI_COMM_WORLD) ;
        
      }
      
    }

  
  if(Loci::MPI_rank == 0) H5Gclose(group_id) ;
  
  
  
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

void colorMatrix(Map &cl, Map &cr, multiMap &face2node) {
    
    entitySet  faces,cells ;
   std::vector<entitySet> ptn = getDist(faces,cells,
                                    cl,cr,face2node);
    entitySet loc_faces = faces & ptn[MPI_rank] ;
    entitySet geom_cells = cells & ptn[MPI_rank] ;
    entitySet negs = interval(UNIVERSE_MIN,-1) ;
    entitySet boundary_faces = cr.preimage(negs).first ;
    entitySet interior_faces = loc_faces - boundary_faces ;

    using std::pair ;
   std:: vector<pair<Entity,Entity> > cellmap(interior_faces.size()*2) ;
    int cnt = 0 ;
    FORALL(interior_faces,fc) {
      cellmap[cnt++] = pair<Entity,Entity>(cl[fc],cr[fc]) ;
      cellmap[cnt++] = pair<Entity,Entity>(cr[fc],cl[fc]) ;
    } ENDFORALL ;
    multiMap c2c ;
    Loci::distributed_inverseMap(c2c,cellmap,cells,cells,ptn) ;
    int ncells = cells.size() ;

    store<int> ctmp ;
    ctmp.allocate(geom_cells) ;
    FORALL(geom_cells,cc) {
      ctmp[cc] = -1 ;
    } ENDFORALL ;

    int col = ncells*Loci::MPI_rank ;
    
   std:: vector<int> visited ;
    entitySet left_out = geom_cells ;
    int lo_p = geom_cells.Min() ;
    while(left_out != EMPTY) {
     std:: vector<int> work ;
      work.push_back(left_out.Min()) ;
      while(work.size() != 0) {
	std::vector<int> working ;
	for(size_t i=0;i<work.size();++i) {
          int cc = work[i] ;
          if(ctmp[cc] == -1) {
            ctmp[cc] = col++ ;
            visited.push_back(cc) ;
            for(const int *pi = c2c.begin(cc);pi!=c2c.end(cc);++pi)
              if(geom_cells.inSet(*pi) && ctmp[*pi] == -1) {
                working.push_back(*pi) ;
              }
          }
	}
        work.swap(working) ;
      }
      left_out = EMPTY ;
      entitySet candidates = geom_cells & interval(lo_p,UNIVERSE_MAX) ;
      FORALL(candidates,cc) {
        if(ctmp[cc] == -1) {
          left_out += cc ;
          break ;
        }
      } ENDFORALL ;
      if(left_out != EMPTY)
        lo_p = left_out.Min() ;
    }

    dstore<int> color ;
    FORALL(geom_cells,cc) {
      color[cc] = ctmp[cc];
    } ENDFORALL ;

    entitySet clone_cells = cl.image(interior_faces)
      + cr.image(interior_faces) ;
    clone_cells -= geom_cells ;
    Loci::storeRepP cp_sp = color.Rep() ;
    Loci::fill_clone(cp_sp, clone_cells, ptn) ;

    FORALL(interior_faces,fc) {
      int color_l = color[cl[fc]] ;
      int color_r = color[cr[fc]] ;
      if(color_l == color_r) 
        cerr << "color equal == " << color_l << endl ;
      if(color_l == -1 |+ color_r == -1)
        cerr << "matrix coloring internal error" << endl ;
                                                              
      if(color_l > color_r) {
        // change face orientation to match matrix coloring
        std::swap(cl[fc],cr[fc]) ;
        int i = 0 ;
        int j = face2node[fc].size() - 1;
       	while(i < j) {
          std::swap(face2node[fc][i],face2node[fc][j]) ;
          i++ ;
          j-- ;
        } 
      }
    } ENDFORALL ;
    
  }

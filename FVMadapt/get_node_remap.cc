#include <Loci.h>
#include <map>
#include "defines.h"
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
#include <algorithm>
using std::sort ;
using std::unique ;

//#include "dist_tools.h"
using std::cout ;


using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using Loci::MPI_processes;
void
par_sort(vector<pair<int,int> >& data, MPI_Comm comm){
  // first get the processor id and total number of processors
  int my_id, num_procs ;
  MPI_Comm_size(comm, &num_procs) ;
  MPI_Comm_rank(comm, &my_id) ;
  if(num_procs <= 1)
    return ;                  // single process, no need to proceed
  // get the number of local elements
  int local_size = data.size() ;
  // then select num_procs-1 equally spaced elements as splitters
  int* splitters = new int[num_procs] ;
  int even_space = local_size / (num_procs-1) ;
  int start_idx = even_space / 2 ;
  int space_idx = start_idx ;
  for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
    splitters[i] = data[space_idx].first ;
  // gather the splitters to all processors as samples
  int sample_size = num_procs * (num_procs-1) ;
  int* samples = new int[sample_size] ;
  MPI_Allgather(splitters, num_procs-1, MPI_INT,
                samples, num_procs-1, MPI_INT, comm) ;
  // now we've obtained all the samples, first we sort them
  sort(samples, samples+sample_size) ;
  // select new splitters in the sorted samples
  even_space = sample_size / (num_procs-1) ;
  start_idx = even_space / 2 ;
  space_idx = start_idx ;
  for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
    splitters[i] = samples[space_idx] ;
  // the last one set as maximum possible integer
  splitters[num_procs-1] = std::numeric_limits<int>::max() ;
  
  // now we can assign local elements to buckets (processors)
  // according to the new splitters. first we will compute
  // the size of each bucket and communicate them first
  int* scounts = new int[num_procs] ;
  for(int i=0;i<num_procs;++i)
    scounts[i] = 0 ;
  { // using a block just to make the definition of "i" and "j" local
    int i, j ;
    for(j=i=0;i<local_size;++i) {
      if(data[i].first < splitters[j])
        scounts[j]++ ;
      else {
        ++j ;
        while(data[i].first >= splitters[j]) {
          scounts[j] = 0 ;
          ++j ;
          }
        scounts[j]++ ;
      }
    }
    }
  // but since one local element contains two integers (a pair of int),
    // we will need to double the size
  for(int i=0;i<num_procs;++i)
      scounts[i] *= 2 ;
  // now we compute the sending displacement for each bucket
  int* sdispls = new int[num_procs] ;
  sdispls[0] = 0 ;
  for(int i=1;i<num_procs;++i)
      sdispls[i] = sdispls[i-1] + scounts[i-1] ;
  // communicate this information to all processors so that each will
  // know how many elements are expected from every other processor
  int* rcounts = new int[num_procs] ;
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm) ;
  // then based on the received info. we will need to compute the
  // receive displacement
  int* rdispls = new int[num_procs] ;
  rdispls[0] = 0 ;
  for(int i=1;i<num_procs;++i)
    rdispls[i] = rdispls[i-1] + rcounts[i-1] ;
  // then we will need to pack the elements in local into
  // a buffer and communicate them
  int* local_pairs = new int[local_size*2] ;
  int count = 0 ;
  for(int i=0;i<local_size;++i) {
    local_pairs[count++] = data[i].first ;
    local_pairs[count++] = data[i].second ;
  }
  // then we allocate buffer for new local elements
  int new_local_size = rdispls[num_procs-1] + rcounts[num_procs-1] ;
  int* sorted_pairs = new int[new_local_size] ;
  // finally we communicate local_pairs to each processor
  MPI_Alltoallv(local_pairs, scounts, sdispls, MPI_INT,
                sorted_pairs, rcounts, rdispls, MPI_INT, comm) ;
  // release buffers
  delete[] splitters ;
  delete[] samples ;
  delete[] scounts ;
  delete[] sdispls ;
  delete[] rcounts ;
  delete[] rdispls ;
  delete[] local_pairs ;
  // finally we unpack the buffer into a vector of pairs
  data.resize(new_local_size/2) ;
  int data_idx = 0 ;
  for(int i=0;i<new_local_size;i+=2,data_idx++)
    data[data_idx] = pair<int,int>(sorted_pairs[i],sorted_pairs[i+1]) ;
  // release the final buffer
  delete[] sorted_pairs ;
  // finally we sort the new local vector
  sort(data.begin(), data.end()) ;
}



int
global_sum(int l) {
  int g ;
  MPI_Allreduce(&l, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
  return g ;
}


vector<entitySet>
gather_all_entitySet(const entitySet& eset){
  int local_size = eset.size() ;
  int global_size = global_sum(local_size) ;
  // compute receive counts from all processors
  int* recv_counts = new int[Loci::MPI_processes] ;
  MPI_Allgather(&local_size, 1, MPI_INT,
                recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;
  // then compute receive displacement
  int* recv_displs = new int[Loci::MPI_processes] ;
  recv_displs[0] = 0 ;
  for(int i=1;i<MPI_processes;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
      // pack the local eset into an array
      int* local_eset = new int[local_size] ;
      int count = 0 ;
      for(entitySet::const_iterator ei=eset.begin();
          ei!=eset.end();++ei,++count)
        local_eset[count] = *ei ;
      // allocate the entire array for all data from all processors
      int* global_eset = new int[global_size] ;
      // communicate to obtain all esets from every processors
      MPI_Allgatherv(local_eset, local_size, MPI_INT,
                     global_eset, recv_counts, recv_displs,
                     MPI_INT, MPI_COMM_WORLD) ;
      delete[] local_eset ;
      delete[] recv_counts ;
      // unpack the raw buffer into a vector<entitySet>
      vector<entitySet> ret(MPI_processes) ;
      int k = 0 ;
      for(int i=0;i<MPI_processes;++i) {
        int limit ;
        if(i == MPI_processes-1)
          limit = global_size ;
        else
          limit = recv_displs[i+1] ;
        for(;k<limit;++k)
          ret[i] += global_eset[k] ;
      }
      delete[] recv_displs ;
      delete[] global_eset ;

      return ret ;
}
void
parallel_balance_pair_vector(vector<pair<int,int> >& vp,
                             MPI_Comm comm) {
  int num_procs = 0 ;
  MPI_Comm_size(comm,&num_procs) ;
  
  // we still use an all-to-all personalized communication
  // algorithm to balance the element numbers on processes.
  // we pick (p-1) equally spaced element as the splitters
  // and then re-split the global vector sequence to balance
  // the number of elements on processes.
  
  int vp_size = vp.size() ;
  int global_vp_size = 0 ;
  MPI_Allreduce(&vp_size, &global_vp_size,
                1, MPI_INT, MPI_SUM, comm) ;
  
  int space = global_vp_size / num_procs ;
  // compute a global range for the elements on each process
  int global_end = 0 ;
  MPI_Scan(&vp_size, &global_end, 1, MPI_INT, MPI_SUM, comm) ;
  int global_start = global_end - vp_size ;
  
  vector<int> splitters(num_procs) ;
  // splitters are just global index number
  splitters[0] = space ;
  for(int i=1;i<num_procs-1;++i)
    splitters[i] = splitters[i-1] + space ;
  splitters[num_procs-1] = global_vp_size ;
  
  // split and communicate the vector of particles
  vector<int> send_counts(num_procs, 0) ;
  int part_start = global_start ;
  for(int idx=0;idx<num_procs;++idx) {
    if(part_start == global_end)
      break ;
    if(splitters[idx] > part_start) {
      int part_end ;
      if(splitters[idx] < global_end)
        part_end = splitters[idx] ;
      else
        part_end = global_end ;
      send_counts[idx] = part_end - part_start ;
      part_start = part_end ;
    }
  }
  
  for(size_t i=0;i<send_counts.size();++i)
    send_counts[i] *= 2 ;
  
  vector<int> send_displs(num_procs) ;
  send_displs[0] = 0 ;
  for(int i=1;i<num_procs;++i)
    send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
  
  vector<int> recv_counts(num_procs) ;
  MPI_Alltoall(&send_counts[0], 1, MPI_INT,
               &recv_counts[0], 1, MPI_INT, comm) ;
  
  vector<int> recv_displs(num_procs) ;
  recv_displs[0] = 0 ;
  for(int i=1;i<num_procs;++i)
    recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
  
  int total_recv_size = recv_displs[num_procs-1] +
    recv_counts[num_procs-1] ;
  
  // prepare send and recv buffer
  vector<int> send_buf(vp_size*2) ;
  int count = 0 ;
  for(int i=0;i<vp_size;++i) {
    send_buf[count++] = vp[i].first ;
    send_buf[count++] = vp[i].second ;
  }
  // release vp buffer to save some memory because we no longer need it
  vector<pair<int,int> >().swap(vp) ;
  // prepare recv buffer
  vector<int> recv_buf(total_recv_size) ;
  
  MPI_Alltoallv(&send_buf[0], &send_counts[0],
                &send_displs[0], MPI_INT,
                &recv_buf[0], &recv_counts[0],
                &recv_displs[0], MPI_INT, comm) ;
  // finally extract the data to fill the pair vector
  // release send_buf first to save some memory
  vector<int>().swap(send_buf) ;
  vp.resize(total_recv_size/2) ;
  count = 0 ;
  for(int i=0;i<total_recv_size;i+=2,count++)
    vp[count] = pair<int,int>(recv_buf[i], recv_buf[i+1]) ;
}





namespace Loci{
 void
  createEdgesParallel(fact_db &facts) {
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet faces = face2node.domain() ;

    // Loop over faces and create list of edges (with duplicates)
    vector<pair<Entity,Entity> > emap ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity e1 = face2node[*ei][i] ;
        Entity e2 = face2node[*ei][i+1] ;
        emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
      }
      Entity e1 = face2node[*ei][0] ;
      Entity e2 = face2node[*ei][sz-1] ;
      emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
    }

    // before we do the parallel sorting, we perform a check
    // to see if every process at least has one data element in
    // the "emap", if not, then the parallel sample sort would fail
    // and we pre-balance the "emap" on every process before the
    // sorting
    if(GLOBAL_OR(emap.empty())) {
      parallel_balance_pair_vector(emap, MPI_COMM_WORLD) ;
    }
    // Sort edges and remove duplicates
    sort(emap.begin(),emap.end()) ;
    vector<pair<Entity,Entity> >::iterator uend ;
    uend = unique(emap.begin(), emap.end()) ;
    emap.erase(uend, emap.end()) ;
    // then sort emap in parallel
    // but we check again to see if every process has at least one
    // element, if not, that means that the total element number is
    // less than the total number of processes, we split the communicator
    // so that only those do have elements would participate in the
    // parallel sample sorting
    if(GLOBAL_OR(emap.empty())) {
      MPI_Comm sub_comm ;
      int color = emap.empty() ;
      MPI_Comm_split(MPI_COMM_WORLD, color, MPI_rank, &sub_comm) ;
      if(!emap.empty())
        par_sort(emap, sub_comm) ;
      MPI_Comm_free(&sub_comm) ;
    } else {
      par_sort(emap, MPI_COMM_WORLD) ;
    }
    // remove duplicates again in the new sorted vector
    uend = unique(emap.begin(), emap.end()) ;
    emap.erase(uend, emap.end()) ;
#ifdef BOUNDARY_DUPLICATE_DETECT
    if(MPI_processes > 1) {
      // then we will need to remove duplicates along the boundaries
      // we send the first element in the vector to the left neighbor
      // processor (my_id - 1) and each processor compares its last
      // element with the received element. if they are the same,
      // then the processor will remove its last element

      // HOWEVER if the parallel sort was done using the sample sort
      // algorithm, then this step is not necessary. Because in the
      // sample sort, elements are partitioned to processors according
      // to sample splitters, it is therefore guaranteed that no
      // duplicates will be crossing the processor boundaries.
      int sendbuf[2] ;
      int recvbuf[2] ;
      if(!emap.empty()) {
        sendbuf[0] = emap[0].first ;
        sendbuf[1] = emap[0].second ;
      } else {
        // if there is no local data, we set the send buffer
        // to be the maximum integer so that we don't have problems
        // in the later comparing stage
        sendbuf[0] = std::numeric_limits<int>::max() ;
        sendbuf[1] = std::numeric_limits<int>::max() ;
      }
      MPI_Status status ;
      if(MPI_rank == 0) {
        // rank 0 only receives from 1, no sending needed
        MPI_Recv(recvbuf, 2, MPI_INT,
                 1/*source*/, 0/*msg tag*/,
                 MPI_COMM_WORLD, &status) ;
      } else if(MPI_rank == MPI_processes-1) {
        // the last processes only sends to the second last processes,
        // no receiving is needed
        MPI_Send(sendbuf, 2, MPI_INT,
                 MPI_rank-1/*dest*/, 0/*msg tag*/, MPI_COMM_WORLD) ;
      } else {
        // others will send to MPI_rank-1 and receive from MPI_rank+1
        MPI_Sendrecv(sendbuf, 2, MPI_INT, MPI_rank-1/*dest*/,0/*msg tag*/,
                     recvbuf, 2, MPI_INT, MPI_rank+1/*source*/,0/*tag*/,
                     MPI_COMM_WORLD, &status) ;
      }
      // then compare the results with last element in local emap
      if( (MPI_rank != MPI_processes-1) && (!emap.empty())){
        const pair<Entity,Entity>& last = emap.back() ;
        if( (recvbuf[0] == last.first) &&
            (recvbuf[1] == last.second)) {
          emap.pop_back() ;
        }
      }
    } // end if(MPI_Processes > 1)
#endif
    
    // Allocate entities for new edges
    int num_edges = emap.size() ;
    entitySet edges = facts.get_distributed_alloc(num_edges).first ;

    //added by Qiuhan, for a constraint edges
    constraint edges_tag;
    *edges_tag = edges;
    facts.create_fact("edges", edges_tag);
    
    // Copy edge nodes into a MapVec
    MapVec<2> edge ;
    edge.allocate(edges) ;
    vector<pair<Entity,Entity> >::iterator pi = emap.begin() ;
    for(entitySet::const_iterator ei=edges.begin();
        ei!=edges.end();++ei,++pi) {
      edge[*ei][0] = pi->first ;
      edge[*ei][1] = pi->second ;
    }

    // Add edge2node data structure to fact databse
    // facts.create_fact("edge2node",edge) ;

    // Now create face2edge data-structure
    // We need to create a lower node to edge mapping to facilitate the
    // searches.  First get map from edge to lower node
    Map el ; // Lower edge map
    el.allocate(edges) ;
    for(entitySet::const_iterator ei=edges.begin();
        ei!=edges.end();++ei,++pi) {
      el[*ei] = edge[*ei][0] ;
    }

    // Now invert this map to get nodes-> edges that have this as a first entry
    multiMap n2e ;
    // Get nodes
    // Get mapping from nodes to edges from lower numbered node
    
    // note inorder to use the distributed_inverseMap, we need
    // to provide a vector of entitySet partitions. for this 
    // case, it is NOT the node (pos.domain()) distribution,
    // instead it is the el Map image distribution
    entitySet el_image = el.image(el.domain()) ;
    vector<entitySet> el_image_partitions =
      gather_all_entitySet(el_image) ;
    distributed_inverseMap(n2e, el, el_image, edges, el_image_partitions) ;

    // Now create face2edge map with same size as face2node
    multiMap face2edge ;
    store<int> count ;
    count.allocate(faces) ;
    for(entitySet::const_iterator ei = faces.begin();
        ei!=faces.end();++ei) 
      count[*ei] = face2node[*ei].size() ;
    face2edge.allocate(count) ;

    // before computing the face2edge map, we will need to gather
    // necessary info among all processors since the edge map is
    // distributed across all the processors. we need to retrieve
    // those that are needed from other processors.

    // we will first need to figure out the set of edges we need
    // but are not on the local processor

    // but we need to access the n2e map in the counting and it
    // is possible that the local n2e map does not have enough
    // data we are looking for, therefore we need to expand it
    // first to include possible clone regions
    entitySet nodes_accessed ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        nodes_accessed += e1 ;
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      nodes_accessed += e1 ;
    }
    // we then expand the n2e map
    entitySet nodes_out_domain = nodes_accessed - n2e.domain() ;    
    n2e.setRep(MapRepP(n2e.Rep())->expand(nodes_out_domain,
                                          el_image_partitions)) ;
    // okay, then we are going to expand the edge map
    // first count all the edges we need
    entitySet edges_accessed ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        for(int j=0;j<n2e[e1].size();++j) {
          int e = n2e[e1][j] ;
          edges_accessed += e ;
        }
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      for(int j=0;j<n2e[e1].size();++j) {
        int e = n2e[e1][j] ;
        edges_accessed += e ;
      }
    }
    vector<entitySet> edge_partitions = gather_all_entitySet(edge.domain()) ;
    entitySet edges_out_domain = edges_accessed - edge.domain() ;
    // but since there is no expand method implemented for
    // MapVec at this time, we will just do a hack to convert
    // the MapVec to a multiMap to reuse the expand code.
    multiMap edge2 ;
    store<int> edge2_count ;
    entitySet edge_dom = edge.domain() ;
    edge2_count.allocate(edge_dom) ;
    for(entitySet::const_iterator ei=edge_dom.begin();
        ei!=edge_dom.end();++ei)
      edge2_count[*ei] = 2 ;
    edge2.allocate(edge2_count) ;
    for(entitySet::const_iterator ei=edge_dom.begin();
        ei!=edge_dom.end();++ei) {
      edge2[*ei][0] = edge[*ei][0] ;
      edge2[*ei][1] = edge[*ei][1] ;
    }
    edge2.setRep(MapRepP(edge2.Rep())->expand(edges_out_domain,
                                              edge_partitions)) ;
    // we are now ready for the face2edge map

    // Now loop over faces, for each face search for matching edge and
    // store in the new face2edge structure
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      // Loop over edges of the face
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        Entity e2 = max(t1,t2) ;
        face2edge[*ei][i] = -1 ;
        // search for matching edge
        for(int j=0;j<n2e[e1].size();++j) {
          int e = n2e[e1][j] ;
          if(edge2[e][0] == e1 && edge2[e][1] == e2) {
            face2edge[*ei][i] = e ;
            break ;
          }
        }
        if(face2edge[*ei][i] == -1)
          cerr << "ERROR: not able to find edge for face " << *ei << endl ;
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      Entity e2 = max(t1,t2) ;
      face2edge[*ei][sz-1] = -1 ;
      for(int j=0;j<n2e[e1].size();++j) {
        int e = n2e[e1][j] ;
        if(edge2[e][0] == e1 && edge2[e][1] == e2) {
          face2edge[*ei][sz-1] = e ;
          break ;
        }
      }
      if(face2edge[*ei][sz-1] == -1)
        cerr << "ERROR: not able to find edge for face " << *ei << endl ;
      
    }
    // Add face2edge to the fact database
    facts.create_fact("face2edge",face2edge) ;
    //sort edge2node according to fileNumbering

  
    if(MPI_processes > 1){    
      //create Map node_l2f
       Map node_l2f; 
      entitySet nodes;
      
      FORALL(edge.domain(), e){
        nodes += edge[e][0];
        nodes += edge[e][1];
      }ENDFORALL;
    
      
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      
      dMap g2f ;
      g2f = df->g2f.Rep() ;
       
      entitySet localNodes = nodes&init_ptn[MPI_rank] ;
      node_l2f.allocate(localNodes);
      FORALL(localNodes, d){
        node_l2f[d] = g2f[d];
      }ENDFORALL;
       
      entitySet out_of_dom = nodes - localNodes;
      vector<entitySet> tmp_ptn = gather_all_entitySet(localNodes);
      node_l2f.setRep(MapRepP(node_l2f.Rep())->expand(out_of_dom, tmp_ptn)) ;
     
      
      //end of create Map
       
      FORALL(edge.domain(), e){
        if(node_l2f[edge[e][0] ]> node_l2f[edge[e][1]]){
          std:: swap(edge[e][0], edge[e][1]);
        }
               
      }ENDFORALL;
    }
      
    

    multiMap edge3 ;
    store<int> edge3_count ;
    entitySet edge3_dom = edge.domain() ;
    edge3_count.allocate(edge3_dom) ;
    for(entitySet::const_iterator ei=edge3_dom.begin();
        ei!=edge3_dom.end();++ei)
      edge3_count[*ei] = 2 ;
    edge3.allocate(edge3_count) ;
    for(entitySet::const_iterator ei=edge3_dom.begin();
        ei!=edge3_dom.end();++ei) {
      edge3[*ei][0] = edge[*ei][0] ;
      edge3[*ei][1] = edge[*ei][1] ;
    }
    
    // Add edge3node data structure to fact databse
    facts.create_fact("edge2node",edge3) ;

  
    
    
 } // end of createEdgesPar
 





}




class fileNumRule: public pointwise_rule {
  store<int> fileNumX ;
public:
  fileNumRule() {
    name_store("fileNumber(X)",fileNumX) ;
    output("fileNumber(X)") ;
    constraint("X") ;
    disable_threading() ;
  }
  void compute(const sequence &seq) {
   if(Loci::MPI_processes == 1) {
    
     for(sequence::const_iterator si=seq.begin();si!= seq.end();++si){
       fileNumX[*si] = *si ;
       // if(*si < 0) cout << "negative file Number " << *si << endl; 
     }
     return;
   }
   fact_db::distribute_infoP df = Loci::exec_current_fact_db->get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    dMap g2f ;
    g2f = df->g2f.Rep() ;

    for(sequence::const_iterator si=seq.begin();si!= seq.end();++si){
      fileNumX[*si] = g2f[l2g[*si]] ;
      //if(fileNumX[*si] < 0) cout << "negative file Number " << fileNumX[*si] << endl; 
    }
  }
  
} ;

register_rule<fileNumRule> register_fileNumRule ;






// namespace Loci{
//   entitySet findBoundingSet(entitySet dom);
  
//   storeRepP  my_get_node_remap(entitySet nodes, entitySet faces, entitySet edges) {
    
//     entitySet dom = nodes + faces + edges;
//     if(MPI_processes == 1) {
//       int minNode = nodes.Min() ;
//       Map nm ;
     
//       nm.allocate(dom) ;
      
//       FORALL(nodes,nd) {
//         nm[nd] = nd - minNode +1;;
//        } ENDFORALL ;
      
//       FORALL(faces,nd) {
//         nm[nd] = nd;
//       } ENDFORALL ;
      
//       FORALL(edges,nd) {
//         nm[nd] = nd;
//       } ENDFORALL ;
      
//       return nm.Rep() ;
//     }
    
//     std::vector<entitySet> init_ptn = Loci::exec_current_fact_db->get_init_ptn() ;
//     fact_db::distribute_infoP df = Loci::exec_current_fact_db->get_distribute_info() ;
//     Map l2g ;
//     l2g = df->l2g.Rep() ;
//     dMap g2f ;
//     g2f = df->g2f.Rep() ;
    
//     entitySet gnodes = l2g.image(nodes&l2g.domain()) ;
//     entitySet gset = findBoundingSet(gnodes) ;
//     int minNode = gset.Min();
    
//     entitySet gelements =l2g.image(dom&l2g.domain()) ; 
//     Map newnum ;
//     newnum.allocate(dom) ;
    
//     // Expand g2f to include clone regions
//     entitySet out_of_dom = gelements - init_ptn[MPI_rank] ;
//     g2f.setRep(MapRepP(g2f.Rep())->expand(out_of_dom, init_ptn)) ;
    
//     FORALL(nodes,i) {
//       newnum[i] = g2f[l2g[i]] - minNode +1;
//     } ENDFORALL ;
    
//     FORALL(faces,i) {
//       newnum[i] = g2f[l2g[i]];
//     } ENDFORALL ;
//     FORALL(edges,i) {
//       newnum[i] = g2f[l2g[i]];
//     } ENDFORALL ;
//     return newnum.Rep() ;
//   }
// }













// // Note that this is a unit_rule, since this is the
// // only way we can have stores as input to a rule outputting blackboxes.
// class set_cell_remap_unit : public unit_rule {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//   blackbox<Loci::storeRepP > node_remap ;
// public:

//       // Define input and output.
//   set_cell_remap_unit() {
//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("node_remap", node_remap);
    
//     input("(lower, upper, boundary_map)->face2node->pos");
//     input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     output("node_remap");
//     constraint("geom_cells") ;
//     disable_threading() ;
//   }
  
//   // Do the set-up.
//   void compute(const sequence & seq) {
//    *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
//   }
  
// };
// register_rule<set_cell_remap_unit> register_set_cell_remap_unit ;

//  // Empty apply rule required by Loci. The data type and operator do not
//   // matter since nothing is done by this rule. Keep the same inputs and
//   // outputs as the unit rule, even though we don't have to.
// class set_cell_remap_apply : public apply_rule<blackbox<Map>,
//                              Loci::NullOp<Map> > {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//   blackbox<Loci::storeRepP > node_remap ;
// public:
  
//   // Define input and output.
//   set_cell_remap_apply() {
//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("node_remap", node_remap);
//     input("(lower, upper, boundary_map)->face2node->pos");
//     input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     output("node_remap");
//     constraint("geom_cells") ;
//     disable_threading() ;
//   }
  
//   // Do nothing.
//   void compute(const sequence & seq) {}
//   } ;

//   register_rule<set_cell_remap_apply> registerset_cell_remap_apply ;


// class set_interior_face_remap_unit : public unit_rule {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//    const_Map cl;
//   const_Map cr;
//   blackbox<Loci::storeRepP > node_remap ;
// public:

//       // Define input and output.
//   set_interior_face_remap_unit() {
//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("iface_remap", node_remap);
//     name_store("cl", cl);
//     name_store("cr", cr);
//     input("(cl, cr)->(lower, upper, boundary_map)->face2node->pos");
//     input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("iface_remap");
//     constraint("interior_faces") ;
//     disable_threading() ;
//   }
  
//   // Do the set-up.
//   void compute(const sequence & seq) {
//    *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
//   }
  
// };
// register_rule<set_interior_face_remap_unit> register_set_interior_face_remap_unit ;

//  // Empty apply rule required by Loci. The data type and operator do not
//   // matter since nothing is done by this rule. Keep the same inputs and
//   // outputs as the unit rule, even though we don't have to.
// class set_interior_face_remap_apply : public apply_rule<blackbox<Map>,
//                              Loci::NullOp<Map> > {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//    const_Map cl;
//   const_Map cr;
//   blackbox<Loci::storeRepP > node_remap ;
// public:
  
//   // Define input and output.
//   set_interior_face_remap_apply() {

//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("iface_remap", node_remap);
//     name_store("cl", cl);
//     name_store("cr", cr);
//     input("(cl, cr)->(lower, upper, boundary_map)->face2node->pos");
//     input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("iface_remap");
//     constraint("interior_faces") ;
//     disable_threading() ;
//   }
    
  
//   // Do nothing.
//   void compute(const sequence & seq) {}
//   } ;

// register_rule<set_interior_face_remap_apply> register_set_interior_face_remap_apply ;

// class set_boundary_face_remap_unit : public unit_rule {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//   const_Map cl;
//   blackbox<Loci::storeRepP > node_remap ;
// public:

//       // Define input and output.
//   set_boundary_face_remap_unit() {
//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("bface_remap", node_remap);
//     name_store("cl", cl);
  
//     input("cl->(lower, upper, boundary_map)->face2node->pos");
//     input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("bface_remap");
//     constraint("boundary_faces") ;
//     disable_threading() ;
//   }
  
//   // Do the set-up.
//   void compute(const sequence & seq) {
//    *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
//   }
  
// };
// register_rule<set_boundary_face_remap_unit> register_set_boundary_face_remap_unit ;

//  // Empty apply rule required by Loci. The data type and operator do not
//   // matter since nothing is done by this rule. Keep the same inputs and
//   // outputs as the unit rule, even though we don't have to.
// class set_boundary_face_remap_apply : public apply_rule<blackbox<Map>,
//                              Loci::NullOp<Map> > {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//    const_Map cl;
 
//   blackbox<Loci::storeRepP > node_remap ;
// public:
  
//   // Define input and output.
//   set_boundary_face_remap_apply() {

//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("bface_remap", node_remap);
//     name_store("cl", cl);
 
//     input("cl->(lower, upper, boundary_map)->face2node->pos");
//     input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("bface_remap");
//     constraint("boundary_faces") ;
//     disable_threading() ;
//   }
    
  
//   // Do nothing.
//   void compute(const sequence & seq) {}
//   } ;

// register_rule<set_boundary_face_remap_apply> register_set_boundary_face_remap_apply ;

// //for fl
// class set_face_remap_unit : public unit_rule {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//   const_Map cl;
//   blackbox<Loci::storeRepP > node_remap ;
// public:

//       // Define input and output.
//   set_face_remap_unit() {
//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("face_remap", node_remap);
//     name_store("cl", cl);
  
//     input("cl->(lower, upper, boundary_map)->face2node->pos");
//     input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("face_remap");
//     constraint("faces") ;
//     disable_threading() ;
//   }
  
//   // Do the set-up.
//   void compute(const sequence & seq) {
//    *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
//   }
  
// };
// register_rule<set_face_remap_unit> register_set_face_remap_unit ;

//  // Empty apply rule required by Loci. The data type and operator do not
//   // matter since nothing is done by this rule. Keep the same inputs and
//   // outputs as the unit rule, even though we don't have to.
// class set_face_remap_apply : public apply_rule<blackbox<Map>,
//                              Loci::NullOp<Map> > {
// private:
//   const_store<vect3d> pos;
//   const_multiMap upper;
//   const_multiMap lower;
//   const_multiMap boundary_map;
  
//   const_multiMap edge2node;
//   const_multiMap face2edge;
//   const_multiMap face2node;
//    const_Map cl;
 
//   blackbox<Loci::storeRepP > node_remap ;
// public:
  
//   // Define input and output.
//   set_face_remap_apply() {

//     name_store("pos", pos);
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("face2edge", face2edge);
//     name_store("edge2node", edge2node);
//     name_store("face_remap", node_remap);
//     name_store("cl", cl);
 
//     input("cl->(lower, upper, boundary_map)->face2node->pos");
//     input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
//     input("face2node->pos");
//     output("face_remap");
//     constraint("faces") ;
//     disable_threading() ;
//   }
    
  
//   // Do nothing.
//   void compute(const sequence & seq) {}
//   } ;

// register_rule<set_face_remap_apply> register_set_face_remap_apply ;

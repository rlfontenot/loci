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
#include <Tools/tools.h>
#include <map>
#include "pnn.h"
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
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

#include "dist_tools.h"
using std::cout ;


#include <gconstraint.h>
#include <gstore.h>
#include <gmap.h>
#include <gmultimap.h>
#include <gfact_db.h>
#include <gmapvec.h>
#include <distribute_long.h>



#define vect3d vector3d<double>
namespace Loci{
  
  void createLowerUpper(gfact_db &facts) {
    gConstraint faces, geom_cells,interior_faces,boundary_faces ;
    faces = facts.get_variable("faces") ;
    geom_cells = facts.get_variable("geom_cells") ;
    interior_faces = facts.get_variable("interior_faces") ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    gEntitySet bfaces = *boundary_faces ;
    gEntitySet ifaces = *interior_faces ;
 
    gStoreRepP pfacesP = facts.get_variable("periodicFaces") ;
    if(pfacesP != 0) {
      gConstraint periodicFaces ;
      periodicFaces = pfacesP ;
      bfaces -= *periodicFaces ;
      ifaces += *periodicFaces ;
    }
 
    gEntitySet global_interior_faces = g_all_collect_entitySet<gEntity>(ifaces) ;
    gEntitySet global_boundary_faces = g_all_collect_entitySet<gEntity>(bfaces) ;
 
    gMap cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    gEntitySet global_geom_cells ; 
    //init_ptn can not use key_ptn of gKeySpace, because key_ptn is in file numbering
    std::vector<gEntitySet> init_ptn = g_all_collect_vectors<gEntity>(*geom_cells) ;
    global_geom_cells = g_all_collect_entitySet<gEntity>(*geom_cells) ;
    gMultiMap lower,upper,boundary_map ;
    upper = cl.distributed_inverse(global_geom_cells, global_interior_faces, init_ptn);
    lower = cr.distributed_inverse(global_geom_cells, global_interior_faces, init_ptn);
    boundary_map = cl.distributed_inverse(global_geom_cells, global_boundary_faces, init_ptn);
        
    facts.create_fact("lower",lower, lower.get_domain_space(), lower.get_image_space()) ;
    facts.create_fact("upper",upper, upper.get_domain_space(), upper.get_image_space()) ;
    facts.create_fact("boundary_map",boundary_map,boundary_map.get_domain_space(), boundary_map.get_image_space()) ;
  }

  namespace{
    // this is a general routine that balances the pair vector
    // on each process, it redistributes the pair vector
    // among the processes such that each process holds
    // roughly the same number of elements, it maintains the
    // original global elements ordering in the redistribution
    void
      parallel_balance_pair_vector(vector<pair<gEntity,gEntity> >& vp,
                                   MPI_Comm comm) {
      int num_procs = 0 ;
      MPI_Comm_size(comm,&num_procs) ;

      // we still use an all-to-all personalized communication
      // algorithm to balance the element numbers on processes.
      // we pick (p-1) equally spaced element as the splitters
      // and then re-split the global vector sequence to balance
      // the number of elements on processes.

      gEntity vp_size = vp.size() ;
      gEntity global_vp_size = 0 ;
      MPI_Allreduce(&vp_size, &global_vp_size,
                    1, MPI_GENTITY_TYPE, MPI_SUM, comm) ;

      gEntity space = global_vp_size / num_procs ;
      // compute a global range for the elements on each process
      gEntity global_end = 0 ;
      MPI_Scan(&vp_size, &global_end, 1, MPI_GENTITY_TYPE, MPI_SUM, comm) ;
      gEntity global_start = global_end - vp_size ;

      vector<gEntity> splitters(num_procs) ;
      // splitters are just global index number
      splitters[0] = space ;
      for(int i=1;i<num_procs-1;++i)
        splitters[i] = splitters[i-1] + space ;
      splitters[num_procs-1] = global_vp_size ;

      // split and communicate the vector of particles
      vector<gEntity> send_counts(num_procs, 0) ;
      gEntity part_start = global_start ;
      for(int idx=0;idx<num_procs;++idx) {
        if(part_start == global_end)
          break ;
        if(splitters[idx] > part_start) {
          gEntity part_end ;
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

      vector<gEntity> send_displs(num_procs) ;
      send_displs[0] = 0 ;
      for(gEntity i=1;i<num_procs;++i)
        send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

      vector<gEntity> recv_counts(num_procs) ;
      MPI_Alltoall(&send_counts[0], 1, MPI_GENTITY_TYPE,
                   &recv_counts[0], 1, MPI_GENTITY_TYPE, comm) ;

      vector<gEntity> recv_displs(num_procs) ;
      recv_displs[0] = 0 ;
      for(gEntity i=1;i<num_procs;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

      gEntity total_recv_size = recv_displs[num_procs-1] +
        recv_counts[num_procs-1] ;

      // prepare send and recv buffer
      vector<gEntity> send_buf(vp_size*2) ;
      gEntity count = 0 ;
      for(gEntity i=0;i<vp_size;++i) {
        send_buf[count++] = vp[i].first ;
        send_buf[count++] = vp[i].second ;
      }
      // release vp buffer to save some memory because we no longer need it
      vector<pair<gEntity,gEntity> >().swap(vp) ;
      // prepare recv buffer
      vector<gEntity> recv_buf(total_recv_size) ;

      MPI_Alltoallv(&send_buf[0], &send_counts[0],
                    &send_displs[0], MPI_GENTITY_TYPE,
                    &recv_buf[0], &recv_counts[0],
                    &recv_displs[0], MPI_GENTITY_TYPE, comm) ;
      // finally extract the data to fill the pair vector
      // release send_buf first to save some memory
      vector<gEntity>().swap(send_buf) ;
      vp.resize(total_recv_size/2) ;
      count = 0 ;
      for(int i=0;i<total_recv_size;i+=2,count++)
        vp[count] = pair<gEntity,gEntity>(recv_buf[i], recv_buf[i+1]) ;
    }
  
    // a parallel sample sort for vector<pair<gEntity, gEntity> >
    // the passed in vector is the local SORTED data. NOTE:
    // the precondition to this routine is that the passed
    // in vector is sorted!!!
    // after sorting, this function puts the new sorted pairs
    // that are local to a processor in the data argument.
    void
      par_sort(vector<pair<gEntity,gEntity> >& data, MPI_Comm comm) {
      // first get the processor id and total number of processors
      int my_id, num_procs ;
      MPI_Comm_size(comm, &num_procs) ;
      MPI_Comm_rank(comm, &my_id) ;
      if(num_procs <= 1)
        return ;                  // single process, no need to proceed
      // get the number of local elements
      int local_size = data.size() ;
      // then select num_procs-1 equally spaced elements as splitters
      gEntity* splitters = new gEntity[num_procs] ;
      int even_space = local_size / (num_procs-1) ;
      int start_idx = even_space / 2 ;
      int space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = data[space_idx].first ;
      // gather the splitters to all processors as samples
      int sample_size = num_procs * (num_procs-1) ;
      gEntity* samples = new gEntity[sample_size] ;
      MPI_Allgather(splitters, num_procs-1, MPI_GENTITY_TYPE,
                    samples, num_procs-1, MPI_GENTITY_TYPE, comm) ;
      // now we've obtained all the samples, first we sort them
      sort(samples, samples+sample_size) ;
      // select new splitters in the sorted samples
      even_space = sample_size / (num_procs-1) ;
      start_idx = even_space / 2 ;
      space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = samples[space_idx] ;
      // the last one set as maximum possible gEntityeger
      splitters[num_procs-1] = std::numeric_limits<gEntity>::max() ;

      // now we can assign local elements to buckets (processors)
      // according to the new splitters. first we will compute
      // the size of each bucket and communicate them first
      int* scounts = new gEntity[num_procs] ;
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
      // but since one local element contains two gEntityegers (a pair of gEntity),
      // we will need to double the size
      for(int i=0;i<num_procs;++i)
        scounts[i] *= 2 ;
      // now we compute the sending displacement for each bucket
      int* sdispls = new gEntity[num_procs] ;
      sdispls[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        sdispls[i] = sdispls[i-1] + scounts[i-1] ;
      // communicate this information to all processors so that each will
      // know how many elements are expected from every other processor
      int* rcounts = new gEntity[num_procs] ;
      MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm) ;
      // then based on the received info. we will need to compute the
      // receive displacement
      int* rdispls = new gEntity[num_procs] ;
      rdispls[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        rdispls[i] = rdispls[i-1] + rcounts[i-1] ;
      // then we will need to pack the elements in local gEntityo
      // a buffer and communicate them
      gEntity* local_pairs = new gEntity[local_size*2] ;
      int count = 0 ;
      for(int i=0;i<local_size;++i) {
        local_pairs[count++] = data[i].first ;
        local_pairs[count++] = data[i].second ;
      }
      // then we allocate buffer for new local elements
      int new_local_size = rdispls[num_procs-1] + rcounts[num_procs-1] ;
      gEntity* sorted_pairs = new gEntity[new_local_size] ;
      // finally we communicate local_pairs to each processor
      MPI_Alltoallv(local_pairs, scounts, sdispls, MPI_GENTITY_TYPE,
                    sorted_pairs, rcounts, rdispls, MPI_GENTITY_TYPE, comm) ;
      // release buffers
      delete[] splitters ;
      delete[] samples ;
      delete[] scounts ;
      delete[] sdispls ;
      delete[] rcounts ;
      delete[] rdispls ;
      delete[] local_pairs ;
      // finally we unpack the buffer gEntityo a vector of pairs
      data.resize(new_local_size/2) ;
      int data_idx = 0 ;
      for(int i=0;i<new_local_size;i+=2,data_idx++)
        data[data_idx] = pair<gEntity,gEntity>(sorted_pairs[i],sorted_pairs[i+1]) ;
      // release the final buffer
      delete[] sorted_pairs ;
      // finally we sort the new local vector
      sort(data.begin(), data.end()) ;
    }
    void
      par_sort2(vector<pair<pair<gEntity,gEntity>, gEntity> >& data, MPI_Comm comm){
      // first get the processor id and total number of processors
      int my_id, num_procs ;
      MPI_Comm_size(comm, &num_procs) ;
      MPI_Comm_rank(comm, &my_id) ;
      if(num_procs <= 1)
        return ;                  // single process, no need to proceed
      // get the number of local elements
      int local_size = data.size() ;
      // then select num_procs-1 equally spaced elements as splitters
      pair<gEntity, gEntity> *splitters = new pair<gEntity, gEntity>[num_procs] ;
      int even_space = local_size / (num_procs-1) ;
      int start_idx = even_space / 2 ;
      int space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = data[space_idx].first ;
      // gather the splitters to all processors as samples
      gEntity sample_size = num_procs * (num_procs-1) ;
      pair<gEntity, gEntity>* samples = new pair<gEntity, gEntity>[sample_size] ;
      MPI_Allgather(splitters, (num_procs-1)*2, MPI_GENTITY_TYPE,
                    samples, (num_procs-1)*2, MPI_GENTITY_TYPE, comm) ;
      // now we've obtained all the samples, first we sort them
      sort(samples, samples+sample_size) ;
      // select new splitters in the sorted samples
      even_space = sample_size / (num_procs-1) ;
      start_idx = even_space / 2 ;
      space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = samples[space_idx] ;
      // the last one set as maximum possible gEntityeger
      gEntity maxnumber = std::numeric_limits<gEntity>::max();
      splitters[num_procs-1] =pair<gEntity, gEntity>(maxnumber, maxnumber);
                               
  
      // now we can assign local elements to buckets (processors)
      // according to the new splitters. first we will compute
      // the size of each bucket and communicate them first
      int *scounts = new gEntity[num_procs] ;
      for(gEntity i=0;i<num_procs;++i)
        scounts[i] = 0;
      { // using a block just to make the definition of "i" and "j" local
        gEntity i, j ;
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
      // but since one local element contains two gEntityegers (a pair of gEntity),
      // we will need to double the size
      for(int i=0;i<num_procs;++i)
        scounts[i] *= 3 ;
      // now we compute the sending displacement for each bucket
      int* sdispls = new gEntity[num_procs] ;
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
      // then we will need to pack the elements in local gEntityo
      // a buffer and communicate them
      gEntity* local_pairs = new gEntity[local_size*3] ;
      int count = 0 ;
      for(int i=0;i<local_size;++i) {
        local_pairs[count++] = data[i].first.first ;
        local_pairs[count++] = data[i].first.second;
        local_pairs[count++] = data[i].second;
      }
      // then we allocate buffer for new local elements
      int new_local_size = rdispls[num_procs-1] + rcounts[num_procs-1] ;
      gEntity* sorted_pairs = new gEntity[new_local_size] ;
      // finally we communicate local_pairs to each processor
      MPI_Alltoallv(local_pairs, scounts, sdispls, MPI_GENTITY_TYPE,
                    sorted_pairs, rcounts, rdispls, MPI_GENTITY_TYPE, comm) ;
      // release buffers
      delete[] splitters ;
      delete[] samples ;
      delete[] scounts ;
      delete[] sdispls ;
      delete[] rcounts ;
      delete[] rdispls ;
      delete[] local_pairs ;
      // finally we unpack the buffer gEntityo a vector of pairs
      data.resize(new_local_size/3) ;
      int data_idx = 0 ;
      for(int i=0;i<new_local_size;i+=3,data_idx++)
        data[data_idx] = pair<pair<gEntity,gEntity>, gEntity>(pair<gEntity, gEntity>(sorted_pairs[i],sorted_pairs[i+1]), sorted_pairs[i+2]) ;
      // release the final buffer
      delete[] sorted_pairs ;
      // finally we sort the new local vector
      sort(data.begin(), data.end()) ;
    }


    void
      parallel_balance_pair2_vector(vector<pair<pair<gEntity,gEntity>, gEntity> >& vp,
                                    MPI_Comm comm) {
      int num_procs = 0 ;
      MPI_Comm_size(comm,&num_procs) ;
  
      // we still use an all-to-all personalized communication
      // algorithm to balance the element numbers on processes.
      // we pick (p-1) equally spaced element as the splitters
      // and then re-split the global vector sequence to balance
      // the number of elements on processes.
  
      gEntity vp_size = vp.size() ;
      gEntity global_vp_size = 0 ;
      MPI_Allreduce(&vp_size, &global_vp_size,
                    1, MPI_GENTITY_TYPE, MPI_SUM, comm) ;
  
      gEntity space = global_vp_size / num_procs ;
      // compute a global range for the elements on each process
      gEntity global_end = 0 ;
      MPI_Scan(&vp_size, &global_end, 1, MPI_GENTITY_TYPE, MPI_SUM, comm) ;
      gEntity global_start = global_end - vp_size ;
  
      vector<gEntity> splitters(num_procs) ;
      // splitters are just global index number
      splitters[0] = space ;
      for(int i=1;i<num_procs-1;++i)
        splitters[i] = splitters[i-1] + space ;
      splitters[num_procs-1] = global_vp_size ;
  
      // split and communicate the vector of particles
      vector<gEntity> send_counts(num_procs, 0) ;
      gEntity part_start = global_start ;
      for(int idx=0;idx<num_procs;++idx) {
        if(part_start == global_end)
          break ;
        if(splitters[idx] > part_start) {
          gEntity part_end ;
          if(splitters[idx] < global_end)
            part_end = splitters[idx] ;
          else
            part_end = global_end ;
          send_counts[idx] = part_end - part_start ;
          part_start = part_end ;
        }
      }
  
      for(size_t i=0;i<send_counts.size();++i)
        send_counts[i] *= 3 ;
  
      vector<gEntity> send_displs(num_procs) ;
      send_displs[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
  
      vector<gEntity> recv_counts(num_procs) ;
      MPI_Alltoall(&send_counts[0], 1, MPI_GENTITY_TYPE,
                   &recv_counts[0], 1, MPI_GENTITY_TYPE, comm) ;
  
      vector<gEntity> recv_displs(num_procs) ;
      recv_displs[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
  
      gEntity total_recv_size = recv_displs[num_procs-1] +
        recv_counts[num_procs-1] ;
  
      // prepare send and recv buffer
      vector<gEntity> send_buf(vp_size*3) ;
      gEntity count = 0 ;
      for(gEntity i=0;i<vp_size;++i) {
        send_buf[count++] = vp[i].first.first ;
        send_buf[count++] = vp[i].first.second ;
        send_buf[count++] = vp[i].second;
      }
      // release vp buffer to save some memory because we no longer need it
      vector<pair<pair<gEntity,gEntity>, gEntity>  >().swap(vp) ;
      // prepare recv buffer
      vector<gEntity> recv_buf(total_recv_size) ;
  
      MPI_Alltoallv(&send_buf[0], &send_counts[0],
                    &send_displs[0], MPI_GENTITY_TYPE,
                    &recv_buf[0], &recv_counts[0],
                    &recv_displs[0], MPI_GENTITY_TYPE, comm) ;
      // finally extract the data to fill the pair vector
      // release send_buf first to save some memory
      vector<gEntity>().swap(send_buf) ;
      vp.resize(total_recv_size/3) ;
      count = 0 ;
      for(gEntity i=0;i<total_recv_size;i+=3,count++)
        vp[count] = pair<pair<gEntity,gEntity>, gEntity>(pair<gEntity, gEntity>(recv_buf[i], recv_buf[i+1]), recv_buf[i+2]) ;
    }
    // end of unnamed namespace
  }

  
  void
  createEdgesPar(gfact_db &facts) {
    gMultiMap face2node ;
    face2node = facts.get_variable("face2node") ;
   
    // Loop over faces and create list of edges (with duplicates)
    vector<pair<gEntity,gEntity> > emap ;
    fatal(!face2node.sorted());
    {
      gMultiMap::const_iterator previous = face2node.begin();
      gMultiMap::const_iterator itr = face2node.begin();
      gMultiMap::const_iterator start = face2node.begin();
      itr++;
      for(; itr != face2node.end(); itr++){ //start with the second element
        previous = itr;
        previous--; //previous is always the one in front of itr 
        if(itr->first == previous->first){ //if itr and previous belong to the same face, create an edge 
          gEntity e1 = previous->second ;
          gEntity e2 = itr->second ;
          emap.push_back(pair<gEntity,gEntity>(min(e1,e2),max(e1,e2))) ;
        }else{ //else connect the last node(previous) to the first node(start), and update start
          gEntity e1 = previous->second ;
          gEntity e2 = start->second ;
          emap.push_back(pair<gEntity,gEntity>(min(e1,e2),max(e1,e2))) ;
          start = itr;
        }
      }
      //the last edge
      previous = itr;
      previous--; 
      gEntity e1 = previous->second ;
      gEntity e2 = start->second ;
      emap.push_back(pair<gEntity,gEntity>(min(e1,e2),max(e1,e2))) ;
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
    vector<pair<gEntity,gEntity> >::iterator uend ;
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
      gEntity sendbuf[2] ;
      gEntity recvbuf[2] ;
      if(!emap.empty()) {
        sendbuf[0] = emap[0].first ;
        sendbuf[1] = emap[0].second ;
      } else {
        // if there is no local data, we set the send buffer
        // to be the maximum integer so that we don't have problems
        // in the later comparing stage
        sendbuf[0] = std::numeric_limits<gEntity>::max() ;
        sendbuf[1] = std::numeric_limits<gEntity>::max() ;
      }
      MPI_Status status ;
      if(MPI_rank == 0) {
        // rank 0 only receives from 1, no sending needed
        MPI_Recv(recvbuf, 2, MPI_GENTITY_TYPE,
                 1/*source*/, 0/*msg tag*/,
                 MPI_COMM_WORLD, &status) ;
      } else if(MPI_rank == MPI_processes-1) {
        // the last processes only sends to the second last processes,
        // no receiving is needed
        MPI_Send(sendbuf, 2, MPI_GENTITY_TYPE,
                 MPI_rank-1/*dest*/, 0/*msg tag*/, MPI_COMM_WORLD) ;
      } else {
        // others will send to MPI_rank-1 and receive from MPI_rank+1
        MPI_Sendrecv(sendbuf, 2, MPI_GENTITY_TYPE, MPI_rank-1/*dest*/,0/*msg tag*/,
                     recvbuf, 2, MPI_GENTITY_TYPE, MPI_rank+1/*source*/,0/*tag*/,
                     MPI_COMM_WORLD, &status) ;
      }
      // then compare the results with last element in local emap
      if( (MPI_rank != MPI_processes-1) && (!emap.empty())){
        const pair<gEntity,gEntity>& last = emap.back() ;
        if( (recvbuf[0] == last.first) &&
            (recvbuf[1] == last.second)) {
          emap.pop_back() ;
        }
      }
    } // end if(MPI_Processes > 1)
#endif
     
    //create EdgeSpace and register it
    string casename; //cheat here, should find it somewhere in facts
    gKeySpaceP edge_space = new gKeySpace();
    edge_space->register_space("EdgeSpace", casename); 
    //get key manager and generate the global keys for edges
    gKeyManagerP key_manager = facts.get_key_manager();
    int num_edges = emap.size() ;
    gEntitySet edges = edge_space->generate_key(key_manager,num_edges);
    
    //create constraint edges
    gConstraint edges_tag;
    *edges_tag = edges;
    facts.create_fact("edges", edges_tag, edge_space);
   
    // Copy edge nodes into a MapVec
    gMapVec<2> edge ;
    vector<pair<gEntity,gEntity> >::iterator pi = emap.begin() ; 
    for(gEntitySet::const_iterator ei=edges.begin();
        ei!=edges.end();++ei,++pi) {
      edge.insert(*ei, pi->first);
      edge.insert(*ei, pi->second);
    }
    {// add block here to make local variables dissepear at the end  
      // Now create face2edge data-structure
      // We need to create a lower node to edge mapping to facilitate the
      // searches.  First get map from edge to lower node
      gMap el ; // Lower edge map
      for(gMapVec<2>::const_iterator itr = edge.begin(); itr!= edge.end(); itr++){ 
        el.insert(itr->first, itr->second);
        itr++;
      }
    
    
      // Now invert this map to get nodes-> edges that have this as a first entry
      gMultiMap n2e ;
      // Get nodes
      // Get mapping from nodes to edges from lower numbered node
    
      // note inorder to use the distributed_inverseMap, we need
      // to provide a vector of gEntitySet partitions. for this 
      // case, it is NOT the node (pos.domain()) distribution,
      // instead it is the el Map image distribution
      gEntitySet el_image = el.image() ;
      vector<gEntitySet> el_image_partitions = g_all_collect_vectors<gEntity>(el_image);
      n2e = el.distributed_inverse(el_image, edges, el_image_partitions) ;
      el.clear();
      
      // Now create face2edge map with same size as face2node

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
      gEntitySet nodes_accessed ;
      {
        gMultiMap::const_iterator previous = face2node.begin();
        gMultiMap::const_iterator itr = face2node.begin();
        gMultiMap::const_iterator start = face2node.begin();
        itr++;
        for(; itr != face2node.end(); itr++){
          previous = itr;
          previous--;
          if(itr->first == previous->first){
            gEntity e1 = previous->second ;
            gEntity e2 = itr->second ;
            gEntity n1 = min(e1,e2) ;
            nodes_accessed += n1;
          }else{
            gEntity e1 = previous->second ;
            gEntity e2 = start->second ;
            gEntity n1 = min(e1,e2) ;
            nodes_accessed += n1;
            start = itr;
          }
        }
        //the last edge
        previous = itr;
        previous--; 
        gEntity e1 = previous->second ;
        gEntity e2 = start->second ;
        gEntity n1 = min(e1,e2) ;
        nodes_accessed += n1;
      }


    
      // we then expand the n2e map
      n2e = n2e.expand(nodes_accessed, el_image_partitions) ;
      // okay, then we are going to expand the edge map
      // first count all the edges we need
      gEntitySet edges_accessed = n2e.image() ;
      n2e.clear();
      
      vector<gEntitySet> edge_partitions =  g_all_collect_vectors<gEntity>(edges) ;
      gMapVec<2> expanded_edge;
      expanded_edge = edge.expand(edges_accessed, edge_partitions) ;
      // we are now ready for the face2edge map

      // Now loop over faces, for each face search for matching edge and
      // store in the new face2edge structure
      vector<f2e_comp> f2e_vec;
      {
        short ind = 0;
        gMultiMap::const_iterator previous = face2node.begin();
        gMultiMap::const_iterator itr = face2node.begin();
        gMultiMap::const_iterator start = face2node.begin();
        itr++;
        for(; itr != face2node.end(); itr++){
          previous = itr;
          previous--;
          if(itr->first == previous->first){
            gEntity e1 = previous->second ;
            gEntity e2 = itr->second ;
            gEntity n1 = min(e1,e2) ;
            gEntity n2 = max(e1,e2) ;
            f2e_vec.push_back(f2e_comp(itr->first, n1, n2, ind++));
          }else{
            gEntity e1 = previous->second ;
            gEntity e2 = start->second ;
            gEntity n1 = min(e1,e2) ;
            gEntity n2 = max(e1,e2) ;
            f2e_vec.push_back(f2e_comp(previous->first, n1, n2, ind));
            ind = 0;
            start = itr;
          }
        }
        //the last edge
        previous = itr;
        previous--; 
        gEntity e1 = previous->second ;
        gEntity e2 = start->second ;
        gEntity n1 = min(e1,e2) ;
        gEntity n2 = max(e1,e2) ;
        f2e_vec.push_back(f2e_comp(previous->first, n1, n2, ind)); 
      }
      sort(f2e_vec.begin(), f2e_vec.end(), f2e_field_sort_e);

      vector<e2n_comp> e2n_vec;
      for(gMapVec<2>::const_iterator itr = expanded_edge.begin(); itr != expanded_edge.end(); itr++){ 
        gEntity i1 = itr->first;
        gEntity i2 = itr->second;
        itr++;
        gEntity i3 = itr->second;
        e2n_vec.push_back(e2n_comp(i1, i2, i3));
      }
      sort(e2n_vec.begin(), e2n_vec.end(), e2n_field_sort_n);

      vector<ms_comp<gEntity> > new_vec;
      //equiJoin
      vector<f2e_comp>::const_iterator itr2 = f2e_vec.begin();
      for(vector<e2n_comp>::const_iterator itr1 = e2n_vec.begin(); itr1!= e2n_vec.end(); itr1++){
        while(itr2!= f2e_vec.end() && (itr2->n1 < itr1->n1 || (itr2->n1 == itr1->n1 &&  itr2->n2 < itr1->n2))) itr2++; 
        while(itr2!= f2e_vec.end() && (itr2->n1 == itr1->n1 &&  itr2->n2 == itr1->n2)) {
          new_vec.push_back(ms_comp<gEntity>(itr2->f, itr1->e, itr2->ind));
          itr2++;
        }
      }
      vector<f2e_comp>().swap(f2e_vec);//free up memory
      vector<e2n_comp>().swap(e2n_vec);//free up memory
      sort(new_vec.begin(), new_vec.end(),field_sort_dom<gEntity>);

      gMultiMap face2edge;
      for(vector<ms_comp<gEntity> >::const_iterator itr = new_vec.begin(); itr != new_vec.end(); itr++){
        face2edge.insert(itr->dom, itr->img);
      }
      // Add face2edge to the fact database
      facts.create_fact("face2edge",face2edge) ;
    }

   
    //sort edge2node according to fileNumbering
    if(MPI_processes > 1){    
      gKeySpaceP node_space = gKeySpace::get_space("NodeSpace", casename);
      gMapVec<2> edge_file;
      gMap node_g2f;
      node_g2f = node_space->get_g2f_map();
      
      edge_file = edge.recompose(node_space->get_g2f_map());
      
 
      //make inside edge and edge_file, edge direction is from node with lower file numbe to node with higher file number
      gMapVec<2>::iterator itr2 = edge.begin();            
      for(gMapVec<2>::iterator itr = edge_file.begin(); itr != edge_file.end(); itr++, itr2++){
        gMapVec<2>::iterator next = itr;
        next++;
        gMapVec<2>::iterator next2 = itr2;
        next2++;
        if(itr->second > next->second){
          std::swap(itr->second, next->second);
          std::swap(itr2->second, next2->second);
        }
        itr++, itr2++;
      }
    
      //then update keyspace so that the file number of edges is consistent with the file number of nodes

      //give each edge a file number
      vector<pair<pair<gEntity, gEntity> , gEntity> > edge_file_vec(num_edges);
      int eindex = 0;
      for(gMapVec<2>::iterator itr = edge_file.begin(); itr != edge_file.end(); itr++){
        gMapVec<2>::iterator next = itr;
        next++;
        edge_file_vec[eindex++] = pair<pair<gEntity, gEntity>, gEntity>(pair<gEntity, gEntity>(itr->second,next->second), itr->first);
        itr++;
      }
     
      if(GLOBAL_OR(edge_file_vec.empty())) {
        parallel_balance_pair2_vector(edge_file_vec, MPI_COMM_WORLD) ;
      }
      // Sort edges and remove duplicates
      //unless parallel_balance_pair2_vector create duplicates, there should be no duplicates
      sort(edge_file_vec.begin(),edge_file_vec.end()) ;

      vector<pair<pair<Entity,Entity>, Entity> >::iterator uend2 ;
      uend2 = unique(edge_file_vec.begin(), edge_file_vec.end()) ;
      edge_file_vec.erase(uend2, edge_file_vec.end()) ;
           
      // then sort emap in parallel
      // but we check again to see if every process has at least one
      // element, if not, that means that the total element number is
      // less than the total number of processes, we split the communicator
      // so that only those do have elements would participate in the
      // parallel sample sorting
      if(GLOBAL_OR(edge_file_vec.empty())) {
        MPI_Comm sub_comm ;
        int color = edge_file_vec.empty() ;
        MPI_Comm_split(MPI_COMM_WORLD, color, MPI_rank, &sub_comm) ;
        if(!edge_file_vec.empty())
          par_sort2(edge_file_vec, sub_comm) ;
        MPI_Comm_free(&sub_comm) ;
      } else {
        par_sort2(edge_file_vec, MPI_COMM_WORLD) ;
      }
      // remove duplicates again in the new sorted vector
      uend2 = unique(edge_file_vec.begin(), edge_file_vec.end()) ;
      edge_file_vec.erase(uend2, edge_file_vec.end()) ;
     
      
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
        gEntity sendbuf[3] ;
        gEntity recvbuf[3] ;
        if(!edge_file_vec.empty()) {
          sendbuf[0] = edge_file_vec[0].first.first ;
          sendbuf[1] = edge_file_vec[0].first.second ;
          sendbuf[2] = edge_file_vec[0].second;
        } else {
          // if there is no local data, we set the send buffer
          // to be the maximum integer so that we don't have problems
          // in the later comparing stage
          sendbuf[0] = std::numeric_limits<int>::max() ;
          sendbuf[1] = std::numeric_limits<int>::max() ;
          sendbuf[2] = std::numeric_limits<int>::max() ;
        }
        MPI_Status status ;
        if(MPI_rank == 0) {
          // rank 0 only receives from 1, no sending needed
          MPI_Recv(recvbuf, 3, MPI_GENTITY_TYPE,
                   1/*source*/, 0/*msg tag*/,
                   MPI_COMM_WORLD, &status) ;
        } else if(MPI_rank == MPI_processes-1) {
          // the last processes only sends to the second last processes,
          // no receiving is needed
          MPI_Send(sendbuf, 3, MPI_GENTITY,
                   MPI_rank-1/*dest*/, 0/*msg tag*/, MPI_COMM_WORLD) ;
        } else {
          // others will send to MPI_rank-1 and receive from MPI_rank+1
          MPI_Sendrecv(sendbuf, 3, MPI_GENTITY_TYPE, MPI_rank-1/*dest*/,0/*msg tag*/,
                       recvbuf, 3, MPI_GENTITY_TYPE, MPI_rank+1/*source*/,0/*tag*/,
                       MPI_COMM_WORLD, &status) ;
        }
        // then compare the results with last element in local emap
        if( (MPI_rank != MPI_processes-1) && (!edge_file_vec.empty())){
          const pair<pair<gEntity,gEntity>, gEntity>& last = edge_file_vec.back() ;
          if( (recvbuf[0] == last.first.first) &&
              (recvbuf[1] == last.first.second)&&
              (recvbuf[2] == last.secon)) {
            edge_file_vec.pop_back() ;
          }
        }
      } // end if(MPI_Processes > 1)
#endif
      
    
      int local_num_edge = edge_file_vec.size();
      vector<int> edge_sizes(MPI_processes);
      MPI_Allgather(&local_num_edge,1,MPI_INT,&edge_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    
      gEntity file_num_offset = 0;
      for(int i = 0; i < MPI_rank; i++){
        file_num_offset += edge_sizes[i];
      }
    
      gMap edge_global2file;
      edge_global2file.reserve(edge_file_vec.size());
      gEntity index = file_num_offset;
      for(int i = 0; i < local_num_edge; i++){
        edge_global2file.insert(edge_file_vec[i].second, index);
        index++;
      }
      edge_global2file.local_sort();
      //expand the map and update keyspace
      {
        gEntitySet domain_global = edge_global2file.domain();
        vector<gEntitySet> ptn = g_all_collect_vectors<gEntity>(domain_global);
        edge_space->set_g2f_map(edge_global2file.expand(edges, ptn));
      } 
    }//end of if(MPI_processes > 1)

    // Add edge2node to fact databse
    facts.create_fact("edge2node",edge) ;
   
  } // end of createEdgesPar
  
  int classify_cell(const gEntitySet& faces, const_gMultiMap &face2node) {
   
    int num_triangles = 0 ;
    int num_quads = 0 ;
    int num_others = 0 ;
    gEntity triangle_nodes[3][2] ;
    GFORALL(faces, fc){
      std::pair<const_gMultiMap::const_iterator, const_gMultiMap::const_iterator> r = face2node.range(fc);
      int count = distance(r.first, r.second);
      if(count == 3) {
        if(num_triangles < 2) {
          const_gMultiMap::const_iterator itr = r.first;
          triangle_nodes[0][num_triangles] = itr->second ; itr++;
          triangle_nodes[1][num_triangles] = itr->second ; itr++; 
          triangle_nodes[2][num_triangles] = itr->second ; itr++;
        }
        num_triangles++ ;
      } else if(count == 4)
        num_quads++ ;
      else
        num_others++ ;
    }ENDGFORALL;

   
    bool prism_test = false ;

    if((num_triangles == 2) && (num_quads == 3) && (num_others == 0)) {
      prism_test = true ;
      for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
          if(triangle_nodes[i][0] == triangle_nodes[j][1])
            prism_test = false ;
    }

    
    bool hex_test = false ;
    if( (num_triangles == 0) && (num_quads == 6) && (num_others == 0)) {
      gEntitySet::const_iterator fi = faces.begin();
      const gEntity ef = *fi ;// the first face
      std::pair<const_gMultiMap::const_iterator, const_gMultiMap::const_iterator> r1 = face2node.range(ef);
      int count = 0 ;
      fi++;
      for(;fi!=faces.end(); fi++) {//check the other 5 faces
        gEntity ej = *fi ;
        bool find = false;
        std::pair<const_gMultiMap::const_iterator, const_gMultiMap::const_iterator> r2 = face2node.range(ej);
        for(const_gMultiMap::const_iterator ef_itr = r1.first; ef_itr != r1.second; ef_itr++){
          for(const_gMultiMap::const_iterator ej_itr = r2.first; ej_itr != r2.second; ej_itr++){   
            if(ef_itr->second == ej_itr->second){
              find = true ;
              break;
            }
          }
          if(find) break;
        }
        if(find)
          count++ ;
      }
      
      if(count == 4)
        hex_test = true ;
    }
    // new classification code
    if( (num_triangles == 4) && (num_quads == 0) && (num_others == 0)) {
      return 0 ;
    } else if( hex_test ) {
      return 1 ;
    } else if( prism_test ) {
      return 2 ;
    } else if( (num_triangles == 4) && (num_quads == 1) && (num_others == 0)) {
      return 3 ;
    }
    return 4 ;
    
  }

}

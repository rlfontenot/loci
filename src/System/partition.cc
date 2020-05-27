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

#include <partition.h>
#include <distribute_io.h>
#include <gparameter.h>
#include <fact_db.h>
#include <gconstraint.h>
#ifdef HAS_MALLINFO
#include <malloc.h>
#endif

#include <Tools/tools.h>
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <parmetis.h>
#include <LociGridReaders.h> //real_t is defined here

#if REALTYPEWIDTH == 32
typedef float metisreal_t ;
#else
typedef double metisreal_t ;
#endif

namespace Loci{
  extern bool load_cell_weights ;
  extern string cell_weight_file ;
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
  

  enum mappingDirection{INWARD, OUTWARD};
  namespace {
    inline bool fieldSort2(const std::pair<gEntity,gEntity> &p1,
                           const std::pair<gEntity,gEntity> &p2) {
      return p1.second < p2.second ;
    }
  }

  //Given send_ptn, retuen recv_ptn
  //or given recv_ptn, return send_ptn
  vector<gEntitySet> transposePtn(const vector<gEntitySet> &ptn, MPI_Comm comm) {
    
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
    gEntity *send_store = new gEntity[size_send] ;
    gEntity *recv_store = new gEntity[size_recv] ;
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

    MPI_Datatype MPI_T_type = MPI_traits<gEntity>::get_MPI_type() ;
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_T_type,
		  recv_store, &recv_sz[0], recv_displacement, MPI_T_type,
		  comm) ;

    vector<gEntitySet> ptn_t(MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<recv_sz[i]/2;++j) {
        gEntity i1 = recv_store[recv_displacement[i]+j*2]  ;
        gEntity i2 = recv_store[recv_displacement[i]+j*2+1] ;
        ptn_t[i] += gInterval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return ptn_t ;
  }

  
  
  //assume input is sorted
  void get_selfmap(const gMultiMap& input,
                   const vector<gEntitySet>&dom_ptn,
                   gMultiMap& result){

    gMultiMap selfmap;
    gMultiMap::const_iterator first = input.begin();
    while(first != input.end()){
      gMultiMap::const_iterator last = first;
      while(last != input.end() && last->first == first->first) last++;
      gEntitySet image;
      for(gMultiMap::const_iterator itr = first; itr != last; itr++)image += itr->second;
      GFORALL(image, ii){
        GFORALL(image, jj){
          if(ii!=jj) selfmap.insert(ii, jj);
        }ENDGFORALL;
      }ENDGFORALL;
      first = last;
    }
    selfmap.local_sort();
    gStoreRepP rep = selfmap.split_redistribute(dom_ptn);
    result.setRep(rep);
  }
  
  //merge cl, cr into a gMultiMap
  //For boundary faces, only merge cl
  gStoreRepP merge_cl_cr( fact_db& facts){
    gMap cl, cr;
    gMultiMap result;
    cl = facts.get_gfact("cl");
    cr = facts.get_gfact("cr");
    fatal(cl.domain() != cr.domain()); 
    //all boundary surfaces
    gConstraint bcSurf;
    bcSurf = facts.get_gfact("bcSurf");
    gEntitySet refSet = g_all_collect_entitySet<gEntity>(*bcSurf);
   
    for( gMap::const_iterator itr1 = cl.begin(); itr1 != cl.end(); itr1++){
      result.insert(itr1->first, itr1->second);
    }
    for( gMap::const_iterator itr2 = cr.begin(); itr2 != cr.end(); itr2++){
      if(!refSet.inSet(itr2->second)) {
        result.insert(itr2->first, itr2->second);
      }
    }
    result.local_sort();
    result.set_domain_space(cl.get_domain_space());
    result.set_image_space(cl.get_image_space());
    return result.Rep();
  }
  


  //given the selfmap and the initial partition
  //return the owner(process) of entity in the local domain
  void newMetisPartition_raw( gMultiMap& selfmap,
                              const vector<gEntitySet> &init_ptn,//of the space to partition
                              vector<idx_t>& part
                              ) { //assignment of each entity to process  
    size_t map_size = init_ptn[Loci::MPI_rank].size() ;
    part.resize(map_size) ;
    vector<idx_t> xadj(map_size+1, 0) ;
    idx_t edgecut ;
    vector<idx_t> vdist(Loci::MPI_processes + 1) ;
    gEntity cmin = init_ptn[0].Min();
    for(int i = 0; i < Loci::MPI_processes; i++) {
      cmin = min(init_ptn[i].Min(), cmin);
    }
    
    // check init_ptn to be consistent with parmetis partitioning
    for(int i=0;i<Loci::MPI_processes;++i) {
      if(init_ptn[i].size() == 0 ||
         init_ptn[i].size() != size_t(init_ptn[i].Max()-init_ptn[i].Min()+1)) {
        cerr << "invalid local cell set, p=" << i
             << ", init_ptn=" << init_ptn[i] << endl ;
      }
      if(i>0 && init_ptn[i-1].Max()+1 != init_ptn[i].Min()) {
        cerr << "gap between processors in cell numbering" << endl ;
      }
    }
         
    //turn multimap into two vectors: adjncy and xadj(offset in adjncy)
    edgecut = 0 ;
    size_t tot = selfmap.size(); 
    vector<idx_t> adjncy(tot) ;
  
    idx_t count = 0 ;
    size_t ind = 0;
    xadj[ind] = count;
    gMultiMap::const_iterator itr = selfmap.begin();
    gEntity previous_value = itr->first;
    for(; itr != selfmap.end(); itr++) {
      adjncy[count] = itr->second-cmin ;
      if(itr->first != previous_value){
        xadj[++ind] = count;
        previous_value = itr->first;
      }
      count++;
    }
    xadj[++ind] = count;

    selfmap.setRep(gMultiMap().Rep()) ;// Free up memory from multiMap
       
    //for each process, the start index
    vdist[0] = 0 ;
    for(int i = 1; i <= Loci::MPI_processes; ++i)
      vdist[i] = vdist[i-1] + init_ptn[i-1].size() ;
    idx_t top = vdist[Loci::MPI_processes] ;
    
    bool trouble = false ;
    for(size_t i=0;i<tot;++i)
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
        hid_t file_id = hdf5OpenFile(cell_weight_file.c_str(),
                                     H5F_ACC_RDONLY, H5P_DEFAULT) ;
        if(file_id < 0) {
          std::cerr << "...file reading failed..., Aborting" << std::endl ;
          Loci::Abort() ;
        }
        
        // read
        gEntitySet dom = init_ptn[Loci::MPI_rank] ;
	// Needs fixing
        store<int> cell_weights ;

        readContainerRAW(file_id,"cell weight", cell_weights.Rep(),
                         MPI_COMM_WORLD) ;

	gEntitySet cwdom ;
	entitySet cwdomb = cell_weights.domain() ;
	size_t nivls = cwdomb.num_intervals() ;
	for(size_t i=0;i<nivls;++i)
	  cwdom += gInterval(cwdomb[i].first,cwdomb[i].second) ;
	
        if(cwdom != init_ptn[Loci::MPI_rank]) {
          cerr << "cell weights partition inconsistent!" << endl ;
          Loci::Abort() ;
        }
        
        Loci::hdf5CloseFile(file_id) ;
       
        // compute necessary ParMETIS data-structure
        wgtflag = 2 ;           // weights on the vertices only
        idx_t ncon = 2 ;          // number of weights per vertex
        size_t tpwgts_len = ncon*nparts ;
        vector<metisreal_t> tpwgts(tpwgts_len) ;

        for(size_t i=0;i<tpwgts_len;++i)
          tpwgts[i] = 1.0 / double(nparts) ;
        
        vector<metisreal_t> ubvec(ncon) ;
        for(idx_t i=0;i<ncon;++i)
          ubvec[i] = 1.05 ;     // as recommended by the ParMETIS manual

        // now construct the vertex weights
        vector<idx_t> vwgt(ncon*map_size) ;
        size_t cnt = 0 ;
        for(gEntitySet::const_iterator
              ei=init_ptn[Loci::MPI_rank].begin();
            ei!=init_ptn[Loci::MPI_rank].end();++ei,cnt+=ncon) {
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
        size_t tpwgts_len = ncon*nparts ;
        vector<metisreal_t> tpwgts(tpwgts_len) ;
        for(size_t i=0;i<tpwgts_len;++i)
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
      size_t tpwgts_len = ncon*nparts ;
      vector<metisreal_t> tpwgts(tpwgts_len) ;
      for(size_t i=0;i<tpwgts_len;++i)
        tpwgts[i] = 1.0 / double(nparts) ;
      
      metisreal_t ubvec = 1.05 ;
      wgtflag = 0 ;
      ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],NULL,NULL,
                           &wgtflag,&numflag,&ncon,&nparts,
                           &tpwgts[0],&ubvec,&options,&edgecut,
                           &part[0],&mc) ;
      
   
    }
    
    Loci::debugout << " Parmetis Edge cut   " <<  edgecut << endl ;
    return;
  }
  

  inline bool sortFirstORB(const pair<gEntity,pair<int,int> > &p1,
                           const pair<gEntity,pair<int,int> > &p2) {
    return p1.first < p2.first ;
  }

  //send_ptn: send_ptn of space1
  //f2n: and the mapping from space1 to space2
  //output: the vector of pair<entity in space2, pair<process it belongs to, 1> >
  //both input sorted, output not sorted, and duplication exists
  vector<pair<gEntity,pair<int,int> > > get_scratchPad(const vector<gEntitySet>& send_ptn, const gMultiMap& f2n){
    gStore<int> procmap;
    for(unsigned int i=0;i<send_ptn.size();++i) { 
      GFORALL(send_ptn[i],fc) {
        procmap.insert(fc, i);
      } ENDGFORALL ;
    }
    procmap.local_sort();
    gEntitySet fdom = f2n.domain();
    fatal(procmap.size() != fdom.size());
    vector<pair<gEntity,pair<int,int> > >   scratchPad(f2n.size()) ;
    
    int spad_ind = 0;
    gMultiMap::const_iterator itr2 = f2n.begin();
    gStore<int>::const_iterator itr1 = procmap.begin();
    GFORALL(fdom, fi){
      int proc = itr1->second;
      itr1++;
      while(itr2 != f2n.end() && itr2->first < fi)itr2++;
      while(itr2 != f2n.end() && itr2->first == fi){
        scratchPad[spad_ind++] = std::pair<gEntity,pair<int,int> > (itr2->second, std::pair<int, int>(proc, 1));
        itr2++;
      }
      
    }ENDGFORALL;     
    return scratchPad;
  }
  
  //promap: for each entity in space1, the process it belongs to
  //f2n: and the mapping from space1 to space2
  //output: the vector of pair<entity in space2, pair<process it belongs to, 1> >
  //both input sorted, output not sorted, and duplication exists
  vector<pair<gEntity,pair<int,int> > >  get_scratchPad(const vector<int>& procmap, const gMultiMap& f2n){
    gEntitySet fdom = f2n.domain();
    fatal(fdom.size() != procmap.size());
    vector<pair<gEntity,pair<int,int> > >   scratchPad(f2n.size()) ; 
    gMultiMap::const_iterator itr2 = f2n.begin();
    int procmap_ind = 0;
    int spad_ind = 0;
    GFORALL(fdom, fi){
      int proc = procmap[procmap_ind++];
      while(itr2 != f2n.end() && itr2->first < fi)itr2++;
      while(itr2 != f2n.end() && itr2->first == fi){
        scratchPad[spad_ind++] = std::pair<gEntity,pair<int,int> > (itr2->second, std::pair<int, int>(proc, 1));
        itr2++;
      }
    }ENDGFORALL;
    return scratchPad;    
  }


  //scratchPad: the vector of pair<entity in a space, pair<process it belongs to, 1> >
  //ptn: the original key partition of the space
  //out_ptn: the resulted send_ptn of local domain
  void assignOwner(vector<pair<gEntity,pair<int,int> > > &scratchPad,
                   const vector<gEntitySet>& ptn,
                   vector<gEntitySet> &out_ptn,
                   MPI_Comm comm) {

    std::sort(scratchPad.begin(),scratchPad.end()) ;
    
    //After processing, for each entity, how many times it fall onto to a process
    // Note, we are assuming scratch pad is non-zero size
    vector<pair<gEntity,pair<int,int> > >::iterator i1,i2 ;
    if(scratchPad.size() > 0) {
      i1 = scratchPad.begin() ;
      i2 = i1 +1;
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
    vector<pair<gEntity,pair<int,int> > > nsplits(MPI_processes-1) ;
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
      gEntity nd = scratchPad[x].first ;
      out_ptn[pmax] += nd ;
      x = y ;
    }
    // Check to see if any sets are left out.
    gEntitySet assigned ;
    for(int i=0;i<MPI_processes;++i)
      assigned += out_ptn[i] ;
    gEntitySet unassigned = ptn[MPI_rank] - assigned ;
    out_ptn[MPI_rank] += unassigned ; // allocate unassigned entities to
    // original processor
  }


 
  //this function returns if var_name is in_var of space
  //i.e., if the image space of var_name is space
  //var_name should be the name of a map
  bool is_in_var(const string& var_name, fact_db &facts, gKeySpaceP space){
    gStoreRepP st = facts.get_gvariable(var_name);
    gMapRepP mp = gMapRepP(st);
    if(mp==0){
      cerr<<"ERROR: " << var_name << " is not a map  in function  is_in_var()" << endl;
      return false;
    }
    gKeySpaceP ptr = mp->get_image_space();
    if(space == ptr) return true;
    return false;
  }
              
  void affinity_partition(fact_db &facts, const vector<int>& procmap,parStrategy& s){
    bool need_reverse2 = !is_in_var(s.map2, facts, s.space2);
   
    //from space1 to space2
    gMultiMap inward_map;
    if(s.map2 == "cl"||s.map2== "cr"){
      //merge cl and cr, for boundary face, only cl is inserted
      inward_map = merge_cl_cr(facts);
    }else{
      inward_map =facts.get_gfact(s.map2);
    }
   
    //if needed, reverse map
    
    if(need_reverse2){
      inward_map.setRep(inward_map.distributed_inverse((inward_map.get_image_space())->get_keys(),
                                                       (inward_map.get_domain_space())->get_keys(),
                                                       (inward_map.get_image_space())-> get_key_ptn()));
    }
       
    //equiJoin operation
    std::vector<pair<gEntity, pair<int,int> > >  scratchPad;
    scratchPad = get_scratchPad(procmap, inward_map);

   
    vector<gEntitySet> old_ptn = s.space2->get_key_ptn();
    MPI_Comm comm =  s.space2->get_mpi_comm();
    vector<gEntitySet> new_ptn2( s.space2->get_np());
    assignOwner(scratchPad,old_ptn, new_ptn2, comm) ;
    s.space2->set_send_recv_ptn( new_ptn2);
    
    gKeySpaceP from_space;
    if(s.from_space1){
      from_space = s.space1;
    }else{
      from_space = s.space2;
    }
       
    bool need_reverse3 = !is_in_var(s.map3, facts, s.space3);
    
   
    if(s.map3 == "cl" || s.map3 == "cr"){
      //merge cl and cr, for boundary face, only cl is inserted
      inward_map = merge_cl_cr(facts);
    }else{
      inward_map =facts.get_gfact(s.map3);
    }
    
    if(need_reverse3){
      inward_map.setRep(inward_map.distributed_inverse((inward_map.get_image_space())->get_keys(),
                                                       (inward_map.get_domain_space())->get_keys(),
                                                       (inward_map.get_image_space())-> get_key_ptn()));
    }
        
    if(s.from_space1){
      scratchPad = get_scratchPad(procmap, inward_map); 
    }else{
      scratchPad = get_scratchPad(new_ptn2, inward_map);
    }
    {
      vector<gEntitySet> old_ptn = s.space3->get_key_ptn();
      MPI_Comm comm =  s.space3->get_mpi_comm();
      vector<gEntitySet> new_ptn3( s.space3->get_np());
      assignOwner(scratchPad,old_ptn, new_ptn3, comm) ;
      s.space3->set_send_recv_ptn( new_ptn3);
    }
    
    return ;
  }

     
  //return entity2process vector     
  void primary_partition_metis(fact_db& facts, vector<int>& e2p, parStrategy& s){
    int np = s.space1->get_np();
    vector<gEntitySet> init_ptn = s.space1->get_key_ptn(); 
    bool need_reverse1 = !is_in_var(s.map1, facts, s.space1);
    gMultiMap inward_map;
    if(s.map1 == "cl" || s.map1 == "cr"){//merge cl and cr
      inward_map = merge_cl_cr(facts);
    }else{
      inward_map =facts.get_gfact(s.map1);
    }
    //if needed, reverse map
    if(need_reverse1){
      inward_map.setRep(inward_map.distributed_inverse((inward_map.get_image_space())->get_keys(),
                                                       (inward_map.get_domain_space())->get_keys(),
                                                       init_ptn));
    }
    

    gMultiMap self_map;
    get_selfmap(inward_map, init_ptn, self_map);

    inward_map.clear();
    newMetisPartition_raw(self_map, init_ptn, e2p);
    

    gEntitySet local_keys = s.space1->get_my_keys();
    gEntitySet::const_iterator ei = local_keys.begin();
    vector<gEntitySet> ptn(np);
    for(unsigned int i  = 0; i < e2p.size(); i++){
      ptn[e2p[i]] += *ei;
      ei++;
    }
    s.space1->set_send_recv_ptn(ptn);
  }

 
          
  //return entity2process vector     
  void primary_partition_orb(fact_db& facts, vector<int>& e2p,parStrategy& s ){
    int np = s.space1->get_np();
    MPI_Comm comm = s.space1->get_mpi_comm();
    gEntitySet domain = s.space1->get_my_keys();
    gStore<vector3d<real_t> > pos;
    pos = facts.get_gfact("pos");
    vector<vector3d<float> > center;
    if(s.map1=="pos"){
      //copy pos into center
      center.resize(pos.size());
      size_t ind = 0;
      for(gStore<vector3d<real_t> >::const_iterator itr1 = pos.begin(); itr1 != pos.end(); itr1++){
        center[ind++] = vector3d<float>(itr1->second.x, itr1->second.y, itr1->second.z);
      }
    }else if(s.map1=="facecenter"){
      gStoreRepP face2node;
      face2node = facts.get_gfact("face2node");
      gMultiStore<vector3d<real_t> > fpos; 
      fpos = pos.recompose(face2node, comm);
      gStore<vector3d<real_t> > fcenter;
      fcenter =  fpos.get_simple_center(fpos.domain());
      for(gStore<vector3d<real_t> >::const_iterator itr = fcenter.begin();
          itr != fcenter.end(); itr++){
        center.push_back(vector3d<float>((itr->second).x, (itr->second).y, (itr->second).z));
      }
    }else if(s.map1=="cellcenter"){
      gMultiMap inward_map;
      inward_map = merge_cl_cr(facts); //face2cell
      inward_map.setRep(inward_map.distributed_inverse((inward_map.get_image_space())->get_keys(),
                                                       (inward_map.get_domain_space())->get_keys(),
                                                       (inward_map.get_image_space())-> get_key_ptn()));//cell2face
      gMultiMap face2node;
      face2node = facts.get_gfact("face2node");
      gStoreRepP cell2node = inward_map.recompose(face2node, comm); //cell2node
      gMultiStore<vector3d<real_t> > cpos; 
      cpos = pos.recompose(cell2node, comm);
      gStore<vector3d<real_t> > ccenter;
      ccenter = cpos.get_simple_center(cpos.domain());
      for(gStore<vector3d<real_t> >::const_iterator itr = ccenter.begin();
          itr != ccenter.end(); itr++){
        center.push_back(vector3d<float>((itr->second).x, (itr->second).y, (itr->second).z));
      }
    }

   
    // perform ORB partition
    ORBPartition(center,e2p,comm) ;
    vector<gEntitySet> ptn(np);
    gEntitySet local_keys = s.space1->get_my_keys();
    gEntitySet::const_iterator ei = local_keys.begin();
    for(unsigned int i  = 0; i < e2p.size(); i++){
      ptn[e2p[i]] += *ei;
      ei++;
    }
    s.space1->set_send_recv_ptn(ptn);
  }
}

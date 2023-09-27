//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include "loci_globs.h"
#include "sched_tools.h"
#include "dist_tools.h"
#include "param_rule.h"
#include "thread.h"
#include <Tools/except.h>
#include <constraint.h>
#include <new>
#include "comp_tools.h"
#include <gpurep.h>
#include "gpuMap.h"

using std::bad_alloc ;
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ios ;
using std::ofstream ;
using std::istream ;
using std::ostream ;

///////////////////////////////////
#include <sstream>
#include <algorithm>
#include <sys/time.h> // for gettimeofday function
///////////////////////////////////

namespace Loci {
  std::vector<GPUstoreAllocateInfo> GPUstoreAllocateData ;
  std::vector<int> GPUstoreAllocateFreeList ;

  int getGPUStoreAllocateID() {
    // allocate slot in storeAllocateData
    int id = GPUstoreAllocateData.size() ;
    if(!GPUstoreAllocateFreeList.empty()) {
      id = GPUstoreAllocateFreeList.back() ;
      GPUstoreAllocateFreeList.pop_back() ;
    } else {
      GPUstoreAllocateData.push_back(GPUstoreAllocateInfo()) ;
    }
    GPUstoreAllocateData[id].alloc_ptr1 = 0 ;
    GPUstoreAllocateData[id].alloc_ptr2 = 0 ;
    GPUstoreAllocateData[id].base_ptr = 0 ;
    GPUstoreAllocateData[id].base_offset = 0 ;
    GPUstoreAllocateData[id].size = 0 ;
    GPUstoreAllocateData[id].allocated_size = 0 ;
    GPUstoreAllocateData[id].allocated = true ;
    GPUstoreAllocateData[id].allocset = EMPTY ;
    return id ;
  }
    
  void releaseGPUStoreAllocateID(int id) {
    GPUstoreAllocateData[id].alloc_ptr1 = 0 ;
    GPUstoreAllocateData[id].alloc_ptr2 = 0 ;
    GPUstoreAllocateData[id].base_ptr = 0 ;
    GPUstoreAllocateData[id].base_offset = 0 ;
    GPUstoreAllocateData[id].size = 0 ;
    GPUstoreAllocateData[id].allocated_size = 0 ;
    GPUstoreAllocateData[id].allocated = false ;
    GPUstoreAllocateData[id].allocset = EMPTY ;
    GPUstoreAllocateFreeList.push_back(id) ;
  }

  using std::pair ;
  using std::make_pair ;
  
  void gpuMapRepI::allocate(const entitySet &ptn) {
    if(alloc_id < 0)
      alloc_id = getGPUStoreAllocateID() ;

    GPUstoreAllocateData[alloc_id].template allocBasic<Entity>(ptn,1) ;
    store_domain = GPUstoreAllocateData[alloc_id].allocset ;
    base_ptr = (((Entity *) GPUstoreAllocateData[alloc_id].base_ptr) -
		GPUstoreAllocateData[alloc_id].base_offset) ;
      
    dispatch_notify() ;
    return ;
  }

  gpuMapRepI::~gpuMapRepI() {
    if(alloc_id>=0) {
      GPUstoreAllocateData[alloc_id].template release<Entity>() ;
      releaseGPUStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
  }

  storeRep *gpuMapRepI::new_store(const entitySet &p) const {
    return new gpuMapRepI(p)  ;
  }
  storeRep *gpuMapRepI::new_store(const entitySet &p, const int* count) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a Map " << endl ;
    return sp ;
  }
  storeRepP gpuMapRepI::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    pair<entitySet,entitySet> mappimage = preimage(rm.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = dm.image(newdomain) ;
    Map s ;
    s.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(s.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(dm,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(rm,mapimage) ;
    return s.Rep() ;
  }

  storeRepP gpuMapRepI::remap(const dMap &m) const {
    cerr << "Map shouldn't use remap!" << endl ;
    return MapRemap(m,m) ;
  }

  void gpuMapRepI::compose(const dMap &m, const entitySet &context) {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = m[base_ptr[i]] ;
    } ENDFORALL ;
  }

  void gpuMapRepI::copy(storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }

  void gpuMapRepI::gather(const dMap &m, storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ; 
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  void gpuMapRepI::scatter(const dMap &m,storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);

    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }
 
  int gpuMapRepI::pack_size(const entitySet &e) {
    fatal((e - domain()) != EMPTY);
    int size ;
    size = sizeof(Entity) * e.size() ;
    return(size) ;
  }
  int gpuMapRepI::estimated_pack_size(const entitySet &e) {
   
    return e.size()*sizeof(Entity) ;
  }
  int gpuMapRepI::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    int size = sizeof(Entity) * packed.size() ;
    return size ;
  }
  
  void gpuMapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      MPI_Pack( &base_ptr[begin], t, MPI_INT, outbuf, outcount, 
                &position, MPI_COMM_WORLD) ;
    }
  }
  
  void gpuMapRepI::pack(void *outbuf, int &position,
                     int &outcount, const entitySet &eset, const Map& remap) 
  {
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      int* img = new int[t] ;
      for(int k=0;k<t;++k)
        img[k] = remap[base_ptr[begin+k]] ;
      MPI_Pack(img, t, MPI_INT, outbuf, outcount, &position, MPI_COMM_WORLD) ;
      delete[] img ;
    }
  }
  
  void gpuMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) {

    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                      1 , MPI_INT, MPI_COMM_WORLD) ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                    t, MPI_INT, MPI_COMM_WORLD) ;
      }
    }
  }

  void gpuMapRepI::unpack(void *inbuf, int &position,
                       int &insize, const sequence &seq, const dMap& remap) {

    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                      1 , MPI_INT, MPI_COMM_WORLD) ;
        // remap
        for(Loci::int_type indx=seq[i].first;indx!=stop-1;--indx)
          base_ptr[indx] = remap[base_ptr[indx]] ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                    t, MPI_INT, MPI_COMM_WORLD) ;
        // remap
        for(int k=0;k<t;++k)
          base_ptr[indx+k] = remap[base_ptr[indx+k]] ;
      }
    }
  }

  entitySet gpuMapRepI::domain() const {
    return store_domain ;
  }

  entitySet gpuMapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(size_t i=0;i<d.num_intervals();++i)
        codomain += image_section(base_ptr+d[i].first,
                                  base_ptr+(d[i].second+1)) ;
    } else {
      std::vector<int> img(d.size()) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(size_t i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          *ins = base_ptr[j] ;
          ++ins ;
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
      
    return codomain ;
  }

  pair<entitySet,entitySet>
  gpuMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domain ;
    FORALL(store_domain,i) {
      if(codomain.inSet(base_ptr[i]))
        domain += i ;
    } ENDFORALL ;
    return make_pair(domain,domain) ;
  }
  
  storeRepP gpuMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
    storeRepP sp ;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(MPI_processes), send_clone(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = out_of_dom & ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      std::sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      std::sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(store_domain.inSet(*vi))
	  (map_entities[i])[*vi] = base_ptr[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	hm[recv_map[j]] = recv_map[j+1];
      }
    }
    Map dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    
    entitySet tm = store_domain + out_of_dom ;
    dm.allocate(tm) ;
    for(ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      dm[*ei] = base_ptr[*ei] ;
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) 
      dm[hmi->first] = hmi->second ;
    sp = dm.Rep() ;
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return sp ;
  }

  storeRepP gpuMapRepI::freeze() {
    return getRep() ;
  }
  
  storeRepP gpuMapRepI::thaw() {
    
    dMap dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    FORALL(store_domain,i) {
      dm[i] = base_ptr[i] ;
    } ENDFORALL ;
    return dm.Rep() ;
  }
  storeRepP gpuMapRepI::get_map() {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = 1 ;
    } ENDFORALL ;
    multiMap result ;
    result.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      result.begin(i)[0] = base_ptr[i] ;
    } ENDFORALL ;
    return result.Rep() ;
  }
    
  std::ostream &gpuMapRepI::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << std::endl ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  std::istream &gpuMapRepI::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      s >> base_ptr[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }
  DatatypeP gpuMapRepI::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }

  frame_info gpuMapRepI::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
   }

  void gpuMapRepI::copyFrom(const storeRepP &p, entitySet set) {
    int setivals = set.num_intervals() ;
    Map m ;
    m.setRep(p) ;
    Entity *gpu_base_ptr = get_base_ptr() ;
    for(int i=0;i<setivals;++i) {
      int start = set[i].first ;
      int end = set[i].second ;
      int sz = end-start+1 ;
      cudaError_t err = cudaMemcpy(gpu_base_ptr+start,&m[start],sizeof(int)*sz,
			       cudaMemcpyHostToDevice) ;
      if(err!= cudaSuccess) {
	cerr << "cudaMemcpy failed in gpuMapRepI::copyFrom" << endl ;
	Loci::Abort() ;
      }
    }
  }

  store_type gpuMapRepI::RepType() const  {
    return GPUMAP ;
  }
  
  void gpuMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset){
    warn(true) ; 
  } 

#ifdef H5_HAVE_PARALLEL 
  void gpuMapRepI::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset, hid_t xfer_plist_id){
    warn(true) ; 
  } 
#endif
  void gpuMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const{
    warn(true) ;
  } 

#ifdef H5_HAVE_PARALLEL 
  void gpuMapRepI::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset, hid_t xfer_plist_id) const{
    warn(true) ;
  } 
#endif  
  gpuMap::~gpuMap() {}

  void gpuMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }    

  const_gpuMap::~const_gpuMap() {}

  void const_gpuMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  store_instance::instance_type const_gpuMap::access() const
  { return READ_ONLY ; }
    
  namespace {
    rule create_rule(variable sv, variable tv, string qualifier) {
      ostringstream oss ;
      oss << "source(" << sv << ')' ;
      oss << ",target(" << tv << ')' ;
      oss << ",qualifier(" << qualifier << ')' ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
  }

  inline variable makeGPUVAR(variable v) {
    variable::info vinfo = v.get_info() ;
    string vname = string("__GPU__") + vinfo.name ;
    vinfo.name = vname ;
    return variable(vinfo) ;
  }
  rule_db rename_gpu_containers(fact_db  &facts,
				const rule_db &rdb)  {

    rule_db gpu_rdb ;

    variableSet gpuInputs ;
    variableSet gpuOutputs ;
    variableSet gpuMaps ;
    variableSet inputs ;
    variableSet outputs ;
    ruleSet rset = rdb.all_rules() ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end();++rsi) {
      rule_implP rp = rsi->get_rule_implP() ;
      if(rp == 0) {
	// Not a rule with an implementation so no GPU container
	// types
	gpu_rdb.add_rule(*rsi) ;
        continue ;
      }
      variableSet targets = rsi->targets() ;
      variableSet sources = rsi->sources() ;
      bool hasGPUVar = false ;
      bool hasCPUVar = false ;
      map<variable,variable> vm ;
      for(variableSet::const_iterator vsi = targets.begin(); vsi !=
	    targets.end(); ++vsi) {
	storeRepP sp = rp->get_store(*vsi) ;
	if(isGPU(sp)) {
	    gpuOutputs += *vsi ;
	  hasGPUVar = true ;
	  vm[*vsi] = makeGPUVAR(*vsi) ;
	} else {
	  outputs += *vsi ;
	  hasCPUVar = true ;
	  vm[*vsi] = *vsi ;
	}
      }
      for(variableSet::const_iterator vsi = sources.begin(); vsi !=
	    sources.end(); ++vsi) {
	storeRepP sp = rp->get_store(*vsi) ;
	if(isGPU(sp)) {
	  if(isMAP(sp))
	    gpuMaps += *vsi ;
	  else
	    gpuInputs += *vsi ;
	  hasGPUVar = true ;
	  vm[*vsi] = makeGPUVAR(*vsi) ;
	} else {
	  inputs += *vsi ;
	  hasCPUVar = true ;
	  vm[*vsi] = *vsi ;
	}
      }
      if(hasGPUVar) {
	if(hasCPUVar) {
	  cerr << "WARNING: rule " << *rsi << " contains both cpu and gpu containers! ---------" << endl ;
	}
	rp->rename_vars(vm) ;
	gpu_rdb.add_rule(rule(rp)) ;
      } else {
	gpu_rdb.add_rule(*rsi) ;
      }
      
    }
   

    cout << "gpuMaps = " << gpuMaps << endl ;
    //    cout << "inputs = " << inputs << endl ;
    //    cout << "outputs = " << outputs << endl ;
    //    cout << "gpuInputs = " << gpuInputs << endl ;
    //    cout << "gpuOutputs = " << gpuOutputs << endl ;

    variableSet cpu2gpu = gpuInputs ;
    cpu2gpu -= gpuOutputs ;
    for(variableSet::const_iterator vsi = cpu2gpu.begin(); vsi !=
	  cpu2gpu.end(); ++vsi) {
      gpu_rdb.add_rule(create_rule(*vsi,makeGPUVAR(*vsi),"cpu2gpu")) ;
    }
    for(variableSet::const_iterator vsi = gpuMaps.begin(); vsi !=
	  gpuMaps.end(); ++vsi) {
      gpu_rdb.add_rule(create_rule(*vsi,makeGPUVAR(*vsi),"map2gpu")) ;
    }
    variableSet gpu2cpu = gpuOutputs ;
    gpuOutputs -= gpuInputs ;
    for(variableSet::const_iterator vsi = gpu2cpu.begin(); vsi !=
	  gpu2cpu.end(); ++vsi) {
      gpu_rdb.add_rule(create_rule(makeGPUVAR(*vsi),*vsi,"gpu2cpu")) ;
    }    

    variableSet overlap = cpu2gpu ;
    overlap &= gpu2cpu ;
    if(overlap != EMPTY) {
      cerr << "warning, loops formed in interactions between gpu and cpu kernels" << endl ;
      cerr << "offending variables is " << overlap << endl ;
    }
    return gpu_rdb ;
  }

  void gpu2cpu_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    //    existential_rule_analysis(r,facts, scheds) ;
    variable vin = *r.sources().begin() ;
    variable vout = *r.targets().begin() ;

    entitySet dom = scheds.variable_existence(vin) ;
    scheds.set_existential_info(vout,r,dom) ;
  }

  void gpu2cpu_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    variable vin = *r.sources().begin() ;
    variable vout = *r.targets().begin() ;
    entitySet exec_seq = scheds.get_variable_request(r,vout) ;
    scheds.variable_request(vin,exec_seq) ;
    //entitySet exec_seq = process_rule_requests(r,facts, scheds) ;
    scheds.update_exec_seq(r, exec_seq);
  }

  executeP gpu2cpu_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(r) ;
    variable vgpu = *r.sources().begin() ;
    storeRepP p = facts.get_variable(vgpu)->getRep() ;
    gpuRepP gp = gpuRepP(p) ; 
    variable vcpu = *r.targets().begin() ;
    storeRepP cp = facts.get_variable(vcpu) ;
    executeP execute = executeP(new execute_gpu2cpu_copy(r,gp,cp,exec_seq)) ;
    return execute;
  }

  void cpu2gpu_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    variable vin = *r.sources().begin() ;
    variable vout = *r.targets().begin() ;

    entitySet dom = scheds.variable_existence(vin) ;
    scheds.set_existential_info(vout,r,dom) ;
    //    existential_rule_analysis(r,facts, scheds) ;
  }

  void cpu2gpu_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    variable vin = *r.sources().begin() ;
    variable vout = *r.targets().begin() ;
    entitySet exec_seq = scheds.get_variable_request(r,vout) ;
    scheds.variable_request(vin,exec_seq) ;
    //    entitySet exec_seq = process_rule_requests(r,facts, scheds) ;
    scheds.update_exec_seq(r, exec_seq);
  }

  executeP cpu2gpu_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(r) ;
    variable vgpu = *r.targets().begin() ;
    storeRepP p = facts.get_variable(vgpu)->getRep() ;
    gpuRepP gp = gpuRepP(p) ; 
    variable vcpu = *r.sources().begin() ;
    storeRepP cp = facts.get_variable(vcpu);
    
    executeP execute = executeP(new execute_cpu2gpu_copy(r,gp,cp,exec_seq)) ;
    return execute;
  }

  void map2gpu_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    //    existential_rule_analysis(r,facts, scheds) ;
    scheds.synonym_variable(*r.sources().begin(),*r.targets().begin()) ;
  }

  void map2gpu_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    //    entitySet exec_seq = process_rule_requests(r,facts, scheds) ;
    //    scheds.update_exec_seq(r, exec_seq);
  }

  executeP map2gpu_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(r) ;
    variable vgpu = *r.targets().begin() ;
    MapRepP p = MapRepP(facts.get_variable(vgpu)->getRep()) ;
    gpuMapRepP gp = gpuMapRepP(p) ; 
    variable vcpu = *r.sources().begin() ;
    storeRepP sp = facts.get_variable(vcpu) ;
    exec_seq = sp->domain() ;
    executeP execute = executeP(new execute_map2gpu_copy(r,gp,sp,exec_seq)) ;
    return execute;
  }
  
  void execute_gpu2cpu_copy::execute(fact_db &facts, sched_db &scheds) {
    gpuvar->copyTo(cpuvar,copyset) ;
  }

  void execute_gpu2cpu_copy::Print(ostream &s) const {
    printIndent(s) ;
    s << r << " over sequence " ;
    if(verbose || copyset.num_intervals() < 4) {
      s << copyset << endl ;
    } else {
      s << "[ ... ], l=" << copyset.size() << endl ;
    }
  }

  void execute_gpu2cpu_copy::dataCollate(collectData &data_collector) const {
    //    ostringstream oss ;
    //    oss << "rule: "<<rule_tag ;
    //
    //    data_collector.accumulateTime(timer,EXEC_COMPUTATION,osstr()) ;
  }

  void execute_cpu2gpu_copy::execute(fact_db &facts, sched_db &scheds) {
    gpuvar->copyFrom(cpuvar,copyset) ;
  }

  void execute_cpu2gpu_copy::Print(ostream &s) const {
    printIndent(s) ;
    s << r << " over sequence " ;
    if(verbose || copyset.num_intervals() < 4) {
      s << copyset << endl ;
    } else {
      s << "[ ... ], l=" << copyset.size() << endl ;
    }
  }

  void execute_cpu2gpu_copy::dataCollate(collectData &data_collector) const {
    //    ostringstream oss ;
    //    oss << "rule: "<<rule_tag ;
    //
    //    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  void execute_map2gpu_copy::execute(fact_db &facts, sched_db &scheds) {
    cerr << "copy map r=" << r << " set=" << copyset << endl ;
    gpuvar->allocate(copyset) ;
    gpuvar->copyFrom(cpuvar,copyset) ;
  }

  void execute_map2gpu_copy::Print(ostream &s) const {
    printIndent(s) ;
    s << r << " over sequence " ;
    if(verbose || copyset.num_intervals() < 4) {
      s << copyset << endl ;
    } else {
      s << "[ ... ], l=" << copyset.size() << endl ;
    }
  }

  void execute_map2gpu_copy::dataCollate(collectData &data_collector) const {
    //    ostringstream oss ;
    //    oss << "rule: "<<rule_tag ;
    //
    //    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }
  

}

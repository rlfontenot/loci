#include <Map.h>
#include <DMultiMap.h>
#include <multiMap.h>
#include <Tools/stream.h>
#include <set> 
namespace Loci 
{

  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;
  storeRepP dmultiMapRepI::thaw() {
    return getRep() ;
  }
  storeRepP dmultiMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
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
      sort(copy[i].begin(), copy[i].end()) ;
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
      sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector<HASH_MAP(int, std::vector<int>) > map_entities(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      for(HASH_MAP(int, std::vector<int>)::iterator hi = map_entities[i].begin(); hi != map_entities[i].end(); ++hi)
	send_count[i] += hi->second.size() ; 
    }
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += send_count[i] ;
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(HASH_MAP(int, std::vector<int> )::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second.size() ;
	++size_send ;
	for(std::vector<int>::const_iterator vi = miv->second.begin(); vi != miv->second.end(); ++vi) { 
	  send_map[size_send] = *vi ;
	  ++size_send ;
	}
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
    HASH_MAP(int, std::vector<int> ) hm ;
    std::vector<int> ss ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	int count = recv_map[j+1] ;
	if(count)
	  for(int k = 0; k < count; ++k)
	    hm[recv_map[j]].push_back(recv_map[j+k+2]);
	else
	  hm[recv_map[j]] = ss ;
	j += count + 1 ;
      }
    }
    std::vector<int> tmp_vec ;
    for(HASH_MAP(int, std::vector<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
      if(hmi->second.size()) 
	for(std::vector<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	  attrib_data[hmi->first].push_back(*si) ;
      else
	attrib_data[hmi->first] = tmp_vec ;
    dmultiMap dmul ;
    for(HASH_MAP(int, std::vector<int> )::const_iterator hi = attrib_data.begin(); hi != attrib_data.end(); ++hi)
      dmul[hi->first] = hi->second ;
    storeRepP sp = dmul.Rep() ;
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
  //------------------------------------------------------------------
  
  void dmultiMapRepI::allocate(const entitySet &eset) 
  {

    entitySet redundant, newSet;
    entitySet :: const_iterator  ci;

    redundant = domain() -  eset;
    newSet    = eset - domain();

    for( ci = redundant.begin(); ci != redundant.end(); ++ci)
         attrib_data.erase(*ci);

    std::vector<int>   emptyVec;
    for( ci = newSet.begin(); ci != newSet.end(); ++ci)
         attrib_data[*ci] = emptyVec;

    dispatch_notify() ;
  }
  
  //------------------------------------------------------------------
  
  void dmultiMapRepI::allocate(const store<int> &sizes) 
  {
    entitySet ptn = sizes.domain() ;
    entitySet :: const_iterator  ci;
    for( ci = ptn.begin(); ci != ptn.end(); ++ci) {
      std::vector<int>   newVec(sizes[*ci]) ;
      attrib_data[*ci] = newVec;
    }
    dispatch_notify() ;
  }
  
  //**************************************************************************/

  dmultiMapRepI::~dmultiMapRepI() 
  {
    attrib_data.clear();
  }
  
  //**************************************************************************/
  
  storeRep *dmultiMapRepI::new_store(const entitySet &p) const 
  {
    return new dmultiMapRepI()  ;
  }
  storeRep *dmultiMapRepI::new_store(const entitySet &p, const int* cnt) const 
  {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new dmultiMapRepI(count)  ;
  }
  //**************************************************************************/
  
  storeRepP dmultiMapRepI::remap(const dMap &m) const {
    dmultiMap s ;
    
    //-------------------------------------------------------------------------
    // Select only those entities from domain of "m" (which is input) which 
    // are part of multimap. Since, we have mapping defined only on those
    // entities.
    //-------------------------------------------------------------------------
    
    entitySet newdomain = m.domain() & domain() ;
    
    //-------------------------------------------------------------------------
    // Get the preimage of entire map. We will get two entity set. The
    // first is the intersection, and the second  is union of entities.
    //-------------------------------------------------------------------------
    pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    
    //
    // Get all entities which are both in preimage and newdomain
    //
    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    
    multiMap   newmap; 
    newmap = MapRepP(s.Rep())->get_map() ;
    return newmap.Rep() ;  
   
    // return s.Rep() ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::compose(const dMap &m, const entitySet &context) 
  {
    vector<int>    vec;
    HASH_MAP(int,vector<int> ) ::const_iterator ci;
    
    //-------------------------------------------------------------------------
    // All the entities in the context should be present in the domain. ie. A->B
    // should be valid.
    //-------------------------------------------------------------------------
    fatal((context-domain()) != EMPTY) ;
    
    //-------------------------------------------------------------------------
    // All in the entities in B should be part of Map. Also. i.e. B->C should
    // be valid.
    //-------------------------------------------------------------------------
    fatal((image(context)-m.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      ci = attrib_data.find(i);
      if( ci != attrib_data.end() ) {
        vec =  ci->second;
        for(size_t j = 0; j < vec.size(); j++) {
          attrib_data[i][j] =   m[vec[j]];
        }
      }
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiMap s(st) ;
    vector<int>    newVec;
    
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      for(size_t j = 0; j < newVec.size(); j++)
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::gather(const dMap &m, storeRepP &st, const entitySet  &context) 
  {
    const_dmultiMap s(st) ;
    vector<int>    newVec;
    
    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[m[i]];
      for(size_t j = 0; j < newVec.size(); j++) 
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::scatter(const dMap &m, storeRepP &st, const entitySet  &context) 
  {
    const_dmultiMap s(st) ;
    vector<int>    newVec;
    FORALL(context,i) {
      attrib_data[m[i]].clear();
      newVec  =   s[i];
      for(size_t j = 0; j < newVec.size(); j++) 
        attrib_data[m[i]].push_back( newVec[j] );
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  int dmultiMapRepI::pack_size(const  entitySet &e ) 
  {
    int size = 0 ;
    FORALL(e,i) {
      size  +=  attrib_data[i].size();
    } ENDFORALL ;
    
    return( size*sizeof(int) + e.size()*sizeof(int) ) ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::pack( void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
     int vsize;
    entitySet :: const_iterator ci;
    std::vector<int>   newVec;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize  = attrib_data[*ci].size();
      MPI_Pack( &vsize, 1, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
      MPI_Pack( &attrib_data[*ci][0], vsize, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
    }

  }
  
  
  //**************************************************************************/

  void dmultiMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) 
  {
    sequence:: const_iterator ci;
    std::vector<int>   newVec;

    int vsize;
    for( ci = seq.begin(); ci != seq.end(); ++ci){
         MPI_Unpack( inbuf, insize, &position, &vsize,
                     1, MPI_INT, MPI_COMM_WORLD) ;
         newVec.resize(vsize);
         MPI_Unpack( inbuf, insize, &position, &newVec[0],
                      vsize, MPI_INT, MPI_COMM_WORLD) ;
         attrib_data[*ci] = newVec;
    }

  }   
      
  //**************************************************************************/
    
  entitySet dmultiMapRepI::domain() const 
  {
    HASH_MAP(int,vector<int> ) :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;
    
    sort( vec.begin(), vec.end() );

    for(size_t i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];
    
    return storeDomain ;
  }

  //**************************************************************************/
  
  entitySet dmultiMapRepI::image(const entitySet &domain) const 
  {
    entitySet     codomain ;
    vector<int>   mapvec;
    HASH_MAP(int,vector<int> )  :: const_iterator   ai;
    
    entitySet :: const_iterator  ei;
    
    for( ei = domain.begin(); ei != domain.end(); ++ei){
      ai = attrib_data.find(*ei);
      if( ai != attrib_data.end() ) {
        mapvec = ai->second;
        for(size_t i = 0; i < mapvec.size(); i++)
          codomain +=   mapvec[i];
      }
    }
    return codomain ;
  }

  //**************************************************************************/
 
  pair<entitySet,entitySet>
  dmultiMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domaini,domainu ;
    HASH_MAP(int,vector<int> )::const_iterator   ai ;
    FORALL(domain(),i) {
      bool vali = true ;
      ai = attrib_data.find(i);
      warn(ai == attrib_data.end()) ;
      bool valu = (!ai->second.size()) ; //begin(i) == end(i)?true:false ;
      for(std::vector<int>::const_iterator vi = ai->second.begin(); vi != ai->second.end(); ++vi) {
        bool in_set = codomain.inSet(*vi) ;
	vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return make_pair(domaini,domainu) ;
  }

  //**************************************************************************/

  multiMap dmultiMapRepI::get_map() 
  {
    multiMap   newmap;
    entitySet::const_iterator  ei;
    store<int> count ;
    entitySet dom = domain() ;
    count.allocate(dom) ;
    for(ei = dom.begin(); ei != dom.end(); ++ei)
      count[*ei] = attrib_data[*ei].size() ;
    newmap.allocate(count) ;
    for(ei = dom.begin(); ei != dom.end(); ++ei) 
      for(int i = 0; i < count[*ei]; i++)
	newmap[*ei][i] = attrib_data[*ei][i];
    
    return newmap;
  }
  
  //**************************************************************************/
    
  ostream &dmultiMapRepI::Print(ostream &s) const 
  {
    HASH_MAP(int,vector<int>)::const_iterator ci;
    vector<int>   newVec;

    s << '{' << domain() << endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
        newVec    = ci->second;
        s << newVec.size() << endl;
      }
    } ENDFORALL ;
    
    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
        newVec    = ci->second;
        for(size_t i = 0; i < newVec.size(); i++)
          s << newVec[i] << "    ";
        s << endl;
      }
    } ENDFORALL ;
    
    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/

  istream &dmultiMapRepI::Input(istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;
    store<int> sizes ;

    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    allocate(sizes) ;
        
    vector<int>   emptyVec;
    int           newEntity;

    FORALL(e,ii) {
      attrib_data[ii] = emptyVec;
      for( int i = 0; i < sizes[ii]; i++){
        s >> newEntity;
        attrib_data[ii].push_back( newEntity ); 
      }
    } ENDFORALL ;
            
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }
    return s ;
  }
  DatatypeP dmultiMapRepI::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  
  frame_info dmultiMapRepI::read_frame_info(hid_t group_id) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  frame_info dmultiMapRepI::write_frame_info(hid_t group_id) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  //**************************************************************************/

  void dmultiMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
    /*
    entitySet::const_iterator ei;
    hsize_t       dimension;
    entitySet     eset;	
    std::vector<int>   vec;


    hid_t vDatatype = H5T_NATIVE_INT;
    Loci::HDF5_ReadDomain( group_id, eset );

    store<int> sizes;
    sizes.allocate(eset);

    dimension  = eset.size();
    hid_t v1Dataset   = H5Dopen(group_id,"ContainerSize");
    hid_t v1Dataspace = H5Dget_space(v1Dataset);

    std::vector<int>  data(dimension);
    H5Dread(v1Dataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    H5Dclose( v1Dataset  );
    H5Sclose( v1Dataspace);

    int indx=0;
    for(ei= eset.begin(); ei != eset.end(); ++ei)
      sizes[*ei] = data[indx++];

    entitySet  ecommon = eset;
    store<int> user_sizes;
    user_sizes.allocate(ecommon);

    for(ei= ecommon.begin(); ei != ecommon.end(); ++ei)
      user_sizes[*ei] =  sizes[*ei]; 

    allocate(user_sizes);

    hid_t v2Dataset   = H5Dopen(group_id,"MultiMap");
    hid_t v2Dataspace = H5Dget_space(v2Dataset);
    H5Sget_simple_extent_dims (v2Dataspace, &dimension, NULL);

    data.resize(dimension);
    H5Dread(v2Dataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( v2Dataset  );
    H5Sclose( v2Dataspace);

    //-----------------------------------------------------------------------
    // For each hyperslab, read the data as contiguous chunck and assign
    // values to the multimap.
    //-----------------------------------------------------------------------

    int vecsize;
    int num_intervals = ecommon.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i<num_intervals;i++) it[i] = ecommon[i];

    indx = 0;
    for(int i=0;i<num_intervals;i++){
      for(int j=it[i].first;j<=it[i].second;j++) {
        vecsize = sizes[j];
        attrib_data[j].clear();
        for( int k = 0; k < vecsize; k++)
          attrib_data[j].push_back( data[indx++] );
      }
    }
    delete [] it;
    */
  }

  //**************************************************************************/

  void dmultiMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const
  {
    warn(true) ;
    /*
    int rank = 1;
    entitySet  :: const_iterator ci;
    HASH_MAP(int,vector<int>)::const_iterator iter;

    entitySet   eset(usr_eset&domain());

    int arraySize = eset.size();
 
    if( arraySize < 1) return;

    Loci::HDF5_WriteDomain(group_id, eset);

    std::vector<int> container, data, vec;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
        iter = attrib_data.find(*ci);
        if( iter == attrib_data.end() ) continue;
        vec  = iter->second;
        container.push_back(vec.size());
        data.insert( data.end(), vec.begin(), vec.end() );
    }

    hsize_t dimension = arraySize;
    hid_t v1Datatype  = H5T_NATIVE_INT;
    hid_t v1Dataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t v1Dataset   = H5Dcreate(group_id, "ContainerSize", v1Datatype,
                                  v1Dataspace, H5P_DEFAULT);
    H5Dwrite(v1Dataset, v1Datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);
    H5Dclose( v1Dataset  );
    H5Sclose( v1Dataspace);

    dimension   = data.size();

    hid_t v2Datatype  = H5T_NATIVE_INT;
    hid_t v2Dataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t v2Dataset   = H5Dcreate(group_id, "MultiMap", v2Datatype,
                                  v2Dataspace, H5P_DEFAULT);
    H5Dwrite(v2Dataset, v2Datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &data[0]);

    H5Dclose( v2Dataset  );
    H5Sclose( v2Dataspace);
    */
  } 

  //**************************************************************************/
  
  dmultiMap::~dmultiMap() {}

  //**************************************************************************/

  void dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  const_dmultiMap::~const_dmultiMap() { }

  //**************************************************************************/

  void const_dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  store_instance::instance_type const_dmultiMap::access() const
  { return READ_ONLY ; }

  //**************************************************************************/


  void inverseMap(dmultiMap &result, const dMap &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) {
    entitySet preloop = input_preimage & input_map.domain() ;
    std::vector<int> tmp_vec ;
    FORALL(input_image,i) {
      result[i] = tmp_vec ;
    } ENDFORALL ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
	result[elem].push_back(i) ;
    } ENDFORALL ;
  }
  void inverseMap(dmultiMap &result, const Map &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) {
    entitySet preloop = input_preimage & input_map.domain() ;
    std::vector<int> tmp_vec ;
    FORALL(input_image,i) {
      result[i] = tmp_vec ;
    } ENDFORALL ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
	result[elem].push_back(i) ;
    } ENDFORALL ;
  }
}


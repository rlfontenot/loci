
#include <DMap.h>
#include <multiMap.h>
#include <Tools/stream.h>
#include <hdf5_readwrite.h>
namespace Loci {
  
  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;
  //**************************************************************************/
  storeRepP dMapRepI::thaw() {
    return getRep() ;
  } 
  storeRepP dMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
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
    
    std::vector<HASH_MAP(int, int) > map_entities(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
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
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	hm[recv_map[j]] = recv_map[j+1];
	j++ ;
      }
    }
    
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) 
      attrib_data[hmi->first] = hmi->second ;
    
    dMap dm ;
    for(HASH_MAP(int, int)::const_iterator hi = attrib_data.begin();  hi != attrib_data.end(); ++hi)
      dm[hi->first] = hi->second ;
    storeRepP sp = dm.Rep() ;
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
  
  //****************************************************************************/
  void dMapRepI::allocate(const entitySet &eset)
  {
    entitySet redundant, newSet;
    entitySet :: const_iterator   ci;

    redundant = domain() -  eset;
    newSet    = eset - domain();
    
    for( ci = redundant.begin(); ci != redundant.end(); ++ci)
         attrib_data.erase(*ci);
    
    for( ci = newSet.begin(); ci != newSet.end(); ++ci)
      attrib_data[*ci] = 0;
    
    dispatch_notify() ;
  }
  
  //**************************************************************************/
  
  dMapRepI::~dMapRepI() 
  { 
    attrib_data.clear(); 
  }
  
  //**************************************************************************/

  storeRep *dMapRepI::new_store(const entitySet &p) const 
  {
    return new dMapRepI(p)  ;
  }
  //**************************************************************************/
  
  storeRepP dMapRepI::remap(const dMap &newmap) const 
  {
    dMap s ;
    entitySet newdomain = newmap.domain() & domain() ;
    pair<entitySet,entitySet> mappimage = preimage(newmap.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = newmap.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(newmap,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(newmap,mapimage) ;
    
    Map m ;
    m.allocate(s.domain()) ;
    FORALL(s.domain(), i) {
      m[i] = s[i] ;
    } ENDFORALL ;
    return m.Rep() ;
    
    //return s.Rep() ;
  }
  
  //**************************************************************************/
  
  void dMapRepI::compose(const dMap &newmap, const entitySet &context) 
  {
    fatal((context-domain()) != EMPTY) ;
    fatal((image(context)-newmap.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      attrib_data[i] = newmap[attrib_data[i]] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  
  void dMapRepI::copy(storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;
    
    fatal((context-domain()) != EMPTY) ;

    fatal((context-s.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/

  void dMapRepI::gather(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;
    
    fatal((m.image(context) - s.domain()) != EMPTY) ; 
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]] ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  void dMapRepI::scatter(const dMap &m,storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    
    FORALL(context,i) {
      attrib_data[m[i]] = s[i] ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  int dMapRepI::pack_size(const entitySet &e) {
    int size ;
    size = sizeof(int) * e.size() ;
    return(size) ;
  }
  
  //**************************************************************************/
  
  void dMapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
    entitySet :: const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      MPI_Pack( &attrib_data[*ci], 1, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
  }
  
  //**************************************************************************/

  void dMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) 
  {
    sequence:: const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci)
      MPI_Unpack( inbuf, insize, &position, &attrib_data[*ci],
		  1, MPI_INT, MPI_COMM_WORLD) ;
  }
  
  //**************************************************************************/
  
  entitySet dMapRepI::domain() const 
  {
    HASH_MAP(int,int ) :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;
    
    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;
    
    sort( vec.begin(), vec.end() );
    
    for(unsigned  int i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];
    
    return storeDomain ;
  }

  //**************************************************************************/
  
  entitySet dMapRepI::image(const entitySet &domain) const 
  {
    
    entitySet codomain ;
    entitySet :: const_iterator  ei;
    HASH_MAP(int,int) ::  const_iterator   ai;   
  
    for( ei = domain.begin(); ei != domain.end(); ++ei){
      ai = attrib_data.find(*ei);
      if( ai != attrib_data.end() ) {
	codomain +=   ai->second;
      }
    }
    return codomain ;
  }
  
  //**************************************************************************/
  
  pair<entitySet,entitySet>
  dMapRepI::preimage(const entitySet &codomain) const  
  {
    entitySet domain ;
    HASH_MAP(int,int) :: const_iterator  ci;
    for(ci = attrib_data.begin(); ci != attrib_data.end(); ++ci)
      if(codomain.inSet( ci->second) )  domain  += ci->first;
    return make_pair(domain,domain);
  }
  
  //**************************************************************************/
  
  multiMap dMapRepI::get_map() 
  {
    multiMap result ;
    store<int> sizes ;
    HASH_MAP(int,int ) :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;
    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;
  
    sort( vec.begin(), vec.end() );
    for(unsigned int i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];
    sizes.allocate(storeDomain) ;
    FORALL(storeDomain,i) {
      sizes[i] = 1 ;
    }ENDFORALL ;
    result.allocate(sizes) ;
    ci =  attrib_data.begin() ;
    FORALL(storeDomain,i) {
      result.begin(i)[0] = ci->second ;
      ++ci ;
    }ENDFORALL ;
    return result ;
  }

  //**************************************************************************/

  ostream &dMapRepI::Print(ostream &s) const 
  {
    HASH_MAP(int,int ) :: const_iterator ci;

    s << '{' << domain() << endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
        s << ci->second << endl ;
      }
    } ENDFORALL ;

    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/

  istream &dMapRepI::Input(istream &s) 
  {

    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      s >> attrib_data[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  //**************************************************************************/

  void dMapRepI::readhdf5(hid_t group_id, entitySet &usr_eset)
  {
    hsize_t       dimension;
    entitySet     eset;	
    vector<int>   vec;

    HDF5_ReadDomain( group_id, eset );
    hid_t vDatatype   = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dopen(group_id,"Map");
    hid_t vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims (vDataspace, &dimension, NULL);

    int *data = new int[dimension];
    H5Dread(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    entitySet  ecommon = eset & usr_eset;

    int num_intervals = ecommon.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i<num_intervals;i++) it[i] = ecommon[i];

    int indx = 0;
    for(int i=0;i<num_intervals;i++){
      for(int j=it[i].first;j<=it[i].second;j++) 
        attrib_data[j] = data[indx++];
    }

    H5Dclose( vDataset   );
    H5Sclose( vDataspace );
    delete [] it;
    delete [] data;
  } 

  //**************************************************************************/

  void dMapRepI::writehdf5(hid_t group_id, entitySet& usr_eset) const
  {
    int       rank = 1;
    hsize_t   dimension;

    entitySet eset(usr_eset & domain());

    int arraySize = eset.size();
    if( arraySize < 1) return;

    HDF5_WriteDomain( group_id, eset);

    vector<int> data(arraySize);
    entitySet :: const_iterator   ei;
    HASH_MAP(int,int) :: const_iterator  ci;

    int indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci == attrib_data.end() )  continue;
      data[indx++] =  ci->second;
    }

    dimension       = arraySize;
    hid_t dataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t datatype  = H5T_NATIVE_INT;
    hid_t dataset   = H5Dcreate(group_id, "Map", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Sclose( dataspace );
    H5Dclose( dataset   );
  } 

  //**************************************************************************/

  dMap::~dMap() {}

  //**************************************************************************/

  void dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }    

  //**************************************************************************/

  const_dMap::~const_dMap() {}

  //**************************************************************************/

  void const_dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  store_instance::instance_type const_dMap::access() const
  { return READ_ONLY ; }

//****************************************************************************/
  void inverseMap(multiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      if(input_image.inSet(input_map[i]))
        sizes[input_map[i]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) {
        sizes[elem] -= 1 ;
        FATAL(sizes[elem] < 0) ;
        result[elem][sizes[elem]] = i ;
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }
}

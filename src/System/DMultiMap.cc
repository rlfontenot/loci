#include <Map.h>
#include <DMultiMap.h>
#include <multiMap.h>

#include <Tools/stream.h>

namespace Loci 
{

  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;
  //-----------------------------------------------------------------

  /*
    multiMap MapRepI::get_map() {

    multiMap result ;
    store<int> sizes ;
    sizes.allocate(store_domain) ;

    FORALL(store_domain,i) {
    sizes[i] = 1 ;
    } ENDFORALL ;


    result.allocate(sizes) ;

    FORALL(store_domain,i) {
    result.begin(i)[0] = base_ptr[i] ;
    } ENDFORALL ;
    return result ;
    }
  */

  //------------------------------------------------------------------

  void dmultiMapRepI::allocate(const entitySet &ptn) {
    vector<int>   emptyVec;
    entitySet :: const_iterator  ci;

    for( ci = ptn.begin(); ci != ptn.end(); ++ci)
      attrib_data[*ci] = emptyVec;

    store_domain = ptn ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------

  void dmultiMapRepI::allocate(const store<int> &sizes) 
  {
    entitySet ptn = sizes.domain() ;
    entitySet :: const_iterator  ci;
    int   veclength;

    for( ci = ptn.begin(); ci != ptn.end(); ++ci) {
      veclength = sizes[*ci];
      attrib_data[*ci].reserve(veclength);
    }

    store_domain = ptn ;
    dispatch_notify() ;
  }

  //***************************************************************************

  dmultiMapRepI::~dmultiMapRepI() 
  {
    attrib_data.clear();
  }

  //***************************************************************************

  storeRep *dmultiMapRepI::new_store(const entitySet &p) const 
  {
    warn(true) ;
    return new dmultiMapRepI()  ;
  }

  //***************************************************************************

  storeRepP dmultiMapRepI::remap(const Map &m) const 
  {
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

    return s.Rep() ;
  }

  //***************************************************************************

  void dmultiMapRepI::compose(const Map &m, const entitySet &context) 
  {
    vector<int>    vec;
    hash_map<int,vector<int> > ::const_iterator ci;

    //-------------------------------------------------------------------------
    // All the entities in the contex should be present in the domain. ie. A->B
    // should be valid.
    //-------------------------------------------------------------------------
    fatal((context-store_domain) != EMPTY) ;

    //-------------------------------------------------------------------------
    // All in the entities in B should be part of Map. Also. i.e. B->C should
    // be valid.
    //-------------------------------------------------------------------------
    fatal((image(context)-m.domain()) != EMPTY) ;

    FORALL(context,i) {
      ci = attrib_data.find(i);
      if( ci != attrib_data.end() ) {
        vec =  ci->second;
        for( int j = 0; j < vec.size(); j++) {
          attrib_data[i][j] =   m[vec[j]];
        }
      }
    } ENDFORALL ;

  }

  //***************************************************************************

  void dmultiMapRepI::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiMap s(st) ;
    vector<int>    newVec;

    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[i].reserve( newVec.size() );
      for( int j = 0; j < newVec.size(); j++)
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;

  }

  //***************************************************************************

  void dmultiMapRepI::gather(const Map &m, storeRepP &st, const entitySet  &context) 
  {

    const_dmultiMap s(st) ;
    vector<int>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[m[i]];
      attrib_data[i].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[i].push_back( newVec[i] );
    } ENDFORALL ;
  }

  //***************************************************************************

  void dmultiMapRepI::scatter(const Map &m, storeRepP &st, const entitySet  &context) 
  {
    const_dmultiMap s(st) ;
    vector<int>    newVec;
    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[m[i]].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[m[i]].push_back( newVec[i] );
    } ENDFORALL ;

  }

  //***************************************************************************
  
  int dmultiMapRepI::pack_size(const  entitySet &e ) 
  {
    int size = 0 ;
    FORALL(e,i) {
      size  +=  attrib_data[i].size();
    } ENDFORALL ;
    
    return( size*sizeof(int) + e.size()*sizeof(int) ) ;
  }

  //***************************************************************************

  void dmultiMapRepI::pack( void *ptr, int &loc, int &size, const entitySet &e) 
  {

    int   *buf;
    int   j, numBytes, numentity = e.size();

    entitySet  :: const_iterator ei;

    //-------------------------------------------------------------------------
    // In the multi-map it is necessary to pack the size of each entity too.
    // For example.
    // 1    2 3 4
    // 2    5 6 2
    // 3    1
    // 4    2 4
    //
    // Then we will pack as size+map
    // 3 2 3 4 3 5 6 2 1 1 2 2 4
    //-------------------------------------------------------------------------

    numBytes = pack_size( e );

    buf   = ( int *) malloc( numBytes );
    if( buf == NULL ) {
      cout << "Warning: Cann't allocate memory for packing data " << endl;
      return;
    }

    vector<int>   newVec;

    hash_map<int,vector<int> > :: const_iterator  ci;
 
    int indx = 0;
    for( ei = e.begin(); ei != e.end(); ++ei) {
      ci = attrib_data.find( *ei );
      if( ci != attrib_data.end() ) {
        newVec = ci->second;
        buf[indx++] = newVec.size();
        for( j = 0; j < newVec.size(); j++) 
          buf[indx++] = newVec[j];
      }
    }

    //------------------------------------------------------------------------
    // At present, we are packing everything at once, hoping that size of
    // buffer is not very large :). May have to change later.
    //------------------------------------------------------------------------

    MPI_Pack(buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;

    free( buf );
  }

  //***************************************************************************

  void dmultiMapRepI::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    int                numBytes;
    int                *buf;
    entitySet          eset;

    //-------------------------------------------------------------------------
    // Unfortunately, pack size is known only for entitySet. So create it from
    // sequence.
    //-------------------------------------------------------------------------
  
    sequence  :: const_iterator  si;
    for( si = seq.begin(); si != seq.end(); ++si) 
      eset +=  *si;

    numBytes = pack_size( eset );

    buf   = ( int *) malloc( numBytes );
    if( buf == NULL ) {
      cout << "Warning: Cann't allocate memory for packing data " << endl;
      return;
    }

    MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

    hash_map<int,vector<int> > :: const_iterator  ci;

    //-------------------------------------------------------------------------
    // At present, I am clearing the attrib_data for the entity, where unpack is
    // being done. I have no idea whether it will lead to some problems later.
    // Chaman Singh Verma :   25 th May 2001
    //-------------------------------------------------------------------------
 
    int indx = 0;
    int vecsize;

    for( si = seq.begin(); si != seq.end(); ++si) {
      attrib_data[*si].clear();
      vecsize = buf[indx++];
      for( int i = 0; i < vecsize; i++)
        attrib_data[*si].push_back( buf[indx++] );
    }

    free( buf );
  }   
      
  //**************************************************************************
    
  entitySet dmultiMapRepI::domain() const 
  {

    hash_map<int,vector<int> > :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;

    sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];

    return storeDomain ;
  }

  //***************************************************************************
    
  entitySet dmultiMapRepI::image(const entitySet &domain) const 
  {
    entitySet     codomain ;
    vector<int>   mapvec;
    hash_map<int,vector<int> >  :: const_iterator   ai;

    entitySet :: const_iterator  ei;
    
    for( ei = domain.begin(); ei != domain.end(); ++ei){
      ai = attrib_data.find(*ei);
      if( ai != attrib_data.end() ) {
        mapvec = ai->second;
        for( int i = 0; i < mapvec.size(); i++)
          codomain +=   mapvec[i];
      }
    }
    return codomain ;
  }

  //***************************************************************************

  pair<entitySet,entitySet> dmultiMapRepI::preimage(const entitySet &codomain) const  
  {
    entitySet domaini,domainu ;
    /*
      FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = begin(i) == end(i)?true:false ;
      for(const int *ip = begin(i);ip!= end(i);++ip) {
      bool in_set = codomain.inSet(*ip) ;
      vali = vali && in_set ;
      valu = valu || in_set ;
      }
      if(vali)
      domaini += i ;
      if(valu)
      domainu += i ;
      } ENDFORALL ;
    */
    return make_pair(domaini,domainu) ;
  }

  //***************************************************************************

  multiMap dmultiMapRepI::get_map() 
  {
    multiMap   newmap;

    cout << "Warning: Dynamic get_map not implemented yet " << endl;
    exit(0);

    return newmap;
  }

  //***************************************************************************
    
  ostream &dmultiMapRepI::Print(ostream &s) const 
  {
    hash_map<int,vector<int> > :: const_iterator ci;
    vector<int>   newVec;

    s << '{' << domain() << endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
        newVec    = ci->second;
        s << newVec.size() << "   ";
      }
    } ENDFORALL ;
    s << endl;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
        newVec    = ci->second;
        for(int i = 0; i < newVec.size(); i++)
          s << newVec[i] << "    ";
        s << endl;
      }
    } ENDFORALL ;

    s << '}' << endl ;

    return s ;
  }

  //***************************************************************************

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
      attrib_data[ii].reserve(sizes[ii]);
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

  //***************************************************************************

  void dmultiMapRepI::readhdf5( H5::Group group, entitySet &user_eset) 
  {

    //-------------------------------------------------------------------------
    // Objective   : Reading Multimap Data in HDF5 Format and storing them as
    //               dynamic multimap.
    // Programmer  : Chaman Singh Verma
    // Date        : 29th May, 2001
    // Place       : ERC, Mississippi State University.
    // Status      :
    // Testing     :
    // Limitations :
    // Future Work :
    // Reference   : See the writehdf module for FORMAT.
    //-------------------------------------------------------------------------

    hash_map<int, vector<int> > :: const_iterator  ci;
    entitySet::const_iterator ei;
    int           *data;
    hsize_t       dimension;
    entitySet     eset;	
    vector<int>   vec;

    if( user_eset.size() < 1) {
        cout << "Warning : Reading entity set is empty " << endl;
        return;
    }

    H5::DataType  datatype = H5::PredType::NATIVE_INT;

    HDF5_ReadDomain( group, eset );
      
    try{
      //-----------------------------------------------------------------------
      // Part II : Read size of each entity i.e. how many mapped elements on 
      //           each entity.
      //-----------------------------------------------------------------------
      store<int> sizes;
      sizes.allocate(eset);

      H5::DataSet   sDataset   = group.openDataSet( "ContainerSize");
      H5::DataSpace sDataspace = sDataset.getSpace();

      sDataspace.getSimpleExtentDims( &dimension, NULL);

      data = new int[dimension];

      sDataset.read( data,  H5::PredType::NATIVE_INT );

      int indx=0, numentities = 0;
      for(ei= eset.begin(); ei != eset.end(); ++ei){
        sizes[*ei]   =  data[indx];
        numentities +=  data[indx];
        indx++;
      }

      //-----------------------------------------------------------------------
      // Now we know all the entities and their size, we can allocate memory 
      // for them.
      //-----------------------------------------------------------------------
      entitySet  ecommon = eset & user_eset;
      store<int> user_sizes;
      user_sizes.allocate(ecommon);

      for(ei= eset.begin(); ei != eset.end(); ++ei)
          if( user_eset.inSet(*ei) ) user_sizes[*ei] =  sizes[*ei]; 

      allocate(user_sizes);

      //-----------------------------------------------------------------------
      // Part III: Read the mapped data( in this case entities )
      //-----------------------------------------------------------------------

      H5::DataSet   vDataset   = group.openDataSet( "multimap");
      H5::DataSpace vDataspace = vDataset.getSpace();

      vDataspace.getSimpleExtentDims( &dimension, NULL);    // file dataspace
      
      //declear the variables used by hyperslab

      hssize_t  start_mem[] = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      // memory dataspace
      int rank = 1;
      dimension = numentities;
      H5::DataSpace mDataspace(rank, &dimension);   // memory dataspace

      //-----------------------------------------------------------------------
      // For each hyperslab, read the data as contiguous chunck and assign
      // values to the multimap.
      //-----------------------------------------------------------------------

      int num_intervals = ecommon.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i<num_intervals;i++) it[i] = ecommon[i];

      //-----------------------------------------------------------------------
      // Calculate the offset of the first element of each interval set
      //-----------------------------------------------------------------------
      store<int> offset;
      offset.allocate( eset );
      indx = 0;
      for(ei= eset.begin(); ei != eset.end(); ++ei) {
          offset[*ei] = indx;
          indx  += sizes[*ei];
      }


      for(int i=0;i<num_intervals;i++){
          
          // How many elements in each hyperslab ...
          count[0]  = 0;          
          for(int j=it[i].first;j<=it[i].second;j++) 
              count[0] += sizes[j];
          int *buf = new int[count[0]];

          // Read the HDF5 Data ...
          foffset[0] = offset[it[i].first];
          mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);	
          vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
          vDataset.read( buf, datatype, mDataspace, vDataspace);

 	       // Create multimap ....
          indx = 0;
          int vecsize;
          for(int j=it[i].first;j<=it[i].second;j++) {
              vecsize = sizes[j];
              attrib_data[j].clear();
              for( int k = 0; k < vecsize; k++)
                   attrib_data[j].push_back( buf[indx++] );
          }

          delete [] buf;
      }

      delete [] it;
    }
    catch( H5::HDF5DatasetInterfaceException error )  { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error ) { error.printerror(); }
  }

  //***************************************************************************

  void dmultiMapRepI::writehdf5( H5::Group group,entitySet& en) const
  {

    //-------------------------------------------------------------------------
    // Objective   : Write Multimap Data in HDF5 Format.
    // Programmer  : Chaman Singh Verma
    // Date        : 29th May, 2001
    // Place       : ERC, Mississippi State University.
    // Status      :
    // Testing     :
    // Limitations :
    // Future Work :
    // Example     :
    // Let us assume that we have following multimap data
    // entity         mapped entities        degree
    // 0              1                         1
    // 1              1 2                       2
    // 2              2 3                       2
    // 5              5 6 7                     3
    // 6              6 7 8                     3
    // 10             1 2 3                     3
    // 11             2                         1
    //
    // Then we have 3 intervals set (0,2),(5,6),(10,11)
    // 
    //-------------------------------------------------------------------------

    hash_map<int,vector<int> > :: const_iterator  ci;
    vector<int>   vec;
    int           numentity;

    hsize_t dimf_domain[1];
    hsize_t dimf_range[1];
    hsize_t dimension[1];

    int rank = 1;    // One dimensional array

    int num_intervals = en.num_intervals();
    dimf_domain[0]    = num_intervals*2;

    interval *it   = new interval[num_intervals];
    for(int i=0;i<num_intervals;i++) it[i] = en[i];
    //
    // Calculate total number of mapped data which will be written in the file.
    //
    int   range, degree, sum_degree = 0;
    for(int i=0;i<num_intervals;i++){
      for(int j=it[i].first;j<=it[i].second;j++) {
        ci = attrib_data.find(j);
        if( ci != attrib_data.end() ) {
          vec         = ci->second;
          sum_degree += vec.size();
        }
      }
    }

    numentity = en.size();

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is written.
    hsize_t   count[1];           // how many positions to select from the dataspace

    dimension[0]   = sum_degree;

    H5::DataSpace fdataspace( rank, dimension );     // Dataspace of file
    H5::DataSpace mdataspace( rank, dimension );     // Dataspace of memory 
    H5::DataType  datatype = H5::PredType::NATIVE_INT;

    H5::DataSet   dataset = group.createDataSet( "multimap", datatype, fdataspace );

    //-------------------------------------------------------------------------
    // For each interval, write data in contiguous fashion. There are 3 data 
    // which have to be written. 
    // Write the mapped entities using hyperslab. Each hyperslab consists of
    // contiguous entitySet.
    // From the example.
    // 1 st hyperslab :   (1,1,2,2,3)
    // 2 nd hyperslab :   (5,6,7,6,7,8)
    // 3 rd hyperslab :   (1,2,3,2)
    //------------------------------------------------------------------------

    for(int i= 0; i< num_intervals; i++){
      count[0] = 0;
      for(int j=it[i].first;j<=it[i].second;j++) {
        ci = attrib_data.find(j);
        if( ci != attrib_data.end() ) {
          vec       = ci->second;
          count[0] += vec.size();
        }
      }

      int *buf = new int[count[0]];

      mdataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      fdataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);

      //
      // Copy data in temporary buffer, so that it can be written in single call.
      //
      int indx = 0;
      for(int j=it[i].first;j<=it[i].second;j++) {
        ci = attrib_data.find(j);
        if( ci != attrib_data.end() ) {
          vec  = ci->second;
          for( int k = 0; k < vec.size(); k++)
            buf[indx++] =  vec[k];
        }
      }

      dataset.write( buf, datatype, mdataspace, fdataspace);
      //
      // change offset for the second record, and clear the buffer.
      //
      foffset[0] += count[0]; 
      delete[] buf;
  }

  //---------------------------------------------------------------------------
  // Second Part: Write the interval set.
  // From the example, we have three intervals namely (0.2),(5,6),(11,12)
  // Write them in the file as 1D array (0,2,5,6,11,12). 
  //---------------------------------------------------------------------------

  int *buf =  new int[2*num_intervals];
  for(int i=0;i<num_intervals;i++){
    it[i]      = en[i];
    buf[i*2]   = it[i].first;
    buf[i*2+1] = it[i].second;
  }

  dimension[0] = 2*num_intervals;

  try {
    H5::DataSpace dDataspace( rank, dimension );
    H5::DataSet   dDataset = group.createDataSet( "Domain", 
                              H5::PredType::NATIVE_INT, dDataspace );
    dDataset.write( buf, H5::PredType::NATIVE_INT );
  }

  catch( H5::HDF5DatasetInterfaceException   error ){error.printerror();}
  catch( H5::HDF5DatatypeInterfaceException  error ){error.printerror();} 
  catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

  delete [] buf;
   
  //---------------------------------------------------------------------------
  // Third part: For each entity write down number of mapped entities
  // For each entity, write number of degree, so that we can reconstruct the
  // data later.
  // From the example. Write in 1D array (1,2,2,3,3,3,1)
  //---------------------------------------------------------------------------

  int indx = 0;
  buf  = new int[numentity];
  entitySet :: const_iterator   ei;

  for( ei = en.begin(); ei != en.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        vec   = ci->second;
        buf[indx++] =  vec.size();
      }
  }

  dimension[0] = en.size();
  try {
    H5::DataSpace sDataspace( rank, dimension );
    H5::DataSet   sDataset = group.createDataSet( "ContainerSize", 
                              H5::PredType::NATIVE_INT, sDataspace );
    sDataset.write( buf, H5::PredType::NATIVE_INT );
  }

  catch( H5::HDF5DatasetInterfaceException   error ){error.printerror();}
  catch( H5::HDF5DatatypeInterfaceException  error ){error.printerror();} 
  catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

  delete [] buf;

  //---------------------------------------------------------------------------
  // Clean up :
  //---------------------------------------------------------------------------
  delete [] it;

  } 

  //***************************************************************************
  
  dmultiMap::~dmultiMap() {}

  //**************************************************************************

  void dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //***************************************************************************

  const_dmultiMap::~const_dmultiMap() { }

  //***************************************************************************

  void const_dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //***************************************************************************

  store_instance::instance_type const_dmultiMap::access() const
  { return READ_ONLY ; }

  //***************************************************************************

    
  void inverseMap(const dmultiMap &input, dmultiMap &result)
  {
    vector<int>    vec;
    entitySet      eset;

    eset  =   input.domain();

    entitySet :: const_iterator   ei;

    for( ei = eset.begin(); ei != eset.end(); ++ei)
      result[*ei].clear();

    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      vec = input[*ei];
      for(int j = 0; j < vec.size(); j++)
        result[vec[j]].push_back(*ei);
    }

  }
  //****************************************************************************
}


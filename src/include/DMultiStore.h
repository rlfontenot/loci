#ifndef DMULTISTORE_H
#define DMULTISTORE_H

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>
#include <hdf5CC/H5cpp.h>

#include <storeVec.h>
#include <Map.h>
#include <multiMap.h>

#include <DMultiMap.h>


namespace Loci {
  template<class T> class dmultiStoreRepI : public storeRep {
    entitySet                 store_domain ;
    hash_map<int,vector<T> >  attrib_data;

    void  hdf5read( H5::Group group, DEFAULT_CONVERTER c,      entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void  hdf5write( H5::Group group, DEFAULT_CONVERTER c,      const entitySet &en) const;
    void  hdf5write( H5::Group group, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( H5::Group group, USER_DEFINED_CONVERTER c, const entitySet &en) const;

  public:

    //  Constructor ...
    dmultiStoreRepI() { }

    dmultiStoreRepI(const entitySet &p) { store_domain=p;}

    dmultiStoreRepI(const store<int> &sizes) { allocate(sizes) ; }

    //  Destructors ...
    virtual ~dmultiStoreRepI() ;

    //  Member function ...
    void allocate(const store<int> &sizes) ;

    virtual void allocate(const entitySet &ptn) ;

    virtual storeRep *new_store(const entitySet &p) const ;

    virtual storeRepP remap(const Map &m) const ;

    virtual void copy(storeRepP &st, const entitySet &context) ;

    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;

    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group, entitySet &user_eset) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;

    hash_map<int,vector<T> > *get_attrib_data(){return &attrib_data; }
  } ;

  //***************************************************************************/
  
  template<class T> class dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T>  storeType ;
    hash_map<int, vector<T> > *attrib_data;
  public:
    typedef vector<T> containerType ;
    dmultiStore() {setRep(new storeType) ;}
    dmultiStore(dmultiStore<T> &var) {setRep(var.Rep()) ;}
    dmultiStore(storeRepP &rp) { setRep(rp) ;}
    
    virtual ~dmultiStore() ;
    virtual void notification() ;
    
    dmultiStore<T> & operator=(dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    
    dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }

    void setSizes(const const_dmultiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }

    entitySet domain() const { return Rep()->domain() ; }

    vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  //*************************************************************************/

  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const dmultiStore<T> &m)
  { return m.Print(s) ; }

  //*************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dmultiStore<T> &m)
  { return m.Input(s) ; }

  //**************************************************************************/
 
  template<class T> 
  dmultiStore<T>::~dmultiStore() {}

  //**************************************************************************/
  
  template<class T> 
  void dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //**************************************************************************/
  
  template<class T> class const_dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T> storeType ;
    hash_map<int, vector<T> >   *attrib_data;
  public:
    //    typedef const_Vect<T> containerType ;
    const_dmultiStore() {setRep(new storeType) ;}
    const_dmultiStore(const_dmultiStore<T> &var) {setRep(var.Rep()) ;}
    const_dmultiStore(storeRepP &rp) { setRep(rp) ;}
    
    virtual ~const_dmultiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
    const_dmultiStore<T> & operator=(const dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dmultiStore<T> & operator=(const const_dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    const entitySet domain() const { return Rep()->domain() ; }

    vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //**************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_dmultiStore<T>::access() const
  { return READ_ONLY ; }

  //**************************************************************************/
  
  template<class T> 
  const_dmultiStore<T>::~const_dmultiStore() {}

  // *************************************************************************
  
  template<class T> 
  void const_dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::allocate(const store<int> &sizes) 
  {

    //------------------------------------------------------------------------
    // Objective : reserve the memory of multiStore using Store..
    //------------------------------------------------------------------------

    entitySet eset = sizes.domain() ;
    entitySet :: const_iterator  ci;
    int   veclength;
    T     newObj;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      veclength = sizes[*ci];
      attrib_data[*ci].reserve(veclength);
      for( int i = 0; i < veclength; i++)
        attrib_data[*ci].push_back( newObj );
    }

    store_domain = eset ;
    dispatch_notify() ;
  }
  //***************************************************************************/
  
  template<class T> 
  void dmultiStoreRepI<T>::allocate(const entitySet &ptn) 
  {
    vector<T>   emptyVec;
    entitySet :: const_iterator  ci;

    for( ci = ptn.begin(); ci != ptn.end(); ++ci)
      attrib_data[*ci] = emptyVec;

    store_domain = ptn ;
    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> 
  dmultiStoreRepI<T>::~dmultiStoreRepI() 
  {
    attrib_data.clear();
  }

  //**************************************************************************/

  template<class T> 
  storeRep *dmultiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dmultiStoreRepI<T>(p) ;
  }

  //**************************************************************************/

  template<class T> 
  storeRepP dmultiStoreRepI<T>::remap(const Map &m) const {
    dmultiStore<T> s ;
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;

    return s.Rep() ;
  }

  //**************************************************************************/
  
  template<class T> 
  void dmultiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[i].reserve( newVec.size() );
      for( int j = 0; j < newVec.size(); j++)
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::gather(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[m[i]];
      attrib_data[i].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[i].push_back( newVec[i] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::scatter(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[m[i]].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[m[i]].push_back( newVec[i] );
    } ENDFORALL ;

    dispatch_notify() ;
  }

  //**************************************************************************/
 
  template <class T> 
  int dmultiStoreRepI<T>::pack_size(const entitySet &e ) 
  {
    int size = 0 ;
    FORALL(e,i) {
      size  +=  attrib_data[i].size();
    } ENDFORALL ;
    
    return( size*sizeof(T) + e.size()*sizeof(int) ) ;
  }

  //**************************************************************************/
  
  template <class T> 
  void dmultiStoreRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &eset ) 
  {
    T    *buf;
    int    numBytes, vecsize;
    vector<T>   newVec;
    entitySet  :: const_iterator ei;

    hash_map<int,vector<T> > :: const_iterator  ci;
 
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find( *ei );
      if( ci != attrib_data.end() ) {
        newVec  = ci->second;
        vecsize = newVec.size();
        MPI_Pack( &vecsize, 1, MPI_INT,  ptr, size, &loc, MPI_COMM_WORLD) ;
        buf  =  new T[vecsize];
        for( int j = 0; j < vecsize; j++) 
          buf[j] = newVec[j];
        numBytes = vecsize*sizeof(T);
        MPI_Pack(buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
        delete [] buf;
      }
    }

  }
  
  //**************************************************************************/
  
  template <class T> 
  void dmultiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    T             *buf;
    int            vecsize, numBytes;

    sequence :: const_iterator  si;

    for( si = seq.begin(); si != seq.end(); ++si) {
      MPI_Unpack(ptr, size, &loc, &vecsize, 1, MPI_INT, MPI_COMM_WORLD) ;
      numBytes =  vecsize*sizeof(T);
      buf      =  new T[vecsize];
      MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;
      attrib_data[*si].clear();
      for( int i = 0; i < vecsize; i++)
        attrib_data[*si].push_back( buf[i] );
      delete [] buf;
    }

  }
  
  //**************************************************************************/
 
  template<class T> 
  store_type dmultiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //***************************************************************************/
  
  template<class T> 
  entitySet dmultiStoreRepI<T>::domain() const 
  {
    hash_map<int,vector<T> > :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;

    sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];

    return storeDomain ;
  }

  //***************************************************************************/
  
  template<class T> 
  std::ostream &dmultiStoreRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << endl ;

    hash_map<int,vector<T> >  :: const_iterator ci;
    vector<T>    vec;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci != attrib_data.end() ) {
        vec =  ci->second;
        s << vec.size() << std::endl ;
      }
    } ENDFORALL ;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci != attrib_data.end() ) {
        vec = ci->second;
        for( int i = 0; i < vec.size(); i++)
          s << vec[i] << "  ";
        s << endl;
      }
    } ENDFORALL ;

    s << '}' << std::endl ;
    return s ;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &dmultiStoreRepI<T>::Input(std::istream &s) 
  {

    entitySet  e ;
    char ch ;

    // Read the opening brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    // Read the interval set ....
    s >> e ;

    // Read the size of each entity map.
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;
    allocate(sizes) ;

    // Read the attribute data
    T          val;
    vector<T>  vec;
    FORALL(e,ii) {
      vec.clear();
      for( int i = 0; i < sizes[ii]; i++) {
        s >> val;
        vec.push_back(val);
      }
      attrib_data[ii] = vec;
    } ENDFORALL ;
            
    // Close the bracket ..
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;

  }

  //***************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::writehdf5(H5::Group group, entitySet &en) const 
  {
    typedef typename hdf5_schema_traits<T> ::Schema_Converter 
      schema_converter;
    schema_converter traits_output_type;

    hdf5write(group, traits_output_type, en);

  }

  //***************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::readhdf5( H5::Group group, entitySet &user_eset) 
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    // Read the entitySet available in file ...
    HDF5_ReadDomain(group, eset);

    // Intersection of entityset in file and user defined. Only common entitySet
    // are read from the file ...
    //
    ecommon = eset & user_eset ;   
    allocate( ecommon );

    // Read the common entitities ...
    hdf5read( group, traits_type, eset,  ecommon );

  }

  //***************************************************************************/
  template <class T> 
  void dmultiStoreRepI<T> :: hdf5read( H5::Group group, DEFAULT_CONVERTER c, 
                                       entitySet &eset, entitySet &user_eset)
  {
     cout << "Fatal error: Default read converter not implemented for dynamic multiStore " << endl;
     exit(0);
  }

  //***************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T> :: hdf5read( H5::Group group, IDENTITY_CONVERTER c, 
                                       entitySet &eset, entitySet &user_eset)
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1, size;

    entitySet::const_iterator ci;
    typedef hdf5_schema_traits<T> traits_type;

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------

    H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   bDataset   = group.openDataSet( "ContainerSize");
    H5::DataSpace bDataspace = bDataset.getSpace();

    bDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

    dimension[0]  = eset.size();
    bDataset.read( ibuf, H5::PredType::NATIVE_INT );

    store<int> container;
    container.allocate( eset );

    indx  = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      container[*ci] = ibuf[indx++];

    delete [] ibuf;

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<unsigned>   offset;
    store<int>        bucket;

    offset.allocate( eset );

    bucket.allocate( user_eset );
    for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
      bucket[*ci] = container[*ci];
    allocate( bucket );

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += container[*ci];
    }

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    T   *data;

    dimension[0] = arraySize;
    H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
    H5::DataSpace vDataspace(rank, dimension);

    H5::DataType vDatatype = traits_type::get_type();
    H5::DataSet  vDataset   = group.openDataSet( "VariableData");

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  container[i];

      data = new T[count[0]];

      foffset[0] = offset[it[k].first];

      mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
      vDataset.read( data, vDatatype, mDataspace, vDataspace);

      indx = 0;
      int size;
      for( int i = it[k].first; i <= it[k].second; i++) {
        attrib_data[i].clear();
        size = container[i];
        for( int m = 0; m < size; m++)
          attrib_data[i].push_back( data[indx++] );
      }

      delete[] data;
    }

  }

  //***************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T> :: hdf5write( H5::Group group, DEFAULT_CONVERTER c, 
                                        const entitySet &eset)  const
  {
     cout << "Fatal error: Default write converter not implemented for dynamic multiStore " << endl;
     exit(0);


  }
  //***************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T> :: hdf5write( H5::Group group, IDENTITY_CONVERTER c, 
                                        const entitySet &eset)  const
  {

    int      rank = 1;
    hsize_t  dimension[1];

    vector<T>   newvec;
    entitySet :: const_iterator  ei;
    hash_map<int,vector<T> >:: const_iterator ci;

    //write out the domain   
    HDF5_WriteDomain(group, eset);

    //-----------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-----------------------------------------------------------------------------
    int  *container = new int[eset.size()];
    size_t  arraySize= 0;
    int     count;

    size_t indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec      = ci->second;
        arraySize  += newvec.size();
        container[indx++] = newvec.size();
      }
    }

    //-----------------------------------------------------------------------------
    // Write the Size of each multiStore ....
    //-----------------------------------------------------------------------------
    dimension[0]=  eset.size();

    try {
      H5::DataSpace sDataspace( rank, dimension );
      H5::DataType  sDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   sDataset  = group.createDataSet( "ContainerSize", sDatatype, sDataspace);

      sDataset.write( container, sDatatype );
    }

    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    delete [] container;

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------
    T  *data, *buf;

    data =  new T[arraySize];

    indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec  = ci->second;
        for( int j = 0; j < newvec.size(); j++) 
          data[indx++] = newvec[j];
      }
    }

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------

    typedef hdf5_schema_traits<T> traits_type;

    dimension[0] =  arraySize;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = traits_type::get_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    //-----------------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------------
    delete [] data;
  }

  //***************************************************************************/

  template<class T>
  void  dmultiStoreRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, 
                                        entitySet &eset, entitySet &user_eset)
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet::const_iterator ci;

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------

    H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   bDataset   = group.openDataSet( "ContainerSize");
    H5::DataSpace bDataspace = bDataset.getSpace();

    bDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

    dimension[0]  = eset.size();
    bDataset.read( ibuf, H5::PredType::NATIVE_INT );

    store<int> container;
    container.allocate( eset );

    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      container[*ci] = ibuf[indx++];

    delete [] ibuf;

    //---------------------------------------------------------------------------
    // Size of each sub-container ....
    //---------------------------------------------------------------------------

    H5::DataType  sDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   sDataset   = group.openDataSet( "SubContainerSize");
    H5::DataSpace sDataspace = sDataset.getSpace();

    sDataspace.getSimpleExtentDims( dimension, NULL);
    ibuf = new int[dimension[0]];

    sDataset.read( ibuf, H5::PredType::NATIVE_INT );

    int maxBucketSize = *max_element( ibuf, ibuf + (int)dimension[0] );

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<unsigned>   offset;
    store<int>        bucket;
    dmultiStore<int>  subcontainer;

    offset.allocate( eset );

    bucket.allocate( user_eset );
    for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
      bucket[*ci] = container[*ci];
    
    allocate( bucket );

    arraySize = 0;
    int indx1 = 0, indx2 = 0, size;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      for( int i = 0; i < container[*ci]; i++)  {
        size = ibuf[indx2];
        arraySize  += size;
        subcontainer[*ci].push_back( size );
        indx2++;
      }
    }

    delete [] ibuf;
    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *data, *buf;

    dimension[0] = arraySize;
    H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
    H5::DataSpace vDataspace(rank, dimension);

    H5::DataType  vDatatype  = converter_traits::get_variable_HDF5_type();
    H5::DataSet   vDataset   = group.openDataSet( "VariableData");

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    buf  = new typename converter_traits::memento_type[maxBucketSize];

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++){
        for( int j = 0; j < subcontainer[i].size(); j++)
          count[0] +=  subcontainer[i][j];
      }

      data = new typename converter_traits::memento_type[count[0]];

      foffset[0] = offset[it[k].first];

      mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
      vDataset.read( data, vDatatype, mDataspace, vDataspace);

      indx = 0;
      int size;
      T      newObj;
      for( int i = it[k].first; i <= it[k].second; i++) {
        attrib_data[i].clear();
        for( int j = 0; j < subcontainer[i].size(); j++) {
          attrib_data[i].push_back( newObj );
          Memento<T> memento( attrib_data[i][j] );
          size = subcontainer[i][j];
          for( int m = 0; m < size; m++)
            buf[m] = data[indx++];
          attrib_data[i][j] = memento.setState( buf, size );
        }
      }

      delete[] data;
    }

    delete[] buf;

  }

  //***************************************************************************/
  
  template <class T> 
  void dmultiStoreRepI<T> :: hdf5write( H5::Group group, USER_DEFINED_CONVERTER c, 
                                        const entitySet &eset)  const
  {   

    int rank = 1;
    hsize_t  dimension[1];

    HDF5_WriteDomain(group, eset);

    //-----------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-----------------------------------------------------------------------------

    entitySet :: const_iterator ei;
    size_t       arraySize= 0;
    int          count, stateSize, *storeSize, maxBucketSize;
    vector<int>  bucketSize;    //  Because we don't know in advance the size
    vector<T>    newvec;
    hash_map<int,vector<T> >:: const_iterator ci;

    storeSize = new int[eset.size()];

    size_t indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec = ci->second;
        storeSize[indx++] = newvec.size();
        for( int j = 0; j < newvec.size(); j++) {
          Memento<T> memento( newvec[j] );
          stateSize  = memento.getSize();
          arraySize += stateSize;
          bucketSize.push_back( stateSize );
        }
      }
    }

    maxBucketSize = *max_element( bucketSize.begin(), bucketSize.end() );

    //-----------------------------------------------------------------------------
    // Write the Size of each multiStore ....
    //-----------------------------------------------------------------------------
    dimension[0]=  eset.size();

    try {
      H5::DataSpace sDataspace( rank, dimension );
      H5::DataType  sDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   sDataset  = group.createDataSet( "ContainerSize", sDatatype, sDataspace);

      sDataset.write( storeSize, sDatatype );
    }

    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    //-----------------------------------------------------------------------------
    // Write the size of each bucket...
    //-----------------------------------------------------------------------------
    dimension[0]=  bucketSize.size();
    int  *bucket = new int[bucketSize.size()];

    for( int i=0; i< bucketSize.size();i++)
      bucket[i] = bucketSize[i];
    bucketSize.clear();

    try {
      H5::DataSpace bDataspace( rank, dimension );
      H5::DataType  bDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   bDataset  = group.createDataSet("SubContainerSize", bDatatype, bDataspace);

      bDataset.write( bucket, bDatatype );
    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *data, *buf;

    data =  new typename converter_traits::memento_type[arraySize];
    buf  =  new typename converter_traits::memento_type[maxBucketSize];
    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------

    indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec =  ci->second;
        for( int j = 0; j < newvec.size(); j++) {
          Memento<T> memento( newvec[j] );
          memento.getState(buf, stateSize);
          for( int i = 0; i < stateSize; i++)
            data[indx++] =  buf[i];
        }
      }
    }

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------

    dimension[0] =  arraySize;

    try {
      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );
    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    //-----------------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------------
    delete [] data;
    delete [] buf;

  }

  //***************************************************************************/
}

#endif

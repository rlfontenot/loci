#ifndef DSTORE_H
#define DSTORE_H

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

#include <Map.h>
#include <Tools/intervalSet.h>
#include <algorithm>
#include <functional>
#include <hdf5_readwrite.h>

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif

namespace Loci {

  using std::hash_map ;
  class Map ;

  template<class T> class dstoreRepI : public storeRep {
    hash_map<int,T>      attrib_data;
    mutable entitySet    store_domain ;

    void  hdf5read( hid_t group, DEFAULT_CONVERTER c,      entitySet &en, entitySet &usr);
    void  hdf5read( hid_t group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( hid_t group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void  hdf5write( hid_t group, DEFAULT_CONVERTER c,      const entitySet &en) const;
    void  hdf5write( hid_t group, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( hid_t group, USER_DEFINED_CONVERTER c, const entitySet &en) const;

  public:
    dstoreRepI(){}
    dstoreRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dstoreRepI()  ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( hid_t group, entitySet &user_eset) ;
    virtual void writehdf5( hid_t group,entitySet& en) const ;
    virtual entitySet domain() const;
    hash_map<int,T> *get_attrib_data() { return &attrib_data; }
    const hash_map<int,T> *get_attrib_data() const { return &attrib_data; }
  } ;

  //**********************************************************************************

  template<class T> 
  void dstoreRepI<T>::allocate(const entitySet &ptn)
  {
    entitySet :: const_iterator ci;

    T   newvalue;
    for( ci = ptn.begin(); ci != ptn.end(); ++ci)
      attrib_data[*ci] =   newvalue;
  
    store_domain = ptn ;
    dispatch_notify() ;
  }

  //**********************************************************************************

  template<class T> 
  std::ostream &dstoreRepI<T>::Print(std::ostream &s) const 
  {
    hash_map<int,T> :: const_iterator   ci;

    s << '{' << domain() << std::endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end() ) s << ci->second << std::endl ;
    } ENDFORALL ;

    s << '}' << std::endl ;

    return s ;
  }

  //**********************************************************************************

  template<class T> 
  std::istream &dstoreRepI<T>::Input(std::istream &s) 
  {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
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
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  //***************************************************************************

  template<class T>  
  dstoreRepI<T>::~dstoreRepI<T>() {
    attrib_data.clear();
  }

  //***************************************************************************
    
  template<class T>  
  entitySet dstoreRepI<T>::domain() const 
  {
    hash_map<int,T> :: const_iterator    ci;
    entitySet          storeDomain;
    std::vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci ) 
      vec.push_back( ci->first ) ;

    std::sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++) 
      storeDomain +=  vec[i];

    return storeDomain ;
  }

  //***************************************************************************

  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dstoreRepI<T>(p) ;
  }

  //***************************************************************************

  template<class T> 
  store_type dstoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //**********************************************************************************

  template<class T> class dstore : public store_instance {
    typedef dstoreRepI<T>  storeType ;
    hash_map<int,T>       *attrib_data;
  public:
    typedef T containerType ;
    dstore() { setRep(new storeType); }
    dstore(dstore &var) { setRep(var.Rep()) ; }
    dstore(storeRepP &rp) { setRep(rp) ; }

    virtual ~dstore() ;
    virtual void notification() ;

    dstore<T> & operator=(dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    dstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

    T &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    const T &elem(int indx) const {
      hash_map<int,T>::const_iterator  citer;

      citer = attrib_data->find(indx);

      if( citer != attrib_data->end() )
        return( (*attrib)[indx] );

      cout << "Error: Entity out of bound " << endl;

    }
  
    T &operator[](int indx) { 
      return elem(indx); 
    }

    const T&operator[](int indx) const { return elem(indx); }

  } ;

  //**********************************************************************************

  template<class T> 
  dstore<T>::~dstore<T>() { }

  //***************************************************************************
    
  template<class T> 
  void dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //****************************************************************************

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const dstore<T> &t)
  { return t.Print(s) ; }

  //****************************************************************************

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dstore<T> &t)
  { return t.Input(s) ; }

  //****************************************************************************

  template<class T> class const_dstore : public store_instance {
    typedef dstoreRepI<T> storeType ;
    hash_map<int,T>      *attrib_data;
  public:
    typedef T containerType ;
    const_dstore() { setRep(new storeType) ; }
    const_dstore(dstore<T> &var) { setRep(var.Rep()) ; }
    const_dstore(const_dstore &var) { setRep(var.Rep()) ; }
    const_dstore(storeRepP &rp) { setRep(rp) ; }

    virtual ~const_dstore() ;
    virtual void notification() ;

    virtual instance_type access() const  ;
        
    const_dstore<T> & operator=(const_dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dstore<T> & operator=(dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    entitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

    const T &elem(int indx) const {
      hash_map<int,T> :: const_iterator  citer;

      citer = attrib_data->find(indx);

      if( citer == attrib_data->end() )
	fatal(true) ;
      
      return ( (*attrib_data)[indx] );
      
    } 
    const T&operator[](int indx) const { return elem(indx); }
      
  } ;

  //****************************************************************************

  template<class T> 
  const_dstore<T>::~const_dstore<T>() { }

  //****************************************************************************
    
  template<class T> 
  void const_dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //***************************************************************************

  template<class T> 
  store_instance::instance_type
  const_dstore<T>::access() const 
  { return READ_ONLY; }
        
  //****************************************************************************

  template<class T> 
  storeRepP dstoreRepI<T>::remap(const Map &m) const 
  {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    dstore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    /*
    store<T> static_store ;
    entitySet tmp_dom = domain();
    static_store.allocate(tmp_dom) ;
    for(entitySet::const_iterator ei = tmp_dom.begin(); ei != tmp_dom.end(); ++ei)
    static_store[*ei] = s[*ei] ;
    */
    return s.Rep() ;
  }

  //***************************************************************************

  template<class T> 
  void dstoreRepI<T>::copy(storeRepP &st, const entitySet &context)  
  {
    const_dstore<T> s(st) ;

    fatal(context != EMPTY ) ;
    fatal((context-domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;

  }

  //***************************************************************************

  template<class T> 
  void dstoreRepI<T>::gather(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dstore<T> s(st) ;

    fatal( context != EMPTY ) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]] ;
    } ENDFORALL ;

  }

  //****************************************************************************

  template<class T> 
  void dstoreRepI<T>::scatter(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dstore<T> s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[m[i]] = s[i] ;
    } ENDFORALL ;

  }

  //***************************************************************************

  template<class T> 
  void dstoreRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset)
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    HDF5_ReadDomain(group_id, eset);

    ecommon = eset & user_eset;
    allocate( ecommon );

    hdf5read( group_id, traits_type, eset,  ecommon );
  }

  //**************************************************************************

  template<class T> 
  void dstoreRepI<T>::writehdf5( hid_t group_id, entitySet& en) const
  {
    if( en.size() < 1) {
      cout << " Entity Set empty for writing " << endl;
      return;
    }

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    hdf5write(group_id, traits_output_type, en );
  }

  //***************************************************************************
    
  template <class T> 
  int dstoreRepI<T>::pack_size( const entitySet &e) 
  {
    int size ;
    size = sizeof(T) * e.size() ;
    return(size) ;
  }

  //***************************************************************************

  template <class T> 
  void dstoreRepI<T>::pack(void *ptr, int &loc, int &size,  const entitySet &e )  
  {

    T   *buf;
    int  numBytes, numentity = e.size();

    numBytes = numentity*sizeof(T);
    buf = ( T *) malloc( numBytes );
    entitySet  :: const_iterator ei;

    int indx = 0;
    for( ei = e.begin(); ei != e.end(); ++ei) {
      buf[indx++] =  attrib_data[*ei];
    }
  
    MPI_Pack(buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;

    free( buf );
 
  }

  //****************************************************************************
  
  template <class T> 
  void dstoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    Loci :: int_type   indx, jndx;
    int                numentity, numBytes;
    T                  *buf;

    for(int i = 0; i < seq.num_intervals(); ++i) {
      numentity =  abs(seq[i].second - seq[i].first) + 1; 
      numBytes  =  numentity*sizeof(T);
      buf       =  (T *) malloc( numBytes );
      MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

      jndx = 0;
      if(seq[i].first > seq[i].second) {
        for(indx = seq[i].first; indx >= seq[i].second; --indx) 
          attrib_data[indx] =  buf[jndx++];
      } else {
        for(indx = seq[i].first; indx <= seq[i].second; ++indx)
          attrib_data[indx] =  buf[jndx++];
      }
      free(buf);
    }
  }  

  //***************************************************************************

  template<class T>
  void  dstoreRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                   entitySet &eset, entitySet &user_eset)
  {
    int      rank=1;
    hsize_t  dimension;

    dimension = user_eset.size();

    char *data = new char[dimension];

    hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_CHAR;
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);

    entitySet :: const_iterator ci;

    istringstream iss;
    iss.str( data );

    for( ci = eset.begin(); ci != eset.end(); ++ci)
      iss >> attrib_data[*ci];

    H5Dclose( vDataset   );
    H5Sclose( mDataspace );
    H5Sclose( vDataspace );

  }

  //***************************************************************************

  template<class T>
  void  dstoreRepI<T> :: hdf5read( hid_t group_id, IDENTITY_CONVERTER c, 
                                   entitySet &eset, entitySet &user_eset)
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1, size;

    entitySet::const_iterator ci;
    typedef data_schema_traits<T> traits_type;

    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    //--------------------------------------------------------------------------
    // Calculate the offset 
    //--------------------------------------------------------------------------

    store<unsigned> offset;
    offset.allocate( eset );

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      offset[*ci] = arraySize++;

    //--------------------------------------------------------------------------

    dimension[0] = arraySize;
    hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);   // memory  dataspace
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDatatype  = traits_type::get_hdf5_type();
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    T   *data;
    for( int k = 0; k < num_intervals; k++) {

      count[0] = it[k].second - it[k].first + 1;
   
      data = new T[count[0]];

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, 
               H5P_DEFAULT, data);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) 
        attrib_data[i] = data[indx++];

      delete[] data;
    }

    H5Dclose( vDataset   );
    H5Tclose( vDatatype  );
    H5Sclose( mDataspace );
    H5Sclose( vDataspace );

  }
  //***************************************************************************

  template<class T>
  void  dstoreRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                   entitySet &eset, entitySet &usr_eset)
  {

    hsize_t  dimension[1];
    size_t   indx = 0, arraySize;
    int      bucketID = 0, rank = 1;
    entitySet::const_iterator  ci;

    typedef data_schema_converter_traits<T> converter_traits; 
    /*
      //-------------------------------------------------------------------------
      // Size of each Bucket ....
      //--------------------------------------------------------------------------

      H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
      H5::DataSet   bDataset   = group.openDataSet( "BucketSize");
      H5::DataSpace bDataspace = bDataset.getSpace();

      bDataspace.getSimpleExtentDims( dimension, NULL);
      int *ibuf = new int[dimension[0]];

      if( converter_traits::var_length== 0) {
      dimension[0] = 1;

      int  fixed_size;
      bDataset.read( &fixed_size, H5::PredType::NATIVE_INT );
      for( int i = 0; i < eset.size(); i++) ibuf[i] = fixed_size;

      } else {
      dimension[0]  = eset.size();
      bDataset.read( ibuf, H5::PredType::NATIVE_INT );
      }

      int maxBucketSize = *std::max_element( ibuf, ibuf + (int)dimension[0] );

      store<int> container;
      container.allocate( eset );

      indx      = 0;
      arraySize = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      container[*ci] = ibuf[indx++];
      arraySize     += container[*ci];
      }
      delete [] ibuf;

      //---------------------------------------------------------------------------
      // Calculate the offset of each entity in file ....
      //---------------------------------------------------------------------------
      store<unsigned>   offset;
      offset.allocate( eset );

      indx     = 0;
      bucketID = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = indx;
      indx       += container[*ci];
      }

      //---------------------------------------------------------------------------
      // Read the data now ....
      //---------------------------------------------------------------------------
      int num_intervals = usr_eset.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

      typedef data_schema_converter_traits<T> converter_traits; 
      typename converter_traits::memento_type *data, *buf;

      dimension[0] = arraySize;
      H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
      H5::DataSpace vDataspace(rank, dimension);

      H5::DataType  vDatatype  = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset   = group.openDataSet( "variable");

      hssize_t  start_mem[] = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      buf  = new typename converter_traits::memento_type[maxBucketSize];

      for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  container[i];

      data = new typename converter_traits::memento_type[count[0]];

      foffset[0] = offset[it[k].first];

      mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
      vDataset.read( data, vDatatype, mDataspace, vDataspace);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
      Memento<T> memento( attrib_data[i] );
      for( int k = 0; k < container[i]; k++)
      buf[k] = data[indx++];
      attrib_data[i] = memento.setState( buf, container[i]);
      }

      delete[] data;
      }

      delete[] buf;
    */

  }

  //***************************************************************************
  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, DEFAULT_CONVERTER c, 
                                   const entitySet &eset)  const
  {
    HDF5_WriteDomain(group_id, eset);

    std::ostringstream oss;
    entitySet :: const_iterator  ei;
    hash_map<int,T>:: const_iterator ci;

    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() )
        oss << ci->second << std::endl ;
    }

    hsize_t dimension[1];
    int rank = 1;
  
    std::string memento = oss.str();
    hsize_t size  =  memento.length();

    dimension[0]  =  size+1;

    hid_t dataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t datatype  = H5T_NATIVE_CHAR;
    hid_t dataset   = H5Dcreate(group_id, "VariableData", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, memento.c_str());

  }
  //***************************************************************************
  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, IDENTITY_CONVERTER c, 
                                   const entitySet &eset)  const
  {

    int      rank = 1;
    hsize_t  dimension[1];

    std::vector<T>   newvec;
    entitySet :: const_iterator  ei;
    hash_map<int,T>:: const_iterator ci;

    //write out the domain   
    HDF5_WriteDomain(group_id, eset);

    int arraySize = eset.size();

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------
    T  *data, *buf;

    data =  new T[arraySize];

    int indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) 
        data[indx++]  = ci->second;
    }

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------

    typedef data_schema_traits<T> traits_type;

    dimension[0] =  arraySize;

    hid_t dataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t datatype  = traits_type::get_hdf5_type();
    hid_t dataset   = H5Dcreate(group_id, "VariableData", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    delete [] data;

  }

  //***************************************************************************

  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                   const entitySet &eset)  const
  {   

    //-----------------------------------------------------------------------------
    // Objective : Write store datatype into HDF5 Format which are user defined
    //             datatypes or STL containers. Such datatypes are first written
    //             in memento class, which store data in NATIVE datatypes. This
    //             memento objects is then written into HDF5 format. A user need
    // to should provide interface to convert data into memento class
    //
    //-----------------------------------------------------------------------------
    /*
      //write out the domain   
      HDF5_WriteDomain(group, eset);

      entitySet :: const_iterator ci;

      //-----------------------------------------------------------------------------
      // Get the sum of each object size and maximum size of object in the 
      // container for allocation purpose
      //-----------------------------------------------------------------------------
      hash_map<int,T> :: const_iterator   iter;

      size_t  arraySize= 0;
      int     stateSize, maxStateSize = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
      Memento<T> memento( iter->second );
      stateSize    = memento.getSize();
      arraySize   += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
      }
      }

      typedef data_schema_converter_traits<T> converter_traits; 
      typename converter_traits::memento_type *data, *buf;

      data =  new typename converter_traits::memento_type[arraySize];
      buf  =  new typename converter_traits::memento_type[maxStateSize];

      //-----------------------------------------------------------------------------
      // Collect state data from each object and put into 1D array
      //-----------------------------------------------------------------------------

      size_t indx = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
      Memento<T> memento( iter->second );
      memento.getState( buf, stateSize);
      for( int i = 0; i <  stateSize; i++) 
      data[indx++] =  buf[i];
      }
      }

      //-----------------------------------------------------------------------------
      // Write (variable) Data into HDF5 format
      //-----------------------------------------------------------------------------
      int rank = 1;
      hsize_t  dimension[1];

      dimension[0] =  arraySize;

      try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "variable", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );

      }
      catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
      catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
      catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

      delete [] data;
      delete [] buf;

      //-----------------------------------------------------------------------------
      // Write Container 
      //-----------------------------------------------------------------------------

      if( converter_traits::var_length== 0) {
      dimension[0]    = 1;
      int bucketSize  = converter_traits::var_size;

      H5::DataType  frameDatatype = H5::PredType::NATIVE_INT;
      H5::DataSpace frameDataspace( rank, dimension );
      H5::DataSet   frameDataset = group.createDataSet( "framing", 
      frameDatatype, frameDataspace);
      frameDataset.write( &bucketSize, frameDatatype );

      } else {
      dimension[0] = eset.size();
      int *vbucket = new int[ eset.size() ];

      indx = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
      Memento<T> memento( iter->second );
      vbucket[indx++] = memento.getSize();
      }
      }

      H5::DataType  frameDatatype = H5::PredType::NATIVE_INT;
      H5::DataSpace frameDataspace( rank, dimension );
      H5::DataSet   frameDataset = group.createDataSet( "BucketSize", 
      frameDatatype, frameDataspace);
      frameDataset.write( vbucket, frameDatatype );
      }
   */

    }
  

}

#endif
  

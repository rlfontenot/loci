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

#ifdef ALLOW_DEFAULT_CONVERTER
    void StringVal( const int &entity, std::string &memento);
    int  get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    void packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                    const sequence &seq) ;

    void  hdf5read( hid_t group, DEFAULT_CONVERTER c,      entitySet &en, entitySet &usr);
    void  hdf5write( hid_t group, DEFAULT_CONVERTER c,      const entitySet &en) const;
#endif
    
    void  hdf5read( hid_t group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( hid_t group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void  hdf5write( hid_t group, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( hid_t group, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int   get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void  packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void  unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq) ;
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void  packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;

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

  //*************************************************************************/

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

  //*********************************************************************/

  template<class T> 
  std::ostream &dstoreRepI<T>::Print(std::ostream &s) const 
  {
    hash_map<int,T> :: const_iterator   ci;

    s << '{' << domain() << std::endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci == attrib_data.end() ) continue;
      s << ci->second << std::endl ;
    } ENDFORALL ;

    s << '}' << std::endl ;

    return s ;
  }

  //************************************************************************/

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

  //*************************************************************************/

  template<class T>  
  dstoreRepI<T>::~dstoreRepI<T>() {
    attrib_data.clear();
  }

  //*************************************************************************/
    
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

  //*************************************************************************/

  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dstoreRepI<T>(p) ;
  }

  //*************************************************************************/

  template<class T> 
  store_type dstoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //************************************************************************/

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

  //*************************************************************************/

  template<class T> 
  dstore<T>::~dstore<T>() { }

  //*************************************************************************/
    
  template<class T> 
  void dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const dstore<T> &t)
  { return t.Print(s) ; }

  //*************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dstore<T> &t)
  { return t.Input(s) ; }

  //************************************************************************/

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

  //************************************************************************/

  template<class T> 
  const_dstore<T>::~const_dstore<T>() { }

  //************************************************************************/
    
  template<class T> 
  void const_dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_dstore<T>::access() const 
  { return READ_ONLY; }
        
  //*************************************************************************/

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

  //*************************************************************************/

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

  //*************************************************************************/

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

  //**************************************************************************/

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

  //*************************************************************************/
  template <class T> 
  int dstoreRepI<T>::pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }

  //*******************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  int dstoreRepI<T>::get_mpi_size( DEFAULT_CONVERTER c,
                                  const entitySet &eset)
  {
    IDENTITY_CONVERTER cc;
    return get_mpi_size( cc, eset);

    std::ostringstream oss;
    int  size;

    FORALL(eset,ii){
      oss << attrib_data[ii] << std::endl ;
    }ENDFORALL ;

    std::string memento = oss.str();

    size = memento.length() + eset.size()*sizeof(int);

    return( size );
  }
#endif

  //*******************************************************************/

  template <class T>
  int dstoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                  const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }

  //*******************************************************************/
  template <class T>
  int dstoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                  const entitySet &eset)
  {

    int       size, numBytes;
    entitySet  ecommon;
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    T   obj;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      obj  = attrib_data[*ci];
      typename converter_traits::Converter_Type cvtr(obj);
      size      = cvtr.getSize();
      numBytes += size*sizeof(typename converter_traits::Converter_Base_Type) ;
    }

    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
  }
  //*******************************************************************/
  template <class T> 
  void dstoreRepI<T>::pack( void *outbuf, int &position, int &size, 
                           const entitySet &usr_eset )  
  {
    entitySet ecommon;

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    ecommon = domain()&usr_eset;

    packdata( traits_type, outbuf, position, size, ecommon);

  }
 
  //*******************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  void dstoreRepI<T>::StringVal( const int &entity, std::string &memento)
  {
    std::ostringstream oss;

    oss << base_ptr[entity] << endl;

    memento = oss.str();
  }

  template <class T> 
  void dstoreRepI<T>::packdata(DEFAULT_CONVERTER c, void *outbuf,
                              int &position, int outcount,
                              const entitySet &eset )  
  {
    std::ostringstream oss, os1;
    std::string        memento;
    int   bufSize, currpos = 0;

    entitySet :: const_iterator ci;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      StringVal(*ci, memento );
      bufSize = memento.length();
      MPI_Pack( &bufSize, 1, MPI_INT, outbuf, outcount, 
                &position, MPI_COMM_WORLD) ;
      MPI_Pack( &memento[0], bufSize, MPI_CHAR, outbuf, outcount, 
                &position, MPI_COMM_WORLD) ;
    }
  }
#endif
  //*******************************************************************/
  template <class T> 
  void dstoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                              int &position, int outcount,
                              const entitySet &eset )  
  {

    entitySet :: const_iterator ci;

    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      MPI_Pack( &attrib_data[*ci], sizeof(T), MPI_BYTE, outbuf,outcount, 
                &position, MPI_COMM_WORLD) ;
  }

  //*******************************************************************/

  template <class T>
  void dstoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                               int &position, int outcount, 
                               const entitySet &eset )
  {

    entitySet :: const_iterator ci;
    entitySet  ecommon;

    ecommon = domain()&eset;
    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------
    size_t  incount = 0;
    int     stateSize, maxStateSize = 0; 
    typedef data_schema_traits<T> schema_traits ;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( attrib_data[*ci] );
      stateSize    = cvtr.getSize();
      incount     += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
    }

    typename schema_traits::Converter_Base_Type *inbuf;

    int typesize = sizeof(typename schema_traits::Converter_Base_Type);

    inbuf = new typename schema_traits::Converter_Base_Type[maxStateSize];

    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------

    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( attrib_data[*ci]);
      cvtr.getState( inbuf, stateSize);

      MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
               MPI_COMM_WORLD);

      incount =  stateSize*typesize;
      MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position, 
               MPI_COMM_WORLD) ;
    }
    delete [] inbuf;
  }

  //*******************************************************************/

  template <class T> 
  void dstoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size, seq); 
  }

  //*********************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void dstoreRepI<T>::unpackdata(DEFAULT_CONVERTER c, void *inbuf,
                                int &position,  int insize,
                                const sequence &seq) 
  {
    IDENTITY_CONVERTER cc;
    unpackdata(cc, inbuf, position,  insize, seq);
    return;

    char *outbuf;
    int   outcount;
    sequence:: const_iterator ci;
    entitySet eset(seq);

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      MPI_Unpack( inbuf, insize, &position, &outcount, 1, MPI_INT, 
                  MPI_COMM_WORLD) ;
      outbuf   = new char[outcount];

      MPI_Unpack( inbuf, insize, &position, outbuf, outcount, MPI_CHAR, 
                  MPI_COMM_WORLD) ;

      std::istringstream iss(outbuf);
      iss >> base_ptr[*ci];
      delete [] outbuf;
    }
  }
#endif
  //*********************************************************************/
  template <class T> 
  void dstoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                int &position,  int insize,
                                const sequence &seq) 
  {
    sequence:: const_iterator ci;

    for( ci = seq.begin(); ci != seq.end(); ++ci) 
          MPI_Unpack( inbuf, insize, &position, &attrib_data[*ci],
                      sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
  }
  //*********************************************************************/
  template <class T> 
  void dstoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                 int &position, int insize, const sequence &seq) 
  {

    sequence :: const_iterator ci;

    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-------------------------------------------------------------------
    int  stateSize, outcount;

    typedef  data_schema_traits<T> converter_traits; 
    typename converter_traits::Converter_Base_Type *outbuf;

    int typesize = sizeof(converter_traits::Converter_Base_Type);

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( !store_domain.inSet( *ci ) ) {
        std::cout << "Warning: entity not in entitySet  : " << *ci << endl;
        continue;
      }
      outcount = sizeof(int);
      MPI_Unpack(inbuf, insize, &position, &stateSize, 1, 
                 MPI_INT, MPI_COMM_WORLD) ;

      outbuf = new typename converter_traits::Converter_Base_Type[stateSize];

      outcount = stateSize*typesize;
      MPI_Unpack(inbuf, insize, &position, outbuf, outcount, 
                 MPI_BYTE, MPI_COMM_WORLD) ;

      typename converter_traits::Converter_Type cvtr( attrib_data[*ci] );
      cvtr.setState( outbuf, stateSize);
      delete [] outbuf;
    }

  }
  //*********************************************************************/
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


#ifdef ALLOW_DEFAULT_CONVERTER
  template<class T>
  void  dstoreRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                   entitySet &eset, entitySet &user_eset)
  {
    int      rank=1;
    hsize_t  dimension;

    dimension = user_eset.size();

    std::vector<char> data(dimension);

    hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_CHAR;
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, &data[0]);

    entitySet :: const_iterator ci;

    istringstream iss;
    iss.str( &data[0] );

    for( ci = eset.begin(); ci != eset.end(); ++ci)
      iss >> attrib_data[*ci];

    H5Dclose( vDataset   );
    H5Sclose( mDataspace );
    H5Sclose( vDataspace );

  }
#endif
  //*************************************************************************/

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

    //------------------------------------------------------------------------
    // Calculate the offset 
    //------------------------------------------------------------------------

    store<unsigned> offset;
    offset.allocate( eset );

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      offset[*ci] = arraySize++;

    //------------------------------------------------------------------------
    DatatypeP dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    dimension[0] = arraySize;
    hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);   // memory  dataspace
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    std::vector<T>  data;
    for( int k = 0; k < num_intervals; k++) {
      count[0] = it[k].second - it[k].first + 1;

      if( count[0] > data.size() ) data.resize(count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, 
               H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) 
        attrib_data[i] = data[indx++];
    }

    H5Dclose( vDataset   );
    H5Tclose( vDatatype  );
    H5Sclose( mDataspace );
    H5Sclose( vDataspace );

  }
  //*************************************************************************/

  template<class T>
  void  dstoreRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                   entitySet &eset, entitySet &usr_eset)
  {

    hsize_t  dimension;
    size_t   indx = 0, arraySize;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;
    int      rank = 1;
    entitySet::const_iterator  ci;

    //---------------------------------------------------------------
    // Size of each Bucket ....
    //---------------------------------------------------------------
    vDatatype  = H5T_NATIVE_INT;
    vDataset   = H5Dopen(group_id,"ContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    std::vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);
    H5Dclose(vDataset  );
    H5Sclose(vDataspace);

    store<int> container;
    container.allocate( eset );

    indx      = 0;
    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      container[*ci] = ibuf[indx++];
      arraySize     += container[*ci];
    }

    //-------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //-------------------------------------------------------------------
    store<unsigned>   offset;
    offset.allocate( eset );

    indx     = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = indx;
      indx       += container[*ci];
    }

    //------------------------------------------------------------------
    // Read the data now ....
    //------------------------------------------------------------------
    int num_intervals = usr_eset.num_intervals();
    if( num_intervals == 0) {
      std::cout << "Warning: Number of intervals are zero : " << endl;
      return;
    }

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

    typedef typename data_schema_traits<T>::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    DatatypeP atom_type = traits_type::get_type() ;
    vDatatype = atom_type->get_hdf5_type();

    dimension  = arraySize;
    vDataset   = H5Dopen(group_id,"VariableData");
    vDataspace = H5Dget_space(vDataset);
    mDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    std::vector<dtype> data;

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  container[i];
      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];
      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride,
                          count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride,
                          count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
          typename data_schema_traits<T>::Converter_Type cvtr( attrib_data[i]);
          cvtr.setState( &data[0]+indx, container[i] );
          indx += container[i];
      }
    }

    H5Tclose(vDatatype );
    H5Dclose(vDataset  );
    H5Sclose(vDataspace);
    H5Sclose(mDataspace);
  }
  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::writehdf5( hid_t group_id, entitySet& usr_eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet eset(usr_eset &domain());

    HDF5_WriteDomain(group_id, eset);

    hdf5write(group_id, traits_output_type, eset);
  }
  //*************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, DEFAULT_CONVERTER c, 
                                   const entitySet &eset)  const
  {
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

    dimension        = memento.length()+1;
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_CHAR;
    hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, memento.c_str());

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);

  }
#endif
  //*************************************************************************/
  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, IDENTITY_CONVERTER c, 
                                   const entitySet &eset)  const
  {
    int      rank = 1;
    hsize_t  dimension;

    std::vector<T>   newvec;
    entitySet :: const_iterator  ei;
    hash_map<int,T>:: const_iterator ci;

    int arraySize = eset.size();

    if( arraySize < 1) return;

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------
    std::vector<T>  data(arraySize);

    int indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci == attrib_data.end() )  continue;
      data[indx++]  = ci->second;
    }

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------

    typedef data_schema_traits<T> traits_type;
    DatatypeP dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    dimension        = arraySize;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

  }

  //*************************************************************************/

  template <class T> 
  void dstoreRepI<T> :: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                   const entitySet &eset)  const
  {   
    T   Obj;
    entitySet :: const_iterator ci;
    hid_t   vDataspace, vDataset, vDatatype;

    hash_map<int,T> :: const_iterator   iter;

    //-----------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-----------------------------------------------------------------------------

    typedef data_schema_traits<T> schema_traits ;

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        Obj = iter->second;
        typename schema_traits::Converter_Type cvtr( Obj );
        stateSize    = cvtr.getSize();
        arraySize   += stateSize;
        maxStateSize = max( maxStateSize, stateSize );
      }
    } 

    typedef typename schema_traits::Converter_Base_Type dtype;

    dtype *data, *buf;

    data =  new dtype[arraySize];
    buf  =  new dtype[maxStateSize];

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        Obj = iter->second;
        typename schema_traits::Converter_Type cvtr( Obj );
        cvtr.getState( data+indx, stateSize);
        indx +=stateSize ;
      }
    }

    //--------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //--------------------------------------------------------------------

    typedef data_schema_traits<dtype> traits_type;
    DatatypeP atom_type = traits_type::get_type() ;
    vDatatype = atom_type->get_hdf5_type();

    int rank = 1;
    hsize_t  dimension;

    dimension  =  arraySize;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
                           vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    delete [] data;
    delete [] buf;

    //--------------------------------------------------------------------
    // Write Container 
    //--------------------------------------------------------------------
    dimension    = eset.size();
    int *vbucket = new int[ eset.size() ];

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        Obj = iter->second;
        typename schema_traits::Converter_Type cvtr( Obj );
        vbucket[indx++] =  cvtr.getSize();
      }
    }

    vDatatype  = H5T_NATIVE_INT;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "ContainerSize",
                           vDatatype, vDataspace, H5P_DEFAULT);

    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vbucket);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    delete [] vbucket;

  }
  

}

#endif
  

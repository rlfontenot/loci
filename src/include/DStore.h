#ifndef DSTORE_H
#define DSTORE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
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

#include <Tools/block_hash.h>


namespace Loci {

  class Map ;

  template<class T> class dstoreRepI : public storeRep {
    block_hash<T>  attrib_data;

    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en);

    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

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
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
  public:
    dstoreRepI(){}
    dstoreRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~dstoreRepI()  ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    virtual entitySet domain() const;
    block_hash<T> *get_attrib_data() { return &attrib_data; }
    const block_hash<T> *get_attrib_data() const { return &attrib_data; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ; 

  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::allocate(const entitySet &eset)
  {
    entitySet :: const_iterator ci;

    entitySet dom = domain() ;
    entitySet remove = dom - eset ;
    entitySet add = eset - dom ;

    attrib_data.erase_set(remove) ;

    for( ci = add.begin(); ci != add.end(); ++ci)
      attrib_data.access(*ci) ;
  
    dispatch_notify() ;
  }

  template<class T>
    void dstoreRepI<T>::shift(int_type offset) {
    entitySet new_domain = domain() ;
    new_domain >>= offset ;
    allocate(new_domain) ;
  }
  //*********************************************************************/

  template<class T> 
  std::ostream &dstoreRepI<T>::Print(std::ostream &s) const 
  {
    entitySet :: const_iterator it;
    entitySet dom = domain() ;

    s << '{' << dom << std::endl ;

    for(it = dom.begin();it!=dom.end();++it) {
      Loci::streamoutput(&(attrib_data.elem(*it)),1,s) ;
    }

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
        
    FORALL(e,ii) {
      Loci::streaminput(&attrib_data.access(ii),1,s) ;
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
  }

  //*************************************************************************/
    
  template<class T>  
  entitySet dstoreRepI<T>::domain() const 
  {
    return attrib_data.domain() ;
  }

  //*************************************************************************/

  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dstoreRepI<T>(p) ;
  }
  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p, const int* cnt) const 
    {
      storeRep* sp = 0 ;
      cerr << " This method should not be called for a dstore " << endl ;
      return sp ;
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
    block_hash<T>       *attrib_data;
  public:
    typedef T containerType ;
    dstore() { setRep(new storeType); }
    dstore(const dstore &var) { setRep(var.Rep()) ; }
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
      return attrib_data->access(indx) ;
    }

    const T &elem(int indx) const { return attrib_data->elem(indx) ; }
  
    T &operator[](int indx) { return elem(indx); }
    const T&operator[](int indx) const { return elem(indx); }
    const T&operator()(int indx) const { return elem(indx) ; }
    

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
    block_hash<T>      *attrib_data;
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
      return attrib_data->elem(indx) ;
    } 
    const T&operator[](int indx) const { return elem(indx); }

    const T&operator()(int indx) const { return elem(indx) ; }
      
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
  storeRepP dstoreRepI<T>::remap(const dMap &m) const 
  {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    dstore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    
    store<T> static_store ;
    entitySet tmp_dom = s.domain();
    static_store.allocate(tmp_dom) ;
    for(entitySet::const_iterator ei = tmp_dom.begin(); ei != tmp_dom.end(); ++ei)
      static_store[*ei] = s[*ei] ;
    return static_store.Rep() ;
    
    //return s.Rep() ;
  }

  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::copy(storeRepP &st, const entitySet &context)  
  {
    const_dstore<T> s(st) ;
    fatal((context-domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;

  }

  //*************************************************************************/

  template<class T> 
    void dstoreRepI<T>::gather(const dMap &m, storeRepP &st, const entitySet &context) 
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
  void dstoreRepI<T>::scatter(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dstore<T> s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);

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

    int       size, numBytes=0;
    entitySet  ecommon;
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename converter_traits::Converter_Type cvtr(attrib_data[*ci]);
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
    int  stateSize, incount;
    entitySet :: const_iterator ci;

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  inbuf;

    int typesize = sizeof(dtype);
    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( attrib_data[*ci]);

      stateSize    = cvtr.getSize();
      if( stateSize > static_cast<int>(inbuf.size()) ) inbuf.resize(stateSize);

      cvtr.getState( &inbuf[0], stateSize);
      MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
               MPI_COMM_WORLD);

      incount =  stateSize*typesize;
      MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
               MPI_COMM_WORLD) ;
    }
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
    int  stateSize, outcount;
    sequence :: const_iterator ci;

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  outbuf;

    int typesize = sizeof(dtype);
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      MPI_Unpack(inbuf, insize, &position, &stateSize, 1, 
                 MPI_INT, MPI_COMM_WORLD) ;

      if( stateSize > static_cast<int>(outbuf.size()) ) outbuf.resize(stateSize);
 
      outcount = stateSize*typesize;
      MPI_Unpack(inbuf, insize, &position, &outbuf[0], outcount, 
                 MPI_BYTE, MPI_COMM_WORLD) ;

      typename schema_traits::Converter_Type cvtr( attrib_data[*ci] );
      cvtr.setState( &outbuf[0], stateSize);
    }

  }
  template<class T> 
    frame_info dstoreRepI<T>::read_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return read_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info dstoreRepI<T>::read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
    frame_info dstoreRepI<T>::read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
    frame_info dstoreRepI<T>::write_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return write_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info dstoreRepI<T>::write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
    frame_info dstoreRepI<T>::write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
    DatatypeP dstoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
    DatatypeP dstoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
    DatatypeP dstoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  
  //*********************************************************************/
  template<class T> 
    void dstoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
    /*
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    HDF5_ReadDomain(group_id, eset);

    ecommon = eset;
    allocate( ecommon );

    hdf5read( group_id, traits_type, eset,  ecommon );
    */
  }

  //*************************************************************************/

  template<class T>
  void  dstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    warn(true) ;
    /*
    hsize_t dimension;
    size_t indx = 0, arraySize;
    int    rank = 1;

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

    dimension = arraySize;
    hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);   // memory  dataspace
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
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
    */
  }
  //*************************************************************************/

  template<class T>
    void  dstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
    {
    warn(true) ;
    /*
    hsize_t  dimension;
    size_t   indx = 0, arraySize;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;
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
    */
  }
  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const
  {
    warn(true) ;
    /*
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet eset(usr_eset &domain());

    HDF5_WriteDomain(group_id, eset);

    hdf5write(group_id, traits_output_type, eset);
    */
  }

  //*************************************************************************/
  template <class T> 
  void dstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset)  const
  {
    warn(true) ;
    /*
    int      rank = 1;
    hsize_t  dimension;

    std::vector<T>   newvec;
    entitySet :: const_iterator  ei;

    int arraySize = eset.size();

    if( arraySize < 1) return;

    //------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //------------------------------------------------------------------------
    std::vector<T>  data(arraySize);

    int indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      data[indx++]  = attrib_data.elem(*ei) ;
    }

    //-------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-------------------------------------------------------------------------

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
    */
  }

  //*************************************************************************/

  template <class T> 
  void dstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &eset)  const
    {   
      warn(true) ;
      /*
    T         Obj;
    hid_t     vDataspace, vDataset, vDatatype;
    int       rank = 1;
    hsize_t   dimension;

    entitySet :: const_iterator ci;

    //-------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-------------------------------------------------------------------------

    typedef data_schema_traits<T> schema_traits ;

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      stateSize    = cvtr.getSize();
      arraySize   += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
    } 

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------
    typedef typename schema_traits::Converter_Base_Type dtype;

    dtype *data = new dtype[arraySize];

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      cvtr.getState( data+indx, stateSize);
      indx +=stateSize ;
    }

    //--------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //--------------------------------------------------------------------
    typedef data_schema_traits<dtype> traits_type;
    DatatypeP atom_type = traits_type::get_type() ;
    vDatatype = atom_type->get_hdf5_type();

    dimension  = arraySize;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
                           vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    delete [] data;

    //--------------------------------------------------------------------
    // Write Container 
    //--------------------------------------------------------------------
    dimension    = eset.size();
    std::vector<int> vbucket(eset.size());

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      vbucket[indx++] =  cvtr.getSize();
    }

    vDatatype  = H5T_NATIVE_INT;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "ContainerSize",
                           vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vbucket[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
      */
  }
  //*************************************************************************/
  
}

#endif
  

#ifndef DSTOREVEC_H
#define DSTOREVEC_H 

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <hdf5_readwrite.h>

#include <store_rep.h>

#include <Tools/lmutex.h>

#include <Map.h>
#include <multiMap.h>

#include <algorithm>

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif

namespace Loci {
  using std::hash_map ;
  using std::vector;

  template<class T> class dstoreVecRepI : public storeRep {

    lmutex                    mutex ;
    entitySet                 store_domain ;
    int                       size;
    hash_map<int,std::vector<T> >  attrib_data;
#ifdef ALLOW_DEFAULT_CONVERTER
    void hdf5write( hid_t group_id, DEFAULT_CONVERTER  c,
                    const entitySet &en ) const;
    void hdf5read( hid_t group, DEFAULT_CONVERTER      c,
                   entitySet &en, entitySet &usr );
    int get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    void packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int &size,
                    const sequence &seq) ;
    void StringVal( const int &entity, const int &ivec, std::string &memento);
#endif
    
    void hdf5read( hid_t group_id, IDENTITY_CONVERTER     c, entitySet &en, entitySet &usr );
    void hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr );

    void hdf5write( hid_t group, IDENTITY_CONVERTER c, const entitySet &en ) const;
    void hdf5write( hid_t group, USER_DEFINED_CONVERTER c, const entitySet &en ) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size, const entitySet &e ) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size, const entitySet &e ) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size, const sequence &seq) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size, const sequence &seq) ;

  public:
    dstoreVecRepI() {}

    dstoreVecRepI(const entitySet &p) { allocate(p) ; }

    virtual ~dstoreVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e ) ;
    virtual void pack(void * ptr, int &loc, int &size, const entitySet &e ) ;
    virtual void unpack(void * ptr, int &loc, int &size, const sequence &seq) ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(  hid_t group, entitySet &eset) ;
    virtual void writehdf5( hid_t group, entitySet &eset) const ;

    virtual void set_elem_size(int sz) ;

    hash_map<int,std::vector<T> > *get_attrib_data() { return &attrib_data; }
    int get_size() const { return size; }
  } ;

  //*************************************************************************/

  template<class T> 
  std::ostream &dstoreVecRepI<T>::Print(std::ostream &s) const
  {

    s << '{' << domain() << std::endl ;

    s << size << std::endl ;

    hash_map<int, std::vector<T> >  :: const_iterator   ci;
    std::vector<T>   newVec;
    
    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci != attrib_data.end() ) {
        newVec = ci->second;
        for( int i = 0; i < newVec.size(); i++)
          s << newVec[i] <<  "   ";
        s << std::endl;
      }
    }ENDFORALL ;

    s << '}' << std::endl ;

    return s ;
  }

  //************************************************************************/

  template<class T> 
  std::istream &dstoreVecRepI<T>::Input(std::istream &s)
  {
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    entitySet e ;

    s >> e ;
    s >> size ;

    std::vector<T>  newVec(size);
    allocate(e) ;

    FORALL(e,ii) {
      for( int i = 0; i < size; i++)
        s  >> newVec[i];
      attrib_data[ii] = newVec;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }


  //************************************************************************/

  template<class T> 
  void dstoreVecRepI<T>::allocate(const entitySet &ptn) 
  {

    entitySet :: const_iterator  ci;

    std::vector<T>   newVec;

    if( size > 1) newVec.reserve( size );

    for( ci = ptn.begin(); ci != ptn.end(); ++ci)
      attrib_data[*ci] = newVec;

    store_domain = ptn ;
    dispatch_notify() ;

  }

  //************************************************************************/

  template<class T> 
  dstoreVecRepI<T>::~dstoreVecRepI<T>() 
  {
    attrib_data.clear();
  }

  //***********************************************************************/

  template<class T>
  storeRep *dstoreVecRepI<T>::new_store(const entitySet &p) const 
  {
    return new dstoreVecRepI<T>(p)  ;
  }

  //*************************************************************************/

  template<class T> 
  store_type dstoreVecRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //*************************************************************************/

  template<class T> 
  entitySet dstoreVecRepI<T>::domain() const 
  {
    hash_map<int,std::vector<T> > :: const_iterator    ci;
    entitySet          storeDomain;
    std::vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;

    std::sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++) 
      storeDomain +=  vec[i];

    dispatch_notify() ;

    return storeDomain ;
  }

  //*************************************************************************/

  template<class T> 
  void dstoreVecRepI<T>::set_elem_size(int sz) 
  {

    mutex.lock() ;

    fatal(sz<1);

    size = sz ;
    allocate(store_domain) ;

    mutex.unlock() ;

  }

  //*************************************************************************/
      
  template<class T> class dstoreVec : public store_instance {
    typedef dstoreVecRepI<T>     storeType ;
    hash_map<int, std::vector<T> > *attrib_data;
    int                          size;
  public:
    typedef std::vector<T> containerType ;
    dstoreVec() {setRep(new storeType) ;}
    dstoreVec(dstoreVec<T> &var) {setRep(var.Rep()) ;}
    dstoreVec(storeRepP &rp) { setRep(rp) ;}

    virtual ~dstoreVec() ;

    virtual void notification() ;

    dstoreVec<T> & operator=(dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    dstoreVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      Rep()->set_elem_size(size) ;
    }

    int getVecSize() {
      return( Rep()->get_size() );
    }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    std::vector<T> &operator[](int indx) {
      return( elem(indx) );
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //*************************************************************************/

  template<class T> 
  dstoreVec<T>::~dstoreVec<T>() { }

  //*************************************************************************/

  template<class T> 
  void dstoreVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) 
      attrib_data = p->get_attrib_data() ;

    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const dstoreVec<T> &t)
  { return t.Print(s) ; }

  //*************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dstoreVec<T> &t)
  { return t.Input(s) ; }

  //*************************************************************************/

  template<class T> class const_dstoreVec : public store_instance {
    typedef dstoreVecRepI<T> storeType ;
    hash_map<int, std::vector<T> >  *attrib_data;
    int size ;
  public:
    typedef std::vector<T> containerType ;
    const_dstoreVec() { setRep(new storeType) ; }
    const_dstoreVec(const_dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec(dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dstoreVec() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_dstoreVec<T> & operator=(dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_dstoreVec<T> & operator=(const_dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dstoreVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }

    entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> elem(int indx) const {
      std::vector<T>   newVec;
      hash_map<int, std::vector<T> > :: const_iterator   ci;
      ci = attrib_data->find(indx);
      if( ci != attrib_data->end()){
        newVec = ci->second;
        return( newVec );
      } else {
        cout << "Error: Out of range data " << endl;
        exit(0);
      }
      return newVec;
    }

    int vecSize() const { return size ; }
		       
    std::vector<T>  operator[](int indx) const {
      return( elem(indx) );
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  } ;

  //*************************************************************************/

  template<class T> 
  const_dstoreVec<T>::~const_dstoreVec<T>() { }

  //*************************************************************************/

  template<class T> 
  void const_dstoreVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      attrib_data = p->get_attrib_data() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  //************************************************************************/

  template<class T> 
  store_instance::instance_type const_dstoreVec<T>::access() const
  { return READ_ONLY ; }

  //************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_dstoreVec<T> &t)
  { return t.Print(s) ; }

  //*************************************************************************/

  template <class T> 
  storeRepP dstoreVecRepI<T>::remap(const Map &m) const {
    dstoreVec<T> s ;

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;

    return s.Rep() ;
  }

  //*************************************************************************/

  template <class T> 
  void dstoreVecRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_dstoreVec<T> s(st) ;

    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i].clear();
      for(int j=0;j<sz;++j)
        attrib_data[i].push_back(s[i][j]) ;
    } ENDFORALL ;

  }

  //*************************************************************************/

  template <class T> 
  void dstoreVecRepI<T>::gather(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dstoreVec<T> s(st) ;

    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]];
    } ENDFORALL ;

  }

  //*************************************************************************/

  template <class T> 
  void dstoreVecRepI<T>::scatter(const Map &m, storeRepP &st, const entitySet &context)
  {
    const_dstoreVec<T> s(st) ;

    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[m[i]] = s[i];
    } ENDFORALL ;
  }

  //*************************************************************************/
  //**************************************************************************/

  template <class T> 
  int dstoreVecRepI<T>::pack_size( const entitySet &eset) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }


  //*************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  void dstoreVecRepI<T>::StringVal( const int &entity, const int &ivec, 
                                    std::string &memento)
  {
    std::ostringstream oss;
    std::hash_map<int,std::vector<T> > :: const_iterator iter;

    iter = attrib_data.find(entity);
    if( iter != attrib_data.end() )
      oss << attrib_data[entity][ivec] << endl;

    memento = oss.str();
  }

  template <class T> 
  int dstoreVecRepI<T>::get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset) 
  {
    std::ostringstream oss;
    entitySet::const_iterator ci;
    std::vector<T>   newvec;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++)
        oss << attrib_data[*ci][ivec] << " ";
      oss << endl;
    }

    std::string memento = oss.str();
    return ( memento.length());


  }
#endif
  
  //*************************************************************************/

  template <class T> 
  int dstoreVecRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset) 
  {
    return (sizeof(T)*eset.size()*size);
  }
  
  //*************************************************************************/
  template <class T> 
  int dstoreVecRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset) 
  {
    entitySet                :: const_iterator ci;
    std::hash_map<int,std::vector<T> > :: const_iterator iter;
    int       arraySize =0, numContainers = 0;
    std::vector<T> newVec;

    typedef data_schema_traits<T> schema_traits; 

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newVec = iter->second;
        for( int ivec = 0; ivec < size; ivec++){
             typename schema_traits::Converter_Type cvtr( newVec[ivec] );
             arraySize += cvtr.getSize();
        }
      }
    }
    numContainers =  size*eset.size();

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) +
           numContainers*sizeof(int) );
  }
  //*************************************************************************/
  template <class T> 
  void dstoreVecRepI<T>::pack(void * ptr, int &loc, int &size, const entitySet &eset ) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    packdata( traits_type, ptr, loc, size, eset);
  }

  //*************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  void dstoreVecRepI<T>::packdata( DEFAULT_CONVERTER c, void *outbuf,
                                   int &position, int outcount,
                                   const entitySet &eset )
  {
    std::string memento;
    int  bufSize;

    FORALL(eset,ii){
      for( int ivec = 0; ivec < size; ivec++){
        StringVal( ii, ivec, memento );
        bufSize = memento.length();
        MPI_Pack( &bufSize, 1, MPI_INT, outbuf, outcount,
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack( &memento[0], bufSize, MPI_BYTE, outbuf, outcount,
                  &position, MPI_COMM_WORLD) ;
      }
    }ENDFORALL ;
  }
#endif

  //*************************************************************************/
  template <class T> 
  void dstoreVecRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position, 
                                   int outcount, const entitySet &eset ) 
  {
    int incount = size*eset.size();

    std::vector<T>   inbuf;
    entitySet :: const_iterator   ci;
    hash_map<int,std::vector<T> >::const_iterator iter;

    for( ci = eset.begin(); ci != eset.end(); ++ci){
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        inbuf = iter->second;
        MPI_Pack( &inbuf[0], size*sizeof(T), MPI_BYTE, outbuf, outcount, &position, 
                  MPI_COMM_WORLD) ;
      } 
    }
  }

  //************************************************************************/
  template <class T> 
  void dstoreVecRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, int &position, 
                                   int outcount, const entitySet &eset ) 
  {
    entitySet :: const_iterator   ci;

    hash_map<int,std::vector<T> >::const_iterator iter;
    typedef data_schema_traits<T> schema_traits; 

    typedef typename schema_traits::Converter_Base_Type dtype;

    //---------------------------------------------------------------------------
    // Get the maximum size of container 
    //---------------------------------------------------------------------------

    vector<T>  newVec;
    int stateSize, maxStateSize=0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newVec = iter->second;
        for( int ivec = 0; ivec < size; ivec++){
             typename schema_traits::Converter_Type cvtr( newVec[ivec] );
             stateSize  = cvtr.getSize();
             maxStateSize = max( maxStateSize, stateSize);
        }
      }
    }

    int typesize = sizeof(dtype);

    vector<dtype> inbuf(maxStateSize);

    int incount;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newVec = iter->second;
        for( int ivec = 0; ivec < size; ivec++){
             typename schema_traits::Converter_Type cvtr( newVec[ivec] );
             cvtr.getState( &inbuf[0], stateSize);
             maxStateSize = max( maxStateSize, stateSize);

             MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
                      MPI_COMM_WORLD);

             incount =  stateSize*typesize;
             MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
                      MPI_COMM_WORLD) ;
        }
      }
    }

  }

  //************************************************************************/

  template <class T> 
  void dstoreVecRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata( traits_type, ptr, loc, size, seq);

  }

  //************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER

  template <class T>
  void dstoreVecRepI<T>::unpackdata( DEFAULT_CONVERTER c, void *inbuf,
                                     int &position,  int &insize,
                                     const sequence &seq)
  {

    int   outcount, offset;
    sequence:: const_iterator ci;
    entitySet eset(seq);

    vector<char> outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        MPI_Unpack( inbuf, insize, &position, &outcount, 1,
                    MPI_INT, MPI_COMM_WORLD) ;
        if( outcount > outbuf.size() ) outbuf.resize(outcount);

        MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount,
                    MPI_BYTE, MPI_COMM_WORLD) ;
        istringstream iss(&outbuf[0]);
        iss >> attrib_data[*ci][ivec];
      }
    }

  }
#endif
  //************************************************************************/
  
  template <class T> 
  void dstoreVecRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position, 
                                     int &insize, const sequence &seq)
  {

    int   outcount;
    std::vector<T>  outbuf(size);
    sequence :: const_iterator ci;

    if( size == 0) {
      cout << "Error: Size of vector is unknown " << endl;
      return;
    }

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      attrib_data[*ci].resize(size);
      outcount = size*sizeof(T);
      MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount, 
                  MPI_BYTE, MPI_COMM_WORLD) ;
      for( int i = 0; i < size; i++) 
        attrib_data[*ci][i] = outbuf[i];
    }

  }

  //************************************************************************/
  template <class T> 
  void dstoreVecRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, int &position, 
                                     int &insize, const sequence &seq)
  {
    sequence :: const_iterator ci;

    //------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------------
    hash_map<int,vector<T> >::const_iterator   iter;
    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);
    int  stateSize, outcount;

    if( size == 0) {
      cout << "Error: Size of vector is unknown " << endl;
      return;
    }

    std::vector<dtype> outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( attrib_data[*ci].size() < size) {
        attrib_data[*ci].resize(size);
      }
      for( int ivec = 0; ivec < size; ivec++) {
        outcount = sizeof(int);
        MPI_Unpack(inbuf, insize, &position, &stateSize, 1, MPI_INT, 
                   MPI_COMM_WORLD) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);

        outcount = stateSize*typesize;
        MPI_Unpack(inbuf, insize, &position, &outbuf[0], outcount, MPI_BYTE, 
                   MPI_COMM_WORLD) ;

        typename schema_traits::Converter_Type  cvtr(attrib_data[*ci][ivec]);
        cvtr.setState( &outbuf[0], stateSize);
      }
    }

  }

  //************************************************************************/

  template<class T> 
  void dstoreVecRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet  eset, ecommon;

    Loci::HDF5_ReadDomain(group_id, eset);
    Loci::HDF5_ReadVecSize(group_id, &size);

    ecommon = eset & user_eset;

    allocate( ecommon );
    hdf5read( group_id, traits_type, eset, ecommon);
  }

  //************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                     entitySet &eset, entitySet &user_eset )
  {

    hsize_t  dimension;

    hid_t vDatatype  = H5Tcopy(H5T_NATIVE_CHAR);
    hid_t vDataset   = H5Dopen(group_id,"VariableData");
    hid_t vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    vector<char> ibuf(dimension);
    H5Dread(vDataset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    int vsize = get_size();

    entitySet :: const_iterator ci;
    istringstream iss(ibuf);
    for( ci = eset.begin(); ci != eset.end(); ++ci){
      attrib_data[*ci].resize(vsize);
      for(int i=0;i<vsize;++i,++p)
        iss >> attrib_data[*ci][vsize];
    }

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );
  }

#endif
  //************************************************************************/

  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( hid_t group_id, IDENTITY_CONVERTER c, 
                                     entitySet &eset, entitySet &user_eset )
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet::const_iterator ci;

    //------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //------------------------------------------------------------------------
    store<unsigned> offset;
    offset.allocate(eset);

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += size;
    }

    //------------------------------------------------------------------------
    // Read the data now ....
    //------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    T   *data;

    dimension[0] = arraySize;

    typedef data_schema_traits<T> traits_type;

    hid_t mDataspace = H5Screate_simple(rank, dimension, NULL); 
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    DatatypeP dtype = traits_type::get_type() ;
    hid_t vDatatype  = dtype->get_hdf5_type();
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  size;

      data = new T[count[0]];

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        attrib_data[i].clear();
        for( int m = 0; m < size; m++) 
          attrib_data[i].push_back( data[indx++] );
      }

      delete[] data;
    }
  }

  //************************************************************************/

  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                     entitySet &eset, entitySet &user_eset )
  {

    hsize_t  dimension;
    hid_t    vDataspace, vDataset, vDatatype;

    entitySet::const_iterator ci;

    int rank  = 1;
    int vsize = get_size();

    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"SubContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    int maxStateSize = *std::max_element( ibuf.begin(), ibuf.end() );

    //---------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------

    store< unsigned int >   offset;
    dmultiStore<int>  subcontainer;
    offset.allocate( eset );

    int vecsize;
    int arraySize = 0;
    int indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      for( int i = 0; i < size; i++)  {
        vecsize    =  ibuf[indx++];
        arraySize  += vecsize;
        subcontainer[*ci].push_back( vecsize );
      }
    }

    //---------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type  dtype;

    typedef data_schema_traits<dtype> traits_type;
    DatatypeP atom_type = traits_type::get_type() ;
    vDatatype = atom_type->get_hdf5_type();

/*
    vector<dtype>  data(arraySize);

    dimension  = arraySize;
    vDataset   = H5Dopen(group_id,"VariableData");
    vDataspace = H5Dget_space(vDataset);
    mDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, &data[0]);

    //-----------------------------------------------------------------------
    // Fill the objects ....
    //-----------------------------------------------------------------------
    vector<dtype> buf(maxStateSize);
   
    entitySet::const_iterator ci;
    hash_map<int, std::vector<T> > ::const_iterator iter;
    T       newObj;

    size_t indx = 0;
    int    offset, bucketSize, bucketID = 0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      attrib_data[*ci].resize(size);
      for( int ivec = 0; ivec < size; ivec++) {
          typename data_schema_traits<T>::Converter_Type cvtr( attrib[i][ivec] );
          bucsize = subcontainer[i][j];
          cvtr.setState( data+indx, bucsize );
          indx += bucsize ;
      }
    }
*/

  }

  //*************************************************************************/
  template<class T> 
  void dstoreVecRepI<T>::writehdf5( hid_t group_id, entitySet& usr_eset) const
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet eset( usr_eset & domain());

    Loci::HDF5_WriteDomain(group_id, eset);
    Loci::HDF5_WriteVecSize( group_id, size);
   
    hdf5write(group_id, traits_output_type, eset );
  }
  //*************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( hid_t group_id, DEFAULT_CONVERTER c, 
                                     const entitySet &eset ) const
  {
    hsize_t dimension[1];
    int rank = 1;

    std::ostringstream oss;
    entitySet::const_iterator ci;
    hash_map<int, std::vector<T> > ::const_iterator iter;
    std::vector<T>   newvec;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newvec = iter->second;
        for( int i = 0; i < size; i++)
          oss << newvec[i] << " " ;
      }
    }

    std::string memento = oss.str();
    hsize_t arraySize  =  memento.length();

    dimension[0]  =  arraySize+1;
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDatatype  = H5Tcopy(H5T_NATIVE_CHAR);
    hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, H5P_DEFAULT);

    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, memento.c_str());

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

  }
#endif
  //*************************************************************************/

  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( hid_t group_id, IDENTITY_CONVERTER c, 
                                     const entitySet &eset ) const
  {

    hsize_t   dimension[1];
    int       rank = 1;

    //------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //------------------------------------------------------------------------

    int arraySize  = size*eset.size();

    T *data =  new T[arraySize];

    entitySet::const_iterator ci;
    hash_map<int, std::vector<T> > ::const_iterator iter;
    std::vector<T>   newvec;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newvec = iter->second;
        for( int i = 0; i < size; i++)
          data[indx++] =  newvec[i];
      }
    }

    //------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //------------------------------------------------------------------------
    typedef data_schema_traits<T> traits_type;

    dimension[0]=  arraySize;

    hid_t dataspace = H5Screate_simple(rank, dimension, NULL);
    DatatypeP dtype = traits_type::get_type() ;
    hid_t datatype  = dtype->get_hdf5_type();
    hid_t dataset   = H5Dcreate(group_id, "VariableData", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    delete [] data;

  }

  //*************************************************************************/

  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, 
                                     const entitySet &eset ) const
  {
    hsize_t   dimension;
    hid_t     vDataspace, vDataset, vDatatype;

    int       rank = 1;

    //------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //------------------------------------------------------------------------

    entitySet :: const_iterator ci;
    int   bucketID;

    int *vbucket = new int[size*eset.size()];

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    hash_map<int, std::vector<T> > ::const_iterator iter;
    typedef data_schema_traits<T> schema_traits ;

    std::vector<T>   newVec;

    bucketID = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        newVec = iter->second;
        for( int ivec = 0; ivec < size; ivec++){
             typename schema_traits::Converter_Type cvtr( newVec[ivec] );
             stateSize           = cvtr.getSize();
             vbucket[bucketID++] = stateSize;
             arraySize          += stateSize;
             maxStateSize        = max( stateSize, maxStateSize);
        }
      }
    }

    //-----------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype *data ;

    data =  new dtype[arraySize];

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
           iter = attrib_data.find(*ci);
           if( iter != attrib_data.end() ) {
               newVec = iter->second;
               for( int ivec = 0; ivec < size; ivec++){
                    typename schema_traits::Converter_Type cvtr( newVec[ivec] );
                    cvtr.getState( data+indx, stateSize);
                    indx +=stateSize ;

               }
           }
    }

    //-------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-------------------------------------------------------------------------
    typedef data_schema_traits<dtype> traits_type;
    DatatypeP atom_type = traits_type::get_type() ;
    vDatatype = atom_type->get_hdf5_type();

    dimension =  arraySize;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace,
                           H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    //-----------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------
    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    delete [] data;

    dimension =  size*eset.size();

    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5Tcopy( H5T_NATIVE_INT);
    vDataset   = H5Dcreate(group_id, "SubContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vbucket);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );
    delete [] vbucket;

  }


}
#endif

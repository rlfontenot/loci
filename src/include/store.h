#ifndef STORE_H
#define STORE_H

#include <Map.h>
#include <store_rep.h>
#include <data_traits.h>
#include <sstream>
#include <hdf5_memento.h>
#include <hdf5_readwrite.h>
#include <mpi.h>
#include <string.h>

namespace Loci {
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;
  
  class Map ;

  template<class T> class storeRepI : public storeRep {
    T *alloc_pointer ;
    T *base_ptr ;
    entitySet store_domain ;
    bool istat ;

#ifdef ALLOW_DEFAULT_CONVERTER
    void StringVal( const int &entity, std::string &memento);
    int  get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    void packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                    const sequence &seq) ;
    void hdf5read( hid_t group, DEFAULT_CONVERTER c,  entitySet &en, entitySet &usr);
    void hdf5write( hid_t group_id, DEFAULT_CONVERTER g,   const entitySet &en) const;
#endif

    void hdf5read( hid_t group, IDENTITY_CONVERTER c, entitySet &en, entitySet &usr);
    void hdf5write( hid_t group_id, IDENTITY_CONVERTER g,  const entitySet &en) const;
    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq) ;

    void hdf5read( hid_t group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);
    void hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, const entitySet &en) const;
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;

  public:
    storeRepI() { alloc_pointer = 0 ; base_ptr = 0;  istat = 0; }
    storeRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ;istat = 0; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~storeRepI()  ;
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
    virtual void readhdf5( hid_t group, entitySet &en) ;
    virtual void writehdf5( hid_t group_id,entitySet& en) const ;
    bool is_static() {return istat; }
    virtual entitySet domain() const ;
    T * get_base_ptr() const { return base_ptr ; }
    
  } ;
  
  template<class T> void storeRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ; int size = ptn.Max()-top+1 ;
      alloc_pointer = new T[size] ;
      base_ptr = alloc_pointer - top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  template<class T> std::ostream &storeRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << std::endl ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }
    

  template<class T> std::istream &storeRepI<T>::Input(std::istream &s) {
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
      s >> base_ptr[ii] ;
    } ENDFORALL ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //*************************************************************************/

  template<class T>  storeRepI<T>::~storeRepI<T>() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }
    
  template<class T>  entitySet storeRepI<T>::domain() const {
    return store_domain ;
  }

  template<class T>
  storeRep *storeRepI<T>::new_store(const entitySet &p) const {
    return new storeRepI<T>(p)  ;
  }

  template<class T> store_type storeRepI<T>::RepType() const {
    return STORE ;
  }

  template<class T> class store : public store_instance {
    typedef storeRepI<T> storeType ;
    T* base_ptr ;
  public:
    typedef T containerType ;
    store() { setRep(new storeType); }
    store(store &var) { setRep(var.Rep()) ; }
    store(storeRepP rp) { setRep(rp) ; }
    
    virtual ~store() ;
    virtual void notification() ;

    store<T> & operator=(store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    
    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
    //    operator storeRepP() { return Rep() ; }
    T &elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
    const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
    
    T &operator[](int indx) { return elem(indx); }
    const T&operator[](int indx) const { return elem(indx); }
    
  } ;

  template<class T> store<T>::~store<T>() { }
    
  template<class T> void store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }
 
  template<class T> inline std::ostream & operator<<(std::ostream &s, const store<T> &t)
  { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, store<T> &t)
  { return t.Input(s) ; }



  template<class T> class const_store : public store_instance {
    typedef storeRepI<T> storeType ;
    const T * base_ptr ;
  public:
    typedef T containerType ;
    const_store() { setRep(new storeType) ; }
    const_store(store<T> &var) { setRep(var.Rep()) ; }
    const_store(const_store &var) { setRep(var.Rep()) ; }
    const_store(storeRepP rp) { setRep(rp) ; }

    virtual ~const_store() ;
    virtual void notification() ;


    virtual instance_type access() const  ;
        
    const_store<T> & operator=(const_store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_store<T> & operator=(store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    bool is_static() {return 0 ; }
   

    const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
            
    const T&operator[](int indx) const { return elem(indx); }

  } ;

  template<class T> const_store<T>::~const_store<T>() { }
    
  template<class T> void const_store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
  const_store<T>::access() const { return READ_ONLY; }
        
  template<class T> storeRepP storeRepI<T>::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    store<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
      
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T> void storeRepI<T>::copy(storeRepP &st, const entitySet &context)  {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr ==0)) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }
  template<class T> void storeRepI<T>::gather(const Map &m, storeRepP &st,
                                              const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T> void storeRepI<T>::scatter(const Map &m, storeRepP &st,
                                               const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }
  
  //*******************************************************************/
  template <class T> int storeRepI<T>::pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }

  //*******************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  int storeRepI<T>::get_mpi_size( DEFAULT_CONVERTER c,
                                  const entitySet &eset)
  {
    std::ostringstream oss;
    int  size;

    FORALL(eset,ii){
      oss << base_ptr[ii] << std::endl ;
    }ENDFORALL ;

    std::string memento = oss.str();

    size = memento.length() + eset.size()*sizeof(int);

    return( size );
  }
#endif

  //*******************************************************************/

  template <class T>
  int storeRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                  const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }

  //*******************************************************************/
  template <class T>
  int storeRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                  const entitySet &eset)
  {

    int       size, numBytes;
    entitySet  ecommon;
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    ecommon = eset;

    T   obj;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      obj  = base_ptr[*ci];
      typename converter_traits::Converter_Type cvtr(obj);
      size      = cvtr.getSize();
      numBytes += size*sizeof(typename converter_traits::Converter_Base_Type) ;
    }

    numBytes  += ecommon.size()*sizeof(int);
    return(numBytes) ;
  }
  //*******************************************************************/

  template <class T> 
  void storeRepI<T>::pack( void *outbuf, int &position, int &size, 
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
  void storeRepI<T>::StringVal( const int &entity, std::string &memento)
  {
    std::ostringstream oss;

    oss << base_ptr[entity] << endl;

    memento = oss.str();
  }

  template <class T> 
  void storeRepI<T>::packdata(DEFAULT_CONVERTER c, void *outbuf,
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
  void storeRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                              int &position, int outcount,
                              const entitySet &eset )  
  {
    for( int i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      MPI_Pack( &base_ptr[begin], t*sizeof(T), MPI_BYTE, outbuf,outcount, 
                &position, MPI_COMM_WORLD) ;
    }
  }

  //*******************************************************************/

  template <class T>
  void storeRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
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
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
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
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci]);
      cvtr.getState( inbuf, stateSize);

      incount =  sizeof(int);
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
  void storeRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size, seq); 
  }

  //*********************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void storeRepI<T>::unpackdata(DEFAULT_CONVERTER c, void *inbuf,
                                int &position,  int insize,
                                const sequence &seq) 
  {
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
  void storeRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                int &position,  int insize,
                                const sequence &seq) 
  {

    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                      sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                    t*sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      }
    }
  }
  //*********************************************************************/
  template <class T> 
  void storeRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
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

      typename converter_traits::Converter_Type cvtr( base_ptr[*ci] );
      cvtr.setState( outbuf, stateSize);
      delete [] outbuf;
    }

  }

  //**********************************************************************/
  template<class T> 
  void storeRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset)
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    HDF5_ReadDomain(group_id, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );

    hdf5read( group_id, traits_type, eset, ecommon);

  }

  //**************************************************************************/

  template<class T> 
  void storeRepI<T>::writehdf5( hid_t group_id, entitySet& eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet     ecommon;
    ecommon = store_domain & eset;

    if( ecommon.size() < 1) return; 
    HDF5_WriteDomain(group_id, ecommon);

    hdf5write(group_id, traits_output_type, ecommon );

  }

  //**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void storeRepI<T> :: hdf5write( hid_t group_id, DEFAULT_CONVERTER g,
                                  const entitySet &en) const
  {

    int rank = 1;
    hsize_t dimension[1];
    std::ostringstream oss;

    FORALL(en,ii) {
      oss << base_ptr[ii] << std::endl ;
    }ENDFORALL ;

    std::string memento = oss.str();
    hsize_t size  =  memento.length();
    dimension[0]  =  size+1;

    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_CHAR;
    hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
                                 vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             memento.c_str());

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);

  }
#endif

  //**********************************************************************/
  template <class T> 
  void storeRepI<T> :: hdf5write( hid_t group_id, IDENTITY_CONVERTER c, 
                                  const entitySet &eset)  const
  {

    int arraySize =  eset.size(); 

    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    int rank = 1;
    hsize_t  dimension = eset.size();

    entitySet :: const_iterator ci;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);

    T *data = new T[eset.size()];
    int indx =0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
         data[indx++] = base_ptr[*ci];

    hid_t cparms   = H5Pcreate (H5P_DATASET_CREATE);
    hid_t vDataset = H5Dcreate(group_id, "VariableData", vDatatype,
                               vDataspace, cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    delete [] data;

  }
  //*********************************************************************/

  template <class T> 
  void storeRepI<T>::hdf5write( hid_t group_id,USER_DEFINED_CONVERTER g, 
                                const entitySet &eset) const
  {   
    entitySet :: const_iterator ci;
    hid_t   vDataspace, vDataset, vDatatype;

    //---------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //---------------------------------------------------------------
    typedef data_schema_traits<T> schema_traits ;

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
      stateSize    = cvtr.getSize();
      arraySize   += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
    }

    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> data(arraySize);

    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
      cvtr.getState( &data[0]+indx, stateSize);
      indx +=stateSize ;
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
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    //--------------------------------------------------------------------
    // Write Container 
    //--------------------------------------------------------------------
    dimension    = eset.size();
    std::vector<int> vbucket(eset.size());

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
      vbucket[indx++] =  cvtr.getSize();
    }

    vDatatype  = H5T_NATIVE_INT;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "ContainerSize",
                           vDatatype, vDataspace, H5P_DEFAULT);

    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vbucket[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
  };


  //**********************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template<class T>
  void  storeRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                  entitySet &eset, entitySet &usr_eset)
  {

    hsize_t  dimension;
    hid_t vDatatype  = H5T_NATIVE_CHAR;
    hid_t vDataset   = H5Dopen(group_id,"VariableData");
    hid_t vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims (vDataspace, &dimension, NULL);

    vector<char> cbuf(dimension);

    H5Dread(vDataset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    entitySet :: const_iterator ci;

    std::istringstream iss(&ibuf[0]);

    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      iss >> base_ptr[*ci];

    H5Dclose(vDataset  );
    H5Sclose(vDataspace);

  }
#endif

  //*********************************************************************/

  template<class T>
  void  storeRepI<T> :: hdf5read(hid_t group_id, IDENTITY_CONVERTER c, 
                                 entitySet &eset, entitySet &usr_eset)
  {
    int      rank = 1, indx = 0;
    hid_t    mDataspace, vDataspace, vDataset, vDatatype;
    hsize_t  dimension;
    entitySet::const_iterator  ci;

    store<int>  offset;
    offset.allocate( eset );

    indx     = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      offset[*ci] = indx++;

    //--------------------------------------------------------------------
    // Read the data now ....
    //--------------------------------------------------------------------
    int num_intervals = usr_eset.num_intervals();

    if( num_intervals == 0) {
      std::cout << "Warning: Number of intervals are zero : " << endl;
      return;
    }

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    vDatatype = dtype->get_hdf5_type();


    dimension  = eset.size();
    mDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    T  *data;
    int preallocated = 0;
    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] += 1;

      if( count[0] > preallocated) {
          delete [] data;
          preallocated = count[0];
          data = new T[preallocated];
      }

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride,
                          count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride,
                          count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace,
               H5P_DEFAULT, data);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) 
        base_ptr[i] = data[indx++];
    }

    if( preallocated) delete [] data;

    H5Sclose( mDataspace );
    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

  }
  //*********************************************************************/

  template<class T>
  void  storeRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                  entitySet &eset, entitySet &usr_eset)
  {

    hsize_t  dimension;
    size_t   indx = 0, arraySize;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;
    entitySet::const_iterator  ci;

    //---------------------------------------------------------------
    // Size of each Bucket ....
    //---------------------------------------------------------------
    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"ContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    std::vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);
    H5Tclose(vDatatype );
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
        typename data_schema_traits<T>::Converter_Type cvtr( base_ptr[i]);
        cvtr.setState( &data[0]+indx, container[i] );
        indx += container[i];
      }
    }

    H5Tclose(vDatatype );
    H5Dclose(vDataset  );
    H5Sclose(vDataspace);
    H5Sclose(mDataspace);

  }
  //*********************************************************************/
}

#endif

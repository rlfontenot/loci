#ifndef STORE_H
#define STORE_H

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

#include <hdf5CC/H5cpp.h>
#include <hdf5_traits.h>
#include <Map.h>
#include <hdf5_write_template.h>
#include <Tools/intervalSet.h>
#include <algorithm>
#include <functional>

#include <mpi.h>
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

    int get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, 
                  int size, const entitySet &e) ;
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, 
                  int size, const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, 
                  int size, const entitySet &e) ;

    void unpackdata(DEFAULT_CONVERTER c,     void *ptr, int &loc, 
                    int size, const sequence &seq) ;
    void unpackdata(IDENTITY_CONVERTER c,    void *ptr, int &loc, 
                    int size, const sequence &seq) ;
    void unpackdata(USER_DEFINED_CONVERTER c,void *ptr, int &loc, 
                    int size, const sequence &seq) ;

    void  hdf5read( H5::Group group, DEFAULT_CONVERTER c,      entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void hdf5write( H5::Group group, DEFAULT_CONVERTER g,      const entitySet &en) const;
    void hdf5write( H5::Group group, IDENTITY_CONVERTER g,     const entitySet &en) const;
    void hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, const entitySet &en) const;


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
    virtual void readhdf5( H5::Group group, entitySet &en) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
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

  //**************************************************************************

  template<class T> 
  void storeRepI<T>::readhdf5( H5::Group group, entitySet &user_eset)
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    HDF5_ReadDomain(group, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );

    hdf5read( group, traits_type, eset, ecommon);
  }

  //**************************************************************************

  template<class T> 
  void storeRepI<T>::writehdf5( H5::Group group,entitySet& en) const
  {
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    hdf5write(group, traits_output_type, en );
  }

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
    store(storeRepP &rp) { setRep(rp) ; }
    
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
    const_store(storeRepP &rp) { setRep(rp) ; }

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

  //*******************************************************************
  
  template <class T> int storeRepI<T>::pack_size( const entitySet &eset) 
  {
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }
  
  //*******************************************************************
  template <class T>
  int storeRepI<T>::get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }

  //*******************************************************************

  template <class T>
  int storeRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }

  //*******************************************************************
  template <class T>
  int storeRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {

    int size, numBytes = 0;
    entitySet  ecommon;

    entitySet :: const_iterator ci;

    typedef hdf5_schema_converter_traits<T> converter_traits;

    ecommon =  domain() & eset;

    T   obj;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
         obj  = base_ptr[*ci];
         Memento<T> memento( obj );
         size      = memento.getSize();
         numBytes += size*sizeof(typename converter_traits::memento_type );
    }
    numBytes  += ecommon.size()*sizeof(int);
    return(numBytes) ;
  }
  //*******************************************************************

  template <class T> 
  void storeRepI<T>::pack( void *outbuf, int &position, int &size, 
                            const entitySet &eset )  
  {
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    packdata( traits_type, outbuf, position, size, eset);

  }
  //*******************************************************************
  template <class T> 
  void storeRepI<T>::packdata(DEFAULT_CONVERTER c, void *outbuf, int &position, 
                              int outcount,  const entitySet &eset )  
  {
     IDENTITY_CONVERTER ic;
     packdata(ic, outbuf, position, outcount, eset);
  }
  //*******************************************************************

  template <class T> 
  void storeRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf, int &position, 
                              int outcount,  const entitySet &eset )  
  {
    entitySet :: const_iterator ei;
    entitySet  ecommon;

    ecommon = eset & store_domain;

    T *inbuf =  new T[ecommon.size()];

    int indx = 0;
    for( ei = ecommon.begin(); ei != ecommon.end(); ++ei)
         inbuf[indx++] = base_ptr[*ei];

    int incount =  ecommon.size()*sizeof(T);
    MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position, 
             MPI_COMM_WORLD) ;

    delete [] inbuf;
  }

  //*******************************************************************

  template <class T>
  void storeRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                               int &position, int outcount, 
                               const entitySet &eset )

  {

    entitySet :: const_iterator ci;
    entitySet  ecommon;

    ecommon = store_domain & eset;
//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the
// container for allocation purpose
//-----------------------------------------------------------------------------
    size_t  incount = 0;
    int     stateSize, maxStateSize = 0;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
         Memento<T> memento( base_ptr[*ci] );
         stateSize    = memento.getSize();
         incount     += stateSize;
         maxStateSize = max( maxStateSize, stateSize );
    }

    typedef hdf5_schema_converter_traits<T> converter_traits;
    typename converter_traits::memento_type *inbuf;

    int typesize = sizeof(typename converter_traits::memento_type);

    inbuf = new typename converter_traits::memento_type[maxStateSize];

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
         Memento<T> memento( base_ptr[*ci]);
         memento.getState( inbuf, stateSize);

         incount =  sizeof(int);
         MPI_Pack(&stateSize,incount, MPI_BYTE, outbuf, outcount,&position,
                  MPI_COMM_WORLD);

         incount =  stateSize*typesize;
         MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position, 
                  MPI_COMM_WORLD) ;
    }

    delete [] inbuf;
  }

  //*******************************************************************

  template <class T> 
  void storeRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size, seq); 
  }

  //****************************************************************************

  template <class T> 
  void storeRepI<T>::unpackdata(DEFAULT_CONVERTER c, void *inbuf, int &position, 
                                 int insize, const sequence &seq) 
  {
     IDENTITY_CONVERTER ic;
     unpackdata(ic, inbuf, position, insize, seq);
  }
  //****************************************************************************

  template <class T> 
  void storeRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf, int &position, 
                                 int insize, const sequence &seq) 
  {

    int outcount = seq.size();
    T  *outbuf = new T[outcount];
    MPI_Unpack(inbuf, insize, &position, outbuf, outcount*sizeof(T), 
               MPI_BYTE, MPI_COMM_WORLD) ;

    int indx, jndx;
    for(int i = 0; i < seq.num_intervals(); ++i) {
        jndx = 0;
        if(seq[i].first > seq[i].second) {
           for(indx = seq[i].first; indx >= seq[i].second; --indx)
               base_ptr[indx] =  outbuf[jndx++];
        } else {
           for(indx = seq[i].first; indx <= seq[i].second; ++indx)
               base_ptr[indx] =  outbuf[jndx++];
        }
    }

    delete [] outbuf;
  }

  //****************************************************************************
  template <class T> 
  void storeRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                 int &position, int insize, const sequence &seq) 
  {

    sequence :: const_iterator ci;

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------
    int  stateSize, outcount;

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *outbuf;
    int typesize = sizeof(typename converter_traits::memento_type);

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
        if( !store_domain.inSet( *ci ) ) {
             cout << "Warning: entity not in entitySet  : " << *ci << endl;
             continue;
        }
        outcount = sizeof(int);
        MPI_Unpack(inbuf, insize, &position, &stateSize, outcount, 
                   MPI_BYTE, MPI_COMM_WORLD) ;

        outbuf = new typename converter_traits::memento_type[stateSize];

        outcount = stateSize*typesize;
        MPI_Unpack(inbuf, insize, &position, outbuf, outcount, 
                   MPI_BYTE, MPI_COMM_WORLD) ;

        Memento<T> memento( base_ptr[*ci] );
        base_ptr[*ci] = memento.setState( outbuf, stateSize);
        delete [] outbuf;
    }
  }
  //****************************************************************************
  template <class T> 
  void storeRepI<T> :: hdf5write( H5::Group group, DEFAULT_CONVERTER g,
                                  const entitySet &en) const
  {

    int rank = 1;

    hsize_t dimension[1];
    std::ostringstream oss;
    oss << '{' << en << std::endl ;

    FORALL(en,ii) {
      oss << base_ptr[ii] << std::endl ;
    }ENDFORALL ;

    oss << '}' << std::endl ;
   
    std::string memento = oss.str();
    hsize_t size  =  memento.length();
    dimension[0]  =  size+1;

    try {
      H5::DataSpace dataspace( rank, dimension );
      H5::DataSet dataset = group.createDataSet( "store", 
                                                 H5::PredType::NATIVE_CHAR,
	                                              dataspace);
      dataset.write( memento.c_str(), H5::PredType::NATIVE_CHAR );
    }

    catch( H5::HDF5DatasetInterfaceException error  ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error ) { error.printerror(); }
   
  }

  //****************************************************************************

  template <class T> 
  void storeRepI<T> :: hdf5write( H5::Group group, IDENTITY_CONVERTER g, 
                                  const entitySet &eset) const
  {
    entitySet  ecommon;

    ecommon = store_domain & eset;

    //write out the domain   
    HDF5_WriteDomain(group, eset);

    int arraySize =  eset.size(); 

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------
    T    *data;

 	 data =  new T[arraySize];

    entitySet :: const_iterator ci;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) 
         data[indx++] =  base_ptr[*ci];

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------
    typedef hdf5_schema_traits<T> traits_type;

    int rank = 1;
    hsize_t  dimension[1];

    dimension[0] =  arraySize;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = traits_type::get_type();
      H5::DataSet   vDataset  = group.createDataSet("VariableData", 
                                                    vDatatype, vDataspace);
      vDataset.write( data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

//-----------------------------------------------------------------------------
// Clean up
//-----------------------------------------------------------------------------
    delete [] data;
  };

  //****************************************************************************

  template <class T> 
  void storeRepI<T> :: hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, 
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


    //write out the domain   
    HDF5_WriteDomain(group, eset);

    entitySet :: const_iterator ci;

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         Memento<T> memento( base_ptr[*ci] );
         stateSize    = memento.getSize();
         arraySize   += stateSize;
         maxStateSize = max( maxStateSize, stateSize );
    }

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *data, *buf;

 	 data =  new typename converter_traits::memento_type[arraySize];
 	 buf  =  new typename converter_traits::memento_type[maxStateSize];

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         Memento<T> memento( base_ptr[*ci] );
         memento.getState( buf, stateSize);
         for( int i = 0; i <  stateSize; i++) {
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
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", 
                                                     vDatatype, vDataspace);
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
    typedef hdf5_schema_traits<T>   schema_traits;

    dimension[0]    = eset.size();
    int *vbucket =  new int[ eset.size() ];

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         Memento<T> memento( base_ptr[*ci] );
         vbucket[indx++] = memento.getSize();
    }

    H5::DataType  frameDatatype = H5::PredType::NATIVE_INT;
    H5::DataSpace frameDataspace( rank, dimension );
    H5::DataSet   frameDataset = group.createDataSet( "ContainerSize", 
                                                       frameDatatype, 
                                                       frameDataspace);
    frameDataset.write( vbucket, frameDatatype );

  };


  //*****************************************************************************
  template<class T>
  void  storeRepI<T> :: hdf5read( H5::Group group, DEFAULT_CONVERTER c, 
                                  entitySet &eset, entitySet &usr_eset)
  {
     cout << " Warning: Default converter for store not implemented " << endl;
  }

  //*****************************************************************************

  template<class T>
  void  storeRepI<T> :: hdf5read( H5::Group group, IDENTITY_CONVERTER c, 
                                  entitySet &eset, entitySet &usr_eset)
  {

    int      rank = 1, indx = 0;
    hsize_t  dimension[1];
    entitySet::const_iterator  ci;

   store<int>  offset;
   offset.allocate( eset );

   indx     = 0;
   for( ci = eset.begin(); ci != eset.end(); ++ci)
        offset[*ci] = indx++;

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------
   int num_intervals = usr_eset.num_intervals();

   if( num_intervals == 0) {
       cout << "Warning: Number of intervals are zero : " << endl;
       return;
   }

   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

   typedef hdf5_schema_traits<T> traits_type;

   T  *data;

   dimension[0] = eset.size();
   H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
   H5::DataSpace vDataspace(rank, dimension);

   H5::DataType  vDatatype = traits_type::get_type();
	H5::DataSet   vDataset   = group.openDataSet( "VariableData");

   hssize_t  start_mem[] = {0};  // determines the starting coordinates.
   hsize_t   stride[]    = {1};  // which elements are to be selected.
   hsize_t   block[]     = {1};  // size of element block;
   hssize_t  foffset[]   = {0};  // location (in file) where data is read.
   hsize_t   count[]     = {0};  // how many positions to select from the dataspace

   for( int k = 0; k < num_intervals; k++) {
        count[0] = 0;
        for( int i = it[k].first; i <= it[k].second; i++)
             count[0] += 1;

        data = new T[count[0]];

        foffset[0] = offset[it[k].first];

        mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
        vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        vDataset.read( data, vDatatype, mDataspace, vDataspace);

        indx = 0;
        for( int i = it[k].first; i <= it[k].second; i++) 
             base_ptr[i] = data[indx++];

        delete[] data;
   }

  }
  //*****************************************************************************

  template<class T>
  void  storeRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, 
                                  entitySet &eset, entitySet &usr_eset)
  {
    hsize_t  dimension[1];
    size_t   indx = 0, arraySize;
    int      rank = 1;
    entitySet::const_iterator  ci;

    //-------------------------------------------------------------------------
    // Size of each Bucket ....
    //--------------------------------------------------------------------------

    H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   bDataset   = group.openDataSet( "ContainerSize");
    H5::DataSpace bDataspace = bDataset.getSpace();

    bDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

	 bDataset.read( ibuf, H5::PredType::NATIVE_INT );

    int maxBucketSize = *max_element( ibuf, ibuf + (int)dimension[0] );

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
   for( ci = eset.begin(); ci != eset.end(); ++ci) {
        offset[*ci] = indx;
        indx       += container[*ci];
   }

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------
   int num_intervals = usr_eset.num_intervals();
   if( num_intervals == 0) {
       cout << "Warning: Number of intervals are zero : " << endl;
       return;
   }
 
   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

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
        for( int i = it[k].first; i <= it[k].second; i++)
             count[0] +=  container[i];

        data = new typename converter_traits::memento_type[count[0]];

        foffset[0] = offset[it[k].first];

        mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
        vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        vDataset.read( data, vDatatype, mDataspace, vDataspace);

        indx = 0;
        for( int i = it[k].first; i <= it[k].second; i++) {
             Memento<T> memento( base_ptr[i] );
             for( int k = 0; k < container[i]; k++) {
                  buf[k] = data[indx++];
             }
             base_ptr[i] = memento.setState( buf, container[i]);
        }

        delete[] data;
   }

   delete[] buf;


  }
  //***************************************************************************


  
}

#endif

#ifndef STOREVEC_H
#define STOREVEC_H 

#include <istream>
#include <ostream>

#include <algorithm>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>
#include <store.h>
#include <DMultiStore.h>

#include <Tools/lmutex.h>
#include <hdf5CC/H5cpp.h>

#include <Map.h>

namespace Loci {

  //**************************************************************************/

  template <class T> struct Scalar {
    T val ;
    Scalar(T v) : val(v) { }
  } ;

  template <class T> Scalar<T> mk_Scalar(T v) { return Scalar<T>(v) ;} 

  //**************************************************************************/
  
  template <class T> class Vect ;
  
  template <class T> class const_Vect {
  public:
    friend class Vect<T> ;
    const T *ptr ;
#ifdef BOUNDS_CHECK
    int size ;
#endif
  public:

    const_Vect(const T *p ,int sz) {
      ptr = p ;
#ifdef BOUNDS_CHECK
      size = sz ;
#endif
    }

    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    operator const T *() const {
      return ptr ;
    }
  } ;

  //**************************************************************************/

  template <class T> class Vect {
  public:
    T *ptr ;
    int size ;
  public:
    Vect() {};
    void setSize( int s ) {
         size = s;
    }

    Vect(T *p ,int sz) {
      ptr = p ;
      size = sz ;
    }

    void operator=(const Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    void operator=(const const_Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    template <class S> void operator=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ = s.val ;
    }

    template <class S> void operator+=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ += s.val ;
    }
      
    template <class S> void operator*=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ *= s.val ;
    }
      
    template <class S> void operator-=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ -= s.val ;
    }
      
    template <class S> void operator/=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ /= s.val ;
    }

    template <class S> void operator+=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    template <class S> void operator-=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    operator T*() {
      return ptr ;
    }

    operator const T *() const {
      return ptr ;
    }
  } ;

  //**************************************************************************/

  template<class T> class storeVecRepI : public storeRep {
    entitySet    store_domain ;
    T           *alloc_pointer, *base_ptr ;
    int          size ;
    lmutex       mutex ;
    bool istat ;
    void hdf5read( H5::Group group, DEFAULT_CONVERTER      c, entitySet &en, entitySet &usr);
    void hdf5read( H5::Group group, IDENTITY_CONVERTER     c, entitySet &en, entitySet &usr);
    void hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void hdf5write( H5::Group group, DEFAULT_CONVERTER g,      const entitySet &en) const; 
    void hdf5write( H5::Group group, IDENTITY_CONVERTER g,     const entitySet &en) const;
    void hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, const entitySet &en) const;

  public:
    storeVecRepI() 
    { alloc_pointer= 0 ; base_ptr = 0 ; size=0; istat = 1 ; }

    storeVecRepI(const entitySet &p) 
    { size = 0; alloc_pointer=0 ; allocate(p) ; istat = 1 ; }
    
    virtual ~storeVecRepI() ;
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
    virtual void readhdf5( H5::Group group, entitySet &en) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    virtual void set_elem_size(int sz) ;
    bool is_static() {return istat; }
    T *get_base_ptr() const { return base_ptr ; }
    int get_size() const { return size ; }
  } ;

  //**************************************************************************/

   template<class T> 
	std::ostream &storeVecRepI<T>::Print(std::ostream &s) const
   {

      s << '{' << domain() << std::endl ;
      s << size << std::endl ;
    
      FORALL(domain(),ii) {
        T * p = base_ptr + ii*size ;
        for(int i=0;i<size;++i,++p)
          s << *p << " " ;
        s << std::endl ;
      }ENDFORALL ;
      s << '}' << std::endl ;

      return s ;
   }

  //*************************************************************************/

  template<class T> 
  std::istream &storeVecRepI<T>::Input(std::istream &s)
  {

      cout << " Commented for time being " << endl;
      exit(0);

  //-------------------------------------------------------------------------
  // Objective : Read the storeVec from the input stream.
  //-------------------------------------------------------------------------
  /*
      char ch ;
    
	// Look for the opening brackets ...
      do ch = s.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
        return s ;
      }

      entitySet e ;
      int sz ;

      s >> e ;               // Read the entitySet intervals.
      s >> sz ;              // Read the size of the vector.

      set_elem_size(sz) ;
      allocate(e) ;

      FORALL(e,ii) {
        T * p = base_ptr + ii*size ;
        for(int i=0;i<size;++i,++p)
          s >> *p ;
      } ENDFORALL ;
    
	// Look for the closing brackets ...
      do ch = s.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
      }
   */
	  
      return s ;
  }

  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::readhdf5( H5::Group group, entitySet &user_eset)
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    //--------------------------------------------------------------------------
    // Read the vector size ...
    //--------------------------------------------------------------------------
    hsize_t dimension[1];

    dimension[0] = 1;

    H5::DataType  datatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   dataset   = group.openDataSet( "VecSize");
    H5::DataSpace dataspace = dataset.getSpace();

    dataspace.getSimpleExtentDims( dimension, NULL);

    dataset.read( &size, H5::PredType::NATIVE_INT );

    //************************************************************************/
    // Only ecommon should be read, but I didn't touch this code
    //************************************************************************/

    entitySet   eset, ecommon;

    HDF5_ReadDomain(group, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );
    hdf5read( group, traits_type, eset, ecommon);

  }

  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::writehdf5( H5::Group group,entitySet &eset) const
  {
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    hdf5write(group, traits_output_type, eset);
  }

  //**************************************************************************/


  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::allocate(const entitySet &ptn) 
  {

  //---------------------------------------------------------------------------
  // Allocation reclaims all previously hold memeory 
  //---------------------------------------------------------------------------

    if(alloc_pointer) delete[] alloc_pointer ;

    alloc_pointer = 0 ;
    base_ptr      = 0 ;

  //---------------------------------------------------------------------------
  // Get the minimum and maximum entity ID from the entitySet and allocate 
  // memory of the size = ( max-min+1). Notice that, if the entityset 
  // contains the entities with ID quite sparse, it will create lots of 
  // unused block of memory. 
  //---------------------------------------------------------------------------

    if(size != 0) {
      fatal(size < 1) ;
      if(ptn != EMPTY) {
        int top = ptn.Min() ; int sza = (ptn.Max()-top+1)*size ;
        alloc_pointer = new T[sza] ;
        base_ptr = alloc_pointer - top*size ;
      }
    }

  //---------------------------------------------------------------------------
  // Domain equals to entitySet provided by the argument.
  //--------------------------------------------------------------------------

    store_domain = ptn ;
  //--------------------------------------------------------------------------
  // Let everybody know about the change in memeory location.
  //--------------------------------------------------------------------------

    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> 
  storeVecRepI<T>::~storeVecRepI<T>() 
  {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  //**************************************************************************/

  template<class T>
  storeRep *storeVecRepI<T>::new_store(const entitySet &p) const 
  {
    storeRep *sp = new storeVecRepI<T>(p)  ;
    sp->set_elem_size(size) ;
    return sp ;
  }

  //*************************************************************************/

  template<class T> 
  store_type storeVecRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //**************************************************************************/

  template<class T> 
  entitySet storeVecRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::set_elem_size(int sz) 
  {
    mutex.lock() ;

    //------------------------------------------------------------------------
    // Change the size of vector held. It will reclaim the memory used before
    // his call. and allocate new one.
    //------------------------------------------------------------------------

    if(size != sz) {
      if(size != 0) {
	cout << " sz = " << sz << "   size =  " << size << endl ;
        warn(size != sz) ;
      }
      size = sz ;
      fatal(sz<1) ;
      allocate(store_domain) ;
    }
    
    mutex.unlock() ;
  }

  //**************************************************************************/
      
  template<class T> class storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    storeVec() {setRep(new storeType) ;}
    storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    storeVec(storeRepP &rp) { setRep(rp) ;}

    virtual ~storeVec() ;
    virtual void notification() ;

    storeVec<T> & operator=(storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    storeVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      Rep()->set_elem_size(size) ;
    }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size ; }

    const entitySet domain() const { return Rep()->domain() ; }

    Vect<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr+(indx*size),size) ; }
    Vect<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr+(indx*size),size) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //**************************************************************************/

  template<class T> 
  storeVec<T>::~storeVec<T>() { }
  
  //*************************************************************************/

  template<class T> 
  void storeVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;

    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size = p->get_size() ;
    }

    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const storeVec<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, storeVec<T> &t)
    { return t.Input(s) ; }

  //**************************************************************************/

  template<class T> class const_storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_storeVec() { setRep(new storeType) ; }
    const_storeVec(const_storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_storeVec() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_storeVec<T> & operator=(storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_storeVec<T> & operator=(const_storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_storeVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }

    int vecSize() const { return size ; }

    const entitySet domain() const { return Rep()->domain() ; }

    const_Vect<T> elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Vect<T>(base_ptr+(indx*size),size) ; }
    const_Vect<T> operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Vect<T>(base_ptr+(indx*size),size) ;
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

  } ;

  //**************************************************************************/

  template<class T> 
  const_storeVec<T>::~const_storeVec<T>() { }

  //**************************************************************************/

  template<class T> 
  void const_storeVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> store_instance::instance_type
    const_storeVec<T>::access() const
    { return READ_ONLY ; }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_storeVec<T> &t)
    { return t.Print(s) ; }

  //**************************************************************************/

  template <class T> 
  storeRepP storeVecRepI<T>::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    storeVec<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  //**************************************************************************/

  template <class T> 
  void storeVecRepI<T>::copy(storeRepP &st, const entitySet &context) {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + i*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }

  //**************************************************************************/

  template <class T> 
  void storeVecRepI<T>::gather(const Map &m, storeRepP &st,
                               const entitySet &context) 
  {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal(base_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + i*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[m[i]][j] ;
    } ENDFORALL ;
  }

  //*************************************************************************/

  template <class T> 
  void storeVecRepI<T>::scatter(const Map &m, storeRepP &st,
                                const entitySet &context) 
  {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal(base_ptr == 0) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + m[i]*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }

  //**************************************************************************/
  template <class T> 
    int storeVecRepI<T>::pack_size( const entitySet &e) {
    int size, M ;
    M = get_size() ;
    size = (sizeof(T) * e.size() * M ) + sizeof(int) ;
    return(size) ;
  }

  //*************************************************************************/
  template <class T> 
    void storeVecRepI<T>::pack( void * ptr, int &loc, int &size, 
				const entitySet &e ) 
    {
      int M = get_size() ;
      MPI_Pack(&M, sizeof(int), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
      for(int i = 0; i < e.num_intervals(); ++i) {
	Loci::int_type indx1 = e[i].first ;
	Loci::int_type stop = e[i].second ;
	T *p = base_ptr + M * indx1 ;
	int t = (stop - indx1 + 1) * M ;
	MPI_Pack(p, t * sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
      }
    }
  
  //**************************************************************************/
  template <class T> void storeVecRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    int init_size = get_size() ;
    int M ;
    MPI_Unpack(ptr, size, &loc, &M, sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
    
    if(init_size < M) {
      fatal(init_size != 0) ;
      set_elem_size(M) ;
    }
    
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
	const Loci::int_type indx1 = seq[i].first ;
	const Loci::int_type stop =  seq[i].second ;
	for(Loci::int_type indx = indx1; indx != stop-1; --indx) {
	  T *p = base_ptr + M * indx ;
	  MPI_Unpack(ptr, size, &loc, p, sizeof(T) * M, MPI_BYTE, MPI_COMM_WORLD) ;
	}
      }
      else {
	Loci::int_type indx1 = seq[i].first ;
	Loci::int_type stop = seq[i].second ;
	T *p = base_ptr + M * indx1 ;
	int t = (stop - indx1 + 1) * M ;
	MPI_Unpack(ptr, size, &loc, p, t * sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ; 
      }
    }
  }
  
  //*************************************************************************/

  template <class T> 
  void storeVecRepI<T>::hdf5write( H5::Group group, DEFAULT_CONVERTER g, 
                                   const entitySet &en) const 
  {
      
    int rank=1;
    hsize_t dimf_store[1];
    std::ostringstream oss;
	
    oss << '{' << en << std::endl ;
    oss << size << std::endl ;
	
    FORALL(en,ii) {
      T* p = base_ptr + ii*size ;
      for(int i=0;i<size;++i,++p)
        oss << *p << " " ;
      oss << std::endl ;
    }ENDFORALL ;
	
    oss << '}' << std::endl ;
	
    std::string memento = oss.str();
    hsize_t m_size      = memento.length();
    dimf_store[0]       = m_size+1;

    try{
      H5::DataSpace dataspace( rank, dimf_store );
      H5::DataSet dataset = group.createDataSet( "storeVec",
                                                 H5::PredType::NATIVE_CHAR, 
                                                 dataspace);
      dataset.write( memento.c_str(), H5::PredType::NATIVE_CHAR );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
      
  };
  
  //************************************************************************/
  
  template <class T>  
  void storeVecRepI<T>:: hdf5write( H5::Group group, IDENTITY_CONVERTER g,
                                    const entitySet &eset ) const
  {
    hsize_t dimension[] = {1};
    int rank = 1;

    HDF5_WriteDomain(group, eset);

    //-----------------------------------------------------------------------
    // write the Vector size
    //-----------------------------------------------------------------------
    dimension[0]=  1;
    try{
      H5::DataSpace sdataspace( rank, dimension );
      H5::DataSet   sdataset = group.createDataSet( "VecSize",
                                                    H5::PredType::NATIVE_INT, 
                                                    sdataspace );
      sdataset.write( &size, H5::PredType::NATIVE_INT );
    }
    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}

    //------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //------------------------------------------------------------------------
    
    entitySet :: const_iterator ci;
    int   offset;

    int arraySize =  size*eset.size();
    
    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------

    T  *data;
    data =  new T[arraySize];
    
    size_t indx= 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        offset = (*ci)*size + ivec;
        data[indx++] = base_ptr[offset];
      }
    }
    
    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------
    typedef hdf5_schema_traits<T> traits_type;
    
    rank = 1;
    dimension[0] =  arraySize;
    
    try {
    
      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = traits_type::get_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", vDatatype, 
                                                      vDataspace);
    
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
  
  //*************************************************************************/
  
  template <class T>  
  void storeVecRepI<T> :: hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, 
                                     const entitySet &eset ) const
  {
    hsize_t   dimension[1];
    int       rank = 1;
      
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

    //-----------------------------------------------------------------------------
    // write the Vector size
    //-----------------------------------------------------------------------------
    int vecsize = size;
    dimension[0]=  1;
    try{
      H5::DataSpace sdataspace( rank, dimension );
      H5::DataSet   sdataset = group.createDataSet( "VecSize",
                                                    H5::PredType::NATIVE_INT, 
                                                    sdataspace );
      sdataset.write( &vecsize, H5::PredType::NATIVE_INT );
    }
    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}
      
    //-----------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-----------------------------------------------------------------------------
    entitySet :: const_iterator ci;
    int   offset, bucketID;
      
    int *vbucket = new int[size*eset.size()];
      
    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;
      
    bucketID = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        offset = (*ci)*size + ivec;
        Memento<T> memento( base_ptr[offset] );
        stateSize           = memento.getSize();
        vbucket[bucketID++] = stateSize;
        arraySize          += stateSize;
        maxStateSize        = max( stateSize, maxStateSize);
      }
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
      for( int ivec = 0; ivec < size; ivec++){
        offset = (*ci)*size + ivec;
        Memento<T> memento( base_ptr[offset] );
        memento.getState( buf, stateSize);
        for( int i = 0; i < stateSize; i++)
          data[indx++] =  buf[i];
      }
    }

    //-----------------------------------------------------------------------------
    // Write size of each container ...
    //-----------------------------------------------------------------------------
    dimension[0]=  size*eset.size();
      
    try {
      H5::DataSpace fDataspace( rank, dimension );
      H5::DataType  fDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   fDataset  = group.createDataSet( "SubContainerSize", 
                                                     fDatatype, fDataspace);
	
      fDataset.write( vbucket, fDatatype );
    }
      
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }
      
    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------
    dimension[0]=  arraySize;

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
      
    //-----------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------
      
    delete [] data;
    delete [] buf;
      
  };
  
  //***************************************************************************/
  
  template <class T> 
  void storeVecRepI<T> :: hdf5read( H5::Group group, DEFAULT_CONVERTER c,
                                    entitySet &en, entitySet &usr )
  {
      
    char ch;

    try{
      H5::DataSet dataset_store = group.openDataSet( "store");
      H5::DataSpace dataspace_store = dataset_store.getSpace();
	
      hsize_t dims_store[1];
      dataspace_store.getSimpleExtentDims( dims_store, NULL);
	
      char* memento = new char[dims_store[0]];
      dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
	
      std::istringstream iss(memento);
      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
	
      entitySet e ;
      iss >> e ;
      iss >> size ;
	
      FORALL(e,ii) {
	      T * p = base_ptr + ii*size ;
	      for(int i=0;i<size;++i,++p)
             iss >> *p;
      } ENDFORALL ;
	
      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
	
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
  };
  
  //**************************************************************************/
  
  template <class T> 
  void storeVecRepI<T>::hdf5read(H5::Group group, IDENTITY_CONVERTER convert, 
                                 entitySet &eset, entitySet &user_eset)
  {

    hsize_t dimension[1];
    int  indx=0, rank=1;

    entitySet::const_iterator ci;

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<unsigned> offset;
    offset.allocate(eset);

    int arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += size;
    }

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    T   *data;

    dimension[0] = size*eset.size();
    H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
    H5::DataSpace vDataspace(rank, dimension);

    typedef hdf5_schema_traits<T> traits_type;
    H5::DataType vDatatype = traits_type::get_type();
    H5::DataSet  vDataset   = group.openDataSet( "VariableData");

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    int voffset;
    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  size;

      data = new T[count[0]];

      foffset[0] = offset[it[k].first];

      mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
      vDataset.read( data, vDatatype, mDataspace, vDataspace);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int ivec = 0; ivec < size; ivec++){
          voffset           = i*size+ivec;
          base_ptr[voffset] = data[indx++];
        }
      }
      delete[] data;
    }
    
  };
  
  //**************************************************************************/

  template <class T> 
  void storeVecRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
  {
      
    hsize_t dimension[1];
    int indx = 0, arraySize;
    int    rank = 1, vecsize;
    
    entitySet::const_iterator ci;
    
    typedef hdf5_schema_converter_traits<T> converter_traits; 

    //--------------------------------------------------------------------------
    // Size of each sub-container ....
    //--------------------------------------------------------------------------
    
    H5::DataType  sDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   sDataset   = group.openDataSet( "SubContainerSize");
    H5::DataSpace sDataspace = sDataset.getSpace();
    
    sDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

    sDataset.read( ibuf, H5::PredType::NATIVE_INT );

    int maxBucketSize = *std::max_element( ibuf, ibuf + (int)dimension[0] );

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store< unsigned int >   offset;
    dmultiStore<int>  subcontainer;
    offset.allocate( eset );

    arraySize = 0;
    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      for( int i = 0; i < size; i++)  {
        vecsize    =  ibuf[indx++];
        arraySize  += vecsize;
        subcontainer[*ci].push_back( vecsize );
      }
    }
    delete [] ibuf;

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
   
    int num_intervals = user_eset.num_intervals();

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

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
      for( int i = it[k].first; i <= it[k].second; i++){
        for( int j = 0; j < size; j++)
          count[0] +=  subcontainer[i][j];
      }

      data = new typename converter_traits::memento_type[count[0]];

      foffset[0] = offset[it[k].first];

      mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
      vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
      vDataset.read( data, vDatatype, mDataspace, vDataspace);

      indx = 0;
      int bucsize;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int j = 0; j < size; j++) {
          Memento<T> memento( base_ptr[i*size+j] );
          bucsize = subcontainer[i][j];
          for( int m = 0; m < bucsize; m++) 
            buf[m] = data[indx++];
          base_ptr[i*size+j] = memento.setState( buf, bucsize );
        }
      }
      delete[] data;
    }

    delete[] buf;

  }; 
  
  



}

#endif

#ifndef STOREVEC_H
#define STOREVEC_H 

#include <istream>
#include <ostream>
#include <hdf5_readwrite.h>
#include <DMultiStore.h>
#include <store.h>
#include <vector>

namespace Loci {
  extern double total_memory_usage ;
  extern std::ofstream debugout ;
  //*******************************************************************/
  
  template <class T> struct Scalar {
    T val ;
    Scalar(T v) : val(v) { }
  } ;

  template <class T> Scalar<T> mk_Scalar(T v) { return Scalar<T>(v) ;} 

  //*******************************************************************/
  
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

  //******************************************************************/

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

  //******************************************************************/

  template<class T> class storeVecRepI : public storeRep {
    entitySet    store_domain ;
    T           *alloc_pointer, *base_ptr ;
    int          size ;
    lmutex       mutex ;
    bool istat ;

    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void hdf5read( hid_t group_id, IDENTITY_CONVERTER     c, entitySet &en, entitySet &usr);
    void hdf5write( hid_t group_id, IDENTITY_CONVERTER g,     const entitySet &en) const;
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size,
                    const sequence &seq) ;

    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);
    void hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, const entitySet &en) const;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size,
                    const sequence &seq) ;

  public:
    storeVecRepI() 
    { alloc_pointer= 0 ; base_ptr = 0 ; size=0; istat = 1 ; }

    storeVecRepI(const entitySet &p) 
    { size = 0; alloc_pointer=0 ; allocate(p) ; istat = 1 ; }
    
    virtual ~storeVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e ) ;
    virtual void pack(void * ptr, int &loc, int &size, const entitySet &e ) ;
    virtual void unpack(void * ptr, int &loc, int &size, const sequence &seq) ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( hid_t group, entitySet &en) ;
    virtual void writehdf5( hid_t group, entitySet& en) const ;
    virtual void set_elem_size(int sz) ;
    bool is_static() {return istat; }
    T *get_base_ptr() const { return base_ptr ; }
    int get_size() const { return size ; }
  } ;

  //******************************************************************/

  template<class T> 
  std::ostream &storeVecRepI<T>::Print(std::ostream &s) const
  {

    s << '{' << domain() << std::endl ;
    s << size << std::endl ;
    
    FORALL(domain(),ii) {
      T * p = base_ptr + ii*size ;
      Loci::streamoutput(p,size,s) ;
    }ENDFORALL ;
    s << '}' << std::endl ;

    return s ;
  }

  //******************************************************************/

  template<class T> 
  std::istream &storeVecRepI<T>::Input(std::istream &s)
  {

    //------------------------------------------------------------------
    // Objective : Read the storeVec from the input stream.
    //------------------------------------------------------------------
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
      for(int i=0;i<size;++i)
        p[i] = T() ;
      
      Loci::streaminput(p,size,s) ;
    } ENDFORALL ;
    
    // Look for the closing brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
	  
    return s ;
  }

  //*****************************************************************/
  template<class T> 
  void storeVecRepI<T>::allocate(const entitySet &ptn) {
  
    /*
    entitySet common = store_domain & ptn ;
    T *tmp_base_ptr, *tmp_alloc_pointer ;
    tmp_alloc_pointer = 0 ;
    tmp_base_ptr = 0 ;
    if(store_domain != EMPTY) {
      if(size > 0) {
	if(common != EMPTY) {
	  int top = ptn.Min() ; 
	  int sza = (ptn.Max() - top + 1) * size ;
	  tmp_alloc_pointer = new T[sza] ; 
	  tmp_base_ptr = tmp_alloc_pointer - top*size ;
	  FORALL(common,i) { 
	    T *p1 = tmp_base_ptr + i*size ;
	    T* p2 = base_ptr + i*size ;
	    for(int j = 0; j < size; ++j)
	      p1[j] = p2[j] ;
	  } ENDFORALL ; 
	}
      }
    }
    else {
      if(size != 0) {
	fatal(size < 1) ;
	if(ptn != EMPTY) {
	  int top = ptn.Min() ; 
	  int sza = (ptn.Max()-top+1)*size ;
	  tmp_alloc_pointer = new T[sza] ;
	  tmp_base_ptr = alloc_pointer - top*size ;
	}
      }
    }
    
    //------------------------------------------------------------------
    // Allocation reclaims all previously hold memeory 
    //------------------------------------------------------------------
    //int p = 0 ;
    
    if(alloc_pointer) 
      delete[] alloc_pointer ;
    
    //-----------------------------------------------------------------
    // Get the minimum and maximum entity ID from the entitySet and
    //allocate 
    // memory of the size = ( max-min+1). Notice that, if the entityset 
    // contains the entities with ID quite sparse, it will create lots of 
    // unused block of memory. 
    //------------------------------------------------------------------
    alloc_pointer = tmp_alloc_pointer ;
    base_ptr = tmp_base_ptr ;
    */
    /*
      if(size != 0) {
      fatal(size < 1) ;
      if(common == EMPTY) {
      if(ptn != EMPTY) {
      int top = ptn.Min() ; 
      int sza = (ptn.Max() - top + 1) * size ;
      alloc_pointer = new T[sza] ;
      base_ptr = alloc_pointer - top*size ;
      } 
      }
      }
    */
    /*
      if(!p) {
      debugout << " ******************************************" << endl ;
      debugout << "Size = " << size << "  Domain = " << ptn << endl ;
      total_memory_usage += double(sza*sizeof(T)) ;
      debugout << " Allocated  " << double(sza*sizeof(T)) << " bytes :   Total Memory Usage is now " << total_memory_usage << endl ;
      debugout << " ******************************************" << endl ;
      }
    */
    if(alloc_pointer) delete[] alloc_pointer ;
    
    alloc_pointer = 0 ;
    base_ptr      = 0 ;

    if(size != 0) {
      fatal(size < 1) ;
      if(ptn != EMPTY) {
        int top = ptn.Min() ; int sza = (ptn.Max()-top+1)*size ;
        alloc_pointer = new T[sza] ;
        base_ptr = alloc_pointer - top*size ;
      }
    }
    
    store_domain = ptn ;
    dispatch_notify() ;
  }
  
  
  //*******************************************************************/

  template<class T> 
  storeVecRepI<T>::~storeVecRepI<T>() 
  {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  //*******************************************************************/

  template<class T>
  storeRep *storeVecRepI<T>::new_store(const entitySet &p) const 
  {
    storeRep *sp = new storeVecRepI<T>(p)  ;
    sp->set_elem_size(size) ;
    return sp ;
  }

  //*******************************************************************/

  template<class T> 
  store_type storeVecRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //******************************************************************/

  template<class T> 
  entitySet storeVecRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //*****************************************************************/

  template<class T> 
  void storeVecRepI<T>::set_elem_size(int sz) 
  {
    mutex.lock() ;

    //----------------------------------------------------------------
    // Change the size of vector held. It will reclaim the memory used before
    // his call. and allocate new one.
    //----------------------------------------------------------------
    
    if(size != sz) {
      if(size > 0) {
        if(sz < 1000*size)
	  size = sz ;
      }
      else if (size == 0)
	if(sz < 1000)
	  size = sz ;
      allocate(store_domain) ;
    }
    mutex.unlock() ;
  }
  
  //*******************************************************************/
      
  template<class T> class storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    storeVec() {setRep(new storeType) ;}
    storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    storeVec(storeRepP rp) { setRep(rp) ;}

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

  //*******************************************************************/

  template<class T> 
  storeVec<T>::~storeVec<T>() { }
  
  //******************************************************************/

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

  //******************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const storeVec<T> &t)
  { return t.Print(s) ; }

  //*******************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, storeVec<T> &t)
  { return t.Input(s) ; }

  //*******************************************************************/

  template<class T> class const_storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_storeVec() { setRep(new storeType) ; }
    const_storeVec(const_storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeRepP rp) { setRep(rp) ; }
    
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

  //*******************************************************************/

  template<class T> 
  const_storeVec<T>::~const_storeVec<T>() { }

  //*******************************************************************/

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

  //******************************************************************/

  template<class T> store_instance::instance_type
  const_storeVec<T>::access() const
  { return READ_ONLY ; }

  //******************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_storeVec<T> &t)
  { return t.Print(s) ; }

  //************************************************************************/

  template <class T> 
  storeRepP storeVecRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    storeVec<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  //*************************************************************************/

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
  void storeVecRepI<T>::gather(const dMap &m, storeRepP &st,
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
  void storeVecRepI<T>::scatter(const dMap &m, storeRepP &st,
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
  inline int storeVecRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    int size, M ;
    M = get_size() ;
    size = (sizeof(T) * eset.size() * M ) + sizeof(int) ;
    return(size) ;

  }

  //**************************************************************************/

  template <class T>
  int storeVecRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {

    int     arraySize =0, numContainers = 0;
    entitySet::const_iterator ci;

    typedef data_schema_traits<T> converter_traits;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename converter_traits::Converter_Type cvtr(base_ptr[(*ci)*size+ivec]) ;
        arraySize += cvtr.getSize() ;
      }
    }

    numContainers =  size*eset.size();

    return(arraySize*sizeof(typename converter_traits::Converter_Base_Type) +
           (numContainers+1)*sizeof(int));
  }

  //**************************************************************************/

  template <class T>
  int storeVecRepI<T>::pack_size( const entitySet &eset)
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }

  //*************************************************************************/

  template <class T>
  inline void storeVecRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                  int outcount, const entitySet &eset )
  {

    const int M = get_size() ;
    for(int i = 0; i < eset.num_intervals(); ++i) {
      Loci::int_type indx1 = eset[i].first ;
      Loci::int_type stop  = eset[i].second ;
      T *p = base_ptr + M * indx1 ;
      const int t = (stop - indx1 + 1) * M ;
      MPI_Pack( p, t*sizeof(T), MPI_BYTE, outbuf, outcount, &position, 
                MPI_COMM_WORLD) ;
    }
  }

  //**************************************************************************/

  template <class T> 
  inline void storeVecRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                  int &position, int outcount, 
                                  const entitySet &eset ) 
  {
    entitySet::const_iterator ci;

    warn(ecommon-store_domain != EMPTY) ;

    //-------------------------------------------------------------------------
    // Get the maximum size of container 
    //-------------------------------------------------------------------------
    int stateSize, maxStateSize=0;
    
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( base_ptr[(*ci)*size+ivec] );
        stateSize = cvtr.getSize();
        maxStateSize = max( maxStateSize, stateSize);
      }
    }

    typedef data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);
    std::vector<dtype> inbuf(maxStateSize);

    int incount;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( base_ptr[(*ci)*size+ivec] );
        cvtr.getState( &inbuf[0], stateSize);

        MPI_Pack(&stateSize,1, MPI_INT, outbuf, outcount,&position,
                 MPI_COMM_WORLD);

        incount =  stateSize*typesize;
        MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
                 MPI_COMM_WORLD) ;
      }
    }
  }

  //**************************************************************************/

  template <class T>
  void storeVecRepI<T>::pack(void *outbuf, int &position, int &outcount,
                             const entitySet &eset )
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    int M = get_size() ;
    MPI_Pack( &M, 1, MPI_INT, outbuf, outcount, &position, 
	      MPI_COMM_WORLD) ;
    
    packdata( traits_type, outbuf, position, outcount, eset);
  }

  //**************************************************************************/
  template <class T> 
  inline void storeVecRepI<T>::unpackdata( IDENTITY_CONVERTER c,
                                           void *inbuf,
                                           int &position,
                                           int &insize,
                                           const sequence &seq)
  {
    const int M = get_size() ;
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type indx1 = seq[i].first ;
        const Loci::int_type stop =  seq[i].second ;
        for(Loci::int_type indx = indx1; indx != stop-1; --indx) {
          T *p = base_ptr + M * indx ;
          MPI_Unpack( inbuf, insize, &position, p, sizeof(T) * M, 
                      MPI_BYTE, MPI_COMM_WORLD) ;
        }
      } else {
        Loci::int_type indx1 = seq[i].first ;
        Loci::int_type stop = seq[i].second ;
        T *p = base_ptr + M * indx1 ;
        const int t = (stop - indx1 + 1) * M ;
        MPI_Unpack( inbuf, insize, &position, p, t * sizeof(T), 
                    MPI_BYTE, MPI_COMM_WORLD) ;
      }
    }

  }

  //***********************************************************************/
  template <class T> 
  void storeVecRepI<T>::unpackdata( USER_DEFINED_CONVERTER c,
                                    void *inbuf, 
                                    int &position,
                                    int &insize,
                                    const sequence &seq)
  {

    sequence :: const_iterator ci;
    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    std::vector<dtype> outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++) {
        MPI_Unpack( inbuf, insize, &position, &stateSize, 1, 
                    MPI_INT, MPI_COMM_WORLD) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);

        outcount = stateSize*typesize;
        MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount, 
                    MPI_BYTE, MPI_COMM_WORLD) ;

        typename schema_traits::Converter_Type  cvtr( base_ptr[(*ci)*size+ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }
  }

  //**************************************************************************/
  template <class T> 
  void storeVecRepI<T>::unpack(void *inbuf, int &position, int &insize, const sequence &seq)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

#ifdef DEBUG
    entitySet ecommon, ediff,eset(seq);

    ediff = eset - domain();
    if( ediff.size() > 0) { 
      std::cout << "Error:Entities not part of domain " << ediff <<endl;
      abort();
    }
#endif
    
    int init_size = get_size() ;
    int M ;
    MPI_Unpack(inbuf, insize, &position, &M, 1, MPI_INT, MPI_COMM_WORLD) ;
    warn(M<=0) ;
    if(M > init_size) {
      set_elem_size(M) ;
    }
    unpackdata( traits_type, inbuf, position, insize, seq);
  }
  
  //*******************************************************************/

  template<class T> 
  void storeVecRepI<T>::writehdf5( hid_t group_id, entitySet &usr_eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet eset(usr_eset&domain());

    if( eset.size() < 1) return;
    int vsize = get_size();

    Loci::HDF5_WriteVecSize(group_id, vsize);
    Loci::HDF5_WriteDomain(group_id,  eset);

    hdf5write(group_id, traits_output_type, eset);

  }

  //************************************************************************/
  
  template <class T>  
  void storeVecRepI<T>:: hdf5write( hid_t group_id, IDENTITY_CONVERTER g,
                                    const entitySet &eset ) const
  {

    hsize_t dimension;
    int rank = 1;
    int vsize = get_size();
    
    //------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //------------------------------------------------------------------------
    
    entitySet :: const_iterator ci;

    int arraySize =  vsize*eset.size();

    //------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //------------------------------------------------------------------------

    T  *data;
    data =  new T[arraySize];
    
    size_t indx= 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < vsize; ivec++)
        data[indx++] = base_ptr[(*ci)*vsize+ivec];
    }

    //------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //------------------------------------------------------------------------
    typedef data_schema_traits<T> traits_type;
    DatatypeP dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    rank      = 1;
    dimension = arraySize;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, 
                                 vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );
    delete [] data;

  }
  
  //*************************************************************************/
  
  template <class T>  
  void storeVecRepI<T> :: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, 
                                     const entitySet &eset ) const
  {
    hsize_t   dimension;
    hid_t     vDataspace, vDataset, vDatatype;

    int rank  = 1;
    int vsize = get_size();
      
    //-----------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-----------------------------------------------------------------------

    entitySet :: const_iterator ci;
    int   bucketID;
      
    std::vector<int> vbucket(size*eset.size());
      
    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    typedef data_schema_traits<T> schema_traits ;
    
    bucketID = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)*size+ivec] );
        stateSize           = cvtr.getSize();
        vbucket[bucketID++] = stateSize;
        arraySize          += stateSize;
        maxStateSize        = max( stateSize, maxStateSize);
      }
    }


    //-------------------------------------------------------------------------
    // Write size of each container ...
    //-------------------------------------------------------------------------
    dimension =  size*eset.size();

    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5Tcopy( H5T_NATIVE_INT);
    vDataset   = H5Dcreate(group_id, "SubContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vbucket[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    //-------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------------
      
    typedef typename schema_traits::Converter_Base_Type dtype;

    dtype *data ;
      
    data =  new dtype[arraySize];

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename schema_traits::Converter_Type cvtr( base_ptr[(*ci)*size+ivec] );
        cvtr.getState( data+indx, stateSize);
        indx +=stateSize ;
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
      
  }
  
  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    int     vec_size;
    entitySet   eset, ecommon;

    Loci::HDF5_ReadVecSize(group_id, &vec_size);

    set_elem_size(vec_size) ;

    Loci::HDF5_ReadDomain(group_id, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );
    hdf5read( group_id, traits_type, eset, ecommon);

  }
  //**************************************************************************/

  template <class T> 
  void storeVecRepI<T>::hdf5read(hid_t group_id, IDENTITY_CONVERTER convert, 
                                 entitySet &eset, entitySet &user_eset)
  {

    hsize_t dimension;
    int  indx=0, rank=1;

    entitySet::const_iterator ci;

    //------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //------------------------------------------------------------------------
    store<unsigned> offset;
    offset.allocate(eset);
    int vsize = get_size();

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci){
      offset[*ci] = indx;
      indx       += vsize;
    }
    int arraySize = vsize*eset.size();

    //------------------------------------------------------------------------
    // Read the data now ....
    //------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    dimension = arraySize;

    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    std::vector<T>  data;
    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  size;

      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, &data[0]);
 
      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int ivec = 0; ivec < vsize; ivec++)
          base_ptr[i*vsize+ivec] = data[indx++];
      }
    }

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Sclose( mDataspace);
    H5Tclose( vDatatype );
  }
  
  //************************************************************************/

  template <class T> 
  void storeVecRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
  {
    hsize_t  dimension;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;
    
    int indx = 0, arraySize;
    int    rank = 1, vecsize;

    entitySet::const_iterator ci;

    //---------------------------------------------------------------------

    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"SubContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    std::vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    int maxBucketSize = *std::max_element( ibuf.begin(), ibuf.end() );

    //---------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------
    
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


    //---------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    typedef data_schema_traits<T> converter_traits;
    typename converter_traits::Converter_Base_Type  *data, *buf, dtype;
    DatatypeP base_type =
      data_schema_traits<typename converter_traits::Converter_Base_Type>::get_type() ;
    
    vDatatype = base_type->get_hdf5_type();

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


    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++){
        for( int j = 0; j < size; j++)
          count[0] +=  subcontainer[i][j];
      }

      data = new typename data_schema_traits<T>::Converter_Base_Type[count[0]];

      foffset[0] = offset[it[k].first];
      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, data);

      indx = 0;
      int bucsize;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int j = 0; j < size; j++) {
          typename data_schema_traits<T>::Converter_Type cvtr( base_ptr[i*size+j] );
          bucsize = subcontainer[i][j];
          cvtr.setState( data+indx, bucsize );
          indx += bucsize ;
        }
      }
      delete[] data;
    }
  }
  //******************************************************************/


}

#endif

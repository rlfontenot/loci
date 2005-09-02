#ifndef STOREVEC_H
#define STOREVEC_H 

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <hdf5_readwrite.h>
#include <DMultiStore.h>
#include <store.h>
#include <vector>
#include <dist_internal.h>

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
      fatal(idx >= size || idx < 0) ;
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
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    void operator=(const const_Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
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
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
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
    
    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
    void packdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int &size,
                    const sequence &seq) ;

    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en) const;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size,
                    const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
  public:
    storeVecRepI() 
      { alloc_pointer= 0 ; base_ptr = 0 ; size=0 ; }
    
    storeVecRepI(const entitySet &p) 
      { size = 0; alloc_pointer=0 ; allocate(p) ;}
    
    virtual ~storeVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &eset, const int* p) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
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
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    virtual void set_elem_size(int sz) ;
    T *get_base_ptr() const { return base_ptr ; }
    int get_size() const { return size ; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
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
  
  
  template<class T>
    void storeVecRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    base_ptr -= offset*size ;
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
  
  template<class T>
  storeRep *storeVecRepI<T>::new_store(const entitySet &p, const int* cnt) const 
    {
      storeRep* sp = 0;
      cerr << " This method should not be called for a storeVec " << endl ;
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
    storeVec(const storeVec<T> &var) {setRep(var.Rep()) ;}
    storeVec<T> & operator=(const storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

  public:
    typedef Vect<T> containerType ;
    storeVec() {setRep(new storeType) ;}
    storeVec(storeRepP rp) { setRep(rp) ;}

    virtual ~storeVec() ;
    virtual void notification() ;

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
    const_storeVec(const const_storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(const storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec<T> & operator=(const storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_storeVec<T> & operator=(const const_storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

  public:
    typedef const_Vect<T> containerType ;
    const_storeVec() { setRep(new storeType) ; }
    const_storeVec(storeRepP rp) { setRep(rp) ; }
    
    virtual ~const_storeVec() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
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

  template<class T>
  storeRepP storeVecRepI<T>::freeze() {
    return getRep() ;
  }

  // This method is also incorrect for reasons
  // see the comments before storeRepI<T>::thaw()
  // in store.h
  template<class T>
  storeRepP storeVecRepI<T>::thaw() {
    return getRep() ;
//     dstoreVec<T> ds ;
//     ds.setVecSize(size) ;
//     for(entitySet::const_iterator ei=store_domain.begin();
//         ei!=store_domain.end();++ei) {
//       T* b = base_ptr+( (*ei)*size) ;
//       T* e = b + size ; 
//       ds[*ei] = std::vector<T>(b,e) ;
//     }
//     return ds.Rep() ;
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
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
	
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

    fatal((eset - domain()) != EMPTY);

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
    warn(M == 0);
    fatal(M<0) ;
    if(M > init_size) {
      set_elem_size(M) ;
    }
    unpackdata( traits_type, inbuf, position, insize, seq);
  }
  
  template<class T> 
    frame_info storeVecRepI<T>::read_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return read_frame_info(group_id, schema_converter()) ;
  }
  template<class T>  
    frame_info storeVecRepI<T>::read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    int is_stat = 0 ;
    int sz = 0 ;
    if(Loci::MPI_rank == 0) {
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dopen(group_id, "is_stat") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dopen(group_id, "vec_size") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
      H5Dclose(dataset) ;
    }
    int dim[2] ;
    dim[0] = is_stat ;
    dim[1] = sz ;
    MPI_Bcast(&dim, 2, MPI_INT, 0, MPI_COMM_WORLD) ;
    return frame_info(dim[0], dim[1]);
  }
  template<class T> 
    frame_info storeVecRepI<T>::read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    hid_t datatype = H5T_NATIVE_INT ;
    hid_t dataset ;
    int is_stat = 0 ;
    int sz = 0 ;
    frame_info fi ;
    if(Loci::MPI_rank == 0) {
      dataset = H5Dopen(group_id, "is_stat") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dopen(group_id, "vec_size") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
      H5Dclose(dataset) ;
    }
    int dim[2] ;
    dim[0] = is_stat ;
    dim[1] = sz ;
    MPI_Bcast(&dim, 2, MPI_INT, 0, MPI_COMM_WORLD) ;
    fi.is_stat = dim[0] ;
    fi.size = dim[1] ;
    std::vector<int> vint ;
    int dom_size = domain().size() * fi.size ;
    read_vector_int(group_id, "second_level", vint,dom_size) ;
    fi.second_level = vint ; 
    return fi ;
  }
  
  template<class T> 
    frame_info storeVecRepI<T>::write_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return write_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info storeVecRepI<T>::write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = get_size() ;
    if(Loci::MPI_rank == 0 ) {
      hsize_t dimension = 1 ;
      int rank = 1 ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
    return fi ;
  }
  
  template<class T> 
    frame_info storeVecRepI<T>::write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = get_size() ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) 
      for(int ivec = 0; ivec < fi.size; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)*fi.size+ivec] );
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
    hsize_t dimension = 0 ;
    hid_t dataspace ;
    hid_t datatype = H5T_NATIVE_INT ;
    int rank = 1 ;
    if(MPI_rank == 0) {
      dimension = 1 ;
      dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ; 
    }
    write_vector_int(group_id, "second_level", fi.second_level) ;
    return fi ;
  }
  
  template<class T> 
    DatatypeP storeVecRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
    DatatypeP storeVecRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
    DatatypeP storeVecRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //*******************************************************************/

  template<class T> 
  void storeVecRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset) ;
  }
  
  //************************************************************************/
  
  template <class T>  
    void storeVecRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
    {
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	T* tmp_array = new T[dimension] ;
	size_t tmp = 0 ;
	int qs = get_size() ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
	  for(int ivec = 0; ivec < qs; ivec++){
	    tmp_array[tmp++] = base_ptr[(*si)*qs+ivec] ;
	  }
	}
	H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  
  //*************************************************************************/
  
  template <class T>  
    void storeVecRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
    { 
      typedef data_schema_traits<T> schema_traits ;
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	typedef typename schema_traits::Converter_Base_Type dtype;
	dtype* tmp_array = new dtype[dimension] ;
	size_t tmp = 0 ;
	int stateSize = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
	  for(int ivec = 0; ivec < size; ivec++){
	    typename schema_traits::Converter_Type cvtr(base_ptr[(*si)*size+ivec]);
	    cvtr.getState(tmp_array+tmp, stateSize) ;
	    tmp +=stateSize ;
	  }
	}
	H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  
  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
    {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
    }
  //**************************************************************************/
  
  template <class T> 
    void storeVecRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER convert, frame_info &fi, entitySet &eset)
    {
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	T* tmp_array = new T[dimension] ;
	size_t tmp = 0 ;
	int qs = fi.size ;
	hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_array) ;
	if(err < 0) {
	  cerr << "H5Dread() failed" << endl ;
	}

	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
	  for(int ivec = 0; ivec < qs; ivec++) {
	    base_ptr[(*si)*qs+ivec] = tmp_array[tmp++] ;
	  }
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  
  //************************************************************************/

  template <class T> 
    void storeVecRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
    {
      typedef data_schema_traits<T> schema_traits ;
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	std::vector<int> vint = fi.second_level ;
	typedef typename schema_traits::Converter_Base_Type dtype;
	dtype* tmp_array = new dtype[dimension] ;    
	hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_array) ;
	if(err < 0) {
	  cerr << "H5Dread() failed" << endl ;
	}
	int qs = fi.size ;
	size_t tmp = 0 ;
	int bucsize ;
	size_t indx = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) 
	  for(int ivec = 0; ivec < qs; ivec++) {
	    typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[(*si)*qs+ivec]);
	    bucsize = vint[indx++] ;
	    cvtr.setState(tmp_array+tmp, bucsize) ;
	    tmp += bucsize ;
	  }
      
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  //******************************************************************/
  

}
 
#endif

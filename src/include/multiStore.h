#ifndef MULTISTORE_H
#define MULTISTORE_H

#include <istream>
#include <ostream>

extern "C" {
#include <hdf5.h>
}

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>

#include <Map.h>
#include <multiMap.h>
#include <DMultiStore.h>

namespace Loci {

  template<class T> class multiStoreRepI : public storeRep {
    entitySet store_domain ;
    T **index ;
    T *alloc_pointer ;
    T **base_ptr ;
    int size ;
    lmutex mutex ;
    bool istat ;

#ifdef ALLOW_DEFAULT_CONVERTER
    void  hdf5read( hid_t group_id, DEFAULT_CONVERTER c,      entitySet &en,
                    entitySet &usr);
    void  hdf5write( hid_t group_id, DEFAULT_CONVERTER c,      const entitySet &en) const;
    int   get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    void  packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                     const sequence &seq );
    void  StringVal( const int &entity, const int &ivec, std::string &memento);

#endif
    void  hdf5read( hid_t group_id, IDENTITY_CONVERTER c,     entitySet &en,
                    entitySet &usr);
    void  hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, entitySet &en,
                    entitySet &usr);
    void  hdf5write( hid_t group_id, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( hid_t group_id, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq );
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq);

  public:
    multiStoreRepI()
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; istat = 1; }

    multiStoreRepI(const entitySet &p)
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; store_domain=p; istat = 1;}

    multiStoreRepI(const store<int> &sizes) {
      index = 0 ; alloc_pointer=0 ; base_ptr = 0; allocate(sizes) ;istat = 1; }

    void allocate(const store<int> &sizes) ;
    void multialloc(const store<int> &count, T ***index, T **alloc_pointer, T ***base_ptr) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~multiStoreRepI() ;
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
    virtual void readhdf5( hid_t group_id, entitySet &en) ;
    virtual void writehdf5( hid_t group_id, entitySet& en) const ;

    bool is_static() { return istat ; } 
    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
  } ;
  
  //*************************************************************************/
  
  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~multiStore() ;
    virtual void notification() ;
    
    multiStore<T> & operator=(multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    
    multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void setSizes(const const_multiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    const entitySet domain() const { return Rep()->domain() ; }

    Vect<T> elem(int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }

    Vect<T> operator[](int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }
    
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int index) const { return end(index)-begin(index); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  //*************************************************************************/

  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const multiStore<T> &m)
  { return m.Print(s) ; }

  //************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, multiStore<T> &m)
  { return m.Input(s) ; }
 
  //************************************************************************/

  template<class T> 
  multiStore<T>::~multiStore() {}

  //************************************************************************/
  
  template<class T> 
  void multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //************************************************************************/
  
  template<class T> class const_multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~const_multiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
    const_multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(const const_multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    const entitySet domain() const { return Rep()->domain() ; }

    containerType elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    containerType operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }

    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int index) const { return end(index)-begin(index); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //*************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_multiStore<T>::access() const
  { return READ_ONLY ; }

  //*************************************************************************/
  
  template<class T> 
  const_multiStore<T>::~const_multiStore() {}

  //*************************************************************************/
  
  template<class T> 
  void const_multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::allocate(const store<int> &sizes) 
  {
    //-------------------------------------------------------------------------
    // Objective: Allocate memeory for multiStore data. This call reclaims 
    // all previously held memory
    //-------------------------------------------------------------------------
    // Assign new entitySet ...
    entitySet ptn = sizes.domain() ;
    store_domain  = ptn ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    if(index) delete[] index ;

    index = 0 ;
    int sz = 0 ;
    if(ptn != EMPTY) {
      int top  = ptn.Min() ;
      int len  = ptn.Max() - top + 2 ;
      index    = new T *[len] ;
      base_ptr = index - top ;

      FORALL(ptn,i) {
        sz += sizes[i] ;
      } ENDFORALL ;

      alloc_pointer = new T[sz+1] ;
      sz = 0 ;
      for(int ivl=0;ivl< ptn.num_intervals(); ++ivl) {
        int i       = ptn[ivl].first ;
        base_ptr[i] = alloc_pointer + sz ;
        while(i<=ptn[ivl].second) {
          sz += sizes[i] ;
          ++i ;
          base_ptr[i] = alloc_pointer + sz ;
        }
      }

    }
    dispatch_notify();
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::multialloc(const store<int> &count, T ***index, 
                                     T **alloc_pointer, T ***base_ptr ) {
    entitySet ptn = count.domain() ;
    int top = ptn.Min() ;
    int len = ptn.Max() - top + 2 ;
    T **new_index = new T *[len] ;
    T **new_base_ptr = new_index - top ;
    int sz = 0 ;
    
    FORALL(ptn, i) {
      sz += count[i] ;
    } ENDFORALL ;
    
    T *new_alloc_pointer = new T[sz + 1] ;
    sz = 0 ;
    
    for(int ivl = 0; ivl < ptn.num_intervals(); ++ivl) {
      int i = ptn[ivl].first ;
      new_base_ptr[i] = new_alloc_pointer + sz ;
      while(i <= ptn[ivl].second) {
	sz += count[i] ;
	++i ;
	new_base_ptr[i] = new_alloc_pointer + sz ;
      }
    }
    
    *index = new_index ;
    *alloc_pointer = new_alloc_pointer ;
    *base_ptr = new_base_ptr ;
    
  }

  //*************************************************************************/
   
  template<class T> 
  void multiStoreRepI<T>::setSizes(const const_multiMap &mm)
  {
    //------------------------------------------------------------------------
    // Objective : Set the degree of each entity specified by the map..
    //------------------------------------------------------------------------

    mutex.lock() ;

    if(alloc_pointer != 0 && base_ptr[store_domain.Min()] == base_ptr[store_domain.Max()]) {
      delete[] index ;
      delete[] alloc_pointer ;
      index = 0 ;
      alloc_pointer = 0 ;
    }

    if(alloc_pointer != 0) {
      entitySet map_set = mm.domain() & store_domain ;
      entitySet problem ;
      FORALL(map_set,i) {
        if((end(i)-begin(i))<(mm.end(i)-mm.begin(i)))
          problem += i ;
      } ENDFORALL ;

      if(problem != EMPTY) {
        std::cerr << "reallocation of multiStore required for entities"
                  << problem << endl
                  << "Currently this reallocation isn't implemented."
                  << endl ;
      }
    } else {
      store<int> sizes ;
      sizes.allocate(store_domain) ;
      FORALL(store_domain,i) {
        sizes[i] = 0 ;
      } ENDFORALL ;
      entitySet map_set = mm.domain() & store_domain ;
      FORALL(map_set,i) {
        sizes[i] = (mm.end(i) - mm.begin(i)) ;
      } ENDFORALL ;
      allocate(sizes) ;
    }
    mutex.unlock() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::allocate(const entitySet &ptn) 
  {
    //------------------------------------------------------------------------
    // Objective : allocate memory specified by the entitySet. Allocation
    // doesn't resize the memory, therefore reclaims previously held memory.
    //------------------------------------------------------------------------
    
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
    
    alloc_pointer = 0 ;
    index         = 0 ;
    base_ptr      = 0 ;

    store_domain  = ptn ;

    //-------------------------------------------------------------------------
    // Initialize degree of each entity to zero 
    //-------------------------------------------------------------------------

    store<int> count ;
    count.allocate(ptn) ;
    
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    
    allocate(count) ;

    //-------------------------------------------------------------------------
    // Notify all observers ...
    //-------------------------------------------------------------------------

    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  multiStoreRepI<T>::~multiStoreRepI() 
  {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
  }

  //*************************************************************************/

  template<class T> 
  storeRep *multiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    store<int> count ;
    count.allocate(p) ;
    entitySet ent = p - domain() ;
    
    for(entitySet::const_iterator ei = p.begin(); ei != p.end(); ++ei)
      count[*ei] = base_ptr[*ei+1] - base_ptr[*ei] ;
    
    for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei)
      count[*ei] = 0 ;
    
    return new multiStoreRepI<T>(count) ;
  }

  //*************************************************************************/

  template<class T> 
  storeRepP multiStoreRepI<T>::remap(const Map &m) const 
  {

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    multiStore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_multiStore<T> s(st) ;
    fatal(alloc_pointer == 0) ;
    fatal((context - domain()) != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    store<int> count ;
    count.allocate(domain()) ;
    
    FORALL(domain() - context, i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context, i) {
      count[i] = s.end(i) - s.begin(i) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::gather(const Map &m, storeRepP &st,
                                 const entitySet &context) 
  {
    store<int> count ;
    const_multiStore<T> s(st) ;
    count.allocate(domain()) ;
    
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[i] = s.end(m[i])-s.begin(m[i]) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;

    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    
    FORALL(domain()-context,i) {
      for(int j = 0; j < count[i]; ++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j = 0; j < count[i]; ++j)
        new_base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::scatter(const Map &m, storeRepP &st,
                                  const entitySet &context) 
  {
    
    store<int> count;
    
    const_multiStore<T> s(st) ;
    count.allocate(domain());
    
    FORALL(domain()-m.image(context),i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[m[i]] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    
    FORALL(domain() - m.image(context),i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[m[i]];++j) {
        new_base_ptr[m[i]][j] = s[i][j] ;
      }
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    
    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  store_type multiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //*************************************************************************/
  
  template<class T> 
  entitySet multiStoreRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //*************************************************************************/
  
  template<class T> 
  std::ostream &multiStoreRepI<T>::Print(std::ostream &s) const 
  {
    //-------------------------------------------------------------------------
    // Objective : Print the multiStore data in the output stream.
    //-------------------------------------------------------------------------

    s << '{' << domain() << endl ;

    //-------------------------------------------------------------------------
    // Write the size of each entity
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Write the data of each entity in the domain 
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      for(const T *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << ' ' ;
      s << std::endl ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // close the bracket ...
    //-------------------------------------------------------------------------

    s << '}' << std::endl ;

    return s ;
  }

  //*************************************************************************/

  template<class T> 
  std::istream &multiStoreRepI<T>::Input(std::istream &s) 
  {
    //-------------------------------------------------------------------------
    // Objective : Read the multiStore data from the input stream
    //-------------------------------------------------------------------------

    entitySet e ;
    char ch ;

    // Look for opening bracket ....
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    //-------------------------------------------------------------------------
    // Read the entitySet intervals ....
    //-------------------------------------------------------------------------
    
    s >> e ;

    //-------------------------------------------------------------------------
    // Read the size of each entity in the set ...
    //-------------------------------------------------------------------------
    
    store<int> sizes ;
    sizes.allocate(e) ;

    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // read the data
    //-------------------------------------------------------------------------
    
    allocate(sizes) ;
    FORALL(e,ii) {
      for(T *ip = begin(ii);ip!=end(ii);++ip)
        s >> *ip  ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Look for closing brackets
    //-------------------------------------------------------------------------
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //**************************************************************************/

  template <class T> 
  int multiStoreRepI<T>::pack_size(const entitySet &eset ) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet ecommon;
    ecommon = eset & domain();

    return get_mpi_size( traits_type, ecommon );
  }
  //**************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  void multiStoreRepI<T>::StringVal( const int &entity, const int &ivec, std::string &memento)
  {
    std::ostringstream oss;

    oss << base_ptr[entity][ivec] << endl;

    memento = oss.str();
  }
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(DEFAULT_CONVERTER c, const entitySet &eset ) 
  {
    int size;
    std::ostringstream oss;

    FORALL(eset,i) {
      size = base_ptr[i+1] - base_ptr[i] ;
      for(int j=0; j<size; j++)
        oss << base_ptr[i][j] << " ";
      oss << endl;
    } ENDFORALL ;

    std::string memento = oss.str();
    return(memento.length() );
  }
#endif
  //**************************************************************************/

  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset ) 
  {

    int size = 0 ;
    FORALL(eset,i) {
      size += end(i)- begin(i) ;
    } ENDFORALL ;

    size *= sizeof(T) ;
    size += eset.size()*sizeof(int) ;
    return(size) ;
  }

  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset ) 
  {
    int        arraySize =0, vsize, numContainers = 0;

    entitySet  :: const_iterator ci;
    typedef data_schema_traits<T> schema_traits; 

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize = end(*ci) - begin(*ci) ;
      numContainers  += vsize;
      for( int ivec = 0; ivec < vsize; ivec++){
          typename schema_traits::Converter_Type cvtr(base_ptr[*ci][ivec] );
          arraySize += cvtr.getSize() ;
      }
    }

    numContainers += eset.size();
    return( arraySize*sizeof(typename schema_traits::Converter_Base_Type) +
            numContainers*sizeof(int) );
  }
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::pack( void *outbuf, int &position, int &outcount, 
                                const entitySet &eset ) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    FORALL(eset,ii){
      size = end(ii) - begin(ii) ;
      MPI_Pack( &size, 1, MPI_INT, outbuf, outcount, &position, MPI_COMM_WORLD) ;
    }ENDFORALL ;

    packdata( traits_type, outbuf, position, outcount, eset);
  }
  //**************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void multiStoreRepI<T>::packdata( DEFAULT_CONVERTER c, void *outbuf,
                                    int &position,  int outcount,
                                    const entitySet &eset ) 
  {
    int  size;
    int  bufSize;
    std::ostringstream oss;
    std::string memento;

    FORALL(eset,ii){
      size = end(ii) - begin(ii) ;
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
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf,
                                    int &position,  int outcount,
                                    const entitySet &eset ) 
  {
    int  vsize, incount;
    FORALL(eset,i) {
      vsize   = end(i)-begin(i);
      incount = vsize*sizeof(T);
      MPI_Pack( &base_ptr[i][0], incount, MPI_BYTE, outbuf, outcount, &position, 
                MPI_COMM_WORLD) ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                    int &position, int outcount,
                                    const entitySet &eset ) 
  {

    entitySet::const_iterator ci;

    //-------------------------------------------------------------------------
    // Get the maximum size of container
    //-------------------------------------------------------------------------
    int vecsize, stateSize, maxStateSize=0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vecsize = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vecsize; ivec++){
        typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*ci][ivec]);
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
      vecsize = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vecsize; ivec++){
        typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*ci][ivec]);
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
void multiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, 
                               const sequence &seq) 
{
  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_type;
  if(base_ptr == 0) return ;

  entitySet new_dom = domain() | entitySet(seq) ;
  entitySet ent = domain() - entitySet(seq);

  store<int> ecount ;
  ecount.allocate(new_dom) ;

  bool conflict = 0 ;
  for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
    MPI_Unpack(ptr, size, &loc, &ecount[*si], 1, MPI_INT, MPI_COMM_WORLD) ;
    if(ecount[*si] != (base_ptr[*si+1] - base_ptr[*si])) conflict = 1 ;
  }

  if(conflict) {
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ; 
    multialloc(ecount, &new_index, &new_alloc_pointer, &new_base_ptr) ;
      
    for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei) {
      for(int j = 0 ; j < ecount[*ei]; ++j) 
        new_base_ptr[*ei][j] = base_ptr[*ei][j] ;
    }
      
    if(alloc_pointer) delete [] alloc_pointer ;
      
  
    alloc_pointer = new_alloc_pointer;
      
    if(index) delete[] index ;
      
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  unpackdata( traits_type, ptr, loc, size, seq); 
}

//**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T>::unpackdata( DEFAULT_CONVERTER c, void *inbuf, int &position, 
                                    int insize, const sequence &seq) 
{
  sequence:: const_iterator ci;
  char *outbuf;
  int   outcount, offset, vsize;
  entitySet eset(seq);

  for( ci = seq.begin(); ci != seq.end(); ++ci) {
    vsize = end(*ci)-begin(*ci);
    for( int ivec = 0; ivec < vsize; ivec++){
      MPI_Unpack( inbuf, insize, &position, &outcount, 1, MPI_INT, 
                  MPI_COMM_WORLD) ;
      outbuf   = new char[outcount];

      MPI_Unpack( inbuf, insize, &position, outbuf, outcount, MPI_BYTE, 
                  MPI_COMM_WORLD) ;
      std::istringstream iss(outbuf);
      iss >> base_ptr[*ci][ivec];
      delete [] outbuf;
    }
  }
}
#endif

//**************************************************************************/
template <class T> 
void multiStoreRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position,
                                    int insize, const sequence &seq) 
{
  int   outcount, vsize;
  sequence :: const_iterator si;

  for(si = seq.begin(); si != seq.end(); ++si) {
    vsize    = end(*si)-begin(*si);
    outcount = vsize*sizeof(T);
    MPI_Unpack(inbuf, insize, &position, &base_ptr[*si][0], outcount, MPI_BYTE,
               MPI_COMM_WORLD) ;
  }
}
  
//**************************************************************************/
template <class T> 
void multiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                    int &position, int insize,
                                    const sequence &seq) 
{
    int  stateSize, outcount, vsize;
    sequence :: const_iterator ci;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    std::vector<dtype> outbuf;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      vsize  = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vsize; ivec++) {
        MPI_Unpack( inbuf, insize, &position, &stateSize, 1,
                    MPI_INT, MPI_COMM_WORLD) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);

        outcount = stateSize*typesize;
        MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount,
                    MPI_BYTE, MPI_COMM_WORLD) ;

        typename schema_traits::Converter_Type  cvtr( base_ptr[*ci][ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }
}
//**************************************************************************/
  
template<class T> 
void multiStoreRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset) 
{
  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_type;

  entitySet eset;
  HDF5_ReadDomain(group_id, eset);

  entitySet dom = eset & user_eset ;


  //-------------------------------------------------------------------------
  // Size of each main container....
  //--------------------------------------------------------------------------
  hsize_t  dimension;

  hid_t vDatatype  = H5T_NATIVE_INT;
  hid_t vDataset   = H5Dopen(group_id,"ContainerSize");
  hid_t vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims (vDataspace, &dimension, NULL);

  std::vector<int> ibuf(dimension);
  H5Dread(vDataset,vDatatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  store<int> container;
  container.allocate( eset );

  size_t indx      = 0;
  entitySet::const_iterator ci;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[*ci] = ibuf[indx++];

  allocate( container);
  hdf5read( group_id, traits_type, eset, dom );

}

//**************************************************************************/

template<class T> 
void multiStoreRepI<T>::writehdf5(hid_t group_id, entitySet &usr_eset) const 
{

  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_output_type;

  entitySet eset = domain()&usr_eset ;

  if( eset.size() < 1) return;

  //write out the domain
  HDF5_WriteDomain(group_id, eset);

  hdf5write(group_id, traits_output_type, eset);

}

//**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
{
  hsize_t  dimension;
  hid_t vDatatype  = H5T_NATIVE_CHAR;
  hid_t vDataset   = H5Dopen(group_id,"VariableData");
  hid_t vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  std::vector<char> cbuf(dimension);

  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

  entitySet :: const_iterator ci;

  std::istringstream iss( &cbuf[0] );

  int count;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      iss >> base_ptr[*ci][j];
  }
  H5Dclose(vDataset);
  H5Sclose(vDataspace);
}
#endif
//**************************************************************************/

template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, IDENTITY_CONVERTER c, 
                                    entitySet &eset, entitySet &usr_eset )
{

  hsize_t dimension;
  hid_t   vDatatype, vDataset, vDataspace, mDataspace;
  size_t indx = 0, arraySize;
  int    rank = 1;

  entitySet::const_iterator ci;

  //-------------------------------------------------------------------------
  // Size of each main container....
  //--------------------------------------------------------------------------
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dopen(group_id,"ContainerSize");
  vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  int *ibuf = new int[dimension];
  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

  store<int> container;
  container.allocate( eset );

  indx  = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[*ci] = ibuf[indx++];

  delete [] ibuf;
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  //---------------------------------------------------------------------------
  // Calculate the offset of each entity in file ....
  //---------------------------------------------------------------------------
  store<int>        bucket;
  bucket.allocate( usr_eset );

  for( ci = usr_eset.begin(); ci != usr_eset.end(); ++ci) 
    bucket[*ci] = container[*ci];
  allocate( bucket );

  store<unsigned>   offset;
  offset.allocate( eset );

  arraySize = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    offset[*ci] = arraySize;
    arraySize  += container[*ci];
  }

  //---------------------------------------------------------------------------
  // Read the data now ....
  //---------------------------------------------------------------------------
  int num_intervals = usr_eset.num_intervals();
  interval *it = new interval[num_intervals];

  for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

  dimension  = arraySize;
  mDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataset   = H5Dopen(group_id,"VariableData");

  DatatypeP dtype = data_schema_traits<T>::get_type();
  vDatatype = dtype->get_hdf5_type();

  std::vector<T> data;

  hssize_t  start_mem[] = {0};  // determines the starting coordinates.
  hsize_t   stride[]    = {1};  // which elements are to be selected.
  hsize_t   block[]     = {1};  // size of element block;
  hssize_t  foffset[]   = {0};  // location (in file) where data is read.
  hsize_t   count[]     = {0};  // how many positions to select from the dataspace

  for( int k = 0; k < num_intervals; k++) {
    count[0] = 0;
    for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  container[i];

    if( count[0] > data.size() ) data.resize( count[0] );

    foffset[0] = offset[it[k].first];

    H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
    H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
    H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, &data[0]);

    indx = 0;
    for( int i = it[k].first; i <= it[k].second; i++) {
      for( int j = 0; j < container[i]; j++) 
        base_ptr[i][j] = data[indx++];
    }
  }
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);
  H5Sclose(mDataspace);
}

//**************************************************************************/

template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                    entitySet &eset, entitySet &usr_eset)
{
    hsize_t  dimension;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;

    int       indx=0, rank=1, vecsize, vsize;
    unsigned  arraySize;

    entitySet::const_iterator ci;

    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"ContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    std::vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------
    store<int> container;
    container.allocate( eset );
    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
       container[*ci] = ibuf[indx++];

    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"SubContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    ibuf.resize(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    //---------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------
    store< unsigned int>   offset;
    dmultiStore<int>  subContainer;
    offset.allocate( eset );

    arraySize = 0;
    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      vsize       = container[*ci];
      for( int i = 0; i < vsize; i++)  {
        vecsize    =  ibuf[indx++];
        arraySize  += vecsize;
        subContainer[*ci].push_back( vecsize );
      }
    }
    //---------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------
    int num_intervals = usr_eset.num_intervals();

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

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
        vsize = container[i];
        for( int j = 0; j < vsize; j++)
          count[0] +=  subContainer[i][j];
      }

      data = new typename data_schema_traits<T>::Converter_Base_Type[count[0]];

      foffset[0] = offset[it[k].first];
      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, data);

      indx = 0;
      int bucsize;
      for( int i = it[k].first; i <= it[k].second; i++) {
        vsize = container[i];
        for( int j = 0; j < vsize; j++) {
          typename data_schema_traits<T>::Converter_Type cvtr( base_ptr[i][j] );
          bucsize = subContainer[i][j];
          cvtr.setState( data+indx, bucsize );
          indx += bucsize ;
        }
      }
      delete[] data;
    }

}; 

//**************************************************************************/

template <class T> 
void multiStoreRepI<T>:: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, 
                                    const entitySet &eset)  const
{   
  entitySet :: const_iterator ci;
  hid_t     vDataset, vDataspace, vDatatype;

  typedef data_schema_traits<T> schema_traits;
    
  //------------------------------------------------------------------------
  // Write the container size:
  //------------------------------------------------------------------------
  int vsize, rank = 1, arraySize=0, stateSize;
  hsize_t  dimension;
  hid_t cparms     = H5Pcreate (H5P_DATASET_CREATE);

    
  std::vector<int>  container(eset.size());
  std::vector<int> bucketSize;

  size_t indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci){
    vsize           = end(*ci) - begin(*ci);
    container[indx++] = vsize;
    for( int ivec = 0; ivec < vsize; ivec++) {
       typename schema_traits::Converter_Type cvtr( base_ptr[*ci][ivec] );
       bucketSize.push_back(cvtr.getSize());
       arraySize  += cvtr.getSize();
    }
  }

  dimension  = eset.size();
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5T_NATIVE_INT;
  vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);

  dimension  = bucketSize.size();
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5T_NATIVE_INT;
  vDataset   = H5Dcreate(group_id, "SubContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
   H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bucketSize[0]);

   H5Dclose( vDataset  );
   H5Sclose( vDataspace);

  //-------------------------------------------------------------------------
  // Collect state data from each object and put into 1D array
  //-------------------------------------------------------------------------

   typedef typename schema_traits::Converter_Base_Type dtype;
   std::vector<dtype> data(arraySize);

   indx = 0;
   for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize  = end(*ci) - begin(*ci);
      for( int ivec = 0; ivec < vsize; ivec++){
        typename schema_traits::Converter_Type cvtr( base_ptr[*ci][ivec] );
        cvtr.getState( &data[0]+indx, stateSize);
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
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Tclose(vDatatype);
    H5Dclose(vDataset);
    H5Sclose(vDataspace);
};

//*************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T>::hdf5write( hid_t group_id, DEFAULT_CONVERTER g, 
                                   const entitySet &eset) const
{
  entitySet :: const_iterator ci;
  hid_t     vDataset, vDataspace, vDatatype;
    
  //------------------------------------------------------------------------
  // Write the container size:
  //------------------------------------------------------------------------
  int rank = 1;
  hsize_t  dimension;
  hid_t cparms     = H5Pcreate (H5P_DATASET_CREATE);

  dimension=  eset.size();
    
  std::vector<int>  container(eset.size());

  size_t indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[indx++] = end(*ci) - begin(*ci);

  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5T_NATIVE_INT;
  vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);

  std::ostringstream oss;
  int count;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      oss << base_ptr[*ci][j];
  }
  
  std::string memento = oss.str();
  hsize_t arraySize   =  memento.length();

  dimension  = arraySize+1;
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5T_NATIVE_CHAR;
  vDataset   = H5Dcreate( group_id, "VariableData", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, memento.c_str());

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);

}
#endif
//*************************************************************************/

template <class T> 
void multiStoreRepI<T>::hdf5write( hid_t group_id, IDENTITY_CONVERTER g, 
                                   const entitySet &eset) const
{

  hid_t  vDatatype, vDataset, vDataspace;

  entitySet :: const_iterator ci;

  //---------------------------------------------------------------------------
  // Get the sum of each object size and maximum size of object in the 
  // container for allocation purpose
  //---------------------------------------------------------------------------
  int     count, newsize;

  std::vector<int>  container(eset.size());

  size_t indx = 0, arraySize = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    newsize    = end(*ci) - begin(*ci);
    arraySize  += newsize;
    container[indx++] = newsize;
  }
         
  //------------------------------------------------------------------------
  // Write the Size of each multiStore ....
  //------------------------------------------------------------------------
  int rank = 1;
  hsize_t  dimension;
  hid_t cparms  = H5Pcreate (H5P_DATASET_CREATE);

  dimension  = eset.size();
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5T_NATIVE_INT;
  vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);

  //--------------------------------------------------------------------------
  // Collect state data from each object and put into 1D array
  //--------------------------------------------------------------------------

  std::vector<T>  data(arraySize);

  indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      data[indx++] = base_ptr[*ci][j];
  }

  //--------------------------------------------------------------------------
  // Write (variable) Data into HDF5 format
  //--------------------------------------------------------------------------
  typedef data_schema_traits<T> traits_type;
  DatatypeP dtype = traits_type::get_type();
  vDatatype = dtype->get_hdf5_type();

  dimension        = arraySize;
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataset   = H5Dcreate( group_id, "VariableData", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Tclose( vDatatype );
  
}; 

}

#endif


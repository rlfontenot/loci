#ifndef MULTISTORE_H
#define MULTISTORE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/tools.h>

#include <data_traits.h>
#include <store_rep.h>

#include <Tools/lmutex.h>

#include <Map.h>
#include <multiMap.h>
#include <dist_internal.h>

namespace Loci {

  template<class T> class multiStoreRepI : public storeRep {
    entitySet store_domain ;
    T **index ;
    T *alloc_pointer ;
    T **base_ptr ;
    int size ;
    lmutex mutex ;
    
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info&fi, entitySet &en);
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq );
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq);
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
  public:
    multiStoreRepI()
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0;}

    multiStoreRepI(const entitySet &p)
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; store_domain=p;}

    multiStoreRepI(const store<int> &sizes) {
      index = 0 ; alloc_pointer=0 ; base_ptr = 0; allocate(sizes) ; }

    void allocate(const store<int> &sizes) ;
    virtual void shift(int_type offset) ;
    void multialloc(const store<int> &count, T ***index, T **alloc_pointer, T ***base_ptr) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~multiStoreRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;

    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ;
  
  //*************************************************************************/
  
  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
    multiStore(const multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~multiStore() ;
    virtual void notification() ;
    
    
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
    const_multiStore(const const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(const multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(const const_multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~const_multiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
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

  template<class T>
    void multiStoreRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    allocate(store_domain) ;
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
  template<class T> 
    storeRep *multiStoreRepI<T>::new_store(const entitySet &p, const int* cnt) const {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new multiStoreRepI<T>(count) ;
  }
  //*************************************************************************/

  template<class T> 
  storeRepP multiStoreRepI<T>::remap(const dMap &m) const 
  {

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    multiStore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T>
  storeRepP multiStoreRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP multiStoreRepI<T>::thaw() {
    dmultiStore<T> ds ;
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei) {
      std::vector<T> v ;
      ds[*ei] = v ;
      int t = end(*ei) - begin(*ei) ;
      for(int i=0;i<t;++i)
        ds[*ei].push_back(base_ptr[*ei][i]) ;
    }
    return ds.Rep() ;
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
  void multiStoreRepI<T>::gather(const dMap &m, storeRepP &st,
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
  void multiStoreRepI<T>::scatter(const dMap &m, storeRepP &st,
                                  const entitySet &context) 
  {
    store<int> count;
    
    const_multiStore<T> s(st) ;
    count.allocate(domain());

    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    
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
      Loci::streamoutput(begin(ii),end(ii)-begin(ii),s) ;
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
        *ip = T() ;
      Loci::streaminput(begin(ii),end(ii)-begin(ii),s) ;
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
  inline int multiStoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset ) 
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
  int multiStoreRepI<T>::pack_size(const entitySet &eset ) 
  {
    fatal((eset - domain()) != EMPTY);
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    warn(eset-domain() != EMPTY) ;
    return get_mpi_size( traits_type, eset );
  }

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
  template<class T> 
    frame_info multiStoreRepI<T>::read_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return read_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info multiStoreRepI<T>::read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
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
    frame_info fi ;
    fi.is_stat = dim[0] ;
    fi.size = dim[1] ;
    int dom_size = domain().size() ;
    read_vector_int(group_id, "first_level", fi.first_level, dom_size) ;
    return fi ;
  }
  template<class T> 
    frame_info multiStoreRepI<T>::read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
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
    int dom_size = domain().size() ;
    read_vector_int(group_id, "first_level", fi.first_level, dom_size) ;
    int total_size = fi.first_level.size() ;
    int dims = 0 ;
    for(int i = 0; i < total_size; ++i)
      dims += (fi.first_level)[i] ;
    read_multi_vector_int(group_id, "second_level", dims, fi.second_level) ;
    return fi ;
  }
  
  template<class T> 
    frame_info multiStoreRepI<T>::write_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return write_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info multiStoreRepI<T>::write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 0 ;
    entitySet dom = domain() ;
    int newsize = 0 ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) {
      newsize = end(*ci) - begin(*ci);
      fi.first_level.push_back(newsize) ;
    }
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
    write_vector_int(group_id, "first_level", fi.first_level) ;
    return fi ;
  }
  template<class T> 
    frame_info multiStoreRepI<T>::write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    int vsize ;
    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 0 ;
    std::vector<int> fl(dom.size()) ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) {
      vsize = end(*ci) - begin(*ci) ;
      fi.first_level.push_back(vsize) ;
      for(int ivec = 0; ivec < vsize; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)][ivec]) ;
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
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
    write_vector_int(group_id, "first_level", fi.first_level) ;
    write_vector_int(group_id, "second_level", fi.second_level) ;
    return fi ;
  }
  template<class T> 
    DatatypeP multiStoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
    DatatypeP multiStoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
    DatatypeP multiStoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  
  template<class T> 
    void multiStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) 
    {
      typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
      schema_converter traits_output_type;
      hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, user_eset) ;
    }

  //**************************************************************************/

  template <class T> 
    void multiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
    {
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	T* tmp_array = new T[dimension] ;
	size_t tmp = 0 ;
	hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_array) ;
	if(err < 0) {
	  cerr << "H5Dread() failed" << endl ;
	}
	int loc = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
	  for(int ivec = 0; ivec < (fi.first_level)[loc]; ivec++) {
	    base_ptr[*si][ivec] = tmp_array[tmp++] ;
	    loc++ ;
	  }
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  
  //**************************************************************************/
  
  template <class T> 
    void multiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER, frame_info &fi, entitySet &eset)
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
	size_t tmp = 0 ;
	int bucsize ;
	size_t indx = 0 ;
	int loc = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
	  for(int ivec = 0; ivec < (fi.first_level)[loc]; ivec++) {
	    typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si][ivec]);
	    bucsize = vint[indx++] ;
	    cvtr.setState(tmp_array+tmp, bucsize) ;
	    tmp += bucsize ;
	  }
	  loc++ ;
	}
      
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }  
  //**************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset) ;
  }
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
    {
      if(dimension != 0) {
	storeRepP qrep = getRep() ;
	int rank = 1 ;
	DatatypeP dp = qrep->getType() ;
	hid_t datatype = dp->get_hdf5_type() ;
	hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	T* tmp_array = new T[dimension] ;
	size_t tmp = 0 ;
	int newsize = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
	  newsize = end(*si) - begin(*si);
	  for(int ivec = 0; ivec < newsize; ivec++){
	    tmp_array[tmp++] = base_ptr[*si][ivec] ;
	  }
	}
	H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }
    }
  //*************************************************************************/

  template <class T> 
    void multiStoreRepI<T>:: hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset)  const
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
	int newsize = 0 ;
	for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
	  newsize = end(*si) - begin(*si);
	  for(int ivec = 0; ivec < newsize; ivec++){
	    typename schema_traits::Converter_Type cvtr(base_ptr[*si][ivec]);
	    cvtr.getState(tmp_array+tmp, stateSize) ;
	    tmp +=stateSize ;
	  }
	}
	H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
	H5Sclose(memspace) ;
	delete [] tmp_array ;
      }   
    }

  //*************************************************************************/

}

#endif


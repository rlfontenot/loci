#ifndef MULTISTORE_H
#define MULTISTORE_H

#include <istream>
#include <ostream>

#include <algorithm>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>
#include <hdf5CC/H5cpp.h>

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

    void  hdf5read( H5::Group group, DEFAULT_CONVERTER c,      entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void  hdf5write( H5::Group group, DEFAULT_CONVERTER c,      const entitySet &en) const;
    void  hdf5write( H5::Group group, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( H5::Group group, USER_DEFINED_CONVERTER c, const entitySet &en) const;
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
    virtual void readhdf5( H5::Group group, entitySet &en) ;
    virtual void writehdf5( H5::Group group, entitySet& en) const ;

    bool is_static() { return istat ; } 
    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
  } ;
  
  //***************************************************************************/
  
  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore(storeRepP &rp) { setRep(rp) ;}
    
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
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  //***************************************************************************/

  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const multiStore<T> &m)
    { return m.Print(s) ; }

  //***************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, multiStore<T> &m)
    { return m.Input(s) ; }
 
  //***************************************************************************/

  template<class T> 
  multiStore<T>::~multiStore() {}

  //***************************************************************************/
  
  template<class T> 
  void multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //***************************************************************************/
  
  template<class T> class const_multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(storeRepP &rp) { setRep(rp) ;}
    
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

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //***************************************************************************/

  template<class T> 
  store_instance::instance_type
    const_multiStore<T>::access() const
    { return READ_ONLY ; }

  //***************************************************************************/
  
  template<class T> 
  const_multiStore<T>::~const_multiStore() {}

  //***************************************************************************/
  
  template<class T> 
  void const_multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //***************************************************************************/

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

  //***************************************************************************/

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

  //***************************************************************************/
   
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

  //***************************************************************************/
  
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

  //***************************************************************************/

  template<class T> 
  multiStoreRepI<T>::~multiStoreRepI() 
  {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
  }

  //***************************************************************************/

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

  //***************************************************************************/

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

  //***************************************************************************/
  
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

  //***************************************************************************/
  
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

  //***************************************************************************/

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

  //***************************************************************************/
  template <class T> int multiStoreRepI<T>::pack_size(const entitySet &e ) {
    int size = 0 ;
    FORALL(e,i) {
      size += base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    size *= sizeof(T) ;
    size += e.size() * sizeof(int) ;
    return(size) ;
  }
  
  //***************************************************************************/
  
  template <class T> 
  void multiStoreRepI<T>::pack(void * ptr, int &loc, int &size, const entitySet &e ) 
  {
/*
    //-------------------------------------------------------------------------
    // Objective : Pack the multiStore data into a 1D array, which is passed to
    // MPI for packing the data.
    //-------------------------------------------------------------------------
    
  */
    store<int> count ;
    count.allocate(e) ;
    FORALL(e,i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    FORALL(e,i) {
      MPI_Pack(&count[i], sizeof(int), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ; 
      MPI_Pack(&base_ptr[i][0], count[i] * sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
    } ENDFORALL ;
  }
  
  //***************************************************************************/
  
  template <class T> 
  void multiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    if(base_ptr == 0)
      return ;
    store<int> count ;
    bool conflict = 0 ;
    int temploc = loc ;
    entitySet new_dom = domain() | entitySet(seq) ;
    entitySet ent = domain() - entitySet(seq);
    count.allocate(new_dom) ;
    for(entitySet::const_iterator ei = domain().begin(); ei != domain().end(); ++ei)
      count[*ei] = base_ptr[*ei+1] - base_ptr[*ei] ;
    for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
      MPI_Unpack(ptr, size, &loc, &count[*si],sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
      if(count[*si] != (base_ptr[*si+1] - base_ptr[*si]))
	conflict = 1 ;
      loc += count[*si] * sizeof(T) ;
    }
    if(conflict) {
      T **new_index ;
      T *new_alloc_pointer ;
      T **new_base_ptr ; 
      multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
      for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei) {
	for(int j = 0 ; j < count[*ei]; ++j) 
	  new_base_ptr[*ei][j] = base_ptr[*ei][j] ;
      }
      if(alloc_pointer) delete [] alloc_pointer ;
      alloc_pointer = new_alloc_pointer;
      if(index) delete[] index ;
      index = new_index ;
      base_ptr = new_base_ptr ;
      dispatch_notify() ;
    }
    loc = temploc ;
    
    for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
      loc += sizeof(int) ;
      MPI_Unpack(ptr, size, &loc, &base_ptr[*si][0],count[*si]*sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
    }
  }
  
  //***************************************************************************/
  
  template<class T> 
  store_type multiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //***************************************************************************/
  
  template<class T> 
  entitySet multiStoreRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //***************************************************************************/
  
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

  //***************************************************************************/

  template<class T> 
  std::istream &multiStoreRepI<T>::Input(std::istream &s) 
  {
    cout << " Commented for the time being " << endl;
    exit(0);

    /*
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
    */
    return s ;
  }

  //**************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::readhdf5( H5::Group group, entitySet &user_eset) 
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset;
    HDF5_ReadDomain(group, eset);

    entitySet dom = eset & user_eset ;
    allocate( dom );
    
    hdf5read( group, traits_type, eset, dom );

  }

  //**************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::writehdf5(H5::Group group, entitySet &en) const 
  {

   typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
   schema_converter traits_output_type;

   hdf5write(group, traits_output_type, en );

  }

  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T> :: hdf5read( H5::Group group, DEFAULT_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
  {
      cout << "Fatal: Default read converter  not implemented for multiStore" << endl;
      exit(0);

  }
  //**************************************************************************/

  template <class T> 
  void multiStoreRepI<T> :: hdf5read( H5::Group group, IDENTITY_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet::const_iterator ci;

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------

    H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   bDataset   = group.openDataSet( "ContainerSize");
    H5::DataSpace bDataspace = bDataset.getSpace();

    bDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

    dimension[0]  = eset.size();
	 bDataset.read( ibuf, H5::PredType::NATIVE_INT );

    store<int> container;
    container.allocate( eset );

    indx  = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
        container[*ci] = ibuf[indx++];

    delete [] ibuf;

   //---------------------------------------------------------------------------
   // Calculate the offset of each entity in file ....
   //---------------------------------------------------------------------------
   store<int>        bucket;
   bucket.allocate( user_eset );
   for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
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
   int num_intervals = user_eset.num_intervals();
   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

   T  *data;

   dimension[0] = arraySize;
   H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
   H5::DataSpace vDataspace(rank, dimension);

   typedef hdf5_schema_traits<T> traits_type;
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
                count[0] +=  container[i];

        data = new T[count[0]];

        foffset[0] = offset[it[k].first];

        mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
        vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        vDataset.read( data, vDatatype, mDataspace, vDataspace);

        indx = 0;
        for( int i = it[k].first; i <= it[k].second; i++) {
             for( int j = 0; j < container[i]; j++) 
                  base_ptr[i][j] = data[indx++];
        }

        delete[] data;
   }
}

  //**************************************************************************/

  template <class T> 
  void multiStoreRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet::const_iterator ci;

    typedef hdf5_schema_converter_traits<T> converter_traits; 

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------

    H5::DataType  bDatatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   bDataset   = group.openDataSet( "ContainerSize");
    H5::DataSpace bDataspace = bDataset.getSpace();

    bDataspace.getSimpleExtentDims( dimension, NULL);
    int *ibuf = new int[dimension[0]];

    dimension[0]  = eset.size();
	 bDataset.read( ibuf, H5::PredType::NATIVE_INT );

    store<int> container;
    container.allocate( eset );

    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
        container[*ci] = ibuf[indx++];

    delete [] ibuf;

    //---------------------------------------------------------------------------
    // Size of each sub-container ....
    //---------------------------------------------------------------------------

    H5::DataType  sDatatype  = H5::PredType::NATIVE_INT;
	 H5::DataSet   sDataset   = group.openDataSet( "SubContainerSize");
	 H5::DataSpace sDataspace = sDataset.getSpace();

	 sDataspace.getSimpleExtentDims( dimension, NULL);
	 ibuf = new int[dimension[0]];

	 sDataset.read( ibuf, H5::PredType::NATIVE_INT );

    int maxBucketSize = *std::max_element( ibuf, ibuf + (int)dimension[0] );

   //---------------------------------------------------------------------------
   // Calculate the offset of each entity in file ....
   //---------------------------------------------------------------------------
   store<unsigned>   offset;
   store<int>        bucket;
   dmultiStore<int>  subcontainer;

   offset.allocate( eset );

   bucket.allocate( user_eset );
   for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
        bucket[*ci] = container[*ci];
   allocate( bucket );

   arraySize = 0;
   int indx1 = 0, indx2 = 0;
   for( ci = eset.begin(); ci != eset.end(); ++ci) {
        offset[*ci] = arraySize;
        for( int i = 0; i < container[*ci]; i++)  {
             size = ibuf[indx2];
             arraySize  += size;
             subcontainer[*ci].push_back( size );
             indx2++;
        }
   }

   delete [] ibuf;
   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------
   int num_intervals = user_eset.num_intervals();
   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

   typedef hdf5_schema_converter_traits<T> converter_traits; 
   converter_traits::memento_type *data, *buf;

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
             for( int j = 0; j < subcontainer[i].size(); j++)
                  count[0] +=  subcontainer[i][j];
        }

        data = new typename converter_traits::memento_type[count[0]];

        foffset[0] = offset[it[k].first];

        mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
        vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        vDataset.read( data, vDatatype, mDataspace, vDataspace);

        indx = 0;
        int size;
        for( int i = it[k].first; i <= it[k].second; i++) {
             for( int j = 0; j < subcontainer[i].size(); j++) {
                  Memento<T> memento( base_ptr[i][j] );
                  size = subcontainer[i][j];
                  for( int m = 0; m < size; m++)
                       buf[m] = data[indx++];
                  base_ptr[i][j] = memento.setState( buf, size );
             }
        }

        delete[] data;
   }

   delete[] buf;

  }; 

  //**************************************************************************/

  template <class T> 
  void multiStoreRepI<T>:: hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, 
                                       const entitySet &eset)  const
  {   

    int rank = 1;
    hsize_t  dimension[1];

    //write out the domain   
    HDF5_WriteDomain(group, eset);

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------

    entitySet :: const_iterator ci;
    size_t       arraySize= 0;
    int          count, stateSize, *storeSize, maxBucketSize;
    std::vector<int>  bucketSize;    //  Because we don't know in advance the size

    storeSize = new int[eset.size()];

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         count  = end(*ci) - begin(*ci);
         storeSize[indx++] = count;

         for( int j = 0; j < count; j++) {
             Memento<T> memento( base_ptr[*ci][j] );
             stateSize  = memento.getSize();
             arraySize += stateSize;
             bucketSize.push_back( stateSize );
         }
    }

    maxBucketSize = *max_element( bucketSize.begin(), bucketSize.end() );

//-----------------------------------------------------------------------------
// Write the Size of each multiStore ....
//-----------------------------------------------------------------------------
    dimension[0]=  eset.size();

    try {
      H5::DataSpace sDataspace( rank, dimension );
      H5::DataType  sDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   sDataset  = group.createDataSet( "ContainerSize", sDatatype, sDataspace);

      sDataset.write( storeSize, sDatatype );
    }

    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

//-----------------------------------------------------------------------------
// Write the size of each bucket...
//-----------------------------------------------------------------------------
    dimension[0]=  bucketSize.size();
    int  *bucket = new int[bucketSize.size()];

    for( int i=0; i< bucketSize.size();i++)
         bucket[i] = bucketSize[i];
    bucketSize.clear();

    try {
      H5::DataSpace bDataspace( rank, dimension );
      H5::DataType  bDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   bDataset  = group.createDataSet( "SubContainerSize", bDatatype, bDataspace);

      bDataset.write( bucket, bDatatype );
    }

    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    converter_traits::memento_type *data, *buf;

 	 data =  new typename converter_traits::memento_type[arraySize];
 	 buf  =  new typename converter_traits::memento_type[maxBucketSize];
//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         count  = end(*ci) - begin(*ci);
         for( int j = 0; j < count; j++) {
              Memento<T> memento( base_ptr[*ci][j] );
              memento.getState(buf, stateSize);
              for( int i = 0; i < stateSize; i++)
                   data[indx++] =  buf[i];
         }
    }

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------

    dimension[0] =  arraySize;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "variable", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

//-----------------------------------------------------------------------------
// Clean up
//-----------------------------------------------------------------------------
    delete [] data;
    delete [] buf;

  };

  //*************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::hdf5write( H5::Group group, DEFAULT_CONVERTER g, 
                                     const entitySet &eset) const
  {
    cout << "Fatal: default converter write not implemented yet for multiStore" << endl;
    exit(0);
  }
  //*************************************************************************/

  template <class T> 
  void multiStoreRepI<T>::hdf5write( H5::Group group, IDENTITY_CONVERTER g, 
                                     const entitySet &eset) const
  {

/*
    //write out the domain   
    domain_hdf5write(group, eset);

    entitySet :: const_iterator ci;

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------
    size_t  arraySize= 0;
    int     count;

    for( ci = eset.begin(); ci != eset.end(); ++ci)
         arraySize  += end(*ci) - begin(*ci);

    T  *data, *buf;

 	 data =  new T[arraySize];
//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         count  = end(*ci) - begin(*ci);
         for( int j = 0; j < count; j++) 
              data[indx++] = base_ptr[*ci][j];
    }

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
      H5::DataSet   vDataset  = group.createDataSet( "variable", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

//-----------------------------------------------------------------------------
// Clean up
//-----------------------------------------------------------------------------
    delete [] data;
*/

  }; 

  //*************************************************************************/
}

#endif


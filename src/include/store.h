#ifndef STORE_H
#define STORE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Map.h>
#include <store_rep.h>
#include <data_traits.h>
#include <sstream>
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
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &usr);
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void packdata(IDENTITY_CONVERTER, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;

    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en) const;
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
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
    storeRepI() { alloc_pointer = 0 ; base_ptr = 0; }
    storeRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ;}
    virtual void allocate(const entitySet &ptn) ;
    virtual ~storeRepI()  ;
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
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    virtual entitySet domain() const ;
    T * get_base_ptr() const { return base_ptr ; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ;

  template<class T> void storeRepI<T>::allocate(const entitySet &eset) {
    // if the pass in is EMPTY, we delete the previous allocated memory
    // this equals to free the memory
    if( eset == EMPTY ) {
      delete [] alloc_pointer ;
      alloc_pointer = 0 ; base_ptr = 0;
      store_domain = eset ;
      dispatch_notify() ;
      return ;
    }

    int_type old_range_min = store_domain.Min() ;
    int_type old_range_max = store_domain.Max() ;
    int_type new_range_min = eset.Min() ;
    int_type new_range_max = eset.Max() ;

    // if the old range and the new range are equal, nothing
    // needs to be done, just return
    if( (old_range_min == new_range_min) &&
        (old_range_max == new_range_max)) {
      store_domain = eset ;
      return ;
    }

    // is there any overlap between the old and the new domain?
    // we copy the contents in the overlap region to the new
    // allocated storage
    entitySet ecommon = store_domain & eset ;
    
    T* tmp_alloc_pointer = new T[new_range_max - new_range_min + 1] ;
    T* tmp_base_ptr = tmp_alloc_pointer - new_range_min ;

    // if ecommon == EMPTY, then nothing is done in the loop
    FORALL(ecommon,i) {
      tmp_base_ptr[i] = base_ptr[i] ;
    } ENDFORALL ;
    
    delete [] alloc_pointer ;
    alloc_pointer = tmp_alloc_pointer ;
    base_ptr = tmp_base_ptr ;
      
    store_domain = eset ;
    dispatch_notify() ;
    return ;
  }


  template<class T> std::ostream &storeRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      Loci::streamoutput(&base_ptr[ii],1,s) ;
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
      base_ptr[ii] = T() ;
      Loci::streaminput(&base_ptr[ii],1,s) ;
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
  template<class T>
  storeRep *storeRepI<T>::new_store(const entitySet &p, const int* cnt) const {
    storeRep* sp ;
    cerr << " This method should not be called for a store " << endl ;
    return sp ;
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
    store(const store &var) { setRep(var.Rep()) ; }
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
        
  template<class T> storeRepP storeRepI<T>::remap(const dMap &m) const {
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
  template<class T> void storeRepI<T>::gather(const dMap &m, storeRepP &st,
                                              const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T> void storeRepI<T>::scatter(const dMap &m, storeRepP &st,
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

  template <class T>
  inline int storeRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                  const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }

  //*******************************************************************/
  template <class T>
  inline int storeRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                  const entitySet &eset)
  {

    int       size ;
    int numBytes = 0 ;
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
  template <class T> int storeRepI<T>::pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }


  //*******************************************************************/
  template <class T> 
  inline void storeRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
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
  inline void storeRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
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
  void storeRepI<T>::pack( void *outbuf, int &position, int &size, 
                           const entitySet &usr_eset )  
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    warn(usr_eset-domain() != EMPTY) ;

    packdata( schema_converter(), outbuf, position, size, usr_eset);

  }
 
  //*********************************************************************/
  template <class T> 
  inline void storeRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
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
  inline void storeRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
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


  //*******************************************************************/

  template <class T> 
  void storeRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata(schema_converter(), ptr, loc, size, seq); 
  }

  //**********************************************************************/
 
  
  template<class T> 
    frame_info storeRepI<T>::read_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return read_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info storeRepI<T>::read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
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
    frame_info storeRepI<T>::read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    hsize_t dimension = 0 ;
    hid_t dataspace ;
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
      dataset = H5Dopen(group_id, "second_level") ;
      dataspace = H5Dget_space(dataset) ;
      H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    }
    int dim[2] ;
    dim[0] = is_stat ;
    dim[1] = sz ;
    MPI_Bcast(&dim, 2, MPI_INT, 0, MPI_COMM_WORLD) ;
    fi.is_stat = dim[0] ;
    fi.size = dim[1] ;
    std::vector<int> vint ;
    read_vector_int(group_id, "second_level", vint) ;
    fi.second_level = vint ; 
    return fi ;
  }
  
  template<class T> 
    frame_info storeRepI<T>::write_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return write_frame_info(group_id, schema_converter()) ;
  }
  template<class T> 
    frame_info storeRepI<T>::write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
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
    frame_info storeRepI<T>::write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
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
    DatatypeP storeRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
    DatatypeP storeRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
    DatatypeP storeRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  
  template<class T> 
    void storeRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }

  //**********************************************************************/
  template <class T> 
  void storeRepI<T> :: hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset)  const
    {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
	tmp_array[tmp++] = base_ptr[*si] ;
      
      
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  //*********************************************************************/

  template <class T> 
  void storeRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
    {   
      typedef data_schema_traits<T> schema_traits ;
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
	typename schema_traits::Converter_Type cvtr(base_ptr[*si]);
	cvtr.getState(tmp_array+tmp, stateSize) ;
	tmp +=stateSize ;
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  
  
  //*********************************************************************/
 template<class T> 
  void storeRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }
  template<class T>
    void  storeRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
    {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			    H5P_DEFAULT, tmp_array) ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
	base_ptr[*si] = tmp_array[tmp++] ;
      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  //*********************************************************************/

  template<class T>
    void  storeRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits ;
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
    size_t tmp = 0 ;
    int bucsize ;
    size_t indx = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) { 
      typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si]);
      bucsize = vint[indx++] ;
      cvtr.setState(tmp_array+tmp, bucsize) ;
      tmp += bucsize ;
    }
    H5Sclose(memspace) ;
    delete [] tmp_array ;
  }
  //*********************************************************************/
}

#endif

#ifndef DMULTISTORE_H
#define DMULTISTORE_H

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <hdf5_readwrite.h>

#include <store_rep.h>
#include <Tools/lmutex.h>
#include <storeVec.h>
#include <Map.h>
#include <multiMap.h>

#include <DMultiMap.h>

#include <algorithm>

namespace Loci {
  template<class T> class dmultiStoreRepI : public storeRep {
    entitySet                 store_domain ;
    hash_map<int,std::vector<T> >  attrib_data;

    void  hdf5read( hid_t group, IDENTITY_CONVERTER c,     entitySet &en, entitySet &usr);
    void  hdf5read( hid_t group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr);

    void  hdf5write( hid_t group, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( hid_t group, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size, const entitySet &e ) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size, const entitySet &e ) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size, const sequence &seq) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size, const sequence &seq) ;

  public:

    //  Constructor ...
    dmultiStoreRepI() { }

    dmultiStoreRepI(const entitySet &p) { store_domain=p;}

    dmultiStoreRepI(const store<int> &sizes) { allocate(sizes) ; }

    //  Destructors ...
    virtual ~dmultiStoreRepI() ;

    //  Member function ...
    void allocate(const store<int> &sizes) ;

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
    virtual void readhdf5( hid_t group, entitySet &user_eset) ;
    virtual void writehdf5( hid_t group, entitySet& en) const ;

    hash_map<int,std::vector<T> > *get_attrib_data(){return &attrib_data; }
  } ;

  //***************************************************************************/
  
  template<class T> class dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T>  storeType ;
    hash_map<int, std::vector<T> > *attrib_data;
  public:
    typedef std::vector<T> containerType ;
    dmultiStore() {setRep(new storeType) ;}
    dmultiStore(dmultiStore<T> &var) {setRep(var.Rep()) ;}
    dmultiStore(storeRepP &rp) { setRep(rp) ;}
    
    virtual ~dmultiStore() ;
    virtual void notification() ;
    
    dmultiStore<T> & operator=(dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    
    dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }

    void setSizes(const const_dmultiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }

    entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    std::vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  //*************************************************************************/

  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const dmultiStore<T> &m)
  { return m.Print(s) ; }

  //*************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dmultiStore<T> &m)
  { return m.Input(s) ; }

  //**************************************************************************/
 
  template<class T> 
  dmultiStore<T>::~dmultiStore() {}

  //**************************************************************************/
  
  template<class T> 
  void dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //**************************************************************************/
  
  template<class T> class const_dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T> storeType ;
    hash_map<int, std::vector<T> >   *attrib_data;
  public:
    //    typedef const_Vect<T> containerType ;
    const_dmultiStore() {setRep(new storeType) ;}
    const_dmultiStore(const_dmultiStore<T> &var) {setRep(var.Rep()) ;}
    const_dmultiStore(storeRepP &rp) { setRep(rp) ;}
    
    virtual ~const_dmultiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
    const_dmultiStore<T> & operator=(const dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dmultiStore<T> & operator=(const const_dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    const entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }

    std::vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //**************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_dmultiStore<T>::access() const
  { return READ_ONLY ; }

  //**************************************************************************/
  
  template<class T> 
  const_dmultiStore<T>::~const_dmultiStore() {}

  // *************************************************************************
  
  template<class T> 
  void const_dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::allocate(const store<int> &sizes) 
  {

    //------------------------------------------------------------------------/
    // Objective : reserve the memory of multiStore using Store..
    //------------------------------------------------------------------------/

    entitySet eset = sizes.domain() ;
    entitySet :: const_iterator  ci;
    int   veclength;
    T     newObj;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      veclength = sizes[*ci];
      attrib_data[*ci].reserve(veclength);
      for( int i = 0; i < veclength; i++)
        attrib_data[*ci].push_back( newObj );
    }

    store_domain = eset ;
    dispatch_notify() ;
  }
  //***************************************************************************/
  
  template<class T> 
  void dmultiStoreRepI<T>::allocate(const entitySet &ptn) 
  {
    std::vector<T>   emptyVec;
    entitySet :: const_iterator  ci;

    for( ci = ptn.begin(); ci != ptn.end(); ++ci)
      attrib_data[*ci] = emptyVec;

    store_domain = ptn ;
    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> 
  dmultiStoreRepI<T>::~dmultiStoreRepI() 
  {
    attrib_data.clear();
  }

  //**************************************************************************/

  template<class T> 
  storeRep *dmultiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dmultiStoreRepI<T>(p) ;
  }

  //**************************************************************************/

  template<class T> 
  storeRepP dmultiStoreRepI<T>::remap(const Map &m) const {
    dmultiStore<T> s ;
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    /*
    multiStore<T> static_mul ;
    store<int> count ;
    entitySet tmp_dom = domain() ;
    count.allocate(tmp_dom) ;
    for(entitySet::const_iterator ei = tmp_dom.begin(); ei != tmp_dom.end(); ++ei)
      count[*ei] = s[*ei].size() ;
    for(entitySet::const_iterator ei = tmp_dom.begin(); ei != tmp_dom.end(); ++ei) {
      int i = 0 ;
      for(std::vector<int>::const_iterator vi = s[*ei].begin(); vi != s[*ei].end(); ++vi) {
	static_mul[*ei][i] = *vi ;
	++i ;
      }
    }
    return static_mul.Rep() ;
    */
    return s.Rep() ;
  }

  //**************************************************************************/
  
  template<class T> 
  void dmultiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    std::vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[i].reserve( newVec.size() );
      for( int j = 0; j < newVec.size(); j++)
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::gather(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    std::vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[m[i]];
      attrib_data[i].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[i].push_back( newVec[i] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::scatter(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    std::vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[m[i]].reserve( newVec.size() );
      for( int i = 0; i < newVec.size(); i++) 
        attrib_data[m[i]].push_back( newVec[i] );
    } ENDFORALL ;

    dispatch_notify() ;
  }

  //**************************************************************************/
 
  template<class T> 
  store_type dmultiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //***************************************************************************/
  
  template<class T> 
  entitySet dmultiStoreRepI<T>::domain() const 
  {
    hash_map<int,std::vector<T> > :: const_iterator    ci;
    entitySet          storeDomain;
    std::vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;

    sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];

    return storeDomain ;
  }

  //***************************************************************************/
  
  template<class T> 
  std::ostream &dmultiStoreRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << endl ;

    hash_map<int,std::vector<T> >  :: const_iterator ci;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci != attrib_data.end() ) continue;
      s << ci->second.size() << std::endl ;
    } ENDFORALL ;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci == attrib_data.end() ) continue;
      const std::vector<T> &vec =  ci->second;
      streamoutput(&vec[0],vec.size(),s) ;
    } ENDFORALL ;

    s << '}' << std::endl ;
    return s ;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &dmultiStoreRepI<T>::Input(std::istream &s) 
  {

    entitySet  e ;
    char ch ;

    // Read the opening brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    // Read the interval set ....
    s >> e ;

    // Read the size of each entity map.
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;
    allocate(sizes) ;

    // Read the attribute data
    std::vector<T>  vec;
    FORALL(e,ii) {
      vec.clear();
      for( int i = 0; i < sizes[ii]; i++) {
        T val ;
        streaminput(&val,1,s) ;
        vec.push_back(val);
      }
      attrib_data[ii] = vec;
    } ENDFORALL ;
            
    // Close the bracket ..
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;

  }

  //***************************************************************************/
  template <class T> 
  int dmultiStoreRepI<T>::pack_size(const entitySet &e ) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }

  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    int size = 0 ;
    FORALL(eset,i) {
      size  +=  attrib_data[i].size();
    } ENDFORALL ;
    
    return( size*sizeof(T) + eset.size()*sizeof(int) ) ;
  }

  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    std::vector<T>   newVec;
    entitySet :: const_iterator  ei;
    hash_map<int,std::vector<T> >:: const_iterator ci;

    int     arraySize =0, numContainers = 0;
    typedef data_schema_traits<T> schema_traits;

    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newVec      = ci->second;
        numContainers += newVec.size();
        for( int ivec = 0; ivec < newVec.size(); ivec++) {
          typename schema_traits::Converter_Type cvtr(newVec[ivec] );
          arraySize += cvtr.getSize() ;
        }
      }
    }

    int  vsize = arraySize*sizeof(typename schema_traits::Converter_Base_Type) + 
      numContainers*sizeof(int)  +   // size of each object
      eset.size()*sizeof(int);       // size of each entity

    return(vsize);
  }

  //**************************************************************************/
  
  template <class T> 
  void dmultiStoreRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &eset ) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
  
    packdata( traits_type, ptr, loc, size, eset);
  }
  
  //**************************************************************************/
  template <class T>
  void dmultiStoreRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                     int outcount, const entitySet &eset )
  {
    int vsize, incount;
    std::vector<T>   inbuf;
    entitySet :: const_iterator   ci;
    hash_map<int,std::vector<T> >::const_iterator iter;

    for( ci = eset.begin(); ci != eset.end(); ++ci){
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        inbuf   = iter->second;
        vsize   = inbuf.size();
        incount = vsize*sizeof(T);
        MPI_Pack( &vsize, 1, MPI_INT, outbuf, outcount, &position, 
                  MPI_COMM_WORLD) ;
        MPI_Pack( &inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
                  MPI_COMM_WORLD) ;
      } 
    }

  }

  //**************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                     int &position, int outcount, 
                                     const entitySet &eset ) 
  {
    entitySet::const_iterator ci;
    int vsize, stateSize;
    
    typedef data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);
    std::vector<dtype> inbuf;

    int incount;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize = attrib_data[*ci].size();
      MPI_Pack(&vsize,1, MPI_INT, outbuf, outcount, &position,
               MPI_COMM_WORLD);
      for( int ivec = 0; ivec < vsize; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( attrib_data[*ci][ivec] );

        stateSize = cvtr.getSize();
        inbuf.resize(stateSize);

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
  void dmultiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet ecommon, ediff,eset(seq);

    ediff = eset - domain();

    if( ediff.size() > 0)
      cout << " Warning: Entities not part of domain, and not unpacked " << ediff <<endl;

    unpackdata( traits_type, ptr, loc, size, seq);
  }
  

  //**************************************************************************/
  template <class T> 
  void dmultiStoreRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position, 
                                       int &insize, const sequence &seq)
  {

    int   outcount, vsize;
    std::vector<T>  outbuf;
    sequence :: const_iterator ci;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      MPI_Unpack( inbuf, insize, &position, &vsize, 1, MPI_INT, 
                  MPI_COMM_WORLD) ;

      if( vsize > outbuf.size() ) outbuf.resize(vsize);

      attrib_data[*ci].resize(vsize);
      outcount = vsize*sizeof(T);
      MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount, 
                  MPI_BYTE, MPI_COMM_WORLD) ;
      for( int ivec = 0; ivec < vsize; ivec++) 
        attrib_data[*ci][ivec] = outbuf[ivec];
    }
  }

  //***********************************************************************/
  template <class T> 
  void dmultiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                       int &position, int &insize, const sequence &seq)
  {

    sequence :: const_iterator ci;
    int  stateSize, outcount, vsize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    std::vector<dtype> outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( !store_domain.inSet( *ci ) ) {
        cout << "Warning: Entity not present in entityset " << *ci << endl;
        continue;
      }
      MPI_Unpack( inbuf, insize, &position, &vsize, 1, MPI_INT, MPI_COMM_WORLD) ;
      attrib_data[*ci].resize(vsize);

      for( int ivec = 0; ivec < vsize; ivec++) {
        MPI_Unpack( inbuf, insize, &position, &stateSize, 1, 
                    MPI_INT, MPI_COMM_WORLD) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);

        outcount = stateSize*typesize;
        MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount, 
                    MPI_BYTE, MPI_COMM_WORLD) ;

        typename schema_traits::Converter_Type  cvtr( attrib_data[*ci][ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }

  }
 
  //***************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet eset, ecommon;

    // Read the entitySet available in file ...
    HDF5_ReadDomain(group_id, eset);

    // Intersection of entityset in file and user defined. Only common entitySet
    // are read from the file ...
    //
    ecommon = eset & user_eset ;   
    allocate( ecommon );

    // Read the common entitities ...
    hdf5read( group_id, traits_type, eset,  ecommon );
  }

  //*************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T> :: hdf5read( hid_t group_id, IDENTITY_CONVERTER c, 
                                       entitySet &eset, entitySet &user_eset)
  {
    hsize_t  dimension;
    hid_t    vDataset, vDataspace, vDatatype, mDataspace;

    int indx = 0, arraySize, vsize;
    int    rank = 1, vecsize;

    entitySet::const_iterator ci;

    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dopen(group_id,"ContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    std::vector<int> ibuf(dimension);
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

    //---------------------------------------------------------------------
    // Size of each main container....
    //---------------------------------------------------------------------

    store<int> container;
    container.allocate( eset );

    indx      = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      container[*ci] = ibuf[indx++];

    container.Print(cout);

    //-----------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //-----------------------------------------------------------------------
    store<unsigned>   offset;
    store<int>        bucket;

    offset.allocate( eset );

    bucket.allocate( user_eset );
    for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
      bucket[*ci] = container[*ci];
    allocate( bucket );

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += container[*ci];
    }
    //------------------------------------------------------------------------
    // Read the data now ....
    //------------------------------------------------------------------------
    int num_intervals = user_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

    dimension = arraySize;

    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    vDatatype = dtype->get_hdf5_type();

    mDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dopen( group_id, "VariableData");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    int voffset;

    std::vector<T>  data;

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  container[i];

      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        vsize = container[i];
        attrib_data[i].resize(vsize);
        for( int ivec = 0; ivec < vsize; ivec++)
          attrib_data[i][ivec] = data[indx++];
      }
    }

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Sclose( mDataspace);
    H5Tclose( vDatatype );

  }

  //*************************************************************************/
  template<class T>
  void  dmultiStoreRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                        entitySet &eset, entitySet &user_eset)
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

    //-----------------------------------------------------------------------
    // Size of each main container....
    //------------------------------------------------------------------------
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
        attrib_data[i].resize(vsize);
        for( int j = 0; j < vsize; j++) {
          typename data_schema_traits<T>::Converter_Type cvtr( attrib_data[i][j] );
          bucsize = subContainer[i][j];
          cvtr.setState( data+indx, bucsize );
          indx += bucsize ;
        }
      }

      delete[] data;
    }
  }
  //*************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::writehdf5(hid_t group_id, entitySet &eset) const 
  {
    typedef typename data_schema_traits<T> ::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    Loci::HDF5_WriteDomain(group_id, eset);

    hdf5write(group_id, traits_output_type, eset);

  }

  //*************************************************************************/
  template <class T> 
  void dmultiStoreRepI<T> :: hdf5write( hid_t group_id, IDENTITY_CONVERTER c, 
                                        const entitySet &eset)  const
  {
    int      rank = 1;
    hsize_t  dimension;
    hid_t    vDataspace, vDataset, vDatatype;

    std::vector<T>   newvec;
    entitySet :: const_iterator  ei;
    hash_map<int,std::vector<T> >:: const_iterator ci;
    std::vector<int> container(eset.size());
    size_t  arraySize= 0;
    int     count;

    size_t indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec      = ci->second;
        arraySize  += newvec.size();
        container[indx++] = newvec.size();
      }
    }

    dimension =  eset.size();
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5Tcopy( H5T_NATIVE_INT);
    vDataset   = H5Dcreate(group_id, "ContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    //------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //------------------------------------------------------------------------
    std::vector<T>  data(arraySize);

    indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newvec  = ci->second;
        for( int j = 0; j < newvec.size(); j++) 
          data[indx++] = newvec[j];
      }
    }

    typedef data_schema_traits<T> traits_type;
    DatatypeP dtype = traits_type::get_type();
    vDatatype = dtype->get_hdf5_type();

    rank      = 1;
    dimension = arraySize;

    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
                           vDataspace, H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

  }

  //***************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T> :: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                        const entitySet &eset)  const
  {   
    int      rank = 1, stateSize;
    hsize_t  dimension;
    hid_t    vDataspace, vDataset, vDatatype;
    typedef data_schema_traits<T> schema_traits;

    std::vector<T>   newVec;
    entitySet :: const_iterator  ei;
    hash_map<int,std::vector<T> >:: const_iterator ci;

    //------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //-------------------------------------------------------------------------
    std::vector<int> container(eset.size());
    std::vector<int> bucketSize;
    
    size_t  arraySize= 0;
    int     count;

    size_t indx = 0;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newVec  = ci->second;
        container[indx++] = newVec.size();
        for( int ivec = 0; ivec < newVec.size(); ivec++) {
          typename schema_traits::Converter_Type cvtr(newVec[ivec] );
          bucketSize.push_back(cvtr.getSize());
          arraySize  += cvtr.getSize();
        }
      }
    }

    dimension  = eset.size();
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5T_NATIVE_INT;
    vDataset   = H5Dcreate(group_id, "ContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
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
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci == attrib_data.end() ) continue;
      newVec =  ci->second;
      for( int ivec = 0; ivec < newVec.size(); ivec++) {
        typename schema_traits::Converter_Type cvtr( newVec[ivec] );
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

    //-----------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------
    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );
  }

  //************************************************************************/
}

#endif

#ifndef PARAMETER_H
#define PARAMETER_H

#include <mpi.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>
#include <data_traits.h>

namespace Loci {
  
  template<class T> class paramRepI : public storeRep {
    entitySet store_domain ;
    T attrib_data ;

    void hdf5read(hid_t group, IDENTITY_CONVERTER g );
    void hdf5read(hid_t group, USER_DEFINED_CONVERTER g );

    void hdf5write( hid_t group, IDENTITY_CONVERTER g,    const entitySet &en) const;
    void hdf5write( hid_t group, USER_DEFINED_CONVERTER g,const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size );
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size );

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size);
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size);

  public:
    paramRepI() { store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ; }
    paramRepI(const entitySet &p) { store_domain = p ;}
    virtual void allocate(const entitySet &p)  ;
    virtual ~paramRepI() ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq)  ;
    
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( hid_t group, entitySet &en) ;
    virtual void writehdf5( hid_t group,entitySet& en) const ;
    T *get_param() { return &attrib_data ; }
  } ;

  //**************************************************************************/

  template<class T> void paramRepI<T>::allocate(const entitySet &p) {
    store_domain = p ;
    dispatch_notify();
  }

  //**************************************************************************/

  template<class T> paramRepI<T>::~paramRepI<T>() {}

  //**************************************************************************/
    
  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p) const 
  {
    return new paramRepI<T>(p) ;
  }

  //**************************************************************************/

  template<class T> 
  store_type paramRepI<T>::RepType() const 
  {
    return PARAMETER ;
  }

  //**************************************************************************/

  template<class T> entitySet paramRepI<T>::domain() const {
    return store_domain ;
  }

  //**************************************************************************/
        
  template<class T> 
  std::ostream &paramRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << std::endl ;
    Loci::streamoutput(&attrib_data,1,s) ;
    s << '}' << std::endl ;
    return s ;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &paramRepI<T>::Input(std::istream &s) 
  {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      s.putback(ch) ;
      e = ~EMPTY ;
      allocate(e) ;
      attrib_data = T() ;
      Loci::streaminput(&attrib_data,1,s) ;
      return s ;
    }
        
    s >> e ;
    allocate(e) ;
        
    attrib_data = T() ;
    Loci::streaminput(&attrib_data,1,s) ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading parameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  //**************************************************************************/
  template<class T> 
  void paramRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    entitySet eset, ecommon;
    Loci::HDF5_ReadDomain(group_id, eset);

    ecommon = eset & user_eset;
    allocate( ecommon );
    hdf5read(group_id, schema_converter() );
  }

  //**************************************************************************/

  template<class T> 
  void paramRepI<T>::writehdf5( hid_t group_id, entitySet &eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    Loci::HDF5_WriteDomain(group_id, eset);

    hdf5write(group_id, schema_converter(), eset);

  }

  //**************************************************************************/

  template<class T> class param : public store_instance {
    typedef paramRepI<T> paramType ;
    T * data ;
  public:
    typedef T containerType ;
    param() { setRep(new paramType) ; }
    param(param &var) { setRep(var.Rep()) ; }
    param(storeRepP rp) { setRep(rp); }

    virtual ~param() ;

    param & operator=(param &p) {setRep(p.Rep()) ; return *this ; }

    param & operator=(storeRepP p) {setRep(p) ; return *this ; }
    param & operator=(const T &v) { *data = v ; return *this ; }

    virtual void notification() ;
    
    T * operator->() { return data ; }
    const T * operator->() const { return data ; }
    
    T * operator&() { return data ; }
    const T * operator &() const { return data ; }

    T &operator*() { return *data ; }
    const T &operator*() const { return *data ; }

    T &operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    void set_entitySet(const entitySet &ptn) {Rep()->allocate(ptn); }

    entitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  //**************************************************************************/

  template<class T> param<T>::~param() {}
    
  //**************************************************************************/

  template<class T> 
  void param<T>::notification()
  {  
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const param<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, param<T> &t)
  { return t.Input(s) ; }

  //**************************************************************************/

  template<class T> 
  class const_param : public store_instance {
    typedef T containerType ;
    typedef paramRepI<T> paramType ;
    const T * data ;
  public:
    const_param() { setRep(new paramType) ; }
    const_param(const_param<T> &var) { setRep(var.Rep()) ; }
    const_param(param<T> &var) { setRep(var.Rep()) ; }
    const_param(storeRepP rp) { setRep(rp); }
    
    virtual ~const_param() ;

    const_param & operator=(const_param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(storeRepP p)
    { setRep(p) ; return *this ;}

    virtual void notification() ;
    virtual instance_type access() const ;
        
    const T * operator->() const { return data ; }
    
    const T * operator &() const { return data ; }

    const T &operator*() const { return *data ; }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    entitySet domain() const { return Rep()->domain(); }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  //**************************************************************************/

  template<class T> const_param<T>::~const_param() {}

  //**************************************************************************/

  template<class T> 
  void const_param<T>::notification() 
  {  
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }
    
  //**************************************************************************/

  template<class T> 
  storeRepP paramRepI<T>::remap(const Map &m) const 
  {
    param<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    return r.Rep() ;
  }

  //**************************************************************************/

  template<class T> 
  void paramRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    param<T> p(st) ;
    attrib_data = *p ;
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> 
  void paramRepI<T>::gather(const Map &m, storeRepP &st,
                            const entitySet &context) 
  {
    param<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T> 
  void paramRepI<T>::scatter(const Map &m, storeRepP &st,
                             const entitySet &context) 
  {

    fatal((context - store_domain) != EMPTY) ;
    store_domain = m.image(context) ;

    param<T> p(st) ;
    attrib_data = *p ;
  }

  //**************************************************************************/
 
  template <class T> 
  int paramRepI<T>::pack_size( const entitySet &eset) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    return get_mpi_size( schema_converter(), eset );
  }

  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    return( sizeof(T) ) ;
  }

  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits;

    typename schema_traits::Converter_Type cvtr(attrib_data);
    int arraySize = cvtr.getSize() ;

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) );
  }
  //**************************************************************************/

  template <class T> 
  void paramRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &e ) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    packdata( schema_converter(), ptr, loc, size);
  }

  //**************************************************************************/
  template <class T>
  void paramRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                               int outcount )
  {
    MPI_Pack( &attrib_data, sizeof(T), MPI_BYTE, outbuf, outcount, &position, 
              MPI_COMM_WORLD) ;
  }

  //**************************************************************************/

  template <class T> 
  void paramRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                               int &position, int outcount )
  {
    int stateSize;
    
    typedef data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    typename data_schema_traits<T>::Converter_Type cvtr( attrib_data );

    std::vector<dtype> inbuf(cvtr.getSize());
    cvtr.getState( &inbuf[0], stateSize);

    MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
             MPI_COMM_WORLD);
    int incount =  stateSize*typesize;
    MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
             MPI_COMM_WORLD) ;

  }
  //**************************************************************************/

  template <class T> 
  void paramRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata( schema_converter(), ptr, loc, size);
  }  


  //**************************************************************************/
  template <class T> 
  void paramRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position, 
                                 int &insize)
  {

    /*
      typedef data_schema_traits<T> traits_type;
      DatatypeP    atom_type = traits_type::get_type();
      MPI_Datatype datatype  = atom_type->get_mpi_type();
    */
    MPI_Unpack( inbuf, insize, &position, &attrib_data, sizeof(T), 
                MPI_BYTE, MPI_COMM_WORLD) ;

  }

  //***********************************************************************/
  template <class T> 
  void paramRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                 int &position, int &insize)
  {

    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    MPI_Unpack( inbuf, 1, &position, &stateSize, 1, 
                MPI_INT, MPI_COMM_WORLD) ;

    std::vector<dtype> outbuf(stateSize);

    outcount = stateSize*sizeof(dtype);
    MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount, 
                MPI_BYTE, MPI_COMM_WORLD) ;

    typename schema_traits::Converter_Type  cvtr( attrib_data );
    cvtr.setState( &outbuf[0], stateSize);

  }

  //***********************************************************************/
  
  template<class T> store_instance::instance_type
  const_param<T>::access() const
  { return READ_ONLY; }

  //**************************************************************************/
    
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_param<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/
  template <class T> 
  void paramRepI<T> :: hdf5write( hid_t group_id, IDENTITY_CONVERTER g,
                                  const entitySet &eset ) const
  {
    int arraySize =  eset.size(); 

    typedef data_schema_traits<T> traits_type;

    DatatypeP  dtype = traits_type::get_type();
    hid_t vDatatype = dtype->get_hdf5_type();

    int      rank = 1;
    hsize_t  dimension = 1;

    entitySet :: const_iterator ci;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);

    hid_t cparms   = H5Pcreate (H5P_DATASET_CREATE);
    hid_t vDataset = H5Dcreate(group_id, "VariableData", vDatatype,
                               vDataspace, cparms);
    T data;
    data = attrib_data;

    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

  }

  //*********************************************************************/

  template <class T> 
  void paramRepI<T> :: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, 
                                  const entitySet &eset) const
  {

    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;

    dtype *data ;

    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------
    T Obj;
    Obj = attrib_data;

    typename schema_traits::Converter_Type cvtr(Obj);
    int stateSize  = cvtr.getSize();

    data =  new dtype[stateSize];

    cvtr.getState( data, stateSize);

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------

    typedef data_schema_traits<dtype> traits_type;

    DatatypeP atom_type = traits_type::get_type() ;
    hid_t vDatatype = atom_type->get_hdf5_type();

    int rank = 1;
    hsize_t dimension = stateSize;
    hid_t vDataspace  = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDataset    = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace,
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
  template <class T> 
  void paramRepI<T> :: hdf5read(hid_t group_id, IDENTITY_CONVERTER g )
  { 

    hsize_t  dimension;

    typedef data_schema_traits<T> traits_type;
    DatatypeP  dtype = traits_type::get_type();

    hid_t vDatatype  = dtype->get_hdf5_type();
    hid_t vDataset   = H5Dopen(group_id,"VariableData");
    hid_t vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    H5Dread(vDataset, vDatatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &attrib_data);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

  }

  //*************************************************************************/

  template <class T> 
  void paramRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER g )
  { 

    hsize_t   dimension;

    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;

    dtype *data;

    DatatypeP atom_type = data_schema_traits<dtype>::get_type() ;
    
    hid_t vDatatype = atom_type->get_hdf5_type();

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    hid_t vDataset   = H5Dopen(group_id,"VariableData");
    hid_t vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    size_t  stateSize;
    stateSize =  dimension;
    data      = new dtype[stateSize];

    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, data);

    entitySet::const_iterator ci;
    typename data_schema_traits<T>::Converter_Type cvtr( attrib_data);
    cvtr.setState( data, stateSize );

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );

    delete [] data;

  }

  //***************************************************************************

}

#endif
    

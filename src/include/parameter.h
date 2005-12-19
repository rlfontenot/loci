#ifndef PARAMETER_H
#define PARAMETER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>


#include <mpi.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <distribute.h>
#include <data_traits.h>

namespace Loci {

  template<class T> class paramRepI : public storeRep {
    entitySet store_domain ;
    T attrib_data ;

    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, frame_info &fi, const entitySet &en);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const entitySet &en);

    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g,const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size );
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size );

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size);
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size);
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) ;
    frame_info write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) ;
  public:
    paramRepI() { store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ; }
    paramRepI(const entitySet &p) { store_domain = p ;}
    virtual void allocate(const entitySet &p)  ;
    virtual void shift(int_type offset) ;
    virtual ~paramRepI() ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
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
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq)  ;

    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    T *get_param() { return &attrib_data ; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ;

  //**************************************************************************/

  template<class T> void paramRepI<T>::allocate(const entitySet &p) {
    store_domain = p ;
    dispatch_notify();
  }

  template<class T> void paramRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> paramRepI<T>::~paramRepI() {}

  //**************************************************************************/

  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p) const
  {
    return new paramRepI<T>(p) ;
  }

  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a parameter " << endl ;
    return sp ;
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
    entitySet dom = domain() ;
    if(dom == ~EMPTY) {
      Loci::streamoutput(&attrib_data,1,s) ;
    } else {
      s << '{' << domain() << std::endl ;
      Loci::streamoutput(&attrib_data,1,s) ;
      s << '}' << std::endl ;
    }
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
  template<class T>
  frame_info paramRepI<T>::read_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return read_frame_info(group_id, schema_converter()) ;
  }

  template<class T>
  frame_info paramRepI<T>::read_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
    int is_stat = 0 ;
    int sz = 0 ;
    hid_t datatype = H5T_NATIVE_INT ;
    hid_t dataset = H5Dopen(group_id, "is_stat") ;
    H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
    H5Dclose(dataset) ;
    dataset = H5Dopen(group_id, "vec_size") ;
    H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
    H5Dclose(dataset) ;
    return frame_info(is_stat, sz);
  }
  template<class T>
  frame_info paramRepI<T>::read_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    hsize_t dimension = 0 ;
    hid_t dataspace ;
    hid_t datatype = H5T_NATIVE_INT ;
    hid_t dataset ;
    int is_stat = 0 ;
    int sz = 0 ;
    frame_info fi ;
    dataset = H5Dopen(group_id, "is_stat") ;
    H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
    H5Dclose(dataset) ;
    dataset = H5Dopen(group_id, "vec_size") ;
    H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
    H5Dclose(dataset) ;
    dataset = H5Dopen(group_id, "second_level") ;
    dataspace = H5Dget_space(dataset) ;
    H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    fi.is_stat = is_stat ;
    fi.size = sz ;
    std::vector<int> vint ;
    int tmp ;
    H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT, &tmp) ;
    vint.push_back(tmp) ;
    fi.second_level = vint ;
    H5Dclose(dataset) ;
    H5Sclose(dataspace) ;
    return fi ;
  }
  template<class T>
  frame_info paramRepI<T>::write_frame_info(hid_t group_id) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return write_frame_info(group_id, schema_converter()) ;
  }

  template<class T>
  frame_info paramRepI<T>::write_frame_info(hid_t group_id, IDENTITY_CONVERTER g) {
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
  frame_info paramRepI<T>::write_frame_info(hid_t group_id, USER_DEFINED_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    typename schema_traits::Converter_Type cvtr(attrib_data);
    stateSize = cvtr.getSize();
    fi.second_level.push_back(stateSize) ;

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
    rank = 1 ;
    dimension = 1 ;
    dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    hid_t dataset = H5Dcreate(group_id,"second_level",H5T_NATIVE_INT,
                              dataspace,H5P_DEFAULT) ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &fi.second_level[0]) ;
    H5Sclose(memspace) ;
    H5Dclose(dataset) ;
    H5Sclose(dataspace) ;
    
    return fi ;
  }


  template<class T>
  DatatypeP paramRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T>
  DatatypeP paramRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T>
  DatatypeP paramRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  template<class T>
  void paramRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }

  //**************************************************************************/

  template<class T>
  void paramRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }

  //**************************************************************************/

  template<class T> class param : public store_instance {
    typedef paramRepI<T> paramType ;
    T * data ;
  public:
    typedef T containerType ;
    param() { setRep(new paramType) ; }
    param(const param &var) {
      *data = *(var.data) ;
      set_entitySet(var.domain());
    }
    param(storeRepP rp) { setRep(rp); }

    virtual ~param() ;

    param & operator=(const param &p) {
      *data = *(p.data) ;
      set_entitySet(p.domain());
      return *this ;
    }
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
    const T * restrict data ;
    const_param(const const_param<T> &var) {setRep(var.Rep()) ;}
    const_param(param<T> &var) { setRep(var.Rep()) ; }
  public:
    const_param() { setRep(new paramType) ; }
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

    const T * restrict operator->() const { return data ; }

    const T * restrict operator &() const { return data ; }

    const T & restrict operator*() const { return *data ; }

    const T & restrict operator[](int indx) const {
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
  storeRepP paramRepI<T>::remap(const dMap &m) const
  {
    param<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    return r.Rep() ;
  }

  template<class T>
  storeRepP paramRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP paramRepI<T>::thaw() {
    return getRep() ;
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
  void paramRepI<T>::gather(const dMap &m, storeRepP &st,
                            const entitySet &context)
  {
    param<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T>
  void paramRepI<T>::scatter(const dMap &m, storeRepP &st,
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

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) + sizeof(int));
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

    MPI_Unpack( inbuf, insize, &position, &stateSize, 1,
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
  void paramRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;

      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &attrib_data) ;
      H5Sclose(memspace) ;
    }
  }

  //*********************************************************************/

  template <class T>
  void paramRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
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
      int stateSize = 0 ;
      T tmp = attrib_data;
      typename schema_traits::Converter_Type cvtr(tmp);
      cvtr.getState(tmp_array, stateSize) ;
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  }

  //**************************************************************************/
  template <class T>
  void paramRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, frame_info &fi, const entitySet &en)
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
      }
      attrib_data = tmp_array[0] ;

      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  }

  //*************************************************************************/

  template <class T>
  void paramRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const entitySet &en)
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
      typename data_schema_traits<T>::Converter_Type cvtr(attrib_data);
      int bucsize = vint[0] ;
      cvtr.setState(tmp_array, bucsize) ;
      H5Sclose(memspace) ;
      delete [] tmp_array ;
    }
  }

  //***************************************************************************

}

#endif

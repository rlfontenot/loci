#ifndef PARAMETER_H
#define PARAMETER_H

#include <mpi.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

#include <hdf5CC/H5cpp.h>
#include <hdf5_traits.h>
#include <hdf5_memento.h>

#include <hdf5_readwrite.h>


namespace Loci {
  
  template<class T> class paramRepI : public storeRep {
    entitySet store_domain ;
    T attrib_data ;

    void hdf5read(H5::Group group, DEFAULT_CONVERTER  g );
    void hdf5read(H5::Group group, IDENTITY_CONVERTER g );
    void hdf5read(H5::Group group, USER_DEFINED_CONVERTER g );

    void hdf5write( H5::Group group, DEFAULT_CONVERTER g,     const entitySet &en) const;
    void hdf5write( H5::Group group, IDENTITY_CONVERTER g,    const entitySet &en) const;
    void hdf5write( H5::Group group, USER_DEFINED_CONVERTER g,const entitySet &en) const;

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
    virtual void readhdf5( H5::Group group, entitySet &en) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    T *get_param() { return &attrib_data ; }
  } ;

  //****************************************************************************

  template<class T> void paramRepI<T>::allocate(const entitySet &p) {
    store_domain = p ;
    dispatch_notify();
  }

  //****************************************************************************

  template<class T> paramRepI<T>::~paramRepI<T>() {}

  //****************************************************************************
    
  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p) const 
  {
    return new paramRepI<T>(p) ;
  }

  //****************************************************************************

  template<class T> 
  store_type paramRepI<T>::RepType() const 
  {
    return PARAMETER ;
  }

  //****************************************************************************

  template<class T> entitySet paramRepI<T>::domain() const {
    return store_domain ;
  }

  //****************************************************************************
        
  template<class T> 
  std::ostream &paramRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << std::endl ;
    s << attrib_data << std::endl ;
    s << '}' << std::endl ;

    return s ;
  }

  //****************************************************************************

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
      s>>attrib_data ;
      return s ;
    }
        
    s >> e ;
    allocate(e) ;
        
    s >> attrib_data ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading parameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  //****************************************************************************
  template<class T> 
  void paramRepI<T>::readhdf5( H5::Group group, entitySet &user_eset)
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    entitySet eset, ecommon;
    HDF5_ReadDomain(group, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );

    hdf5read(group, traits_output_type );

  }

  //****************************************************************************

  template<class T> 
  void paramRepI<T>::writehdf5( H5::Group group,entitySet &en) const
  {

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    hdf5write(group, traits_output_type, en);

  }

  //****************************************************************************

  template<class T> class param : public store_instance {
    typedef paramRepI<T> paramType ;
    T * data ;
  public:
    typedef T containerType ;
    param() { setRep(new paramType) ; }
    param(param &var) { setRep(var.Rep()) ; }
    param(storeRepP &rp) { setRep(rp); }

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

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  //****************************************************************************

  template<class T> param<T>::~param() {}
    
  //****************************************************************************

  template<class T> 
  void param<T>::notification()
  {  
     NPTR<paramType> p(Rep());
     if(p!=0) data = p->get_param() ;
     warn(p==0);
  }

  //****************************************************************************

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const param<T> &t)
  { return t.Print(s) ; }

  //****************************************************************************

  template<class T> 
  inline std::istream & operator>>(std::istream &s, param<T> &t)
  { return t.Input(s) ; }

  //****************************************************************************

  template<class T> 
  class const_param : public store_instance {
    typedef T containerType ;
    typedef paramRepI<T> paramType ;
    const T * data ;
  public:
    const_param() { setRep(new paramType) ; }
    const_param(const_param<T> &var) { setRep(var.Rep()) ; }
    const_param(param<T> &var) { setRep(var.Rep()) ; }
    const_param(storeRepP &rp) { setRep(rp); }
    
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
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  //****************************************************************************

  template<class T> const_param<T>::~const_param() {}

  //****************************************************************************

  template<class T> 
  void const_param<T>::notification() 
  {  
     NPTR<paramType> p(Rep());
     if(p!=0) data = p->get_param() ;
     warn(p==0);
  }
    
  //****************************************************************************

  template<class T> 
  storeRepP paramRepI<T>::remap(const Map &m) const 
  {
    param<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    return r.Rep() ;
  }

  //****************************************************************************

  template<class T> 
  void paramRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    param<T> p(st) ;
    attrib_data = *p ;
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
    dispatch_notify() ;
  }

  //****************************************************************************

  template<class T> 
  void paramRepI<T>::gather(const Map &m, storeRepP &st,
                            const entitySet &context) 
  {
    param<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = context ;
  }

  //****************************************************************************

  template<class T> 
  void paramRepI<T>::scatter(const Map &m, storeRepP &st,
                             const entitySet &context) 
  {

    param<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = m.image(context) ;
  }

  //****************************************************************************
 
  template <class T> 
  int paramRepI<T>::pack_size( const entitySet &e) 
  {
    int size ;
    size = sizeof(T);
    return(size) ;
  }

  //****************************************************************************

  template <class T> 
  void paramRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &e ) 
  {
    MPI_Pack(&attrib_data, sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
  }

  //****************************************************************************

  template <class T> 
  void paramRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)  {
    MPI_Unpack(ptr, size, &loc, &attrib_data, sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
  }  

  //****************************************************************************
  
  template<class T> store_instance::instance_type
    const_param<T>::access() const
  { return READ_ONLY; }

  //****************************************************************************
    
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_param<T> &t)
  { return t.Print(s) ; }

  //****************************************************************************

  template <class T> 
  void paramRepI<T> :: hdf5write( H5::Group group, DEFAULT_CONVERTER g,
                                  const entitySet &eset ) const
  {
    //write out the domain   
    HDF5_WriteDomain(group, eset);

    typedef hdf5_schema_traits<T> traits_type;
    std::ostringstream oss;
    oss<< attrib_data;
    std::string memento = oss.str();
    hsize_t size = memento.length();

    int rank = 1;
    hsize_t dimension[1];
    dimension[0] =  size+1;

    int num_intervals = eset.num_intervals();
    interval *it = new interval[num_intervals];


    try{
      H5::DataType  datatype = traits_type::get_type();
      H5::DataSpace dataspace( rank, dimension );
      H5::DataSet   dataset = group.createDataSet( "param", datatype, dataspace);
      dataset.write( memento.c_str(), datatype );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

    delete [] it;
  }


  //****************************************************************************
  
  template <class T> 
  void paramRepI<T> :: hdf5write( H5::Group group, IDENTITY_CONVERTER g,
                                  const entitySet &eset ) const
  {

    //write out the domain   
    HDF5_WriteDomain(group, eset);

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------
    typedef hdf5_schema_traits<T> traits_type;

    int rank = 1;
    hsize_t  dimension[1];

    dimension[0] =  1;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = traits_type::get_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", vDatatype, vDataspace);
      vDataset.write( &attrib_data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

  };

  //***********************************************************************

  template <class T> 
  void paramRepI<T> :: hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, 
                                  const entitySet &eset) const
  {

    //write out the domain   
    HDF5_WriteDomain(group, eset);

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------
    typedef hdf5_schema_converter_traits<T> converter_traits; 
    converter_traits::memento_type *data, *buf;

    Memento<T> memento( &attrib_data );
    int arraySize = memento.getSize();
    
 	 data =  new typename converter_traits::memento_type[arraySize];
 	 buf  =  new typename converter_traits::memento_type[arraySize];

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    size_t indx = 0;
    int    stateSize;
    memento.getState(data, stateSize);

    for( int i = 0; i < stateSize; i++)
         cout << data[i] << endl;

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------
    int rank = 1;
    hsize_t  dimension[1];

    dimension[0] =  stateSize;
    cout << stateSize << endl;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", vDatatype, vDataspace);

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

  //***********************************************************************

  template <class T> 
  void paramRepI<T> :: hdf5read(H5::Group group,DEFAULT_CONVERTER g )
  {

    entitySet    eset;

    typedef hdf5_schema_traits <T> traits_type;
    entitySet num;	

    //get domain data
    try{
      H5::DataSet dataset_domain      = group.openDataSet( "domain");
      H5::DataSpace dataspace_domain  = dataset_domain.getSpace();

      hsize_t dims_domain[1];
      dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
      int *data_domain = new int[dims_domain[0]];
      dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );

      for(int i=0;i<dims_domain[0];i++){
        num |=interval(data_domain[i],data_domain[i+1]);
        i++;
      }

      //get param data
      H5::DataType datatype         =  traits_type::get_type();
      H5::DataSet dataset_param     =  group.openDataSet( "param");
      H5::DataSpace dataspace_param =  dataset_param.getSpace();

      hsize_t dims_param[1];
      dataspace_param.getSimpleExtentDims( dims_param, NULL);

      char* data_param = new char[dims_param[0]];
      dataset_param.read( data_param, datatype );
      std::istringstream iss(data_param);

      iss >> attrib_data;
      delete [] data_param;
      delete [] data_domain;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

  };

  //***************************************************************************

  template <class T> 
  void paramRepI<T> :: hdf5read(H5::Group group,IDENTITY_CONVERTER g )
  { 

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------

    typedef hdf5_schema_traits <T> traits_type;
    H5::DataType  vDatatype = traits_type::get_type();

    try {
	   H5::DataSet   vDataset   = group.openDataSet( "VariableData");

      vDataset.read( &attrib_data, vDatatype );
    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }
  };

  //***************************************************************************

  template <class T> 
  void paramRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER g )
  { 

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------
   hsize_t   dimension[1];

   typedef hdf5_schema_converter_traits<T> converter_traits; 
   converter_traits::memento_type *data, *buf;

   H5::DataType  vDatatype  = converter_traits::get_variable_HDF5_type();
	H5::DataSet   vdataset   = group.openDataSet( "VariableData");
	H5::DataSpace vdataspace = vdataset.getSpace();

	vdataspace.getSimpleExtentDims( dimension, NULL);

   size_t  stateSize;
   stateSize =  dimension[0];
   data      = new typename converter_traits::memento_type[stateSize];
   buf       = new typename converter_traits::memento_type[stateSize];

	vdataset.read( data, vDatatype);

   entitySet::const_iterator ci;

   Memento<T> memento( &attrib_data );
   for( int i = 0; i < stateSize; i++) 
        buf[i] = data[i];
   memento.setState( buf, stateSize);

   delete [] data;
   delete [] buf;

  };

  //***************************************************************************

}

#endif
    

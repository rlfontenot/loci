#ifndef DSTOREVEC_H
#define DSTOREVEC_H 

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>
#include <hdf5CC/H5cpp.h>

#include <Map.h>
#include <multiMap.h>

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif


namespace Loci {

  using std::hash_map ;

  template<class T> class dstoreVecRepI : public storeRep {

  lmutex                    mutex ;
  entitySet                 store_domain ;
  int                       size;
  hash_map<int,std::vector<T> >  attrib_data;
  
  void hdf5read( H5::Group group, DEFAULT_CONVERTER      c, entitySet &en, entitySet &usr );
  void hdf5read( H5::Group group, IDENTITY_CONVERTER     c, entitySet &en, entitySet &usr );
  void hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, entitySet &en, entitySet &usr );

  void hdf5write( H5::Group group, DEFAULT_CONVERTER  c, const entitySet &en ) const;
  void hdf5write( H5::Group group, IDENTITY_CONVERTER c, const entitySet &en ) const;
  void hdf5write( H5::Group group, USER_DEFINED_CONVERTER c, const entitySet &en ) const;


  public:
    dstoreVecRepI() {}

    dstoreVecRepI(const entitySet &p) { allocate(p) ; }

    virtual ~dstoreVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e ) ;
    virtual void pack(void * ptr, int &loc, int &size, const entitySet &e ) ;
    virtual void unpack(void * ptr, int &loc, int &size, const sequence &seq) ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group, entitySet &user_eset) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;

    virtual void set_elem_size(int sz) ;

    hash_map<int,std::vector<T> > *get_attrib_data() { return &attrib_data; }
    int get_size() const { return size; }
  } ;

  //****************************************************************************

   template<class T> 
	std::ostream &dstoreVecRepI<T>::Print(std::ostream &s) const
   {

      s << '{' << domain() << std::endl ;

      s << size << std::endl ;

      hash_map<int, std::vector<T> >  :: const_iterator   ci;
      std::vector<T>   newVec;
    
      FORALL(domain(),ii) {
         ci =  attrib_data.find(ii);
         if( ci != attrib_data.end() ) {
             newVec = ci->second;
             for( int i = 0; i < newVec.size(); i++)
                  s << newVec[i] <<  "   ";
             s << std::endl;
         }
      }ENDFORALL ;

      s << '}' << std::endl ;

      return s ;
   }

  //***************************************************************************

  template<class T> 
  std::istream &dstoreVecRepI<T>::Input(std::istream &s)
  {
      char ch ;

      do ch = s.get(); while(ch==' ' || ch=='\n') ;

      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
        return s ;
      }
      entitySet e ;
      int sz ;

      s >> e ;
      s >> size ;

      std::vector<T>  newVec(size);
      allocate(e) ;

      FORALL(e,ii) {
        for( int i = 0; i < size; i++)
             s  >> newVec[i];
        attrib_data[ii] = newVec;
      } ENDFORALL ;
    
      do ch = s.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
      }

      return s ;
    }


  //****************************************************************************

template<class T> 
void dstoreVecRepI<T>::readhdf5( H5::Group group, entitySet &user_eset)
{
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    //--------------------------------------------------------------------------
    // Read the vector size ...
    //--------------------------------------------------------------------------
    hsize_t dimension[1];

    dimension[0] = 1;

    H5::DataType  datatype  = H5::PredType::NATIVE_INT;
    H5::DataSet   dataset   = group.openDataSet( "VecSize");
    H5::DataSpace dataspace = dataset.getSpace();

    dataspace.getSimpleExtentDims( dimension, NULL);

    dataset.read( &size, H5::PredType::NATIVE_INT );

    //--------------------------------------------------------------------------
    entitySet  eset, ecommon;

    HDF5_ReadDomain(group, eset);

    ecommon = eset & user_eset;

    allocate( ecommon );
    hdf5read( group, traits_type, eset, ecommon);
}

//******************************************************************************

template<class T> 
void dstoreVecRepI<T>::writehdf5( H5::Group group,entitySet& en) const
{

    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;

    hdf5write(group, traits_output_type, en );
}

//******************************************************************************

template<class T> 
void dstoreVecRepI<T>::allocate(const entitySet &ptn) 
{

	 entitySet :: const_iterator  ci;

    std::vector<T>   newVec;

    if( size > 1) newVec.reserve( size );

	 for( ci = ptn.begin(); ci != ptn.end(); ++ci)
	     attrib_data[*ci] = newVec;

    store_domain = ptn ;
    dispatch_notify() ;

  }

  //****************************************************************************

  template<class T> 
  dstoreVecRepI<T>::~dstoreVecRepI<T>() 
  {
    attrib_data.clear();
  }

  //****************************************************************************

  template<class T>
  storeRep *dstoreVecRepI<T>::new_store(const entitySet &p) const 
  {
       return new dstoreVecRepI<T>(p)  ;
  }

  //****************************************************************************

  template<class T> 
  store_type dstoreVecRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //****************************************************************************

  template<class T> 
  entitySet dstoreVecRepI<T>::domain() const 
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

  //***************************************************************************

  template<class T> 
  void dstoreVecRepI<T>::set_elem_size(int sz) 
  {

    mutex.lock() ;

    if(size != sz) {
      if(size != 0) warn(size != sz) ;
    
      size = sz ;
      fatal(sz < 1) ;
      allocate(store_domain) ;
    }

    mutex.unlock() ;

  }

  //****************************************************************************
      
  template<class T> class dstoreVec : public store_instance {
    typedef dstoreVecRepI<T>    storeType ;
	 hash_map<int, std::vector<T> >  *attrib_data;
    int                         size;
  public:
    typedef std::vector<T> containerType ;
    dstoreVec() {setRep(new storeType) ;}
    dstoreVec(dstoreVec<T> &var) {setRep(var.Rep()) ;}
    dstoreVec(storeRepP &rp) { setRep(rp) ;}

    virtual ~dstoreVec() ;

    virtual void notification() ;

    dstoreVec<T> & operator=(dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    dstoreVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      Rep()->set_elem_size(size) ;
    }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> &elem(int indx) {
	   return( (*attrib_data)[indx] );
    }

    std::vector<T> &operator[](int indx) {
      return( elem(indx) );
	 }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //****************************************************************************

  template<class T> 
  dstoreVec<T>::~dstoreVec<T>() { }

  //****************************************************************************

  template<class T> 
  void dstoreVec<T>::notification() 
  {
      NPTR<storeType> p(Rep()) ;
      if(p!=0) 
         attrib_data = p->get_attrib_data() ;

      warn(p == 0) ;
  }

  //****************************************************************************

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const dstoreVec<T> &t)
    { return t.Print(s) ; }

  //****************************************************************************

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dstoreVec<T> &t)
    { return t.Input(s) ; }

  //****************************************************************************

  template<class T> class const_dstoreVec : public store_instance {
    typedef dstoreVecRepI<T> storeType ;
	 hash_map<int, std::vector<T> >  *attrib_data;
    int size ;
  public:
    typedef std::vector<T> containerType ;
    const_dstoreVec() { setRep(new storeType) ; }
    const_dstoreVec(const_dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec(dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dstoreVec() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_dstoreVec<T> & operator=(dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_dstoreVec<T> & operator=(const_dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_dstoreVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }

    entitySet domain() const { return Rep()->domain() ; }

    std::vector<T> elem(int indx) const {
	     std::vector<T>   newVec;
	     hash_map<int, std::vector<T> > :: const_iterator   ci;
		  ci = attrib_data->find(indx);
		  if( ci != attrib_data->end()){
		      newVec = ci->second;
			   return( newVec );
	      } else {
		    cout << "Error: Out of range data " << endl;
          exit(0);
        }
		  return newVec;
    }
		       
    std::vector<T>  operator[](int indx) const {
	     return( elem(indx) );
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  } ;

  //***************************************************************************

  template<class T> 
  const_dstoreVec<T>::~const_dstoreVec<T>() { }

  //****************************************************************************

  template<class T> 
  void const_dstoreVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      attrib_data = p->get_attrib_data() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  //****************************************************************************

  template<class T> 
  store_instance::instance_type const_dstoreVec<T>::access() const
  { return READ_ONLY ; }

  //****************************************************************************

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_dstoreVec<T> &t)
  { return t.Print(s) ; }

  //****************************************************************************

  template <class T> 
  storeRepP dstoreVecRepI<T>::remap(const Map &m) const {
    dstoreVec<T> s ;

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;

    return s.Rep() ;
  }

  //***************************************************************************

  template <class T> 
  void dstoreVecRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
   /*
    const_storeVec<T> s(st) ;

    int sz = s.size() ;
    set_elem_size(sz) ;
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      T *p = base_ptr + i*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
   */

  }

  //***************************************************************************

  template <class T> 
  void dstoreVecRepI<T>::gather(const Map &m, storeRepP &st, const entitySet &context) 
  {
    const_dstoreVec<T> s(st) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]];
    } ENDFORALL ;

  }

  //****************************************************************************

  template <class T> 
  void dstoreVecRepI<T>::scatter(const Map &m, storeRepP &st, const entitySet &context)
  {
    const_dstoreVec<T> s(st) ;

    FORALL(context,i) {
      attrib_data[m[i]] = s[i];
    } ENDFORALL ;
  }

  //****************************************************************************

  template <class T> 
  int dstoreVecRepI<T>::pack_size( const entitySet &e) 
  {
    return (sizeof(T)*e.size()*size);
  }
  
  //****************************************************************************

  template <class T> 
  void dstoreVecRepI<T>::pack(void * ptr, int &loc, int &size, const entitySet &e ) 
  {
      
/*
     T     *buf = NULL;
     int   numBytes, numentity = e.size();

     numBytes =  size*numentity*sizeof(T);
     buf      =  (T *)malloc( numBytes);
     fatal( buf == NULL );
         
     entitySet :: const_iterator   ci;
     hash_map<int,std::vector<T> >  iter;

     int indx = 0;
     for( ci = e.begin(); ci != e.end(); ++ci){
         iter = attrib_data.find(ii);
         if( iter != attrib_data.end() ) {
             newVec = iter->second;
             for( j = 0; j < size; j++)
                  buf[indx++] =  newVec[j];
         }
     }

     MPI_Pack( buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD);

     free( buf );
*/
  }

  //***************************************************************************
  
  template <class T> 
  void dstoreVecRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {
/*
    Loci :: int_type   indx, jndx;
    int                numentity, numBytes;
    T                  *buf;

    for(int i = 0; i < seq.num_intervals(); ++i) {
        numentity =  abs(seq[i].second - seq[i].first) + 1; 
        numBytes  =  size*numentity*sizeof(T);
        buf       =  (T *) malloc( numBytes );
        MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

        jndx = 0;
        if(seq[i].first > seq[i].second) {
	        for(indx = seq[i].first; indx >= seq[i].second; --indx) {
               attrib_data[indx] =  buf[size*jndx];
               jndx++;
           }
        } else {
	        for(indx = seq[i].first; indx <= seq[i].second; ++indx){
               attrib_data[indx] =  buf[size*jndx];
               jndx++;
           }
        }
        free(buf);
    }
*/
  }

  //***************************************************************************
  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( H5::Group group, DEFAULT_CONVERTER c, 
                                     entitySet &eset, entitySet &user_eset )
  {
   cout << "Fatal: Default read converter not implemeneted for dynamic storeVec " << endl;
   exit(0);
  }

  //***************************************************************************

  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( H5::Group group, IDENTITY_CONVERTER c, 
                                     entitySet &eset, entitySet &user_eset )
  {

    hsize_t dimension[1];
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet::const_iterator ci;

   //---------------------------------------------------------------------------
   // Calculate the offset of each entity in file ....
   //---------------------------------------------------------------------------
   store<unsigned> offset;
   offset.allocate(eset);

   arraySize = 0;
   for( ci = eset.begin(); ci != eset.end(); ++ci) {
        offset[*ci] = arraySize;
        arraySize  += size;
   }

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------
   int num_intervals = user_eset.num_intervals();
   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

   T   *data;

   dimension[0] = arraySize;
   H5::DataSpace mDataspace(rank, dimension);   // memory  dataspace
   H5::DataSpace vDataspace(rank, dimension);

   typedef hdf5_schema_traits<T> traits_type;
   H5::DataType vDatatype = traits_type::get_type();
	H5::DataSet  vDataset   = group.openDataSet( "VariableData");

   hssize_t  start_mem[] = {0};  // determines the starting coordinates.
   hsize_t   stride[]    = {1};  // which elements are to be selected.
   hsize_t   block[]     = {1};  // size of element block;
   hssize_t  foffset[]   = {0};  // location (in file) where data is read.
   hsize_t   count[]     = {0};  // how many positions to select from the dataspace

   for( int k = 0; k < num_intervals; k++) {
        count[0] = 0;
        for( int i = it[k].first; i <= it[k].second; i++)
             count[0] +=  size;

        data = new T[count[0]];

        foffset[0] = offset[it[k].first];

        mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);
        vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        vDataset.read( data, vDatatype, mDataspace, vDataspace);

        indx = 0;
        for( int i = it[k].first; i <= it[k].second; i++) {
             attrib_data[i].clear();
             for( int m = 0; m < size; m++) 
                 attrib_data[i].push_back( data[indx++] );
        }

        delete[] data;
   }

  }

  //***************************************************************************

  template <class T> 
  void dstoreVecRepI<T> :: hdf5read( H5::Group group, USER_DEFINED_CONVERTER c, 
                                     entitySet &en, entitySet &user_eset )
  {

    hsize_t dimension[1];

    //--------------------------------------------------------------------------
    // Size of each container ....
    //--------------------------------------------------------------------------

    dimension[0] = 1;

    H5::DataType  fdatatype  = H5::PredType::NATIVE_INT;
	   H5::DataSet   fdataset   = group.openDataSet( "ContainerSize");
	   H5::DataSpace fdataspace = fdataset.getSpace();

	   fdataspace.getSimpleExtentDims( dimension, NULL);
	   int *vbucket = new int[dimension[0]];

	   fdataset.read( vbucket, H5::PredType::NATIVE_INT );
    int maxStateSize = *max_element( vbucket, vbucket + (int)dimension[0] );

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------

   typedef hdf5_schema_converter_traits<T> converter_traits; 
   typename converter_traits::memento_type *data, *buf;

   H5::DataType  vDatatype  = converter_traits::get_variable_HDF5_type();
   H5::DataSet   vdataset   = group.openDataSet( "variable");
   H5::DataSpace vdataspace = vdataset.getSpace();

   vdataspace.getSimpleExtentDims( dimension, NULL);

   size_t arraySize = dimension[0];
   data   = new typename converter_traits::memento_type[arraySize];

   vdataset.read( data, vDatatype);

   //---------------------------------------------------------------------------
   // Fill the objects ....
   //---------------------------------------------------------------------------
   
   buf    = new typename converter_traits::memento_type[maxStateSize];

   entitySet::const_iterator ci;
   hash_map<int, std::vector<T> > ::const_iterator iter;
   T       newObj;

   Memento<T> memento( newObj );
 
   size_t indx = 0;
   int    offset, bucketSize, bucketID = 0;
   for( ci = en.begin(); ci != en.end(); ++ci) {
        for( int ivec = 0; ivec < size; ivec++) {
             bucketSize = vbucket[bucketID];
             for( int i = 0; i < bucketSize; i++) 
                  buf[i] = data[indx++];
             newObj = memento.setState( buf, bucketSize);
             attrib_data[*ci].push_back(newObj);
             bucketID++;
        }
   }

  }

  //****************************************************************************
  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( H5::Group group, DEFAULT_CONVERTER c, 
                                     const entitySet &eset ) const
  {
    cout << "Fatal: Default converter not implemented for dynamic storevec " << endl;
    exit(0);
  }
  //****************************************************************************

  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( H5::Group group, IDENTITY_CONVERTER c, 
                                     const entitySet &eset ) const
  {

    hsize_t   dimension[1];
    int       rank = 1;

    //write out the domain   
    HDF5_WriteDomain(group, eset);

//-----------------------------------------------------------------------------
// write the Vector size
//-----------------------------------------------------------------------------
   dimension[0]=  1;
   try{
      H5::DataSpace sdataspace( rank, dimension );
      H5::DataSet   sdataset = group.createDataSet( "VecSize",
                                   H5::PredType::NATIVE_INT, sdataspace );
      sdataset.write( &size, H5::PredType::NATIVE_INT );
    }
    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    int arraySize  = size*eset.size();

    T *data =  new T[arraySize];

    entitySet::const_iterator ci;
    hash_map<int, std::vector<T> > ::const_iterator iter;
    std::vector<T>   newvec;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         iter = attrib_data.find(*ci);
         if( iter != attrib_data.end() ) {
             newvec = iter->second;
             for( int i = 0; i < size; i++)
                  data[indx++] =  newvec[i];
        }
    }

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------
    typedef hdf5_schema_traits<T> traits_type;

    dimension[0]=  arraySize;

    try {
      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = traits_type::get_type();
      H5::DataSet   vDataset  = group.createDataSet( "VariableData", 
                                                      vDatatype, vDataspace);

      vDataset.write( data, vDatatype );
    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    //-----------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------

    delete [] data;

  }

  //****************************************************************************

  template <class T>  
  void dstoreVecRepI<T> ::hdf5write( H5::Group group, USER_DEFINED_CONVERTER g, 
                                     const entitySet &eset ) const
  {
    hsize_t   dimension[1];
    int       rank = 1;

//-----------------------------------------------------------------------------
// Objective : Write store datatype into HDF5 Format which are user defined
//             datatypes or STL containers. Such datatypes are first written
//             in memento class, which store data in NATIVE datatypes. This
//             memento objects is then written into HDF5 format. A user need
// to should provide interface to convert data into memento class
//
//-----------------------------------------------------------------------------

    //write out the domain   
    HDF5_WriteDomain(group, eset);

//-----------------------------------------------------------------------------
// write the Vector size
//-----------------------------------------------------------------------------
   dimension[0]=  1;
   try{
      H5::DataSpace sdataspace( rank, dimension );
      H5::DataSet   sdataset = group.createDataSet( "VecSize",
                                   H5::PredType::NATIVE_INT, sdataspace );
      sdataset.write( &size, H5::PredType::NATIVE_INT );
    }
    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}

//-----------------------------------------------------------------------------
// Get the sum of each object size and maximum size of object in the 
// container for allocation purpose
//-----------------------------------------------------------------------------
    entitySet :: const_iterator ci;
    int   offset, bucketID;

    int *vbucket = new int[size*eset.size()];

    size_t  arraySize= 0;
    int     stateSize, maxStateSize = 0;

    hash_map<int, std::vector<T> > ::const_iterator iter;
    std::vector<T>   newVec;

    bucketID = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         iter = attrib_data.find(*ci);
         if( iter != attrib_data.end() ) {
             newVec = iter->second;
             for( int ivec = 0; ivec < size; ivec++){
                  Memento<T> memento( newVec[ivec] );
                  stateSize           = memento.getSize();
                  vbucket[bucketID++] = stateSize;
                  arraySize          += stateSize;
                  maxStateSize        = max( stateSize, maxStateSize);
             }
         }
    }

    typedef hdf5_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *data, *buf;

    data =  new typename converter_traits::memento_type[arraySize];
    buf  =  new typename converter_traits::memento_type[maxStateSize];

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         iter = attrib_data.find(*ci);
         if( iter != attrib_data.end() ) {
             newVec = iter->second;
             for( int ivec = 0; ivec < size; ivec++){
                  Memento<T> memento( newVec[ivec] );
                  memento.getState( buf, stateSize);
                  for( int i = 0; i < stateSize; i++)
                       data[indx++] =  buf[i];
            }
        }
    }

//-----------------------------------------------------------------------------
// Write size of each container ...
//-----------------------------------------------------------------------------
    dimension[0]=  size*eset.size();

    try {

      H5::DataSpace fDataspace( rank, dimension );
      H5::DataType  fDatatype = H5::PredType::NATIVE_INT;
      H5::DataSet   fDataset  = group.createDataSet( "ContainerSize", fDatatype, fDataspace);

      fDataset.write( vbucket, fDatatype );

    }

    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

//-----------------------------------------------------------------------------
// Write (variable) Data into HDF5 format
//-----------------------------------------------------------------------------
    dimension[0]=  arraySize;

    try {

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = converter_traits::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "variable", vDatatype, vDataspace);

      vDataset.write( data, vDatatype );

    }
    catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

    //-----------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------

    delete [] data;
    delete [] buf;

  }

//***************************************************************************

}
#endif

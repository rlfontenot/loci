#ifndef DMAPVEC_H
#define DMAPVEC_H

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <hdf5CC/H5cpp.h>
#include <store.h>
#include <multiMap.h>
#include <Loci_types.h>

namespace Loci {
  template <unsigned int M> class dMapVecRepI : public MapRep {
  public:
    typedef Array<int,M>   VEC;
  private:
    entitySet            store_domain ;
 	 hash_map<int,VEC>    attrib_data;
  public:
    dMapVecRepI() {}
    dMapVecRepI(const entitySet &p) { allocate(p); }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dMapVecRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void compose(const Map &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq)  ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;

    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;

    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group, entitySet &user_eset) ;
    virtual void writehdf5(H5::Group group,entitySet& en) const ;

    hash_map<int,VEC> *get_attrib_data() { return &attrib_data; }
  private:
    int* get_hdf5_data(H5::Group group,const char* datasetname) ;
    void put_hdf5_data(H5::Group group, int* data, const char* datasetname,hsize_t* dimf) const ;
  } ;

//-----------------------------------------------------------------------------

template<unsigned int M> 
void dMapVecRepI<M>::readhdf5( H5::Group group, entitySet &user_eset)
{
    //-------------------------------------------------------------------------
    // Objective   : Reading mapVec Data in HDF5 Format and storing them as
    //               dynamic mapVec
    // Programmer  : Chaman Singh Verma
    // Date        : 29th May, 2001
    // Place       : ERC, Mississippi State University.
    // Status      :
    // Testing     :
    // Limitations :
    // Future Work :
    // Reference   : See the writehdf module for FORMAT.
    //-------------------------------------------------------------------------

    hash_map<int, VEC> :: const_iterator  ci;
    VEC     vec;
    int     numentities;

    H5::DataType  datatype = H5::PredType::NATIVE_INT;

    try{
      //-----------------------------------------------------------------------
      // Part I : read interval sets...
      //-----------------------------------------------------------------------
      H5::DataSet dataset_domain = group.openDataSet( "domain");
      H5::DataSpace dataspace_domain = dataset_domain.getSpace();

      hsize_t dims_domain[1];
      dataspace_domain.getSimpleExtentDims( dims_domain, NULL);

      int *data_domain = new int[dims_domain[0]];

      dataset_domain.read( data_domain,  H5::PredType::NATIVE_INT );

      entitySet eset;	
      for(int i=0;i< dims_domain[0];i++){
        eset |= interval(data_domain[i], data_domain[i+1]);
        i++;
      }
      delete [] data_domain;
 
      //-----------------------------------------------------------------------
      // Part II : Read size of each entity i.e. how many mapped elements on 
      //           each entity.
      //-----------------------------------------------------------------------
      int *range = get_hdf5_data(group,"range");

      if( range[0] != M ) {
          cout << "Error: Reading Invalid mapvec  : " << endl;
          return;
      }
	
      //-----------------------------------------------------------------------
      // Part III: Read the mapped data( in this case entities )
      //-----------------------------------------------------------------------
      hsize_t dims_map[1];

      H5::DataSet dataset_map = group.openDataSet( "mapvec");
      H5::DataSpace fdataspace = dataset_map.getSpace();

      fdataspace.getSimpleExtentDims( dims_map, NULL);    // file dataspace
      
      entitySet en = domain();
      int num_intervals=en.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i<num_intervals;i++) it[i] = en[i];

      //declear the variables used by hyperslab

      hsize_t dimension[1];
      dimension[0] = numentities;
	
      hssize_t  start_mem[] = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace
	
      int rank = 1;
      H5::DataSpace mdataspace(rank,dimension);   // memory dataspace

      //-----------------------------------------------------------------------
      // For each hyperslab, read the data as contiguous chunck and assign
      // values to the multimap.
      //-----------------------------------------------------------------------

      for(int i=0;i<num_intervals;i++){
        count[0]  = 0;          // How many elements in each hyperslab
        for(int j=it[i].first;j<=it[i].second;j++) 
            count[0] += M;

        int *buf = new int[count[0]];
        mdataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);	
        fdataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
        dataset_map.read( buf, datatype, mdataspace, fdataspace);

        int indx = 0;
        int vecsize;
        for(int j=it[i].first;j<=it[i].second;j++) {
          for( int k = 0; k < M; k++)
            attrib_data[j][k] = buf[indx++];
        }
        foffset[0] += count[0];
        delete [] buf;
      }
      delete [] it;
  }
  catch( H5::HDF5DatasetInterfaceException error )  { error.printerror(); }
  catch( H5::HDF5DataspaceInterfaceException error) { error.printerror(); }
  catch( H5::HDF5DatatypeInterfaceException error ) { error.printerror(); }
}

//------------------------------------------------------------------------
    
template<unsigned int M> 
void dMapVecRepI<M>::writehdf5(H5::Group group,entitySet& en) const 
{
    //-------------------------------------------------------------------------
    // Objective   : Write Multimap Data in HDF5 Format.
    // Programmer  : Chaman Singh Verma
    // Date        : 29th May, 2001
    // Place       : ERC, Mississippi State University.
    // Status      :
    // Testing     :
    // Limitations :
    // Future Work :
    // Example     :
    // Let us assume that we have following 
    // entity         mapped entities        
    // 0              1 2 3                    
    // 1              1 2 4                      
    // 2              2 3 5                   
    // 5              5 6 7                     
    // 6              6 7 8                     
    // 10             1 2 3                     
    // 11             2 6 7
    //
    // Then we have 3 intervals set (0,2),(5,6),(10,11)
    // 
    //-------------------------------------------------------------------------

    hash_map<int,VEC> :: const_iterator  ci;
    VEC     vec;
    int     numentity;

    hsize_t dimf_domain[1];
    hsize_t dimf_range[1];
    hsize_t dimension[1];

    int rank = 1;    // One dimensional array

    int num_intervals = en.num_intervals();
    dimf_domain[0]    = num_intervals*2;

    interval *it   = new interval[num_intervals];
    for(int i=0;i<num_intervals;i++) it[i] = en[i];
    //
    // Calculate total number of mapped data which will be written in the file.
    //
    int   range, degree;
    
    numentity = en.size();

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is written.
    hsize_t   count[1];           // how many positions to select from the dataspace

    dimension[0]   = M*numentity;

    H5::DataSpace fdataspace( rank, dimension );     // Dataspace of file
    H5::DataSpace mdataspace( rank, dimension );     // Dataspace of memory 
    H5::DataType  datatype = H5::PredType::NATIVE_INT;

    H5::DataSet   dataset = group.createDataSet( "mapvec", datatype, fdataspace );

    //---------------------------------------------------------------------------
    // For each interval, write data in contiguous fashion. There are 3 data which
    // have to be written. 
    // Write the mapped entities using hyperslab. Each hyperslab consists of
    // contiguous entitySet.
    // From the example.
    // 1 st hyperslab :   (1,1,2,2,3)
    // 2 nd hyperslab :   (5,6,7,6,7)
    // 3 rd hyperslab :   (1,2,3,2,9)
    //---------------------------------------------------------------------------

    for(int i= 0; i< num_intervals; i++){
      count[0] = 0;
      for(int j=it[i].first;j<=it[i].second;j++) {
        ci = attrib_data.find(j);
        if( ci != attrib_data.end() ) count[0] += M;
      }

      int *buf = new int[count[0]];

      mdataspace.selectHyperslab( H5S_SELECT_SET,  count, start_mem, stride, block);
      fdataspace.selectHyperslab( H5S_SELECT_SET,  count, foffset,   stride, block);

      //
      // Copy data in temporary buffer, so that it can be written in single call.
      //
      int indx = 0;
      for(int j=it[i].first;j<=it[i].second;j++) {
        ci = attrib_data.find(j);
        if( ci != attrib_data.end() ) {
          vec  = ci->second;
          for(int k = 0; k < M; k++)
              buf[indx++] =  vec[k];
        }
      }

      dataset.write( buf, datatype, mdataspace, fdataspace);
      //
      // change offset for the second record, and clear the buffer.
      //
      foffset[0] += count[0]; 
      delete[] buf;
  }

  //---------------------------------------------------------------------------
  // Second Part: Write the interval set.
  // From the example, we have three intervals namely (0.2),(5,6),(11,12)
  // Write them in the file as 1D array (0,2,5,6,11,12). 
  //---------------------------------------------------------------------------

  int *buf =  new int[2*num_intervals];
  for(int i=0;i < num_intervals; i++){
    it[i]      = en[i];
    buf[i*2]   = it[i].first;
    buf[i*2+1] = it[i].second;
  }
  put_hdf5_data(group, buf, "domain", dimf_domain );  // domain data
  delete [] buf;
   
  //---------------------------------------------------------------------------
  // Third part: For each entity write down number of mapped entities
  //---------------------------------------------------------------------------

  dimf_range[0]  = 1;

  buf    =  new int[1];
  buf[0] = M;

  put_hdf5_data(group, buf,  "range",  dimf_range  ); 
  delete [] buf;

  //---------------------------------------------------------------------------
  // Clean up :
  //---------------------------------------------------------------------------
  delete [] it;
}

//------------------------------------------------------------------------

template<unsigned int M> 
int* dMapVecRepI<M>::get_hdf5_data(H5::Group group, const char* datasetname)
{
   H5::DataSet dataset = group.openDataSet( datasetname);
   H5::DataSpace dataspace = dataset.getSpace(); 
   hsize_t dims [1];
   dataspace.getSimpleExtentDims(dims , NULL);
   int *data  = new int[dims [0]];
   dataset.read(data , H5::PredType::NATIVE_INT );
   return data;
}

//------------------------------------------------------------------------

template<unsigned int M> 
void dMapVecRepI<M>::put_hdf5_data(H5::Group group, int* data, 
                                   const char* datasetname,hsize_t* dimf) const
{
  int RANK=1;
  try {
     H5::DataSpace dataspace( RANK, dimf );
     H5::DataSet dataset = group.createDataSet( datasetname, 
		                         H5::PredType::NATIVE_INT, dataspace );
     dataset.write( data, H5::PredType::NATIVE_INT );
  }
  catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
  catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
  catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
}

//------------------------------------------------------------------------

template<unsigned int M> 
void dMapVecRepI<M>::allocate(const entitySet &ptn) 
{
   VEC   newvalue;
   entitySet :: const_iterator  ci;

   for( ci = ptn.begin(); ci != ptn.end(); ++ci)
     attrib_data[*ci] = newvalue;
		 
   store_domain = ptn ;
   dispatch_notify() ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
dMapVecRepI<M>::~dMapVecRepI<M>() 
{
   attrib_data.clear();
}

//------------------------------------------------------------------------

template<unsigned int M> 
multiMap dMapVecRepI<M>::get_map()  
{
    multiMap result ;
    cout << "Error: get_map in MapVec Not yet implemented " << endl;
    exit(0);
/*
    store<int> sizes ;
    sizes.allocate(store_domain) ;

    FORALL(store_domain,i) {
      sizes[i] = M ;
    } ENDFORALL ;

    result.allocate(sizes) ;

    FORALL(store_domain,i) {
      for(int j=0;j<M;++j) 
        result.begin(i)[j] = attrib_data[i][j] ;
    } ENDFORALL ;
*/

    return result ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
storeRep *dMapVecRepI<M>::new_store(const entitySet &p) const 
{
    return new dMapVecRepI<M>(p) ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
entitySet dMapVecRepI<M>::domain() const 
{
    hash_map<int,VEC> :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci ) 
         vec.push_back( ci->first ) ;

    sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++) 
         storeDomain +=  vec[i];

    return storeDomain ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
entitySet dMapVecRepI<M>::image(const entitySet &domain) const 
{
    entitySet codomain ;

	/*
    entitySet d = domain & store_domain ;

    for(int i=0;i<d.num_intervals();++i)
      codomain += Loci::image_section(&base_ptr[d[i].first][0],
				      &base_ptr[d[i].second+1][0]) ;

    */
    return codomain ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
std::pair<entitySet,entitySet> dMapVecRepI<M>::preimage(const entitySet &codomain) const 
{
    entitySet domaini ;
    entitySet domainu ;
/*
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false ;
      for(int j=0;j<M;++j) {
        bool in_set = codomain.inSet(base_ptr[i][j]) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
*/
    return std::make_pair(domaini,domainu) ;
}

//------------------------------------------------------------------------
  
template<unsigned int M> 
std::ostream &dMapVecRepI<M>::Print(std::ostream &s) const 
{

    hash_map<int,VEC>  :: const_iterator  ci;

    s << '{' << domain() << std::endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) s << ci->second;
      s << std::endl ;
    } ENDFORALL ;

    s << '}' << std::endl ;

    return s ;
}

//------------------------------------------------------------------------

template<unsigned int M> 
std::istream &dMapVecRepI<M>::Input(std::istream &s) 
{
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      for(int j=0;j<M;++j)
        s >> attrib_data[ii][j] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
}

//------------------------------------------------------------------------

template<unsigned int M> class const_dMapVec ;

//------------------------------------------------------------------------
    
template<unsigned int M> class dMapVec : public store_instance 
{
    friend  class const_dMapVec<M> ;
    typedef dMapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    hash_map<int, VEC>     *attrib_data;
 public:
    dMapVec() { setRep(new MapVecType) ; }
    dMapVec(const dMapVec<M> &var) { setRep(var.Rep()) ; }
    dMapVec(storeRepP &rp) { setRep(rp) ; }

    virtual ~dMapVec() ;

    virtual void notification() ;

    // -----------------------------------------------------------------
        
    dMapVec & operator=(const dMapVec<M> &str) { 
      setRep(str.Rep()) ; 
      return *this ;
    }
    // -----------------------------------------------------------------

    dMapVec &operator=(storeRepP p) { 
      setRep(p) ; 
      return *this ;
    }

    // -----------------------------------------------------------------
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    // -----------------------------------------------------------------

    entitySet domain() const { 
      return Rep()->domain() ; 
    }

    // -----------------------------------------------------------------
    //    operator storeRepP() { return Rep() ; }
    // -----------------------------------------------------------------

    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    // -----------------------------------------------------------------
 
    VEC &elem(int indx) {
      return (*attrib_data)[indx]; 
    }

    // -----------------------------------------------------------------

    const VEC &const_elem(int indx)  const { 
      hash_map<int,VEC> :: const_iterator ci;
      ci = attrib_data.find(index);
      if( ci != attrib_data->end()) return ci->second();
    }

    // -----------------------------------------------------------------

    VEC &operator[](int indx) { 
      return elem(indx); 
    }

    // -----------------------------------------------------------------

    const VEC &operator[](int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }

    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

//  entitySet image(const entitySet &dom) const;

  } ;

  //-----------------------------------------------------------------------

  template<unsigned int M>  
  dMapVec<M>::~dMapVec<M>() { }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  inline std::ostream & operator<<(std::ostream &s, const dMapVec<M> &m)
  {
    return m.Print(s) ;
  }
    
  //------------------------------------------------------------------------

  template<unsigned int M> 
  inline std::istream & operator>>(std::istream &s, dMapVec<M> &m) 
  {
    return m.Input(s) ;
  }

  //------------------------------------------------------------------------

/*
  template<unsigned int M> 
  entitySet image(const entitySet &dom) const 
  {
      cout << "Warn: Not implemented " << endl;
      return( entitySet );
      return MapRepP(Rep())->image(dom) ;
  }
*/

  //------------------------------------------------------------------------

  template<unsigned int M> class const_dMapVec : public store_instance {
    typedef dMapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    hash_map<int,VEC>  *attrib_data;
  public:
    const_dMapVec() { setRep(new MapVecType) ; }
    const_dMapVec(const const_dMapVec<M> &var) { setRep(var.Rep()) ; } 
    const_dMapVec(const dMapVec<M> &var) { setRep(var.Rep()) ; }

    virtual ~const_dMapVec() ;

    virtual void notification() ;

    virtual instance_type access() const ;

    const_dMapVec & operator=(const const_dMapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dMapVec & operator=(const dMapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dMapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }

    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    const VEC &const_elem(int indx)  const {
      hash_map<int,VEC> :: const_iterator ci;
	  
	  ci = attrib_data->find(indx);
	  if( ci != attrib_data->end() ) return( ci->second );
    }

    const VEC &operator[](int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  //------------------------------------------------------------------------

  template<unsigned int M>  
  const_dMapVec<M>::~const_dMapVec() { }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void const_dMapVec<M>::notification()
  {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  store_instance::instance_type
  const_dMapVec<M>::access() const { return READ_ONLY; }
    
  //------------------------------------------------------------------------

  template<unsigned int M> 
  inline std::ostream & operator<<(std::ostream &s,const const_dMapVec<M> &m)
  {
    return m.Print(s) ;
  }
    
  //------------------------------------------------------------------------

  template<unsigned int M> 
  storeRepP dMapVecRepI<M>::remap(const Map &m) const
  {
    entitySet newdomain = m.domain() & domain() ;
    std::pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;

    MapVec<M> s ;
    s.allocate(mapimage) ;

    storeRepP my_store = getRep() ;
      
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    
    return s.Rep() ;
  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVecRepI<M>::compose(const Map &m, const entitySet &context)
  {
    /*
      fatal((context-store_domain) != EMPTY) ;
      fatal((image(context)-m.domain()) != EMPTY) ;

      FORALL(context,i) {
           for(int j=0;j<M;++j)
           base_ptr[i][j] = m[base_ptr[i][j]] ;
      } ENDFORALL ;
    */
  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVecRepI<M>::copy(storeRepP &st, const entitySet &context)
  {
  /*
      const_dMapVec<M> s(st) ;

      FORALL(context,i) {
          attrib_data[i] =  s[i];
      } ENDFORALL ;
  */

  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVecRepI<M>::gather(const Map &m, storeRepP &st, const entitySet &context)
  {
  /*
      const_dMapVec<M> s(st) ;

      FORALL(context,i) {
	      attrib_data[i] = s[m[i]];
      } ENDFORALL ;
  */

  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVecRepI<M>::scatter(const Map &m, storeRepP &st, const entitySet &context)
  {
  /*
      const_dMapVec<M> s(st) ;

      FORALL(context,i) {
         attrib_data[m[i]] = s[i];
      } ENDFORALL ;
  */

  }

  //------------------------------------------------------------------------
  
  template <unsigned int M> 
  int dMapVecRepI<M>::pack_size( const entitySet &e)
  {
    int size ;
    size = sizeof(int) * e.size() * M ;
    return size ;
  }

  //------------------------------------------------------------------------

  template <unsigned int M> 
  void dMapVecRepI<M>::pack(void *ptr, int &loc, int &size, const entitySet &e )
  {
     int          *buf, j, numBytes, numentity = e.size();
     entitySet :: const_iterator ei;

     numBytes = pack_size( e );

     buf = ( int *) malloc( numBytes );

     if( buf == NULL ) {
         cout << "Warning: Cann't allocate memory for packing data " << endl;
         return;
     }

     Array<int,M>   array;
     int indx = 0;

     for( ei = e.begin(); ei != e.end(); ++ei) {
       array  = attrib_data[*ei];
       for( j = 0; j < M; j++) 
            buf[M*indx+j] =  array[j];
       indx++;
     }
     MPI_Pack(buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
     free( buf );
  }

  //------------------------------------------------------------------------

  template <unsigned int M> 
  void dMapVecRepI<M>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {
  /*
    Loci :: int_type   indx, jndx;
    int                numentity, numBytes;
    int                *buf;

    for(int i = 0; i < seq.num_intervals(); ++i) {
        numentity =  abs(seq[i].second - seq[i].first) + 1; 
        numBytes  =  M*numentity*sizeof(T);
        buf       =  (int *) malloc( numBytes );
        MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

        jndx = 0;
        if(seq[i].first > seq[i].second) {
	        for(indx = seq[i].first; indx >= seq[i].second; --indx) {
               attrib_data[indx] =  buf[vecSize*jndx];
               jndx++;
           }
        } else {
	        for(indx = seq[i].first; indx <= seq[i].second; ++indx){
               attrib_data[indx] =  buf[vecSize*jndx];
               jndx++;
           }
        }
        free(buf);
    }
  */
  }  

  //------------------------------------------------------------------------

  /*
    template<unsigned int M> void inverseMap (multiMap &result,
    const dMapVec<M> &input_map,
    const entitySet &input_image,
    const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;

    FORALL(input_image,i) {
    sizes[i] = 0 ;
    } ENDFORALL ;

    entitySet preloop = input_preimage & input_map.domain() ;

    FORALL(preloop,i) {
    for(int k=0;k<M;++k)
    if(input_image.inSet(input_map[i][k]))
    sizes[input_map[i][k]] += 1 ;
    } ENDFORALL ;

    result.allocate(sizes) ;

    FORALL(preloop,i) {
    for(int k=0;k<M;++k) {
    int elem = input_map[i][k] ;
    if(input_image.inSet(elem)) {
    sizes[elem] -= 1 ;
    FATAL(sizes[elem] < 0) ;
    result[elem][sizes[elem]] = i ;
    }

    }
    } ENDFORALL ;

    #ifdef DEBUG
    FORALL(input_image,i) {
    FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
    #endif
    }      
  */
}
//------------------------------------------------------------------------

#endif


#include <DMap.h>
#include <multiMap.h>
#include <Tools/stream.h>
#include <hash_map.h>

namespace Loci {

  using std::pair ;
  using std::make_pair ;

//********************************************************************
        
void dMapRepI::allocate(const entitySet &ptn)
{
  entitySet :: const_iterator   ci;

  for( ci = ptn.begin(); ci != ptn.end(); ++ci)
    attrib_data[*ci] = 0;

  store_domain = ptn ;
		
  dispatch_notify() ;
}

//********************************************************************

dMapRepI::~dMapRepI() 
{
  attrib_data.clear();
}

//********************************************************************

storeRep *dMapRepI::new_store(const entitySet &p) const 
{
  return new dMapRepI(p)  ;
}
//********************************************************************
  
storeRepP dMapRepI::remap(const Map &newmap) const 
{
  dMap s ;

  entitySet newdomain = newmap.domain() & domain() ;

  pair<entitySet,entitySet> mappimage = preimage(newmap.domain()) ;
  newdomain &= mappimage.first ;
  entitySet mapimage = newmap.image(newdomain) ;
  s.allocate(mapimage) ;
  storeRepP my_store = getRep() ;
  s.Rep()->scatter(newmap,my_store,newdomain) ;
  MapRepP(s.Rep())->compose(newmap,mapimage) ;

  return s.Rep() ;
}

//********************************************************************

void dMapRepI::compose(const Map &newmap, const entitySet &context) 
{
  fatal((context-store_domain) != EMPTY) ;
  fatal((image(context)-newmap.domain()) != EMPTY) ;

  FORALL(context,i) {
    attrib_data[i] = newmap[attrib_data[i]] ;
  } ENDFORALL ;

}

//********************************************************************


void dMapRepI::copy(storeRepP &st, const entitySet &context) 
{
    const_dMap s(st) ;

    fatal((context-domain()) != EMPTY) ;

    fatal((context-s.domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;

}

//********************************************************************

void dMapRepI::gather(const Map &m, storeRepP &st, const entitySet &context) 
{
  const_Map s(st) ;

  fatal((m.image(context) - s.domain()) != EMPTY) ; 
  fatal((context - domain()) != EMPTY) ;

  FORALL(context,i) {
      attrib_data[i] = s[m[i]] ;
  } ENDFORALL ;
}

//********************************************************************

void dMapRepI::scatter(const Map &m,storeRepP &st, const entitySet &context) 
{
  const_Map s(st) ;

  fatal((context - s.domain()) != EMPTY) ;
  fatal((m.image(context) - domain()) != EMPTY) ;

  FORALL(context,i) {
      attrib_data[m[i]] = s[i] ;
  } ENDFORALL ;
}

//********************************************************************

  int dMapRepI::pack_size(const entitySet &e) {
    int size ;
    size = sizeof(int) * e.size() ;
    return(size) ;
  }

//********************************************************************

void dMapRepI::pack(void *ptr, int &loc, int &size, const entitySet &e) 
{
    int  *data;
    
    data = ( int *) malloc( e.size()*sizeof(int));
    if( data == NULL ) {
        cout << "Error: Cann't allocate memory for packing data " << endl;
        return;
    }
    
    FORALL(e, i) {
        data[i] = attrib_data[i];
    } ENDFORALL ;

    MPI_Pack(data, e.size()*sizeof(int), MPI_BYTE, ptr, size, &loc, 
             MPI_COMM_WORLD) ;
   
    free ( data );
}

//********************************************************************

void dMapRepI::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
{

    Loci :: int_type   indx, jndx;
    int                numentity, numBytes;
    int                *buf;

    for(int i = 0; i < seq.num_intervals(); ++i) {

        numentity =  abs(seq[i].second - seq[i].first) + 1; 
        numBytes  =  numentity*sizeof(int);
        buf       =  (int *) malloc( numBytes );

        MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

        jndx = 0;
        if(seq[i].first > seq[i].second) {
	        for(indx = seq[i].first; indx >= seq[i].second; --indx) {
               attrib_data[indx] =  buf[jndx++];
           }
        } else {
	        for(indx = seq[i].first; indx <= seq[i].second; ++indx){
               attrib_data[indx] =  buf[jndx++];
           }
        }
        free(buf);
    }
}

//********************************************************************

entitySet dMapRepI::domain() const 
{
    hash_map<int,int > :: const_iterator    ci;
    entitySet          storeDomain;
    vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
         vec.push_back( ci->first ) ;

    sort( vec.begin(), vec.end() );

    for( int i = 0; i < vec.size(); i++)
         storeDomain +=  vec[i];

    return storeDomain ;
}

//********************************************************************

entitySet dMapRepI::image(const entitySet &domain) const 
{

    entitySet codomain ;
    entitySet :: const_iterator  ei;
    hash_map<int,int> ::  const_iterator   ai;   
    
    for( ei = domain.begin(); ei != domain.end(); ++ei){
         ai = attrib_data.find(*ei);
         if( ai != attrib_data.end() )
             codomain +=   ai->second;
    }
    return codomain ;
}

//********************************************************************

pair<entitySet,entitySet>
dMapRepI::preimage(const entitySet &codomain) const  
{
    entitySet domain ;
    hash_map<int,int> :: const_iterator  ci;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci)
         if(codomain.inSet( ci->second) )  domain  += ci->first;
    
    return make_pair(domain,domain);
}

//********************************************************************

multiMap dMapRepI::get_map() 
{
    multiMap result ;
/*

    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      result.begin(i)[0] = base_ptr[i] ;
    } ENDFORALL ;
*/

    return result ;
}

//********************************************************************

ostream &dMapRepI::Print(ostream &s) const 
{
    hash_map<int,int > :: const_iterator ci;

    s << '{' << domain() << endl ;

    FORALL(domain(),ii) {
      ci = attrib_data.find(ii);
      if( ci != attrib_data.end()) {
          s << ci->second << endl ;
      }
    } ENDFORALL ;

    s << '}' << endl ;
    return s ;
}

//********************************************************************

istream &dMapRepI::Input(istream &s) 
{

    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      s >> attrib_data[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }

    return s ;
}

//*********************************************************************************

void dMapRepI::readhdf5(H5::Group group, entitySet &usr_eset)
{

    hash_map<int, int > :: const_iterator  ci;
    entitySet::const_iterator ei;
    hsize_t       dimension;
    entitySet     eset;	
    vector<int>   vec;

    H5::DataType  datatype = H5::PredType::NATIVE_INT;

    HDF5_ReadDomain( group, eset );
      
    try{
      //-----------------------------------------------------------------------
      // Part III: Read the mapped data( in this case entities )
      //-----------------------------------------------------------------------

      H5::DataSet   vDataset   = group.openDataSet( "Map");
      H5::DataSpace vDataspace = vDataset.getSpace();

      vDataspace.getSimpleExtentDims( &dimension, NULL);    // file dataspace
      
      //declear the variables used by hyperslab

      hssize_t  start_mem[] = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      // memory dataspace
      int rank = 1;
      dimension = eset.size();;
      H5::DataSpace mDataspace(rank, &dimension);   // memory dataspace

      //-----------------------------------------------------------------------
      // For each hyperslab, read the data as contiguous chunck and assign
      // values to the multimap.
      //-----------------------------------------------------------------------

      entitySet  ecommon = eset & usr_eset;

      int num_intervals = ecommon.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i<num_intervals;i++) it[i] = ecommon[i];

      //-----------------------------------------------------------------------
      // Calculate the offset of the first element of each interval set
      //-----------------------------------------------------------------------
      store<int> offset;
      offset.allocate( eset );

      int indx = 0;
      for(ei= eset.begin(); ei != eset.end(); ++ei) {
          offset[*ei] = indx;
          indx  += 1;
      }

      for(int i=0;i<num_intervals;i++){
          
          // How many elements in each hyperslab ...
          count[0]  = it[i].second - it[i].first + 1;          
          
          int *buf = new int[count[0]];

          // Read the HDF5 Data ...
          foffset[0] = offset[it[i].first];
          mDataspace.selectHyperslab(H5S_SELECT_SET, count, start_mem, stride, block);	
          vDataspace.selectHyperslab(H5S_SELECT_SET, count, foffset,   stride, block);
          vDataset.read( buf, datatype, mDataspace, vDataspace);

 	       // Create Map ....
          indx = 0;
          for(int j=it[i].first;j<=it[i].second;j++) 
              attrib_data[j] = buf[indx++];

          delete [] buf;
      }

      delete [] it;
    }
    catch( H5::HDF5DatasetInterfaceException error )  { error.printerror(); }
    catch( H5::HDF5DataspaceInterfaceException error) { error.printerror(); }
    catch( H5::HDF5DatatypeInterfaceException error ) { error.printerror(); }
} 


//********************************************************************************

void dMapRepI::writehdf5(H5::Group group, entitySet& eset) const
{
  int       rank = 1;
  hsize_t   dimension;


  HDF5_WriteDomain( group, eset);
   
  //---------------------------------------------------------------------------
  // Third part: For each entity write down number of mapped entities
  // For each entity, write number of degree, so that we can reconstruct the
  // data later.
  // From the example. Write in 1D array (1,2,2,3,3,3,1)
  //---------------------------------------------------------------------------

  int indx = 0;
  int *data = new int[eset.size()];
  entitySet :: const_iterator   ei;
  hash_map<int,int> :: const_iterator  ci;

  for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) 
        data[indx++] =  ci->second;
  }

  dimension = eset.size();
  try {
    H5::DataSpace sDataspace( rank, &dimension );
    H5::DataSet   sDataset = group.createDataSet( "Map", 
                              H5::PredType::NATIVE_INT, sDataspace );
    sDataset.write( data, H5::PredType::NATIVE_INT );
  }

  catch( H5::HDF5DatasetInterfaceException   error ){error.printerror();}
  catch( H5::HDF5DatatypeInterfaceException  error ){error.printerror();} 
  catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

  delete [] data;
} 

//******************************************************************************

  dMap::~dMap() {}

//******************************************************************************

  void dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }    

//******************************************************************************

  const_dMap::~const_dMap() {}

//******************************************************************************

  void const_dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
       attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

//**************************************************************************

  store_instance::instance_type const_dMap::access() const
    { return READ_ONLY ; }

//**************************************************************************

}

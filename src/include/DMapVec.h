#ifndef DMAPVEC_H
#define DMAPVEC_H

#include <istream>
#include <ostream>

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif

#include <Tools/debug.h>
#include <Map_rep.h>
#include <store.h>
#include <multiMap.h>
#include <Loci_types.h>

namespace Loci {
  using std::hash_map ;

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
    virtual void readhdf5( hid_t group, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group,entitySet &en) const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    hash_map<int,VEC> *get_attrib_data() { return &attrib_data; }
  } ;

//-----------------------------------------------------------------------------

template<unsigned int M> 
void dMapVecRepI<M>::readhdf5( hid_t group_id, entitySet &user_eset)
{
    VEC         vec;
    int         size, numentities, rank = 1;
    hsize_t     dimension[1];
    entitySet   eset;
    hash_map<int, VEC> :: const_iterator  ci;

    HDF5_ReadDomain( group_id, eset );

    //--------------------------------------------------------------------------
    // Read the vector size ...
    //--------------------------------------------------------------------------
    dimension[0] = 1;

    hid_t vDataset   = H5Dopen(group_id, "VecSize");
    H5Dread(vDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);

    if( size != M ) {
          cout << "Error: Reading Invalid mapvec  : " << endl;
          return;
    }
    entitySet ecommon;

    ecommon = eset & user_eset;

    allocate( ecommon );

   //---------------------------------------------------------------------------
   // Calculate the offset of each entity in file ....
   //---------------------------------------------------------------------------
   cout << ecommon << endl;
   store<unsigned> offset;
   offset.allocate(eset);
   entitySet :: const_iterator ei;

   int arraySize = 0;
   for( ei = eset.begin(); ei != eset.end(); ++ei) {
        offset[*ei] = arraySize;
        arraySize  += size;
   }

   //---------------------------------------------------------------------------
   // Read the data now ....
   //---------------------------------------------------------------------------

   int num_intervals = ecommon.num_intervals();
   interval *it = new interval[num_intervals];

   for(int i=0;i< num_intervals;i++) it[i] = ecommon[i];

   dimension[0] = arraySize;

   hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);
   hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
   hid_t vDatatype   = H5Tcopy(H5T_NATIVE_INT);

   vDataset   = H5Dopen( group_id, "MapVec");

   hssize_t  start[]     = {0};  // determines the starting coordinates.
   hsize_t   stride[]    = {1};  // which elements are to be selected.
   hsize_t   block[]     = {1};  // size of element block;
   hssize_t  foffset[]   = {0};  // location (in file) where data is read.
   hsize_t   count[]     = {0};  // how many positions to select from the dataspace

   int   *data;

   for( int k = 0; k < num_intervals; k++) {
        count[0] = 0;
        for( int i = it[k].first; i <= it[k].second; i++) 
             count[0] +=  size;

        data = new int[count[0]];

        foffset[0] = offset[it[k].first];

        H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
        H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
        H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);

        int indx = 0;
        for( int i = it[k].first; i <= it[k].second; i++) {
             for( int m = 0; m < size; m++)
                 attrib_data[i][m] = data[indx++];
        }

        delete[] data;
   }

   H5Dclose( vDataset   );
   H5Tclose( vDatatype  );
   H5Sclose( mDataspace );
   H5Sclose( vDataspace );
}
//------------------------------------------------------------------------
    
template<unsigned int M> 
void dMapVecRepI<M>::writehdf5(hid_t group_id, entitySet &eset) const 
{

    hsize_t   dimension[1];
    int       size, rank = 1;

    //write out the domain
    HDF5_WriteDomain(group_id, eset);

    dimension[0]=  1;
    size        =  M;

    hid_t dataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t datatype  = H5Tcopy(H5T_NATIVE_INT);
    hid_t dataset   = H5Dcreate(group_id, "VecSize", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);

//-----------------------------------------------------------------------------
// Collect state data from each object and put into 1D array
//-----------------------------------------------------------------------------

    int arraySize  = size*eset.size();

    int *data =  new int[arraySize];

    entitySet::const_iterator ci;
    hash_map<int, VEC> ::const_iterator iter;
    VEC   newvec;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
         iter = attrib_data.find(*ci);
         if( iter != attrib_data.end() ) {
             newvec = iter->second;
             for( int i = 0; i < size; i++)
                  data[indx++] =  newvec[i];
        }
    }

    dimension[0] = arraySize;
    dataspace = H5Screate_simple(rank, dimension, NULL);
    dataset   = H5Dcreate(group_id, "MapVec", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose( dataset   );
    H5Tclose( datatype  );
    H5Sclose( dataspace );

    delete [] data;
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
storeRepP dMapVecRepI<M>::expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
  int *recv_count = new int[MPI_processes] ;
  int *send_count = new int[MPI_processes] ;
  int *send_displacement = new int[MPI_processes] ;
  int *recv_displacement = new int[MPI_processes] ;
  entitySet::const_iterator ei ;
  std::vector<int>::const_iterator vi ;
  int size_send = 0 ;
  std::vector<std::vector<int> > copy(MPI_processes), send_clone(MPI_processes) ;
  for(int i = 0; i < MPI_processes; ++i) {
    entitySet tmp = out_of_dom & ptn[i] ;
    for(ei = tmp.begin(); ei != tmp.end(); ++ei)
      copy[i].push_back(*ei) ;
    sort(copy[i].begin(), copy[i].end()) ;
    send_count[i] = copy[i].size() ;
    size_send += send_count[i] ; 
  }
  int *send_buf = new int[size_send] ;
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
	       MPI_COMM_WORLD) ; 
  size_send = 0 ;
  for(int i = 0; i < MPI_processes; ++i)
    size_send += recv_count[i] ;
  
  int *recv_buf = new int[size_send] ;
  size_send = 0 ; 
  for(int i = 0; i < MPI_processes; ++i)
    for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
      send_buf[size_send] = *vi ;
      ++size_send ;
    }
  send_displacement[0] = 0 ;
  recv_displacement[0] = 0 ;
  for(int i = 1; i < MPI_processes; ++i) {
    send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
    recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
  }
  MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		recv_buf, recv_count, recv_displacement, MPI_INT,
		MPI_COMM_WORLD) ;  
  for(int i = 0; i < MPI_processes; ++i) {
    for(int j = recv_displacement[i]; j <
	  recv_displacement[i]+recv_count[i]; ++j) 
      send_clone[i].push_back(recv_buf[j]) ;
    sort(send_clone[i].begin(), send_clone[i].end()) ;
  }
  std::vector<hash_map<int, VEC > > map_entities(MPI_processes) ;
  for(int i = 0; i < MPI_processes; ++i) 
    for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
      if(attrib_data.find(*vi) != attrib_data.end())
	(map_entities[i])[*vi] = attrib_data[*vi] ;
  
  for(int i = 0; i < MPI_processes; ++i) {
    send_count[i] = 2 * map_entities[i].size() ;
    for(hash_map<int, VEC >::iterator hi = map_entities[i].begin(); hi != map_entities[i].end(); ++hi)
      send_count[i] += hi->second.size() ; 
  }
  size_send = 0 ;
  for(int i = 0; i < MPI_processes; ++i)
    size_send += send_count[i] ;
  int *send_map = new int[size_send] ;
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
	       MPI_COMM_WORLD) ; 
  size_send = 0 ;
  for(int i = 0; i < MPI_processes; ++i)
    size_send += recv_count[i] ;
  int *recv_map = new int[size_send] ;
  size_send = 0 ;
  for(int i = 0; i < MPI_processes; ++i) 
    for(hash_map<int, std::vector<int> >::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
      send_map[size_send] = miv->first ;
      ++size_send ;
      send_map[size_send] = miv->second.size() ;
      ++size_send ;
      for(int j = 0; j < M; ++j) { 
	send_map[size_send] = (miv->second)[j] ;
	++size_send ;
      }
    }
  send_displacement[0] = 0 ;
  recv_displacement[0] = 0 ;
  for(int i = 1; i < MPI_processes; ++i) {
    send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
    recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
  }
  MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		recv_map, recv_count, recv_displacement, MPI_INT,
		MPI_COMM_WORLD) ;  
  hash_map<int, std::set<int> > hm ;
  std::set<int> ss ;
  for(int i = 0; i < MPI_processes; ++i) {
    for(int j = recv_displacement[i]; j <
	  recv_displacement[i]+recv_count[i]-1; ++j) {
      int count = recv_map[j+1] ;
      if(count)
	for(int k = 0; k < count; ++k)
	  hm[recv_map[j]].insert(recv_map[j+k+2]);
      else
	hm[recv_map[j]] = ss ;
      j += count + 1 ;
    }
  }
  std::vector<int> tmp_vec ;
  for(hash_map<int, std::set<int> >::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
    if(hmi->second.size()) 
      for(std::set<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si)
	attrib_data[hmi->first].push_back(*si) ;
    else
      attrib_data[hmi->first] = tmp_vec ;
  dmultiMap dmul ;
  for(hash_map<int, std::vector<int> >::const_iterator hi = attrib_data.begin(); hi != attrib_data.end(); ++hi)
    dmul[hi->first] = hi->second ;
  storeRepP sp = dmul.Rep() ;
  delete [] send_buf ;
  delete [] recv_buf ;
  delete [] send_map ;
  delete [] recv_map ;
  delete [] recv_count ;
  delete [] send_count ;
  delete [] send_displacement ;
  delete [] recv_displacement ; 
  return sp ;
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
    std::vector<int>        vec;

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

    entitySet d = domain & store_domain ;

    entitySet codomain ;

    entitySet :: const_iterator ei;
   
    unsigned int i;
    VEC    aArray;
    hash_map<int,VEC> :: const_iterator iter;
    for( ei = d.begin(); ei != d.end(); ++ei){
         iter = attrib_data.find( *ei );
         if( iter != attrib_data.end() ) {
             aArray = iter->second;
             for( i = 0; i < M; i++) 
                 codomain +=  aArray[i];
         }
    }
    return codomain;
}

//------------------------------------------------------------------------

template<unsigned int M> 
std::pair<entitySet,entitySet> dMapVecRepI<M>::preimage(const entitySet &codomain) const 
{
    entitySet domaini ;   // intersection 
    entitySet domainu ;   // union
 	 hash_map<int,VEC> :: const_iterator ci;
    VEC  aVec;

    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false ;
      ci = attrib_data.find( i );
      if( ci != attrib_data.end() ) {
          aVec = ci->second;
          for(int j=0;j<M;++j) {
              bool in_set = codomain.inSet(aVec[j]) ;
              vali = vali && in_set ;
              valu = valu || in_set ;
          }
          if(vali)
             domaini += i ;
          if(valu)
            domainu += i ;
      }
    } ENDFORALL ;

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

    dMapVec<M> s ;
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


    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;

    entitySet dom = m.domain() ;

    FORALL(context,i) {
      for(int j=0;j<M;++j) {
	       if(dom.inSet(attrib_data[i][j]))
	           attrib_data[i][j] = m[attrib_data[i][j]] ;
	       else
	           attrib_data[i][j] = -1 ;
      }
    } ENDFORALL ;

  }

  //------------------------------------------------------------------------

  template<unsigned int M> 
  void dMapVecRepI<M>::copy(storeRepP &st, const entitySet &context)
  {
/*
      const_dMapVec<M> s(st) ;

      fatal((context-domain()) != EMPTY) ;
      fatal((context-s.domain()) != EMPTY) ;

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

      fatal((m.image(context) - s.domain()) != EMPTY) ;
      fatal((context - domain()) != EMPTY) ;

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

      fatal((context - s.domain()) != EMPTY) ;
      fatal((m.image(context) - domain()) != EMPTY) ;

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

     int indx = 0;
     for( ei = e.begin(); ei != e.end(); ++ei) {
       for( j = 0; j < M; j++) 
            buf[M*indx+j] =  attrib_data[*ei][j];
       indx++;
     }

     MPI_Pack(buf, numBytes, MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
     free( buf );
  }

  //------------------------------------------------------------------------

  template <unsigned int M> 
  void dMapVecRepI<M>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {
    Loci :: int_type   indx;
    int                numentity, numBytes;
    int                *buf;
    int  vecSize = M;

    for(int i = 0; i < seq.num_intervals(); ++i) {
        numentity =  abs(seq[i].second - seq[i].first) + 1; 
        numBytes  =  M*numentity*sizeof(int);
        buf       =  (int *) malloc( numBytes );
        MPI_Unpack(ptr, size, &loc, buf, numBytes, MPI_BYTE, MPI_COMM_WORLD) ;

        if(seq[i].first > seq[i].second) {
	        for(indx = seq[i].first; indx >= seq[i].second; --indx) {
                for( int j = 0; j < M; j++) 
                     attrib_data[indx][j] =  buf[M*indx+j];
           }
        } else {
	        for(indx = seq[i].first; indx <= seq[i].second; ++indx){
                for( int j = 0; j < M; j++) 
                     attrib_data[indx][j] =  buf[M*indx+j];
           }
        }
        free(buf);
    }
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

#ifndef MAPVEC_H
#define MAPVEC_H

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <store.h>
#include <multiMap.h>
#include <hdf5_readwrite.h>

namespace Loci {

  template <int M> class MapVecRepI : public MapRep {
  public:
    typedef int VEC[M] ;
  private:
    entitySet store_domain ;
    VEC *alloc_pointer ;
    VEC *base_ptr ;
  public:
    MapVecRepI() { alloc_pointer = 0 ; base_ptr = 0 ; }
    MapVecRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~MapVecRepI() ;
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
    virtual void readhdf5( hid_t group, entitySet &en) ;
    virtual void writehdf5(hid_t group,entitySet& en) const ;
    VEC * get_base_ptr() const { return base_ptr ; }
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  } ;
  
  template<int M> void MapVecRepI<M>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int size = ptn.Max()-top+1 ;
      alloc_pointer = new VEC[size] ;
      base_ptr = alloc_pointer-top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  template<int M> MapVecRepI<M>::~MapVecRepI<M>() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  template<int M> multiMap MapVecRepI<M>::get_map()  {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = M ;
    } ENDFORALL ;
    multiMap result ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      for(int j=0;j<M;++j) 
        result.begin(i)[j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return result ;
  }

  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p) const {
    return new MapVecRepI<M>(p) ;
  }

  template<int M> entitySet MapVecRepI<M>::domain() const {
    return store_domain ;
  }

  template<int M> entitySet MapVecRepI<M>::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(int i=0;i<d.num_intervals();++i)
        codomain += Loci::image_section(&base_ptr[d[i].first][0],
                                        &base_ptr[d[i].second+1][0]) ;
    } else {
      std::vector<int> img(d.size()*M) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(int i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          for(int k=0;k<M;++k) {
            *ins = base_ptr[j][k] ;
            ++ins ;
          }
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
    return codomain ;
  }

  template<int M> std::pair<entitySet,entitySet>
  MapVecRepI<M>::preimage(const entitySet &codomain) const {
    entitySet domaini ;
    entitySet domainu ;
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
    return std::make_pair(domaini,domainu) ;
  }

  
  template<int M> std::ostream &MapVecRepI<M>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii][0] ;
      for(int j=1;j<M;++j)
        s << " " << base_ptr[ii][j] ;
      s << std::endl ;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }

  template<int M> std::istream &MapVecRepI<M>::Input(std::istream &s) {
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
        s >> base_ptr[ii][j] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  template<int M> class const_MapVec ;
    
  template<int M> class MapVec : public store_instance {
    friend  class const_MapVec<M> ;
    typedef MapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    VEC * base_ptr ;
  public:
    MapVec() { setRep(new MapVecType) ; }
    MapVec(const MapVec<M> &var) { setRep(var.Rep()) ; }
    MapVec(storeRepP rp) { setRep(rp) ; }

    virtual ~MapVec() ;

    virtual void notification() ;
        
    MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }
    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
    //    operator storeRepP() { return Rep() ; }
    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    VEC &elem(int indx) { fatal(base_ptr==NULL); 
    fatal(!((Rep()->domain()).inSet(indx))) ;
    return base_ptr[indx]; }
    const VEC &const_elem(int indx)  const { fatal(base_ptr==NULL); 
    fatal(!((Rep()->domain()).inSet(indx))) ;
    return base_ptr[indx]; }
    VEC &operator[](int indx) { return elem(indx); }
    const VEC &operator[](int indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  template<int M>  MapVec<M>::~MapVec<M>() { }

  template<int M> void MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  template<int M> inline std::ostream & operator<<(std::ostream &s, const MapVec<M> &m){
    return m.Print(s) ;
  }
    
  template<int M> inline std::istream & operator>>(std::istream &s, MapVec<M> &m) {
    return m.Input(s) ;
  }

  template<int M> class const_MapVec : public store_instance {
    typedef MapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    const VEC * base_ptr ;
  public:
    const_MapVec() { setRep(new MapVecType) ; }
    const_MapVec(const const_MapVec<M> &var) { setRep(var.Rep()) ; } 
    const_MapVec(const MapVec<M> &var) { setRep(var.Rep()) ; }
    const_MapVec(storeRepP rp) { setRep(rp) ; }

    virtual ~const_MapVec() ;

    virtual void notification() ;

    virtual instance_type access() const ;

    const_MapVec & operator=(const const_MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    const entitySet domain() const { return Rep()->domain() ; }
    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }

    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    const VEC &const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const VEC &operator[](int indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  template<int M>  const_MapVec<M>::~const_MapVec() { }

  template<int M> void const_MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  template<int M> store_instance::instance_type
  const_MapVec<M>::access() const { return READ_ONLY; }
    

  template<int M> inline std::ostream & operator<<(std::ostream &s,
                                                   const const_MapVec<M> &m){
    return m.Print(s) ;
  }
    
  template<int M> storeRepP MapVecRepI<M>::remap(const Map &m) const {
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
  template<int M> void MapVecRepI<M>::compose(const Map &m,
                                              const entitySet &context)  {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    entitySet dom = m.domain() ;
    FORALL(context,i) {
      for(int j=0;j<M;++j) {
	if(dom.inSet(base_ptr[i][j]))
          base_ptr[i][j] = m[base_ptr[i][j]] ;
	else
	  base_ptr[i][j] = -1 ;
      }
    } ENDFORALL ;
  }
  template<int M> void MapVecRepI<M>::copy(storeRepP &st,
                                           const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
  }
  template<int M> void MapVecRepI<M>::gather(const Map &m, storeRepP &st,
                                             const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal(base_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;
  }
  template<int M> void MapVecRepI<M>::scatter(const Map &m, storeRepP &st,
                                              const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((base_ptr == 0) &&(context != EMPTY)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[m[i]][j] = s[i][j] ;
    } ENDFORALL ;
  }
  
  template <int M> int MapVecRepI<M>::pack_size( const entitySet &eset) {
    int size ;
    size = sizeof(int)*eset.size()*M  + sizeof(int);
    return size ;
  }

  template <int M> void MapVecRepI<M>::pack(void *outbuf, int &position, 
                 int &outcount, const entitySet &eset ) 
  {
    int init_size = M;
    MPI_Pack( &init_size, 1,  MPI_INT, outbuf, outcount, &position, 
              MPI_COMM_WORLD) ;

    entitySet::const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      MPI_Pack( base_ptr[*ci], M, MPI_INT, outbuf, outcount, &position, 
                MPI_COMM_WORLD);
  }

  template <int M> void MapVecRepI<M>::unpack(void *inbuf, int &position, int &insize, 
                  const sequence &seq) {

    int init_size;

    MPI_Unpack(inbuf, insize, &position, &init_size, 1, MPI_INT, MPI_COMM_WORLD) ;

    if(init_size != M) {
       cout << " Invalid MapVec container for unpack data " << endl;
       abort();
    }

    sequence::const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) 
         MPI_Unpack( inbuf, insize, &position, base_ptr[*ci], M, MPI_INT, 
                     MPI_COMM_WORLD) ;
  } 
  
  template<int M> storeRepP MapVecRepI<M>::expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    storeRepP sp ;
    warn(true) ;
    return sp ;
  }     
  template<int M> void inverseMap (multiMap &result,
                                   const MapVec<M> &input_map,
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

template<int M> 
void MapVecRepI<M>::readhdf5( hid_t group_id, entitySet &eset)
{
  hsize_t dimension[1];
  int  indx=0, rank=1, size;
 
  Loci::HDF5_ReadDomain(group_id, eset);
  Loci::HDF5_ReadVecSize(group_id, &size);

  allocate( eset );

  entitySet::const_iterator ci;

  //---------------------------------------------------------------------------
  // Calculate the offset of each entity in file ....
  //---------------------------------------------------------------------------
  store<unsigned> offset;
  offset.allocate(eset);

  int arraySize = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    offset[*ci] = arraySize;
    arraySize  += size;
  }

  //---------------------------------------------------------------------------
  // Read the data now ....
  //---------------------------------------------------------------------------
  int num_intervals = eset.num_intervals();
  interval *it = new interval[num_intervals];

  for(int i=0;i< num_intervals;i++) it[i] = eset[i];

  int   *data;
  dimension[0] = size*eset.size();

  hid_t vDatatype  = H5T_NATIVE_INT;
  hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);
  hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
  hid_t vDataset   = H5Dopen( group_id, "MapVec");

  hssize_t  start[]     = {0};  // determines the starting coordinates.
  hsize_t   stride[]    = {1};  // which elements are to be selected.
  hsize_t   block[]     = {1};  // size of element block;
  hssize_t  foffset[]   = {0};  // location (in file) where data is read.
  hsize_t   count[]     = {0};  // how many positions to select from the dataspace

  int voffset;
  for( int k = 0; k < num_intervals; k++) {
    count[0] = 0;
    for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  size;

    data = new int[count[0]];

    foffset[0] = offset[it[k].first];

    H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
    H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
    H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);
 
    indx = 0;
    for( int i = it[k].first; i <= it[k].second; i++) {
      for( int ivec = 0; ivec < size; ivec++){
        voffset           = i*size+ivec;
        base_ptr[i][ivec] = data[indx++];
      }
    }
    delete[] data;
  }

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Sclose( mDataspace);

}

//*************************************************************************
    
template<int M> 
void MapVecRepI<M>::writehdf5(hid_t group_id,entitySet& usr_eset) const 
{
  int     rank = 1;
  int     voffset;
  hsize_t dimension;

  entitySet eset(usr_eset &domain());
  int arraySize = M*eset.size();
  if( arraySize < 1) return;

  Loci::HDF5_WriteDomain( group_id, eset);
  Loci::HDF5_WriteVecSize(group_id, M);

  std::vector<int> data(arraySize);
  entitySet::const_iterator ci;

  size_t indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    for( int ivec = 0; ivec < M; ivec++)
      data[indx++] =  base_ptr[*ci][ivec];
  }

  dimension        = arraySize;
  hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
  hid_t vDatatype  = H5T_NATIVE_INT;
  hid_t vDataset   = H5Dcreate(group_id, "MapVec", vDatatype, vDataspace, 
                               H5P_DEFAULT);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
}

}

#endif

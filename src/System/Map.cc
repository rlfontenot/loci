#include <Map.h>
#include <multiMap.h>
#include <iostream>

namespace Loci {

  using std::pair ;
  using std::make_pair ;

  MapRep::~MapRep() {}
  
  store_type MapRep::RepType() const { return MAP ; }
  
  void MapRepI::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int size = ptn.Max()-top+1 ;
      alloc_pointer = new(Entity[size]) ;
      base_ptr = alloc_pointer - top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }


  MapRepI::~MapRepI() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  storeRep *MapRepI::new_store(const entitySet &p) const {
    return new MapRepI(p)  ;
  }
  
  storeRepP MapRepI::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;
    Map s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    return s.Rep() ;
  }

  void MapRepI::compose(const Map &m, const entitySet &context) {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = m[base_ptr[i]] ;
    } ENDFORALL ;
  }

  void MapRepI::copy(storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }

  void MapRepI::gather(const Map &m, storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ; 
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  void MapRepI::scatter(const Map &m,storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }
 
  int MapRepI::pack_size(const entitySet &e) {
    int size ;
    size = sizeof(int) * e.size() ;
    return(size) ;
  }
  
  void MapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
    for( int i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      MPI_Pack( &base_ptr[begin], t, MPI_INT, outbuf, outcount, 
                &position, MPI_COMM_WORLD) ;
    }
  }
  
  void MapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) {

    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                      1 , MPI_INT, MPI_COMM_WORLD) ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                    t, MPI_INT, MPI_COMM_WORLD) ;
      }
    }
  }

  entitySet MapRepI::domain() const {
    return store_domain ;
  }

  entitySet image_section(const Entity *start, const Entity *end) {
    if(start == end)
      return EMPTY ;
    int mx,mn ;
    mx = *start ;
    mn = *start ;
    for(const Entity *i=start;i!=end;++i) {
      mx = max(mx,*i) ;
      mn = min(mn,*i) ;
    }                                         
    int sz = mx-mn+1 ;
    int sz2 = end-start ;
    if(sz>2*sz2) {
      // If the image is too sparse then we are better off sorting
      // using standard sort
      entitySet dom ;
      std::vector<int> img(sz2) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(const Entity *i=start;i!=end;++i) {
        *ins = *i ;
        ++ins ;
      }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        dom += *ins ;
      return dom ;
    }
    std::vector<bool> bits(sz) ;
    for(int i=0;i<sz;++i)
      bits[i] = false ;
    for(const int *i=start;i!=end;++i)
      bits[*i-mn] = true ;

    WARN(!bits[0]);

    entitySet result ;
    interval iv(mn,mn) ;
    for(int i=0;i<sz;++i)
      if(!bits[i]) {
        iv.second = i+mn-1 ;
        result += iv ;
        for(;i<sz;++i)
          if(bits[i])
            break ;
        iv.first = i+mn ;
      }

    WARN(!bits[sz-1]) ;
    
    iv.second = mx ;
    result += iv ;
    return result ;
  }
      
        
    
  entitySet MapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(int i=0;i<d.num_intervals();++i)
        codomain += image_section(base_ptr+d[i].first,
                                  base_ptr+(d[i].second+1)) ;
    } else {
      std::vector<int> img(d.size()) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(int i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          *ins = base_ptr[j] ;
          ++ins ;
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
      
    return codomain ;
  }

  pair<entitySet,entitySet>
  MapRepI::preimage(const entitySet &codomain) const  {
    entitySet domain ;
    FORALL(store_domain,i) {
      if(codomain.inSet(base_ptr[i]))
        domain += i ;
    } ENDFORALL ;
    return make_pair(domain,domain) ;
  }
  
  storeRepP MapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    storeRepP sp ;
    warn(true) ;
    return sp ;
  }
  multiMap MapRepI::get_map() {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = 1 ;
    } ENDFORALL ;
    multiMap result ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      result.begin(i)[0] = base_ptr[i] ;
    } ENDFORALL ;
    return result ;
  }
    
  std::ostream &MapRepI::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << std::endl ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  std::istream &MapRepI::Input(std::istream &s) {
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
      s >> base_ptr[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  void MapRepI::readhdf5(hid_t group_id, entitySet &usr_eset){

    int      rank = 1, indx = 0;
    hid_t    mDataspace, vDataspace, vDataset, vDatatype;
    hsize_t  dimension;
    entitySet  eset;
    entitySet::const_iterator  ci;

    HDF5_ReadDomain(group_id, eset);

    store<int>  offset;
    offset.allocate( eset );

    indx     = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      offset[*ci] = indx++;

    //--------------------------------------------------------------------
    // Read the data now ....
    //--------------------------------------------------------------------
    int num_intervals = usr_eset.num_intervals();

    if( num_intervals == 0) {
      std::cout << "Warning: Number of intervals are zero : " << endl;
      return;
    }

    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

    vDatatype = H5T_NATIVE_INT;

    std::vector<int>  data;

    dimension  = eset.size();
    mDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dopen( group_id, "Map");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] += 1;

      if( count[0] > data.size()) data.resize(count[0]);

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride,
                          count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride,
                          count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace,
               H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) 
        base_ptr[i] = data[indx++];
    }

    H5Sclose( mDataspace );
    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

  } 

  void MapRepI::writehdf5(hid_t group_id,entitySet &usr_eset) const{

    entitySet eset(usr_eset&domain());

    HDF5_WriteDomain(group_id, eset);

    int arraySize =  eset.size(); 
    if( arraySize < 1) return;

    hid_t vDatatype = H5T_NATIVE_INT;

    int rank = 1;
    hsize_t  dimension = arraySize;

    entitySet :: const_iterator ci;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);

    std::vector<int> data(arraySize);
    int indx =0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      data[indx++] = base_ptr[*ci];

    hid_t cparms   = H5Pcreate (H5P_DATASET_CREATE);
    hid_t vDataset = H5Dcreate(group_id, "Map", vDatatype,
                               vDataspace, cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);

  } 

  Map::~Map() {}

  void Map::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }    

  const_Map::~const_Map() {}

  void const_Map::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  store_instance::instance_type const_Map::access() const
  { return READ_ONLY ; }
    
  void multiMapRepI::allocate(const entitySet &ptn) {
    store<int> count ;
    count.allocate(ptn) ;
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    allocate(count) ;
  }
  storeRepP multiMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    storeRepP sp ;
    warn(true) ;
    return sp ;
  }
  void multiMapRepI::allocate(const store<int> &sizes) {
    int sz = 0 ;
    entitySet ptn = sizes.domain() ;
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    if(index) delete[] index ;
    index = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int len = ptn.Max() - top + 2 ;
      index = new int *[len] ;
      base_ptr = index - top ;
      FORALL(ptn,i) {
        sz += sizes[i] ;
      } ENDFORALL ;
      alloc_pointer = new int[sz+1] ;
      sz = 0 ;
      for(int ivl=0;ivl<ptn.num_intervals();++ivl) {
        int i = ptn[ivl].first ;
        base_ptr[i] = alloc_pointer + sz ;
        while(i<=ptn[ivl].second) {
          sz += sizes[i] ;
          ++i ;
          base_ptr[i] = alloc_pointer + sz ;
        }
      }
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }


  multiMapRepI::~multiMapRepI() {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
  }

  storeRep *multiMapRepI::new_store(const entitySet &p) const {
    warn(true) ;
    return new multiMapRepI()  ;
  }

  storeRepP multiMapRepI::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    //    pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    //    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;
    multiMap s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    return s.Rep() ;
  } 

  void multiMapRepI::compose(const Map &m, const entitySet &context) {
    fatal(alloc_pointer == 0) ;
    fatal((context-store_domain) != EMPTY) ;
    //fatal((image(context)-m.domain()) != EMPTY) ;
    entitySet dom = m.domain() ;
    FORALL(context,i) {
      for(int *ii = base_ptr[i];ii!=base_ptr[i+1];++ii) {
        if(dom.inSet(*ii))
          *ii = m[*ii] ;
        else
          *ii = -1 ;
      }
    } ENDFORALL ;
  }

  void multialloc(const store<int> &count, int ***index, int **alloc_pointer,
                  int ***base_ptr) {
    entitySet ptn = count.domain() ;
    int top = ptn.Min() ;
    int len = ptn.Max() - top + 2 ;
    int **new_index = new int *[len] ;
    int **new_base_ptr = new_index-top ;
    int sz = 0 ;
    FORALL(ptn,i) {
      sz += count[i] ;
    } ENDFORALL ;
    int *new_alloc_pointer = new int[sz+1] ;
    sz = 0 ;
    for(int ivl=0;ivl<ptn.num_intervals();++ivl) {
      int i = ptn[ivl].first ;
      new_base_ptr[i] = new_alloc_pointer + sz ;
      while(i<=ptn[ivl].second) {
        sz += count[i] ;
        ++i ;
        new_base_ptr[i] = new_alloc_pointer + sz ;
      }
    }
    *index = new_index ;
    *alloc_pointer = new_alloc_pointer ;
    *base_ptr = new_base_ptr ;
  }
    
  void multiMapRepI::copy(storeRepP &st, const entitySet &context) {
    const_multiMap s(st) ;
    fatal(alloc_pointer == 0) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;

    store<int> count ;
    count.allocate(domain()) ;
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[i] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    
    int **new_index ;
    int *new_alloc_pointer ;
    int **new_base_ptr ;

    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  void multiMapRepI::gather(const Map &m, storeRepP &st,
                            const entitySet  &context) {
    store<int> count ;
    const_multiMap s(st) ;
    count.allocate(domain()) ;
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[i] = s.end(m[i])-s.begin(m[i]) ;
    } ENDFORALL ;
    int **new_index ;
    int *new_alloc_pointer ;
    int **new_base_ptr ;

    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
 
  void multiMapRepI::scatter(const Map &m, storeRepP &st,
                             const entitySet  &context) {
    store<int> count ;
    const_multiMap s(st) ;
    count.allocate(domain()) ;
    FORALL(domain()-m.image(context),i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[m[i]] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    int **new_index ;
    int *new_alloc_pointer ;
    int **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-m.image(context),i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    FORALL(context,i) {
      for(int j=0;j<count[m[i]];++j) {
        new_base_ptr[m[i]][j] = s[i][j] ;
      }
    } ENDFORALL ;
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
  
  int multiMapRepI::pack_size(const  entitySet &eset ) {

    int size = 0 ;
    FORALL(eset,i) {
      size += end(i) - begin(i);
    } ENDFORALL ;
    
    return( (size+eset.size())*sizeof(int) ) ;
  }

  void multiMapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) {

    int vsize;
    entitySet :: const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize    = end(*ci) - begin(*ci);
      MPI_Pack( &vsize, 1, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
      MPI_Pack( &base_ptr[*ci], vsize, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
    }

  }
  
  void multiMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) {

    int vsize;
    sequence:: const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
          MPI_Unpack( inbuf, insize, &position, &vsize,
                      1, MPI_INT, MPI_COMM_WORLD) ;
          MPI_Unpack( inbuf, insize, &position, &base_ptr[*ci],
                      vsize, MPI_INT, MPI_COMM_WORLD) ;
    }
  }   
    
  entitySet multiMapRepI::domain() const {
    return store_domain ;
  }
    
  entitySet multiMapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(int i=0;i<d.num_intervals();++i)
        codomain += image_section(begin(d[i].first),end(d[i].second)) ;
    } else {
      int sz = 0 ;
      for(int i=0;i<d.num_intervals();++i)
        sz += end(d[i].second)-begin(d[i].first) ;

      std::vector<int> img(sz) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(int i=0;i<d.num_intervals();++i)
        for(const Entity *j=begin(d[i].first);j!=end(d[i].second);++j) {
          *ins = *j ;
          ++ins ;
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
    return codomain ;
  }

  pair<entitySet,entitySet>
  multiMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domaini,domainu ;
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = begin(i) == end(i)?true:false ;
      for(const int *ip = begin(i);ip!= end(i);++ip) {
        bool in_set = codomain.inSet(*ip) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return make_pair(domaini,domainu) ;
  }

  multiMap multiMapRepI::get_map() {
    return multiMap(storeRepP(this)) ;
  }
    
  std::ostream &multiMapRepI::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;
    FORALL(domain(),ii) {
      for(const int *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << " " ;
      s << std::endl;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  std::istream &multiMapRepI::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    allocate(sizes) ;
        
    FORALL(e,ii) {
      for(int *ip = begin(ii);ip!=end(ii);++ip)
        s >> *ip  ;
    } ENDFORALL ;
            
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  void multiMapRepI::readhdf5( hid_t group_id, entitySet &usr_eset) {

    hsize_t dimension;
    hid_t   vDatatype, vDataset, vDataspace, mDataspace;
    size_t indx = 0, arraySize;
    int    rank = 1;

    entitySet eset;
    entitySet::const_iterator ci;

    HDF5_ReadDomain(group_id, eset);

    //-------------------------------------------------------------------------
    // Size of each main container....
    //--------------------------------------------------------------------------
    vDatatype  = H5T_NATIVE_INT;
    vDataset   = H5Dopen(group_id,"ContainerSize");
    vDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    int *ibuf = new int[dimension];
    H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

    store<int> container;
    container.allocate( eset );

    indx  = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      container[*ci] = ibuf[indx++];

    delete [] ibuf;
    H5Dclose(vDataset);
    H5Sclose(vDataspace);

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<int>        bucket;
    bucket.allocate( usr_eset );

    for( ci = usr_eset.begin(); ci != usr_eset.end(); ++ci) 
      bucket[*ci] = container[*ci];
    allocate( bucket );

    store<unsigned>   offset;
    offset.allocate( eset );

    arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += container[*ci];
    }

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = usr_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

    dimension  = arraySize;
    mDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dopen(group_id,"Map");

    vDatatype = H5T_NATIVE_INT;

    std::vector<int> data;

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  container[i];

      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int j = 0; j < container[i]; j++) 
          base_ptr[i][j] = data[indx++];
      }
    }
    H5Dclose(vDataset);
    H5Sclose(vDataspace);
    H5Sclose(mDataspace);

  }

  void multiMapRepI::writehdf5( hid_t group_id, entitySet& eset) const{
    HDF5_WriteDomain(group_id, eset);

    hid_t  vDatatype, vDataset, vDataspace;

    entitySet :: const_iterator ci;

    //---------------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the 
    // container for allocation purpose
    //---------------------------------------------------------------------------
    int     count, newsize;

    std::vector<int>  container(eset.size());

    size_t indx = 0, arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      newsize    = end(*ci) - begin(*ci);
      arraySize  += newsize;
      container[indx++] = newsize;
    }
         
    //------------------------------------------------------------------------
    // Write the Size of each multiStore ....
    //------------------------------------------------------------------------
    int rank = 1;
    hsize_t  dimension;
    hid_t cparms  = H5Pcreate (H5P_DATASET_CREATE);

    dimension  = eset.size();
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5T_NATIVE_INT;
    vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                            cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);

    //--------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //--------------------------------------------------------------------------
    std::vector<int>  data(arraySize);

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      count  = end(*ci) - begin(*ci);
      for( int j = 0; j < count; j++) 
        data[indx++] = base_ptr[*ci][j];
    }

    //--------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //--------------------------------------------------------------------------
    dimension  = arraySize;
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5T_NATIVE_INT;
    vDataset   = H5Dcreate( group_id, "Map", vDatatype, vDataspace,
                            cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);

  } 
  
  multiMap::~multiMap() {}

  void multiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  const_multiMap::~const_multiMap() { }

  void const_multiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  store_instance::instance_type const_multiMap::access() const
  { return READ_ONLY ; }
    
  void inverseMap(multiMap &result, const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;

    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      if(input_image.inSet(input_map[i]))
        sizes[input_map[i]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) {
        sizes[elem] -= 1 ;
        FATAL(sizes[elem] < 0) ;
        result[elem][sizes[elem]] = i ;
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }


  void inverseMap(multiMap &result, const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;

    FORALL(preloop,i) {
      for(const int *mi=input_map.begin(i);mi!=input_map.end(i);++mi)
        if(input_image.inSet(*mi))
          sizes[*mi] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(const int *mi=input_map.begin(i);mi!=input_map.end(i);++mi) {
        int elem = *mi ;
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
   
}

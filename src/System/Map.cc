#include <Map.h>
#include <multiMap.h>

#include <Tools/stream.h>

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
  
  void MapRepI::pack(void * ptr, int &loc, int &size, const entitySet &e) {
    warn(true) ; 
    /*
      FORALL(e, i) {
      MPI_Pack(base_ptr[i], sizeof(int), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
      } ENDFORALL ;
      
    */
      
  }
  
  void MapRepI::unpack(void * ptr, int &loc, int &size, const sequence &seq) {
    warn(true) ;
    /*
      loc = 0 ;
    */ 
      
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
    if(sz>3*sz2) {
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
    
  ostream &MapRepI::Print(ostream &s) const {
    s << '{' << domain() << endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << endl ;
    }ENDFORALL ;
    s << '}' << endl ;
    return s ;
  }


  istream &MapRepI::Input(istream &s) {
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
      s >> base_ptr[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  void MapRepI::readhdf5(H5::Group group, entitySet &eset){
      hsize_t dims_map[1];
      
      try{
	H5::DataSet dataset_domain = group.openDataSet( "domain");//get domain data
	H5::DataSpace dataspace_domain = dataset_domain.getSpace();
	hsize_t dims_domain[1];
	dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
	int *data_domain = new int[dims_domain[0]];
	dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
	entitySet num;	
	for(int i=0;i<dims_domain[0];i++){
	  num |=interval(data_domain[i],data_domain[i+1]);
	  i++;
	}
	allocate(num);

	H5::DataSet dataset_map = group.openDataSet( "map");//get map data  
	H5::DataSpace dataspace_map = dataset_map.getSpace();
	dataspace_map.getSimpleExtentDims( dims_map, NULL);
	int RANK = dataspace_map.getSimpleExtentNdims();
	int *data_map = new int[dims_map[0]]; 
	dataset_map.read( data_map, H5::PredType::NATIVE_INT );
	//set the base_ptr
	entitySet en=domain();
	int num_intervals=en.num_intervals();
	interval *it = new interval[num_intervals];
	//int x=0;//get the map data 
	int bound=0;
	for(int i=0;i<num_intervals;i++){
	  it[i]=en[i];
	}
      if(en.Min()<0&&en.Max()>0){
	bound=en.Max()-en.Min()+1;
      }
      else if(en.Min()<0){
	bound=abs(en.Min());
      }
      else
	bound=en.Max();
	//for negative domain, we do the coord. transformation
      if(en.Min()<0){
	for(int i=0;i<num_intervals;i++){
	  it[i].first=it[i].first+abs(en.Min());
	  it[i].second=it[i].second+abs(en.Min());
	}
      }
	//declear the variables used by hyperslab
	hsize_t dim_mem[1];
	dim_mem[0]=bound;
	hssize_t start_mem[1];
	hsize_t stride_mem[1];
	hsize_t count_mem[1];
	hsize_t block_mem[1];
	hssize_t start_file[1];
	start_file[0]=0;
	stride_mem[0]=1;	    
	block_mem[0]=1;

	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	for(int i=0;i<num_intervals;i++){
	    start_mem[0]=it[i].first;
	    count_mem[0]=it[i].second-it[i].first+1;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    if(en.Min()<0){//negative domain
	      dataset_map.read(base_ptr+en.Min(),H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	    }
	    else//positive domain
	      dataset_map.read(base_ptr,H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	}
      } 
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
    } 

  void MapRepI::writehdf5(H5::Group group,entitySet& en) const{
    //entitySet en=domain();
      int dim=en.size();
      //int *data_map = new int[dim];
      hsize_t dimf[1];
      dimf[0]=dim;
      int RANK=1;
      hsize_t dimf_domain[1];

      int num_intervals=en.num_intervals();
      dimf_domain[0]=num_intervals*2;
      interval *it = new interval[num_intervals];
      int bound=0;
      int *data_domain = new int[num_intervals*2];//get the domain data
      for(int i=0;i<num_intervals;i++){
        it[i]=en[i];
	data_domain[i*2]=it[i].first;
	data_domain[i*2+1]=it[i].second;
      }
      if(en.Min()<0&&en.Max()>0){
	bound=en.Max()-en.Min()+1;
      }
      else if(en.Min()<0){
	bound=abs(en.Min());
      }
      else
	bound=en.Max();

      //for negative domain, we do the coord. transformation
      if(en.Min()<0){
	for(int i=0;i<num_intervals;i++){
	  it[i].first=it[i].first+abs(en.Min());
	  it[i].second=it[i].second+abs(en.Min());
	}
      }
       //declear the variables used by hyperslab
      hsize_t dim_mem[1];
      dim_mem[0]=bound;
      hssize_t start_mem[1];
      hsize_t stride_mem[1];
      hsize_t count_mem[1];
      hsize_t block_mem[1];
      hssize_t start_file[1];
      start_file[0]=0;
      stride_mem[0]=1;	    
      block_mem[0]=1;

      try{
	//create the domain data
	H5::DataSpace dataspace_domain( RANK, dimf_domain );
	H5::DataSet dataset_domain = group.createDataSet( "domain", H5::PredType::NATIVE_INT, dataspace_domain );
	dataset_domain.write( data_domain, H5::PredType::NATIVE_INT );
	//create the map data
	H5::DataSpace dataspace_map( RANK, dimf );
	H5::DataSet dataset_map = group.createDataSet( "map", H5::PredType::NATIVE_INT, dataspace_map );
	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);

	for(int i=0;i<num_intervals;i++){
	    start_mem[0]=it[i].first;
	    count_mem[0]=it[i].second-it[i].first+1;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    if(en.Min()<0){//negative domain
	      dataset_map.write(base_ptr+en.Min(),H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	    }
	    else//positive domain
	      dataset_map.write(base_ptr,H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	}
	//dataset_map.write( data_map, H5::PredType::NATIVE_INT );
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();} 
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
  
  int multiMapRepI::pack_size(const  entitySet &e ) {
    int size = 0 ;
    store<int> count ;
    count.allocate(e) ;
    FORALL(e,i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
      size += count[i] ;
    } ENDFORALL ;
    
    return( size * sizeof(int)) ;
  }
  void multiMapRepI::pack(void * ptr, int &loc, int &size, const entitySet &e) {
    warn(true) ;
  }
  
  void multiMapRepI::unpack(void * ptr, int &loc, int &size, const sequence &seq) {
    warn(true) ;
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
    
  ostream &multiMapRepI::Print(ostream &s) const {
    s << '{' << domain() << endl ;
    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << endl ;
    } ENDFORALL ;
    FORALL(domain(),ii) {
      for(const int *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << " " ;
      s << endl;
    } ENDFORALL ;
    s << '}' << endl ;
    return s ;
  }


  istream &multiMapRepI::Input(istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
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
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }
    return s ;
  }

    void multiMapRepI::readhdf5( H5::Group group, entitySet &en) {
      try{
	//get domain data
	H5::DataSet dataset_domain = group.openDataSet( "domain");
	H5::DataSpace dataspace_domain = dataset_domain.getSpace();
	hsize_t dims_domain[1];
	dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
	int *data_domain = new int[dims_domain[0]];
	dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
	entitySet num;	
	for(int i=0;i<dims_domain[0];i++){
	  num |=interval(data_domain[i],data_domain[i+1]);
	  i++;
	}

	store<int> sizes;
	sizes.allocate(num);
	
	int *range = get_hdf5_data(group,"range");//get range data
	
	entitySet::const_iterator ii;
	int ij=0;
	int bound=0;
	for(ii=num.begin();ii!=num.end();++ii){
	  sizes[*ii]=range[ij];
	  bound+=range[ij];
	  ij++;
	}
	allocate(sizes);
	
	//get map data 
	hsize_t dims_map[1];
	H5::DataSet dataset_map = group.openDataSet( "multimap");
	H5::DataSpace dataspace_map = dataset_map.getSpace();
	dataspace_map.getSimpleExtentDims( dims_map, NULL);
	int RANK = dataspace_map.getSimpleExtentNdims();

	//set the base_ptr
	entitySet en=domain();
	int num_intervals=en.num_intervals();
	interval *it = new interval[num_intervals];
	for(int i=0;i<num_intervals;i++){
	  it[i]=en[i];
	}

	//declear the variables used by hyperslab
	hsize_t dim_mem[1];
	dim_mem[0]=bound;
	hssize_t start_mem[1];
	hsize_t stride_mem[1];
	hsize_t count_mem[1];
	hsize_t block_mem[1];
	hssize_t start_file[1];
	start_file[0]=0;
	stride_mem[0]=1;	    
	block_mem[0]=1;
	start_mem[0]=0;
	count_mem[0]=0;
	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	for(int i=0;i<num_intervals;i++){
	  for(int j=it[i].first;j<=it[i].second;j++){
	    count_mem[0]+=end(j)-begin(j);
	  }
	  dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	  dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	  start_file[0]=start_file[0]+count_mem[0];//for next interval
	  dataset_map.read(base_ptr[it[i].first],H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	  count_mem[0]=0;
	}

      //reclaim memory
      delete [] it;
      delete [] data_domain;
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
    }

    void multiMapRepI::writehdf5( H5::Group group,entitySet& en) const{
      //entitySet en=domain();
      hsize_t dimf_domain[1];
      hsize_t dimf_range[1];
      hsize_t dimf_map[1];

      int RANK=1;
      int num_intervals=en.num_intervals();
      dimf_domain[0]=num_intervals*2;
      interval *it = new interval[num_intervals];
      int *data_domain = new int[num_intervals*2];//get the domain data
      for(int i=0;i<num_intervals;i++){
        it[i]=en[i];
	data_domain[i*2]=it[i].first;
	data_domain[i*2+1]=it[i].second;
      }
      int x=0,y=0;//get the map data
      int range;
      int bound=0;
      for(int i=0;i<num_intervals;i++){
	for(int j=it[i].first;j<=it[i].second;j++){
	  range=end(j)-begin(j);
	  bound+=range;
	  y++;
	  x+=range;
	}
      }

      //declare the variables used by hyperslab
      hsize_t dim_mem[1];
      dim_mem[0]=bound;
      hssize_t start_mem[1];
      hsize_t stride_mem[1];
      hsize_t count_mem[1];
      hsize_t block_mem[1];
      hssize_t start_file[1];
      start_file[0]=0;
      stride_mem[0]=1;	    
      block_mem[0]=1;
      start_mem[0]=0;
      dimf_map[0]=x;
      H5::DataSpace dataspace_map( RANK, dimf_map );
      H5::DataSet dataset_map = group.createDataSet( "multimap", H5::PredType::NATIVE_INT, dataspace_map );
      //get data using HDF5 hyperslab	
      H5::DataSpace dataspace_memory(RANK,dim_mem);
      int *data_map = new int[x];
      int *data_range = new int[y];
      dimf_range[0]=y;
      count_mem[0]=0;
      y=0;//reset y for data_range
      for(int i=0;i<num_intervals;i++){
	for(int j=it[i].first;j<=it[i].second;j++){
	  range=end(j)-begin(j);
	  data_range[y]=range;
	  y++;
	  count_mem[0]+=range;
	}
	dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);
	dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	start_file[0]=start_file[0]+count_mem[0];//for next interval
	dataset_map.write(base_ptr[it[i].first],H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	count_mem[0]=0;
      }
      put_hdf5_data(group,data_domain,"domain",dimf_domain);//domain data
      put_hdf5_data(group,data_range,"range",dimf_range);//range data

      //reclaim memory
      delete [] it;
      delete [] data_domain;
    } 
  
  int* multiMapRepI::get_hdf5_data(H5::Group group,const char* datasetname){
      H5::DataSet dataset = group.openDataSet( datasetname);
      H5::DataSpace dataspace = dataset.getSpace();
      hsize_t dims [1];
      dataspace.getSimpleExtentDims(dims , NULL);
      int *data  = new int[dims [0]];
      dataset.read(data , H5::PredType::NATIVE_INT );
      return data;
    }
  
    void multiMapRepI::put_hdf5_data(H5::Group group, int* data, const char* datasetname,hsize_t* dimf) const{
       int RANK=1;
      try{
	H5::DataSpace dataspace( RANK, dimf );
	H5::DataSet dataset = group.createDataSet( datasetname, H5::PredType::NATIVE_INT, dataspace );
	dataset.write( data, H5::PredType::NATIVE_INT );
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
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

#ifndef MAPVEC_H
#define MAPVEC_H

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <hdf5CC/H5cpp.h>
#include <store.h>
#include <multiMap.h>

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
    virtual void readhdf5( H5::Group group, entitySet &en) ;
    virtual void writehdf5(H5::Group group,entitySet& en) const ;
    VEC * get_base_ptr() const { return base_ptr ; }
  private:
    int* get_hdf5_data(H5::Group group,const char* datasetname) ;
    void put_hdf5_data(H5::Group group, int* data, const char* datasetname,hsize_t* dimf) const ;
  } ;

 template<int M> void MapVecRepI<M>::readhdf5( H5::Group group, entitySet &eset){
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
	allocate(num); 
	
	int bound=M;
	int *range = get_hdf5_data(group,"range");//get range data

	//get map data 
	hsize_t dims_map[1];
	H5::DataSet dataset_map = group.openDataSet( "mapvec");
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
	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);

	for(int i=0;i<num_intervals;i++){
	  for(int j=it[i].first;j<=it[i].second;j++){
	    count_mem[0]=bound;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    dataset_map.read(base_ptr[j],H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	  }
	}      

	//reclaim memory
	delete [] it;
	delete [] data_domain;
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
    }
    
  template<int M> void MapVecRepI<M>::writehdf5(H5::Group group,entitySet& en) const {
    //entitySet en=domain();
      hsize_t dimf_domain[1];
      hsize_t dimf_range[1];
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
      int range=M,x=0;
      for(int i=0;i<num_intervals;i++){
	for(int j=it[i].first;j<=it[i].second;j++){
	    x=x+M;
	}
      }
      hsize_t dimf_map[1];
      dimf_map[0]=x;
      int *data_map = new int[x];
      int data_range[] = {range};
      dimf_range[0]=1;
      int bound=range;

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
      H5::DataSet dataset_map = group.createDataSet( "mapvec", H5::PredType::NATIVE_INT, dataspace_map );
      //get data using HDF5 hyperslab	
      H5::DataSpace dataspace_memory(RANK,dim_mem);

       x=0;//reset x
      for(int i=0;i<num_intervals;i++){
	for(int j=it[i].first;j<=it[i].second;j++){
	 count_mem[0]=range;
	  dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	  dataspace_map.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	  start_file[0]=start_file[0]+count_mem[0];//for next interval
	  dataset_map.write(base_ptr[j],H5::PredType::NATIVE_INT,dataspace_memory,dataspace_map);
	}
      }
      put_hdf5_data(group,data_domain,"domain",dimf_domain);//domain data
      put_hdf5_data(group,data_range,"range",dimf_range);//range data

      //reclaim memory
      delete [] it;
      delete [] data_domain;
      delete [] data_map;
    }

  template<int M> int* MapVecRepI<M>::get_hdf5_data(H5::Group group,const char* datasetname){
      H5::DataSet dataset = group.openDataSet( datasetname);
      H5::DataSpace dataspace = dataset.getSpace(); 
      hsize_t dims [1];
      dataspace.getSimpleExtentDims(dims , NULL);
      int *data  = new int[dims [0]];
      dataset.read(data , H5::PredType::NATIVE_INT );
      return data;
    }
    template<int M> void MapVecRepI<M>::put_hdf5_data(H5::Group group, int* data, const char* datasetname,hsize_t* dimf) const{
       int RANK=1;
       try{
	 H5::DataSpace dataspace( RANK, dimf );
	 H5::DataSet dataset = group.createDataSet( datasetname, H5::PredType::NATIVE_INT, dataspace );
	 dataset.write( data, H5::PredType::NATIVE_INT );
       }
       catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
       catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
       catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
     }

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
    MapVec(storeRepP &rp) { setRep(rp) ; }

    virtual ~MapVec() ;

    virtual void notification() ;
        
    MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    const entitySet domain() const { return Rep()->domain() ; }
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
    const_MapVec(storeRepP &rp) { setRep(rp) ; }

    virtual ~const_MapVec() ;

    virtual void notification() ;

    virtual instance_type access() const ;

    const_MapVec & operator=(const const_MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    const entitySet domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }
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
  
  template <int M> int MapVecRepI<M>::pack_size( const entitySet &e) {
    int size ;
    size = sizeof(int) * e.size() * M ;
    return size ;
  }

  template <int M> void MapVecRepI<M>::pack(void * ptr, int &loc, int &size, const entitySet &e ) {
    warn(true) ;
  }
  template <int M> void MapVecRepI<M>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    warn(true) ;
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
}

#endif

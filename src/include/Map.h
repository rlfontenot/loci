#ifndef MAP_H_
#define MAP_H_

#include "debug.h"
#include "Map_rep.h"
#include "store.h"
#include <algorithm>

namespace Loci {

  class MapRepI : public MapRep {
    entitySet store_domain ;
    int *alloc_pointer ;
    int *base_ptr ;
  public:
    MapRepI() { alloc_pointer = 0 ; base_ptr = 0 ; }
    MapRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~MapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual const entitySet &domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    int * get_base_ptr() const { return base_ptr ; }
  } ;
      
  class Map : public store_instance {
    friend class const_Map ;
    typedef MapRepI MapType ;
    int* base_ptr ;
  public:
    Map() { setRep(new MapType) ;}
    Map(Map &var) { setRep(var.Rep()) ; }
    Map(const entitySet &ptn) { setRep(new MapType(ptn)) ; }

    virtual ~Map() ;

    virtual void notification() ;

    Map & operator=(Map &str) { setRep(str.Rep()) ; return *this ;}
    Map & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    const entitySet &domain() const { return Rep()->domain() ; }
    operator storeRepP() { return Rep() ; }
    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    int &elem(int indx) { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const int &const_elem(int indx)  const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    int &operator[](int indx) { return elem(indx); }
    const int &operator[](int indx) const { return const_elem(indx) ; }
    //        operator int*() { return base_ptr; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;


  inline std::ostream & operator<<(std::ostream &s, const Map &m)
    { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, Map &m)
    { return m.Input(s) ; }

  class const_Map : public store_instance {
    typedef MapRepI MapType ;
    const int* base_ptr ;
  public:
    const_Map()
    { setRep(new MapType); }
    const_Map(const_Map &var) {setRep(var.Rep()) ; }
    const_Map(Map &var) {setRep(var.Rep()); }

    virtual ~const_Map() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_Map & operator=(Map &str) { setRep(str.Rep()) ; return *this ;}
    const_Map & operator=(const_Map &str)
    { setRep(str.Rep()) ; return *this ;}
    const_Map & operator=(storeRepP p) { setRep(p) ; return *this ;}

    const entitySet &domain() const { return Rep()->domain(); }
    operator storeRepP() { return Rep() ; }
    operator MapRepP()
    { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    const int &const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const int &operator[](int indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

        
  inline std::ostream & operator<<(std::ostream &s, const const_Map &m)
    { return m.Print(s) ; }


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
    virtual const entitySet &domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    VEC * get_base_ptr() const { return base_ptr ; }
  } ;

  template<int M> void MapVecRepI<M>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int size = ptn.Max()-top+1 ;
      alloc_pointer = new(VEC[size]) ;
      base_ptr = alloc_pointer-top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  template<int M> MapVecRepI<M>::~MapVecRepI<M>() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p) const {
    return new MapVecRepI<M>(p) ;
  }

  template<int M> const entitySet &MapVecRepI<M>::domain() const {
    return store_domain ;
  }

  template<int M> entitySet MapVecRepI<M>::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    FORALL(d,i) {
      for(int j=0;j<M;++j) 
        codomain += base_ptr[i][j] ;
    } ENDFORALL ;
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
      for(int j=0;j<M;++j)
        s << base_ptr[ii][j] << std::endl ;
    }ENDFORALL ;
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
    MapVec(MapVec<M> &var) { setRep(var.Rep()) ; }
    MapVec(const entitySet &ptn) { setRep(new MapVecType(ptn)) ; }

    virtual ~MapVec() ;

    virtual void notification() ;
        
    MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    const entitySet &domain() const { return Rep()->domain() ; }
    operator storeRepP() { return Rep() ; }
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
    const_MapVec(const_MapVec<M> &var) { setRep(var.Rep()) ; } 
    const_MapVec(MapVec<M> &var) { setRep(var.Rep()) ; }

    virtual ~const_MapVec() ;

    virtual void notification() ;

    virtual instance_type access() const ;

    const_MapVec & operator=(const const_MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    const entitySet &domain() const { return Rep()->domain() ; }
    operator storeRepP() { return Rep() ; }
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
    

  class multiMapRepI : public MapRep {
    entitySet store_domain ;
    int **index ;
    int *alloc_pointer ;
    int **base_ptr ;
  public:
    multiMapRepI() { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; }
    multiMapRepI(const store<int> &sizes) {
      index = 0 ;
      alloc_pointer = 0 ;
      base_ptr = 0 ;
      allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    virtual void allocate(const entitySet &ptn) ;
    virtual ~multiMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual const entitySet &domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    int ** get_base_ptr() const { return base_ptr ; }
    int *begin(int indx) { return base_ptr[indx] ; }
    int *end(int indx) { return base_ptr[indx+1] ; }
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
  } ;
      
  class multiMap : public store_instance {
    friend class const_multiMap ;
    typedef multiMapRepI MapType ;
    int **base_ptr ;
  public:
    multiMap() { setRep(new MapType) ; }
        
    multiMap(const store<int> &sizes) { setRep( new MapType(sizes) ); }
        
    multiMap(const multiMap &var) { setRep(var.Rep()) ; }

    multiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~multiMap() ;
    virtual void notification() ;

    multiMap & operator=(const multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    multiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; }

    const entitySet &domain() const { return Rep()->domain() ; }

    operator storeRepP() { return Rep() ; }
    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    int *elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const int *const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    int *operator[](int indx) { return elem(indx); }
    const int *operator[](int indx) const { return const_elem(indx) ; }
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    int *begin(int indx) { return base_ptr[indx] ; }
    int *end(int indx) { return base_ptr[indx+1] ; }
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;


  inline std::ostream & operator<<(std::ostream &s, const multiMap &m)
    { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, multiMap &m)
    { return m.Input(s) ; }

  class const_multiMap : public store_instance {
    typedef multiMapRepI MapType ;
    const int * const * base_ptr ;
  public:
    const_multiMap() { setRep(new MapType) ; }
        
    const_multiMap(const_multiMap &var) {  setRep(var.Rep()) ; }

    const_multiMap(multiMap &var) { setRep(var.Rep()) ; }

    virtual ~const_multiMap() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_multiMap & operator=(const_multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_multiMap & operator=(multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_multiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}

    const entitySet &domain() const { return Rep()->domain(); }
    operator storeRepP() { return Rep() ; }
    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    const int *const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const int *operator[](int indx) const { return const_elem(indx) ; }
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;


  inline std::ostream & operator<<(std::ostream &s, const const_multiMap &m)
    { return m.Print(s) ; }

  void inverseMap(multiMap &result,
                  const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;

  void inverseMap(multiMap &result,
                  const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;



  template<int M> multiMap MapVecRepI<M>::get_map()  {
    store<int> sizes(store_domain) ;
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

  template<int M> void inverseMap (multiMap &result,
                                   const MapVec<M> &input_map,
                                   const entitySet &input_image,
                                   const entitySet &input_preimage) {
    store<int> sizes(input_image) ;
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

#include "Map.h"

#include <stream.h>

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
      alloc_pointer = new(int[size]) ;
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

  const entitySet &MapRepI::domain() const {
    return store_domain ;
  }
    
  entitySet MapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    FORALL(d,i) {
      codomain += base_ptr[i] ;
    } ENDFORALL ;
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
    store<int> sizes(store_domain) ;
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
    if(ptn != EMPTY) {
      cerr << "Warning: multiMapRepI::allocate(const entitySet &) : "
           << endl ;
      cerr << "Generic allocation can not be applied to multiMaps"
           << endl ;
      cerr << "allocate set = " << ptn << endl ;
    }
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
    return new multiMapRepI(p)  ;
  }

  const entitySet &multiMapRepI::domain() const {
    return store_domain ;
  }
    
  entitySet multiMapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    FORALL(d,i) {
      for(const int *ip = begin(i);ip!=end(i);++ip)
        codomain += *ip ;
    } ENDFORALL ;
    return codomain ;
  }

  pair<entitySet,entitySet>
    multiMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domaini,domainu ;
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false ;
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
    store<int> sizes(e) ;
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
    store<int> sizes(input_image) ;

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
    store<int> sizes(input_image) ;
    
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

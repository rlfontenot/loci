#ifndef LOCI_MAP_H_
#define LOCI_MAP_H_

#include <Tools/debug.h>
#include <Map_rep.h>
#include <hdf5CC/H5cpp.h>

namespace Loci {

  entitySet image_section(const int *start, const int *end) ;

  class Map ;
  class multiMap ;
  
  class MapRepI : public MapRep {
    entitySet store_domain ;
    Entity *alloc_pointer ;
    Entity *base_ptr ;
  public:
    MapRepI() { alloc_pointer = 0 ; base_ptr = 0 ; }
    MapRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~MapRepI() ;
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
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(H5::Group group, entitySet &en) ;
    virtual void writehdf5(H5::Group group,entitySet& en) const ;
    int * get_base_ptr() const { return base_ptr ; }
  } ;
      
  class Map : public store_instance {
    friend class const_Map ;
    typedef MapRepI MapType ;
    Entity* base_ptr ;
  public:
    Map() { setRep(new MapType) ;}
    Map(const Map &var) { setRep(var.Rep()) ; }
    Map(storeRepP &rp) { setRep(rp) ; }

    virtual ~Map() ;

    virtual void notification() ;

    Map & operator=(const Map &str) { setRep(str.Rep()) ; return *this ;}
    Map & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    Entity &elem(Entity indx) { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const Entity &const_elem(Entity indx)  const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    Entity &operator[](Entity indx) { return elem(indx); }
    const Entity &operator[](Entity indx) const { return const_elem(indx) ; }
    //        operator int*() { return base_ptr; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const Map &m)
    { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, Map &m)
    { return m.Input(s) ; }

  class const_Map : public store_instance {
    typedef MapRepI MapType ;
    const Entity* base_ptr ;
  public:
    const_Map()
    { setRep(new MapType); }
    const_Map(const const_Map &var) {setRep(var.Rep()) ; }
    const_Map(const Map &var) {setRep(var.Rep()); }
    const_Map(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_Map() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_Map & operator=(const Map &str) { setRep(str.Rep()) ; return *this ;}
    const_Map & operator=(const const_Map &str)
    { setRep(str.Rep()) ; return *this ;}
    const_Map & operator=(storeRepP p) { setRep(p) ; return *this ;}

    entitySet domain() const { return Rep()->domain(); }
    operator MapRepP()
    { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    const Entity &const_elem(Entity indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const Entity &operator[](Entity indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

        
  inline std::ostream & operator<<(std::ostream &s, const const_Map &m)
  { return m.Print(s) ; }

  const int IMAGE_THRESHOLD = 4 ;
}


#endif

#ifndef MAPVEC_DEF_H
#define MAPVEC_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <store.h>
#include <multiMap.h>
#include <DMultiMap.h>
#include <hdf5_readwrite.h>
#include <distribute.h>

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
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq)  ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
    preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    VEC * get_base_ptr() const { return base_ptr; } 
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ;

  template<int M> class const_MapVec ;
  template<int M> class MapVec : public store_instance {
    friend  class const_MapVec<M> ;
    typedef MapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    VEC * base_ptr ;
    MapVec(const MapVec<M> &var) { setRep(var.Rep()) ; }
    MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    MapVec() { setRep(new MapVecType) ; }
    MapVec(storeRepP rp) { setRep(rp) ; }
    virtual ~MapVec() {}
    virtual void notification() ;
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

  template<int M> class const_MapVec : public store_instance {
    typedef MapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    const VEC * base_ptr ;
    const_MapVec(const const_MapVec<M> &var) { setRep(var.Rep()) ; } 
    const_MapVec(const MapVec<M> &var) { setRep(var.Rep()) ; }
    const_MapVec & operator=(const const_MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_MapVec & operator=(const MapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    const_MapVec() { setRep(new MapVecType) ; }
    const_MapVec(storeRepP rp) { setRep(rp) ; }
    virtual ~const_MapVec() {}
    virtual void notification() ;
    virtual instance_type access() const ;
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
  
} // end of namespace Loci

#endif

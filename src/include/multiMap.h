#ifndef MULTIMAP_H
#define MULTIMAP_H

#include <Tools/debug.h>
#include <Map_rep.h>

#include <Map.h>
#include <store.h>

namespace Loci {
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
    virtual storeRepP remap(const Map &m) const ;
    virtual void compose(const Map &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context)  ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
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
    virtual void readhdf5( hid_t group, entitySet &en) ;
    virtual void writehdf5( hid_t group,entitySet& en) const ;
    int ** get_base_ptr() const { return base_ptr ; }
    int *begin(int indx) { return base_ptr[indx] ; }
    int *end(int indx) { return base_ptr[indx+1] ; }
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
  private:
    int* get_hdf5_data(hid_t group, const char* datasetname) ;
    void put_hdf5_data(hid_t group, int* data,  const char* datasetname, 
                       hsize_t* dimf) const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual storeRepP thaw() ; 
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
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; }

    entitySet domain() const { return Rep()->domain() ; }

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
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
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
    
    const_multiMap(storeRepP rp) { setRep(rp) ; }
    
    virtual ~const_multiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_multiMap & operator=(const_multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_multiMap & operator=(multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_multiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }
    //    operator storeRepP() { return Rep() ; }
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


/*
  inline std::ostream & operator<<(std::ostream &s, const const_multiMap &m)
    { return m.Print(s) ; }
*/

  void inverseMap(multiMap &result,
                  const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(multiMap &result,
                  const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;



}

#endif

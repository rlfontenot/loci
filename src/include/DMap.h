#ifndef DMAP_H_
#define DMAP_H_

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <hdf5CC/H5cpp.h>

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif


namespace Loci {
  using std::hash_map ;

  entitySet image_section(const int *start, const int *end) ;
  class Map ;
  class multiMap ;
  class dMap;
  
  class dMapRepI : public MapRep {
    entitySet store_domain ;
    hash_map<int,int> attrib_data;
  public:
    dMapRepI() { }
    dMapRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dMapRepI() ;
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
    virtual void readhdf5(H5::Group group, entitySet &user_eset) ;
    virtual void writehdf5(H5::Group group,entitySet& en) const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    hash_map<int,int> *get_attrib_data() { return &attrib_data; }
  } ;
      
  class dMap : public store_instance {
    friend class const_dMap ;
    typedef dMapRepI MapType ;
    hash_map<int,int> *attrib_data;
  public:
    dMap() { setRep(new MapType) ;}
    dMap(const dMap &var) { setRep(var.Rep()) ; }
    dMap(storeRepP &rp) { setRep(rp) ; }

    virtual ~dMap() ;

    virtual void notification() ;
    
    dMap & operator=(const dMap &str) { setRep(str.Rep()) ; return *this ;}
    dMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }

    int &elem(int indx) { 
      return (*attrib_data)[indx]; 
    }

    const int &const_elem(int indx)  const 
    { 
      return (*attrib_data)[indx]; 
    }

    int &operator[](int indx) { return elem(indx); }

    const int &operator[](int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const dMap &m)
    { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, dMap &m)
    { return m.Input(s) ; }

  class const_dMap : public store_instance {
    typedef dMapRepI MapType ;
    hash_map<int,int>  *attrib_data;
  public:
    const_dMap()
    { setRep(new MapType); }
    const_dMap(const const_dMap &var) {setRep(var.Rep()) ; }
    const_dMap(const dMap &var) {setRep(var.Rep()); }
    const_dMap(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dMap() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_dMap & operator=(const dMap &str) { setRep(str.Rep()) ; return *this ;}

    const_dMap & operator=(const const_dMap &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dMap & operator=(storeRepP p) { setRep(p) ; return *this ;}

    entitySet domain() const { return Rep()->domain(); }

    operator MapRepP()
    { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    const int &const_elem(int indx)  const {
      return (*attrib_data)[indx]; 
    }

    const int &operator[](int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

        
  inline std::ostream & operator<<(std::ostream &s, const const_dMap &m)
    { return m.Print(s) ; }
 
  void inverseMap(multiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage);
}

#endif

#ifndef DMAP_H_
#define DMAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>

#include<Tools/hash_map.h>


namespace Loci {

  entitySet image_section(const int *start, const int *end) ;
  class Map ;
  class multiMap ;
  class dMap;
  
  class dMapRepI : public MapRep {
    HASH_MAP(int,int) attrib_data;
  public:
    dMapRepI() { }
    dMapRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
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
    virtual void readhdf5(hid_t group, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group,entitySet& en) const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual storeRepP thaw() ;
    virtual HASH_MAP(int,int) *get_attrib_data() { return &attrib_data; }
  } ;
      
  class dMap : public store_instance {
    friend class const_dMap ;
    typedef dMapRepI MapType ;
    HASH_MAP(int,int) *attrib_data;
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
    HASH_MAP(int,int)  *attrib_data;
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

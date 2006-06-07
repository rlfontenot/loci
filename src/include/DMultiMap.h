#ifndef DMULTIMAP_H
#define DMULTIMAP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <DMap.h>
#include <vector>

#include <Tools/hash_map.h>

#include <Map.h>
#include <store.h>
namespace Loci {

  class dmultiMapRepI : public MapRep {
  HASH_MAP(int, std::vector<int> )  attrib_data;
  public:
  dmultiMapRepI() { }
    dmultiMapRepI(const store<int> &sizes) { allocate(sizes) ; }
  void allocate(const store<int> &sizes) ;
  virtual void allocate(const entitySet &ptn) ;
  virtual ~dmultiMapRepI() ;
  virtual storeRep *new_store(const entitySet &p) const ;
  virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
  virtual storeRepP remap(const dMap &m) const ;
  virtual void compose(const dMap &m, const entitySet &context) ;
  virtual void copy(storeRepP &st, const entitySet &context)  ;
  virtual void gather(const dMap &m, storeRepP &st,
                      const entitySet &context)  ;
  virtual void scatter(const dMap &m, storeRepP &st,
                       const entitySet &context) ;
  
  virtual int pack_size(const entitySet &e) ;
  virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
  virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
  
  virtual entitySet domain() const ;
  
  virtual entitySet image(const entitySet &domain) const ;
  virtual std::pair<entitySet,entitySet>
  preimage(const entitySet &codomain) const ;
  virtual storeRepP get_map() ;
  virtual std::ostream &Print(std::ostream &s) const ;
  virtual std::istream &Input(std::istream &s) ;
  virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
  virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
  virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  virtual storeRepP freeze() ;
  virtual storeRepP thaw() ;
  HASH_MAP(int,std::vector<int> ) *get_attrib_data() {return &attrib_data;}
  virtual DatatypeP getType() ;
  virtual frame_info read_frame_info(hid_t group_id) ;
  virtual frame_info write_frame_info(hid_t group_id) ;
  } ;
  
  //***************************************************************************
  
  class dmultiMap : public store_instance {
    friend class const_dmultiMap ;
    typedef dmultiMapRepI MapType ;
    HASH_MAP(int, std::vector<int> )   *attrib_data;
    dmultiMap(const dmultiMap &var) { setRep(var.Rep()) ; }
    dmultiMap & operator=(const dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    dmultiMap() { setRep(new MapType) ; }
        
    dmultiMap(const store<int> &sizes) { setRep( new MapType(sizes) ); }
        

    dmultiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~dmultiMap() ;

    virtual void notification() ;

    dmultiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; 
    }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; 
    }

    std::vector<int> &elem(int indx) {
      return (*attrib_data)[indx]; 
    }

    const std::vector<int> &const_elem(int indx)  const 
    {
      HASH_MAP(int, std::vector<int> ) :: const_iterator   ci;
     
      ci = attrib_data->find(indx);
      if( ci == attrib_data->end()) {
          std::cerr << "Error: out of range entity " << std::endl;
          exit(0);
      }
      
      return ci->second ;
    }

    std::vector<int> &operator[](int indx) { return elem(indx); }
    
    const std::vector<int> &operator[](int indx) const 
      { return const_elem(indx) ; }
    
    int num_elems(int indx) const 
      {
	HASH_MAP(int, std::vector<int> ) :: const_iterator   ci;
	std::vector<int>     newVec;
	
	ci = attrib_data->find(indx);
	if( ci != attrib_data->end())
          return( (ci->second).size() );
	
	return(0);
	
      }

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }

    std::istream &Input(std::istream &s) 
    { return Rep()->Input(s) ; }
  } ;
  
  //***************************************************************************

  inline std::ostream & operator<<(std::ostream &s, const dmultiMap &m)
    { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, dmultiMap &m)
    { return m.Input(s) ; }

  //***************************************************************************

  class const_dmultiMap : public store_instance {
    typedef dmultiMapRepI      MapType ;
    HASH_MAP(int,std::vector<int> ) *attrib_data;
    const_dmultiMap(const const_dmultiMap &var) {  setRep(var.Rep()) ; }
    const_dmultiMap(const dmultiMap &var) { setRep(var.Rep()) ; }
    const_dmultiMap & operator=(const const_dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dmultiMap & operator=(const dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    const_dmultiMap() { setRep(new MapType) ; }
    
    
    
    const_dmultiMap(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dmultiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_dmultiMap & operator=(storeRepP p) 
    { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }

    operator MapRepP() 
    {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; 
    }

    const std::vector<int> &const_elem(int indx)  const {
      HASH_MAP(int,std::vector<int> ) :: const_iterator   ci;
      ci = attrib_data->find(indx);
      if( ci == attrib_data->end() ) {
	cerr << " Trying to access the dmultiMap out of bounds " << endl ;
	exit(0) ;
      }
      return( ci->second );
    }

    const std::vector<int> operator[](int indx) const 
      { return const_elem(indx) ; }
    
    int num_elems(int indx) const 
      {
       HASH_MAP(int,std::vector<int> ) :: const_iterator   ci;
       ci = attrib_data->find(indx);
       if( ci != attrib_data->end() )
	 return( (ci->second).size() );
       return(0);
      }

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }
  } ;

  //***************************************************************************

  inline std::ostream & operator<< (std::ostream &s, 
                                    const const_dmultiMap &m)
  { return m.Print(s) ; }

  //***************************************************************************

  void inverseMap(dmultiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const Map &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const const_dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const const_Map &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const dmultiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  
}

#endif

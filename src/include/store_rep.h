#ifndef STORE_REP_H
#define STORE_REP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/nptr.h>
#include <entitySet.h>
#include <istream>
#include <ostream>

#include <data_traits.h>


namespace Loci {
  enum store_type { STORE, PARAMETER, MAP, CONSTRAINT, BLACKBOX } ;

  class Map ;
  class dMap ;
  class storeRep ;

  typedef NPTR<storeRep> storeRepP ;
  typedef const_NPTR<storeRep> const_storeRepP ;
  
  struct frame_info {
    int is_stat ; 
    int size ;
    std::vector<int> first_level ;
    std::vector<int> second_level ;
    frame_info() {
      is_stat = 0 ;
      size = 0 ;
    }
    frame_info(int a , int b) {
      is_stat = a ;
      size = b ;
      
    }
    frame_info(int a , int b, std::vector<int> c, std::vector<int> d) {
      is_stat = a ;
      size = b ;
      
      first_level = c ;
      second_level = d ;
    } 
    frame_info(const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size)
	first_level = fi.first_level ;
      if(is_stat) 
	second_level = fi.second_level ;
    }
    
    frame_info &operator = (const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size) 
	first_level = fi.first_level ;
      if(is_stat)
	second_level = fi.second_level ;
      
      return *this ;
    }
  } ;
  class storeRep : public NPTR_type {
  public:
    virtual ~storeRep() ;
    virtual void allocate(const entitySet &p) = 0 ;
    virtual void shift(int_type offset) = 0 ;
    virtual void set_elem_size(int sz) ;
    virtual storeRep *new_store(const entitySet &p) const = 0 ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const = 0 ;
    virtual storeRepP remap(const dMap &m) const = 0 ;
    virtual void copy(storeRepP &st, const entitySet &context) = 0 ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) = 0 ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) = 0 ;
    virtual int pack_size(const entitySet &e) = 0;
    virtual void pack(void *ptr, int &loc, int &size,  const entitySet &e) = 0 ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) = 0 ;
    virtual store_type RepType() const = 0 ;
    virtual std::ostream &Print(std::ostream &s) const = 0 ;
    virtual std::istream &Input(std::istream &s) = 0 ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) = 0;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &en) const = 0;
    virtual entitySet domain() const = 0 ;
    virtual storeRepP getRep() ;
    virtual storeRepP getRep() const ;
    virtual  DatatypeP getType() = 0 ;
    virtual frame_info read_frame_info(hid_t group_id) = 0 ;
    virtual frame_info write_frame_info(hid_t group_id) = 0 ;
  } ;


  class store_instance : public eventNotify {
    friend class store_ref ;
    storeRepP rep ;
  public:
    enum instance_type { READ_WRITE, READ_ONLY } ;
    void setRep(const storeRepP &p)
      { rep = p ; rep.set_notify(this); notification() ; }
    storeRepP Rep() { return rep->getRep(); }
    storeRepP Rep() const { return rep->getRep(); }
    virtual instance_type access() const ;
  } ;
  
  class store_ref : public store_instance, public storeRep {
  public:
    store_ref() {}
    store_ref(storeRepP &p) { setRep(p); }
    virtual ~store_ref() ;
    
    store_ref &operator=(storeRepP &p) {
      setRep(p) ;
      return *this ;
    }
        
    store_ref &operator=(const store_ref &sr) {
      setRep(sr.rep) ;
      return *this ;
    }

    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) {
      Rep()->readhdf5(group_id, dataspace, dataset, dimension, name, fi, en);
    };
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const {
      Rep()->writehdf5(group_id, dataspace, dataset, dimension, name, en);
    };
    virtual entitySet domain() const ;
    virtual storeRepP getRep() ;
    virtual storeRepP getRep() const ;
    virtual void notification() ;
    virtual DatatypeP getType() {
     return  Rep()->getType() ;
    }
    virtual frame_info read_frame_info(hid_t group_id) {
      return Rep()->read_frame_info(group_id) ;
    }
    virtual frame_info write_frame_info(hid_t group_id) {
      return write_frame_info(group_id) ;
    }
  } ;
  typedef NPTR<store_ref> store_refP ;
    
}

#endif

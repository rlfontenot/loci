#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

namespace Loci {
    
  class constraintRep : public storeRep {
    entitySet constraint_set ;
  public:
    constraintRep() ;
    constraintRep(const entitySet &p) ;
    virtual ~constraintRep() ;
    virtual void allocate(const entitySet &p) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    entitySet *get_constraint() { return &constraint_set ; }
    virtual DatatypeP getType() ;
    virtual frame_info read_frame_info(hid_t group_id) ;
    virtual frame_info write_frame_info(hid_t group_id) ;
  } ;

  class constraint : public store_instance {
    typedef constraintRep constraintType ;
      entitySet *data ;
  public:
    constraint() ;
    constraint(const constraint &var) ;
    constraint(const storeRepP &rp) { setRep(rp) ; }
    virtual ~constraint() ;

    constraint & operator=(const constraint &p)
    { *data = *(p.data) ; return *this ;}
    constraint & operator=(const storeRepP &p)
    { setRep(p) ;  return *this ;}
    constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    entitySet &operator*() { return *data ; }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, constraint &t)
    { return t.Input(s) ; }


  // a duplicate of the constraint class
  class Constraint : public store_instance {
    typedef constraintRep constraintType ;
    entitySet *data ;
  public:
    Constraint() ;
    Constraint(const Constraint &var) ;
    Constraint(const storeRepP &rp) { setRep(rp) ; }
    virtual ~Constraint() ;

    Constraint & operator=(const Constraint &p)
    { setRep(p.Rep()) ; return *this ;}
    Constraint & operator=(const storeRepP &p)
    { setRep(p) ;  return *this ;}
    Constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    entitySet &operator*() { return *data ; }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const Constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, Constraint &t)
    { return t.Input(s) ; }
  
}

#endif
    

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

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
    virtual storeRep *new_store(const entitySet &p) const ;
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
    virtual void readhdf5(  hid_t group, entitySet &en) ;
    virtual void writehdf5( hid_t group,entitySet& en) const ;
    entitySet *get_constraint() { return &constraint_set ; }
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
    { setRep(p.Rep()) ; return *this ;}
    constraint & operator=(const storeRepP &p)
    { setRep(p) ;  return *this ;}
    constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    //    entitySet * operator&() { return data ; }
    //    const entitySet * operator &() const { return data ; }
    entitySet &operator*() { return *data ; }
    
    //    operator storeRepP() { return Rep() ; }

    //    operator entitySet() { return *data ; }
    //    operator entitySet() const { return *data ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, constraint &t)
    { return t.Input(s) ; }

}

#endif
    

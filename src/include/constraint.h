#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <Tools/debug.h>
#include <store_rep.h>
#include <Tools/stream.h>
#include <hdf5CC/H5cpp.h>

namespace Loci {
    
  class constraintRep : public storeRep {
    entitySet constraint ;
  public:
    constraintRep() ;
    constraintRep(const entitySet &p) ;
    virtual ~constraintRep() ;
    virtual void allocate(const entitySet &p) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual store_type RepType() const ;
    virtual const entitySet &domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    entitySet *get_constraint() { return &constraint ; }
  } ;

  class constraint : public store_instance {
    typedef constraintRep constraintType ;
    entitySet *data ;
  public:
    constraint() ;
    constraint(constraint &var) ;
    constraint(const entitySet &ptn) ;
    constraint(storeRepP &rp) { setRep(rp) ; }
    virtual ~constraint() ;

    constraint & operator=(constraint &p)
    { setRep(p.Rep()) ; return *this ;}
    constraint & operator=(storeRepP p)
    { setRep(p) ;  return *this ;}
    constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    entitySet * operator&() { return data ; }
    const entitySet * operator &() const { return data ; }

    entitySet &operator*() { return *data ; }
    
    operator storeRepP() { return Rep() ; }

    operator entitySet() { return *data ; }
    operator entitySet() const { return *data ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, constraint &t)
    { return t.Input(s) ; }

}

#endif
    

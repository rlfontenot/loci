#ifndef ENTITYSET_H
#define ENTITYSET_H 1

#include <Tools/intervalSet.h>


namespace Loci {

  typedef int int_type ;

#ifdef ENTITY

  class Entity {
  private:
    int id;

  public:
    Entity() {}
    Entity(int i) {id = i;}
    Entity(const Entity &e) {id = e.id;} /* copy constructor */
    int getid() const {return id;}
    std::ostream & Print (std::ostream &s) const { 
      s << id ;
      return s;}
    std::istream & Read (std::istream &s) {
      s >> id;
      return s;
    }
    bool operator==(const Entity &e) const {
      return (id == e.getid());
    }
    operator int () const {return getid();}
  };
  
  inline std::ostream &operator<<(std::ostream &s, const Entity &e) {
    return e.Print(s);
  }
  
  inline std::istream &operator>>(std::istream &s, Entity &e) {
    return e.Read(s);
  }
  
  
  inline int getEntityIdentity(const Entity e){return e.getid();}  

  typedef intervalSet entitySet ;
 
#else
  typedef Loci::int_type Entity ;
  inline int getEntityIdentity(const int i){return i;}
  typedef intervalSet entitySet ;
  template<class T> inline entitySet create_entitySet(T start,T end) {
      return create_intervalSet(start,end) ;
  }
#endif

  template<class Op> inline void do_loop(const entitySet &iset, Op f) {
    for(int i=0;i<iset.num_entityIntervals();++i) {
      const Loci::int_type stop = iset[i].second ;
      for(Loci::int_type indx=iset[i].first;indx<=stop;++indx)
        f(indx) ;
    }
  }

  template<class Op> inline void do_loop(const sequence &seq, Op f) {
    for(int i=0;i<seq.num_entityIntervals();++i) {
      const bool dir = seq[i].first>seq[i].second?false:true ;
      const Loci::int_type stop = seq[i].second + (dir?1:-1) ;
      for(Loci::int_type indx=seq[i].first;indx!=stop;dir?++indx:--indx) 
        f(indx) ;
    }
  }
  
  template<class T> inline void do_loop(const sequence &seq, T *cp,
                                        void(T::*pmf)(Entity) ) {
    for(int i=0;i<seq.num_entityIntervals();++i) {
      const bool dir = seq[i].first>seq[i].second?false:true ;
      const Loci::int_type stop = seq[i].second + (dir?1:-1) ;
      for(Loci::int_type indx=seq[i].first;indx!=stop;dir?++indx:--indx) 
        (cp->*pmf)(indx) ;
    }
  }

  template<class T> inline void do_loop(const sequence &seq, T *cp) {
    for(int i=0;i<seq.num_entityIntervals();++i) {
      const bool dir = seq[i].first>seq[i].second?false:true ;
      const Loci::int_type stop = seq[i].second + (dir?1:-1) ;
      for(Loci::int_type indx=seq[i].first;indx!=stop;dir?++indx:--indx) 
        cp->calculate(indx) ;
    }
  }


}

#endif

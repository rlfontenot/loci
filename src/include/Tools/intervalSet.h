#ifndef INTERVALSET_H
#define INTERVALSET_H 1

#include <Tools/debug.h>
#include <Tools/Handle.h>
#include <Tools/tools.h>

#include <limits>
#include <vector>
#include <iostream>

namespace Loci {

  typedef int int_type ;

  const int_type UNIVERSE_MAX = std::numeric_limits<int_type>::max() - 1 ;
  const int_type UNIVERSE_MIN = std::numeric_limits<int_type>::min() + 1 ;

  typedef std::pair<int_type,int_type> interval ;

  std::ostream & operator<<(std::ostream &s, const interval &i) ;
  std::istream & operator>>(std::istream &s, interval &i) ;

  class intervalSet ;
  class sequence ;
    
  extern const intervalSet EMPTY ;

  class intervalSet {
  public:
    class intervalSetIterator ;
    typedef std::vector<interval> intervalSetRep ;
  private:

    friend class sequence ;
    

    Handle<intervalSetRep> Rep ;

    struct rep_holder {
      Handle<intervalSetRep> Rep ;
    } ;
    
    static rep_holder *rhp ;
    static Handle<intervalSetRep> & getEmptyRep() {
      if(rhp==0)
        rhp = new rep_holder ;
      return rhp->Rep ;
    } ;
  public:
    class intervalSetIterator {
      intervalSetRep::const_iterator current_interval,end_interval ;
      int_type current_value ;
      void increment() {
        if(++current_value>current_interval->second) {
          ++current_interval ;
          if(current_interval!=end_interval)
            current_value = current_interval->first ;
          else
            current_value = UNIVERSE_MAX ;
        }
      }                
    public:
      intervalSetIterator() {}
        //      {current_interval=0;end_interval=0;current_value=0;}
      intervalSetIterator(intervalSetRep::const_iterator r,
                          intervalSetRep::const_iterator end,int_type e) :
        current_interval(r),end_interval(end),current_value(e) {}
      int operator*() const { return current_value; }
      intervalSetIterator &operator++() {
        increment() ;
        return *this ;
      }
      intervalSetIterator operator++(int ) {
        intervalSetIterator tmp = *this ;
        increment() ;
        return tmp ;
      }
      bool operator==(const intervalSetIterator &i) const {
        return current_value == i.current_value ; }
      bool operator!=(const intervalSetIterator &i) const {
        return current_value != i.current_value ; }
    } ;


    intervalSet() : Rep(getEmptyRep()) {} ;
    intervalSet(const interval &ivl) { interval i(min(ivl.first,ivl.second),
                                                  max(ivl.first,ivl.second));
    Rep->push_back(i) ; }
    intervalSet(const intervalSet &ptn): Rep(ptn.Rep) {}
    explicit intervalSet(const sequence &seq) ;
    ~intervalSet() {}
    
    typedef intervalSetIterator const_iterator ;

    const_iterator begin() const {
      if(Rep->size() != 0)
        return intervalSetIterator(Rep->begin(),Rep->end(),(*Rep)[0].first) ;
      else
        return end() ;
    }
    const_iterator end() const {
      return intervalSetIterator(Rep->end(),Rep->end(),UNIVERSE_MAX) ;
    }
    

    bool inSet(int_type indx) const ;

    intervalSet & operator=(const intervalSet &ptn)  
    { Rep = ptn.Rep ; return *this;}
    intervalSet & operator=(const interval &ivl) 
    { Rep.New() ; Rep->push_back(ivl) ; return *this ; }

    bool Equal(const intervalSet &ptn) const
    { return *(Rep) == *(ptn.Rep) ; }

    bool less_than(const intervalSet &ptn) const
    { return (*Rep) < *(ptn.Rep) ;}
    bool greater_than(const intervalSet &ptn) const
    {return (*Rep) > *(ptn.Rep) ;}
    
    intervalSet & operator>>=(int rotval) ;
    intervalSet & operator<<=(int rotval) ;

    intervalSet operator<<(int rotval) const ;
    intervalSet operator>>(int rotval) const ;

    int size() const {
      int size = 0 ;
      std::vector<interval>::const_iterator i ;
      for(i = Rep->begin();i!= Rep->end();++i)
        size += i->second-i->first + 1 ;
      return size ;  }

    int num_intervals() const { return Rep->size() ; }
    const interval & operator[](int_type indx) const 
    { fatal(indx<0); fatal(indx>=num_intervals()) ;
      return (*Rep)[indx]; }

    void Union(const interval &ivl) ;
    void Union(const intervalSet &ptn) ;
    static intervalSet Union(const intervalSet &set1,
                             const intervalSet &set2) ;

    void Intersection(const interval &ivl) ;
    void Intersection(const intervalSet &ptn) ;
    static intervalSet Intersection(const intervalSet &set1,
                                    const intervalSet &set2) ;
    void Complement() ;
    static intervalSet Complement(const intervalSet &set) ;
        
    std::ostream &Print(std::ostream &s) const ;
    void Print() const { Print(std::cout) ; }

    std::istream &Input(std::istream &s) ;

    int Min() const {
      int sz = Rep->size();
      return sz==0?UNIVERSE_MAX:(*Rep)[0].first ;}
    int Max() const {
      int sz = Rep->size() ;
      return sz==0?UNIVERSE_MIN:(*Rep)[sz-1].second ;}
  } ;

  inline intervalSet operator~(const intervalSet &e) {
    return intervalSet::Complement(e) ;
  }

  inline intervalSet & operator+=(intervalSet &e, int_type ival) 
    { e.Union(interval(ival,ival)) ; return e ; }

  inline intervalSet & operator+=(intervalSet &e, const interval &ivl)
    { e.Union(ivl) ; return e; }

  inline intervalSet & operator+=(intervalSet &e, const intervalSet &ptn)
    { e.Union(ptn) ; return e ; }

  inline intervalSet & operator|=(intervalSet &e, int_type ival)
    { e.Union(interval(ival,ival)) ; return e ; }

  inline intervalSet & operator|=(intervalSet &e, const interval &ivl)
    { e.Union(ivl) ; return e ; }

  inline intervalSet & operator|=(intervalSet &e, const intervalSet &ptn)
    { e.Union(ptn) ; return e ; }

  inline intervalSet & operator&=(intervalSet &e, int_type ival)
    { e.Intersection(interval(ival,ival)) ;
      return e; }

  inline intervalSet & operator&=(intervalSet &e, const interval &ivl)
    { e.Intersection(ivl) ; return e ; }

  inline intervalSet & operator&=(intervalSet &e, const intervalSet &ptn)
    { e.Intersection(ptn) ; return e ; }

  inline intervalSet operator+(const intervalSet & e, int_type ival) {
    intervalSet retp = e ;
    retp.Union(interval(ival,ival)) ;
    return retp ;
  }

  inline intervalSet operator+(const int &ival, const intervalSet &e) {
    intervalSet retp = e ;
    retp.Union(interval(ival,ival)) ;
    return retp ;
  }

  inline intervalSet operator+(const intervalSet &e, const interval &ivl) {
    intervalSet retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  inline intervalSet operator+(const interval &ivl, const intervalSet &e) {
    intervalSet retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  inline intervalSet operator+(const intervalSet &e, const intervalSet &ptn) {
    intervalSet retp = e ;
    retp.Union(ptn) ;
    return retp ;
  }

  inline intervalSet operator|(const intervalSet &e, int_type ival) {
    intervalSet retp = e ;
    retp.Union(interval(ival,ival)) ;
    return retp ;
  }

  inline intervalSet operator|(int_type ival, const intervalSet &e) {
    intervalSet retp = e ;
    retp.Union(interval(ival,ival)) ;
    return retp ;
  }

  inline intervalSet operator|(const intervalSet &e, const interval &ivl) {
    intervalSet retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  inline intervalSet operator|(const interval &ivl,const intervalSet &e) {
    intervalSet retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  inline intervalSet operator|(const intervalSet &e, const intervalSet &ptn) {
    intervalSet retp = e ;
    retp.Union(ptn) ;
    return retp ;
  }

  inline intervalSet operator&(const intervalSet &e, int_type ival)  {
    return intervalSet::Intersection(e,intervalSet(interval(ival,ival))) ;
  }

  inline intervalSet operator&(int_type ival,const intervalSet &e)  {
    return intervalSet::Intersection(e,intervalSet(interval(ival,ival))) ;
  }

  inline intervalSet operator&(const intervalSet &e, const interval &ivl) {
    return intervalSet::Intersection(e,intervalSet(ivl)) ;
  }

  inline intervalSet operator&(const interval &ivl,const intervalSet &e) {
    return intervalSet::Intersection(e,intervalSet(ivl)) ;
  }

  inline intervalSet operator&(const intervalSet &e, const intervalSet &ptn) {
    return intervalSet::Intersection(e,ptn) ;
  }

  inline intervalSet operator-(const intervalSet &e, int_type &ival) {
    return intervalSet::Intersection(e,~intervalSet(interval(ival,ival))) ;
  }

  inline intervalSet operator-(const intervalSet &e, const interval &ivl) {
    return intervalSet::Intersection(e,~intervalSet(ivl)) ;
  }

  inline intervalSet operator-(const interval &ivl, const intervalSet &e) {
    return intervalSet::Intersection(intervalSet(ivl),~e) ;
  }


  inline intervalSet operator-(const intervalSet &e, const intervalSet &ptn) {
    return intervalSet::Intersection(e,~ptn) ;
  }

  inline intervalSet operator^(const intervalSet &e, const intervalSet &ptn) {
    intervalSet retp = ((e|ptn) & ~(e & ptn)) ;
    return retp ;
  }
  inline intervalSet operator^(const intervalSet &e, int_type ival) {
    intervalSet ivlp = interval(ival,ival) ;
    intervalSet retp = e^ivlp ;
    return retp ;
  }

  inline intervalSet operator^(int_type ival, const intervalSet &e) {
    intervalSet ivlp = interval(ival,ival) ;
    intervalSet retp = e^ivlp ;
    return retp ;
  }

  inline intervalSet operator^(const intervalSet &e, const interval &ivl) {
    intervalSet ivlp = ivl ;
    intervalSet retp = e^ivlp ;
    return retp ;
  }

  inline intervalSet operator^(const interval &ivl, const intervalSet &e) {
    intervalSet ivlp = ivl ;
    intervalSet retp = e^ivlp ;
    return retp ;
  }


  inline bool operator==(const intervalSet &es1, const intervalSet &es2) {
    return es1.Equal(es2) ;
  }

  inline bool operator!=(const intervalSet &es1, const intervalSet &es2) {
    return !es1.Equal(es2) ;
  }

  inline bool operator<(const intervalSet &es1, const intervalSet &es2) {
    return es1.less_than(es2) ;
  }

  inline bool operator>(const intervalSet &es1, const intervalSet &es2) {
    return es1.greater_than(es2) ;
  }

  inline bool operator<=(const intervalSet &es1, const intervalSet &es2) {
    return !(es1.greater_than(es2)) ;
  }

  inline bool operator>=(const intervalSet &es1, const intervalSet &es2) {
    return !(es1.less_than(es2)) ;
  }

  inline intervalSet intervalSet::operator>>(int rotval) const  {
    intervalSet retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  inline intervalSet intervalSet::operator<<(int rotval) const {
    intervalSet retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }


  inline intervalSet & operator^=(intervalSet &e, const intervalSet &ptn) {
    e =  e^ptn ;
    return e ;
  }
    
  inline intervalSet & operator^=(intervalSet &e, const interval &ivl) {
    intervalSet ivlp ;
    ivlp = ivl ;
    e ^= ivlp ;
    return e ;
  }

  inline intervalSet & operator^=(intervalSet &e, int_type ival) {
    intervalSet ivlp ;
    ivlp = interval(ival,ival) ;
    e ^= ivlp ;
    return e ;
  }

  inline intervalSet & operator-=(intervalSet &e, const intervalSet &ptn) {
    intervalSet ptnc = intervalSet::Complement(ptn);
    e &= ptnc ;
    return e ;
  }
    
  inline intervalSet & operator-=(intervalSet &e, const interval &ivl) {
    intervalSet ptnc ;
    ptnc = ivl ;
    ptnc.Complement() ;
    e &= ptnc ;
    return e ;
  }

  inline intervalSet & operator-=(intervalSet &e, int_type ival) {
    intervalSet ptnc;
    ptnc = interval(ival,ival) ;
    ptnc.Complement() ;
    e &= ptnc ;
    return e ;
  }

    
  inline std::ostream & operator<<(std::ostream &s, const intervalSet &e)
    { return e.Print(s) ; }
    
  inline std::istream & operator>>(std::istream &s, intervalSet &e)
    { return e.Input(s) ;}
    
    
  class sequence {

  public:
    typedef intervalSet::intervalSetRep sequenceRep ;
  private:
    
    Handle<sequenceRep> Rep ;

    struct rep_holder {
      Handle<sequenceRep> Rep ;
    } ;
    static rep_holder *rhp ;
    static Handle<sequenceRep> & getEmptyRep() {
      if(rhp==0)
        rhp = new rep_holder ;
      return rhp->Rep ;
    } ;
  public:
    class sequenceIterator {
      sequenceRep::const_iterator current_interval,end_interval ;
      int_type current_value ;
      int dir ;
      void increment() {
        if(current_value == current_interval->second) {
          ++current_interval ;
          if(current_interval!=end_interval) {
            current_value = current_interval->first ;
            dir = current_interval->second>=current_value?1:-1 ;
          } else {
            current_value = UNIVERSE_MAX ;
            dir = 0 ;
          }
        } else
          current_value += dir ;
      }                
    public:
      sequenceIterator() {}
        //      { current_interval=0;end_interval=0;current_value=0; }
      sequenceIterator(sequenceRep::const_iterator r,
                       sequenceRep::const_iterator end,int_type e) :
        current_interval(r),end_interval(end),current_value(e)
      { dir = ((e!=UNIVERSE_MAX)&&(r->second>=r->first))?1:-1 ;}
      int operator*() const { return current_value; }
      sequenceIterator &operator++() {
        increment() ;
        return *this ;
      }
      sequenceIterator operator++(int ) {
        sequenceIterator tmp = *this ;
        increment() ;
        return tmp ;
      }
      bool operator==(const sequenceIterator &i) const {
        return current_value == i.current_value  &&
          current_interval == i.current_interval; }
      bool operator!=(const sequenceIterator &i) const 
      { return !operator==(i) ; }
    } ;

    
    sequence() : Rep(getEmptyRep()) {} ;
    sequence(const sequence &seq): Rep(seq.Rep) {}
    sequence(const interval &ivl) {
      interval i(min(ivl.first,ivl.second),
                 max(ivl.first,ivl.second));
      Rep->push_back(i) ;
    }
    sequence(const intervalSet &ptn) : Rep(ptn.Rep) {}
    ~sequence() {}
    
    typedef sequenceIterator const_iterator ;

    const_iterator begin() const {
      if(Rep->size() != 0)
        return sequenceIterator(Rep->begin(),Rep->end(),(*Rep)[0].first) ;
      return end() ;
    }
    const_iterator end() const {
      return sequenceIterator(Rep->end(),Rep->end(),UNIVERSE_MAX) ;
    }
    

    sequence & operator=(const sequence &ptn)  
    { Rep = ptn.Rep ; return *this;}
    sequence & operator=(const interval &ivl) 
    { Rep.New() ; Rep->push_back(ivl) ; return *this ; }

    bool Equal(const sequence &ptn) const
    { return *(Rep) == *(ptn.Rep) ; }

    bool less_than(const sequence &ptn) const
    { return (*Rep) < *(ptn.Rep) ;}
    bool greater_than(const sequence &ptn) const
    {return (*Rep) > *(ptn.Rep) ;}
    
    sequence & operator>>=(int rotval) ;
    sequence & operator<<=(int rotval) ;

    sequence operator<<(int rotval) const ;
    sequence operator>>(int rotval) const ;

    int size() const {
      int size = 0 ;
      std::vector<interval>::const_iterator i ;
      for(i = Rep->begin();i!= Rep->end();++i)
        size += abs(i->first-i->second) + 1 ;
      return size ;  }
    int num_intervals() const { return Rep->size();  }
    const interval & operator[](int_type indx) const 
    { fatal(indx<0); fatal(indx>=num_intervals()) ;
      return (*Rep)[indx]; }

    void Append(const interval &ivl) ;
    void Append(const sequence &seq) ;
    void Append(const intervalSet &ptn) ;

    sequence &Reverse() ;
    std::ostream &Print(std::ostream &s) const ;
    void Print() const { Print(std::cout) ; }

    std::istream &Input(std::istream &s) ;
  } ;

  inline sequence & operator+=(sequence &e, int_type ival) 
    { e.Append(interval(ival,ival)) ; return e ; }

  inline sequence & operator+=(sequence &e, const interval &ivl)
    { e.Append(ivl) ; return e; }

  inline sequence & operator+=(sequence &e, const sequence &seq)
    { e.Append(seq) ; return e ; }

  inline sequence & operator+=(sequence &e, const intervalSet &ptn)
    { e.Append(ptn) ; return e ; }

  inline sequence operator+(const sequence & e, int_type ival) {
    sequence retp = e ;
    retp.Append(interval(ival,ival)) ;
    return retp ;
  }

  inline sequence operator+(const int &ival, const sequence &e) {
    sequence retp = interval(ival,ival) ;
    retp.Append(e) ;
    return retp ;
  }

  inline sequence operator+(const sequence &e, const interval &ivl) {
    sequence retp = e ;
    retp.Append(ivl) ;
    return retp ;
  }

  inline sequence operator+(const interval &ivl, const sequence &e) {
    sequence retp = ivl ;
    retp.Append(e) ;
    return retp ;
  }

  inline sequence operator+(const sequence &e, const sequence &ptn) {
    sequence retp = e ;
    retp.Append(ptn) ;
    return retp ;
  }

  inline sequence operator+(const sequence &e, const intervalSet &ptn) {
    sequence retp = e ;
    retp.Append(ptn) ;
    return retp ;
  }

  inline bool operator==(const sequence &es1, const sequence &es2) {
    return es1.Equal(es2) ;
  }

  inline bool operator!=(const sequence &es1, const sequence &es2) {
    return !es1.Equal(es2) ;
  }

  inline bool operator<(const sequence &es1, const sequence &es2) {
    return es1.less_than(es2) ;
  }

  inline bool operator>(const sequence &es1, const sequence &es2) {
    return es1.greater_than(es2) ;
  }

  inline bool operator<=(const sequence &es1, const sequence &es2) {
    return !(es1.greater_than(es2)) ;
  }

  inline bool operator>=(const sequence &es1, const sequence &es2) {
    return !(es1.less_than(es2)) ;
  }

  inline sequence sequence::operator>>(int rotval) const  {
    sequence retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  inline sequence sequence::operator<<(int rotval) const {
    sequence retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  inline std::ostream & operator<<(std::ostream &s, const sequence &e)
    { return e.Print(s) ; }
    
  inline std::istream & operator>>(std::istream &s, sequence &e)
    { return e.Input(s) ;}

  inline intervalSet::intervalSet(const sequence &seq) {
    for(int i=0;i<seq.num_intervals();++i)
      Union(seq[i]) ;
  }


  typedef Loci::int_type Entity ;

  template<class Op> inline void do_loop(const intervalSet &iset, Op f) {
    for(int i=0;i<iset.num_intervals();++i) {
      const Loci::int_type stop = iset[i].second ;
      for(Loci::int_type indx=iset[i].first;indx<=stop;++indx)
        f(indx) ;
    }
  }

  template<class Op> inline void do_loop(const sequence &seq, Op f) {
    for(int i=0;i<seq.num_intervals();++i) {
      const bool dir = seq[i].first>seq[i].second?false:true ;
      const Loci::int_type stop = seq[i].second + (dir?1:-1) ;
      for(Loci::int_type indx=seq[i].first;indx!=stop;dir?++indx:--indx) 
        f(indx) ;
    }
  }
  
  template<class T> inline void do_loop(const sequence &seq, T *cp,
                                        void(T::*pmf)(Entity) ) {
    for(int i=0;i<seq.num_intervals();++i) {
      const bool dir = seq[i].first>seq[i].second?false:true ;
      const Loci::int_type stop = seq[i].second + (dir?1:-1) ;
      for(Loci::int_type indx=seq[i].first;indx!=stop;dir?++indx:--indx) 
        (cp->*pmf)(indx) ;
    }
  }

}
  

#define FORALL(var,indx) \
{  const Loci::intervalSet &__p = var ; \
 for(int __ii=0;__ii<__p.num_intervals();++__ii) \
  for(Loci::int_type indx=__p[__ii].first;indx<=__p[__ii].second;++indx)
      
#define ENDFORALL }

#ifdef OLD
#define DOLOOP(seq,indx) \
{ const Loci::sequence &__s = seq ; \
  for(int __ii=0;__ii<seq.num_intervals();++__ii) { \
   const bool __dir = __s[__ii].first>__s[__ii].second?false:true ; \
   const Loci::int_type __stop = __s[__ii].second + (__dir?1:-1) ; \
   for(Loci::int_type indx=__s[__ii].first;indx!=__stop ;__dir?++indx:--indx)

#define ENDDO }}
#endif
#endif

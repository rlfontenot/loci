#ifndef STORE_H
#define STORE_H

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

#include <hdf5CC/H5cpp.h>
//#include <hdf5_traits.h>
#include <hdf5_write_template.h>
#include <Map.h>
#include <Tools/intervalSet.h>
#include <algorithm>
#include <functional>

namespace Loci {
#ifdef VERBOSE
  extern ofstream debugout[] ;
#endif
  extern int MPI_processes;
  extern int MPI_rank ;
  extern int num_threads ;

  class Map ;

  template<class T> class storeRepI : public storeRep {
    T *alloc_pointer ;
    T *base_ptr ;
    entitySet store_domain ;
  public:
    storeRepI() { alloc_pointer = 0 ; base_ptr = 0 ; }
    storeRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~storeRepI()  ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    virtual const entitySet &domain() const ;
    T * get_base_ptr() const { return base_ptr ; }
  } ;

  template<class T> void storeRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ; int size = ptn.Max()-top+1 ;
      alloc_pointer = new T[size] ;
      base_ptr = alloc_pointer - top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  template<class T> std::ostream &storeRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << std::endl ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }
    


  template<class T> std::istream &storeRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;
        
    FORALL(e,ii) {
      s >> base_ptr[ii] ;
    } ENDFORALL ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  template<class T> void storeRepI<T>::readhdf5( H5::Group group){
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    entitySet en=get_store_domain(group,traits_output_type);
    allocate(en);
    //cout<<"read "<<en<<endl;
    store_hdf5read(group,traits_output_type,base_ptr,en);
  }

  template<class T> void storeRepI<T>::writehdf5( H5::Group group,entitySet& en) const{
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    //entitySet en=domain();
    //cout<<"write "<<en<<endl;
    store_hdf5write(group,traits_output_type,base_ptr,en);
  }

  template<class T>  storeRepI<T>::~storeRepI<T>() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }
    
  template<class T>  const entitySet &storeRepI<T>::domain() const {
    return store_domain ;
  }

  template<class T>
    storeRep *storeRepI<T>::new_store(const entitySet &p) const {
    return new storeRepI<T>(p)  ;
  }

  template<class T> store_type storeRepI<T>::RepType() const {
    return STORE ;
  }

  template<class T> class store : public store_instance {
    typedef storeRepI<T> storeType ;
    T* base_ptr ;
  public:
    typedef T containerType ;
    store() { setRep(new storeType); }
    store(store &var) { setRep(var.Rep()) ; }
    store(storeRepP &rp) { setRep(rp) ; }

    virtual ~store() ;
    virtual void notification() ;

    store<T> & operator=(store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

    //    operator storeRepP() { return Rep() ; }
    T &elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
    const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
            
    T &operator[](int indx) { return elem(indx); }
    const T&operator[](int indx) const { return elem(indx); }

  } ;

  template<class T> store<T>::~store<T>() { }
    
  template<class T> void store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const store<T> &t)
    { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, store<T> &t)
    { return t.Input(s) ; }



  template<class T> class const_store : public store_instance {
    typedef storeRepI<T> storeType ;
    const T * base_ptr ;
  public:
    typedef T containerType ;
    const_store() { setRep(new storeType) ; }
    const_store(store<T> &var) { setRep(var.Rep()) ; }
    const_store(const_store &var) { setRep(var.Rep()) ; }
    const_store(storeRepP &rp) { setRep(rp) ; }

    virtual ~const_store() ;
    virtual void notification() ;


    virtual instance_type access() const  ;
        
    const_store<T> & operator=(const_store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_store<T> & operator=(store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

    //    operator storeRepP() { return Rep() ; }

    const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!Rep()->domain().inSet(indx)) ;
#endif 
      return base_ptr[indx]; }
            
    const T&operator[](int indx) const { return elem(indx); }

  } ;

  template<class T> const_store<T>::~const_store<T>() { }
    
  template<class T> void const_store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
    const_store<T>::access() const { return READ_ONLY; }
        
  template<class T> storeRepP storeRepI<T>::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    store<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
      
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T> void storeRepI<T>::copy(storeRepP &st, const entitySet &context)  {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr ==0)) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }
  template<class T> void storeRepI<T>::gather(const Map &m, storeRepP &st,
                                              const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T> void storeRepI<T>::scatter(const Map &m, storeRepP &st,
                                              const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }
  
  template <class T> int storeRepI<T>::pack_size( const entitySet &e) {
    int size ;
    size = sizeof(T) * e.size() ;
    return(size) ;
  }
  
  
  template <class T> void storeRepI<T>::pack(void * ptr, int &loc, int &size,  const entitySet &e )  {
    for(int i = 0; i < e.num_intervals(); i++) {
      const Loci::int_type begin = e[i].first ;
      int t = e[i].second - e[i].first + 1 ;  
      MPI_Pack(&base_ptr[begin], t * sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
    }
#ifdef VERBOSE
    debugout[MPI_rank] << "packing("<<e<<")"<<endl ;
    for(entitySet::const_iterator ei = e.begin(); ei != e.end(); ++ei)
      debugout[MPI_rank] << "   packing   " << base_ptr[*ei]
                         << "   into   " << *ei  << endl ;
#endif
  }
  
  template <class T> void storeRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
	const Loci::int_type stop = seq[i].second ;
	for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx) 
	  MPI_Unpack(ptr, size, &loc, &base_ptr[indx], sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      } else {
	Loci::int_type indx = seq[i].first ;
	int t = seq[i].second - seq[i].first + 1 ; 
	MPI_Unpack(ptr, size, &loc, &base_ptr[indx], t * sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ; 
      }
    }
#ifdef VERBOSE
    debugout[MPI_rank] << "unpack(" << seq << ")" << endl ;
    for(sequence::const_iterator si = seq.begin(); si != seq.end(); si++) 
      debugout[MPI_rank] << "   unpacking   " << base_ptr[*si]
                         <<"    into   " << *si << endl ;
#endif
  }

}

#ifdef GXX_FIXES

// These functions are required when using G++, not because they are actually
// used, but because G++'s instantiation mechanism generates references to
// them.

static inline ostream& operator << (ostream & s, vector<int> &) {
    cerr << "unimplemented operator<<(ostream & s, vector<int> &)" << endl;
    abort();
    return s;
}

static inline istream& operator >> (istream & s, vector<int> &) {
    cerr << "unimplemented operator>>(istream & s, vector<int> &)" << endl;
    abort();
    return s;
}

#endif // GXX_FIXES

#endif

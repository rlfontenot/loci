#ifndef STOREVEC_H
#define STOREVEC_H 

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>
#include <Tools/stream.h>
#include <store_rep.h>

#include <Tools/lmutex.h>
#include <hdf5CC/H5cpp.h>

#include <Map.h>
using std::pair ;
using std::make_pair ;


namespace Loci {
  template <class T> struct Scalar {
    T val ;
    Scalar(T v) : val(v) { }
  } ;

  template <class T> Scalar<T> mk_Scalar(T v) { return Scalar<T>(v) ;} 
  
  template <class T> class Vect ;
  
  template <class T> class const_Vect {
  public:
    friend class Vect<T> ;
    const T *ptr ;
#ifdef BOUNDS_CHECK
    int size ;
#endif
  public:
    const_Vect(const T *p ,int sz) {
      ptr = p ;
#ifdef BOUNDS_CHECK
      size = sz ;
#endif
    }
    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    operator const T *() const {
      return ptr ;
    }
  } ;


  template <class T> class Vect {
  public:
    T *ptr ;
    int size ;
  public:
    Vect(T *p ,int sz) {
      ptr = p ;
      size = sz ;
    }
    void operator=(const Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }
    void operator=(const const_Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }

    template <class S> void operator=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ = s.val ;
    }

    template <class S> void operator+=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ += s.val ;
    }
      
    template <class S> void operator*=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ *= s.val ;
    }
      
    template <class S> void operator-=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ -= s.val ;
    }
      
    template <class S> void operator/=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ /= s.val ;
    }

    template <class S> void operator+=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    template <class S> void operator-=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    operator T*() {
      return ptr ;
    }
    operator const T *() const {
      return ptr ;
    }
  } ;

  template <class T> class Mat ;
  
  template <class T> class const_Mat_partial {
    const T *ptr ;
    int size ;
    int i ;
  public:
    const_Mat_partial(const T *p,int sz, int ii) :ptr(p),size(sz),i(ii) {}
    const T &operator[](int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
  } ;
  typedef unsigned char pivot_type ;

  
  template <class T> class const_Mat {
  public:
    const T *ptr ;
    int size ;
    friend class Mat<T> ;
  public:
    const_Mat(const T *p ,int sz) : ptr(p),size(sz){}

    const_Mat_partial<T> operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return const_Mat_partial<T>(ptr,size,idx) ;
    }

    template<class S> void solve_lu(const S *b, S *x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class T1,class T2> void solve_lu(const_Vect<T1> b, T2 *x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class T1,class T2> void solve_lu(const T1 *b, Vect<T2> x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class S> void solve_lu_pivot(const S *b, S *x,const pivot_type *pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        S xi = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        S xi = x[i] ;
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    template<class T1,class T2> void solve_lu_pivot(const_Vect<T1> b, Vect<T2> x,const_Vect<pivot_type> pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        T2 xi = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }

      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        T2 xi = x[i] ;
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    template<class T1,class T2> void solve_lu_pivot(const T1 *b, Vect<T2> x,const pivot_type *pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        T2 xi = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        T2 xi = x[i] ;
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    template<class Tin,class Tout>
    void dotprod_accum(const Tin *vin, Tout *vout) const {
      const T *Aij = ptr;
      for(int j=0;j<size;++j) {
        const Tin in = vin[j] ;
        for(int i=0;i<size;++i,Aij++)
          vout[i] += (*Aij)*in ;
      }
    } 
    template<class Tin,class Tout>
    void dotprod_accum(const Vect<Tin> &vin, Tout *vout) const {
      const T *Aij = ptr;
      for(int j=0;j<size;++j) {
        const Tin in = vin[j] ;
        for(int i=0;i<size;++i,Aij++)
          vout[i] += (*Aij)*in ;
      }
    }
    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> &vin, Tout *vout) const {
      const T *Aij = ptr;
      for(int j=0;j<size;++j) {
        const Tin in = vin[j] ;
        for(int i=0;i<size;++i,Aij++)
          vout[i] += (*Aij)*in ;
      }
    }
 } ;

  template <class T> class Mat_partial {
    T *ptr ;
    int size ;
    int i ;
  public:
    Mat_partial(T *p,int sz, int ii) : ptr(p),size(sz),i(ii) {}
    T &operator[](int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
  } ;

  
  template <class T> class Mat {
  public:
    T *ptr ;
    int size ;
  public:
    Mat(T *p ,int sz) : ptr(p),size(sz) {}
    void operator=(const Mat<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ = *p2++ ;
    }
    void operator=(const const_Mat<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ = *p2++ ;
    }

    template <class S> void operator=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ = s.val ;
    }

    template <class S> void operator+=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ += s.val ;
    }
      
    template <class S> void operator*=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ *= s.val ;
    }
      
    template <class S> void operator-=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ -= s.val ;
    }
      
    template <class S> void operator/=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ /= s.val ;
    }
    
    template <class S> void operator+=(const Mat<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Mat<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Mat<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ -= *p2++ ;
    }

    template <class S> void operator-=(const const_Mat<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
      
      for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
          *p1++ -= *p2++ ;
    }
    
    Mat_partial<T> operator[](int idx) {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }
    Mat_partial<T> operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    void decompose_lu() {
      // GAXPY from of LU decomposition algorithm 
      T *Aj = ptr ;
      for(int j=0;j<size;++j,Aj += size) {
        const T *Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }

        const T Ajjr = 1./Aj[j] ;
        for(int i=j+1;i<size;++i)
          Aj[i] *= Ajjr ;
      }
    }
    
    void decompose_lu_pivot(pivot_type *pivot) {
      pivot_type piv[256] ;  // Maximum matrix size for pivoting
      for(int i=0;i<size;++i)
        pivot[i] = i ;
      T *Aj = ptr ;
      for(int j=0;j<size;++j,Aj+= size) {
        for(int k=0;k<j;++k)
          if(k!=piv[k])
            std::swap(Aj[k],Aj[piv[k]]) ;
        T *Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        int mu = j ;
        for(int k=j+1;k<size;++k)
          if(abs(Aj[mu]) < abs(Aj[k]))
            mu = k ;
        piv[j] = mu ;
        if(j!= mu)
          std::swap(pivot[j],pivot[mu]) ;
        Ak = ptr ;
        for(int k=0;k<j+1;++k,Ak += size)
          if(j != piv[j])
            std::swap(Ak[j],Ak[piv[j]]) ;
        if(Aj[j] != 0) {
          T ajjr = 1./Aj[j] ;
          for(int i=j+1;i<size;++i)
            Aj[i] *= ajjr ;
        }
      }
    }

    template<class S> void solve_lu(const S *b, S *x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;

      }
    }

    template<class S> void solve_lu(const_Vect<S> &b, S *x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class S> void solve_lu(const S *b, Vect<S> &x) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      S *y = b ;
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class S> void solve_lu_pivot(const S *b, S *x,const pivot_type *pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class S> void solve_lu_pivot(const_Vect<S> &b, S *x,const pivot_type *pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step

      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    template<class S> void solve_lu_pivot(const S *b, Vect<S> &x,const pivot_type *pivot) const {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T *Aj = ptr ;
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *Aj = Ai + size ;
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
  } ;

  
  template<class T> class storeVecRepI : public storeRep {
    entitySet store_domain ;
    T *alloc_pointer ;
    T *base_ptr ;
    int size ;
    lmutex mutex ;
  public:
    storeVecRepI() {
      alloc_pointer = 0 ; base_ptr = 0 ; size=0; }
    storeVecRepI(const entitySet &p) {
      size = 0; alloc_pointer=0 ; allocate(p) ; }
    virtual ~storeVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet &e ) ;
    virtual void pack(void * ptr, int &loc, int &size, const entitySet &e ) ;
    virtual void unpack(void * ptr, int &loc, int &size, const sequence &seq) ;
    virtual store_type RepType() const ;
    virtual const entitySet &domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    virtual void set_elem_size(int sz) ;

    T * get_base_ptr() const { return base_ptr ; }
    int get_size() const { return size ; }
  } ;

    template<class T> std::ostream &storeVecRepI<T>::Print(std::ostream &s) const
    {
      s << '{' << domain() << std::endl ;
      s << size << std::endl ;
    
      FORALL(domain(),ii) {
        T * p = base_ptr + ii*size ;
        for(int i=0;i<size;++i,++p)
          s << *p << " " ;
        s << std::endl ;
      }ENDFORALL ;
      s << '}' << std::endl ;
      return s ;
    }


  template<class T> std::istream &storeVecRepI<T>::Input(std::istream &s)
    {
      char ch ;
    
      do ch = s.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
        return s ;
      }
      entitySet e ;
      int sz ;
      s >> e ;
      s >> sz ;
      set_elem_size(sz) ;
      allocate(e) ;

      FORALL(e,ii) {
        T * p = base_ptr + ii*size ;
        for(int i=0;i<size;++i,++p)
          s >> *p ;
      } ENDFORALL ;
    
      do ch = s.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        s.putback(ch) ;
      }
      return s ;
    }

  template<class T> void storeVecRepI<T>::readhdf5( H5::Group group){
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    entitySet en=get_store_domain(group,traits_output_type);
    int sz=get_storeVec_size(group,traits_output_type);

    set_elem_size(sz);
    allocate(en);
    cout<<"Vec read "<<en<<endl;
    storeVec_hdf5read(group,traits_output_type,base_ptr,en,sz);
  }

  template<class T> void storeVecRepI<T>::writehdf5( H5::Group group,entitySet& en) const{
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    cout<<"Vec write "<<en<<endl;
    storeVec_hdf5write(group,traits_output_type,base_ptr,en,size);
  }

  template<class T> void storeVecRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(size != 0) {
      fatal(size<1) ;
      if(ptn != EMPTY) {
        int top = ptn.Min() ; int sza = (ptn.Max()-top+1)*size ;
        alloc_pointer = new T[sza] ;
        base_ptr = alloc_pointer - top*size ;
      }
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  template<class T> storeVecRepI<T>::~storeVecRepI<T>() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  template<class T>
    storeRep *storeVecRepI<T>::new_store(const entitySet &p) const {
    return new storeVecRepI<T>(p)  ;
  }

  template<class T> store_type storeVecRepI<T>::RepType() const {
    return STORE ;
  }

  template<class T> const entitySet &storeVecRepI<T>::domain() const {
    return store_domain ;
  }

  template<class T> void storeVecRepI<T>::set_elem_size(int sz) {
    //    bmutex l(mutex) ;
    mutex.lock() ;
    if(size != sz) {
      if(size != 0)
        warn(size != sz) ;
    
      size = sz ;
      fatal(sz<1) ;
      allocate(store_domain) ;
    }
    mutex.unlock() ;
  }
      
  template<class T> class storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    storeVec() {setRep(new storeType) ;}
    storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    storeVec(storeRepP &rp) { setRep(rp) ;}

    virtual ~storeVec() ;
    virtual void notification() ;

    storeVec<T> & operator=(storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    storeVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      Rep()->set_elem_size(size) ;
    }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size ; }

    const entitySet &domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }
    Vect<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr+(indx*size),size) ; }
    Vect<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr+(indx*size),size) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  template<class T> storeVec<T>::~storeVec<T>() { }

  template<class T> void storeVec<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const storeVec<T> &t)
    { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, storeVec<T> &t)
    { return t.Input(s) ; }


  template<class T> class const_storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_storeVec() { setRep(new storeType) ; }
    const_storeVec(const_storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_storeVec() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_storeVec<T> & operator=(storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_storeVec<T> & operator=(const_storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_storeVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }

    int vecSize() const { return size ; }

    const entitySet& domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }
    const_Vect<T> elem(int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Vect<T>(base_ptr+(indx*size),size) ; }
    const_Vect<T> operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Vect<T>(base_ptr+(indx*size),size) ;
    }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

  } ;

  template<class T> const_storeVec<T>::~const_storeVec<T>() { }

  template<class T> void const_storeVec<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
    const_storeVec<T>::access() const
    { return READ_ONLY ; }


  template<class T> inline std::ostream & operator<<(std::ostream &s, const const_storeVec<T> &t)
    { return t.Print(s) ; }


  template <class T> storeRepP storeVecRepI<T>::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    storeVec<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template <class T> void storeVecRepI<T>::copy(storeRepP &st,
                                                const entitySet &context) {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + i*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }

  template <class T> void storeVecRepI<T>::gather(const Map &m, storeRepP &st,
                                                  const entitySet &context) {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal(base_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + i*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[m[i]][j] ;
    } ENDFORALL ;
  }

  template <class T> void storeVecRepI<T>::scatter(const Map &m, storeRepP &st,
                                                   const entitySet &context) {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    fatal(base_ptr == 0) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = base_ptr + m[i]*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }
  template <class T> int storeVecRepI<T>::pack_size( const entitySet &e) {
    int size, M ;
    M = get_size() ;
    size = sizeof(T) * e.size() * M  ;
    return(size) ;
  }
  
  template <class T> void storeVecRepI<T>::pack(void * ptr, int &loc, int &size, const entitySet &e ) {
    int M = get_size() ;
    for(int i = 0; i < e.num_intervals(); ++i) {
      Loci::int_type indx1 = e[i].first ;
      Loci::int_type stop = e[i].second ;
      T *p = base_ptr + M * indx1 ;
      int t = (stop - indx1 + 1) * M ;
      MPI_Pack(p, t * sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
    }
  }
  
  template <class T> void storeVecRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    int M = get_size() ;
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
	const Loci::int_type indx1 = seq[i].first ;
	const Loci::int_type stop = seq[i].second ;
	for(Loci::int_type indx = indx1; indx != stop-1; --indx) {
	  T *p = base_ptr + M * indx ;
	  MPI_Unpack(ptr, size, &loc, p, sizeof(T) * M, MPI_BYTE, MPI_COMM_WORLD) ;
	}
      }
      else {
	Loci::int_type indx1 = seq[i].first ;
	Loci::int_type stop = seq[i].second ;
	T *p = base_ptr + M * indx1 ;
	int t = (stop - indx1 + 1) * M ;
	MPI_Unpack(ptr, size, &loc, p, t * sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ; 
	}
      }
  }
  
  template<class T> class storeMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size_tot ;
    int size_dim ;
  public:
    typedef Mat<T> containerType ;
    storeMat() {setRep(new storeType) ;}
    storeMat(storeMat &var) {setRep(var.Rep()) ; }
    storeMat(storeRepP &rp) {setRep(rp) ; }

    virtual ~storeMat() ;
    virtual void notification() ;

    storeMat<T> & operator=(storeMat<T> &str)
      { setRep(str.Rep()) ; return *this ;}

    storeMat<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      size_dim = size ;
      size_tot = size*size ;
      Rep()->set_elem_size(size_tot) ; 
    }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size_dim; }
    const entitySet &domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }
    Mat<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    Mat<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  template<class T> storeMat<T>::~storeMat<T>() { }

  template<class T> void storeMat<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const storeMat<T> &t)
    { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, storeMat<T> &t)
    { return t.Input(s) ; }

  template<class T> class const_storeMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* base_ptr ;
    int size_tot ;
    int size_dim ;
  public:
    typedef const_Mat<T> containerType ;
    const_storeMat() { setRep(new storeType) ; }

    const_storeMat(const_storeMat<T> &var) { setRep(var.Rep()) ; }
    const_storeMat(storeMat<T> &var) { setRep(var.Rep()) ; }
    const_storeMat(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_storeMat() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_storeMat<T> & operator=(const_storeMat<T> &str)
      { setRep(str.Rep()) ; return *this ;}
    const_storeMat<T> & operator=(storeMat<T> &str)
      { setRep(str.Rep()) ; return *this ;}

    const_storeMat<T> & operator=(storeRepP p)
      { setRep(p) ; return *this ; }

    int vecSize() const { return size_dim; }
    const entitySet &domain() { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }

    const_Mat<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    const_Mat<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  } ;

  template<class T> const_storeMat<T>::~const_storeMat<T>() { }

  template<class T> void const_storeMat<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
    const_storeMat<T>::access() const
    { return READ_ONLY; }
        

  template<class T> inline std::ostream &
    operator<<(std::ostream &s, const const_storeMat<T> &t)
    { return t.Print(s) ; }

  template<class T> class multiStoreRepI : public storeRep {
    entitySet store_domain ;
    T **index ;
    T *alloc_pointer ;
    T **base_ptr ;
    int size ;
    lmutex mutex ;
  public:
    multiStoreRepI()
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; }
    multiStoreRepI(const entitySet &p)
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; store_domain=p;}
    multiStoreRepI(const store<int> &sizes) {
      index = 0 ; alloc_pointer=0 ; base_ptr = 0; allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    void multialloc(const store<int> &count, T ***index, T **alloc_pointer, T ***base_ptr) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~multiStoreRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual const entitySet &domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;

    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
  } ;

  
  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore(storeRepP &rp) { setRep(rp) ;}

    virtual ~multiStore() ;
    virtual void notification() ;

    multiStore<T> & operator=(multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void setSizes(const const_multiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    const entitySet &domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }

    Vect<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    Vect<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  template<class T> multiStore<T>::~multiStore() {}
  
  template<class T> void multiStore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }
  
  template<class T> class const_multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(storeRepP &rp) { setRep(rp) ;}

    virtual ~const_multiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
    const_multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(const const_multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    const entitySet &domain() const { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }

    containerType elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    containerType operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }

    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  template<class T> store_instance::instance_type
    const_multiStore<T>::access() const
    { return READ_ONLY ; }
  
  template<class T> const_multiStore<T>::~const_multiStore() {}
  
  template<class T> void const_multiStore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> void multiStoreRepI<T>::allocate(const store<int> &sizes) {
    int sz = 0 ;
    entitySet ptn = sizes.domain() ;
    store_domain = ptn ;
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    if(index) delete[] index ;
    index = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int len = ptn.Max() - top + 2 ;
      index = new T *[len] ;
      base_ptr = index - top ;
      FORALL(ptn,i) {
        sz += sizes[i] ;
      } ENDFORALL ;
      alloc_pointer = new T[sz+1] ;
      sz = 0 ;
      for(int ivl=0;ivl<ptn.num_intervals();++ivl) {
        int i = ptn[ivl].first ;
        base_ptr[i] = alloc_pointer + sz ;
        while(i<=ptn[ivl].second) {
          sz += sizes[i] ;
          ++i ;
          base_ptr[i] = alloc_pointer + sz ;
        }
      }
    }
    dispatch_notify();
  }

  
  template<class T> void multiStoreRepI<T>::multialloc(const store<int> &count, T ***index, 
						       T **alloc_pointer, T ***base_ptr ) {
    entitySet ptn = count.domain() ;
    int top = ptn.Min() ;
    int len = ptn.Max() - top + 2 ;
    T **new_index = new T *[len] ;
    T **new_base_ptr = new_index - top ;
    int sz = 0 ;
    FORALL(ptn, i) {
      sz += count[i] ;
    } ENDFORALL ;
    T *new_alloc_pointer = new T[sz + 1] ;
    sz = 0 ;
    for(int ivl = 0; ivl < ptn.num_intervals(); ++ivl) {
      int i = ptn[ivl].first ;
      new_base_ptr[i] = new_alloc_pointer + sz ;
      while(i <= ptn[ivl].second) {
	sz += count[i] ;
	++i ;
	new_base_ptr[i] = new_alloc_pointer + sz ;
      }
    }
    
    *index = new_index ;
    *alloc_pointer = new_alloc_pointer ;
    *base_ptr = new_base_ptr ;
    
  }
   
  template<class T> void multiStoreRepI<T>::setSizes(const const_multiMap &mm){
    //    bmutex l(mutex) ;
    mutex.lock() ;
    if(alloc_pointer != 0) {
      entitySet map_set = mm.domain() & store_domain ;
      entitySet problem ;
      FORALL(map_set,i) {
        if((end(i)-begin(i))>(mm.end(i)-mm.begin(i)))
          problem += i ;
      } ENDFORALL ;

      if(problem != EMPTY) {
        std::cerr << "reallocation of multiStore required for entities"
                  << problem << endl
                  << "Currently this reallocation isn't implemented."
                  << endl ;
      }
    } else {
      store<int> sizes ;
      sizes.allocate(store_domain) ;
      FORALL(store_domain,i) {
        sizes[i] = 0 ;
      } ENDFORALL ;
      entitySet map_set = mm.domain() & store_domain ;
      FORALL(map_set,i) {
        sizes[i] = (mm.end(i) - mm.begin(i)) ;
      } ENDFORALL ;
      allocate(sizes) ;
    }
    mutex.unlock() ;
  }
  
  template<class T> void multiStoreRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
    alloc_pointer = 0 ;
    index = 0 ;
    base_ptr = 0 ;
    store_domain = ptn ;
    store<int> count ;
    count.allocate(ptn) ;
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    allocate(count) ;
    dispatch_notify() ;
  }

  template<class T> multiStoreRepI<T>::~multiStoreRepI() {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
  }

  template<class T> storeRep *multiStoreRepI<T>::new_store(const entitySet &p)
    const {
    return new multiStoreRepI<T>(p) ;
  }

  template<class T> storeRepP multiStoreRepI<T>::remap(const Map &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    multiStore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }
  
  template<class T> void multiStoreRepI<T>::copy(storeRepP &st,
                                                 const entitySet &context) {
    const_multiStore<T> s(st) ;
    fatal(alloc_pointer == 0) ;
    fatal((context - domain()) != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    store<int> count ;
    count.allocate(domain()) ;
    FORALL(domain() - context, i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context, i) {
      count[i] = s.end(i) - s.begin(i) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
  
  template<class T> void multiStoreRepI<T>::gather(const Map &m, storeRepP &st,
                                                  const entitySet &context) {
    store<int> count ;
    const_multiStore<T> s(st) ;
    count.allocate(domain()) ;
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[i] = s.end(m[i])-s.begin(m[i]) ;
    } ENDFORALL ;
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;

    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j = 0; j < count[i]; ++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j = 0; j < count[i]; ++j)
        new_base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
    
  }
  template<class T> void multiStoreRepI<T>::scatter(const Map &m, storeRepP &st,
                                                  const entitySet &context) {
    
    store<int> count ;
    const_multiStore<T> s(st) ;
    count.allocate(domain()) ;
    FORALL(domain()-m.image(context),i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[m[i]] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain() - m.image(context),i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[m[i]];++j) {
        new_base_ptr[m[i]][j] = s[i][j] ;
      }
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
 
  template <class T> int multiStoreRepI<T>::pack_size(const entitySet &e ) {
    int size = 0 ;
    store<int> count ;
    count.allocate(e) ;
    FORALL(e,i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
      size += count[i] ;
    } ENDFORALL ;
    entitySet::const_iterator ei = e.begin() ;
    return(size * sizeof(T)) ;
  }
  
  template <class T> void multiStoreRepI<T>::pack(void * ptr, int &loc, int &size, const entitySet &e ) {
    
    store<int> count ;
    count.allocate(e) ;
    FORALL(e,i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    
    for(int i = 0; i < e.num_intervals(); ++i) {
      Loci::int_type indx1 = e[i].first ;
      Loci::int_type stop = e[i].second ;
      int size1 = 0 ;
      for(Loci::int_type indx = indx1; indx != stop+1; ++indx)
	size1 += count[indx] ;
      MPI_Pack(&base_ptr[indx1][0], size1 * sizeof(T), MPI_BYTE, ptr, size, &loc, MPI_COMM_WORLD) ;
    }
  }
  
  
  template <class T> void multiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    
    entitySet e ;
    store<int> count, count1 ;
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) 
	e += interval(seq[i].second, seq[i].first) ;
      else 
	e += interval(seq[i].first, seq[i].second) ;
    }
    count1.allocate(e) ;
    count.allocate(e) ;
    for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
      count1[*si] = base_ptr[*si+1] - base_ptr[*si] ;
    
    entitySet::const_iterator ei = e.begin() ;
    
    for( sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
      count[*ei] = count1[*si] ;
      ++ei ;
    }
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) 
      for(int j = 0 ; j < count[*si]; ++j) 
	new_base_ptr[*si][j] = 0 ;
    
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
    
    for(int i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
	Loci::int_type stop = seq[i].second ;
	for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
	  MPI_Unpack(ptr, size, &loc, &base_ptr[indx][0], count[indx]*sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      }
      else {
	Loci::int_type indx1 = seq[i].first ;
	Loci::int_type stop = seq[i].second ;
	int sz = 0 ;
	for(Loci::int_type indx = seq[i].first; indx != stop+1; ++indx)
	  sz += count[indx] ;
	MPI_Unpack(ptr, size, &loc, &base_ptr[indx1][0], sz * sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ; 
      }
    }
    dispatch_notify() ;
  }
  
 
  template<class T> store_type multiStoreRepI<T>::RepType() const {
    return STORE ;
  }
  
  template<class T> const entitySet &multiStoreRepI<T>::domain() const {
    return store_domain ;
  }
  
  template<class T> std::ostream &multiStoreRepI<T>::Print(std::ostream &s)
    const {
    s << '{' << domain() << endl ;
    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;
    FORALL(domain(),ii) {
      for(const T *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << ' ' ;
      s << std::endl ;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }

  template<class T> std::istream &multiStoreRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    allocate(sizes) ;
        
    FORALL(e,ii) {
      for(T *ip = begin(ii);ip!=end(ii);++ip)
        s >> *ip  ;
    } ENDFORALL ;
            
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  template<class T> void multiStoreRepI<T>::readhdf5( H5::Group group) {
    std::cerr << "readhdf5 not implemented" << std::endl ;
  }

  template<class T> void multiStoreRepI<T>::writehdf5(H5::Group group,
                                                      entitySet &en) const {
    std::cerr << "writehdf5 not implemented" << std::endl ;
  }
  
}

#endif

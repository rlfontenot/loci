#ifndef STOREVEC_H
#define STOREVEC_H 1

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#ifdef PTHREADS
#include <pthread.h>
#endif

#ifdef PTHREADS
extern pthread_mutex_t access_mutex ;
#endif

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
      warn(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
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

    template<class T1,class T2> void solve_lu_pivot(const_Vect<T1> b, Vect<T2> x,const_Vect<pivot_type> pivot) const {
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

    template<class T1,class T2> void solve_lu_pivot(const T1 *b, Vect<T2> x,const pivot_type *pivot) const {
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
        for(int k=0;k<j;++k,Ak += size)
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Aj[k] ;
        Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size)
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Aj[k] ;
        
        for(int i=j+1;i<size;++i)
          Aj[i] /= Aj[j] ;
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
        for(int k=0;k<j;++k,Ak += size)
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Aj[k] ;
        Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size)
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Aj[k] ;
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
  public:
    storeVecRepI() {
      alloc_pointer = 0 ; base_ptr = 0 ; size=0; }
    storeVecRepI(const entitySet &p) {
      alloc_pointer=0 ; allocate(p) ; }
    virtual ~storeVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual store_type RepType() const ;
    virtual const entitySet &domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
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


  template<class T> void storeVecRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(size != 0) {
      fatal(size<1) ;
      if(ptn != EMPTY) {
        int top = ptn.Min() ; int sza = (ptn.Max()-top+1)*size ;
        alloc_pointer = new(T[sza]) ;
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
    if(size != sz) {
      if(size != 0)
        warn(size != sz) ;
    
      size = sz ;
      fatal(sz<1) ;
      allocate(store_domain) ;
    }
  }
      
  template<class T> class storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    storeVec() {
      setRep(new storeType) ;
    }
    storeVec(storeVec<T> &var) {
      setRep(var.Rep()) ;
    }
    storeVec(const entitySet &ptn) {
      setRep(new storeType(ptn)) ;
    }

    virtual ~storeVec() ;
    virtual void notification() ;

    storeVec<T> & operator=(storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    storeVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
#ifdef PTHREADS
    int err = pthread_mutex_lock(&access_mutex) ;
    fatal(err==EDEADLK || err==EINVAL) ;
#endif
      Rep()->set_elem_size(size) ;
#ifdef PTHREADS
    err = pthread_mutex_unlock(&access_mutex) ;
    fatal(err==EPERM) ;
#endif
    }
    void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size ; }

    const entitySet &domain() const { return Rep()->domain() ; }
    operator storeRepP() { return Rep() ; }
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
        
    const_storeVec(const_storeVec<T> &var) { 
      setRep(var.Rep()) ;
    }
        
    const_storeVec(storeVec<T> &var) { 
      setRep(var.Rep()) ;
    }
        
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
    operator storeRepP() { return Rep() ; }
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


  template<class T> class storeMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size_tot ;
    int size_dim ;
  public:
    typedef Mat<T> containerType ;
    storeMat() {setRep(new storeType) ;}
    storeMat(storeMat &var) {setRep(var.Rep()) ; }
    storeMat(const entitySet &ptn) {setRep(new storeType(ptn)) ; }

    virtual ~storeMat() ;
    virtual void notification() ;

    storeMat<T> & operator=(storeMat<T> &str)
      { setRep(str.Rep()) ; return *this ;}

    storeMat<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
#ifdef PTHREADS
    int err = pthread_mutex_lock(&access_mutex) ;
    fatal(err==EDEADLK || err==EINVAL) ;
#endif
      size_dim = size ;
    size_tot = size*size ;
    Rep()->set_elem_size(size_tot) ; 
#ifdef PTHREADS
    err = pthread_mutex_unlock(&access_mutex) ;
    fatal(err==EPERM) ;
#endif
    }
    void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size_dim; }
    const entitySet &domain() const { return Rep()->domain() ; }
    operator storeRepP() { return Rep() ; }
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
    operator storeRepP() { return Rep() ; }

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



}

#endif

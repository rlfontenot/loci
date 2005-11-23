#ifndef MATRIX_H
#define MATRIX_H 

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

using std::pair ;
using std::make_pair ;

namespace Loci {

  template <class T> class Mat ;

  //**************************************************************************
  
  template <class T> class const_Mat_partial {
    const T * restrict ptr ;
    int size ;
    int i ;
  public:
    const_Mat_partial(const T *p,int sz, int ii) :ptr(p),size(sz),i(ii) {}
    const T & restrict operator[](int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
  } ;

  //**************************************************************************

  typedef unsigned char pivot_type ;
  
  template <class T> class const_Mat {
  public:
    const T *restrict ptr ;
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

    //************************************************************************

    template<class S> 
    void solve_lu(const S * restrict b, S *restrict x) const restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    //*************************************************************************

    template<class T1,class T2> 
    void solve_lu(const_Vect<T1> b, T2 *restrict x) const  restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T * restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu(const T1 * restrict b, Vect<T2> xin) const restrict {
      T2 * restrict x = xin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T * restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T * restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    //************************************************************************

    template<class S> 
    void solve_lu_pivot(const S * restrict b, S *restrict x,const pivot_type * restrict pivot) const restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        S xi = b[pivot[i]] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }
      // Do back solver Ux = y
      const T * restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
        S xi = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot( const_Vect<T1> b, Vect<T2> xin, 
                         const_Vect<pivot_type> pivot) const restrict {
      T2 * restrict x = xin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        T2 xi = b[pivot[i]] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }

      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T * restrict Aj = Ai + size ;
        T2 xi = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot(const T1 * restrict b, Vect<T2> xin,const pivot_type *restrict pivot) const restrict {
      T2 * restrict x = xin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        T2 xi = b[pivot[i]] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }
      // Do back solver Ux = y
      const T * restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T * restrict Aj = Ai + size ;
        T2 xi = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Tin *restrict vin, Tout * restrict vout) const restrict {
      const T * restrict Aij = ptr;
      for(int j=0;j<size;++j) {
        const Tin in = vin[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=0;i<size;++i,Aij++)
          vout[i] += (*Aij)*in ;
      }
    } 

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Vect<Tin> &vin, Tout *vout) const restrict {
      const T * restrict Aij = ptr;
      const Tin * restrict vi = vin ;
      Tout * restrict vo = vout ;
      for(int j=0;j<size;++j) {
        const Tin in = vi[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=0;i<size;++i,Aij++)
          vo[i] += (*Aij)*in ;
      }
    }

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> &vin, Tout *vout) const restrict {
      const T * restrict Aij = ptr;
      const Tin * restrict vi = vin ;
      Tout * restrict vo = vout ;
      for(int j=0;j<size;++j) {
        const Tin in = vi[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=0;i<size;++i,Aij++)
          vo[i] += (*Aij)*in ;
      }
    }

  } ;

  //************************************************************************

  template <class T> class Mat_partial {
    T *restrict ptr ;
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

  //**************************************************************************
  
  template <class T> class Mat {
  public:
    T *ptr ;
    int size ;
  public:
    Mat(T *p ,int sz) : ptr(p),size(sz) {}

    //------------------------------------------------------------------------
    void operator=(const Mat<T> &t) {
      T *restrict p1 = ptr ;
      const T *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ = *p2++ ;
    }
    //------------------------------------------------------------------------
    void operator=(const const_Mat<T> &t) {
      T *restrict p1 = ptr ;
      const T *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ = *p2++ ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator=(const Scalar<S> &restrict s) {
      T *restrict p1 = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ = s.val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const Scalar<S> &restrict s) {
      T *restrict p1 = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ += s.val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator*=(const Scalar<S> &restrict s) {
      T *restrict p1 = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ *= s.val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const Scalar<S> &restrict s) {
      T *restrict p1 = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ -= s.val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator/=(const Scalar<S> &restrict s) {
      T *restrict p1 = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ /= s.val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const Mat<S> &t) {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ += *p2++ ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const const_Mat<S> &t) {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ += *p2++ ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const Mat<S> &t) {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ -= *p2++ ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const const_Mat<S> &t) {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<size*size;++i)
        *p1++ -= *p2++ ;
    }
    //------------------------------------------------------------------------
    
    Mat_partial<T> operator[](int idx) restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    //------------------------------------------------------------------------
    Mat_partial<T> operator[](int idx) const restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    //------------------------------------------------------------------------

    void decompose_lu() restrict {
      // GAXPY from of LU decomposition algorithm 
      T *restrict Aj = ptr ;
      for(int j=0;j<size;++j,Aj += size) {
        const T *restrict Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }

        const T Ajjr = 1./Aj[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=j+1;i<size;++i)
          Aj[i] *= Ajjr ;
      }
    }
    //------------------------------------------------------------------------
    
    void decompose_lu_pivot(pivot_type *restrict pivot) restrict {
      pivot_type piv[256] ;  // Maximum matrix size for pivoting
      for(int i=0;i<size;++i)
        pivot[i] = i ;
      T *restrict Aj = ptr ;
      for(int j=0;j<size;++j,Aj+= size) {
        for(int k=0;k<j;++k)
          if(k!=piv[k])
            std::swap(Aj[k],Aj[piv[k]]) ;
        T *restrict Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
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
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j+1;++k,Ak += size)
          if(j != piv[j])
            std::swap(Ak[j],Ak[piv[j]]) ;
        if(Aj[j] != 0) {
          T ajjr = 1./Aj[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=j+1;i<size;++i)
            Aj[i] *= ajjr ;
        }
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu(const S *restrict b, S *restrict x) const restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T *Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T * restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu(const_Vect<S> &bin, S *restrict x) const restrict {
      S *restrict b = bin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T * restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu(const S *restrict b, Vect<S> &xin) const restrict {
      S *restrict x = xin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[i] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      S *y = b ;
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu_pivot(const S * restrict b, S *restrict x,const pivot_type *restrict pivot) const restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T *restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu_pivot(const_Vect<S> &bin, S *restrict x,const pivot_type *restrict pivot) const restrict {
      // Perform forward solve Ly = b, note b becomes y after this step
      S *restrict b = bin ;
      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T *restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }
    //------------------------------------------------------------------------

    template<class S> 
    void solve_lu_pivot(const S *restrict b, Vect<S> &xin,const pivot_type *restrict pivot) const restrict {
      S *restrict x = xin ;
      // Perform forward solve Ly = b, note b becomes y after this step
      for(int i=0;i<size;++i) {
        x[i] = b[pivot[i]] ;
        const T * restrict Aj = ptr ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
      }
      // Do back solver Ux = y
      const T *restrict Ai = ptr + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T *restrict Aj = Ai + size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          x[i] -= Aj[i]*x[j] ;
        x[i] = x[i]/Ai[i] ;
      }
    }

    //------------------------------------------------------------------------

  } ;

}

#endif

  

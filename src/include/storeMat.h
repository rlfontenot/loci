#ifndef STOREMAT_H
#define STOREMAT_H 

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>
#include <matrix.h>

#include <hdf5CC/H5cpp.h>

namespace Loci {
 
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
    const entitySet domain() const { return Rep()->domain() ; }

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

  //**************************************************************************/

  template<class T> 
  storeMat<T>::~storeMat<T>() { }

  //**************************************************************************/

  template<class T> 
  void storeMat<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const storeMat<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, storeMat<T> &t)
  { return t.Input(s) ; }

  //*************************************************************************/

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
    const entitySet domain() { return Rep()->domain() ; }
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

  //***************************************************************************/

  template<class T> 
  const_storeMat<T>::~const_storeMat<T>() { }

  //***************************************************************************/

  template<class T> 
  void const_storeMat<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  //***************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_storeMat<T>::access() const
  { return READ_ONLY; }
        
  //***************************************************************************/

  template<class T> inline std::ostream &
  operator<<(std::ostream &s, const const_storeMat<T> &t)
  { return t.Print(s) ; }

  //***************************************************************************/

}

#endif

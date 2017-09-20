//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef MULTISTORE_DEF_H
#define MULTISTORE_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <Tools/debug.h>
#include <Tools/tools.h>
#include <data_traits.h>
#include <store_rep.h>
#include <Tools/lmutex.h>
#include <DMap.h>
#include <multiMap.h>
#include <dist_internal.h>

namespace Loci {

  template<class T> class multiStoreRepI : public storeRep {
    entitySet store_domain ;
    T **index ;
    T *alloc_pointer ;
    T **base_ptr ;
    int size ;
    lmutex mutex ;
    
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info&fi, entitySet &en);
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset);
    
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq );
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq);
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    multiStoreRepI()
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0;}

    multiStoreRepI(const entitySet &p)
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; store_domain=p;}

    multiStoreRepI(const store<int> &sizes) {
      index = 0 ; alloc_pointer=0 ; base_ptr = 0; allocate(sizes) ; }

    void allocate(const store<int> &sizes) ;
    virtual void shift(int_type offset) ;
    void multialloc(const store<int> &count, T ***index, T **alloc_pointer, T ***base_ptr) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~multiStoreRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;

    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;


  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
    multiStore(const multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(storeRepP rp) { setRep(rp) ;}
    virtual ~multiStore() {}
    virtual void notification() ;
    multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void allocate(const_store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      store<int> csizes(sizes.Rep()) ;
      fatal(p==0) ;
      p->allocate(csizes) ;
    }
    void setSizes(const const_multiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    const entitySet domain() const { return Rep()->domain() ; }
    Vect<T> elem(int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }
    Vect<T> operator[](int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int index) const { return end(index)-begin(index); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;

  template<class T> class const_multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
    const_multiStore(const const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(const multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_multiStore<T> & operator=(const const_multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(storeRepP rp) { setRep(rp) ;}
    virtual ~const_multiStore() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    const entitySet domain() const { return Rep()->domain() ; }
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
    int vec_size(int index) const { return end(index)-begin(index); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;
  
} // end of namespace Loci

#endif

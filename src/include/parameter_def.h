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
#ifndef PARAMETER_DEF_H
#define PARAMETER_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>


#include <mpi.h>
#include <Tools/debug.h>
#include <store_rep.h>
//#include <distribute.h>
#include <data_traits.h>
//#include <gstore_rep.h>
#include <gkeyspace.h>

namespace Loci {
  
  template<class T> class paramRepI : public storeRep {
    entitySet store_domain ;
    T *alloc_ptr ;
    T *base_ptr ;
    
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, frame_info &fi, const entitySet &en);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const entitySet &en);

    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g,const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size );
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size );

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size);
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size);
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    paramRepI() {
      alloc_ptr=0; base_ptr = 0 ;
      store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ;
      allocate(store_domain) ;
    }
    
    paramRepI(const entitySet &p) {
      alloc_ptr=0 ; base_ptr=0; store_domain=p ;
      allocate(store_domain) ;
    }
    virtual void allocate(const entitySet &p)  ;
    virtual void erase(const entitySet& rm) ;
    virtual void guarantee_domain(const entitySet& include) ;
    //    virtual gStoreRepP copy2gstore()const;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute_omd(const std::vector<entitySet>& dom_ptn,
                     const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual void shift(int_type offset) ;
    virtual ~paramRepI() ;
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void fast_copy(storeRepP& st, const entitySet& context);
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq)  ;

    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    T *get_param() { return base_ptr ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;

 
  //**************************************************************************/

  template<class T> class param : public store_instance {
    typedef paramRepI<T> paramType ;
    T * data ;
  private:
    void connect_data(){
      NPTR<paramType> p(Rep()) ;
      if(p != 0) data = p->get_param();
    }
      
  public:
    typedef T containerType ;
    param() { setRep(new paramType) ; connect_data();}
    param(const param &var) {
      setRep(new paramType) ;
      connect_data();
      *data = *(var.data) ;
      set_entitySet(var.domain());
    }
    param(storeRepP rp) { setRep(rp); connect_data();}

    virtual ~param() ;

    param & operator=(const param &p) {
      *data = *(p.data) ;
      set_entitySet(p.domain());
      return *this ;
    }
    param & operator=(storeRepP p) {setRep(p) ; connect_data(); return *this ; }
    param & operator=(const T &v) { *data = v ; return *this ; }

    virtual void notification() ;

    T * operator->() { return data ; }
    const T * operator->() const { return data ; }

    T * operator&() { return data ; }
    const T * operator &() const { return data ; }

    T &operator*() { return *data ; }
    const T &operator*() const { return *data ; }

    T &operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    void set_entitySet(const entitySet &ptn) {Rep()->allocate(ptn); }
    //    virtual gStoreRepP copy2gstore()const{return Rep()->copy2gstore();}
    entitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

 
  //**************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const param<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T>
  inline std::istream & operator>>(std::istream &s, param<T> &t)
  { return t.Input(s) ; }

  //**************************************************************************/

  template<class T>
  class const_param : public store_instance {
    typedef T containerType ;
    typedef paramRepI<T> paramType ;
    const T * restrict data ;
    const_param(const const_param<T> &var) {setRep(var.Rep()) ;}
    const_param(param<T> &var) { setRep(var.Rep()) ; }
  public:
    const_param() { setRep(new paramType) ; }
    const_param(storeRepP rp) { setRep(rp); }

    virtual ~const_param() ;

    const_param & operator=(const_param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(storeRepP p)
    { setRep(p) ; return *this ;}

    virtual void notification() ;
    virtual instance_type access() const ;

    const T * restrict operator->() const { return data ; }

    const T * restrict operator &() const { return data ; }

    const T & restrict operator*() const { return *data ; }

    const T & restrict operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    entitySet domain() const { return Rep()->domain(); }
    //    virtual gStoreRepP copy2gstore()const{return Rep()->copy2gstore();}
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

 
}

#endif

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
#ifndef MULTIMAP_H
#define MULTIMAP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <Map_rep.h>

#include <Map.h>
#include <store.h>

namespace Loci {
  class multiMapRepI : public MapRep {
    entitySet store_domain ;
    int **index ;
    int *alloc_pointer ;
    int **base_ptr ;
  public:
    multiMapRepI() { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; }
    multiMapRepI(const store<int> &sizes) {
      index = 0 ;
      alloc_pointer = 0 ;
      base_ptr = 0 ;
      allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    virtual void allocate(const entitySet &ptn) ;
    virtual ~multiMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context)  ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;

    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    virtual void pack(void *ptr, int &loc,
                      int &size, const entitySet &e, const Map& remap) ;
    virtual void unpack(void *ptr, int &loc,
                        int &size, const sequence &seq, const dMap& remap) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
    int ** get_base_ptr() const { return base_ptr ; }
    int *begin(int indx) { return base_ptr[indx] ; }
    int *end(int indx) { return base_ptr[indx+1] ; }
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  private:
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ; 
  } ;
      
  class multiMap : public store_instance {
    friend class const_multiMap ;
    typedef multiMapRepI MapType ;
    int **base_ptr ;
  public:

    class arrayHelper {
      int *first ;
      int *last ;
    public:
      arrayHelper(int *f, int *l) : first(f), last(l){}
      int size() { return last-first ; }
      int &operator[](int indx) { return first[indx] ; }
      int *begin() { return first ; }
      int *end() { return last; }
    } ;
    class arrayHelper_const {
      const int *first ;
      const int *last ;
    public:
      arrayHelper_const(const int *f, const int *l) : first(f), last(l){}
      int size() { return last-first ; }
      const int &operator[](int indx) { return first[indx] ; }
      const int *begin() { return first ; }
      const int *end() { return last; }
    } ;
        
    // These should be private, as they only perform a shallow copy,
    // which is dangerous.  For now we leave them public because it would
    // be too difficult to fix properly.
    multiMap(const multiMap &var) { setRep(var.Rep()) ; }
    multiMap & operator=(const multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    
    multiMap() { setRep(new MapType) ; }
        
    multiMap(const store<int> &sizes) { setRep( new MapType(sizes) ); }

    multiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~multiMap() ;
    virtual void notification() ;

    multiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    arrayHelper elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper_const const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]) ;
    }

    arrayHelper operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper_const operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]) ;
    }
      
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    int *begin(int indx) { return base_ptr[indx] ; }
    int *end(int indx) { return base_ptr[indx+1] ; }
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;
  
  inline std::ostream & operator<<(std::ostream &s, const multiMap &m)
    { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, multiMap &m)
    { return m.Input(s) ; }

  class const_multiMap : public store_instance {
    typedef multiMapRepI MapType ;
    const int * const * base_ptr ;
    const_multiMap(const_multiMap &var) {  setRep(var.Rep()) ; }
    const_multiMap & operator=(const const_multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_multiMap & operator=(const multiMap &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    class arrayHelper_const {
      const int *first ;
      const int *last ;
    public:
      arrayHelper_const(const int *f, const int *l) : first(f), last(l){}
      int size() { return last-first ; }
      const int &operator[](int indx) { return first[indx] ; }
      const int *begin() { return first ; }
      const int *end() { return last; }
    } ;

    const_multiMap() { setRep(new MapType) ; }
    
    
    const_multiMap(multiMap &var) { setRep(var.Rep()) ; }
    
    const_multiMap(storeRepP rp) { setRep(rp) ; }
    
    virtual ~const_multiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_multiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }
    //    operator storeRepP() { return Rep() ; }
    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    arrayHelper_const const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]); }
    arrayHelper_const operator[](int indx) const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]); }
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;


/*
  inline std::ostream & operator<<(std::ostream &s, const const_multiMap &m)
    { return m.Print(s) ; }
*/

  void inverseMap(multiMap &result,
                  const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(multiMap &result,
                  const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(multiMap &result,
                  const const_Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(multiMap &result,
                  const const_multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;



}

#endif

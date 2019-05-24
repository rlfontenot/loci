//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#ifndef STORE_REP_H
#define STORE_REP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/nptr.h>
#include <entitySet.h>
#include <istream>
#include <ostream>

#include <data_traits.h>

#include <mpi.h>

namespace Loci {
  enum store_type { STORE, PARAMETER, MAP, CONSTRAINT, BLACKBOX } ;

  class Map ;
  class dMap ;
  class storeRep ;

  typedef NPTR<storeRep> storeRepP ;
  typedef const_NPTR<storeRep> const_storeRepP ;
  
  struct frame_info {
    int is_stat ; 
    int size ;
    std::vector<int> first_level ;
    std::vector<int> second_level ;
    frame_info() {
      is_stat = 0 ;
      size = 0 ;
    }
    frame_info(int a , int b) {
      is_stat = a ;
      size = b ;
      
    }
    frame_info(int a , int b, std::vector<int> c, std::vector<int> d) {
      is_stat = a ;
      size = b ;
      
      first_level = c ;
      second_level = d ;
    } 
    frame_info(const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size)
	first_level = fi.first_level ;
      if(is_stat) 
	second_level = fi.second_level ;
    }
    
    frame_info &operator = (const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size) 
	first_level = fi.first_level ;
      if(is_stat)
	second_level = fi.second_level ;
      
      return *this ;
    }
  } ;
  class storeRep : public NPTR_type {
    int domainKeySpace ;
  public:
    storeRep() { domainKeySpace = 0 ; }
    virtual ~storeRep() ;
    virtual void allocate(const entitySet &p) = 0 ;
    // erase part of the domain, useful for dynamic containers,
    // default behavior is doing nothing.
    virtual void erase(const entitySet& rm) {
      std::cerr << "storeRep.erase() is not implemented yet" << std::endl ;
      abort() ;
    }
    // this method is used to invalidate part of the store contents
    // the passed-in domain is the valid parts of the store, everything
    // else inside is not valid any more
    virtual void invalidate(const entitySet& valid) {
      std::cerr << "storeRep.invalidate() is not implemented yet"
                << std::endl ;
      abort() ;
    }
    // this method is used to make sure that the storeRep
    // domain include the specified set, if not, reallocation
    // will be performed to make sure that. dynamic containers
    // actually don't need to do anything, so this is mainly
    // for static containers.
    virtual void guarantee_domain(const entitySet& include) {
      std::cerr << "storeRep.guarantee_domain() is not implemented yet"
                << std::endl ;
      abort() ;
    }
    virtual void shift(int_type offset) = 0 ;
    virtual void set_elem_size(int sz) ;
    virtual void setIsMat(bool im){};//for storeVecRep 
    virtual storeRep *new_store(const entitySet &p) const = 0 ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const = 0 ;
    // the remap method merely renumbers the container
    // according to the passed in map
    virtual storeRepP remap(const dMap &m) const = 0 ;
    // virtual storeRepP remap(const Map& m) const = 0 ;
    // the redistribute takes a vector of entitySets as domain
    // distribution over a group of processes and
    // redistributes the stores according to the domain partition
    // NOTE: should be pure virtual (= 0), but here we default it
    // to do nothing just to be able to compile the code, all
    // containers will need to implement this method later.
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) {
      std::cerr << "storeRep.redistribute() is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    // this redistribute version takes an additional remap
    // argument, upon redistribution, the new store is remapped
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) {
      std::cerr << "storeRep.redistribute() is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0);
    }
    // this redistribute version only remaps the new store domain
    // in the process (in case of maps, the image is not remapped)
    virtual storeRepP
    redistribute_omd(const std::vector<entitySet>& dom_ptn,
                     const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) {
      std::cerr << "storeRep.redistribute_omd() is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    // the freeze method converts a dynamic container to
    // a static one. For already static container,
    // its Rep() is returned    
    virtual storeRepP freeze() = 0 ;
    // the thaw method converts a static container to
    // a dynamic one. For already dynamic container,
    // its Rep() is returned ;
    virtual storeRepP thaw() = 0 ;
    // this version of freeze and thaw will take an entitySet
    // and creats a corresponding version of stores with
    // the passed in entitySet as its domain. in the meanwhile,
    // it also initializes the new store, for those domain
    // within the current store, the data will be copied;
    // for the domains outside of the current store, random
    // values picked from the current store will be used to
    // initialize them
    virtual storeRepP freeze(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    virtual storeRepP thaw(const entitySet& es) const {
      std::cerr << "storeRep.thaw(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    virtual void copy(storeRepP &st, const entitySet &context) = 0 ;
    virtual void fast_copy(storeRepP& st, const entitySet& context)
    { copy(st,context); }       // default behavior
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) = 0 ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) = 0 ;
    virtual int pack_size(const entitySet &e) = 0;
     virtual int estimated_pack_size(const entitySet &e) = 0;
    // this function also sets a domain that can actually be packed
    // in other words, the passed in domain "e" can contain
    // entities outside the store domain, the "packed" sets
    // the entities that can be packed.
    virtual int pack_size(const entitySet& e, entitySet& packed) = 0 ;
    virtual void pack(void *ptr, int &loc, int &size,  const entitySet &e) = 0 ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) = 0 ;
    // this version of pack/unpack uses a remap during the process
    // mainly for maps images to transform to another numbering scheme
    // default behavior is to ignore the remaps
    virtual void pack(void* ptr, int& loc,
                      int& size, const entitySet& e, const Map& remap) {
      pack(ptr,loc,size,e) ;
    }
    virtual void unpack(void* ptr, int& loc,
                        int& size, const sequence& seq, const dMap& remap) {
      unpack(ptr,loc,size,seq) ;
    }
    virtual store_type RepType() const = 0 ;
    virtual std::ostream &Print(std::ostream &s) const = 0 ;
    virtual std::istream &Input(std::istream &s) = 0 ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) = 0;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &en) const = 0;
    virtual entitySet domain() const = 0 ;
    virtual storeRepP getRep() ;
    virtual storeRepP getRep() const ;
    virtual  DatatypeP getType() = 0 ;
    virtual frame_info get_frame_info() = 0 ;
    virtual int getDomainKeySpace() const { return domainKeySpace ; }
    virtual void setDomainKeySpace(int v) { domainKeySpace = v ; }
  } ;


  class store_instance : public eventNotify {
    friend class store_ref ;
    storeRepP rep ;
  public:
    int getDomainKeySpace() const { return rep->getDomainKeySpace() ; }
    void setDomainKeySpace(int v) { rep->setDomainKeySpace(v) ; }
    enum instance_type { READ_WRITE, READ_ONLY } ;
    void setRep(const storeRepP &p)
      { rep = p ; rep.set_notify(this); notification() ; }
    storeRepP Rep() { return rep->getRep(); }
    storeRepP Rep() const { return rep->getRep(); }
    storeRepP getRep() {return rep->getRep() ;}
    storeRepP getRep() const {return rep->getRep() ;}
    virtual instance_type access() const ;
  } ;
  
  class store_ref : public store_instance, public storeRep {
  public:
    store_ref() {}
    store_ref(storeRepP &p) { setRep(p); }
    virtual ~store_ref() ;
    
    store_ref &operator=(storeRepP &p) {
      setRep(p) ;
      return *this ;
    }
        
    store_ref &operator=(const store_ref &sr) {
      setRep(sr.rep) ;
      return *this ;
    }

    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void fast_copy(storeRepP& st, const entitySet& context);
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) {
      Rep()->readhdf5(group_id, dataspace, dataset, dimension, name, fi, en);
    };
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const {
      Rep()->writehdf5(group_id, dataspace, dataset, dimension, name, en);
    };
    virtual entitySet domain() const ;
    virtual storeRepP getRep() ;
    virtual storeRepP getRep() const ;
    virtual void notification() ;
    virtual DatatypeP getType() {
     return  Rep()->getType() ;
    }
    virtual frame_info get_frame_info() {
      return Rep()->get_frame_info() ;
    }
    virtual void erase(const entitySet& rm) ;
    virtual void invalidate(const entitySet& valid) ;
    virtual void guarantee_domain(const entitySet& include) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute_omd(const std::vector<entitySet>& dom_ptn,
                     const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP freeze(const entitySet& es) const ;
    virtual storeRepP thaw(const entitySet& es) const ;
    virtual void pack(void* ptr, int& loc,
                      int& size, const entitySet& e, const Map& remap) ;
    virtual void unpack(void* ptr, int& loc,
                        int& size, const sequence& seq, const dMap& remap) ;
    virtual int getDomainKeySpace() const { return Rep()->getDomainKeySpace() ; }
    virtual void setDomainKeySpace(int v) { Rep()->setDomainKeySpace(v) ; }
    
  } ;
  typedef NPTR<store_ref> store_refP ;
    
}

#endif

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
#ifndef GSTORE_DEF_H
#define GSTORE_DEF_H

// This header file contains the class definition of
// GStore
// Their corresponding template implementation is in
// gstore_impl.h
// The same design applies to all other container classes.

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <data_traits.h>
#include <sstream>
#include <iterator>
#include <vector>


namespace Loci {
  enum gstore_type {GSTORE, GSTOREVEC, GMULTISTORE, GMAP, GMULTIMAP, GUNKNOWN} ;
  typedef int  GEntity;   
  typedef genIntervalSet<GEntity> GEntitySet; 
 
  template<class T> class GStore {
   
  public:
    //storage type
    typedef typename std::vector<std::pair<GEntity,T> > gRep ;

    //iterators, use std iterators so that std functions can be used
    typedef typename std::vector<std::pair<GEntity,T> >::iterator iterator ;
    typedef typename std::vector<std::pair<GEntity,T> >::const_iterator const_iterator ;
    iterator begin(){return Rep->begin();}
    iterator end(){return Rep->end();}
    const_iterator begin() const {return Rep->begin();}
    const_iterator end() const {return Rep->end();}
   
    

    //constructors
    //GStore() : gtype(GUNKNOWN),sorted(false) {} 
    GStore(gstore_type t) : gtype(t),sorted(false) {} 

    //copy constructor, shallow copy
    GStore(const GStore& s):Rep(s.Rep),gtype(s.gtype),sorted(s.sorted){} 

    //output
    std::ostream& Print(std::ostream &s) const;
    
    //inspectors
    gstore_type type() const {return gtype;}
    size_t size() const{return Rep->size();}

    //insert elements into store
    void insert(GEntity e, T val){Rep->push_back(std::pair<GEntity,T>(e,val));}
    void insert(const GEntitySet& seq,  const T* vals){
      size_t idx= 0;
      for(typename GEntitySet::const_iterator itr = seq.begin();
          itr!= seq.end(); itr++){
        Rep->push_back(std::pair<GEntity, T>(*itr, vals[idx++]));
      }
    }

    
    void local_sort();
    void sort();
    
   
    void make_consistent();
   
   
    //binary search function, for random access
    iterator search(iterator first,
                    iterator last,
                    const GEntity& key){
      // if(!sorted)throw(bad_access);
      if(first->first==key) return first;
      iterator it;
     
      typename  std::iterator_traits<iterator>::difference_type count, step;
      count = distance(first,last);
      while (count>0)
        {
          it = first; step=count/2; advance(it,step);
          if (it->first < key) {                
            first=++it;
            count-=step+1;
          }
          else count=step;
        }
      return first;
    }
    
    //binary search function with const access, for random access
    const_iterator search(const_iterator first,
                          const_iterator last,
                          const GEntity& key)const{
      // if(!sorted)throw(bad_access);
      if(first->first==key) return first;
      const_iterator it;
     
      typename  std::iterator_traits<const_iterator>::difference_type count, step;
      count = distance(first,last);
      while (count>0)
        {
          it = first; step=count/2; advance(it,step);
          if (it->first < key) {                
            first=++it;
            count-=step+1;
          }
          else count=step;
        }
      return first;
    }
    
 
    
      
    // std::pair<iterator, iterator>
    //     equal_range(iterator first, iterator last, const GEntity& key){
    //       iterator found = search(first, last, key);
    //       if(found == last){
    //         debugout << "key " << key << " not found " << endl; 
    //         return -1;
    //       }
    //       iterator next = found;
    //       while(next != last && next->first == key){
    //         next++;
    //       }
    //       return --next;
    //     }
 
    
    //copy gstore to traditional Loci store, storevec or multistore
    //partially implemented
    storeRepP copy2store(const entitySet& dom) const{
      // if(!sorted)throw(bad_access);
      if(gtype == GSTORE){
        fatal(size() != dom.size()) ;
        store<T> s;
        s.allocate(dom);
        const_iterator itr = begin();  
        FORALL(dom,i) {
          s[i] = itr->second;
          itr++;
        } ENDFORALL ;
        return s.Rep();
      }else{
        debugout<<" gstore_copy cannot deal with the type " << gtype << endl;
        return NULL; 
      }
    }
    
    //copy gstore to traditional Loci Map or multiMap
    storeRepP copy2map(const entitySet& dom) const{
      // if(!sorted)throw(bad_access);
      if(gtype == GMAP){
        fatal(size() != dom.size()) ;
        Map s;
        s.allocate(dom);
        const_iterator itr = begin();  
        FORALL(dom,i) {
          s[i] = int(itr->second);
          itr++;
        } ENDFORALL ;
        return s.Rep();
      }else if(gtype == GMULTIMAP){
        store<int> counts;
        counts.allocate(dom);
        const_iterator itr = begin();
        const_iterator last = end();
        FORALL(dom,i) {
          const_iterator first = itr;
          int cnt = 0;
          while(itr != last && itr->first == first->first){
            itr++;
            cnt++;
          }
          counts[i] = cnt;
        }ENDFORALL ;
        multiMap s;
        s.allocate(counts);
        itr = begin(); 
        FORALL(dom,i) {
          for(int j = 0; j < counts[i]; j++){
            s[i][j] = int(itr->second);
            itr++;
          }
        } ENDFORALL ;
        return s.Rep();
      }else{
        debugout<<" gstore_copy cannot deal with the type " << gtype << endl;
        return NULL; 
      }
    }
    
 
    
  private:
    Handle<gRep> Rep ;
    gstore_type gtype;
    bool sorted;
  } ;
   
} // end of namespace Loci

#endif

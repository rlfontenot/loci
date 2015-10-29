//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#ifndef GPARAMETER_H
#define GPARAMETER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <mpi.h>
#include <Tools/debug.h>
#include <gstore_rep.h>
#include <distribute.h>
#include <data_traits.h>
#include <distribute_long.h>
#include <partition.h>
#include <parameter.h>
namespace Loci {
  
  template<class T> class gParamRepI : public gStoreRep {
    gEntitySet store_domain ;
    T attrib_data ;
    gKeySpace* domain_space;
  private:
    int get_mpi_size( IDENTITY_CONVERTER c, const gEntitySet &eset)const;
    int get_mpi_size( USER_DEFINED_CONVERTER c, const gEntitySet &eset)const;
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size )const;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size )const;
    void unpackdata(IDENTITY_CONVERTER c,     const void *ptr, int &loc, int size);
    void unpackdata(USER_DEFINED_CONVERTER c, const void *ptr, int &loc, int size);
    DatatypeP getType(IDENTITY_CONVERTER g)const ;
    DatatypeP getType(USER_DEFINED_CONVERTER g)const ;
   
  public:

    gParamRepI():store_domain(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX)),domain_space(0){}
    gParamRepI(const gEntitySet &p):store_domain(p),domain_space(0){}
    virtual void allocate(const gEntitySet &p){  store_domain = p ; }
    virtual ~gParamRepI() {}
    void set_domain_space(gKeySpace* space){domain_space = space;}
    gKeySpace* get_domain_space()const{return domain_space;}
    // this method redistributes the container according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const;

    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const;

    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const;
    //  virtual gStoreRepP
    //         redistribute_omd(const std::vector<gEntitySet>& dom_ptn,
    //                          const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD)  ;

    virtual void shift(int_type offset) ;
    virtual gstore_type RepType() const { return GPARAMETER ;}
    virtual gEntitySet domain() const { return store_domain ;}
    virtual gStoreRepP remap(const gMap &m) const;
    virtual gStoreRepP clone() const;
    virtual void* get_attrib_data() {
      return static_cast<void*>(&attrib_data);
    }
    virtual const void* get_attrib_data() const{
      return  static_cast<const void*>(&attrib_data);
    }
    virtual storeRepP copy2store()const{
      param<T> r ;
      r.set_entitySet(store_domain);
      *r = attrib_data;
      return r.Rep() ;
    }
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;//only pack attrib_data
    virtual void unpack(const void *ptr, int &loc, int size)  ;//only unpack attribute_data

    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual DatatypeP getType() const ;
    
  } ;

 
  //**************************************************************************/


  template<class T> class gParam : public gstore_instance {
    typedef gParamRepI<T> gParamType ;
    T * data ;
  public:
    typedef T containerType ;
    gParam() { setRep(new gParamType) ;
      data = static_cast<T*> (Rep()->get_attrib_data());
    }
    
    //shallow copy
    gParam(const gParam &var) {
      setRep(var.Rep());
      data = static_cast<T*> (Rep()->get_attrib_data()); 
    }
    
    gParam(gStoreRepP rp) {
      setRep(rp);
      data = static_cast<T*> (Rep()->get_attrib_data());   
    }

    virtual ~gParam() {}
    //shallow assignment
    gParam & operator=(const gParam &p) {
      setRep(p.Rep());
      data = static_cast<T*> (Rep()->get_attrib_data()); 
      return *this ;
    }
    
    gParam & operator=(gStoreRepP p) {
      setRep(p) ;
      data = static_cast<T*> (Rep()->get_attrib_data()); 
      return *this ;
    }

    gParam & operator=(const T &v) {
      *data = v ;
      return *this ;
    }
    
    T * operator->() { return data ; }

    const T * operator->() const { return data ; }

    T * operator&() { return data ; }

    const T * operator &() const { return data ; }

    T &operator*() { return *data ; }

    const T &operator*() const { return *data ; }

    void set_domain_space(gKeySpace* space){Rep()->set_domain_space(space);}
    gKeySpace* get_domain_space()const{return Rep()->get_domain_space();}
    
    void set_entitySet(const gEntitySet &ptn) {
      CPTR<gParamType> p(Rep()) ;
      if(p != 0) p->allocate(ptn);
      warn(p==0) ; 
    }

    gEntitySet domain() const { return Rep()->domain(); }

    // this method redistributes the container according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm);}

    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split,remap,comm);}


    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->split_redistribute(dom_ptn, comm);}
    

    // the remap method merely renumbers the container
    // according to the passed in map
    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);} 
    virtual gStoreRepP clone() const{return Rep()->clone();}
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;
  
  

  //**************************************************************************/

  template<class T> void gParamRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
  }
    
  //**************************************************************************/

  template<class T>
  std::ostream &gParamRepI<T>::Print(std::ostream &s) const
  {
    gEntitySet dom = domain() ;
    if(dom == ~EMPTY) {
      Loci::streamoutput(&attrib_data,1,s) ;
    } else {
      s << '{' << domain() << std::endl ;
      Loci::streamoutput(&attrib_data,1,s) ;
      s << '}' << std::endl ;
    }
    return s ;
  }

  //**************************************************************************/

  template<class T>
  std::istream &gParamRepI<T>::Input(std::istream &s)
  {
    gEntitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      s.putback(ch) ;
      e = ~EMPTY ;
      allocate(e) ;
      attrib_data = T() ;
      Loci::streaminput(&attrib_data,1,s) ;
      return s ;
    }

    s >> e ;
    allocate(e) ;

    attrib_data = T() ;
    Loci::streaminput(&attrib_data,1,s) ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading gParameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }
  
  //**************************************************************************/

  template<class T>
  DatatypeP gParamRepI<T>::getType()const {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  
  //**************************************************************************/ 

  template<class T>
  DatatypeP gParamRepI<T>::getType(IDENTITY_CONVERTER g)const {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  
  //**************************************************************************/

  template<class T>
  DatatypeP gParamRepI<T>::getType(USER_DEFINED_CONVERTER g)const {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }

  //**************************************************************************/

  template <class T>
  int gParamRepI<T>::pack_size( const gEntitySet &eset)const
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }
 
  //**************************************************************************/

  template <class T>
  int gParamRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const gEntitySet &eset)const
  {

    return( sizeof(T) ) ;
  }
  
  //   //**************************************************************************/

  template <class T>
  int gParamRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const gEntitySet &eset)const
  {
    typedef data_schema_traits<T> schema_traits;

    typename schema_traits::Converter_Type cvtr(const_cast<T&>(attrib_data));
    int arraySize = cvtr.getSize() ;

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) + sizeof(int));
  }
  
  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::pack(void *ptr, int &loc, int size, const gEntitySet &e )const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    packdata( schema_converter(), ptr, loc, size);
  }

  //**************************************************************************/
  template <class T>
  void gParamRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                int outcount )const
  {
    MPI_Pack( const_cast<T*>(&attrib_data), sizeof(T), MPI_BYTE, outbuf, outcount, &position,
              MPI_COMM_WORLD) ;
  }

  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                                int &position, int outcount )const
  {
    int stateSize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    typename data_schema_traits<T>::Converter_Type cvtr(const_cast<T&>(attrib_data));

    std::vector<dtype> inbuf(cvtr.getSize());
    cvtr.getState( &inbuf[0], stateSize);

    MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
             MPI_COMM_WORLD);
    int incount =  stateSize*typesize;
    MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position,
             MPI_COMM_WORLD) ;

  }
  
  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::unpack(const void *ptr, int &loc, int size)  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata( schema_converter(), ptr, loc, size);
  }
  
  //**************************************************************************/
 
  template <class T>
  gStoreRepP gParamRepI<T>::remap(const gMap &m) const{
    gParam<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    r.set_domain_space(domain_space);
    return r.Rep() ;
  }
  
  //**************************************************************************/
 
  template <class T>
  gStoreRepP gParamRepI<T>::clone() const{
    gParam<T> result;
    *result = attrib_data;
    result.set_entitySet(store_domain);
    result.set_domain_space(domain_space);
    return result.Rep();
  }

  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::unpackdata( IDENTITY_CONVERTER c, const void *inbuf, int &position,
                                  int insize)
  {

    
    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &attrib_data, sizeof(T),
                MPI_BYTE, MPI_COMM_WORLD) ;

  }

  //***********************************************************************/

  template <class T>
  void gParamRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, const void *inbuf,
                                  int &position, int insize)
  {

    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &stateSize, 1,
                MPI_INT, MPI_COMM_WORLD) ;
    std::vector<dtype> outbuf(stateSize);

    outcount = stateSize*sizeof(dtype);
    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &outbuf[0], outcount,
                MPI_BYTE, MPI_COMM_WORLD) ;
    typename schema_traits::Converter_Type  cvtr( attrib_data );
    cvtr.setState( &outbuf[0], stateSize);

  }

  
  //**************************************************************************/

  // template<class T>
  // inline std::ostream & operator<<(std::ostream &s, const const_gparam<T> &t)
  // { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> gStoreRepP gParamRepI<T>::
  split_redistribute(const std::vector<gEntitySet>& dom_ptn, MPI_Comm comm)const {
    // for a gParameter, we just redistribute its domain
    gEntitySet dom = domain() ;
    gEntitySet new_all ;
    for(size_t i=0;i<dom_ptn.size();++i)
      new_all += dom_ptn[i] ;
    gEntitySet out = dom - new_all ;

    std::vector<gEntitySet> old_dist = g_all_collect_vectors<gEntity>(dom, comm) ;
    gEntitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    gEntitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom += out ;

    gParam<T> np ;
    np.set_entitySet(new_dom) ;
    *np = attrib_data ;
    np.set_domain_space(domain_space);
    return np.Rep() ;
  }
  
  //**************************************************************************/

  template<class T> gStoreRepP gParamRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               MPI_Comm comm)const{
    if(store_domain == gEntitySet(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX))){
      gParam<T> np ;
      np.set_entitySet(store_domain) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return np.Rep() ;
    }else{
      std::vector<gEntitySet> recv_split = transposePtn(dom_split, comm);
      gEntitySet new_dom;
      for(size_t i=0;i<recv_split.size();++i)new_dom += recv_split[i];
      gEntitySet old_all = g_all_collect_entitySet<gEntity>(store_domain, comm);
      new_dom = new_dom&old_all;
      gParam<T> np ;
      np.set_entitySet(new_dom) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return  np.Rep() ;
    }
  }
  
  //**************************************************************************/
  
  template<class T> gStoreRepP gParamRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               const gMap& remap, MPI_Comm comm)const{
   
    if(store_domain ==  gEntitySet(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX))){
      gParam<T> np ;
      np.set_entitySet(store_domain) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return np.Rep() ;
    }else{
      gParam<T> np ;
      np = redistribute(dom_split, comm);
      np.set_domain_space(domain_space);
      return np.remap(remap);
    } 
  }

  //**************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const gParam<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T>
  inline std::istream & operator>>(std::istream &s, gParam<T> &t)
  { return t.Input(s) ; }

  
  
}

#endif

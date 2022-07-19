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
#ifndef STORE_IMPL_H
#define STORE_IMPL_H

#include <store_def.h>
#include <DStore_def.h>

namespace Loci {
  template<class T> void storeRepI<T>::allocate(const entitySet &eset) {
    // if the pass in is EMPTY, we delete the previous allocated memory
    // this equals to free the memory
    if( eset == EMPTY ) {
#ifdef PAGE_ALLOCATE
      // Explicitly call destructor
      FORALL(store_domain,ii) {
	base_ptr[ii].~T() ;
      } ENDFORALL ;
      pageRelease(store_domain.Max()-store_domain.Min()+1,alloc_pointer) ;
#else
      delete [] alloc_pointer ;
#endif
      alloc_pointer = 0 ; base_ptr = 0;
      store_domain = eset ;
      dispatch_notify() ;
      return ;
    }

    int_type old_range_min = store_domain.Min() ;
    int_type old_range_max = store_domain.Max() ;
    int_type new_range_min = eset.Min() ;
    int_type new_range_max = eset.Max() ;

    // if the old range and the new range are equal, nothing
    // needs to be done, just return
    if( (old_range_min == new_range_min) &&
        (old_range_max == new_range_max)) {
      store_domain = eset ;
      return ;
    }

    // is there any overlap between the old and the new domain?
    // we copy the contents in the overlap region to the new
    // allocated storage
    entitySet ecommon = store_domain & eset ;

#ifdef PAGE_ALLOCATE
    // Allocate memory in whole pages
    T* tmp_alloc_pointer = 0 ;
    tmp_alloc_pointer = pageAlloc(new_range_max-new_range_min+1,
				  tmp_alloc_pointer) ;
    // Now call constructor
    T* tmp_base_ptr = tmp_alloc_pointer - new_range_min ;
    // Call placement new
    FORALL(eset,ii) {
      new(&tmp_base_ptr[ii]) T() ;
    } ENDFORALL ;
#else
    T* tmp_alloc_pointer = new T[new_range_max - new_range_min + 1] ;
    T* tmp_base_ptr = tmp_alloc_pointer - new_range_min ;
#endif
    // if ecommon == EMPTY, then nothing is done in the loop
    FORALL(ecommon,i) {
      tmp_base_ptr[i] = base_ptr[i] ;
    } ENDFORALL ;


#ifdef PAGE_ALLOCATE
    // Explicitly call destructor
    FORALL(store_domain,ii) {
      base_ptr[ii].~T() ;
    } ENDFORALL ;
    pageRelease(store_domain.Max()-store_domain.Min()+1,alloc_pointer) ;
#else
    delete [] alloc_pointer ;
#endif
    alloc_pointer = tmp_alloc_pointer ;
    base_ptr = tmp_base_ptr ;

    store_domain = eset ;
    dispatch_notify() ;
    return ;
  }

  template<class T> void storeRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    base_ptr -= offset ;
    dispatch_notify() ;
  }

  template<class T> std::ostream &storeRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      Loci::streamoutput(&base_ptr[ii],1,s) ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  template<class T> std::istream &storeRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      base_ptr[ii] = T() ;
      Loci::streaminput(&base_ptr[ii],1,s) ;
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //*************************************************************************/

  template<class T>  storeRepI<T>::~storeRepI() {
    if(alloc_pointer) {
#ifdef PAGE_ALLOCATE
      // Explicitly call destructor
      FORALL(store_domain,ii) {
	base_ptr[ii].~T() ;
      } ENDFORALL ;
      
      pageRelease(store_domain.Max()-store_domain.Min()+1,alloc_pointer) ;
#else
      delete[] alloc_pointer ;
#endif
    }
  }

  template<class T>  entitySet storeRepI<T>::domain() const {
    return store_domain ;
  }

  template<class T>
  storeRep *storeRepI<T>::new_store(const entitySet &p) const {
    return new storeRepI<T>(p)  ;
  }
  template<class T>
  storeRep *storeRepI<T>::new_store(const entitySet &p, const int* cnt) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a store " << endl ;
    return sp ;
  }
  template<class T> store_type storeRepI<T>::RepType() const {
    return STORE ;
  }

  template<class T> void store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const store<T> &t)
  { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, store<T> &t)
  { return t.Input(s) ; }

  template<class T> void const_store<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
  const_store<T>::access() const { return READ_ONLY; }

  template<class T>
  storeRepP storeRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    store<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;

    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T> storeRepP storeRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T> storeRepP storeRepI<T>::thaw() {
    dstore<T> ds ;
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei)
      ds[*ei] = base_ptr[*ei] ;
    return ds.Rep() ;
  }

  template<class T> storeRepP
  storeRepI<T>::thaw(const entitySet& es) const {
    entitySet shared = domain() & es ;
    entitySet out = es - domain() ;

    dstore<T> ds ;
    for(entitySet::const_iterator ei=shared.begin();
        ei!=shared.end();++ei)
      ds[*ei] = base_ptr[*ei] ;

    Entity c = *domain().begin() ;
    for(entitySet::const_iterator ei=out.begin();
        ei!=out.end();++ei)
      ds[*ei] = base_ptr[c] ;
    
    return ds.Rep() ;
  }
  
  template<class T> void storeRepI<T>::copy(storeRepP &st, const entitySet &context)  {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr ==0)) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }
  
  template<class T> void storeRepI<T>::gather(const dMap &m, storeRepP &st,
                                              const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T> void storeRepI<T>::scatter(const dMap &m, storeRepP &st,
                                               const entitySet &context) {
    const_store<T> s(st) ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }

  //*******************************************************************/

  template <class T>
  inline int storeRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                         const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }
  //*******************************************************************/
  template <class T>
  inline int storeRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c,
                                                   const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }
  //*******************************************************************/
  template <class T>
  inline int storeRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                         const entitySet &eset)
  {
    
    int       size ;
    int numBytes = 0 ;
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    fatal((eset - domain()) != EMPTY);

    entitySet sdom = eset & domain() ;
    for( ci = sdom.begin(); ci != sdom.end(); ++ci) {
      typename converter_traits::Converter_Type cvtr(base_ptr[*ci]);
      size      = cvtr.getSize();
      numBytes += size*sizeof(typename converter_traits::Converter_Base_Type) ;
    }

    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
  }
  //*******************************************************************/
  template <class T>
  inline int storeRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c,
                                                   const entitySet &eset)
  {
    int numBytes = 0 ;
    numBytes = eset.size()*50*sizeof(double);
    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
  }

  //*******************************************************************/
  template <class T> int storeRepI<T>::pack_size( const entitySet &eset) {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }

  //*******************************************************************/
  template <class T> int storeRepI<T>::estimated_pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    return get_estimated_mpi_size( schema_converter(), eset );
    
  }
  

  template<class T> int storeRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    typedef typename data_schema_traits<T>::Schema_Converter
      schema_converter ;

    return get_mpi_size(schema_converter(), packed) ;
  }

  //*******************************************************************/
  template <class T>
  inline void storeRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                                     int &position, int outcount,
                                     const entitySet &eset )
  {
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      MPI_Pack( &base_ptr[begin], t*sizeof(T), MPI_BYTE, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
    }
  }

  //*******************************************************************/

  template <class T>
  inline void storeRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                                      int &position, int outcount,
                                      const entitySet &eset )
  {

    entitySet :: const_iterator ci;
    entitySet  ecommon;

    ecommon = domain()&eset;
    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------
    size_t  incount = 0;
    int     stateSize, maxStateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
      stateSize    = cvtr.getSize();
      incount     += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
    }

    typename schema_traits::Converter_Base_Type *inbuf;

    int typesize = sizeof(typename schema_traits::Converter_Base_Type);

    inbuf = new typename schema_traits::Converter_Base_Type[maxStateSize];

    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------

    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci]);
      cvtr.getState( inbuf, stateSize);

      incount =  sizeof(int);
      MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
               MPI_COMM_WORLD);

      incount =  stateSize*typesize;
      MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position,
               MPI_COMM_WORLD) ;
    }
    delete [] inbuf;
  }
  //*******************************************************************/

  template <class T>
  void storeRepI<T>::pack( void *outbuf, int &position, int &size,
                           const entitySet &usr_eset )
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    warn(usr_eset-domain() != EMPTY) ;

    packdata( schema_converter(), outbuf, position, size, usr_eset);

  }

  //*********************************************************************/
  template <class T>
  inline void storeRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                       int &position,  int insize,
                                       const sequence &seq)
  {

    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx)
          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                      sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
                    t*sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      }
    }
  }
  //*********************************************************************/
  template <class T>
  inline void storeRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf,
                                        int &position, int insize, const sequence &seq)
  {

    sequence :: const_iterator ci;

    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------
    int  stateSize, outcount;

    typedef  data_schema_traits<T> converter_traits;
    typename converter_traits::Converter_Base_Type *outbuf;

    int typesize = sizeof(typename converter_traits::Converter_Base_Type);

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( !store_domain.inSet( *ci ) ) {
        std::cout << "Warning: entity not in entitySet  : "
                  << *ci << std::endl;
        continue;
      }
      outcount = sizeof(int);
      MPI_Unpack(inbuf, insize, &position, &stateSize, 1,
                 MPI_INT, MPI_COMM_WORLD) ;

      outbuf = new typename converter_traits::Converter_Base_Type[stateSize];

      outcount = stateSize*typesize;
      MPI_Unpack(inbuf, insize, &position, outbuf, outcount,
                 MPI_BYTE, MPI_COMM_WORLD) ;

      typename converter_traits::Converter_Type cvtr( base_ptr[*ci] );
      cvtr.setState( outbuf, stateSize);
      delete [] outbuf;
    }

  }


  //*******************************************************************/

  template <class T>
  void storeRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata(schema_converter(), ptr, loc, size, seq);
  }

  //**********************************************************************/

  template<class T>
  frame_info storeRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  
  template<class T>
  frame_info storeRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }

  template<class T>
  frame_info storeRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {

    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci)
      for(int ivec = 0; ivec < fi.size; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)*fi.size+ivec] );
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
    return fi ;
  }

  
  template<class T>
  DatatypeP storeRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T>
  DatatypeP storeRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T>
  DatatypeP storeRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/

  template<class T>
  void storeRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }


#ifdef H5_HAVE_PARALLEL
  template<class T>
  void storeRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset, hid_t xfer_plist_id) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5writeP(group_id, dataspace, dataset, dimension, name, traits_output_type, eset, xfer_plist_id) ;
  }
#endif
  
  //**********************************************************************/
  template <class T>
  void storeRepI<T> :: hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset)  const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
        tmp_array[tmp++] = base_ptr[*si] ;


      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template <class T>
  void storeRepI<T> :: hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset, hid_t xfer_plist_id)  const
  {
    //   if(dimension != 0) {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
      tmp_array[tmp++] = base_ptr[*si] ;


    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
    //    }
  }
#endif  
  //*********************************************************************/

  template <class T>
  void storeRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
  {
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      size_t tmp = 0 ;
      int stateSize = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        typename schema_traits::Converter_Type cvtr(base_ptr[*si]);
        cvtr.getState(tmp_array+tmp, stateSize) ;
        tmp +=stateSize ;
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template <class T>
  void storeRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id) const
  {
    typedef data_schema_traits<T> schema_traits ;
    //  if(dimension != 0) {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype* tmp_array = new dtype[dimension] ;
    size_t tmp = 0 ;
    int stateSize = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
      typename schema_traits::Converter_Type cvtr(base_ptr[*si]);
      cvtr.getState(tmp_array+tmp, stateSize) ;
      tmp +=stateSize ;
    }
    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
    //}
  }
#endif
  //*********************************************************************/
  template<class T>
  void storeRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }
#ifdef H5_HAVE_PARALLEL
  template<class T>
  void storeRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5readP(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset, xfer_plist_id) ;
  }
#endif
  
  template<class T>
  void  storeRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        std::string error = "H5Dread: failed" ;
        throw StringError(error) ;
      }
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
        base_ptr[*si] = tmp_array[tmp++] ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }
#ifdef H5_HAVE_PARALLEL
  template<class T>
  void  storeRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    // if(dimension != 0) {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                        xfer_plist_id, tmp_array) ;
    if(err < 0) {
      std::string error = "H5Dread: failed" ;
      throw StringError(error) ;
    }
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
      base_ptr[*si] = tmp_array[tmp++] ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
    //}
  }
#endif
  //*********************************************************************/

  template<class T>
  void  storeRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      std::vector<int> vint = fi.second_level ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			  H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
	std::string error = "H5Dread: failed" ;
	throw StringError(error) ;
      }
      size_t tmp = 0 ;
      int bucsize ;
      size_t indx = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
	typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si]);
	bucsize = vint[indx++] ;
	cvtr.setState(tmp_array+tmp, bucsize) ;
	tmp += bucsize ;
      }
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void  storeRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset,hid_t xfer_plist_id)
  {
    typedef data_schema_traits<T> schema_traits ;
    //   if(dimension != 0) {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    std::vector<int> vint = fi.second_level ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype* tmp_array = new dtype[dimension] ;
    hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                        xfer_plist_id, tmp_array) ;
    if(err < 0) {
      std::string error = "H5Dread: failed" ;
      throw StringError(error) ;
    }
    size_t tmp = 0 ;
    int bucsize ;
    size_t indx = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
      typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si]);
      bucsize = vint[indx++] ;
      cvtr.setState(tmp_array+tmp, bucsize) ;
      tmp += bucsize ;
    }
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
    // }
  }
#endif
  //*********************************************************************/

} // end of namespace Loci

#endif

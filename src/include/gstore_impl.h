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

#ifndef GSTORE_IMPL_H
#define GSTORE_IMPL_H 1
#include <data_traits.h>
#include <Tools/except.h>
#include <parSampleSort.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <gstore_rep.h>
#include <hdf5_readwrite_long.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>
#include <gmap.h>
#include <distribute_long.h>
using std::vector;

namespace Loci {
  
  
  namespace {
    template<class T> inline  bool gstore_porder(const std::pair<gEntity, T> &i1, 
                                                 const std::pair<gEntity, T> &i2) {
      return i1.first < i2.first ;
    }
  }


  template<class T> void gStore<T>::make_consistent(){
    debugout<< "WARNING: void gStore<T>::make_consistent() not implemented yet" << endl;
  }
  
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const gStore<T> &t)
  { return t.Print(s) ; }


  template<class T> void gStoreRepI<T>::local_sort(){
    std::stable_sort(attrib_data.begin(), attrib_data.end(), gstore_porder<T>);
    sorted = true;
  }

  // //parSampleSort not working, because balanceDisttribution will split one entity inside a multiMap
  // // or a multiStore 
  // template<class T> void gStoreRepI<T>::sort(){
  //   cerr<<" WARNING: gStoreRepI<T>::sort() not implemented yet" << endl; 
  // }

  // template<class T> void gStoreRepI<T>::parSplitSort(const std::vector<gEntitySet>& ptn){
  //   cerr<<" WARNING:gStoreRepI<T>::parSplitSort not implemented yet" << endl; 
  // }
 
                         


  template<class T>
  void gStoreRepI<T>::shift(gEntity offset) {
    for(iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){
      itr->first >>= offset ;
    }
  }
  
  template<class T>
  gStoreRepP gStoreRepI<T>::remap(const gMap &m) const{
    fatal(!sorted);
    fatal( !m.sorted());
    fatal(m.domain()-domain() != GEMPTY);
    gStore<T> s;
    typename gRep::const_iterator itr1 = begin();
    gMap::const_iterator itr2 = m.begin();
    for(; itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) s.insert(itr2->second, itr1->second);
    }
    s.local_sort();
    s.set_domain_space(domain_space);
    return s.Rep();
  }

  //*********************************************************************/

  template<class T> 
  std::ostream &gStoreRepI<T>::Print(std::ostream &s) const 
  {
    if(attrib_data.empty()){
      s <<"{}" << endl ;
      return s ;
    }
    s << '{' << std::endl ;
    typename gRep::const_iterator itr;
    for( itr = attrib_data.begin();
         itr != attrib_data.end(); itr++){
      s << itr->first << " -> " << itr->second << std::endl ;
    }
    
    s << '}' << std::endl ; 
    return s ;
  }

  //************************************************************************/

  template<class T> 
  std::istream &gStoreRepI<T>::Input(std::istream &s) 
  {
    debugout<< "WARNING: gStoreRepI<T>::Input() not implemented yet " << endl;
    // entitySet e ;
    //     char ch ;

    //     do ch = s.get(); while(ch==' ' || ch=='\n') ;

    //     if(ch != '{') {
    //       cerr << "Incorrect Format while reading store" << std::endl ;
    //       s.putback(ch) ;
    //       return s ;
    //     }

    //     s >> e ;
        
    //     FORALL(e,ii) {
    //       Loci::streaminput(&attrib_data.access(ii),1,s) ;
    //     } ENDFORALL ;
        
    //     do ch = s.get(); while(ch==' ' || ch=='\n') ;

    //     if(ch != '}') {
    //       cerr << "Incorrect Format while reading store" << std::endl ;
    //       s.putback(ch) ;
    //     }

    return s ;
  }

  
  

  //*************************************************************************/
  
  template<class T> 
  gstore_type gStoreRepI<T>::RepType() const 
  {
    return GSTORE ;
  }

  

  
  //For stores,recompose will compose a store whose domain is the domain of m,
  //whose data is the data of the SECOND field of m.
  //for example, pos.recompose(face2node) will produce the positions  for each face  
  template<class T>
  gStoreRepP  gStoreRepI<T>::recompose(gStoreRepP& m, MPI_Comm comm)const{
    gMap remap(m);
    
    gEntitySet new_dom = remap.image();
    gEntitySet old_dom = domain();
    vector<gEntitySet> ptn = g_all_collect_vectors<gEntity>(old_dom);
    gStore<T> expanded_store;//the value in expanded store is unique in each process
    expanded_store = split_redistribute(ptn, comm);
    expanded_store.local_sort();
    remap.local_sort2();

    gStore<T> result;
    gMap::const_iterator m_itr = remap.begin();
    typename gStore<T>::iterator s_itr = expanded_store.begin();
    while(m_itr!= remap.end() && s_itr != expanded_store.end()){
      while(m_itr!= remap.end() && m_itr->second == s_itr->first){
        result.insert(m_itr->first,s_itr->second);
        m_itr++;
      }
      s_itr++;
    }
    remap.local_sort();
    result.set_domain_space(domain_space);
    return result.Rep();
  }


  
  template<class T>
  gStoreRepP gStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
                                         MPI_Comm comm)const {
    fatal(!sorted);
    int np = 1;
    MPI_Comm_size(comm, &np);
    // compute the send_counts
    std::vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i){
      send_counts[i] = pack_size(dom_split[i]) ;
    }
    // communicate to get the receive counts
    std::vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;
    
    // then pack the buffer
    int tot_send_size = 0 ;
    int tot_recv_size = 0 ;
    for(int i=0;i<np;++i){
      tot_send_size += send_counts[i] ;
      tot_recv_size += recv_counts[i] ;
    }
    
    
    int *send_displs = new int[MPI_processes] ;
    int *recv_displs = new int[MPI_processes] ;
    send_displs[0] = 0 ;
    recv_displs[0] = 0 ;
    for(int i=1;i<np;++i){
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    }
    
    unsigned char* send_buffer = new unsigned char[tot_send_size] ;
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      pack(send_buffer, loc_pack, tot_send_size, dom_split[i]) ;
    
    // we now proceed to communicate the store contents
    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buffer[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, comm) ;

    
    gStore<T> ns ;
    gStoreRepP srp = ns.Rep() ;
    // unpack
    loc_pack = 0 ;
    srp->unpack(recv_buffer, loc_pack, tot_recv_size) ;
    srp->local_sort();
    srp->set_domain_space(domain_space);
    
    delete[] recv_displs ;
    delete[] send_displs ;
    delete[] send_buffer ;
    delete[] recv_buffer ;
    return srp ;
  }
  
  
  //split domain first, then redistribute
  template<class T>
  gStoreRepP gStoreRepI<T>::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                                               MPI_Comm comm)const {
    fatal(!sorted);

    gEntitySet my_dom = domain() ;
    // figure out how the domain is split to send to others
    int np  = 1;
    MPI_Comm_size(comm, &np) ;
    fatal(np != dom_ptn.size()) ;
    std::vector<gEntitySet> dom_split(np) ;
    for(int i=0;i<np;++i) {
      dom_split[i] = my_dom & dom_ptn[i] ;
    }
    return redistribute(dom_split, comm);
  }
  
  template<class T> gStoreRepP
  gStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
                              const gMap& remap, MPI_Comm comm)const {
    fatal(!sorted);

    gStore<T> result;
    result.setRep(redistribute( dom_split,comm));
    result.local_sort();
    result.set_domain_space(domain_space);
    return result.remap(remap);
  }

  //   template<class T> storeRepP
  //   gStoreRepI<T>::redistribute_omd(const std::vector<gEntitySet>& dom_ptn,
  //                                   const dMap& remap, MPI_Comm comm) {
  //   //   // this is a push operation, thus the send, recv are reversed
  // //     std::vector<P2pCommInfo> send, recv ;
  // //     get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
  // //     gStore<T> new_store ;
  // //     fill_store_omd(getRep(), 0, new_store.Rep(), &remap, send, recv, comm) ;
  // //     return new_store.Rep() ;
  //   }

 

  //*************************************************************************/
  template <class T> 
  int gStoreRepI<T>::pack_size( const gEntitySet &eset)const {
   
    fatal(eset - domain() != GEMPTY);
   
    fatal(!sorted);
   
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
   
    return get_mpi_size( traits_type, eset );
  }
  
  

  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  int gStoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                   const gEntitySet &eset)const{
    int result = 0;
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      while(itr1 != attrib_data.end() && itr1->first < ci  )itr1++;
      while(itr1 != attrib_data.end()&& itr1->first == ci ){
        result += (sizeof(gEntity)+sizeof(T));
        itr1++;
      }
    }ENDGFORALL;
    return result;
  }
  

  
  
  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  int gStoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                   const gEntitySet &eset)const
  {
   
    int result = 0;
    int       size, numBytes=0;
    gEntitySet  ecommon;
    typedef data_schema_traits<T> converter_traits;
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      while(itr1 != attrib_data.end() && itr1->first < ci  )itr1++;
      while(itr1 != attrib_data.end() && itr1->first == ci){
        typename converter_traits::Converter_Type cvtr(const_cast<T&>(itr1->second));
        size      = cvtr.getSize();
        numBytes = size*sizeof(typename converter_traits::Converter_Base_Type) ;
        result += (sizeof(gEntity)+sizeof(int)+ numBytes);
        itr1++; 
      }
    }ENDGFORALL;
    return result;
  }
  
  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gStoreRepI<T>::pack( void *outbuf, int &position, int size, 
                            const gEntitySet &usr_eset ) const 
  {
    fatal(usr_eset - domain() != GEMPTY);
    fatal(!sorted);
    gEntitySet ecommon;
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    ecommon = domain()&usr_eset;
    packdata( traits_type, outbuf, position, size, ecommon);
  }
 
  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gStoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                               int &position, int outcount,
                               const gEntitySet &eset )  const
  {
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      
      while(itr1 != attrib_data.end() && itr1->first < ci) itr1++;
      while(itr1 != attrib_data.end() && itr1->first == ci){
        MPI_Pack( const_cast<gEntity*>(&ci), sizeof(gEntity), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack( const_cast<T*>(&(itr1->second)), sizeof(T), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        itr1++;
      }
    }ENDGFORALL;
  }

  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  void gStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                int &position, int outcount, 
                                const gEntitySet &eset )const
  {
    warn(sorted);
    int  stateSize, incount;
    const_iterator itr1 = attrib_data.begin();
    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  inbuf;
    int typesize = sizeof(dtype);
    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------
    for(gEntitySet :: const_iterator ci = eset.begin(); ci != eset.end(); ci++) {
      while(itr1 != attrib_data.end() && itr1->first < *ci) itr1++;
      while(itr1 != attrib_data.end()&& itr1->first == *ci ){
        typename schema_traits::Converter_Type cvtr( const_cast<T&>(itr1->second));
        stateSize    = cvtr.getSize();
        if( stateSize > static_cast<int>(inbuf.size()) ) inbuf.resize(stateSize);
        cvtr.getState( &inbuf[0], stateSize);
        
        MPI_Pack( const_cast<gEntity*>(&(itr1->first)), sizeof(gEntity), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
                 MPI_COMM_WORLD);
        
        incount =  stateSize*typesize;
        MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
                 MPI_COMM_WORLD) ;
        itr1++;
      }
    }
  }

  //*******************************************************************/

  template <class T> 
  void gStoreRepI<T>::unpack(const void *ptr, int &loc, int size ) 
  {
    
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size);
    local_sort();
  }

  //*********************************************************************/
  template <class T> 
  void gStoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, const void *inbuf,
                                 int &position,  int insize
                                 ) 
  {
    gEntity ei;
    T val;
    while(position != insize){
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &ei,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &val,
                  sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      insert(ei, val);
    }
  }
  //*********************************************************************/
  template <class T> 
  void gStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, const void *inbuf, 
                                  int &position, int insize) 
  {
    int  stateSize, outcount;
   

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  outbuf;
    gEntity ei;
    int typesize = sizeof(dtype);
    while(position != insize){
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &ei,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      MPI_Unpack(const_cast<void*>(inbuf), insize, &position, &stateSize, 1, 
                 MPI_INT, MPI_COMM_WORLD) ;
      if( stateSize > static_cast<int>(outbuf.size()) ) outbuf.resize(stateSize);
      outcount = stateSize*typesize;
      MPI_Unpack(const_cast<void*>(inbuf),insize, &position, &outbuf[0], outcount, 
                 MPI_BYTE, MPI_COMM_WORLD) ;
      T temp;
      typename schema_traits::Converter_Type cvtr( temp );
      cvtr.setState( &outbuf[0], stateSize);
      insert(ei, temp);
    }

  }

 
 
  // //*********************************************************************/
  //   template<class T> 
  //   void gStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, gEntitySet &user_eset)
  //   {
  //     warn(true) ;
  //     /*
  //       typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  //       schema_converter traits_type;

  //       gEntitySet eset, ecommon;

  //       HDF5_ReadDomain(group_id, eset);

  //       ecommon = eset;
  //       allocate( ecommon );

  //       hdf5read( group_id, traits_type, eset,  ecommon );
  //     */
  //   }

  //   //*************************************************************************/

  //   template<class T>
  //   void  gStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDGENTITY_CONVERTER c, frame_info &fi, gEntitySet &eset)
  //   {
  //     warn(true) ;
  //     /*
  //       hsize_t dimension;
  //       size_t indx = 0, arraySize;
  //       int    rank = 1;

  //       gEntitySet::const_iterator ci;
  //       typedef data_schema_traits<T> traits_type;

  //       int num_intervals = user_eset.num_intervals();
  //       interval *it = new interval[num_intervals];

  //       for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

  //       //------------------------------------------------------------------------
  //       // Calculate the offset 
  //       //------------------------------------------------------------------------

  //       store<unsigned> offset;
  //       offset.allocate( eset );

  //       arraySize = 0;
  //       for( ci = eset.begin(); ci != eset.end(); ++ci)
  //       offset[*ci] = arraySize++;

  //       //------------------------------------------------------------------------

  //       DatatypeP dtype = traits_type::get_type();
  //       hid_t vDatatype = dtype->get_hdf5_type();

  //       dimension = arraySize;
  //       hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);   // memory  dataspace
  //       hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
  //       hid_t vDataset   = H5Dopen( group_id, "VariableData");

  //       hssize_t  start[]     = {0};  // determines the starting coordinates.
  //       hsize_t   stride[]    = {1};  // which elements are to be selected.
  //       hsize_t   block[]     = {1};  // size of element block;
  //       hssize_t  foffset[]   = {0};  // location (in file) where data is read.
  //       hsize_t   count[]     = {0};  // how many positions to select from the dataspace

  //       std::vector<T>  data;
  //       for( int k = 0; k < num_intervals; k++) {
  //       count[0] = it[k].second - it[k].first + 1;

  //       if( count[0] > data.size() ) data.resize(count[0] );

  //       foffset[0] = offset[it[k].first];

  //       H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
  //       H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
  //       H5Dread( vDataset, vDatatype, mDataspace, vDataspace, 
  //       H5P_DEFAULT, &data[0]);

  //       indx = 0;
  //       for( int i = it[k].first; i <= it[k].second; i++) 
  //       attrib_data[i] = data[indx++];
  //       }

  //       H5Dclose( vDataset   );
  //       H5Tclose( vDatatype  );
  //       H5Sclose( mDataspace );
  //       H5Sclose( vDataspace );
  //     */
  //   }
  //   //*************************************************************************/

  //   template<class T>
  //   void  gStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, gEntitySet &eset)
  //   {
  //     warn(true) ;
  //     /*
  //       hsize_t  dimension;
  //       size_t   indx = 0, arraySize;
  //       hid_t    vDataset, vDataspace, vDatatype, mDataspace;
  //       gEntitySet::const_iterator  ci;

  //       //---------------------------------------------------------------
  //       // Size of each Bucket ....
  //       //---------------------------------------------------------------
  //       vDatatype  = H5T_NATIVE_INT;
  //       vDataset   = H5Dopen(group_id,"ContainerSize");
  //       vDataspace = H5Dget_space(vDataset);
  //       H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  //       std::vector<int> ibuf(dimension);
  //       H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);
  //       H5Dclose(vDataset  );
  //       H5Sclose(vDataspace);

  //       store<int> container;
  //       container.allocate( eset );

  //       indx      = 0;
  //       arraySize = 0;
  //       for( ci = eset.begin(); ci != eset.end(); ++ci) {
  //       container[*ci] = ibuf[indx++];
  //       arraySize     += container[*ci];
  //       }

  //       //-------------------------------------------------------------------
  //       // Calculate the offset of each gEntity in file ....
  //       //-------------------------------------------------------------------
  //       store<unsigned>   offset;
  //       offset.allocate( eset );

  //       indx     = 0;
  //       for( ci = eset.begin(); ci != eset.end(); ++ci) {
  //       offset[*ci] = indx;
  //       indx       += container[*ci];
  //       }

  //       //------------------------------------------------------------------
  //       // Read the data now ....
  //       //------------------------------------------------------------------
  //       int num_intervals = usr_eset.num_intervals();
  //       if( num_intervals == 0) {
  //       std::cout << "Warning: Number of intervals are zero : " << endl;
  //       return;
  //       }

  //       interval *it = new interval[num_intervals];

  //       for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

  //       typedef typename data_schema_traits<T>::Converter_Base_Type dtype;
  //       typedef data_schema_traits<dtype> traits_type;
  //       DatatypeP atom_type = traits_type::get_type() ;
  //       vDatatype = atom_type->get_hdf5_type();

  //       dimension  = arraySize;
  //       vDataset   = H5Dopen(group_id,"VariableData");
  //       vDataspace = H5Dget_space(vDataset);
  //       mDataspace = H5Dget_space(vDataset);
  //       H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  //       hssize_t  start_mem[] = {0};  // determines the starting coordinates.
  //       hsize_t   stride[]    = {1};  // which elements are to be selected.
  //       hsize_t   block[]     = {1};  // size of element block;
  //       hssize_t  foffset[]   = {0};  // location (in file) where data is read.
  //       hsize_t   count[]     = {0};  // how many positions to select from the dataspace

  //       std::vector<dtype> data;

  //       for( int k = 0; k < num_intervals; k++) {
  //       count[0] = 0;
  //       for( int i = it[k].first; i <= it[k].second; i++)
  //       count[0] +=  container[i];
  //       if( count[0] > data.size() ) data.resize( count[0] );

  //       foffset[0] = offset[it[k].first];
  //       H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride,
  //       count, block);
  //       H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride,
  //       count, block);
  //       H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, &data[0]);

  //       indx = 0;
  //       for( int i = it[k].first; i <= it[k].second; i++) {
  //       typename data_schema_traits<T>::Converter_Type cvtr( attrib_data[i]);
  //       cvtr.setState( &data[0]+indx, container[i] );
  //       indx += container[i];
  //       }
  //       }

  //       H5Tclose(vDatatype );
  //       H5Dclose(vDataset  );
  //       H5Sclose(vDataspace);
  //       H5Sclose(mDataspace);
  //     */
  //   }
  //   //*************************************************************************/

  //   template<class T> 
  //   void gStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, gEntitySet& usr_eset) const
  //   {
  //     warn(true) ;
  //     /*
  //       typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  //       schema_converter traits_output_type;

  //       gEntitySet eset(usr_eset &domain());

  //       HDF5_WriteDomain(group_id, eset);

  //       hdf5write(group_id, traits_output_type, eset);
  //     */
  //   }

  //   //*************************************************************************/
  //   template <class T> 
  //   void gStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDGENTITY_CONVERTER c, const gEntitySet &eset)  const
  //   {
  //     warn(true) ;
  //     /*
  //       int      rank = 1;
  //       hsize_t  dimension;

  //       std::vector<T>   newvec;
  //       gEntitySet :: const_iterator  ei;

  //       int arraySize = eset.size();

  //       if( arraySize < 1) return;

  //       //------------------------------------------------------------------------
  //       // Collect state data from each object and put into 1D array
  //       //------------------------------------------------------------------------
  //       std::vector<T>  data(arraySize);

  //       int indx = 0;
  //       for( ei = eset.begin(); ei != eset.end(); ++ei) {
  //       data[indx++]  = attrib_data.elem(*ei) ;
  //       }

  //       //-------------------------------------------------------------------------
  //       // Write (variable) Data into HDF5 format
  //       //-------------------------------------------------------------------------

  //       typedef data_schema_traits<T> traits_type;
  //       DatatypeP dtype = traits_type::get_type();
  //       hid_t vDatatype = dtype->get_hdf5_type();

  //       dimension        = arraySize;
  //       hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
  //       hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, H5P_DEFAULT);
  //       H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

  //       H5Dclose( vDataset  );
  //       H5Sclose( vDataspace);
  //       H5Tclose( vDatatype );
  //     */
  //   }

  //   //*************************************************************************/

  //   template <class T> 
  //   void gStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const gEntitySet &eset)  const
  //   {   
  //     warn(true) ;
  //     /*
  //       T         Obj;
  //       hid_t     vDataspace, vDataset, vDatatype;
  //       int       rank = 1;
  //       hsize_t   dimension;

  //       gEntitySet :: const_iterator ci;

  //       //-------------------------------------------------------------------------
  //       // Get the sum of each object size and maximum size of object in the 
  //       // container for allocation purpose
  //       //-------------------------------------------------------------------------

  //       typedef data_schema_traits<T> schema_traits ;

  //       size_t  arraySize= 0;
  //       int     stateSize, maxStateSize = 0;

  //       for( ci = eset.begin(); ci != eset.end(); ++ci) {
  //       Obj = attrib_data.elem(*ci) ;
  //       typename schema_traits::Converter_Type cvtr( Obj );
  //       stateSize    = cvtr.getSize();
  //       arraySize   += stateSize;
  //       maxStateSize = max( maxStateSize, stateSize );
  //       } 

  //       //-----------------------------------------------------------------------------
  //       // Collect state data from each object and put into 1D array
  //       //-----------------------------------------------------------------------------
  //       typedef typename schema_traits::Converter_Base_Type dtype;

  //       dtype *data = new dtype[arraySize];

  //       size_t indx = 0;
  //       for( ci = eset.begin(); ci != eset.end(); ++ci) {
  //       Obj = attrib_data.elem(*ci) ;
  //       typename schema_traits::Converter_Type cvtr( Obj );
  //       cvtr.getState( data+indx, stateSize);
  //       indx +=stateSize ;
  //       }

  //       //--------------------------------------------------------------------
  //       // Write (variable) Data into HDF5 format
  //       //--------------------------------------------------------------------
  //       typedef data_schema_traits<dtype> traits_type;
  //       DatatypeP atom_type = traits_type::get_type() ;
  //       vDatatype = atom_type->get_hdf5_type();

  //       dimension  = arraySize;
  //       vDataspace = H5Screate_simple(rank, &dimension, NULL);
  //       vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
  //       vDataspace, H5P_DEFAULT);
  //       H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  //       H5Dclose( vDataset  );
  //       H5Sclose( vDataspace);
  //       H5Tclose( vDatatype );

  //       delete [] data;

  //       //--------------------------------------------------------------------
  //       // Write Container 
  //       //--------------------------------------------------------------------
  //       dimension    = eset.size();
  //       std::vector<int> vbucket(eset.size());

  //       indx = 0;
  //       for( ci = eset.begin(); ci != eset.end(); ++ci) {
  //       Obj = attrib_data.elem(*ci) ;
  //       typename schema_traits::Converter_Type cvtr( Obj );
  //       vbucket[indx++] =  cvtr.getSize();
  //       }

  //       vDatatype  = H5T_NATIVE_INT;
  //       vDataspace = H5Screate_simple(rank, &dimension, NULL);
  //       vDataset   = H5Dcreate(group_id, "ContainerSize",
  //       vDatatype, vDataspace, H5P_DEFAULT);
  //       H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vbucket[0]);

  //       H5Dclose( vDataset  );
  //       H5Sclose( vDataspace);
  //     */
}
//   //*************************************************************************/
#endif

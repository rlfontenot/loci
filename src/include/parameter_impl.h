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
#ifndef PARAMETER_IMPL_H
#define PARAMETER_IMPL_H
#include <string>
namespace Loci {
  //**************************************************************************/

  template<class T> void paramRepI<T>::allocate(const entitySet &p) {
    store_domain = p ;
    dispatch_notify();
  }

  template<class T> void paramRepI<T>::erase(const entitySet& rm) {
    store_domain -= rm ;
    dispatch_notify() ;
  }

  template<class T> void paramRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    dispatch_notify() ;
  }

  template<class T> void paramRepI<T>::
  guarantee_domain(const entitySet& include) {
    store_domain += include ;
    dispatch_notify() ;
  }

  template<class T> gStoreRepP paramRepI<T>:: copy2gstore()const{
    //store domain must be universal key set
    fatal(store_domain != ~EMPTY);
    gParam<T> r ;
    r.set_entitySet(~GEMPTY);
    *r = attrib_data;
    gKeySpaceP space = gKeySpace::get_space("UniverseSpace", "");
    r.set_domain_space(space);
    return r.Rep() ;
  }

  //**************************************************************************/

  template<class T> paramRepI<T>::~paramRepI() {}

  //**************************************************************************/

  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p) const
  {
    return new paramRepI<T>(p) ;
  }

  template<class T>
  storeRep *paramRepI<T>::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a parameter " << endl ;
    return sp ;
  }

  //**************************************************************************/

  template<class T>
  store_type paramRepI<T>::RepType() const
  {
    return PARAMETER ;
  }

  //**************************************************************************/

  template<class T> entitySet paramRepI<T>::domain() const {
    return store_domain ;
  }

  //**************************************************************************/

  template<class T>
  std::ostream &paramRepI<T>::Print(std::ostream &s) const
  {
    entitySet dom = domain() ;
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
  std::istream &paramRepI<T>::Input(std::istream &s)
  {
    entitySet e ;
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
      std::cerr << "Incorrect Format while reading parameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  template<class T>
  frame_info paramRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }

  template<class T>
  frame_info paramRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }
  template<class T>
  frame_info paramRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    typename schema_traits::Converter_Type cvtr(attrib_data);
    stateSize = cvtr.getSize();
    fi.second_level.push_back(stateSize) ;

    return fi ;
  }


  template<class T>
  DatatypeP paramRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T>
  DatatypeP paramRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T>
  DatatypeP paramRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  template<class T>
  void paramRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }

  //**************************************************************************/

  template<class T>
  void paramRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }
  //**************************************************************************/

  template<class T> param<T>::~param() {}

  //**************************************************************************/

  template<class T>
  void param<T>::notification()
  {
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }

  //**************************************************************************/

  template<class T> const_param<T>::~const_param() {}

  //**************************************************************************/

  template<class T>
  void const_param<T>::notification()
  {
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }

  //**************************************************************************/

  template<class T>
  storeRepP paramRepI<T>::remap(const dMap &m) const
  {
    param<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    return r.Rep() ;
  }

  template<class T>
  storeRepP paramRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP paramRepI<T>::thaw() {
    return getRep() ;
  }
  //**************************************************************************/

  template<class T>
  void paramRepI<T>::copy(storeRepP &st, const entitySet &context)
  {
    param<T> p(st) ;
    attrib_data = *p ;
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
    dispatch_notify() ;
  }

  // note this method can only be used when the paramRepI<T> (i.e., *this)
  // is NOT connected to any of the containers since we don't perform any
  // kind of notification
  template<class T>
  void paramRepI<T>::fast_copy(storeRepP &st, const entitySet &context)
  {
    storeRepP true_rep = st->getRep();
    paramRepI<T>* p = dynamic_cast<paramRepI<T>*>(&(*true_rep));
    fatal(p==0);
    attrib_data = p->attrib_data;
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T>
  void paramRepI<T>::gather(const dMap &m, storeRepP &st,
                            const entitySet &context)
  {
    param<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T>
  void paramRepI<T>::scatter(const dMap &m, storeRepP &st,
                             const entitySet &context)
  {

    fatal((context - store_domain) != EMPTY) ;
    store_domain = m.image(context) ;

    param<T> p(st) ;
    attrib_data = *p ;
  }

  //**************************************************************************/

  template <class T>
  int paramRepI<T>::pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }
  //**************************************************************************/
  template <class T>
  int paramRepI<T>::estimated_pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_estimated_mpi_size( schema_converter(), eset );
    
  }
  //**************************************************************************/
  template<class T> int paramRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;    

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size(schema_converter(), packed);
  }
  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {

    return( sizeof(T) ) ;
  }
  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {

    return( sizeof(T) ) ;
  }
  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits;

    typename schema_traits::Converter_Type cvtr(attrib_data);
    int arraySize = cvtr.getSize() ;

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) + sizeof(int));
  }
  //**************************************************************************/

  template <class T>
  int paramRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    return(50*sizeof(double) + sizeof(int));
  } 
  //**************************************************************************/

  template <class T>
  void paramRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &e )
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    packdata( schema_converter(), ptr, loc, size);
  }

  //**************************************************************************/
  template <class T>
  void paramRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                               int outcount )
  {
    MPI_Pack( &attrib_data, sizeof(T), MPI_BYTE, outbuf, outcount, &position,
              MPI_COMM_WORLD) ;
  }

  //**************************************************************************/

  template <class T>
  void paramRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                               int &position, int outcount )
  {
    int stateSize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    typename data_schema_traits<T>::Converter_Type cvtr( attrib_data );

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
  void paramRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata( schema_converter(), ptr, loc, size);
  }


  //**************************************************************************/
  template <class T>
  void paramRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position,
                                 int &insize)
  {

    /*
      typedef data_schema_traits<T> traits_type;
      DatatypeP    atom_type = traits_type::get_type();
      MPI_Datatype datatype  = atom_type->get_mpi_type();
    */
    MPI_Unpack( inbuf, insize, &position, &attrib_data, sizeof(T),
                MPI_BYTE, MPI_COMM_WORLD) ;

  }

  //***********************************************************************/
  template <class T>
  void paramRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf,
                                 int &position, int &insize)
  {

    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    MPI_Unpack( inbuf, insize, &position, &stateSize, 1,
                MPI_INT, MPI_COMM_WORLD) ;
    std::vector<dtype> outbuf(stateSize);

    outcount = stateSize*sizeof(dtype);
    MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount,
                MPI_BYTE, MPI_COMM_WORLD) ;
    typename schema_traits::Converter_Type  cvtr( attrib_data );
    cvtr.setState( &outbuf[0], stateSize);

  }

  //***********************************************************************/

  template<class T> store_instance::instance_type
  const_param<T>::access() const
  { return READ_ONLY; }

  //**************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const const_param<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/
  template <class T>
  void paramRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;

      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &attrib_data) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
    }
  }

  //*********************************************************************/

  template <class T>
  void paramRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
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
      int stateSize = 0 ;
      T tmp = attrib_data;
      typename schema_traits::Converter_Type cvtr(tmp);
      cvtr.getState(tmp_array, stateSize) ;
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //**************************************************************************/
  template <class T>
  void paramRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, frame_info &fi, const entitySet &en)
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
      }
      attrib_data = tmp_array[0] ;

      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //*************************************************************************/

  template <class T>
  void paramRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const entitySet &en)
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
        cerr << "H5Dread() failed" << endl ;
      }
      typename data_schema_traits<T>::Converter_Type cvtr(attrib_data);
      int bucsize = vint[0] ;
      cvtr.setState(tmp_array, bucsize) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  template<class T> storeRepP paramRepI<T>::
  redistribute(const std::vector<entitySet>& dom_ptn, MPI_Comm comm) {
    // for a parameter, we just redistribute its domain
    entitySet dom = domain() ;
    entitySet new_all ;
    for(size_t i=0;i<dom_ptn.size();++i)
      new_all += dom_ptn[i] ;
    entitySet out = dom - new_all ;

    std::vector<entitySet> old_dist = all_collect_vectors(dom, comm) ;
    entitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    entitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom += out ;

    param<T> np ;
    np.set_entitySet(new_dom) ;
    *np = attrib_data ;

    return np.Rep() ;
  }

  template<class T> storeRepP paramRepI<T>::
  redistribute(const std::vector<entitySet>& dom_ptn,
               const dMap& remap, MPI_Comm comm) {
    // for a parameter, we just redistribute its domain
    entitySet dom = domain() ;
    std::vector<entitySet> old_dist = all_collect_vectors(dom, comm) ;
    entitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    entitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom = remap_entitySet(new_dom, remap) ;

    param<T> np ;
    np.set_entitySet(new_dom) ;
    *np = attrib_data ;

    return np.Rep() ;
  }
  template<class T> storeRepP paramRepI<T>::
  redistribute_omd(const std::vector<entitySet>& dom_ptn,
                   const dMap& remap, MPI_Comm comm) {
    return redistribute(dom_ptn,remap,comm) ;
  }
  //***************************************************************************

}
#endif

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
#ifndef MAPVEC_IMPL_H
#define MAPVEC_IMPL_H

#include <MapVec_def.h>
#include <DMapVec_def.h>

namespace Loci {
  
  template<int M> void MapVecRepI<M>::allocate(const entitySet &ptn) {

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    base_ptr = 0 ;
    if(ptn != EMPTY) {
      int top = ptn.Min() ;
      int size = ptn.Max()-top+1 ;
      alloc_pointer = new VEC[size] ;
      base_ptr = alloc_pointer-top ;
    }
    store_domain = ptn ;
    dispatch_notify() ;
  }

  //*************************************************************************/

  template<int M> MapVecRepI<M>::~MapVecRepI() {
    if(alloc_pointer) delete[] alloc_pointer ;
  }

  //*************************************************************************/

  template<int M> storeRepP MapVecRepI<M>::get_map()  {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = M ;
    } ENDFORALL ;
    multiMap result ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      for(int j=0;j<M;++j) 
        result.begin(i)[j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return result.Rep() ;
  }

  //*************************************************************************/

  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p) const {
    return new MapVecRepI<M>(p) ;
  }
  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p, const int* cnt) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a MapVec " << endl ;
    return sp ;
  }
  //*************************************************************************/

  template<int M> entitySet MapVecRepI<M>::domain() const {
    return store_domain ;
  }

  //*************************************************************************/

  template<int M> entitySet MapVecRepI<M>::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(int i=0;i<d.num_intervals();++i)
        codomain += Loci::image_section(&base_ptr[d[i].first][0],
                                        &base_ptr[d[i].second+1][0]) ;
    } else {
      std::vector<int> img(d.size()*M) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(int i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          for(int k=0;k<M;++k) {
            *ins = base_ptr[j][k] ;
            ++ins ;
          }
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
    return codomain ;
  }

  //*************************************************************************/

  template<int M> std::pair<entitySet,entitySet>
  MapVecRepI<M>::preimage(const entitySet &codomain) const {
    entitySet domaini ;
    entitySet domainu ;
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false ;
      for(int j=0;j<M;++j) {
        bool in_set = codomain.inSet(base_ptr[i][j]) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return std::make_pair(domaini,domainu) ;
  }

  //*************************************************************************/
  
  template<int M> std::ostream &MapVecRepI<M>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii][0] ;
      for(int j=1;j<M;++j)
        s << " " << base_ptr[ii][j] ;
      s << std::endl ;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }

  //*************************************************************************/

  template<int M> std::istream &MapVecRepI<M>::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      for(int j=0;j<M;++j)
        s >> base_ptr[ii][j] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //*************************************************************************/
  template<int M> void MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  template<int M> inline std::ostream & operator<<(std::ostream &s, const MapVec<M> &m){
    return m.Print(s) ;
  }
    
  template<int M> inline std::istream & operator>>(std::istream &s, MapVec<M> &m) {
    return m.Input(s) ;
  }
  //*************************************************************************/
  template<int M> void const_MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  //*************************************************************************/

  template<int M> store_instance::instance_type
  const_MapVec<M>::access() const { return READ_ONLY; }
    
  //*************************************************************************/

  template<int M> inline std::ostream & operator<<(std::ostream &s,
                                                   const const_MapVec<M> &m){
    return m.Print(s) ;
  }

  //*************************************************************************/

  template<int M> storeRepP MapVecRepI<M>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    std::pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;
    MapVec<M> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    
    return s.Rep() ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::compose(const dMap &m,
                                              const entitySet &context)  {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    entitySet dom = m.domain() ;
    FORALL(context,i) {
      for(int j=0;j<M;++j) {
	if(dom.inSet(base_ptr[i][j]))
          base_ptr[i][j] = m[base_ptr[i][j]] ;
	else
	  base_ptr[i][j] = -1 ;
      }
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::copy(storeRepP &st,
                                           const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::gather(const dMap &m, storeRepP &st,
                                             const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal(base_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::scatter(const dMap &m, storeRepP &st,
                                              const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((base_ptr == 0) &&(context != EMPTY)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[m[i]][j] = s[i][j] ;
    } ENDFORALL ;
  }

  //*************************************************************************/
  template <int M> int MapVecRepI<M>::pack_size( const entitySet &eset) {
    fatal((eset - domain()) != EMPTY);
    int size ;
    size = sizeof(int)*eset.size()*M  + sizeof(int);
    return size ;
  }

  //*************************************************************************/
  template <int M> void MapVecRepI<M>::pack(void *outbuf, int &position, 
                                            int &outcount, const entitySet &eset ) 
  {
    int init_size = M;
    MPI_Pack( &init_size, 1,  MPI_INT, outbuf, outcount, &position, 
              MPI_COMM_WORLD) ;

    entitySet::const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) 
      MPI_Pack( base_ptr[*ci], M, MPI_INT, outbuf, outcount, &position, 
                MPI_COMM_WORLD);
  }
  //*************************************************************************/

  template <int M> void MapVecRepI<M>::unpack(void *inbuf, int &position, int &insize, 
                                              const sequence &seq) {

    int init_size;

    MPI_Unpack(inbuf, insize, &position, &init_size, 1, MPI_INT, MPI_COMM_WORLD) ;

    if(init_size != M) {
      std::cerr << " Invalid MapVec container for unpack data " << std::endl;
      Loci::Abort();
    }

    sequence::const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) 
      MPI_Unpack( inbuf, insize, &position, base_ptr[*ci], M, MPI_INT, 
                  MPI_COMM_WORLD) ;
  }
  //*************************************************************************/
  
  template<int M> storeRepP
    MapVecRepI<M>::expand(entitySet &out_of_dom,
                          std::vector<entitySet> &init_ptn) {
    // first do a check to make sure that out_of_dom is
    // entirely within the init_ptn partitions
    entitySet total ;
    for(std::vector<entitySet>::const_iterator vi=init_ptn.begin();
        vi!=init_ptn.end();++vi)
      total += *vi ;
    if(total - out_of_dom != EMPTY) {
      cerr << "Error: MapVecRepI::expand failed! "
           << "out_of_dom passed in is not containted "
           << "entirely within init_ptn." << endl ;
      Abort() ;
    }
    // then divide "out_of_dom" to every other process
    std::vector<std::vector<int> > copy_from(MPI_processes) ;
    int* send_count = new int[MPI_processes] ;
    int* recv_count = new int[MPI_processes] ;
    int total_send_size = 0 ;
    
    for(int i=0;i<MPI_processes;++i) {
      entitySet dom = out_of_dom & init_ptn[i] ;
      for(entitySet::const_iterator ei=dom.begin();
          ei!=dom.end();++ei)
        copy_from[i].push_back(*ei) ;

      send_count[i] = copy_from[i].size() ;
      total_send_size += send_count[i] ;
    }
    // notify everyone else of this information so that
    // every one knows what to send to others
    MPI_Alltoall(send_count, 1, MPI_INT,
                 recv_count, 1, MPI_INT, MPI_COMM_WORLD) ;

    // preparing buffer to send the domains to everyone else
    int total_recv_size = 0 ;
    for(int i=0;i<MPI_processes;++i)
      total_recv_size += recv_count[i] ;

    int* send_buf = new int[total_send_size] ;
    int* recv_buf = new int[total_recv_size] ;

    // fill in the send buffer
    int idx_count = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      for(std::vector<int>::const_iterator vi=copy_from[i].begin();
          vi!=copy_from[i].end();++vi) {
        send_buf[idx_count++] = *vi ;
      }
    }

    // prepare the send/recv displacement buffer
    int* send_displs = new int[MPI_processes] ;
    int* recv_displs = new int[MPI_processes] ;

    send_displs[0] = 0 ;
    recv_displs[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displs[i] = send_displs[i-1] + send_count[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_count[i-1] ;
    }

    // ready to communicate the entity domains buffer
    MPI_Alltoallv(send_buf, send_count, send_displs, MPI_INT,
                  recv_buf, recv_count, recv_displs, MPI_INT, MPI_COMM_WORLD) ;

    // extract the entities from the raw buffer
    std::vector<std::vector<int> > send_to_dom(MPI_processes) ;

    for(int i=0;i<MPI_processes;++i) {
      for(int j=recv_displs[i];j<recv_displs[i]+recv_count[i];++j) {
        send_to_dom[i].push_back(recv_buf[j]) ;
      }
    }

    delete[] send_buf ;
    delete[] recv_buf ;
    
    // then pack in the contents in the map for the domains received
    std::vector<std::map<int,std::vector<int> > > send_to(MPI_processes) ;

    for(int i=0;i<MPI_processes;++i) {
      for(std::vector<int>::const_iterator vi=send_to_dom[i].begin();
          vi!=send_to_dom[i].end();++vi) {
        for(unsigned int k=0;k<M;++k) {
          send_to[i][*vi].push_back(base_ptr[*vi][k]) ;
        }
      }
    }

    // recompute send_count/recv_count for sending the actual contents
    total_send_size = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      // the format for the actual send buffer is as:
      // domain entity id: actual contents
      // NOTE: we don't need to send the size of the map contents
      // for an entity because this is a MapVec<M> and we know
      // that the size is M
      for(std::map<int,std::vector<int> >::const_iterator
            mi=send_to[i].begin();mi!=send_to[i].end();++mi) {
        send_count[i] = 1 + mi->second.size() ;
      }
      total_send_size += send_count[i] ;
    }
    
    // communicate the size first
    MPI_Alltoall(send_count, 1, MPI_INT,
                 recv_count, 1, MPI_INT, MPI_COMM_WORLD) ;

    // calculate the recv size
    total_recv_size = 0 ;
    for(int i=0;i<MPI_processes;++i)
      total_recv_size += recv_count[i] ;

    // preparing the buffers to send/recv the actual contents
    send_buf = new int[total_send_size] ;
    recv_buf = new int[total_recv_size] ;

    // fill the send buffer
    idx_count = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      for(std::map<int,std::vector<int> >::const_iterator
            mi=send_to[i].begin();mi!=send_to[i].end();++mi) {
        send_buf[idx_count++] = mi->first ;
        for(unsigned int k=0;k<M;++k) {
          send_buf[idx_count++] = (mi->second)[k] ;
        }
      }
    }

    // calculate the send/recv displacement buffer
    send_displs[0] = 0 ;
    recv_displs[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displs[i] = send_displs[i-1] + send_count[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_count[i-1] ;
    }

    // Okay, we are finally ready to communicate the actual map contents
    MPI_Alltoallv(send_buf, send_count, send_displs, MPI_INT,
                  recv_buf, recv_count, recv_displs, MPI_INT, MPI_COMM_WORLD) ;

    // now we need to extract the contents from the raw buffer
    entitySet new_domain = store_domain + out_of_dom ;
    MapVec<M> new_map ;
    new_map.allocate(new_domain) ;
    // fill in the original contents first
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei) {
      for(unsigned int i=0;i<M;++i)
        new_map[*ei][i] = base_ptr[*ei][i] ;
    }
    // fill the expanded domain
    for(int i=0;i<MPI_processes;++i) {
      int begin = recv_displs[i] ;
      int end = begin + recv_count[i] ;
      int current = begin ;
      while(current != end) {
        int e = recv_buf[current++] ;
        for(unsigned int k=0;k<M;++k)
          new_map[e][k] = recv_buf[current++] ;
      }
    }

    // release all the allocated buffer
    delete[] send_count ;
    delete[] recv_count ;
    delete[] send_displs ;
    delete[] recv_displs ;
    delete[] send_buf ;
    delete[] recv_buf ;

    cerr << "MapVecRepI::expand called!" << endl ;
    return new_map.Rep() ;
    
  } // end of expand

  template<int M> storeRepP MapVecRepI<M>::freeze() {
    return getRep() ;
  }
  
  template<int M > storeRepP MapVecRepI<M>::thaw() {
    dMapVec<M> dm ;
    FORALL(store_domain, i) {
      for(int j=0;j<M;++j)
        dm[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return dm.Rep() ;
  }
  //*************************************************************************/
  
  template<int M> void inverseMap (multiMap &result,
                                   const MapVec<M> &input_map,
                                   const entitySet &input_image,
                                   const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k)
        if(input_image.inSet(input_map[i][k]))
          sizes[input_map[i][k]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) {
          sizes[elem] -= 1 ;
          FATAL(sizes[elem] < 0) ;
          result[elem][sizes[elem]] = i ;
        }
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }      
  template<int M> void inverseMap (dmultiMap &result,
                                   const MapVec<M> &input_map,
                                   const entitySet &input_image,
                                   const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k)
        if(input_image.inSet(input_map[i][k]))
          sizes[input_map[i][k]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) {
          sizes[elem] -= 1 ;
          FATAL(sizes[elem] < 0) ;
          result[elem][sizes[elem]] = i ;
        }
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }      

  
  template<int M> 
    DatatypeP MapVecRepI<M>::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  template<int M> 
    frame_info MapVecRepI<M>::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  
  //*************************************************************************/
  
  template<int M> 
  void MapVecRepI<M>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    warn(true) ;
    /*
    hsize_t dimension[1];
    int  indx=0, rank=1, size;
 
    Loci::HDF5_ReadDomain(group_id, eset);
    Loci::HDF5_ReadVecSize(group_id, &size);

    allocate( eset );

    entitySet::const_iterator ci;

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<unsigned> offset;
    offset.allocate(eset);

    int arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += size;
    }

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = eset[i];

    int   *data;
    dimension[0] = size*eset.size();

    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDataset   = H5Dopen( group_id, "MapVec");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  size;

      data = new int[count[0]];

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);
 
      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int ivec = 0; ivec < size; ivec++)
          base_ptr[i][ivec] = data[indx++];
      }
      delete[] data;
    }

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Sclose( mDataspace);
    */
  }

  //*************************************************************************/
    
  template<int M> 
  void MapVecRepI<M>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const 
  {
    warn(true) ;
    /*
    int     rank = 1;
    hsize_t dimension;

    entitySet eset(usr_eset &domain());
    int arraySize = M*eset.size();
    if( arraySize < 1) return;

    Loci::HDF5_WriteDomain( group_id, eset);
    Loci::HDF5_WriteVecSize(group_id, M);

    std::vector<int> data(arraySize);
    entitySet::const_iterator ci;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < M; ivec++)
        data[indx++] =  base_ptr[*ci][ivec];
    }

    dimension        = arraySize;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dcreate(group_id, "MapVec", vDatatype, vDataspace, 
                                 H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    */
  }
  //*************************************************************************/

} // end of namespace Loci


#endif

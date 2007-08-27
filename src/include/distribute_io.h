#ifndef DISTRIBUTE_IO_H
#define DISTRIBUTE_IO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <store_rep.h>
#include <DMap.h>
#include <fact_db.h>
#include <mpi.h>
#include <store.h>

namespace Loci {

  //  void write_container(hid_t group_id, storeRepP qrep) ;
  //  void read_container(hid_t group_id, storeRepP qrep, entitySet &dom) ;

  //  storeRepP collect_reorder_store(storeRepP &sp, dMap &remap, fact_db &facts) ;
  //  void distribute_reorder_store(storeRepP &new_sp, storeRepP sp_init, dMap &remap, fact_db &facts) ;

  void redistribute_write_container(hid_t file_id, std::string vname,
                                    Loci::storeRepP var, fact_db &facts) ;
  void read_container_redistribute(hid_t file_id, std::string vname,
                                   Loci::storeRepP var, entitySet read_set,
                                   fact_db &facts) ;

  inline hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id) {
    if(Loci::MPI_rank==0)
      return H5Fcreate(name,flags,create_id,access_id) ;
    else
      return 0 ;
  }

  inline hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id) {
    if(Loci::MPI_rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
  }

  inline herr_t hdf5CloseFile(hid_t file_id) {
    if(Loci::MPI_rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }
    
  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var, fact_db &facts) {

    redistribute_write_container(file_id,vname,var,facts) ;
  }
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet, fact_db &facts) {
    read_container_redistribute(file_id,vname,var,readSet, facts) ;
  }

  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::writeContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      redistribute_write_container(file_id,vname,var,
                                   *Loci::exec_current_fact_db) ;
  }
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::readContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      read_container_redistribute(file_id,vname,var,readSet,
                                  *Loci::exec_current_fact_db) ;
  }

  template<class T> void writeUnorderedVector(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v) {
    int local_size = v.size() ;
    std::vector<int> recv_sizes(MPI_processes) ;
    MPI_Gather(&local_size,1,MPI_INT,
               &recv_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;

    if(MPI_rank == 0) {
      int array_size = 0 ;
      for(int i=0;i<MPI_processes;++i)
        array_size += recv_sizes[i] ;
      if(array_size == 0)
        return ;
      int rank = 1 ;
      hsize_t dimension = array_size ;

      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = recv_sizes[0] ;
      hid_t dataset = H5Dcreate(group_id,element_name,dp->get_hdf5_type(),
                                dataspace, H5P_DEFAULT) ;
      if(count != 0) {
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &v[0]) ;
        H5Sclose(memspace) ;
      }
      for(int i=1;i<MPI_processes;++i) {
        start += recv_sizes[i-1] ;
        if(recv_sizes[i] == 0)
          continue ;
        int flag = 0 ;
        MPI_Send(&flag,1,MPI_INT,i,0,MPI_COMM_WORLD) ;
        std::vector<T> rv(recv_sizes[i]) ;
        MPI_Status mstat ;
        MPI_Recv(&rv[0],sizeof(T)*recv_sizes[i],MPI_BYTE,i,0,MPI_COMM_WORLD,
                 &mstat) ;
        count = recv_sizes[i] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 H5P_DEFAULT, &rv[0]) ;
        H5Sclose(memspace) ;
      }

      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      
    } else {
      if(local_size == 0)
        return ;

      int flag = 0;
      MPI_Status mstat ;
      MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&mstat) ;
      MPI_Send(&v[0],sizeof(T)*local_size,MPI_BYTE,0,0,MPI_COMM_WORLD) ;
    }
  }

  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) ;
  
  hid_t createUnorderedFile(const char * filename, entitySet set,
                            fact_db &facts) ;

  inline hid_t createUnorderedFile(const char * filename, entitySet set) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::createUnorderedFile()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
      return -1 ;
    } else
      return createUnorderedFile(filename, set, *Loci::exec_current_fact_db) ;
  }

  void closeUnorderedFile(hid_t file_id) ;

  template<class T> void writeUnorderedStore(hid_t file_id,
                                             const_store<T> &s, entitySet set,
                                             const char *name) {
    std::vector<T> v(set.size()) ;
    int c = 0 ;
    FORALL(set,ii) {
      v[c++] = s[ii] ;
    } ENDFORALL ;
    writeUnorderedVector(file_id,name,v) ;
  }

  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells,
                                 fact_db &facts) ;

  inline
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::parallelWriteGridTopology()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      parallelWriteGridTopology(filename, upperRep, lowerRep, boundary_mapRep,
                                face2nodeRep, refRep, bnamesRep, posRep,
                                localCells, *Loci::exec_current_fact_db) ;
  }
  


    // Utility routine for sample sort
  template <class T,class Cmp> void parGetSplitters(std::vector<T> &splitters,
                                                    const std::vector<T> &input,
                                                    Cmp cmp,
                                                    MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    splitters = std::vector<T>(p-1) ;
    std::vector<T> allsplits(p*(p-1)) ;

    int nlocal = input.size() ;
    if(nlocal < p) {
      std::cerr << "sample sort needs at least p elements per processor"
                << std::endl ;
    }
    for(int i=1;i<p;++i) 
      splitters[i-1] = input[(i*nlocal)/p] ;

    int tsz = sizeof(T) ;
    MPI_Allgather(&splitters[0],(p-1)*tsz,MPI_BYTE,
                  &allsplits[0],(p-1)*tsz,MPI_BYTE,comm) ;
    
    sort(allsplits.begin(),allsplits.end(),cmp) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }


  template <class T, class Cmp>
  void parSplitSort(std::vector<T> &list, std::vector<T> &splitters,
                    Cmp cmp, MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;

    int s=0 ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || cmp(list[i] , splitters[s]) ) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !cmp(list[i],splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(T) ;

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    std::vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(T) ;

    std::vector<T> sorted_pnts(result_size) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    list.swap(sorted_pnts) ;
    sort(list.begin(),list.end()) ;
    return ;
  }

  template <class T, class Cmp>
  void parSampleSort(std::vector<T> &list, Cmp cmp, MPI_Comm comm) {
    // First sort list locally
    sort(list.begin(),list.end(),cmp) ;

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;

    std::vector<T> splitters ;
    parGetSplitters(splitters,list,cmp,comm) ;

    parSplitSort(list,splitters,cmp,comm) ;
  }

}

#endif

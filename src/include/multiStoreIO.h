#ifndef MULTISTOREIO_H
#define MULTISTOREIO_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <dist_internal.h>
#include <distribute_container.h>
#include <distribute.h>
#include <distribute_io.h>
#include <multiStore.h>
#include "execute.h"

#include <mpi.h>
#include <hdf5.h>


#define DXFER_COLLECTIVE_IO 0x1  /* Collective IO*/
#define DXFER_INDEPENDENT_IO 0x2 /* Independent IO collectively */
//#define dxfer_coll_type 0x1 /* define the data transfer type as Collective IO*/

//from partest/testphdf5.h,version 1.10.2
#define FACC_DEFAULT    0x0     /* default */
#define FACC_MPIO       0x1     /* MPIO */
namespace Loci {
  namespace hdf5_const {
    extern const int PPN; //processes per node
    extern const int facc_type; //file access type
    extern const int dxfer_coll_type;/* define the data transfer type as Collective IO or Independent IO collectively*/
  }


  /*
   * Create the appropriate File access property list
   */
  hid_t create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type);
  
  hid_t  create_xfer_plist(int l_xfer_type);


  hid_t hdf5PCreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, size_t file_size_estimate,MPI_Comm comm) ;

  inline hid_t hdf5PCreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id,size_t file_size_estimate) {
    return hdf5PCreateFile(name,flags,create_id,access_id,file_size_estimate,
			   MPI_COMM_WORLD) ;
  }

  inline hid_t hdf5POpenFile(const char *name, unsigned flags, hid_t access_id) {
#ifndef H5_HAVE_PARALLEL
    if(Loci::MPI_rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
#else
      return H5Fopen(name,flags,access_id) ;
#endif    
  }

  inline hid_t hdf5POpenFile(const char *name, unsigned flags, hid_t access_id,
                            MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fopen(name,flags,access_id) ;
    else
      return 0 ;
#else
      return H5Fopen(name,flags,access_id) ;
#endif
  }

  inline herr_t hdf5PCloseFile(hid_t file_id) {
#ifndef H5_HAVE_PARALLEL
    if(Loci::MPI_rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
#else
      return H5Fclose(file_id) ;
#endif
  }

  inline herr_t hdf5PCloseFile(hid_t file_id, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
#else
    return H5Fclose(file_id) ;
#endif
  }

  
  void findMapping(std::vector<int> &local_num,
                   std::vector<int> &procID,
                   entitySet read_set,
                   const std::vector<int> &fileID,
                   fact_db &facts) ;

  void distributeMapMultiStore(std::vector<int> &send_sz,
                               std::vector<int> &recv_sz,
                               std::vector<int> &recv_local_num,
                               const std::vector<int> &local_num,
                               const std::vector<int> &procID) ;

  void sendCounts(std::vector<int>&recv_count,
                  const std::vector<int> &send_sz,
                  const std::vector<int> &recv_sz,
                  const std::vector<int> &recv_local_num,
                  const std::vector<int> &counts,
                  const std::vector<int> &procID) ;

  template<class T> void writeVectorSerial(hid_t group_id,
                                           const char *element_name,
                                           std::vector<T> &v ) {
    hsize_t array_size_combined = v.size() ;
    if(array_size_combined == 0)
      return ;

    int rank = 1 ;
    hsize_t dimension = array_size_combined ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

    typedef data_schema_traits<T> traits_type ;
    DatatypeP dp = traits_type::get_type() ;

    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    hsize_t count = v.size() ;
    hid_t datatype = dp->get_hdf5_type() ;

    hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                              dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

    if(count != 0) {

      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                          &start, &stride, &count, NULL) ;

      hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

      H5Dwrite(dataset,datatype,memspace,dataspace,
               H5P_DEFAULT, &v[0]) ;

      H5Sclose(memspace) ;

    }
    H5Dclose(dataset) ;
    H5Sclose(dataspace) ;
    H5Tclose(datatype) ;
    return;
  }

  inline size_t containerSizeEstimateKb(storeRepP p) {
    entitySet dom = p->domain() ;
    int szkb = 0 ;
    for(size_t i=0;i<dom.num_intervals();++i) {
      int is = dom[i].first ;
      int ie = dom[i].second ;
      if(ie-is>100) {
	entitySet tmp = interval(is,ie) ;
	szkb += (p->pack_size(tmp)>>10) ;
      } else {
	int dlta = (ie-is+1)/10 ;
	while(is <= ie) {
	  entitySet tmp= interval(is,min(is+dlta-1,ie)) ;
	  szkb += (p->pack_size(tmp))>>10 ;
	  is += dlta ;
	}
      }
    }
    int szkbtot = szkb ;
    MPI_Allreduce(&szkb,&szkbtot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    return szkbtot ;
  }
    
  template< class T > inline void writeMultiStore(hid_t file_id,
                                                  std::string vname,
                                                  multiStore<T> &var,
                                                  entitySet write_set,
                                                  fact_db &facts) {
    const_multiStore<T> var_const ;
    var_const = var.Rep() ;
    writeMultiStore(file_id,vname,var_const,write_set,facts) ;
  }

  template< class T > inline void writeMultiStore(hid_t file_id,
                                                  std::string vname,
                                                  const const_multiStore<T> &var,
                                                  entitySet write_set,
                                                  fact_db &facts) {
    hid_t group_id = -1 ;
    if(MPI_rank == 0) {
      group_id = H5Gcreate(file_id,vname.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    }
    std::vector<int> sizes_local ;

    double time_write = 0 ;
    double pre_time = 0 ;
    double time_wait = 0 ;
    entitySet dom = write_set ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    if(df != 0)
      write_set = write_set & df->my_entities ;

    Loci::stopWatch sp ;
    Loci::stopWatch sw ;
    sp.start() ;
    { // Write out book keeping stuff.
      std::vector<int> counts(dom.size()) ;
      int cnt = 0 ;
      std::vector<int> fileids(dom.size()) ;
      if(df == 0) {
        FORALL(dom,ii) {
          counts[cnt] = var.vec_size(ii) ;
          int file_no = ii ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;
      } else {
        Map l2g ;
        l2g = df->l2g.Rep() ;
        dMap g2f ;
        g2f = df->g2f.Rep() ;

        FORALL(dom,ii) {
          counts[cnt] = var.vec_size(ii) ;
          int file_no = g2f[l2g[ii]] ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;
      }
      writeUnorderedVector(group_id,"counts",counts,MPI_COMM_WORLD) ;
      writeUnorderedVector(group_id,"fileID",fileids,MPI_COMM_WORLD) ;
      std::sort(counts.begin(),counts.end()) ;
      std::vector<int>::const_iterator lu = std::unique(counts.begin(),counts.end()) ;
      int nunique = lu-counts.begin() ;
      std::vector<int> tmp(nunique) ;
      sizes_local.swap(tmp) ;
      for(int i=0;i<nunique;++i)
        sizes_local[i] = counts[i] ;
    }

    // schedule blocking for writing
    // Find global sizes array
    sizes_local.push_back(std::numeric_limits<int>::max()) ;
    int iloc = 0 ;
    std::vector<int> sizes ;
    int cmin = -1 ;
    do {
      MPI_Allreduce(&sizes_local[iloc],&cmin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      if(cmin != std::numeric_limits<int>::max())
        sizes.push_back(cmin) ;
      if(sizes_local[iloc] == cmin)
        iloc++ ;
    } while (cmin != std::numeric_limits<int>::max()) ;

    std::vector<int> block_sizes ;
    size_t cnt = 0 ;
    if(sizes[0] == 0)
      cnt++ ;
    const int blk_factor = 256 ;
    int cloc = 0 ;
    int cloc_prev = 0 ;
    while(cnt < sizes.size()) {
      cloc += blk_factor ;
      cloc = min(cloc,sizes[cnt]) ;
      block_sizes.push_back(cloc-cloc_prev) ;
      cloc_prev = cloc ;
      if(cloc == sizes[cnt])
        cnt++ ;
    }
    if(block_sizes.size() == 0) {
      block_sizes.push_back(0) ;
    }
    std::vector<int> num_local_xfers(block_sizes.size(),0) ;

    FORALL(dom,ii) {
      int lsz = var.vec_size(ii) ;
      int tot = 0 ;
      for(size_t i=0;i<block_sizes.size();++i) {
        tot += block_sizes[i] ;
        if(lsz >= tot && block_sizes[i] != 0)
          num_local_xfers[i]++ ;
      }
    } ENDFORALL ;

    std::vector<int> block_data_elems(block_sizes.size(),0) ;
    MPI_Allreduce(&num_local_xfers[0],&block_data_elems[0],block_sizes.size(),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;

    int p = MPI_processes ;
    int r = MPI_rank ;

    // write out write schedule information
    if(r == 0) {
      writeVectorSerial(group_id,"block_schedule",block_sizes) ;
      writeVectorSerial(group_id,"block_sets", block_data_elems) ;
    }
    pre_time += sp.stop() ;
    // Now write out the main data block

    size_t total_size = 0 ;
    for(size_t i=0;i<block_sizes.size();++i) {
      total_size += block_sizes[i]*block_data_elems[i] ;
    }
    int max_local_size = 0 ;
    for(size_t i=0;i<block_sizes.size();++i) {
      int local_size = block_sizes[i]*num_local_xfers[i] ;
      max_local_size = max(max_local_size,local_size) ;
    }
    int max_size = max_local_size ;
    MPI_Allreduce(&max_local_size,&max_size,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(r == 0)
      max_local_size = max_size ;

    std::vector<int> all_local_sizes(block_sizes.size()*p) ;

    MPI_Gather(&num_local_xfers[0],block_sizes.size(),MPI_INT,
               &all_local_sizes[0],block_sizes.size(),MPI_INT,
               0, MPI_COMM_WORLD) ;

    const int nadvance = 3 ;
    std::vector<std::vector<T> > buffer(nadvance) ;
    if(r == 0) {
      for(int i=0;i<min(nadvance,p);++i) {
        std::vector<T> tmp(max_local_size) ;
        buffer[i].swap(tmp) ;
      }
    } else {
      std::vector<T> tmp(max_local_size) ;
      buffer[0].swap(tmp) ;
    }
    //    std::vector<T> buffer(max_local_size) ;

    int rank = 1 ;
    hsize_t dimension = total_size ;
    hid_t dataspace =-1;
    if(r == 0)
      dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    typedef data_schema_traits<T> traits_type ;
    DatatypeP dp = traits_type::get_type() ;
    hid_t datatype =-1;
    if(r == 0)
      datatype = dp->get_hdf5_type() ;
    hid_t dataset =-1;
    if(r == 0) {
      dataset = H5Dcreate(group_id,"data",datatype,
                          dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
    }
    int loc = 0 ;
    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    for(size_t i=0;i<block_sizes.size();++i) {
      int bs = block_sizes[i] ;
      std::vector<MPI_Request> recv_Requests(nadvance) ;
      std::vector<MPI_Request> send_Requests(nadvance) ;
      if(r == 0)  {
        for(int k=1;k<min(nadvance,p);++k) {
          hsize_t count = all_local_sizes[k*block_sizes.size()+i]*bs ;
          if(count != 0) {
            int id = k ;
            static int flag = 0 ;

            MPI_Irecv(&buffer[id][0],sizeof(T)*count,MPI_BYTE,k,0,MPI_COMM_WORLD,
                      &recv_Requests[id]) ;
            //      MPI_Isend(&flag,1,MPI_INT,k,0,MPI_COMM_WORLD,&send_Requests[id]) ;
            MPI_Send(&flag,1,MPI_INT,k,0,MPI_COMM_WORLD) ;

          }
        }
        int bcnt = 0 ;
        FORALL(dom,ii) {
          int lsz = var.vec_size(ii) ;
          if(lsz >= loc+bs) {
            const T *base = var.begin(ii) ;
            for(int j=0;j<bs;++j) {
              buffer[0][bcnt++] = base[loc+j] ;
            }
          }
        } ENDFORALL ;
        hsize_t count = num_local_xfers[i]*bs ;
        if(count != 0) {
          Loci::stopWatch s ;
          s.start() ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start,&stride,&count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&buffer[0][0]) ;
          H5Sclose(memspace) ;
          time_write += s.stop() ;
          start += count ;
        }

        if(nadvance < p) {
          count = all_local_sizes[nadvance*block_sizes.size()+i]*bs ;
          if(count != 0) {
            int id = 0 ;
            static int flag = 0 ;

            MPI_Irecv(&buffer[id][0],sizeof(T)*count,MPI_BYTE,nadvance,0,MPI_COMM_WORLD,
                      &recv_Requests[id]) ;
            //      MPI_Isend(&flag,1,MPI_INT,nadvance,0,MPI_COMM_WORLD,&send_Requests[id]) ;
            MPI_Send(&flag,1,MPI_INT,nadvance,0,MPI_COMM_WORLD) ;
          }
        }
        for(int k=1;k<p;++k) {
          count = all_local_sizes[k*block_sizes.size()+i]*bs ;
          if(count != 0) {
            int id = k%nadvance ;
            MPI_Status mstat ;
            //      MPI_Wait(&send_Requests[id],&mstat) ;
            //      int flag = 0 ;
            sw.start() ;

            MPI_Wait(&recv_Requests[id],&mstat) ;

            time_wait += sw.stop() ;

            Loci::stopWatch s ;
            s.start() ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                &start,&stride,&count, NULL) ;
            hid_t memspace = H5Screate_simple(rank,&count,NULL) ;
            //      MPI_Status mstat ;

            //      MPI_Recv(&buffer[id][0],sizeof(T)*count,MPI_BYTE,k,0,MPI_COMM_WORLD,
            //               &mstat) ;

            H5Dwrite(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&buffer[id][0]) ;
            H5Sclose(memspace) ;
            time_write += s.stop() ;
            start += count ;
          }
	  int k_new = k+nadvance ;
	  if(k_new < p) {
	    int id = k_new%nadvance ;
	    int lcount = all_local_sizes[k_new*block_sizes.size()+i]*bs ;
	    if(lcount != 0) {
	      static int flag = 0 ;

	      
	      MPI_Irecv(&buffer[id][0],sizeof(T)*lcount,MPI_BYTE,k_new,0,
			MPI_COMM_WORLD, &recv_Requests[id]) ;

	      MPI_Send(&flag,1,MPI_INT,k_new,0,MPI_COMM_WORLD) ;
	    }

	  }
	}
      } else {
	if(num_local_xfers[i] != 0) {
	  int sz = num_local_xfers[i]*bs ;
	  int bcnt = 0 ;
	  FORALL(dom,ii) {
	    int lsz = var.vec_size(ii) ;
	    if(lsz >= loc+bs) {
	      const T *base = var.begin(ii) ;
	      for(int j=0;j<bs;++j) {
		buffer[0][bcnt++] = base[loc+j] ;
	      }
	    }
	  } ENDFORALL ;
	  int flag = 0 ;
	  MPI_Status mstat ;
	  MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&mstat) ;
	  MPI_Send(&buffer[0][0],sizeof(T)*sz,MPI_BYTE,0,0,MPI_COMM_WORLD) ;
	}
      }
      loc += bs ;
    }
    if(r == 0)  {
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      H5Gclose(group_id) ;
    }
  }

  //original design is process 0 gather all the data and write it out
  //since new hdf5 version collective data transfer will do the gathering
  //so just perform parallel writing directly
  template<class T> void writeUnorderedVectorP(hid_t group_id,
                                               const char *element_name,
                                               std::vector<T> &v,
                                               MPI_Comm prime_comm) {
    int mpi_size;
    int mpi_rank;
    MPI_Comm_rank(prime_comm, &mpi_rank);
    MPI_Comm_size(prime_comm, &mpi_size);

    //serial version
    if(mpi_size==1){
      hsize_t array_size_combined = v.size() ;
      if(array_size_combined == 0)
        return ;

      int rank = 1 ;
      hsize_t dimension = array_size_combined ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = v.size() ;
      hid_t datatype = dp->get_hdf5_type() ;

      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

      if(count != 0) {
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,datatype,memspace,dataspace,H5P_DEFAULT, &v[0]) ;
        H5Sclose(memspace) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      return;
    }


    //each process figure out the size and the start position to write
    int local_size = v.size() ; //my size to write
    hsize_t rsize = local_size; //individual size for each process in prime_comm
    std::vector<hsize_t> prime_count(mpi_size);
    MPI_Allgather(&rsize,sizeof(hsize_t),MPI_BYTE,
                  &prime_count[0],sizeof(hsize_t),MPI_BYTE,prime_comm) ;

    std::vector<hsize_t> pdispls(mpi_size) ; //the start point of each process in prime_comm
    pdispls[0] = 0 ;
    for(int i = 1; i < mpi_size; i++) {
      pdispls[i] = pdispls[i-1]+prime_count[i-1] ;
    }
    hsize_t array_size = pdispls[mpi_size-1]+prime_count[mpi_size-1] ; //size of the whole dataset

    if(array_size == 0)
      return ;


    //create dataset collectively
    int rank = 1 ;
    hsize_t dimension = array_size ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    typedef data_schema_traits<T> traits_type ;
    DatatypeP dp = traits_type::get_type() ;
    hsize_t start = pdispls[mpi_rank] ;
    hsize_t stride = 1 ;

    //write data
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t dataset = H5Dcreate2(group_id,element_name,datatype,
                               dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

    herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                     &start,&stride,&rsize, NULL) ;
    WARN(ret < 0) ;

    /* create a memory dataspace independently */
    hid_t memspace = H5Screate_simple (rank, &rsize, NULL);
    WARN(memspace < 0) ;

    hid_t xfer_plist = create_xfer_plist(hdf5_const::dxfer_coll_type);

    H5Dwrite(dataset,datatype,memspace,dataspace, xfer_plist,  &v[0]) ;

    H5Pclose(xfer_plist) ;
    H5Sclose(memspace) ;

    H5Sclose(dataspace) ;
    H5Tclose(datatype) ;
    H5Dclose(dataset) ;
  }




  //Since file_id and group_id is open for prime_comm, even just one process
  //write out data in serial, the dataset should be created by all processes in
  //prim_comm, otherwise the file won't be able to close properly
  template<class T> void writeVectorSerialP(hid_t group_id,
                                            const char *element_name,
                                            std::vector<T> &v ,
                                            MPI_Comm prime_comm) {
#ifndef H5_HAVE_PARALLEL
    writeVectorSerial(group_id,element_name,v,prime_comm) ;
#else

    hsize_t array_size_combined = v.size() ;

    if(array_size_combined == 0)
      return ;

    if(MPI_COMM_NULL != prime_comm){
      int rank = 1 ;
      hsize_t dimension = array_size_combined ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = v.size() ;
      hid_t datatype = dp->get_hdf5_type() ;

      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

      int prime_rank = -1;

      MPI_Comm_rank(prime_comm, &prime_rank);
      if(prime_rank==0){
        if(count != 0) {
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,datatype,memspace,dataspace,
                   H5P_DEFAULT, &v[0]) ;
          H5Sclose(memspace) ;
        }
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
    }
    return;
#endif
  }

  template<class T> void readVectorSerial(hid_t group_id,
                                          const char *element_name,
                                          std::vector<T> &v ) {
    hsize_t    dimension;
    hid_t      dataset, dataspace;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
    dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

      hid_t datatype = dp->get_hdf5_type() ;

      {
        std::vector<T> tmp(dimension) ;
        v.swap(tmp) ;
      }

      H5Dread(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT, &v[0]);

      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Tclose(datatype) ;
    } else {
      std::vector<T> tmp ;
      v.swap(tmp) ;
    }
  }



  //all processes in prim_comm open the data set, prime_rank 0 read all the data
  template<class T> void readVectorSerialP(hid_t group_id,
                                           const char *element_name,
                                           std::vector<T> &v,
                                           MPI_Comm prime_comm){
#ifndef H5_HAVE_PARALLEL
    readVectorSerialP(group_id,element_name,v,prime_comm) ;
#else
    hsize_t    dimension;
    hid_t      dataset, dataspace;
    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
    if (MPI_COMM_NULL != prime_comm){
      dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);

      if( dataset > 0 ) {
        int prime_rank = -1;
        MPI_Comm_rank(prime_comm, &prime_rank);
        if(prime_rank==0){
          dataspace  = H5Dget_space(dataset);
          H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
          typedef data_schema_traits<T> traits_type ;
          DatatypeP dp = traits_type::get_type() ;
          hid_t datatype = dp->get_hdf5_type() ;
          {
            std::vector<T> tmp(dimension) ;
            v.swap(tmp) ;
          }
          H5Dread(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT, &v[0]);
          H5Sclose(dataspace);
          H5Tclose(datatype) ;
        }
        H5Dclose(dataset);
      } else {
        std::vector<T> tmp ;
        v.swap(tmp) ;
      }
    }
#endif
  }

  template<class T> void readUnorderedVector(hid_t group_id,
                                             const char *element_name,
                                             std::vector<T> &v,
                                             MPI_Comm comm) {
    int my_rank = 0 ;
    MPI_Comm_rank(comm,&my_rank) ;
    int procs = 1 ;
    MPI_Comm_size(comm,&procs) ;

    if(procs == 1)
      return readVectorSerial(group_id,element_name,v) ;


    hsize_t    dimension=0;
    hid_t      dataset=0, dataspace=0 ;

    if(my_rank == 0) {
      H5Eset_auto (H5E_DEFAULT,NULL, NULL);
      dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
    }
    int check = 0 ;
    if(my_rank == 0 && dataset > 0)
      check = 1 ;
    MPI_Bcast(&check,1,MPI_INT,0,comm) ;

    if( check ) {
      if(my_rank == 0) {
        dataspace  = H5Dget_space(dataset);
        H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
      }
    }
    int array_size = dimension ;
    MPI_Bcast(&array_size,1,MPI_INT,0,comm) ;
    int delta = array_size/procs ;
    int rem = array_size%procs ;
    int loc_array_size = delta + ((my_rank<rem)?1:0);
    {
      std::vector<T> tmp(loc_array_size) ;
      v.swap(tmp) ;
    }

    if( check ) {
      if(my_rank == 0) {
        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;

        hid_t datatype = dp->get_hdf5_type() ;
        int rank = 1 ;
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = loc_array_size ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start,&stride,&count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

        H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &v[0]);
        H5Sclose(memspace) ;
        start += count ;
        for(int i=1;i<procs;++i) {
          count = delta + ((i<rem)?1:0) ;
          std::vector<T> tmp(count) ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start,&stride,&count, NULL) ;
          memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dread(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&tmp[0]) ;
          MPI_Send(&tmp[0],count*sizeof(T),MPI_BYTE,i,0,comm) ;
          H5Sclose(memspace) ;
          start += count ;
        }
        H5Tclose(datatype) ;
      } else {
        MPI_Status mstat ;
        MPI_Recv(&v[0],loc_array_size*sizeof(T),MPI_BYTE,0,0,comm,&mstat) ;
      }
    } else {
      std::vector<T> tmp ;
      v.swap(tmp) ;
    }
    if(my_rank == 0) {
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  template<class T> void readUnorderedVectorP(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v,
                                              MPI_Comm prime_comm
                                              ) {


#ifndef H5_HAVE_PARALLEL
    readUnorderedVector(group_id,element_name,v,prime_comm) ;
#else
    int procs = 1 ;
    MPI_Comm_size(MPI_COMM_WORLD,&procs) ;

    if(procs == 1)
      return readVectorSerial(group_id,element_name,v) ;


    int mpi_rank = -1;
    int mpi_size = -1;
    hsize_t    dimension=0;
    hid_t      dataset=0, dataspace=0 ;
    if (MPI_COMM_NULL != prime_comm){
      MPI_Comm_rank(prime_comm, &mpi_rank);
      MPI_Comm_size(prime_comm, &mpi_size);
      H5Eset_auto (H5E_DEFAULT,NULL, NULL);
      dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
      WARN(dataset<0) ;

      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int array_size = dimension ;
      int delta = array_size/mpi_size ;
      int rem = array_size%mpi_size ;
      int loc_array_size = delta + ((mpi_rank<rem)?1:0);

      {
        std::vector<T> tmp(loc_array_size) ;
        v.swap(tmp) ;
      }

      //parallel read here
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      hid_t datatype = dp->get_hdf5_type() ;
      int rank = 1 ;
      hsize_t start = (mpi_rank< rem)?mpi_rank*(delta+1):(delta+1)*rem+delta*(mpi_rank-rem) ;
      hsize_t stride = 1 ;
      hsize_t count = loc_array_size ;

      herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                       &start,&stride,&count, NULL) ;
      WARN(ret<0) ;
      hid_t memspace = H5Screate_simple(rank, &count, NULL);
      WARN(memspace<0) ;

      H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &v[0]);

      H5Sclose(memspace) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      H5Dclose(dataset) ;
    }
#endif
  }

  template< class T >
    inline void writeMultiStoreP(hid_t file_id, std::string vname,
                                 multiStore<T> &var, entitySet write_set,
                                 fact_db &facts) {
    const_multiStore<T> var_const ;
    var_const = var.Rep() ;
#ifndef H5_HAVE_PARALLEL
    writeMultiStore(file_id,vname,var_const,write_set,facts) ;
#else
    writeMultiStoreP(file_id,vname,var_const,write_set, facts) ;
#endif
  }

  template< class T >
    inline void writeMultiStoreP(hid_t file_id, std::string vname,
                                 const const_multiStore<T> &var,
                                 entitySet write_set, fact_db &facts) {

#ifndef H5_HAVE_PARALLEL
    writeMultiStore(file_id,vname,var,write_set,facts) ; ;
#else
    
    //open the group
    hid_t group_id = -1 ;
    group_id = H5Gcreate(file_id,vname.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    std::vector<int> sizes_local ;//unique sizes on this process, such as 0, 35, 105,...

    double time_write = 0 ;
    double pre_time = 0 ;

    entitySet dom = write_set ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    //just write out the entities belong to me???
    if(df != 0)
      dom  = write_set & df->my_entities ;

    Loci::stopWatch sp ;

    sp.start() ;
    { // Write out book keeping stuff.
      std::vector<int> counts(dom.size()) ;
      int cnt = 0 ;
      std::vector<int> fileids(dom.size()) ;
      if(df == 0) {
        FORALL(dom,ii) {

          counts[cnt] = var.vec_size(ii) ;
          int file_no = ii ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;

      } else {
        Map l2g ;
        l2g = df->l2g.Rep() ;
        dMap g2f ;
        g2f = df->g2f.Rep() ;

        FORALL(dom,ii) {
          counts[cnt] = var.vec_size(ii) ;

          int file_no = g2f[l2g[ii]] ;
          fileids[cnt] = file_no ;
          cnt++ ;
        } ENDFORALL ;

      }

      writeUnorderedVectorP(group_id,"counts",counts,MPI_COMM_WORLD) ;
      writeUnorderedVectorP(group_id,"fileID",fileids,MPI_COMM_WORLD) ;
      std::sort(counts.begin(),counts.end()) ;
      std::vector<int>::const_iterator lu = std::unique(counts.begin(),counts.end()) ;
      int nunique = lu-counts.begin() ;
      std::vector<int> tmp(nunique) ;
      sizes_local.swap(tmp) ;
      for(int i=0;i<nunique;++i)
        sizes_local[i] = counts[i] ;
    }

    // schedule blocking for writing

    // Find global sizes array
    sizes_local.push_back(std::numeric_limits<int>::max()) ;
    int iloc = 0 ;
    std::vector<int> sizes ; //unique sizes for all processors
    int cmin = -1 ;
    do {
      MPI_Allreduce(&sizes_local[iloc],&cmin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      if(cmin != std::numeric_limits<int>::max())
        sizes.push_back(cmin) ;
      if(sizes_local[iloc] == cmin)
        iloc++ ;
    } while (cmin != std::numeric_limits<int>::max()) ;

    std::vector<int> block_sizes ;//for all processes
    size_t cnt = 0 ;
    if(sizes[0] == 0)
      cnt++ ;
    const int blk_factor = 256 ; //actual block size is Min(internal block size, 256)
    int cloc = 0 ;
    int cloc_prev = 0 ;
    while(cnt < sizes.size()) {
      cloc += blk_factor ;
      cloc = min(cloc,sizes[cnt]) ;
      block_sizes.push_back(cloc-cloc_prev) ;
      cloc_prev = cloc ;
      if(cloc == sizes[cnt])
        cnt++ ;
    }
    if(block_sizes.size() == 0) {
      block_sizes.push_back(0) ;
    }
    std::vector<int> num_local_xfers(block_sizes.size(),0) ; //local::for each size, how many elements need xfer

    FORALL(dom,ii) {
      int lsz = var.vec_size(ii) ;
      int tot = 0 ;
      for(size_t i=0;i<block_sizes.size();++i) {
        tot += block_sizes[i] ;
        if(lsz >= tot && block_sizes[i] != 0)
          num_local_xfers[i]++ ;
      }
    } ENDFORALL ;

    std::vector<int> block_data_elems(block_sizes.size(),0) ;//global::for each block, how many elements need xfer
    MPI_Allreduce(&num_local_xfers[0],&block_data_elems[0],block_sizes.size(),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;

    // write out write schedule information
    writeVectorSerialP(group_id,"block_schedule",block_sizes, MPI_COMM_WORLD) ;
    writeVectorSerialP(group_id,"block_sets", block_data_elems, MPI_COMM_WORLD) ;

    pre_time += sp.stop() ;
    // Now write out the main data block


    size_t total_size = 0 ;//for all processes, all blocks
    for(size_t i=0;i<block_sizes.size();++i) {
      total_size += block_sizes[i]*block_data_elems[i] ;
    }



    int max_local_size = 0 ;
    for(size_t i=0;i<block_sizes.size();++i) {
      int local_size = block_sizes[i]*num_local_xfers[i] ;
      max_local_size = max(max_local_size,local_size) ;
    }
    //no data transfer, only know the max_local_size;


    // int max_size = max_local_size ;
    // MPI_Allreduce(&max_local_size,&max_size,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD)
    ;
    // if(r == 0)
    //   max_local_size = max_size ;

    int p = 0;
    int r = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    //allgather local sizes over
    std::vector<int> all_local_sizes(block_sizes.size()*p) ;
    MPI_Allgather(&num_local_xfers[0],block_sizes.size(),MPI_INT,
                  &all_local_sizes[0],block_sizes.size(),MPI_INT,
                  MPI_COMM_WORLD) ;

    hid_t dataspace =-1;
    hid_t dataset =-1;
    int rank = 1 ;

    herr_t ret;
    hid_t datatype =-1; //data type of T

    /* Create a large dataset for all processes  */
    {
      hsize_t dimension = total_size ; //total size for all processes, all blocks
      dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      WARN(dataspace<0) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      datatype = dp->get_hdf5_type() ;
      WARN(datatype<0) ;

      dataset = H5Dcreate2(group_id,"data",datatype,
                           dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

      WARN(dataset<0) ;

      H5Sclose(dataspace);
      H5Tclose(datatype) ;
    }

    //send_buffer  over all blocks
    std::vector<T> send_buffer(max_local_size) ;

    int loc = 0 ;
    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    for(size_t i=0;i<block_sizes.size();++i) { //write out data block by block
      int bs = block_sizes[i] ;

      //prepare send_buffer
      if(num_local_xfers[i] != 0) {
        int bcnt = 0 ;
        FORALL(dom,ii) {
          int lsz = var.vec_size(ii) ;
          if(lsz >= loc+bs) {
            const T *base = var.begin(ii) ;
            for(int j=0;j<bs;++j) {
              send_buffer[bcnt++] = base[loc+j] ;
            }
          }
        } ENDFORALL ;
      }
      loc += bs ;
      {
        int local_size =  all_local_sizes[r*block_sizes.size()+i]*bs;
        //start depend on prime_rank
        int offset = 0;
        for(int pr = 0; pr < r ; pr++){
          offset +=  all_local_sizes[pr*block_sizes.size()+i]*bs;
        }

        Loci::stopWatch s ;
        s.start() ;
        hsize_t prime_start = start + offset;
        //      if(local_size != 0) {
        {
          Loci::stopWatch s ;
          s.start() ;
          dataspace = H5Dget_space (dataset);
	  WARN(dataspace<0) ;
          hsize_t count = local_size;
          ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                    &prime_start,&stride,&count, NULL) ;
	  WARN(ret <0) ;

          /* create a memory dataspace independently */
          hid_t memspace = H5Screate_simple (rank, &count, NULL);
	  WARN(memspace<0) ;

          /* set up the collective transfer properties list */

          typedef data_schema_traits<T> traits_type ;
          DatatypeP dp = traits_type::get_type() ;
          datatype = dp->get_hdf5_type() ;
	  WARN(datatype<0) ;

          int lsz = (count==0)?0:1 ;
          int gsz  = lsz ;

          MPI_Allreduce(&lsz,&gsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
          if(gsz != 0) {
            ret = H5Dwrite(dataset, datatype, memspace, dataspace,
                           H5P_DEFAULT, &send_buffer[0]);
            WARN(ret<0) ;
          }
          H5Sclose(memspace) ;

          H5Sclose(dataspace) ;
          H5Tclose(datatype) ;
        }
        time_write += s.stop() ;
        start += block_sizes[i]*block_data_elems[i]; //how many total elements written for this block
      }
    }

    {
      H5Dclose(dataset) ;
      H5Gclose(group_id) ;
    }
#endif
  }

  template< class T > inline void readMultiStore(hid_t file_id,
                                                 std::string vname,
                                                 multiStore<T> &var,
                                                 entitySet read_set,
                                                 fact_db &facts) {
    hid_t group_id = -1;
    int r = MPI_rank ;
    if(0 == r) {
      group_id = H5Gopen(file_id, vname.c_str(),H5P_DEFAULT) ;
    }

    std::vector<int> counts ;
    std::vector<int> fileID ;
    readUnorderedVector(group_id,"counts",counts,MPI_COMM_WORLD) ;
    readUnorderedVector(group_id,"fileID",fileID,MPI_COMM_WORLD) ;

    // Now we need to create a mapping from fileID to local and processor
    // number
    std::vector<int> local_num(fileID.size()) ;
    std::vector<int> procID(fileID.size()) ;
    findMapping(local_num,procID,read_set,fileID,facts) ;

    const int p = MPI_processes ;
    std::vector<int> send_sz(p,0) ;
    std::vector<int> recv_sz(p,0) ;
    std::vector<int> recv_local_num ;
    distributeMapMultiStore(send_sz,recv_sz,recv_local_num,local_num,procID) ;
    std::vector<int> recv_count ;
    sendCounts(recv_count,send_sz,recv_sz,recv_local_num,counts,procID) ;
    //    distributeMapMultiStore(send_sz,recv_sz,recv_count,counts,procID) ;

    std::vector<int> alloc_set = recv_local_num ;
    std::sort(alloc_set.begin(),alloc_set.end()) ;
    entitySet alloc ;
    for(size_t i=0;i<alloc_set.size();++i)
      alloc += alloc_set[i] ;

    store<int> count_store ;
    count_store.allocate(alloc) ;
    for(size_t i=0;i<recv_local_num.size();++i) {
      count_store[recv_local_num[i]] = recv_count[i] ;
    }
    var.allocate(count_store) ;

    // Now with the background work we can begin reading in the
    // multiStore in segments, but we need to get the schedule ready
    // schedule of the blocking sizes
    std::vector<int> block_schedule ;
    // sizes of the subarrays for each block
    std::vector<int> block_sets ;
    if(0 == r) {
      readVectorSerial(group_id,"block_schedule", block_schedule) ;
      readVectorSerial(group_id,"block_sets", block_sets) ;
    }
    int bsize = block_schedule.size() ;
    MPI_Bcast(&bsize,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(0 != r) {
      std::vector<int> bsched(bsize) ;
      block_schedule.swap(bsched) ;
      std::vector<int> bsets(bsize) ;
      block_sets.swap(bsets) ;
    }
    MPI_Bcast(&block_schedule[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&block_sets[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;

    // Open data array
    hid_t      dataset = 0, dataspace = 0;
    hsize_t dimension = 0 ;
    hid_t datatype = 0 ;
    if( 0 == r ) {
      dataset  = H5Dopen(group_id, "data", H5P_DEFAULT);
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      datatype = dp->get_hdf5_type() ;
    }


    // Loop over block schedule, communicate each block
    int rank = 1 ;
    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    hsize_t read_size = 0 ;
    int blksz = 0 ;
    std::vector<int> file_read_sizes(p,0) ;
    for(size_t i=0;i<block_sets.size();++i) {
      blksz += block_schedule[i] ;
      int lbsz = 0 ;
      std::vector<int> send_sz_blk(p,0) ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz && blksz > 0) {
          lbsz++ ;
          send_sz_blk[procID[j]]++ ;
        }

      MPI_Gather(&lbsz,1,MPI_INT,&file_read_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;
      int tot = 0 ;
      for(int j=0;j<p;++j)
        tot+= file_read_sizes[j] ;
      if(r == 0 && tot != block_sets[i]) {
        std::cerr << "inconsistent, tot = "
                  << tot << " block_sets[" << i << "]=" << block_sets[i]
                  << endl ;
      }

      int msgsz = lbsz*block_schedule[i] ;
      std::vector<T> read_buffer(msgsz) ;
      if(r == 0) {
        read_size = file_read_sizes[0]*block_schedule[i] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &read_size, NULL) ;
        hid_t memspace = H5Screate_simple(rank,&read_size,NULL) ;
        H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT,
                &read_buffer[0]) ;
        H5Sclose(memspace) ;
        start += read_size ;
        int max_send_size = 0 ;
        for(int j=1;j<p;++j)
          max_send_size = max(file_read_sizes[j],max_send_size) ;
        std::vector<T> send_buffer(max_send_size*block_schedule[i]) ;
        for(int j=1;j<p;++j) {
          read_size = file_read_sizes[j]*block_schedule[i] ;
          if(read_size > 0) {
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                &start, &stride, &read_size, NULL) ;
            hid_t memspace = H5Screate_simple(rank,&read_size,NULL) ;
            H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT,
                    &send_buffer[0]) ;
            H5Sclose(memspace) ;
            start += read_size ;

            MPI_Send(&send_buffer[0],read_size*sizeof(T),MPI_BYTE,j,0,
                     MPI_COMM_WORLD) ;
          }
        }

      } else {
        MPI_Status mstat ;
        if(msgsz > 0) {
          MPI_Recv(&read_buffer[0], msgsz*sizeof(T),MPI_BYTE,0,0,
                   MPI_COMM_WORLD,&mstat) ;
        }
      }

      // After reading in the data, send it to the appropriate processor
      //First figure out what needs to be sent

      // copy data to send ordered on processor
      std::vector<T> xmit_buffer(msgsz) ;
      std::vector<int> cnt(p,0) ;
      std::vector<int> xoffset(p+1,0) ;
      for(int j=0;j<p;++j)
        xoffset[j+1] = xoffset[j]+send_sz_blk[j] ;

      int rcnt = 0 ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz) {
          int ps = procID[j] ;
          int loc = xoffset[ps]+cnt[ps] ;
          int bs = block_schedule[i] ;
          for(int k=0;k<bs;++k)
            xmit_buffer[loc*bs+k] = read_buffer[rcnt++] ;
          cnt[ps]++ ;
        }


      std::vector<int> recv_sz_blk(p,0) ;
      int indx = 0 ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int c = count_store[recv_local_num[indx]] ;
          if(c >= blksz && c != 0) {
            recv_sz_blk[j]++ ;
          }
        }
      std::vector<int> roffset(p+1,0) ;
      for(int j=0;j<p;++j)
        roffset[j+1] = roffset[j]+recv_sz_blk[j] ;
      std::vector<T> rdata(roffset[p]*block_schedule[i]) ;
      int nreq = 0 ;
      for(int j=0;j<p;++j) {
        if(send_sz_blk[j] > 0)
          nreq++ ;
        if(recv_sz_blk[j] > 0)
          nreq++ ;
      }
      std::vector<MPI_Request> recv_Requests(nreq) ;
      int req = 0 ;
      for(int j=0;j<p;++j)
        if(recv_sz_blk[j] > 0) {
          MPI_Irecv(&rdata[roffset[j]*block_schedule[i]],
                    recv_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      for(int j=0;j<p;++j)
        if(send_sz_blk[j] > 0) {
          MPI_Isend(&xmit_buffer[xoffset[j]*block_schedule[i]],
                    send_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      std::vector<MPI_Status> statuslist(nreq) ;
      MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;

      // Now copy data into multistore
      indx = 0 ;
      int l = 0 ;
      int bs = block_schedule[i] ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int le = recv_local_num[indx] ;
          if(count_store[le] >= blksz) {
            for(int m=0;m<bs;++m) {
              var[le][blksz+m-bs] = rdata[l++] ;
            }
          }
        }
    }

    if(0 == r) {
      H5Tclose(datatype) ;
      H5Dclose(dataset) ;
      H5Gclose(group_id) ;
    }
  }


  template< class T > inline void readMultiStoreP(hid_t file_id,
                                                  std::string vname,
                                                  multiStore<T> &var,
                                                  entitySet read_set,
                                                  fact_db &facts) {
#ifndef H5_HAVE_PARALLEL
    readMultiStore(file_id,vname,var,read_set,facts) ;
#else
    //group_id is known to the world
    hid_t group_id = -1;
    int r = MPI_rank ;
    //int mpi_rank = MPI_rank;

    group_id = H5Gopen(file_id, vname.c_str(),H5P_DEFAULT) ;

    std::vector<int> counts ;
    std::vector<int> fileID ;
    readUnorderedVectorP(group_id,"counts",counts,MPI_COMM_WORLD) ;
    readUnorderedVectorP(group_id,"fileID",fileID,MPI_COMM_WORLD) ;

    // Now we need to create a mapping from fileID to local and processor
    // number
    std::vector<int> local_num(fileID.size()) ;
    std::vector<int> procID(fileID.size()) ;
    findMapping(local_num,procID,read_set,fileID,facts) ;

    const int p = MPI_processes ;
    std::vector<int> send_sz(p,0) ;
    std::vector<int> recv_sz(p,0) ;
    std::vector<int> recv_local_num ;
    distributeMapMultiStore(send_sz,recv_sz,recv_local_num,local_num,procID) ;
    std::vector<int> recv_count ;
    sendCounts(recv_count,send_sz,recv_sz,recv_local_num,counts,procID) ;
    //    distributeMapMultiStore(send_sz,recv_sz,recv_count,counts,procID) ;

    std::vector<int> alloc_set = recv_local_num ;
    std::sort(alloc_set.begin(),alloc_set.end()) ;
    entitySet alloc ;
    for(size_t i=0;i<alloc_set.size();++i)
      alloc += alloc_set[i] ;

    store<int> count_store ;
    count_store.allocate(alloc) ;
    for(size_t i=0;i<recv_local_num.size();++i) {
      count_store[recv_local_num[i]] = recv_count[i] ;
    }
    var.allocate(count_store) ;

    // Now with the background work done we can begin reading in the
    // multiStore in segments, but we need to get the schedule ready
    // schedule of the blocking sizes
    std::vector<int> block_schedule ;
    // sizes of the subarrays for each block
    std::vector<int> block_sets ;


    readVectorSerialP(group_id,"block_schedule", block_schedule, MPI_COMM_WORLD) ;
    readVectorSerialP(group_id,"block_sets", block_sets, MPI_COMM_WORLD) ;

    int bsize = block_schedule.size() ;
    MPI_Bcast(&bsize,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(0 != r) {
      std::vector<int> bsched(bsize) ;
      block_schedule.swap(bsched) ;
      std::vector<int> bsets(bsize) ;
      block_sets.swap(bsets) ;
    }
    MPI_Bcast(&block_schedule[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;
    MPI_Bcast(&block_sets[0],bsize,MPI_INT,0,MPI_COMM_WORLD) ;

    //prime processes Open data array
    hid_t      dataset = 0, dataspace = 0;
    hsize_t dimension = 0 ;
    hid_t datatype = 0 ;
    {
      dataset  = H5Dopen(group_id, "data", H5P_DEFAULT);
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      datatype = dp->get_hdf5_type() ;
    }

    // Loop over block schedule, communicate each block
    int rank = 1 ;
    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    //    hsize_t read_size = 0 ;
    herr_t ret;
    int blksz = 0 ;
    std::vector<int> file_read_sizes(p,0) ;
    for(size_t i=0;i<block_sets.size();++i) {
      //for each block
      blksz += block_schedule[i] ;
      int lbsz = 0 ;
      std::vector<int> send_sz_blk(p,0) ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz && blksz > 0) {
          lbsz++ ;
          send_sz_blk[procID[j]]++ ;
        }

      MPI_Allgather(&lbsz,1,MPI_INT,&file_read_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;

      int tot = 0 ;
      for(int j=0;j<p;++j)
        tot+= file_read_sizes[j] ;
      if(r == 0 && tot != block_sets[i]) {
        std::cerr << "inconsistent, tot = "
                  << tot << " block_sets[" << i << "]=" << block_sets[i]
                  << endl ;
      }

      int msgsz = lbsz*block_schedule[i] ;

      std::vector<T> read_buffer ;
      if(msgsz > 0) read_buffer.resize(msgsz);

      {
        int offset = 0;
        for(int pr = 0; pr < r ; pr++){
          offset += file_read_sizes[pr]*block_schedule[i] ;
        }

        hsize_t prime_start = start + offset;

        hsize_t count = msgsz;
        ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                  &prime_start,&stride,&count, NULL) ;
        WARN(ret<0) ;
        /* create a memory dataspace independently */
        hid_t memspace = H5Screate_simple (rank, &count, NULL);
        WARN(memspace<0) ;

        /* set up the collective transfer properties list */

        ret = H5Dread(dataset, datatype, memspace, dataspace,
                      H5P_DEFAULT, &read_buffer[0]);

        H5Sclose(memspace) ;

        start += block_sets[i]*block_schedule[i];////how many total elements read for this block
      }

      // After reading in the data, send it to the appropriate processor
      //First figure out what needs to be sent

      // copy data to send ordered on processor
      std::vector<T> xmit_buffer ; //give it a default capacity
      if(msgsz > 0)xmit_buffer.resize(msgsz);

      std::vector<int> cnt(p,0) ;
      std::vector<int> xoffset(p+1,0) ;
      for(int j=0;j<p;++j)
        xoffset[j+1] = xoffset[j]+send_sz_blk[j] ;

      int rcnt = 0 ;
      for(size_t j=0;j<counts.size();++j)
        if(counts[j] >=blksz) {
          int ps = procID[j] ;
          int loc = xoffset[ps]+cnt[ps] ;
          int bs = block_schedule[i] ;
          for(int k=0;k<bs;++k)
            xmit_buffer[loc*bs+k] = read_buffer[rcnt++] ;
          cnt[ps]++ ;
        }

      std::vector<int> recv_sz_blk(p,0) ;
      int indx = 0 ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int c = count_store[recv_local_num[indx]] ;
          if(c >= blksz && c != 0) {
            recv_sz_blk[j]++ ;
          }
        }
      std::vector<int> roffset(p+1,0) ;
      for(int j=0;j<p;++j)
        roffset[j+1] = roffset[j]+recv_sz_blk[j] ;
      std::vector<T> rdata(roffset[p]*block_schedule[i]) ;
      int nreq = 0 ;
      for(int j=0;j<p;++j) {
        if(send_sz_blk[j] > 0)
          nreq++ ;
        if(recv_sz_blk[j] > 0)
          nreq++ ;
      }

      std::vector<MPI_Request> recv_Requests(nreq) ;
      int req = 0 ;
      for(int j=0;j<p;++j)
        if(recv_sz_blk[j] > 0) {
          MPI_Irecv(&rdata[roffset[j]*block_schedule[i]],
                    recv_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      for(int j=0;j<p;++j)
        if(send_sz_blk[j] > 0) {
          MPI_Isend(&xmit_buffer[xoffset[j]*block_schedule[i]],
                    send_sz_blk[j]*block_schedule[i]*sizeof(T),
                    MPI_BYTE,j,3,
                    MPI_COMM_WORLD,&recv_Requests[req]) ;
          req++ ;
        }
      std::vector<MPI_Status> statuslist(nreq) ;
      MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;

      // Now copy data into multistore
      indx = 0 ;
      int l = 0 ;
      int bs = block_schedule[i] ;
      for(int j=0;j<p;++j)
        for(int k=0;k<recv_sz[j];++k,indx++) {
          int le = recv_local_num[indx] ;
          if(count_store[le] >= blksz) {
            for(int m=0;m<bs;++m) {
              var[le][blksz+m-bs] = rdata[l++] ;
            }
          }
        }
    }

    H5Tclose(datatype) ;
    H5Dclose(dataset) ;
    H5Gclose(group_id) ;
#endif
  }
}
#endif

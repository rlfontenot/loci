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
#include <hdf5_readwrite.h>
#include <string>
using std::string;
namespace Loci {


  /*
   * Create the appropriate File access property list
   */
  namespace hdf5_const {
    extern const int PPN(1);
    extern const int facc_type(FACC_MPIO);		/*Test file access type */
    //seems no big diff, collective io error: H5Dwrite() unable to write 
    extern const int dxfer_coll_type(DXFER_COLLECTIVE_IO); /*use collective IO*/
    // extern const int dxfer_coll_type(DXFER_INDEPENDENT_IO); /*use independent IO*/
  }

  int  create_mpi_info(MPI_Info *info){//create a info and set the
    // *info = MPI_INFO_NULL;
    //return 0;
      MPI_Info_create(info);
  
    //MPI predefined hints
    MPI_Info_set(*info, "striping_unit", "8388608") ;
    MPI_Info_set(*info,"stripping_factor","4");
    MPI_Info_set(*info,"cb_nodes","4");
    MPI_Info_set(*info,"cb_buffer_size","4194304");
   
    //new algorithm parameters
    MPI_Info_set(*info,"ind_rd_buffer_size","41943040"); 
    MPI_Info_set(*info,"ind_wr_buffer_size","5242880");
    //platform-specific hints
    MPI_Info_set(*info,"IBM_largeblock_io","true");
    MPI_Info_set(*info,"H5F_ACS_CORE_WRITE_TRACKING_PAGE_SIZE_DEF","524288");
    MPI_Info_set(*info,"romio_ds_read","disable"); //dissble datat sieving in read
    MPI_Info_set(*info,"romio_ds_write","disable"); // disable data sieving in write
    MPI_Info_set(*info,"romio_cb_write","enable"); //enable aggregation

    //Setting the environment variable
    //MPICH_MPIIO_HINTS_DISPLAY=1 to print out available I/O hints and their values
    //aggreation: processors with fast connection perform io for others
    return 0;
  }
  
  hid_t
  create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type)
  {
#ifdef H5_HAVE_PARALLEL
    hid_t ret_pl = -1;
    ret_pl = H5Pcreate (H5P_FILE_ACCESS);
    WARN(ret_pl<0) ;
    if (l_facc_type == FACC_DEFAULT){
      return (ret_pl);
    }
    if (l_facc_type == FACC_MPIO){
      /* set Parallel access with communicator */
      H5Pset_fapl_mpio(ret_pl, comm, info);
     
      
      //unleased collective metadata read/write options to avoid read storm
      //one process read small data and broadcast it to avoid I/O access
    
      // ret = H5Pset_all_coll_metadata_ops(ret_pl, TRUE);
      //WARN(ret < 0);
      //ret = H5Pset_coll_metadata_write(ret_pl, TRUE);
      //WARN(ret<0);
      return(ret_pl);
    }
    return(ret_pl);
#else
    return H5P_DEFAULT ;
#endif
  }

  hid_t
  create_xfer_plist(int l_xfer_type){
#ifdef H5_HAVE_PARALLEL
    hid_t xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    if(l_xfer_type == DXFER_COLLECTIVE_IO)
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    else
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    return xfer_plist;
#else
    return H5P_DEFAULT ;
#endif
  }

  using std::vector ;
  using std::ostringstream ;
 

 
  
  //**************************************************************************/
  void HDF5_WriteDomain(hid_t group_id, const entitySet &en )
  {
    hsize_t dimension=0;

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    hid_t dataspace, dataset;
    hid_t datatype  = H5T_NATIVE_INT;

    if( num_intervals < 1) return;

    dimension = num_intervals*2; //size of 1D Array
    if(dimension == 0) return;
    dataspace = H5Screate_simple(rank, &dimension, NULL);
#ifdef H5_USE_16_API
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT);
#else
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif

    interval *it = new interval[num_intervals];

    int *data = new int[num_intervals*2];      

    for(int i=0;i<num_intervals;i++){
      it[i]       = en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }

    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    delete [] data;
    delete [] it;

    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  
  //**************************************************************************/
  void HDF5_WriteDomainP(hid_t group_id, const entitySet &en , MPI_Comm comm)
  {//parallel io, all processes create dataset, but only process 0 write  

#ifndef H5_HAVE_PARALLEL
    HDF5_WriteDomain(group_id, en );
#else
    hsize_t dimension=0;

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    hid_t dataspace, dataset;
    hid_t datatype  = H5T_NATIVE_INT;
    int prank = 0 ;
    //int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    //MPI_Comm_size(comm,&np) ;

    if( num_intervals < 1) return;

    dimension = num_intervals*2; //size of 1D Array
    if(dimension == 0) return;
    dataspace = H5Screate_simple(rank, &dimension, NULL);
#ifdef H5_USE_16_API
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT);
#else
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif
  if(prank == 0){
    interval *it = new interval[num_intervals];
  
    int *data = new int[num_intervals*2];      

    for(int i=0;i<num_intervals;i++){
      it[i]       = en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }
  
    //no xfer_plist here, must be independently, which is default
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   

    delete [] data;
    delete [] it;
    }

    H5Sclose(dataspace);
    H5Dclose(dataset);
#endif
  }

  //**************************************************************************/

  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g)
  {
    int        indx = 0;
    hsize_t    dimension;
    hid_t      dataset, dataspace;

#ifdef H5_USE_16_API
    dataset    = H5Dopen(group_id, "Map");
#else
    dataset    = H5Dopen(group_id, "Map",H5P_DEFAULT);
#endif
    if( dataset > 0) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];

      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
               H5P_DEFAULT, data);

      l_g.allocate(eset);
      entitySet :: const_iterator ci;
      for( ci = eset.begin(); ci != eset.end(); ++ci) 
        l_g[*ci] = data[indx++];

      delete [] data;
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }
  //**************************************************************************/

  void HDF5_Local2GlobalP( hid_t group_id, entitySet &eset, Map &l_g)
  {
#ifndef H5_HAVE_PARALLEL
    HDF5_Local2Global( group_id, eset, l_g);
#else
    int        indx = 0;
    hsize_t    dimension;
    hid_t      dataset, dataspace;

#ifdef H5_USE_16_API
    dataset    = H5Dopen(group_id, "Map");
#else
    dataset    = H5Dopen(group_id, "Map",H5P_DEFAULT);
#endif
    if( dataset > 0) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];
      hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
               xfer_plist, data);

   
      l_g.allocate(eset);
      entitySet :: const_iterator ci;
      for( ci = eset.begin(); ci != eset.end(); ++ci) 
        l_g[*ci] = data[indx++];

      delete [] data;
    
      H5Pclose(xfer_plist);
      H5Sclose(dataspace);
      H5Dclose(dataset);

    }
#endif
  }

  //**************************************************************************/

  void HDF5_ReadDomain( hid_t group_id, entitySet &eset)
  {
    hsize_t    dimension;
    hid_t      dataset, dataspace;

#ifdef H5_USE_16_API
    H5Eset_auto (NULL, NULL);
#else
    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
#endif
    eset = EMPTY;
#ifdef H5_USE_16_API
    dataset  = H5Dopen(group_id, "Interval Set");
#else
    dataset  = H5Dopen(group_id, "Interval Set",H5P_DEFAULT);
#endif
    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];

      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
               H5P_DEFAULT, data);

      eset = EMPTY;
      for(size_t i=0;i< dimension;i+=2){
        eset |= interval(data[i],data[i+1]);
      }
      delete [] data;
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  //**************************************************************************/

  void HDF5_ReadDomainP( hid_t group_id, entitySet &eset)
  {
#ifndef H5_HAVE_PARALLEL
    HDF5_ReadDomain(group_id, eset);
#else
    hsize_t    dimension;
    hid_t      dataset, dataspace;

#ifdef H5_USE_16_API
    H5Eset_auto (NULL, NULL);
#else
    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
#endif
    eset = EMPTY;
#ifdef H5_USE_16_API
    dataset  = H5Dopen(group_id, "Interval Set");
#else
    dataset  = H5Dopen(group_id, "Interval Set",H5P_DEFAULT);
#endif
    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];
      hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
               xfer_plist, data);

      eset = EMPTY;
      for(size_t i=0;i< dimension;i+=2){
        eset |= interval(data[i],data[i+1]);
      }
      delete [] data;
      H5Pclose(xfer_plist);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
#endif
  }
  //**************************************************************************/

  void HDF5_WriteVecSize(hid_t group_id, const int &size )
  {

    hsize_t dimension = 1;
    int     rank = 1;
    if(dimension == 0) return;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
#ifdef H5_USE_16_API
    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                 H5P_DEFAULT);
#else
    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                 H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif

    H5Dwrite(vDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);

    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

  }

  //**************************************************************************/

  void HDF5_WriteVecSizeP(hid_t group_id, const int &size )
  {
#ifndef H5_HAVE_PARALLEL
    HDF5_WriteVecSize(group_id, size );
#else
    hsize_t dimension = 1;
    int     rank = 1;
    if(dimension == 0) return;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
#ifdef H5_USE_16_API
    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                 H5P_DEFAULT);
#else
    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                 H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif
    hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
    H5Dwrite(vDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, xfer_plist, &size);
    H5Pclose(xfer_plist);
    H5Sclose( vDataspace );
    H5Dclose( vDataset   );
#endif
  }
  //**************************************************************************/


  void HDF5_ReadVecSize(hid_t group_id, int *size )
  {
    int     rank=1;
    hsize_t dimension = 1;

#ifdef H5_USE_16_API
    H5Eset_auto (NULL, NULL);
#else
    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
#endif

    *size = 0;
  
#ifdef H5_USE_16_API
    hid_t vDataset   = H5Dopen( group_id, "VecSize");
#else
    hid_t vDataset   = H5Dopen( group_id, "VecSize", H5P_DEFAULT);
#endif
    if( vDataset > 0) {
      if(dimension == 0) return;
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDatatype  = H5T_NATIVE_INT;
  
      H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, size);

      H5Sclose( vDataspace );
      H5Dclose( vDataset   );
    }

  }
  //**************************************************************************/
  //**************************************************************************/


  void HDF5_ReadVecSizeP(hid_t group_id, int *size )
  {
#ifndef H5_HAVE_PARALLEL
    HDF5_ReadVecSize(group_id, size );
#else
    int     rank=1;
    hsize_t dimension = 1;

#ifdef H5_USE_16_API
    H5Eset_auto (NULL, NULL);
#else
    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
#endif

    *size = 0;
  
#ifdef H5_USE_16_API
    hid_t vDataset   = H5Dopen( group_id, "VecSize");
#else
    hid_t vDataset   = H5Dopen( group_id, "VecSize", H5P_DEFAULT);
#endif
    if( vDataset > 0) {
      if(dimension == 0) return;
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDatatype  = H5T_NATIVE_INT;
      hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
      H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, xfer_plist, size);
      H5Pclose(xfer_plist);
      H5Sclose( vDataspace );
      H5Dclose( vDataset   );
    }
#endif
  }
  
}

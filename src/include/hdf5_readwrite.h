#ifndef HDF15_READWRITE_H
#define HDF15_READWRITE_H

#include <entitySet.h>

namespace Loci {
  
  inline void HDF5_WriteDomain(hid_t group_id, const entitySet &en )
  {

    hsize_t dimension;

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    if( num_intervals < 1) return;

    dimension = num_intervals*2;        //  size of 1D Array 
    interval *it = new interval[num_intervals];

    int *data = new int[num_intervals*2];      

    for(int i=0;i<num_intervals;i++){
      it[i]       = en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }

    hid_t dataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t datatype  = H5T_NATIVE_INT;
    hid_t dataset   = H5Dcreate(group_id, "Domain", datatype, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    delete [] data;
    delete [] it;
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }


  inline void HDF5_ReadDomain( hid_t group_id, entitySet &eset)
  {

    hsize_t    dimension[1];

    hid_t dataset   = H5Dopen( group_id, "Domain");
    hid_t dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims (dataspace, dimension, NULL);

    int *data = new int[dimension[0]];
    H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
             H5P_DEFAULT, data);

    for(int i=0;i< dimension[0];i++){
      eset |= interval(data[i],data[i+1]);
      i++;
    }

    delete [] data;
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }


  inline void HDF5_WriteVecSize(hid_t group_id, const int &size )
  {

    hsize_t dimension = 1;
    int     rank = 1;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                 H5P_DEFAULT);
    H5Dwrite(vDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);

    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

  }


  inline void HDF5_ReadVecSize(hid_t group_id, int *size )
  {
    int     rank=1;
    hsize_t dimension = 1;

    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    hid_t vDataset   = H5Dopen( group_id, "VecSize");
    H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, size);
    H5Sclose( vDataspace );
    H5Tclose( vDatatype  );
    H5Dclose( vDataset   );

  }

}
#endif

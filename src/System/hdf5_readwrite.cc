#include <hdf5_readwrite.h>


namespace Loci {
  
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
  dataspace = H5Screate_simple(rank, &dimension, NULL);
  dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT);

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

  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g)
{
  int        indx = 0;
  hsize_t    dimension;
  hid_t      dataset, dataspace;

  dataset    = H5Dopen(group_id, "Map");
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

  void HDF5_ReadDomain( hid_t group_id, entitySet &eset)
{
  hsize_t    dimension;
  hid_t      dataset, dataspace;

  H5Eset_auto (NULL, NULL);

  eset = EMPTY;
  dataset  = H5Dopen(group_id, "Interval Set");
  if( dataset > 0 ) {
    dataspace  = H5Dget_space(dataset);
    H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

    int *data = new int[dimension];

    H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
             H5P_DEFAULT, data);

    eset = EMPTY;
    for(size_t i=0;i< dimension;i++){
      eset |= interval(data[i],data[i+1]);
      i++;
    }
    delete [] data;
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
}
  //**************************************************************************/

  void HDF5_WriteVecSize(hid_t group_id, const int &size )
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
  //**************************************************************************/


  void HDF5_ReadVecSize(hid_t group_id, int *size )
{
  int     rank=1;
  hsize_t dimension = 1;

  H5Eset_auto (NULL, NULL);

  *size = 0;
  
  hid_t vDataset   = H5Dopen( group_id, "VecSize");

  if( vDataset > 0) {
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
  
    H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, size);

    H5Sclose( vDataspace );
    H5Dclose( vDataset   );
  }

}
  //**************************************************************************/

}

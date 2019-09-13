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

}

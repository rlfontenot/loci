#ifndef HDF5_READWRITE_H
#define HDF5_READWRITE_H

#include <entitySet.h>
#include <hdf5CC/H5cpp.h>

namespace Loci {
  void HDF5_WriteDomain(H5::Group group, const entitySet &en );
  void HDF5_ReadDomain( H5::Group group, entitySet &eset);

  void HDF5_WriteVecSize(H5::Group group, int size );

  void HDF5_ReadVecSize(H5::Group group, int &size );
}

#endif

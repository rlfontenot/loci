#ifndef HDF5_READWRITE_H
#define HDF5_READWRITE_H

#include <entitySet.h>
#include <Map.h>
#include <data_traits.h>

namespace Loci {
  
  void HDF5_WriteDomain(hid_t group_id, const entitySet &en ) ;
  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g) ;
  void HDF5_ReadDomain( hid_t group_id, entitySet &eset) ;
  void HDF5_WriteVecSize(hid_t group_id, const int &size ) ;
  void HDF5_ReadVecSize(hid_t group_id, int *size ) ;

}
#endif

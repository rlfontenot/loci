#include <Loci_Datatypes.h>

#include <iostream>

using std::cerr ;
using std::endl ;

namespace Loci {
  hid_t   AtomicType::get_hdf5_type() const {
    switch( atom ) {
    case FLOAT:
      return (H5Tcopy(H5T_NATIVE_FLOAT));
    case DOUBLE:
      return (H5Tcopy(H5T_NATIVE_DOUBLE));
    case INT:
      return (H5Tcopy(H5T_NATIVE_INT));
    case CHAR:
      return (H5Tcopy(H5T_NATIVE_CHAR));
    default:
      cerr << "Unknown type" << endl;
    }
    return 0;
  }

  hid_t  ArrayType::get_hdf5_type() const {
    hid_t hdf5T = type_data->get_hdf5_type() ;
#ifdef HDF5V1_2
    size_t array_dims[10] ;
#else
    hsize_t array_dims[10];
#endif

    for(int k = 0; k < rank; k++)
      array_dims[k]  = dimension[k];
#ifdef HDF5V1_2
    hid_t vDatatype  = H5Tcreate( H5T_COMPOUND, numBytes);
    H5Tinsert_array(vDatatype, "Array", 0, 
                    rank, array_dims, NULL, hdf5T);
    return vDatatype ;
#else
    return   H5Tarray_create(hdf5T,rank,array_dims,NULL) ;
#endif
    }



  hid_t CompoundType::get_hdf5_type() const
  {
    hid_t vDatatype  = H5Tcreate( H5T_COMPOUND, numBytes);

    int rank = 1;
    hid_t hdf5T ;

    for( int i = 0; i < type_list.size(); i++) {
      hdf5T = type_list[i].type_data->get_hdf5_type() ;
      H5Tinsert(vDatatype,type_list[i].name.c_str(), type_list[i].offset,hdf5T) ;
    }
    return vDatatype ;
  }

}

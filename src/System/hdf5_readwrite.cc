#include <entitySet.h>
#include <hdf5CC/H5cpp.h>

namespace Loci {
  void HDF5_WriteDomain(H5::Group group, const entitySet &en )
  {

  //----------------------------------------------------------------------------
  // Write the domain ( or entitySet) in the HDF format. An entitySet consists 
  //of intervalset (vbegin,vend). which is converted into one dimensional array,
  // This one-dimensionl array is then written in HDF format.
  //---------------------------------------------------------------------------
  
    hsize_t dimension[1];

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    dimension[0] = num_intervals*2;        //  size of 1D Array 
    interval *it = new interval[num_intervals];

    int *data = new int[num_intervals*2];      

    for(int i=0;i<num_intervals;i++){
      it[i]=en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }

    try{
      H5::DataSpace dataspace( rank, dimension );
      H5::DataSet   dataset = group.createDataSet( "Domain", 
                                   H5::PredType::NATIVE_INT, dataspace );
      dataset.write( data, H5::PredType::NATIVE_INT );
      delete [] it;
    }

    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}
  }

  //******************************************************************************

  void HDF5_ReadDomain( H5::Group group, entitySet &eset)
  {

    hsize_t    dimension[1];

    try{
      H5::DataSet   dataset = group.openDataSet( "Domain");
      H5::DataSpace dataspace = dataset.getSpace();

      dataspace.getSimpleExtentDims( dimension, NULL);

      int *data = new int[dimension[0]];

      dataset.read( data, H5::PredType::NATIVE_INT );

      for(int i=0;i< dimension[0];i++){
          eset |= interval(data[i],data[i+1]);
          i++;
      }
      delete [] data;
    }

    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

  }

  //******************************************************************************

  void HDF5_WriteVecSize(H5::Group group, int size )
  {

    hsize_t dimension[] = {1};

    int data[1];

    data[0] = size;

    int rank = 1;     
    try{
      H5::DataSpace dataspace( rank, dimension );
      H5::DataSet   dataset = group.createDataSet( "VecSize", 
                                   H5::PredType::NATIVE_INT, dataspace );
      dataset.write( data, H5::PredType::NATIVE_INT );
    }

    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}
  }

  //******************************************************************************

  void HDF5_ReadVecSize(H5::Group group, int &size )
  {

	  hsize_t dimension[1];
     int     data[1];
 
     dimension[0] = 1;

     H5::DataType  datatype  = H5::PredType::NATIVE_INT;
	  H5::DataSet   dataset   = group.openDataSet( "VecSize");
	  H5::DataSpace dataspace = dataset.getSpace();

	  dataspace.getSimpleExtentDims( dimension, NULL);

	  dataset.read( data, H5::PredType::NATIVE_INT );
     size = data[0];
  }

  //******************************************************************************
}

#ifndef HDF5_WRITE_TEMPLATE_H_
#define HDF5_WRITE_TEMPLATE_H_

#include <hdf5_memento.h>
#include <hdf5_traits.h>
#include <hdf5CC/H5cpp.h>
#include <Tools/stream.h>
#include <typeinfo.h>

namespace Loci {

  inline void HDF5_WriteDomain(H5::Group group, const entitySet& en )
  {

  //----------------------------------------------------------------------------
  // Write the domain ( or entitySet) in the HDF format. An entitySet consists 
  //of intervalset (vbegin,vend). which is converted into one dimensional array,
  // This one-dimensionl array is then written in HDF format.
  //---------------------------------------------------------------------------
  
    hsize_t dimension[1];

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    dimension[0] = num_intervals*2;              //  size of 1D Array 
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
  };

  //******************************************************************************

  inline void HDF5_ReadDomain( H5::Group group, entitySet &eset)
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

  };

  //******************************************************************************

  inline void HDF5_WriteVecSize(H5::Group group, int size ) 
  {

    hsize_t dimension[] = {1};

    int rank = 1;     

    try{
      H5::DataSpace dataspace( rank, dimension );
      H5::DataSet   dataset = group.createDataSet( "VecSize", 
                                   H5::PredType::NATIVE_INT, dataspace );
      dataset.write( &size, H5::PredType::NATIVE_INT );
    }

    catch( H5::HDF5DatasetInterfaceException error  ) {error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error) {error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ) {error.printerror();}
  };

  //******************************************************************************

  inline void HDF5_ReadVecSize(H5::Group group, int &size )
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
  };

  //******************************************************************************


/*
  template <> 
  inline int get_storeVec_size(H5::Group group,DEFAULT_CONVERTER g)
  {
    char ch;
    int size;
    try{
      H5::DataSet dataset_store = group.openDataSet( "store");
      H5::DataSpace dataspace_store = dataset_store.getSpace();
      hsize_t dims_store[1];
      dataspace_store.getSimpleExtentDims( dims_store, NULL);
      char* memento= new char[dims_store[0]];
      dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
      std::istringstream iss(memento);
     
      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
      
      entitySet e ;
      iss >> e ;
      iss >> size ;

      delete [] memento;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

    return size;
  };


  //******************************************************************************

  template <> 
  inline entitySet get_store_domain(H5::Group group,DEFAULT_CONVERTER g)
  {
    char ch;
    entitySet e;

    try{
      H5::DataSet dataset_store = group.openDataSet( "store");
      H5::DataSpace dataspace_store = dataset_store.getSpace();

      hsize_t dims_store[1];

      dataspace_store.getSimpleExtentDims( dims_store, NULL);
      char* data_store= new char[dims_store[0]];
      dataset_store.read( data_store, H5::PredType::NATIVE_CHAR );
      std::istringstream iss(data_store);

      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
      iss >> e ;

      delete [] data_store;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

    return e;
  }


  template <class W> 
  void store_hdf5read(H5::Group group,DEFAULT_CONVERTER g,W* base_ptr,entitySet &en)
  {
    char ch;
    entitySet e;
    try{
      H5::DataSet dataset_store = group.openDataSet( "store");
      H5::DataSpace dataspace_store = dataset_store.getSpace();

      hsize_t dims_store[1];
      dataspace_store.getSimpleExtentDims( dims_store, NULL);

      char* memento= new char[dims_store[0]];

      dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
      std::istringstream iss(memento);

      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }

      iss >> e ;
      FORALL(en,ii) {
        iss >> base_ptr[ii] ;
      } ENDFORALL ;
        
      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
      delete [] memento;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
  }
*/

  //******************************************************************************

}
#endif

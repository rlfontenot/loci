#ifndef HDF5_WRITE_TEMPLATE_H_
#define HDF5_WRITE_TEMPLATE_H_

#include <typeinfo>
#include <hdf5_memento.h>
#include <hdf5_traits.h>
#include <hdf5CC/H5cpp.h>
#include <hdf5_readwrite.h>


namespace Loci {

  //**************************************************************************/


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
*/


  //************************************************************************/

/*
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

  //**************************************************************************/

}
#endif

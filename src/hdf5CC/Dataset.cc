#include "Refcount.h"
#include "Exception.h"
#include "H5object.h"
#include "Dataspace.h"
#include "AbstractDs.h"
#include "Dataset.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

// Default constructor
DataSet::DataSet() : AbstractDs() {}

// Creates a copy of DataSet using an existing id
DataSet::DataSet( const hid_t dataset_id ) : AbstractDs()
{
   id = dataset_id;
}

// Copy constructor makes a copy of the original object by using base
// class' copy constructors
DataSet::DataSet( const DataSet& original ) : AbstractDs( original ) {}

// Gets a copy of the dataspace of this dataset
DataSpace DataSet::getSpace() const
{
   // Calls C function H5Dget_space to get the id of the dataspace
   hid_t dataspace_id = H5Dget_space( id );

   // If the dataspace id is invalid, throw an exception
   if( dataspace_id <= 0 )
   {
      throw HDF5DatasetInterfaceException();
   }
   //create dataspace object using the existing id then return the object
   DataSpace data_space( dataspace_id );
   return( data_space );
}

// This private member function calls the C API to get the identifier 
// of the datatype that is used by this dataset.  It is used
// by the various AbstractDs functions to get the specific datatype.
hid_t DataSet::p_getType() const
{
   hid_t type_id = H5Dget_type( id );
   if( type_id > 0 )
      return( type_id );
   else
   {
      throw HDF5DatasetInterfaceException();
   }
}

// Returns the amount of storage required for a dataset.  
hsize_t DataSet::getStorageSize() const
{
   hsize_t storage_size = H5Dget_storage_size( id );

   if( storage_size >= 0 )
      return( storage_size );
   else
   {
      throw HDF5DatasetInterfaceException();
   }
}

// Reads raw data from the specified dataset into buf, converting from 
// file datatype and dataspace to memory datatype and dataspace. 
void DataSet::read( void * buf, const DataType& mem_type, const DataSpace& mem_space, const DataSpace& file_space, const hid_t& xfer_plist_id ) const
{
   // Obtain identifiers for C API
  hid_t mem_type_id = mem_type.getId();
   hid_t mem_space_id = mem_space.getId();
   hid_t file_space_id = file_space.getId();
   //hid_t xfer_plist_id = xfer_plist.getId();

   herr_t ret_value = H5Dread( id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf );
   if( ret_value < 0 )
   {
      throw HDF5DatasetInterfaceException();
   }
}

// Writes raw data from an application buffer buffer to a dataset, 
// converting from memory datatype and dataspace to file datatype 
// and dataspace.
void DataSet::write( const void * buf, const DataType& mem_type, const DataSpace& mem_space, const DataSpace& file_space, const hid_t& xfer_plist_id ) const
{
   // Obtain identifiers for C API
   hid_t mem_type_id = mem_type.getId();
   hid_t mem_space_id = mem_space.getId();
   hid_t file_space_id = file_space.getId();
   //hid_t xfer_plist_id = xfer_plist.getId();

   herr_t ret_value = H5Dwrite( id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf );
   if( ret_value < 0 )
   {
      throw HDF5DatasetInterfaceException();
   }
}

// ~DataSet makes sure that all the copies of this dataset are destroyed
// before invoking the C API H5Dclose to close the dataset
DataSet::~DataSet()
{
   // if this is the only copy, then call the C routine H5Dclose 
   // to close the dataset
   if( copy_count->getCounter() == 1 )
   {
      herr_t ret_value = H5Dclose( id );
      if( ret_value < 0 )  // raise exception when H5Dclose returns a neg value
      {
         throw HDF5DatasetInterfaceException();
      }
   }  // end of copy_count
}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif

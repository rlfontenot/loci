//H5file.cc
#include "Exception.h"
#include "H5file.h"
#include "Group.h"
#include "Dataset.h"
#include "H5templates.h"

namespace H5 {

  H5File::H5File( const char* filename, unsigned int flags, const hid_t& create_plist_id, const hid_t& access_plist_id ){
    char buf[512];
    sprintf(buf,"%s.hdf5",filename);
    // These bits only set for creation, so if any of them are set,
    // create the file.
    if( flags & (H5F_ACC_EXCL|H5F_ACC_TRUNC|H5F_ACC_DEBUG ))
      {
	id = H5Fcreate( buf, flags, create_plist_id, access_plist_id );
      }
    // Open the file if none of the bits above are set.
    else
      {
	id = H5Fopen( buf, flags, access_plist_id );
      }
    if( id <= 0 )  // throw an exception when open/create fail
      {
	throw HDF5FileInterfaceException();
      }

    // keep track of number of copies of this file to properly close
    copy_count = new ReferenceCounter;
  }

  // Copy constructor makes a copy of the original H5File object.
  H5File::H5File( const H5File& original ){
    copy_count = original.copy_count;  // points to the same copy count
    copy_count->increment();     // increment copy count
    id = original.id;
  }

  // Sets and gets H5File's data member id
  void H5File::setId( hid_t new_file_id ){
    id = new_file_id;
  }
  hid_t H5File::getId() const{
    return id;
  }

  // Creates a new group in this file using the template function provided
  // in FGtemplates.h  
  Group H5File::createGroup( const string name, size_t size_hint ) const{
    try {
      Group group = createGroupT<H5File>( *this, name, size_hint );
      return( group );
    }
    catch( File_GroupException error )
      {
	throw HDF5FileInterfaceException();
      }
  }

  // Opens an existing group in this file using the template function provided
  // in FGtemplates.h  
  Group H5File::openGroup( const string name ) const{
    try {
      Group group = openGroupT<H5File>( *this, name );
      return( group );
    }
    catch( File_GroupException error )
      {
	throw HDF5FileInterfaceException();
      }
  }

  // This destructor calls the C API H5Fclose to close this file
  H5File::~H5File() {
    // If this file does not have any more references, then
    // call the C routine H5Fclose to close it
    if( copy_count->noReference())
      {
	delete copy_count;  // deallocate the copy count for this file
	herr_t ret_value = H5Fclose( id );
	if( ret_value < 0 )
	  {
	    throw HDF5FileInterfaceException();
	  }
      }
  }

  // Creates a dataset in this file using the template function
  // provided in FGtemplates.h  
  DataSet H5File::createDataSet( const string name, const DataType& data_type, const DataSpace& data_space, const hid_t& create_plist_id ) const{
    try {
      DataSet dataset = createDataSetT<H5File>( *this, name, data_type, data_space, create_plist_id );
      return( dataset );
    }
   catch( File_GroupException error )
     {
       throw HDF5FileInterfaceException();
     }
  }
  
  // Opens an existing dataset in this file using the template function
  // provided in FGtemplates.h  
  DataSet H5File::openDataSet( const string name ) const{
    try {
      DataSet dataset = openDataSetT<H5File>( *this, name );
      return( dataset );
    }
    catch( File_GroupException error )
      {
	throw HDF5FileInterfaceException();
      }
  }
}

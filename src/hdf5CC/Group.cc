//Group.cc
#include "Group.h"
#include "Dataspace.h"
#include "H5file.h"
#include "Exception.h"
#include "H5templates.h"

namespace H5 {

  // Default constructor
  Group::Group() : H5Object() {}

  
  // Copy constructor: makes a copy of the original Group object 
  Group::Group( const Group& original ) : H5Object( original ) {}
  
  // Creates a new group in this group using the common function
  // provided in FGtemplates.h.
  Group Group::createGroup( const string name, size_t size_hint ){
    try {
      Group group = createGroupT<Group>( *this, name, size_hint );
      return( group );
    }
    catch( File_GroupException error )
      {
	throw HDF5GroupInterfaceException();
      }
  }

  // Creates a copy of an existing Group using its id
  Group::Group( const hid_t group_id ) : H5Object( group_id ) {};

  // Opens an existing group in this group using the common function
  // provided in FGtemplates.h.
  Group Group::openGroup( const string name ){
    try {
      Group group = openGroupT<Group>( *this, name );
      return( group );
    }
    catch( File_GroupException error )
      {
	throw HDF5GroupInterfaceException();
      }
  }

  // Creates a dataset in this group using the common function
  // provided in FGtemplates.h  
  DataSet Group::createDataSet( const string name, const DataType& data_type, const DataSpace& data_space, const hid_t& create_plist_id ){
    try {
      DataSet dataset = createDataSetT<Group>( *this, name, data_type, data_space, create_plist_id );
      return( dataset );
    }
    catch( File_GroupException error )
      {
	throw HDF5GroupInterfaceException();
      }
  }
  
  // Opens a dataset in this group using the common function
  // provided in FGtemplates.h  
  DataSet Group::openDataSet( const string name ){
    try {
      DataSet dataset = openDataSetT<Group>( *this, name );
      return( dataset );
    }
    catch( File_GroupException error )
      {
	throw HDF5GroupInterfaceException();
      }
  }

  // This destructor calls the C API H5Gclose to close this group
  Group::~Group(){
    // if this is the only copy, then call the C routine H5Gclose
    // to close the group
    if( copy_count->getCounter() == 1 )
      {
	herr_t ret_value = H5Gclose( id );
	if( ret_value < 0 )
	  {
	    throw HDF5GroupInterfaceException();
	  }
      }  // end of copy_count
  }

}

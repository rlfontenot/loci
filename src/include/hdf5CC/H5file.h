//H5FILE.h
#ifndef _H5File_H
#define _H5File_H

#include "Group.h"
extern "C" {
#include <hdf5.h>
}
#include "Refcount.h"
#include "Dataset.h"
namespace H5 {

  class H5File {
  private:
    hid_t id;
    ReferenceCounter* copy_count; // keeps track of # of copies of an object

  public:
    // copy constructor: makes a copy of the original H5File object.
    H5File(const H5File& original );

    // Creates or opens an HDF5 file.
    H5File( const char* name, unsigned int flags,
	    const hid_t& create_plist_id = H5P_DEFAULT,
	    const hid_t& access_plist_id = H5P_DEFAULT );
    
    // Sets and gets H5File's data member
    void setId( hid_t new_file_id );
    hid_t getId() const;
    
    // Creates a new group in this file
    Group createGroup( const string name, size_t size_hint = 0 ) const;
    
    // Opens an existing group in this file
    Group openGroup( const string name ) const;
    
    // Creates a new dataset in this file
    DataSet createDataSet( const string name, const DataType& data_type, const DataSpace& data_space, const hid_t& create_plist_id = H5P_DEFAULT ) const;
    
    // Opens a existing dataset in this file
    DataSet openDataSet( const string name ) const;
    
    ~H5File();
  };//end of H5FILE

}

#endif


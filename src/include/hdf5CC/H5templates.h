#ifndef _FGtemplates_H
#define _FGtemplates_H

// There are a few comments that are common to most of the functions
// defined in this file so they are listed here.
// - when a failure returned by the C API, the template functions will 
//   throw an exception, called File_GroupException, so Group or File can 
//   catch it and throw the appropriate exception to the user's application,
//   i.e., HDF5GroupInterfaceException or HDF5FileInterfaceException.
//

namespace H5 {
// Creates a new group at this location which can be a file or another group.
template <class T>
Group createGroupT( const T& loc, const string name, size_t size_hint )
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   // Call C routine H5Gcreate to create the named group, giving the 
   // location id which can be a file id or a group id
   hid_t group_id = H5Gcreate( loc.getId(), name_C, size_hint );

   // If the group id is valid, create and return the Group object
   if( group_id > 0 )
   {
      Group group( group_id );
      return( group );
   }
   else
   {
      throw File_GroupException();
   }
}

// Opens an existing group in a location which can be a file or another group
template <class T>
Group openGroupT( const T& loc, const string name )
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   // Call C routine H5Gopen to open the named group, giving the 
   // location id which can be a file id or a group id
   hid_t group_id = H5Gopen( loc.getId(), name_C );

   // If the group id is valid, create and return the Group object
   if( group_id > 0 )
   {
      Group group( group_id );
      return( group );
   }
   else
   {
      throw File_GroupException();
   }
}

// Creates a new dataset at this location.
template <class T>
DataSet createDataSetT( const T& loc, const string name, const DataType& data_type, const DataSpace& data_space, const hid_t& create_plist_id )
{
   // Convert the dataset's name in C++ string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   // Obtain identifiers for C API
   hid_t type_id = data_type.getId();
   hid_t space_id = data_space.getId();
   //hid_t create_plist_id = create_plist.getId();

   // Call C routine H5Dcreate to create the named dataset
   hid_t dataset_id = H5Dcreate( loc.getId(), name_C, type_id, space_id, create_plist_id );

   // If the dataset id is valid, create and return the DataSet object
   if( dataset_id > 0 )
   {
      DataSet dataset( dataset_id );
      return( dataset );
   }
   else
   {
      throw File_GroupException();
   }
}

// Opens an existing dataset at this location.
template <class T>
DataSet openDataSetT( const T& loc, const string name )
{
   // Convert the dataset's name in C++ string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   // Call C function H5Dopen to open the specified dataset, giving
   // the location id and the dataset's name 
   hid_t dataset_id = H5Dopen( loc.getId(), name_C );

   // If the dataset id is valid, create and return the DataSet object
   if( dataset_id > 0 )
   {
      DataSet dataset( dataset_id );
      return( dataset );
   }
   else
   {
      throw File_GroupException();
   }
}

}
#endif

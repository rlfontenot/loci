#include "Datatype.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

DataType::DataType( const hid_t existing_id, bool predefined ) : H5Object( existing_id ), is_predtype( predefined ) {}

// Creates a datatype given its class and size
DataType::DataType( const H5T_class_t type_class, size_t size ) : H5Object()
{
   // Call C routine to create the new datatype
   id = H5Tcreate( type_class, size );
   if( id <= 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Default constructor
DataType::DataType() : H5Object() {}

// Copy constructor: makes a copy of this DataType object.
DataType::DataType( const DataType& original ) : H5Object( original  )
{
   is_predtype = original.is_predtype; // copy data member from original
}

// Copies an existing datatype to this datatype object
void DataType::copy( const DataType& like_type )
{
   // If there are more than one copy of this object, decrement 
   // the reference counter, then start a new counter and reuse the 
   // object for another datatype 
   if( copy_count->getCounter() > 1 )
   {
      copy_count->decrement();
      copy_count = new ReferenceCounter;
   }

   // Call C routine H5Tcopy to copy the datatype
   id = H5Tcopy( like_type.getId() );
   if( id <= 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }

   is_predtype = false;	// this is not a predefined datatype
}

// Determines whether two datatypes are the same. ???
bool DataType::operator==(const DataType& compared_type ) const
{
   // Call C routine H5Tequal to determines whether two datatype 
   // identifiers refer to the same datatype
   htri_t ret_value = H5Tequal( id, compared_type.getId() );
   if( ret_value > 0 )
      return true;
   else if( ret_value == 0 )
      return false;
   else
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Operates a user's function on each attribute of an object - commented
// out because it should use the one from H5Object; need to check 
// the parameter list???
//int DataType::iterate( unsigned * idx, H5A_operator_t op, void *op_data )
//{
   // Call C routine H5Aiterate to iterate the object's attributes
   //int ret_value = H5Aiterate( id, idx, op, op_data );
   //if( ret_value >= 0 )
      //return( ret_value );
   //else
   //{
      //throw HDF5DatatypeInterfaceException();
   //}
//}

// Creates a new variable-length datatype - Note: make it inheritance???
//DataType DataType::vlenCreate( const DataType& base_type )
//{
   // Call C routine to create a new VL datatype
   //hid_t type_id = H5Tvlen_create( id );
   //if( type_id > 0 )
      //id = type_id;
   //else
   //{
      //throw HDF5DatatypeInterfaceException();
   //}
//}

// Commits a transient datatype to a file, creating a new named datatype
void DataType::commit( H5Object& loc, const string name ) const
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   hid_t loc_id = loc.getId(); // get location id for C API

   // Call C routine to commit the transient datatype 
   herr_t ret_value = H5Tcommit( loc_id, name_C, id );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Determines whether a datatype is a named type or a transient type. 
bool DataType::committed() const
{
   // Call C function to determine if a datatype is a named one
   htri_t committed = H5Tcommitted( id );
   if( committed > 0 )
      return true;
   else if( committed == 0 )
      return false;
   else
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Finds a conversion function.
H5T_conv_t DataType::find( const DataType& dest, H5T_cdata_t **pcdata ) const
{
   // Call C routine to find the conversion function
   H5T_conv_t func = H5Tfind( id, dest.getId(), pcdata );
   if( func == NULL )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( func );
}

// Converts data from between specified datatypes. 
void DataType::convert( const DataType& dest, size_t nelmts, void *buf, void *background, hid_t& plist_id ) const
{
   // Get identifiers for C API
   hid_t dest_id = dest.getId();
   //hid_t plist_id = plist.getId();

   // Call C routine H5Tconvert to convert the data
   herr_t ret_value;
   ret_value = H5Tconvert( id, dest_id, nelmts, buf, background, plist_id );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Sets the overflow handler to a specified function. 
void DataType::setOverflow( H5T_overflow_t func ) const
{
   // Call C routine H5Tset_overflow to set the overflow handler
   herr_t ret_value = H5Tset_overflow( func );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Returns a pointer to the current global overflow function. 
H5T_overflow_t DataType::getOverflow(void) const
{
   return( H5Tget_overflow());  // C routine
   // NULL can be returned as well
}

// Locks a datatype. 
void DataType::lock() const
{
   // Call C routine to lock the datatype
   herr_t ret_value = H5Tlock( id );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Returns the datatype class identifier. 
H5T_class_t DataType::getClass() const
{
   H5T_class_t type_class = H5Tget_class( id );

   // Return datatype class identifier if successful
   if( type_class == H5T_NO_CLASS )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( type_class );
}

// Returns the size of a datatype. 
size_t DataType::getSize() const
{
   // Call C routine to get the datatype size
   size_t type_size = H5Tget_size( id );
   if( type_size <= 0 ) // Is 0 valid value ???
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( type_size );
}

// Returns the base datatype from which a datatype is derived. 
// - not correct yet
DataType DataType::getSuper() const
{
   // Call C routine to get the base datatype from which the specified
   // datatype is derived. 
   hid_t base_type_id = H5Tget_super( id );

   // If H5Tget_super returns a valid datatype id, create and return
   // the base type, otherwise, raise exception
   if( base_type_id > 0 )
   {
      DataType base_type( base_type_id );
      return( base_type );
   }
   else {}
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Registers the specified conversion function. 
void DataType::registerFunc( H5T_pers_t pers, const string name, const DataType& dest, H5T_conv_t func ) const
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   hid_t dest_id = dest.getId();  // get id of the destination datatype

   // Call C routine H5Tregister to register the conversion function
   herr_t ret_value = H5Tregister( pers, name_C, id, dest_id, func );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Removes a conversion function from all conversion paths. 
void DataType::unregister( H5T_pers_t pers, const string name, const DataType& dest, H5T_conv_t func ) const
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   hid_t dest_id = dest.getId();  // get id of the dest datatype for C API

   // Call C routine H5Tunregister to remove the conversion function 
   herr_t ret_value = H5Tunregister( pers, name_C, id, dest_id, func );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Tags an opaque datatype. 
void DataType::setTag( const string tag ) const
{
   // Convert string to C-string
   const char* tag_C;
   tag_C = tag.c_str();  // name_C refers to the contents of name as a C-str

   // Call C routine H5Tset_tag to tag an opaque datatype. 
   herr_t ret_value = H5Tset_tag( id, tag_C );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Gets the tag associated with an opaque datatype. 
string DataType::getTag() const
{
   char* tag_Cstr = H5Tget_tag( id );

   // if the tag C-string returned is not NULL, convert it to C++ string
   // and return it, otherwise, raise an exception
   if( tag_Cstr != NULL )
   {
      string tag = string( tag_Cstr );
      return( tag );
   }
   else
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// This destructor terminates access to the datatype
DataType::~DataType()
{
   // if this is the only copy, then call the C routine H5Tclose
   // to close the datatype
   if( copy_count->getCounter() == 1 )
   {
      if( is_predtype == false )
      {
         herr_t ret_value = H5Tclose( id );
         if( ret_value < 0 )
         {
            throw HDF5DatatypeInterfaceException();
         }
      }
   }  // end of copy_count
}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif

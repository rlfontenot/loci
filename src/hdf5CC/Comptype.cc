#include <string>

#ifndef H5_NO_NAMESPACE
using namespace std;
#endif

#include "Comptype.h"
#include "Exception.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

// Creates a new compound datatype
CompType::CompType( size_t size ) : DataType( H5T_COMPOUND, size ) {}

// Creates a compound datatype using an existing id
CompType::CompType( const hid_t existing_id ) : DataType( existing_id ) {}

// Default constructor
CompType::CompType() : DataType() {}

// Copy constructor: makes copy of the original CompType object
CompType::CompType( const CompType& original ) : DataType( original ) {}

// Gets the compound datatype of the specified dataset - reimplement this
CompType::CompType( const DataSet& dataset ) : DataType()
{
   // Calls C function H5Dget_type to get the id of the datatype
   id = H5Dget_type( dataset.getId() );

   // If the datatype id is invalid, throw exception
   if( id <= 0 )
   {
      throw HDF5DatasetInterfaceException();
   }
}

// Retrieves the number of members in this compound datatype. 
int CompType::getNmembers() const
{
   int num_members = H5Tget_nmembers( id );
   if( num_members < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( num_members );
}

// Retrieves the name of a member of this compound datatype. 
string CompType::getMemberName( int member_num ) const
{
   char* member_name_C = H5Tget_member_name( id, member_num );
   if( member_name_C == NULL )  // should this be returned also???
   {
      throw HDF5DatatypeInterfaceException();
   }
   string member_name = string( member_name_C );
   return( member_name );
}

// Retrieves the offset of a member of a compound datatype. 
size_t CompType::getMemberOffset( int member_num ) const
{
   size_t offset = H5Tget_member_offset( id, member_num );
   if( offset == 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( offset );
}
#ifdef OLD_HDF5
// Returns the dimensionality of the member. 
int CompType::getMemberDims( int member_num, size_t *dims, int *perm ) const
{
   int num_dims = H5Tget_member_dims( id, member_num, dims, perm );
   if( 0 > num_dims > 4 )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( num_dims );
}

#endif

// Gets the type class of the specified member.
H5T_class_t CompType::getMemberClass( int member_num ) const
{
   // get the member datatype first
   hid_t member_type_id = H5Tget_member_type( id, member_num );
   if( member_type_id <= 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }

   // then get its class
   H5T_class_t member_class = H5Tget_class( member_type_id );
   if( member_class == H5T_NO_CLASS )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( member_class );
}

// This private member function calls the C API to get the identifier 
// of the specified member.  It is used by the getMemberXxxType
// below for the sub-types.
hid_t CompType::p_getMemberType( int member_num ) const
{
   hid_t member_type_id = H5Tget_member_type( id, member_num );
   if( member_type_id > 0 )
      return( member_type_id );
   else
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Returns the datatype of the specified member in this compound datatype. 
DataType CompType::getMemberDataType( int member_num ) const
{
   DataType datatype( p_getMemberType( member_num )); 
   return( datatype );
}

EnumType CompType::getMemberEnumType( int member_num ) const
{
   EnumType enumtype( p_getMemberType( member_num )); 
   return( enumtype );
}

CompType CompType::getMemberCompType( int member_num ) const
{
   CompType comptype( p_getMemberType( member_num )); 
   return( comptype );
}

IntType CompType::getMemberIntType( int member_num ) const
{
   IntType inttype( p_getMemberType( member_num )); 
   return( inttype );
}

FloatType CompType::getMemberFloatType( int member_num ) const
{
   FloatType floatype( p_getMemberType( member_num )); 
   return( floatype );
}

StrType CompType::getMemberStrType( int member_num ) const
{
   StrType strtype( p_getMemberType( member_num )); 
   return( strtype );
}

/* old style of getMemberType - using overloads; new style above 
   returns the appropriate datatypes but has different named functions.
// Returns the datatype of the specified member in this compound datatype. 
// Several overloading of getMemberType are for different datatypes
void CompType::getMemberType( int member_num, EnumType& enumtype ) const
{
   p_getMemberType( member_num, enumtype ); 
}

void CompType::getMemberType( int member_num, CompType& comptype ) const
{
   p_getMemberType( member_num, comptype ); 
}

void CompType::getMemberType( int member_num, IntType& inttype ) const
{
   p_getMemberType( member_num, inttype ); 
}

void CompType::getMemberType( int member_num, FloatType& floatype ) const
{
   p_getMemberType( member_num, floatype ); 
}

void CompType::getMemberType( int member_num, StrType& strtype ) const
{
   p_getMemberType( member_num, strtype ); 
}
// end of overloading of getMemberType
*/

// Adds a new member to a compound datatype
void CompType::insertMember( const string name, size_t offset, const DataType& new_member ) const
{
   // Convert string to C-string
   const char* name_C;
   name_C = name.c_str();  // name_C refers to the contents of name as a C-str

   hid_t new_member_id = new_member.getId();  // get new_member id for C API

   // Call C routine H5Tinsert to add the new member
   herr_t ret_value = H5Tinsert( id, name_C, offset, new_member_id );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}
#ifdef OLD_HDF5
// Adds an array member to this compound datatype.
void CompType::insertMember( const string member_name, size_t offset, int ndims, const size_t* dim, const int* perm, const DataType& new_member ) const
{
   // Convert string to C-string
   const char* name_C;
   name_C = member_name.c_str();  // name_C refers to the contents 
				  // of member_name as a C-str
   hid_t new_member_id = new_member.getId();  // get new_member id for C API

   // Call C routine H5Tinsert_array to add an array new member
   herr_t ret_value = 
	H5Tinsert_array( id, name_C, offset, ndims, dim, perm, new_member_id );

   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}
#endif

// Recursively removes padding from within a compound datatype. 
void CompType::pack() const
{
   // Calls C routine H5Tpack to remove padding
   herr_t ret_value = H5Tpack( id );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// This destructor just invokes the base-class' destructor
CompType::~CompType() {}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif

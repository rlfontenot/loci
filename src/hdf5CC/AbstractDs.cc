#include "Exception.h"
#include "AbstractDs.h"

namespace H5 {

// Default constructor; used by sub-classes; simply invokes base-class
// default constructor
AbstractDs::AbstractDs() : H5Object() {}

// Copy constructor: makes copy of the original object; simply invokes
// base-class copy constructor.
AbstractDs::AbstractDs( const AbstractDs& original ) : H5Object( original ) {}

// Returns the class of the datatype that is used by this dataset
H5T_class_t AbstractDs::getTypeClass() const
{
   // Gets the id of the datatype used by this dataset or attribute.
   // p_getType calls either H5Dget_type or H5Aget_type depending on 
   // which object invokes getTypeClass
   hid_t type_id = p_getType();

   // Gets the class of the datatype and validate it before returning
   H5T_class_t type_class = H5Tget_class( type_id );
   if( type_class != H5T_NO_CLASS )
      return( type_class );
   else
   {
     throw HDF5DatatypeInterfaceException();
   }
}

} // end namespace


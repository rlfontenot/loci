#include "Dataspace.h"
#include "Exception.h"

namespace H5 {

const DataSpace DataSpace::ALL( H5S_ALL );

// Default constructor
DataSpace::DataSpace() {}

// This constructor creates a DataSpace object, given a dataspace type
DataSpace::DataSpace( H5S_class_t type )
{
   id = H5Screate( type );
   if( id <= 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   // starts counting object references
   copy_count = new ReferenceCounter;
}

// Creates a new simple data space and opens it for access.
DataSpace::DataSpace( int rank, const hsize_t * dims, const hsize_t * maxdims)
{
   id = H5Screate_simple( rank, dims, maxdims );
   if( id <= 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   // starts counting object references
   copy_count = new ReferenceCounter;
}

// Uses an existing dataspace identifier to make a DataSpace object
// or uses a default id to create a default dataspace object
DataSpace::DataSpace( const hid_t space_id )
{
   // If a constant is being created, do not keep count.
   if( space_id == H5S_ALL )
      id = H5S_ALL;

   // making an object from an existing id, which must be closed later
   // with H5Sclose, so start counting
   else
   {
      copy_count = new ReferenceCounter;
      id = space_id;
   }
}

// Copy constructor: makes a copy of the original DataSpace object
DataSpace::DataSpace( const DataSpace& original )
{
   copy_count = original.copy_count;  // points to the same copy count
   copy_count->increment();     // increment copy count
   id = original.id;
}

void DataSpace::copy( const DataSpace& like_space )
{
   // If there are more than one copy of this object, decrement
   // the reference counter, then start a new counter and reuse the
   // object for another dataspace
   if( copy_count->getCounter() > 1 )
   {
      copy_count->decrement();
      copy_count = new ReferenceCounter;
   }
   // call C routine to copy the property list
   id = H5Scopy( like_space.getId() );
   if( id <= 0 )
   {
      throw HDF5PropertylistInterfaceException();
   }
   copy_count = new ReferenceCounter;
}

// Sets and gets DataSpace's data member id
void DataSpace::setId( hid_t new_space_id )
{
   id = new_space_id;
}
hid_t DataSpace::getId() const
{
   return id;
}

// Determines whether this dataspace is a simple dataspace.
bool DataSpace::isSimple () const
{
   htri_t simple = H5Sis_simple( id );
   if( simple > 0 )
      return true;
   else if( simple == 0 )
      return false;
   else
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Sets the offset of this simple dataspace.
void DataSpace::offsetSimple ( const hssize_t *offset ) const
{
   herr_t ret_value = H5Soffset_simple( id, offset );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Retrieves dataspace dimension size and maximum size
int DataSpace::getSimpleExtentDims ( hsize_t *dims, hsize_t *maxdims ) const
{
   int ndims = H5Sget_simple_extent_dims( id, dims, maxdims );
   if( ndims < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( ndims );
}

// Determines the dimensionality of a dataspace
int DataSpace::getSimpleExtentNdims () const
{
   int ndims = H5Sget_simple_extent_ndims( id );
   if( ndims < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( ndims );
}

// Determines the number of elements in a dataspace
hsize_t DataSpace::getSimpleExtentNpoints () const
{
   hsize_t num_elements = H5Sget_simple_extent_npoints( id );

   // num_elements = 0 when failure occurs
   if( num_elements > 0 )
      return( num_elements );
   else
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Determine the current class of a dataspace
H5S_class_t DataSpace::getSimpleExtentType () const
{
   H5S_class_t class_name = H5Sget_simple_extent_type( id );
   if( class_name == H5S_NO_CLASS )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( class_name );
}

// Copies the extent of a dataspace
void DataSpace::extentCopy ( DataSpace& dest_space ) const
{
   hid_t dest_space_id = dest_space.getId();
   herr_t ret_value = H5Sextent_copy( dest_space_id, id );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Sets or resets the size of an existing dataspace
void DataSpace::setExtentSimple( int rank, const hsize_t *current_size, const hsize_t *maximum_size ) const
{
   herr_t ret_value;
   ret_value = H5Sset_extent_simple( id, rank, current_size, maximum_size );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Removes the extent from a dataspace
void DataSpace::setExtentNone () const
{
   herr_t ret_value = H5Sset_extent_none( id );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}


// Determines the number of elements in a dataspace selection
hssize_t DataSpace::getSelectNpoints () const
{
   hssize_t num_elements = H5Sget_select_npoints( id );
   if( num_elements < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( num_elements );
}

// Get number of hyperslab blocks
hssize_t DataSpace::getSelectHyperNblocks () const
{
   hssize_t num_blocks = H5Sget_select_hyper_nblocks( id );
   if( num_blocks < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( num_blocks );
}

// Gets the list of hyperslab blocks currently selected
void DataSpace::getSelectHyperBlocklist( hsize_t startblock, hsize_t numblocks, hsize_t *buf ) const
{
   herr_t ret_value;
   ret_value = H5Sget_select_hyper_blocklist( id, startblock, numblocks, buf );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Gets the number of element points in the current selection
hssize_t DataSpace::getSelectElemNpoints () const
{
   hssize_t num_points = H5Sget_select_elem_npoints( id );
   if( num_points < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
   return( num_points );
}

// Gets the list of element points currently selected
void DataSpace::getSelectElemPointlist ( hsize_t startpoint, hsize_t numpoints, hsize_t *buf ) const
{
   herr_t ret_value;
   ret_value = H5Sget_select_elem_pointlist( id, startpoint, numpoints, buf );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Gets the bounding box containing the current selection
void DataSpace::getSelectBounds ( hsize_t *start, hsize_t *end ) const
{
   herr_t ret_value = H5Sget_select_bounds( id, start, end );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Selects array elements to be included in the selection for a dataspace
void DataSpace::selectElements ( H5S_seloper_t op, const size_t num_elements, const hssize_t *coord[ ] ) const
{
   herr_t ret_value;
   ret_value = H5Sselect_elements( id, op, num_elements, coord );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Selects the entire dataspace
void DataSpace::selectAll () const
{
   herr_t ret_value = H5Sselect_all( id );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

//Resets the selection region to include no elements
void DataSpace::selectNone () const
{
  herr_t ret_value = H5Sselect_none( id );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Verifies that the selection is within the extent of the dataspace
bool DataSpace::selectValid () const
{
  htri_t ret_value = H5Sselect_valid( id );
   if( ret_value > 0 )
      return true;
   else if( ret_value == 0 )
      return false;
   else
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Selects a hyperslab region to add to the current selected region
void DataSpace::selectHyperslab( H5S_seloper_t op, const hsize_t *count, const hssize_t *start, const hsize_t *stride, const hsize_t *block ) const
{
   herr_t ret_value;
   ret_value = H5Sselect_hyperslab( id, op, start, stride, count, block );
   if( ret_value < 0 )
   {
      throw HDF5DataspaceInterfaceException();
   }
}

// Destructor closes the HDF5 file
DataSpace::~DataSpace()
{
   // Only H5Sclose this dataspace if it is not a constant
   if( id != H5S_ALL ) // real dataspace, not a constant
   {
      // decrement the number of copies of this object and if this is the
      // only copy, then call the C routine H5Pclose to close the dataset
      if( copy_count->noReference())
      {
         delete copy_count;  // deallocate the copy count for this prop list
         herr_t ret_value = H5Sclose( id );
         if( ret_value < 0 )
         {
            throw HDF5DataspaceInterfaceException();
         }
      }
   }
}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif

#ifndef _H5DataSet_H
#define _H5DataSet_H

#include "Tools/stream.h"
extern "C" {
#include <hdf5.h>
}
#include "AbstractDs.h"
#include "Datatype.h"

namespace H5 {

class DataSet : public AbstractDs {
   private:
        // This function contains the common code that is used by
        // getTypeClass and various API functions getXxxType 
        // defined in AbstractDs for generic datatype and specific
        // sub-types
	virtual hid_t p_getType() const;

   public:
	// Default constructor
	DataSet();

	// Copy constructor
	DataSet( const DataSet& original );

	// Gets the dataspace of this dataset.
	virtual DataSpace getSpace() const;

	// Gets the storage size of this dataset.
	hsize_t getStorageSize() const;

	// Reads the data of this dataset and stores it in the provided buffer.
	// The memory and file dataspaces and the transferring property list
	// can be defaults.
	void read( void* buf, const DataType& mem_type, const DataSpace& mem_space = DataSpace::ALL, const DataSpace& file_space = DataSpace::ALL, const hid_t& xfer_plist = H5P_DEFAULT ) const;

	// Writes the buffered data to this dataset.
	// The memory and file dataspaces and the transferring property list
	// can be defaults.
	void write( const void* buf, const DataType& mem_type, const DataSpace& mem_space = DataSpace::ALL, const DataSpace& file_space = DataSpace::ALL, const hid_t& xfer_plist_id = H5P_DEFAULT ) const;

	virtual ~DataSet();

	// Creates a copy of an existing DataSet using its id
	// (used only by template functions in FGtemplates.h
	// to return a DataSet, will not be published; Note: should use
	// friend template function)
	DataSet( const hid_t dataset_id );

};
}
#endif

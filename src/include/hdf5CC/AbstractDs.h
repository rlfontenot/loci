#ifndef _AbstractDs_H
#define _AbstractDs_H

#include "Dataspace.h"
#include "H5object.h"
#include "Tools/stream.h"
extern "C" {
#include <hdf5.h>
}

namespace H5 {

class AbstractDs : public H5Object {

	// This member function is implemented by DataSet and Attribute
	virtual hid_t p_getType() const = 0;

   protected:

	// Default constructor
	AbstractDs();

   public:
	// Copy constructor
	AbstractDs( const AbstractDs& original );

	// Gets the dataspace of this abstract dataset - pure virtual
	virtual DataSpace getSpace() const = 0;

        // Gets the class of the datatype that is used by this abstract 
	// dataset        
	H5T_class_t getTypeClass() const;

	virtual ~AbstractDs() {}; 
};
}
#endif _AbstractDs_H


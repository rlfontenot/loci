#ifndef _H5Object_H
#define _H5Object_H

extern "C"{
#include <hdf5.h>
}
#include "Refcount.h"

namespace H5 {
class H5Object {
   protected:
	hid_t id;
	ReferenceCounter* copy_count; // used to keep track of the
	                              // number of copies of an object
	// Default constructor
	H5Object();

	// Creates a copy of an existing object giving the object id
	H5Object( const hid_t object_id );

   public:
	// Copy constructor: makes copy of an H5Object object.
	H5Object( const H5Object& original );

	// Sets and gets H5Object's data member
	void setId( hid_t new_id );
	hid_t getId () const;

	virtual ~H5Object();
}; /* end class H5Object */

}
#endif

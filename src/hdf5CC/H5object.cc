//H5object.cc

#include "H5object.h"

namespace H5 {

  // Default constructor
  H5Object::H5Object(){
    // starts counting object references
    copy_count = new ReferenceCounter;
  }

  // Constructs an object from an existing id
  H5Object::H5Object( hid_t object_id ) : id( object_id ){
    // starts counting object references
    copy_count = new ReferenceCounter;
  }
  
  // Copy constructor makes a copy of the original object
  H5Object::H5Object( const H5Object& original ){
    copy_count = original.copy_count;  // points to the same copy count
    copy_count->increment();     // increment copy count
    id = original.id;
  }

  // Sets and gets H5Object's data member id
  void H5Object::setId( hid_t new_object_id ){
    id = new_object_id;
  }
  hid_t H5Object::getId() const{
    return id;
  }

  H5Object::~H5Object(){
    // check if there are anymore references to this object and if it is the
    // only copy, then deallocate the reference counter for this object
    if( copy_count->noReference())
      delete copy_count;
  }

}

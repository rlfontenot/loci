//Group.h
#ifndef _Group_H
#define _Group_H

#include "Tools/stream.h"
extern "C" {
#include <hdf5.h>
}
#include "Dataset.h"
#include "H5object.h"

namespace H5 {

  class Group : public H5Object {
  public:
    // default constructor
    Group();
    
    // Copy constructor: makes a copy of the original object
    Group( const Group& original );

    // Creates a group in this group
    Group createGroup( const std::string name, size_t size_hint = 0 );
    
    // Opens an existing group in this group 
    Group openGroup( const std::string name );
  
    // Creates a dataset in this group
    DataSet createDataSet( const string name, const DataType& data_type, const DataSpace& data_space, const hid_t& create_plist_id = H5P_DEFAULT );
    
    // Opens a dataset in this group
    DataSet openDataSet( const string name );
    
    // Creates a copy of an existing Group using its id
    // (used only by template functions in FGtemplates.h
    // to return a Group; will not be published; maybe, use friend???)
    Group( const hid_t group_id );

    virtual ~Group();

  };
}
#endif

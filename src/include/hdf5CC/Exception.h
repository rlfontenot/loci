#ifndef _H5Exception_H
#define _H5Exception_H

#include "Tools/stream.h"

extern "C" {
#include <hdf5.h>
}

namespace H5 {

#include <string>
  using std::string ;


class HDF5Exception {
   private:
	string detailMessage;
   public:
	HDF5Exception();
	HDF5Exception( string message );
	string getDetailMesg() const;

	void printerror() const;
};

// This exception is privately used in Group and H5File only
class File_GroupException {
   public:
	File_GroupException();
	File_GroupException( string message );
};

class HDF5FileInterfaceException : public HDF5Exception {
   public:
	HDF5FileInterfaceException();
	HDF5FileInterfaceException( string message );
};

class HDF5GroupInterfaceException : public HDF5Exception {
   public:
	HDF5GroupInterfaceException();
	HDF5GroupInterfaceException( string message );
};

class HDF5ObjectHeaderException : public HDF5Exception {
   public:
	HDF5ObjectHeaderException();
	HDF5ObjectHeaderException( string message );
};

class HDF5DataspaceInterfaceException : public HDF5Exception {
   public:
	HDF5DataspaceInterfaceException();
	HDF5DataspaceInterfaceException( string message );
};

class HDF5DatatypeInterfaceException : public HDF5Exception {
   public:
	HDF5DatatypeInterfaceException();
	HDF5DatatypeInterfaceException( string message );
};

class HDF5PropertylistInterfaceException : public HDF5Exception {
   public:
	HDF5PropertylistInterfaceException();
	HDF5PropertylistInterfaceException( string message );
};

class HDF5DatasetInterfaceException : public HDF5Exception {
   public:
	HDF5DatasetInterfaceException();
	HDF5DatasetInterfaceException( string message );
};

class HDF5AttributeInterfaceException : public HDF5Exception {
   public:
	HDF5AttributeInterfaceException();
	HDF5AttributeInterfaceException( string message );
};

class HDF5FunctionArgumentException : public HDF5Exception {
   public:
	HDF5FunctionArgumentException();
	HDF5FunctionArgumentException( string message );
};

class HDF5ReferenceException : public HDF5Exception {
   public:
	HDF5ReferenceException();
	HDF5ReferenceException( string message );
};

class HDF5DataStorageException : public HDF5Exception {
   public:
	HDF5DataStorageException();
	HDF5DataStorageException( string message );
};

class HDF5LibraryInterfaceException : public HDF5Exception {
   public:
	HDF5LibraryInterfaceException();
	HDF5LibraryInterfaceException( string message );
};

}

#endif // _H5Exception_H

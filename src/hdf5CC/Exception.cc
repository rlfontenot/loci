#include "stream.h"
#include <hdf5.h>
#include <string>
#include "Exception.h"

namespace H5 {

HDF5Exception::HDF5Exception() : detailMessage("")
{
   printerror();
}

HDF5Exception::HDF5Exception( string message ) : detailMessage( message )
{
   printerror();
}

string HDF5Exception::getDetailMesg() const
{
   return( detailMessage );
}

// Prints information about an error
void HDF5Exception::printerror() const
{
   cout << "\nC++ API: <ERROR> ";
   cout << this->getDetailMesg() << endl;
}

HDF5FileInterfaceException::HDF5FileInterfaceException():HDF5Exception(){};
HDF5FileInterfaceException::HDF5FileInterfaceException( string message ): 
HDF5Exception( message ){};

HDF5GroupInterfaceException::HDF5GroupInterfaceException():HDF5Exception(){};
HDF5GroupInterfaceException::HDF5GroupInterfaceException( string message ): 
HDF5Exception( message ){};

HDF5ObjectHeaderException::HDF5ObjectHeaderException():HDF5Exception(){};
HDF5ObjectHeaderException::HDF5ObjectHeaderException( string message ): HDF5Exception( message ) {};

HDF5DataspaceInterfaceException::HDF5DataspaceInterfaceException():HDF5Exception(){};
HDF5DataspaceInterfaceException::HDF5DataspaceInterfaceException( string message ): HDF5Exception( message ) {};

HDF5DatatypeInterfaceException::HDF5DatatypeInterfaceException():HDF5Exception(){};
HDF5DatatypeInterfaceException::HDF5DatatypeInterfaceException( string message ): HDF5Exception( message ) {};

HDF5PropertylistInterfaceException::HDF5PropertylistInterfaceException():HDF5Exception(){};
HDF5PropertylistInterfaceException::HDF5PropertylistInterfaceException( string message ): HDF5Exception( message ) {};

HDF5DatasetInterfaceException::HDF5DatasetInterfaceException():HDF5Exception(){};
HDF5DatasetInterfaceException::HDF5DatasetInterfaceException( string message ): HDF5Exception( message ) {};

HDF5AttributeInterfaceException::HDF5AttributeInterfaceException():HDF5Exception(){};
HDF5AttributeInterfaceException::HDF5AttributeInterfaceException( string message ): HDF5Exception( message ) {};

HDF5FunctionArgumentException::HDF5FunctionArgumentException():HDF5Exception(){};
HDF5FunctionArgumentException::HDF5FunctionArgumentException( string message ): HDF5Exception( message ) {};

HDF5ReferenceException::HDF5ReferenceException():HDF5Exception(){};
HDF5ReferenceException::HDF5ReferenceException( string message ):
HDF5Exception( message ) {};

HDF5DataStorageException::HDF5DataStorageException():HDF5Exception(){};
HDF5DataStorageException::HDF5DataStorageException( string message ):
HDF5Exception( message ) {};

HDF5LibraryInterfaceException::HDF5LibraryInterfaceException():HDF5Exception(){};
HDF5LibraryInterfaceException::HDF5LibraryInterfaceException( string message ):
HDF5Exception( message ) {};

File_GroupException::File_GroupException()
{
   // for now, do nothing
}

File_GroupException::File_GroupException( string message )
{
   // for now, do nothing
}

} // end namespace


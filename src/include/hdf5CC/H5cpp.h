#ifndef _H5CPP_H
#define _H5CPP_H

#include "Alltypes.h"
#include "Refcount.h"
#include "Exception.h"
#include "Dataspace.h"
#include "H5object.h"
#include "AbstractDs.h"
#include "Dataset.h"
#include "Group.h"
#include "H5file.h"
extern "C"{
#include <hdf5.h>
#ifdef inline
#undef inline
#endif
}
#endif

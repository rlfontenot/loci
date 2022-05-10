//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef HDF5_READWRITE_H
#define HDF5_READWRITE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <entitySet.h>
#include <Map.h>
#include <data_traits.h>

#define DXFER_COLLECTIVE_IO 0x1  /* Collective IO*/
#define DXFER_INDEPENDENT_IO 0x2 /* Independent IO collectively */
//#define dxfer_coll_type 0x1 /* define the data transfer type as Collective IO*/

//from partest/testphdf5.h,version 1.10.2
#define FACC_DEFAULT    0x0     /* default */
#define FACC_MPIO       0x1     /* MPIO */

namespace Loci {
   extern MPI_Info PHDF5_MPI_Info ;

namespace hdf5_const {
extern const int PPN; //processes per node
extern const int facc_type; //file access type
extern const int dxfer_coll_type;/* define the data transfer type as Collective IO or Independent IO collectively*/
}



//serial io
void HDF5_WriteDomain(hid_t group_id, const entitySet &en ) ;
void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g) ;
void HDF5_ReadDomain( hid_t group_id, entitySet &eset) ;
void HDF5_WriteVecSize(hid_t group_id, const int &size ) ;
void HDF5_ReadVecSize(hid_t group_id, int *size ) ;
//parallel io
  void HDF5_WriteDomainP(hid_t group_id, const entitySet &en, MPI_Comm comm ) ;
void HDF5_Local2GlobalP( hid_t group_id, entitySet &eset, Map &l_g) ;
  void HDF5_ReadDomainP( hid_t group_id, entitySet &eset) ;
void HDF5_WriteVecSizeP(hid_t group_id, const int &size ) ;
void HDF5_ReadVecSizeP(hid_t group_id, int *size ) ;

/*
 * Create the appropriate File access property list
 */
hid_t create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type);
  
hid_t  create_xfer_plist(int l_xfer_type);
int create_mpi_info(MPI_Info *info);

 

}
#endif

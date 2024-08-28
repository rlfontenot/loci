/** ****************************************************************************
 * @file      hdf5_Functions.h
 * @author
 * @brief     This file...
 * @details   This file is part of the Loci Framework.
 *
 *            The Loci Framework is free software: you can redistribute it
 *            and/or modify it under the terms of the Lesser GNU General Public
 *            License as published by the Free Software Foundation, either
 *            version 3 of the License, or (at your option) any later version.
 *
 *            The Loci Framework is distributed in the hope that it will be
 *            useful, but WITHOUT ANY WARRANTY; without even the implied
 *            warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *            See the Lesser GNU General Public License for more details.
 *
 *            You should have received a copy of the Lesser GNU General Public
 *            License along with the Loci Framework.  If not, see
 *            <http://www.gnu.org/licenses>
 * @version   0.2
 * @date
 * @copyright Copyright (c) 2008, Mississippi State University
 * @defgroup  hdf5_Functions hdf5_Functions
 * @ingroup   extract
 ******************************************************************************/
#ifndef HDF5_FUNCTIONS_H
#define HDF5_FUNCTIONS_H

#include <Loci.h>

using std::vector;


/// @addtogroup hdf5_Functions
/// @{

#define MAX_NAME 1024 //!< Max HDF5 object name

/** ****************************************************************************
 * @brief 
 * @tparam T 
 * @param group_id 
 * @param element_name 
 * @param v 
 * @param start_elem 
 * @param block_size 
 ******************************************************************************/
template<class T>
void readElementTypeBlock(hid_t group_id, const char *element_name,
                          vector<T> &v, size_t start_elem, int block_size)
{
  if(block_size > 0) {
#ifdef H5_USE_16_API
    hid_t dataset = H5Dopen(group_id,element_name) ;
#else
    hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;
#endif

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;

    hsize_t count = block_size ;
    hsize_t start = start_elem ;
    hid_t dataspace = H5Dget_space(dataset) ;
    hsize_t stride = 1 ;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
    int rank = 1 ;
    hsize_t dimension = count ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    hid_t datatype = dp->get_hdf5_type() ;
    H5Dread(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&v[0]) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;
  }
} // End of readElementTypeBlock()


/** ****************************************************************************
 * @brief 
 * @tparam T 
 * @param group_id 
 * @param element_name 
 * @param v 
 ******************************************************************************/
template<class T>
void readElementType(hid_t group_id, const char *element_name, vector<T> &v)
{
  if(v.size() > 0) {
#ifdef H5_USE_16_API
    hid_t dataset = H5Dopen(group_id,element_name) ;
#else
    hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;
#endif
    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }
} // End of readElementType()


/** ****************************************************************************
 * @brief
 * @param group_id
 * @param element_name
 * @return size_t
 ******************************************************************************/
inline size_t sizeElementType(hid_t group_id, const char *element_name)
{
#ifdef H5_USE_16_API
  hid_t dataset = H5Dopen(group_id,element_name);
#else
  hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT);
#endif

  if(dataset < 0)
  {
#ifdef H5_USE_16_API
    H5Eclear();
#else
    H5Eclear(H5E_DEFAULT);
#endif
    return 0;
  } // End If(dataset)

  hid_t   dspace = H5Dget_space(dataset);
  hsize_t size   = 0;
  H5Sget_simple_extent_dims(dspace,&size,NULL);
  H5Dclose(dataset);

  return size;
} // End of sizeElementType()


/** ****************************************************************************
 * @brief
 * @param file_id
 * @param varName
 * @return string
 ******************************************************************************/
inline string getVarNameFromFile(hid_t file_id, string varName)
{
  string  name = varName;
  hid_t   grp  = H5Gopen(file_id,"/", H5P_DEFAULT);
  hsize_t nobj;
  H5Gget_num_objs(grp,&nobj);
  if(nobj == 1)
  {
    char    memb_name[MAX_NAME];
    ssize_t len = H5Gget_objname_by_idx(grp, (hsize_t)0, memb_name, (size_t)MAX_NAME);
    name = string(memb_name);
  }

  H5Gclose(grp);
  return name;
} // End of getVarNameFromFile()


/** ****************************************************************************
 * @brief
 * @param file_id
 * @param vname
 * @param var
 * @param readSet
 * @param facts
 ******************************************************************************/
inline void readData(hid_t file_id, std::string vname, Loci::storeRepP  var,
                     entitySet readSet, fact_db &facts)
{
  hid_t   grp = H5Gopen(file_id,"/", H5P_DEFAULT);
  hsize_t nobj;

  H5Gget_num_objs(grp, &nobj);

  if(nobj == 1)
  {
    char    memb_name[MAX_NAME];
    ssize_t len = H5Gget_objname_by_idx(grp, (hsize_t)0, memb_name, (size_t)MAX_NAME);
    string  name(memb_name);

#ifdef VERBOSE
    if(name != vname)
    {
      cerr << "NOTE: reading dataset '" << name << "' instead of '"
           << vname << "'" << endl;
    }
#endif
    H5Gclose(grp);
    Loci::readContainer(file_id,name,var,readSet,facts);
  }else
  {
    H5Gclose(grp);
    Loci::readContainer(file_id,vname,var,readSet,facts);
  }
} // End of readData()

/// @}

#endif
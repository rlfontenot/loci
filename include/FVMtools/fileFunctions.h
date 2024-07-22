/** ****************************************************************************
 * @file      fileFunctions.h
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
 * @defgroup  fileFunctions fileFunctions
 * @ingroup   extract
 ******************************************************************************/
#ifndef FILEFUNCTIONS_H
#define FILEFUNCTIONS_H

#include <fstream>
#include <string>
#include <sys/stat.h>


using std::string;
using std::ifstream;



/** ****************************************************************************
 * @brief Check if the topo file exists and return its filename
 * @param output_dir Location of the output directory.
 * @param casename   Name of the case being extracted.
 * @param iteration  Iteration number being extracted
 * @return string    - Topo filename
 ******************************************************************************/
inline string getTopoFileName(string output_dir, string casename, string iteration)
{
  string      gridtopo     = output_dir+"/" + casename +".topo";
  string      toponamefile = output_dir + "/topo_file." + iteration + "_" + casename;
  struct stat tmpstat;
  if(stat(toponamefile.c_str(), &tmpstat) == 0)
  {
    ifstream tinput(toponamefile.c_str());
    string   name;
    tinput >> name;
    name = output_dir + "/" + name;
    if(stat(name.c_str(), &tmpstat) ==0)
      gridtopo = name;
  } // End If(toponamefile)

  return gridtopo;
} // End of getTopoFileName()

/** ****************************************************************************
 * @brief Check if the 'grid_pos' file exists and return its filename
 * @param output_dir Location of the output directory.
 * @param iteration  Iteration number being extracted
 * @param casename   Name of the case being extracted.
 * @return string    - Grid Position Filename
 ******************************************************************************/
inline string getPosFile(string output_dir, string iteration, string casename)
{
  string      posname = output_dir+"/grid_pos." + iteration + "_" + casename;
  struct stat tmpstat;

  if(stat(posname.c_str(),&tmpstat) != 0)
  {
    posname = output_dir+"/grid_pos." + casename;
  }else if(tmpstat.st_size == 0)
  {
    posname = output_dir+"/grid_pos." + casename;
  }

  return posname;
} // End of getPosFile()

/// @}

#endif
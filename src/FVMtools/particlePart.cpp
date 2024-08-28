/** ****************************************************************************
 * @file      particlePart.cc
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
 * @copyright Copyright (c) 2008-2019, Mississippi State University
 ******************************************************************************/

#include <sys/stat.h>
#include "hdf5_Functions.h"
#include "particlePart.h"

using std::cerr;
using std::endl;
using std::cout;

/** ****************************************************************************
 * @brief
 * @param output_dir
 * @param iteration
 * @param casename
 * @param vars
 * @param maxparticles
 ******************************************************************************/
particlePart::particlePart(string output_dir, string iteration, string casename,
                           vector<string> vars,
                           int maxparticles) {
  error = true;
  partName = "Particles";
  directory = output_dir;
  posfile = output_dir + "/particle_pos."+iteration + "_" + casename;

  struct stat tmpstat;
  if(stat(posfile.c_str(),&tmpstat) != 0) {
    return;
  }
  hid_t file_id = Loci::hdf5OpenFile(posfile.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);
  if(file_id < 0)
    return;
  numParticles = sizeElementType(file_id, "particle position");
  H5Fclose(file_id);
  stride_size = 1;
  if(maxparticles > 0) {
    stride_size = numParticles/maxparticles;
    stride_size = max(stride_size,1);
    int num_blocks = numParticles/stride_size;
    numParticles = num_blocks; //*stride_size;
  }
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i];
    string scalarfile = output_dir+"/"+varname+"_ptsca."
      +iteration+"_"+casename;
    if(stat(scalarfile.c_str(),&tmpstat) == 0) {
      scalarVars[varname] = scalarfile;
    } else {
      string vectorfile = output_dir+"/"+varname+"_ptvec."
        +iteration+"_"+casename;
      if(stat(vectorfile.c_str(),&tmpstat) == 0)
        vectorVars[varname] = vectorfile;
    }
  }
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool particlePart::hasScalarVar(string var) const {
  map<string,string>::const_iterator mi;
  mi = scalarVars.find(var);
  return !(mi == scalarVars.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool particlePart::hasVectorVar(string var) const {
  map<string,string>::const_iterator mi;
  mi = vectorVars.find(var);
  return !(mi == vectorVars.end());
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> particlePart::getScalarVars() const {
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=scalarVars.begin();mi!=scalarVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> particlePart::getVectorVars() const {
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=vectorVars.begin();mi!=vectorVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}

/** ****************************************************************************
 * @brief
 * @param ppos
 ******************************************************************************/
void particlePart::getParticlePositions(vector<vector3d<float> > &ppos) const {
  hid_t file_id = Loci::hdf5OpenFile(posfile.c_str(),
                                     H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    cerr << "unable to open file '" << posfile << "'!" << endl;
    return;
  }
  size_t np = sizeElementType(file_id, "particle position");
  vector<vector3d<float> > tmp(np);
  readElementType(file_id, "particle position", tmp);
  Loci::hdf5CloseFile(file_id);

  if(stride_size == 1)
    ppos.swap(tmp);
  else {
    vector<vector3d<float> > cpy(numParticles);
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = tmp[i*stride_size];
    ppos.swap(cpy);
  }
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param val
 ******************************************************************************/
void particlePart::getParticleScalar(string varname, vector<float> &val) const {
  map<string,string>::const_iterator mi;
  mi = scalarVars.find(varname);
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl;
  }
  size_t np = sizeElementType(file_id, varname.c_str());
  vector<float> scalar(np);
  readElementType(file_id, varname.c_str(), scalar);
  Loci::hdf5CloseFile(file_id);

  if(stride_size == 1)
    val.swap(scalar);
  else {
    vector<float > cpy(numParticles);
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = scalar[i*stride_size];
    val.swap(cpy);
  }
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param val
 ******************************************************************************/
void particlePart::getParticleVector(string varname,
                                     vector<vector3d<float> > &val) const {
  map<string,string>::const_iterator mi;
  mi = vectorVars.find(varname);
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl;
    return;
  }
  size_t np = sizeElementType(file_id, varname.c_str());
  vector<vector3d<float> > tmp(np);
  readElementType(file_id, varname.c_str(), tmp);
  Loci::hdf5CloseFile(file_id);

  if(stride_size == 1)
    val.swap(tmp);
  else {
    vector<vector3d<float> > cpy(numParticles);
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = tmp[i*stride_size];
    val.swap(cpy);
  }
}


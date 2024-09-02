/** ****************************************************************************
 * @file      extract.h
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
 * @defgroup  FVMtools FVMtools
 * @defgroup  extract extract
 * @ingroup   FVMtools
 ******************************************************************************/
#ifndef EXTRACT_H
#define EXTRACT_H

#include "particlePart.h"
#include "surfacePart.h"
#include "volumePart.h"

using std::list;

/// @addtogroup extract
/// @{

/// @brief 
enum var_type
{
  NODAL_SCALAR,
  NODAL_VECTOR,
  NODAL_DERIVED,
  NODAL_MASSFRACTION,
  BOUNDARY_SCALAR,
  BOUNDARY_VECTOR,
  BOUNDARY_DERIVED_SCALAR,
  BOUNDARY_DERIVED_VECTOR,
  PARTICLE_SCALAR,
  PARTICLE_VECTOR,
  UNDEFINED
};

/// @brief 
enum ePlot_type
{
  ASCII,
  TWODGV,
  ENSIGHT,
  CGNS,
  FIELDVIEW,
  TECPLOT,
  VTK,
  VTK_SURFACE,
  VTK64,
  VTK_SURFACE64,
  CUTTINGPLANE,
  SURFACE,
  MEAN,
  COMBINE,
  FCOMBINE,
  NONE
};

/// @brief 
enum view_type
{
  VIEWXY=0,
  VIEWYZ=1,
  VIEWXZ=2,
  VIEWXR=3
};


// convert a string to an integer
inline int str2int(string s) {
  std::stringstream ss;
  ss << s;
  int ret;
  ss >> ret;
  return ret;
}

// determine whether range [b,e) contains only digits
// [b,e) must be a valid range and contains characters
template <class ForwardIter> inline bool
alldigits(ForwardIter b, ForwardIter e) {
  for(;b!=e;++b)
    if(!isdigit(*b))
      return false;
  return true;
}
// determine whether s contains a valid integer (including the sign)
inline bool
valid_int(const std::string& s) {
  if(s.empty())
    return false;
  if( (s[0] == '-') || (s[0] == '+')) {
    if(s.size() > 1)
      return alldigits(s.begin()+1,s.end());
    else
      return false;
  }else
    return alldigits(s.begin(),s.end());
}



struct affineMapping {
  float M[4][4];
  affineMapping();
  void Combine(affineMapping a);
  void translate(vector3d<float> tv);
  void rotateX(float theta);
  void rotateY(float theta);
  void rotateZ(float theta);
  vector3d<float> MapNode(vector3d<float> v);
};

void get_2dgv(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries,
              int view);

void get_surf(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries);

void process_ascii_nodal(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames);

void process_ascii_bndry(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames,
                         vector<string> boundaries);

void process_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter);

void combine_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter,
                  bool do_favre);


template<class T> void writeElementType(hid_t group_id,
                                        const char *element_name,
                                        std::vector<T> &v) {
  hsize_t array_size = v.size();
  if(array_size == 0)
    return;
  int rank = 1;
  hsize_t dimension = array_size;

  hid_t dataspace = H5Screate_simple(rank,&dimension,NULL);

  typedef data_schema_traits<T> traits_type;
  Loci::DatatypeP dp = traits_type::get_type();


  hsize_t start = 0;
  hsize_t stride = 1;
  hsize_t count = v.size();
  hid_t datatype = dp->get_hdf5_type();
  hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                            dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  if(count != 0) {
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                        &start, &stride, &count, NULL);
    hid_t memspace = H5Screate_simple(rank, &count, NULL);
    H5Dwrite(dataset,datatype,memspace,dataspace,
             H5P_DEFAULT, &v[0]);
    H5Sclose(memspace);
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(datatype);
}

  
void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration);

namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]);
  }
}

extern string output_dir;


class postProcessorConvert : public Loci::CPTR_type {
protected:
  vector<surfacePartP> surfacePartList;
  vector<volumePartP> volumePartList;
  vector<particlePartP> particlePartList;
public:
  void addSurfaceParts(const vector<surfacePartP> &list) {
    for(size_t i=0;i<list.size();++i)
      surfacePartList.push_back(list[i]);
  }
  void addVolumePart(volumePartP volpart) {
    volumePartList.push_back(volpart);
  }
  void addParticlePart(particlePartP particles) {
    particlePartList.push_back(particles);
  }
  virtual bool processesVolumeElements() const = 0;
  virtual bool processesSurfaceElements() const = 0;
  virtual bool processesParticleElements() const = 0;

  virtual void exportPostProcessorFiles(string casename, string iteration) const = 0;
};

typedef Loci::CPTR<postProcessorConvert> postProcessorP;

class ensightPartConverter : public postProcessorConvert {
  bool id_required;
public:
  ensightPartConverter(bool input) {id_required = input; };
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};

//#ifdef HAVE_CGNS
class cgnsPartConverter : public postProcessorConvert {
  bool id_required;
public:
  cgnsPartConverter(bool input) {id_required = input; };
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};
//#endif

class tecplotPartConverter : public postProcessorConvert {
public:
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};

class vtkPartConverter : public postProcessorConvert {
  bool bit64;
public:
  vtkPartConverter() {bit64 = false; }
  vtkPartConverter(bool input) {bit64 = input; };
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};

class vtkSurfacePartConverter : public postProcessorConvert {
  bool bit64;
public:
  vtkSurfacePartConverter() {bit64=false; }
  vtkSurfacePartConverter(bool input) { bit64=input; }
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};

class fieldViewPartConverter : public postProcessorConvert {
public:
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};

class cuttingPlanePartConverter : public postProcessorConvert
{
  affineMapping transformMatrix;
  float xShift, yShift, zShift;
public:
  cuttingPlanePartConverter(const affineMapping &m,
                            float xs,float ys,float zs) {
    transformMatrix = m;
    xShift=xs;
    yShift=ys;
    zShift=zs;
  }
  virtual bool processesVolumeElements() const;
  virtual bool processesSurfaceElements() const;
  virtual bool processesParticleElements() const;

  virtual void exportPostProcessorFiles(string casename, string iteration) const;
};


/// @}

#endif

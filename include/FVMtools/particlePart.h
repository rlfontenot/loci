/** ****************************************************************************
 * @file      particlePart.h
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
 * @defgroup  particlePart particlePart
 * @ingroup   extract
 ******************************************************************************/
#ifndef PARTICLEPART_H
#define PARTICLEPART_H

#include "Loci.h"

using std::string;
using std::vector;
using std::map;


/// @addtogroup particlePart
/// @{


class particlePartBase : public Loci::CPTR_type {
protected:
  bool error;
  string partName;
  size_t numParticles;
public:
  bool fail() const { return error; }
  string getPartName() const { return partName; }
  size_t getNumParticles() const { return numParticles; }
  virtual bool hasScalarVar(string var) const = 0;
  virtual bool hasVectorVar(string var) const = 0;
  virtual vector<string> getScalarVars() const = 0;
  virtual vector<string> getVectorVars() const = 0;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const = 0;
  virtual void getParticleScalar(string varname, vector<float> &val) const = 0;
  virtual void getParticleVector(string varname,
                                 vector<vector3d<float> > &val) const  = 0;
};

typedef Loci::CPTR<particlePartBase> particlePartP;

class particlePart: public particlePartBase {
  string directory;
  string posfile;
  map<string,string> scalarVars;
  map<string,string> vectorVars;
  int stride_size;
public:
  particlePart() { error = true; }
  particlePart(string output_dir, string iteration, string casename,
               vector<string> vars, int maxparticles);
  virtual bool hasScalarVar(string var) const;
  virtual bool hasVectorVar(string var) const;
  virtual vector<string> getScalarVars() const;
  virtual vector<string> getVectorVars() const;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const;
  virtual void getParticleScalar(string varname, vector<float> &val) const;
  virtual void getParticleVector(string varname,
                                 vector<vector3d<float> > &val) const;
};

/// @}

#endif

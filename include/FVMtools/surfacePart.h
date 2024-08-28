/** ****************************************************************************
 * @file      surfacePart.h
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
 * @defgroup  surfacePart surfacePart
 * @ingroup   extract
 ******************************************************************************/
#ifndef SURFACEPART_H
#define SURFACEPART_H

#include "Loci.h"

using std::string;
using std::vector;
using std::map;


/// @addtogroup surfacePart
/// @{

class surfacePartBase : public Loci::CPTR_type {
protected:
  bool error;
  string partName;
  int nnodes, nquads, ntrias, ngenf;
public:
  bool fail() const { return error; }
  string getPartName() const { return partName; }
  int getNumNodes() const { return nnodes; }
  int getNumQuads() const { return nquads; }
  int getNumTrias() const { return ntrias; }
  int getNumGenfc() const { return ngenf; }
  virtual bool hasNodalScalarVar(string var) const = 0;
  virtual bool hasNodalVectorVar(string var) const = 0;
  virtual bool hasElementScalarVar(string var) const = 0;
  virtual bool hasElementVectorVar(string var) const = 0;
  virtual vector<string> getNodalScalarVars() const = 0;
  virtual vector<string> getNodalVectorVars() const = 0;
  virtual vector<string> getElementScalarVars() const = 0;
  virtual vector<string> getElementVectorVars() const = 0;
  virtual void getQuads(vector<Array<int,4> > &quads) const = 0;
  virtual void getTrias(vector<Array<int,3> > &trias) const = 0;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const =0;
  virtual void getQuadsIds(vector<int> &quads_ids) const = 0;
  virtual void getTriasIds(vector<int> &trias_ids) const = 0;
  virtual void getGenfIds(vector<int> &genface_ids) const =0;
  
  virtual void getPos(vector<vector3d<float> > &pos) const = 0;
  virtual void getPos(vector<vector3d<double> > &pos) const = 0;
  virtual void getNodalScalar(string varname, vector<float> &vals) const = 0;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0;
  virtual void getElementScalar(string varname, vector<float> &qvals,
                                vector<float> &tvals, 
                                vector<float> &gvals) const =0;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
                                vector<vector3d<float> > &tvals,
                                vector<vector3d<float> > &gvals) const = 0;
};
typedef Loci::CPTR<surfacePartBase> surfacePartP;

class surfacePartDerivedVars : public surfacePartBase {
  surfacePartP shadowPart;
  enum derivedVar_t { VAR_M, VAR_P, VAR_logp, VAR_U, VAR_0, VAR_1, VAR_2,
                      VAR_X, VAR_Y, VAR_Z };
  map<string,derivedVar_t> derivedVars;
  float Pambient;
  void processDerivedVars(const vector<string> &vars);
public:
  surfacePartDerivedVars() { error = true; shadowPart=0;} 
  surfacePartDerivedVars(surfacePartP part, string output_dir,
                         string casename ,
                         string iteration, 
                         vector<string> vars);
  virtual bool hasNodalScalarVar(string var) const;
  virtual bool hasNodalVectorVar(string var) const;
  virtual bool hasElementScalarVar(string var) const;
  virtual bool hasElementVectorVar(string var) const;
  virtual vector<string> getNodalScalarVars() const;
  virtual vector<string> getNodalVectorVars() const;
  virtual vector<string> getElementScalarVars() const;
  virtual vector<string> getElementVectorVars() const;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const;
  virtual void getTrias(vector<Array<int,3> > &trias) const;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const;
  virtual void getQuadsIds(vector<int> &quads_ids) const;
  virtual void getTriasIds(vector<int> &trias_ids) const;
  virtual void getGenfIds(vector<int> &genface_ids) const;
  virtual void getPos(vector<vector3d<float> > &pos) const;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getNodalScalar(string varname, vector<float> &vals) const;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const;
  virtual void getElementScalar(string varname, vector<float> &qvals,
                                vector<float> &tvals, 
                                vector<float> &gvals) const;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
                                vector<vector3d<float> > &tvals,
                                vector<vector3d<float> > &gvals) const;
};
    
      
class surfacePart : public surfacePartBase {
  string directory;
  string topoFile;
  string posFile;
  // maps from variables to file name
  map<string,string> nodalScalarVars;
  map<string,string> nodalVectorVars;
  map<string,string> elementScalarVars;
  map<string,string> elementVectorVars;

  vector<int> quad_ord;
  vector<int> tri_ord;
  vector<int> gen_ord;
  entitySet quadSet,triSet,genSet;
public:
  surfacePart() {error = true;}
  surfacePart(string name, string directory, string iteration,
              vector<string> vars);
  virtual bool hasNodalScalarVar(string var) const;
  virtual bool hasNodalVectorVar(string var) const;
  virtual bool hasElementScalarVar(string var) const;
  virtual bool hasElementVectorVar(string var) const;
  virtual vector<string> getNodalScalarVars() const;
  virtual vector<string> getNodalVectorVars() const;
  virtual vector<string> getElementScalarVars() const;
  virtual vector<string> getElementVectorVars() const;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const;
  virtual void getTrias(vector<Array<int,3> > &trias) const;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const;
  virtual void getQuadsIds(vector<int> &quads_ids) const;
  virtual void getTriasIds(vector<int> &trias_ids) const;
  virtual void getGenfIds(vector<int> &genface_ids) const;
  
  virtual void getPos(vector<vector3d<float> > &pos) const;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getNodalScalar(string varname, vector<float> &vals) const;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const;
  virtual void getElementScalar(string varname, vector<float> &qvals,
                                vector<float> &tvals, 
                                vector<float> &gvals) const;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
                                vector<vector3d<float> > &tvals,
                                vector<vector3d<float> > &gvals) const;
};

class surfacePartCopy : public surfacePartBase {
  vector<Array<int,3> > trifaces;
  vector<Array<int,4> > quadfaces;
  vector<int> nfacenodes, gennodes;
  vector<int> triaIds;
  vector<int> quadIds;
  vector<int> genIds;
  vector<int> nodemap;
  vector<vector3d<double> > pos;
  map<string,vector<float> > nodalScalars;
  map<string,vector<vector3d<float> > > nodalVectors;
  map<string,Array<vector<float>,3> > elementScalars;
  map<string,Array<vector<vector3d<float> >,3> > elementVectors;
public:
  surfacePartCopy() {error = true;}
  surfacePartCopy(string name, vector<Array<int,3> > &triangles,
                  vector<int> &tria_ids,
                  vector<Array<int,4> > &quads,
                  vector<int>& quads_ids,
                  vector<int> &genface2n, vector<int> &gnodes,
                  vector<int>&gen_ids);
  void registerPos(const vector<vector3d<double> > &pos);
  void registerNodalScalar(string name,const vector<float> &val);
  void registerNodalVector(string name,const vector<vector3d<float> > &val);
  void registerElementScalar(string name, 
                             const vector<float> &qval,
                             const vector<float> &tval,
                             const vector<float> &gval);
  void registerElementVector(string name, 
                             const vector<vector3d<float> > &qval,
                             const vector<vector3d<float> > &tval,
                             const vector<vector3d<float> > &gval);
  virtual bool hasNodalScalarVar(string var) const;
  virtual bool hasNodalVectorVar(string var) const;
  virtual bool hasElementScalarVar(string var) const;
  virtual bool hasElementVectorVar(string var) const;
  virtual vector<string> getNodalScalarVars() const;
  virtual vector<string> getNodalVectorVars() const;
  virtual vector<string> getElementScalarVars() const;
  virtual vector<string> getElementVectorVars() const;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const;
  virtual void getTrias(vector<Array<int,3> > &trias) const;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const;
  virtual void getQuadsIds(vector<int> &quads_ids) const;  
  virtual void getTriasIds(vector<int> &trias_ids) const; 
  virtual void getGenfIds(vector<int> &genface_ids) const; 
  
  virtual void getPos(vector<vector3d<float> > &pos) const;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getNodalScalar(string varname, vector<float> &vals) const;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const;
  virtual void getElementScalar(string varname, vector<float> &qvals,
                                vector<float> &tvals, 
                                vector<float> &gvals) const;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
                                vector<vector3d<float> > &tvals,
                                vector<vector3d<float> > &gvals) const;
};


/// @}

#endif

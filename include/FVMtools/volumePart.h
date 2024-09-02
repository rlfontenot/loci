/** ****************************************************************************
 * @file      volumePart.h
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
 * @defgroup  volumePart volumePart
 * @ingroup   extract
 ******************************************************************************/
#ifndef VOLUMEPART_H
#define VOLUMEPART_H

using std::string;
using std::vector;
using std::map;


/// @addtogroup volumePart
/// @{

/** ****************************************************************************
 * @brief Base class for handling volume part of the extract tool.
 ******************************************************************************/
class volumePartBase : public Loci::CPTR_type
{
protected:
  bool   error;       //!<
  string partName;    //!<
  size_t nnodes;      //!<
  size_t ntets;       //!<
  size_t nhexs;       //!<
  size_t nprsm;       //!<
  size_t npyrm;       //!<
  size_t ngenc;       //!<
  size_t ntetsIblank; //!<
  size_t nhexsIblank; //!<
  size_t nprsmIblank; //!<
  size_t npyrmIblank; //!<
  size_t ngencIblank; //!<

public:
  bool fail() const { return error; }
  string getPartName() const { return partName; }
  size_t getNumNodes() const { return nnodes; }

  size_t getNumTets() const { return ntets; }
  size_t getNumHexs() const { return nhexs; }
  size_t getNumPrsm() const { return nprsm; }
  size_t getNumPyrm() const { return npyrm; }
  size_t getNumGenc() const { return ngenc; }
  size_t getNumTetsIblank() const { return ntetsIblank; }
  size_t getNumHexsIblank() const { return nhexsIblank; }
  size_t getNumPrsmIblank() const { return nprsmIblank; }
  size_t getNumPyrmIblank() const { return npyrmIblank; }
  size_t getNumGencIblank() const { return ngencIblank; }

  virtual bool hasNodalScalarVar(string var) const = 0;
  virtual bool hasNodalVectorVar(string var) const = 0;

  virtual void getPos(vector<vector3d<float> >  &pos) const = 0;
  virtual void getPos(vector<vector3d<double> > &pos) const = 0;

  virtual void getTetBlock( vector<Array<int,4> > &tets,  size_t start, size_t size) const = 0;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const = 0;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const = 0;
  virtual void getHexBlock( vector<Array<int,8> > &hexs,  size_t start, size_t size) const = 0;

  virtual void getTetIds( vector<int> &tetids,  size_t start, size_t size) const = 0;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const = 0;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const = 0;
  virtual void getHexIds( vector<int> &hexids,  size_t start, size_t size) const = 0;
  virtual void getGenCell(vector<int> &genCellNfaces, vector<int> &genCellNsides,
                          vector<int> &genCellNodes) const = 0;
  virtual void getGenIds( vector<int> &genids)       const = 0;

  virtual void getNodalScalar(string varname, vector<float>            &vals) const = 0;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0;
  virtual void getNodalIblank(vector<unsigned char> &blank) const = 0;

  virtual vector<string> getNodalScalarVars() const = 0;
  virtual vector<string> getNodalVectorVars() const = 0;
}; // End Class volumePartBase



/// @brief Loci pointer for volumePartBase class
typedef Loci::CPTR<volumePartBase> volumePartP;



/** ****************************************************************************
 * @brief 
 ******************************************************************************/
class volumePart : public volumePartBase
{
  bool   has_iblank; //!<
  string directory;  //!<
  string topoFile;   //!<
  string posFile;    //!<

  // maps from variables to file name
  map<string,string>   nodalScalarVars; //!<
  map<string,string>   nodalVectorVars; //!<
  store<unsigned char> iblank;          //!<

  size_t ntets_orig;  //!<
  size_t nhexs_orig;  //!<
  size_t nprsm_orig;  //!<
  size_t npyrm_orig;  //!<
  size_t ngenc_orig;  //!<

  Loci::entitySet tetsIblanked; //!<
  Loci::entitySet hexsIblanked; //!<
  Loci::entitySet prsmIblanked; //!<
  Loci::entitySet pyrmIblanked; //!<
  Loci::entitySet gencIblanked; //!<
public:
  volumePart() {error = true;}
  volumePart(string output_dir, string iteration, string casename, vector<string> vars);

  virtual bool hasNodalScalarVar(string var) const;
  virtual bool hasNodalVectorVar(string var) const;

  virtual void getPos(vector<vector3d<float> >  &val) const;
  virtual void getPos(vector<vector3d<double> > &pos) const;

  virtual void getTetBlock( vector<Array<int,4> > &tets,  size_t start, size_t size) const;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const;
  virtual void getHexBlock( vector<Array<int,8> > &hexs,  size_t start, size_t size) const;

  virtual void getTetIds( vector<int> &tetids,  size_t start, size_t size) const;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const;
  virtual void getHexIds( vector<int> &hexids,  size_t start, size_t size) const;
  virtual void getGenCell(vector<int> &genCellNfaces, vector<int> &genCellNsides,
                          vector<int> &genCellNodes) const;
  virtual void getGenIds( vector<int> &genids)       const;

  virtual void getNodalScalar(string varname, vector<float>            &vals) const;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const;
  virtual void getNodalIblank(vector<unsigned char> &blank) const;

  virtual vector<string> getNodalScalarVars() const;
  virtual vector<string> getNodalVectorVars() const;
}; // End Class volumePart


/** ****************************************************************************
 * @brief 
 ******************************************************************************/
class volumePartDerivedVars : public volumePartBase
{
  enum derivedVar_t
  {
    VAR_M,    //!<
    VAR_P,    //!<
    VAR_logp, //!<
    VAR_U,    //!<
    VAR_0,    //!<
    VAR_1,    //!<
    VAR_2,    //!<
    VAR_X,    //!<
    VAR_Y,    //!<
    VAR_Z     //!<
  };
  float                    Pambient;    //!<
  volumePartP              shadowPart;  //!<
  map<string,derivedVar_t> derivedVars; //!<

  void processDerivedVars(const vector<string> &vars);
public:
  volumePartDerivedVars() {error = true;}
  volumePartDerivedVars(volumePartP part, string output_dir, string iteration,
                        string casename,  vector<string> vars);
  virtual bool hasNodalScalarVar(string var) const;
  virtual bool hasNodalVectorVar(string var) const;

  virtual void getPos(vector<vector3d<float> >  &val) const;
  virtual void getPos(vector<vector3d<double> > &pos) const;

  virtual void getTetBlock( vector<Array<int,4> > &tets,  size_t start, size_t size) const;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const;
  virtual void getHexBlock( vector<Array<int,8> > &hexs,  size_t start, size_t size) const;

  virtual void getTetIds( vector<int> &tetids,  size_t start, size_t size) const;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const;
  virtual void getHexIds( vector<int> &hexids,  size_t start, size_t size) const;
  virtual void getGenCell(vector<int> &genCellNfaces, vector<int> &genCellNsides,
                          vector<int> &genCellNodes) const;
  virtual void getGenIds( vector<int> &genids)       const;

  virtual void getNodalScalar(string varname, vector<float>            &vals) const;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const;
  virtual void getNodalIblank(vector<unsigned char> &blank) const;

  virtual vector<string> getNodalScalarVars() const;
  virtual vector<string> getNodalVectorVars() const;
}; // End Class volumePartDerivedVars

/// @}

#endif

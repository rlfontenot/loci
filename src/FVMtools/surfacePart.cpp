/** ****************************************************************************
 * @file      surfacePart.cc
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
#include "surfacePart.h"

using std::cerr;
using std::endl;
using std::cout;
using std::ios;
using std::ifstream;


/** ****************************************************************************
 * @brief
 * @param name
 * @param dir
 * @param iteration
 * @param vars
 ******************************************************************************/
surfacePart::surfacePart(string name, string dir, string iteration,
                         vector<string> vars)
{
  partName = name;
  directory = dir;
  nnodes = 0;
  nquads = 0;
  ntrias = 0;
  ngenf = 0;
  error = true;
  string topo_link = dir + "/topo_file."+iteration;
  ifstream topo_links(topo_link.c_str(),ios::in);
  //  if(topo_links.fail()) cerr << "topo_links fail, " << topo_link << endl;
  if(topo_links.fail()) return;

  string topolink;
  topo_links>> topolink;
  topo_links.close();

  posFile = dir + "/pos." + iteration;
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);
  if(file_id < 0) cerr << posFile << " fail" << endl;
  if(file_id < 0) return;

  nnodes = sizeElementType(file_id,"data");
  Loci::hdf5CloseFile(file_id);

  topoFile = dir + "/" + topolink;
  file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                               H5F_ACC_RDONLY,
                               H5P_DEFAULT);
  if(file_id < 0) cerr << topoFile << " fail" << endl;
  if(file_id < 0) return;

  ngenf = sizeElementType(file_id,"nside_sizes");
  nquads = sizeElementType(file_id,"quads");
  ntrias = sizeElementType(file_id,"triangles");

  if(nquads > 0)
    quadSet = interval(0,nquads-1);
  if(ntrias > 0)
    triSet = interval(0,ntrias-1);
  if(ngenf > 0)
    genSet = interval(0,ngenf-1);

  Loci::hdf5CloseFile(file_id);
  bool has_element_data = false;
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i];
    // try scalar
    string svar = dir+"/" + varname+"_sca."+iteration;
    file_id = Loci::hdf5OpenFile(svar.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    bool found_var = false;
    if(file_id >= 0) {
      int nsz = sizeElementType(file_id,"data");
      if(nsz == nnodes) {
        nodalScalarVars[varname] = svar;
        found_var = true;
      }
      Loci::hdf5CloseFile(file_id);
    }

    if(!found_var) {
      svar = dir+"/" + varname+"_vec."+iteration;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      if(file_id >= 0) {
        int nsz = sizeElementType(file_id,"data");
        if(nsz == nnodes) {
          nodalVectorVars[varname] = svar;
          found_var = true;
        }
        Loci::hdf5CloseFile(file_id);
      }
    }
    if(!found_var) {
      svar = dir+"/" + varname+"_bsca."+iteration;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      if(file_id >= 0) {
        int nsz = sizeElementType(file_id,"data");
        if(nsz == (nquads+ntrias+ngenf)) {
          elementScalarVars[varname] = svar;
          found_var = true;
          has_element_data = true;
        }
        Loci::hdf5CloseFile(file_id);
      }
    }
    if(!found_var) {
      svar = dir+"/" + varname+"_bvec."+iteration;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
                                   H5F_ACC_RDONLY,
                                   H5P_DEFAULT);
      if(file_id >= 0) {
        int nsz = sizeElementType(file_id,"data");
        if(nsz == (nquads+ntrias+ngenf)) {
          elementVectorVars[varname] = svar;
          found_var = true;
          has_element_data = true;
        }
        Loci::hdf5CloseFile(file_id);
      }
    }
  }
  vector<unsigned char> iblank;
  string iblank_file = dir +"/iblank."+iteration;
  file_id = Loci::hdf5OpenFile(iblank_file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id >=0) {
    vector<unsigned char> tmp(nquads+ntrias+ngenf);
    iblank.swap(tmp);
    readElementType(file_id,"data",iblank);
    Loci::hdf5CloseFile(file_id);
    for(size_t i=0;i<iblank.size();++i)
      has_element_data = (iblank[i] > 1) || has_element_data;
  }
  if(has_element_data) {
    file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    if(nquads > 0) {
      vector<int> tmp(nquads);
      readElementType(file_id,"quads_ord",tmp);
      quad_ord.swap(tmp);
    }

    if(ntrias > 0) {
      vector<int> tmp(ntrias);
      readElementType(file_id,"triangles_ord",tmp);
      tri_ord.swap(tmp);
    }

    if(ngenf > 0) {
      vector<int> tmp(ngenf);
      readElementType(file_id,"nside_ord",tmp);
      gen_ord.swap(tmp);
    }
    Loci::hdf5CloseFile(file_id);
    if(iblank.size() > 0) {
      // compute iblanked set
      quadSet = EMPTY;
      for(int i=0;i<nquads;++i)
        if(iblank[quad_ord[i]] < 2)
          quadSet += i;
      nquads = quadSet.size();
      triSet = EMPTY;
      for(int i=0;i<ntrias;++i)
        if(iblank[tri_ord[i]] < 2)
          triSet += i;
      ntrias = triSet.size();
      genSet = EMPTY;
      for(int i=0;i<ngenf;++i)
        if(iblank[gen_ord[i]] < 2)
          genSet += i;
      ngenf = genSet.size();
    }
  }

  error = false;
} // End of surfacePart()


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePart::hasNodalScalarVar(string var) const
{
  map<string,string>::const_iterator mi=nodalScalarVars.find(var);
  return (mi != nodalScalarVars.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePart::hasNodalVectorVar(string var) const
{
  map<string,string>::const_iterator mi=nodalVectorVars.find(var);
  return (mi != nodalVectorVars.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePart::hasElementScalarVar(string var) const
{
  map<string,string>::const_iterator mi=elementScalarVars.find(var);
  return (mi != elementScalarVars.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePart::hasElementVectorVar(string var) const
{
  map<string,string>::const_iterator mi=elementVectorVars.find(var);
  return (mi != elementVectorVars.end());
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePart::getNodalScalarVars() const
{
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=nodalScalarVars.begin();mi!=nodalScalarVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePart::getNodalVectorVars() const
{
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=nodalVectorVars.begin();mi!=nodalVectorVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePart::getElementScalarVars() const
{
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=elementScalarVars.begin();mi!=elementScalarVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePart::getElementVectorVars() const {
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=elementVectorVars.begin();mi!=elementVectorVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @param quads
 ******************************************************************************/
void surfacePart::getQuads(vector<Array<int,4> > &quads) const {
  quads.clear();
  if(nquads > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) return;
    int nq = sizeElementType(file_id,"quads");

    vector<Array<int,4> > tmp(nq);
    readElementType(file_id,"quads",tmp);
    quads.resize(nquads);
    int cnt = 0;
    FORALL(quadSet,ii) {
      quads[cnt] = tmp[ii];
      cnt++;
    } ENDFORALL;
    Loci::hdf5CloseFile(file_id);
  }
}


/** ****************************************************************************
 * @brief
 * @param quads_ids
 ******************************************************************************/
void surfacePart::getQuadsIds(vector<int> &quads_ids) const {
  cout << "start getQuadsIds " << endl;
  quads_ids.clear();
  FORALL(quadSet,ii) {
    // quads_ids.push_back(quad_ord[ii]);
    quads_ids.push_back(ii);
  } ENDFORALL;
  cout << "start getQuadsIds " << endl;
}


/** ****************************************************************************
 * @brief
 * @param trias
 ******************************************************************************/
void surfacePart::getTrias(vector<Array<int,3> > &trias) const {
  trias.clear();
  if(ntrias > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) return;
    int nt = sizeElementType(file_id,"triangles");
    vector<Array<int,3> > tmp(nt);
    readElementType(file_id,"triangles",tmp);
    trias.resize(ntrias);
    int cnt = 0;
    FORALL(triSet,ii) {
      trias[cnt] = tmp[ii];
      cnt++;
    } ENDFORALL;
    Loci::hdf5CloseFile(file_id);
  }
}


/** ****************************************************************************
 * @brief
 * @param trias_ids
 ******************************************************************************/
void  surfacePart::getTriasIds(vector<int> &trias_ids) const{
  cout << "start getTriasIds " << endl;
  trias_ids.clear();
  FORALL(triSet,ii) {
    //  trias_ids.push_back(tri_ord[ii]);
    trias_ids.push_back(ii);
  } ENDFORALL;
  cout << "end getTriasIds " << endl;
}


/** ****************************************************************************
 * @brief
 * @param numGenFnodes
 * @param genNodes
 ******************************************************************************/
void surfacePart::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  numGenFnodes.clear();
  genNodes.clear();
  if(ngenf > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) return;

    int ng = sizeElementType(file_id,"nside_sizes");
    vector<int> tmp(ng);
    readElementType(file_id,"nside_sizes",tmp);

    int nsz = sizeElementType(file_id,"nside_nodes");
    vector<int> tmp2(nsz);
    readElementType(file_id,"nside_nodes",tmp2);

    vector<int> sum(ng);
    sum[0] = 0;
    for(int i=1;i<ng;++i)
      sum[i] = sum[i-1]+tmp[i-1];
    FORALL(genSet,ii) {
      numGenFnodes.push_back(tmp[ii]);
      for(int i=0;i<tmp[ii];++i)
        genNodes.push_back(tmp2[sum[ii]+i]);
    } ENDFORALL;
    Loci::hdf5CloseFile(file_id);
  }
}


/** ****************************************************************************
 * @brief
 * @param genface_ids
 ******************************************************************************/
void  surfacePart::getGenfIds(vector<int> &genface_ids) const{
  cout << "start getGenIds " << endl;
  genface_ids.clear();
  FORALL(genSet,ii) {
    //genface_ids.push_back(gen_ord[ii]);
    genface_ids.push_back(ii);
  } ENDFORALL;
  cout << "end getGenIds " << endl;
}


/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void surfacePart::getPos(vector<vector3d<float> > &pos) const {
  vector<vector3d<float> > tmp(nnodes);
  pos.swap(tmp);
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  readElementType(file_id,"data",pos);
  Loci::hdf5CloseFile(file_id);
}


/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void surfacePart::getPos(vector<vector3d<double> > &pos) const {
  vector<vector3d<double> > tmp(nnodes);
  pos.swap(tmp);
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  readElementType(file_id,"data",pos);
  Loci::hdf5CloseFile(file_id);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePart::getNodalScalar(string varname,
                                 vector<float> &vals) const {
  vector<float> tmp(nnodes);
  vals.swap(tmp);
  map<string,string>::const_iterator mi = nodalScalarVars.find(varname);
  if(mi == nodalScalarVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  readElementType(file_id,"data",vals);
  Loci::hdf5CloseFile(file_id);

}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePart::getNodalVector(string varname,
                                 vector<vector3d<float> > &vals) const {
  vector<vector3d<float> > tmp(nnodes);
  vals.swap(tmp);
  map<string,string>::const_iterator mi = nodalVectorVars.find(varname);
  if(mi == nodalScalarVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  readElementType(file_id,"data",vals);
  Loci::hdf5CloseFile(file_id);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePart::getElementScalar(string varname,
                                   vector<float> &qvals,
                                   vector<float> &tvals,
                                   vector<float> &gvals) const {
  { vector<float> tmpq(nquads);  qvals.swap(tmpq); }
  { vector<float> tmpt(ntrias);  tvals.swap(tmpt); }
  { vector<float> tmpg(ngenf);  gvals.swap(tmpg); }

  map<string,string>::const_iterator mi = elementScalarVars.find(varname);
  if(mi == elementScalarVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  vector<float> vals(quad_ord.size()+tri_ord.size()+gen_ord.size());
  readElementType(file_id,"data",vals);
  Loci::hdf5CloseFile(file_id);
  int i=0;
  FORALL(quadSet,ii) {
    qvals[i] = vals[quad_ord[ii]];
    i++;
  } ENDFORALL;
  i=0;
  FORALL(triSet,ii) {
    tvals[i] = vals[tri_ord[ii]];
    i++;
  } ENDFORALL;
  i=0;
  FORALL(genSet,ii) {
    gvals[i] = vals[gen_ord[ii]];
    i++;
  } ENDFORALL;
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePart::getElementVector(string varname,
                                   vector<vector3d<float> > &qvals,
                                   vector<vector3d<float> > &tvals,
                                   vector<vector3d<float> > &gvals) const
{
  { vector<vector3d<float> > tmpq(nquads);  qvals.swap(tmpq); }
  { vector<vector3d<float> > tmpt(ntrias);  tvals.swap(tmpt); }
  { vector<vector3d<float> > tmpg(ngenf);  gvals.swap(tmpg); }

  map<string,string>::const_iterator mi = elementVectorVars.find(varname);
  if(mi == elementVectorVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);

  if(file_id < 0) return;
  vector<vector3d<float> > vals(quad_ord.size()+tri_ord.size()+gen_ord.size());
  readElementType(file_id,"data",vals);
  Loci::hdf5CloseFile(file_id);
  int i=0;
  FORALL(quadSet,ii) {
    qvals[i] = vals[quad_ord[ii]];
    i++;
  } ENDFORALL;
  i=0;
  FORALL(triSet,ii) {
    tvals[i] = vals[tri_ord[ii]];
    i++;
  } ENDFORALL;
  i=0;
  FORALL(genSet,ii) {
    gvals[i] = vals[gen_ord[ii]];
    i++;
  } ENDFORALL;
}


/** ****************************************************************************
 * @brief
 * @param vars
 ******************************************************************************/
void surfacePartDerivedVars::processDerivedVars(const vector<string> &vars)
{
  for(size_t i=0;i<vars.size();++i) {
    if(vars[i] == "m" && !shadowPart->hasNodalScalarVar("m")) {
      if(shadowPart->hasNodalScalarVar("a") &&
         shadowPart->hasNodalVectorVar("v"))
        derivedVars["m"] = VAR_M;
    }
    if(vars[i] == "P" && !shadowPart->hasNodalScalarVar("P")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
        derivedVars["P"] = VAR_P;
      }
    }
    if(vars[i] == "p" && !shadowPart->hasNodalScalarVar("p")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
        derivedVars["p"] = VAR_logp;
      }
    }
    if(vars[i] == "u" && !shadowPart->hasNodalScalarVar("u")) {
      if(shadowPart->hasNodalVectorVar("v")) {
        derivedVars["u"] = VAR_U;
      }
    }
    if(vars[i] == "0" && !shadowPart->hasNodalScalarVar("0")) {
      if(shadowPart->hasNodalVectorVar("v")) {
        derivedVars["0"] = VAR_0;
      }
    }
    if(vars[i] == "1" && !shadowPart->hasNodalScalarVar("1")) {
      if(shadowPart->hasNodalVectorVar("v")) {
        derivedVars["1"] = VAR_1;
      }
    }
    if(vars[i] == "2" && !shadowPart->hasNodalScalarVar("2")) {
      if(shadowPart->hasNodalVectorVar("v")) {
        derivedVars["2"] = VAR_2;
      }
    }
    if(vars[i] == "x")
      derivedVars["x"] = VAR_X;
    if(vars[i] == "y")
      derivedVars["y"] = VAR_Y;
    if(vars[i] == "z")
      derivedVars["z"] = VAR_Z;
  }
}


/** ****************************************************************************
 * @brief
 * @param part
 * @param output_dir
 * @param casename
 * @param iteration
 * @param vars
 ******************************************************************************/
surfacePartDerivedVars::surfacePartDerivedVars(surfacePartP part,
                                               string output_dir,
                                               string casename ,
                                               string iteration,
                                               vector<string> vars) {
  error = part->fail();
  partName = part->getPartName();
  nnodes = part->getNumNodes();
  nquads = part->getNumQuads();
  ntrias = part->getNumTrias();
  ngenf = part->getNumGenfc();
  shadowPart = part;

  string filename = output_dir+"/Pambient_par." + iteration +"_" + casename;
  struct stat tmpstat;
  if(stat(filename.c_str(),&tmpstat) == 0) {
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT);
    Pambient = 0;
    if(file_id >= 0) {
      fact_db facts;
      param<float> Pamb;
      readData(file_id,"Pambient",Pamb.Rep(),EMPTY,facts);
      Loci::hdf5CloseFile(file_id);
      Pambient = *Pamb;
    } else {
      cerr << "unable to open " << filename << endl;
    }
    processDerivedVars(vars);
  }
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartDerivedVars::hasNodalScalarVar(string var) const {
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(var);
  if(mi==derivedVars.end())
    return shadowPart->hasNodalScalarVar(var);
  else
    return true;
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartDerivedVars::hasNodalVectorVar(string var) const {
  return shadowPart->hasNodalVectorVar(var);
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartDerivedVars::hasElementScalarVar(string var) const {
  return shadowPart->hasElementScalarVar(var);
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartDerivedVars::hasElementVectorVar(string var) const {
  return shadowPart->hasElementVectorVar(var);
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartDerivedVars::getNodalScalarVars() const {
  vector<string> tmp = shadowPart->getNodalScalarVars();
  map<string,derivedVar_t>::const_iterator mi;
  for(mi=derivedVars.begin();mi!=derivedVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartDerivedVars::getNodalVectorVars() const {
  return shadowPart->getNodalVectorVars();
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartDerivedVars::getElementScalarVars() const {
  return shadowPart->getElementScalarVars();
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartDerivedVars::getElementVectorVars() const {
  return shadowPart->getElementVectorVars();
}


/** ****************************************************************************
 * @brief
 * @param quads
 ******************************************************************************/
void surfacePartDerivedVars::getQuads(vector<Array<int,4> > &quads) const {
  shadowPart->getQuads(quads);
}


/** ****************************************************************************
 * @brief
 * @param trias
 ******************************************************************************/
void surfacePartDerivedVars::getTrias(vector<Array<int,3> > &trias) const {
  shadowPart->getTrias(trias);
}


/** ****************************************************************************
 * @brief
 * @param numGenFnodes
 * @param genNodes
 ******************************************************************************/
void surfacePartDerivedVars::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  shadowPart->getGenf(numGenFnodes,genNodes);
}


/** ****************************************************************************
 * @brief
 * @param quads_ids
 ******************************************************************************/
void  surfacePartDerivedVars::getQuadsIds(vector<int> &quads_ids) const{
  shadowPart->getQuadsIds(quads_ids);
}


/** ****************************************************************************
 * @brief
 * @param trias_ids
 ******************************************************************************/
void  surfacePartDerivedVars::getTriasIds(vector<int> &trias_ids) const{
  shadowPart->getTriasIds(trias_ids);
}


/** ****************************************************************************
 * @brief
 * @param genface_ids
 ******************************************************************************/
void surfacePartDerivedVars::getGenfIds(vector<int> &genface_ids) const{
  shadowPart->getGenfIds(genface_ids);
}


/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void surfacePartDerivedVars::getPos(vector<vector3d<float> > &pos) const {
  shadowPart->getPos(pos);
}


/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void surfacePartDerivedVars::getPos(vector<vector3d<double> > &pos) const {
  shadowPart->getPos(pos);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePartDerivedVars::getNodalScalar(string varname,
                                            vector<float> &vals) const
{
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(varname);
  if(mi==derivedVars.end())
    shadowPart->getNodalScalar(varname,vals);
  else {
    derivedVar_t vartype = mi->second;
    switch(vartype) {
    case VAR_M:
      {
        vector<float> a;
        vector<vector3d<float> > v;
        shadowPart->getNodalScalar("a",a);
        shadowPart->getNodalVector("v",v);
        vector<float> m(a.size());
        for(size_t i=0;i<a.size();++i)
          m[i] = norm(v[i])/a[i];
        vals.swap(m);
      }
      break;
    case VAR_P:
    case VAR_logp:
      {
        vector<float> pg;
        shadowPart->getNodalScalar("pg",pg);
        vector<float> P(pg.size());
        for(size_t i=0;i<P.size();++i)
          P[i] = (vartype==VAR_logp)?log10(max(pg[i]+Pambient,1e-30f)):
            (pg[i]+Pambient);
        vals.swap(P);
      }
      break;
    case VAR_U:
    case VAR_0:
    case VAR_1:
    case VAR_2:
      {
        vector<vector3d<float> > v;
        shadowPart->getNodalVector("v",v);
        vector<float> tmp(v.size());
        for(size_t i=0;i<v.size();++i) {
          switch(vartype) {
          case VAR_U:
            tmp[i] = norm(v[i]);
            break;
          case VAR_0:
            tmp[i] = v[i].x;
            break;
          case VAR_1:
            tmp[i] = v[i].y;
            break;
          case VAR_2:
            tmp[i] = v[i].z;
            break;
          default:
            tmp[i] = 0;
          }
        }
        vals.swap(tmp);
      }
      break;
    case VAR_X:
    case VAR_Y:
    case VAR_Z:
      {
        vector<vector3d<float> > pos;
        shadowPart->getPos(pos);
        vector<float> tmp(pos.size());
        for(size_t i=0;i<pos.size();++i) {
          switch(vartype) {
          case VAR_X:
            tmp[i] = pos[i].x;
            break;
          case VAR_Y:
            tmp[i] = pos[i].y;
            break;
          case VAR_Z:
            tmp[i] =pos[i].z;
            break;
          default:
            tmp[i] = 0;
          }
        }
        vals.swap(tmp);
      }
      break;
    }
  }
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePartDerivedVars::getNodalVector(string varname,
                                            vector<vector3d<float> > &vals) const {
  shadowPart->getNodalVector(varname,vals);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePartDerivedVars::getElementScalar(string varname,
                                              vector<float> &qvals,
                                              vector<float> &tvals,
                                              vector<float> &gvals) const {
  shadowPart->getElementScalar(varname,qvals,tvals,gvals);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePartDerivedVars::getElementVector(string varname,
                                              vector<vector3d<float> > &qvals,
                                              vector<vector3d<float> > &tvals,
                                              vector<vector3d<float> > &gvals) const {
  shadowPart->getElementVector(varname,qvals,tvals,gvals);
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param triangles
 * @param tria_ids
 * @param quads
 * @param quad_ids
 * @param genface2n
 * @param gnodes
 * @param gen_ids
 ******************************************************************************/
surfacePartCopy::surfacePartCopy(string name,
                                 vector<Array<int,3> > &triangles,
                                 vector<int> &tria_ids,
                                 vector<Array<int,4> > &quads,
                                 vector<int> &quad_ids,
                                 vector<int> &genface2n,
                                 vector<int> &gnodes,
                                 vector<int> &gen_ids) {
  error = false;

  vector<int> node_set;
  for(size_t i=0;i<gnodes.size();++i)
    node_set.push_back(gnodes[i]);

  for(size_t i=0;i<triangles.size();++i) {
    node_set.push_back(triangles[i][0]);
    node_set.push_back(triangles[i][1]);
    node_set.push_back(triangles[i][2]);
  }
  for(size_t i=0;i<quads.size();++i) {
    node_set.push_back(quads[i][0]);
    node_set.push_back(quads[i][1]);
    node_set.push_back(quads[i][2]);
    node_set.push_back(quads[i][3]);
  }
  sort(node_set.begin(),node_set.end());
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end());

  map<int,int> nmap;
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1;
  }
  nodemap = node_set;
  trifaces = triangles;
  quadfaces = quads;
  nfacenodes = genface2n;
  gennodes = gnodes;
  triaIds = tria_ids;
  quadIds = quad_ids;
  genIds = gen_ids;

  for(size_t i=0;i<trifaces.size();++i) {
    trifaces[i][0] = nmap[trifaces[i][0]];
    trifaces[i][1] = nmap[trifaces[i][1]];
    trifaces[i][2] = nmap[trifaces[i][2]];
  }
  for(size_t i=0;i<quadfaces.size();++i) {
    quadfaces[i][0] = nmap[quadfaces[i][0]];
    quadfaces[i][1] = nmap[quadfaces[i][1]];
    quadfaces[i][2] = nmap[quadfaces[i][2]];
    quadfaces[i][3] = nmap[quadfaces[i][3]];
  }
  for(size_t i=0;i<gennodes.size();++i)
    gennodes[i] = nmap[gennodes[i]];

  partName = name;
  nnodes = nodemap.size();
  nquads = quadfaces.size();
  ntrias = trifaces.size();
  ngenf  = nfacenodes.size();
}


/** ****************************************************************************
 * @brief
 * @param posvol
 ******************************************************************************/
void surfacePartCopy::registerPos(const vector<vector3d<double> > &posvol) {
  vector<vector3d<double> > tmp(nodemap.size());
  pos.swap(tmp);
  for(size_t i=0;i<nodemap.size();++i)
    pos[i] = posvol[nodemap[i]-1];
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param val
 ******************************************************************************/
void surfacePartCopy::registerNodalScalar(string name, const vector<float> &val) {
  vector<float> tmp(nodemap.size());
  for(size_t i=0;i<nodemap.size();++i)
    tmp[i] = val[nodemap[i]-1];
  nodalScalars[name] = tmp;
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param val
 ******************************************************************************/
void surfacePartCopy::registerNodalVector(string name, const vector<vector3d<float> > &val) {
  vector<vector3d<float> > tmp(nodemap.size())
;
  for(size_t i=0;i<nodemap.size();++i)
    tmp[i] = val[nodemap[i]-1];
  nodalVectors[name] = tmp;
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param qval
 * @param tval
 * @param gval
 ******************************************************************************/
void surfacePartCopy::registerElementScalar(string name,
                      const vector<float> &qval,
                      const vector<float> &tval,
                      const vector<float> &gval)  {
  Array<vector<float>,3> &tmp = elementScalars[name];
  tmp[0] = qval;
  tmp[1] = tval;
  tmp[2] = gval;
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param qval
 * @param tval
 * @param gval
 ******************************************************************************/
void surfacePartCopy::registerElementVector(string name,
                      const vector<vector3d<float> > &qval,
                      const vector<vector3d<float> > &tval,
                      const vector<vector3d<float> > &gval)  {
  Array<vector<vector3d<float> >,3> &tmp = elementVectors[name];
  tmp[0] = qval;
  tmp[1] = tval;
  tmp[2] = gval;
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartCopy::hasNodalScalarVar(string var) const {
  map<string,vector<float> >::const_iterator mi;
  mi = nodalScalars.find(var);
  return !(mi == nodalScalars.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartCopy::hasNodalVectorVar(string var) const {
  map<string,vector<vector3d<float> > >::const_iterator mi;
  mi = nodalVectors.find(var);
  return !(mi == nodalVectors.end());
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartCopy::hasElementScalarVar(string var) const {
  map<string,Array<vector<float>,3> >::const_iterator mi;
  mi = elementScalars.find(var);
  return mi != elementScalars.end();
}


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool surfacePartCopy::hasElementVectorVar(string var) const {
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi;
  mi = elementVectors.find(var);
  return mi != elementVectors.end();
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartCopy::getNodalScalarVars() const {
  vector<string> tmp;
  map<string,vector<float> >::const_iterator mi;
  for(mi=nodalScalars.begin();mi!=nodalScalars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartCopy::getNodalVectorVars() const {
  vector<string> tmp;
  map<string,vector<vector3d<float> > >::const_iterator mi;
  for(mi=nodalVectors.begin();mi!=nodalVectors.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartCopy::getElementScalarVars() const {
  vector<string> tmp;
  map<string,Array<vector<float>,3> >::const_iterator mi;
  for(mi=elementScalars.begin();mi!=elementScalars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> surfacePartCopy::getElementVectorVars() const {
  vector<string> tmp;
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi;
  for(mi=elementVectors.begin();mi!=elementVectors.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief
 * @param quads
 ******************************************************************************/
void surfacePartCopy::getQuads(vector<Array<int,4> > &quads) const {
  quads = quadfaces;
}

/** ****************************************************************************
 * @brief
 * @param quads_ids
 ******************************************************************************/
void surfacePartCopy::getQuadsIds(vector<int> &quads_ids) const{
  quads_ids = quadIds;
}

/** ****************************************************************************
 * @brief
 * @param trias
 ******************************************************************************/
void surfacePartCopy::getTrias(vector<Array<int,3> > &trias) const {
  trias = trifaces;
}


/** ****************************************************************************
 * @brief
 * @param trias_ids
 ******************************************************************************/
void surfacePartCopy::getTriasIds(vector<int> &trias_ids) const{
  trias_ids = triaIds;
}


/** ****************************************************************************
 * @brief
 * @param numGenFnodes
 * @param genNodes
 ******************************************************************************/
void surfacePartCopy::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  numGenFnodes = nfacenodes;
  genNodes = gennodes;
}


/** ****************************************************************************
 * @brief
 * @param genface_ids
 ******************************************************************************/
void surfacePartCopy::getGenfIds(vector<int> &genface_ids) const{
  genface_ids = genIds;
}


/** ****************************************************************************
 * @brief
 * @param pos_out
 ******************************************************************************/
void surfacePartCopy::getPos(vector<vector3d<float> > &pos_out) const {
  pos_out.resize(pos.size());
  for(size_t i = 0; i < pos.size(); i++){
    pos_out[i].x = pos[i].x;
    pos_out[i].y = pos[i].y;
    pos_out[i].z = pos[i].z;
  }

}


/** ****************************************************************************
 * @brief
 * @param pos_out
 ******************************************************************************/
void surfacePartCopy::getPos(vector<vector3d<double> > &pos_out) const {
  pos_out = pos;
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePartCopy::getNodalScalar(string varname,
                                     vector<float> &vals) const {
  map<string,vector<float> >::const_iterator mi;
  mi = nodalScalars.find(varname);
  if(!(mi == nodalScalars.end()))
    vals = mi->second;
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void surfacePartCopy::getNodalVector(string varname,
                                     vector<vector3d<float> > &vals) const {
  map<string,vector<vector3d<float> > >::const_iterator mi;
  mi = nodalVectors.find(varname);
  if(!(mi == nodalVectors.end()))
    vals = mi->second;
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePartCopy::getElementScalar(string name,
                                       vector<float> &qvals,
                                       vector<float> &tvals,
                                       vector<float> &gvals) const {
  map<string,Array<vector<float>,3> >::const_iterator mi;
  mi = elementScalars.find(name);
  qvals = mi->second[0];
  tvals = mi->second[1];
  gvals = mi->second[2];
}


/** ****************************************************************************
 * @brief
 * @param name
 * @param qvals
 * @param tvals
 * @param gvals
 ******************************************************************************/
void surfacePartCopy::getElementVector(string name,
                                       vector<vector3d<float> > &qvals,
                                       vector<vector3d<float> > &tvals,
                                       vector<vector3d<float> > &gvals) const {
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi;
  mi = elementVectors.find(name);
  qvals = mi->second[0];
  tvals = mi->second[1];
  gvals = mi->second[2];
}

/** ****************************************************************************
 * @file      volumePart.cc
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

#include "fileFunctions.h"
#include "hdf5_Functions.h"
#include "volumePart.h"

using std::cerr;
using std::endl;
using std::cout;


/** ****************************************************************************
 * @brief
 * @param out_dir
 * @param iteration
 * @param casename
 * @param vars
 ******************************************************************************/
volumePart::volumePart(string out_dir, string iteration, string casename,
                       vector<string> vars)
{
  error = true;
  partName = "Volume";

  // Check number of nodes
  //-------------------------------------------------------------------------
  posFile = getPosFile(out_dir,iteration,casename);
  hid_t file_id = H5Fopen(posFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos");
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT);
#endif
  if(elg < 0) {
    H5Fclose(file_id);
    return;
  }
  nnodes = sizeElementType(elg,"data");
  H5Gclose(elg);
  H5Fclose(file_id);


  // Check for iblank information
  //-------------------------------------------------------------------------

  string iblankname = out_dir+"/grid_iblank." + iteration + "_" + casename;
  struct stat tmpstat;
  has_iblank = false;
  if(stat(iblankname.c_str(),&tmpstat)== 0) {
    file_id = Loci::hdf5OpenFile(iblankname.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    if(file_id < 0) {
      return;
    }
    fact_db facts;
    store<unsigned char> iblank_tmp;
    readData(file_id,"iblank",iblank_tmp.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);
    entitySet pdom = interval(1,nnodes);
    iblank.allocate(pdom);
    entitySet dom = iblank_tmp.domain();
    int cnt = 1;
    FORALL(dom,nd) {
      iblank[cnt++] = iblank_tmp[nd];
    } ENDFORALL;
    has_iblank = true;
  }

  // Check for element types in topo file
  //-------------------------------------------------------------------------
  topoFile = getTopoFileName(out_dir, casename, iteration);
  file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  elg = H5Gopen(file_id,"elements");
#else
  elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;

  ntets = sizeElementType(elg,"tetrahedra");
  nhexs = sizeElementType(elg,"hexahedra");
  nprsm = sizeElementType(elg,"prism");
  npyrm = sizeElementType(elg,"pyramid");
  ngenc = sizeElementType(elg,"GeneralCellNfaces");

  size_t ntets_b = ntets;
  size_t nhexs_b = nhexs;
  size_t nprsm_b = nprsm;
  size_t npyrm_b = npyrm;
  size_t ngenc_b = ngenc;
  const int block_size=65536; // Size of blocking factor
  if(has_iblank) {
    // need to adjust the number of elements based on iblanking.
    if(ntets > 0) {
      int nblocks = (ntets-1)/block_size+1;
      int remain = ntets;
      int start = 0;
      int cnt = 0;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain);
        vector<Array<int,4> > tets(size);
        readElementTypeBlock(elg,"tetrahedra",tets,start,size);
        remain -= size;
        start += size;
        for(int i=0;i<size;++i) {
          bool blank = true;
          for(int j=0;j<4;++j)
            if(iblank[tets[i][j]] < 2)
              blank = false;
          if(!blank)
            cnt++;
          else
            tetsIblanked += start-size+i;
        }
      }
      WARN(remain != 0);
      ntets_b = cnt;
      if(ntets-ntets_b > 0)
        cout << ntets-ntets_b << " tetrahedra iblanked" << endl;
    }
    if(nhexs > 0) {
      int nblocks = (nhexs-1)/block_size+1;
      int remain = nhexs;
      int start = 0;
      int cnt = 0;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain);

        vector<Array<int,8> > hexs(size);
        readElementTypeBlock(elg,"hexahedra",hexs,start,size);
        remain -= size;
        start += size;
        for(int i=0;i<size;++i) {
          bool blank = true;
          for(int j=0;j<8;++j)
            if(iblank[hexs[i][j]] < 2)
              blank = false;
          if(!blank)
            cnt++;
          else
            hexsIblanked += start-size+i;
        }
      }
      WARN(remain != 0);
      nhexs_b = cnt;
      if(nhexs-nhexs_b > 0)
        cout << nhexs-nhexs_b << " hexahedra iblanked" << endl;
    }
    if(nprsm > 0) {
      int nblocks = (nprsm-1)/block_size+1;
      int remain = nprsm;
      int start = 0;
      int cnt = 0;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain);
        vector<Array<int,6> > prsm(size);
        readElementTypeBlock(elg,"prism",prsm,start,size);
        remain -= size;
        start += size;
        for(int i=0;i<size;++i) {
          bool blank = true;
          for(int j=0;j<6;++j)
            if(iblank[prsm[i][j]] < 2)
              blank = false;
          if(!blank)
            cnt++;
          else
            prsmIblanked += start-size+i;
        }
      }
      WARN(remain != 0);
      nprsm_b = cnt;
      if(nprsm-nprsm_b > 0)
        cout << nprsm-nprsm_b << " prisms iblanked" << endl;

    }
    if(npyrm > 0) {
      int nblocks = (npyrm-1)/block_size+1;
      int remain = npyrm;
      int start = 0;
      int cnt = 0;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain);
        vector<Array<int,5> > pyrm(size);
        readElementTypeBlock(elg,"pyramid",pyrm,start,size);
        remain -= size;
        start += size;
        for(int i=0;i<size;++i) {
          bool blank = true;
          for(int j=0;j<5;++j)
            if(iblank[pyrm[i][j]] < 2)
              blank = false;
          if(!blank)
            cnt++;
          else
            pyrmIblanked += start-size+i;
        }
      }
      WARN(remain != 0);
      npyrm_b = cnt;
      if(npyrm-npyrm_b > 0)
        cout << npyrm-npyrm_b << " pyramids iblanked" << endl;
    }
    if(ngenc > 0) {
      vector<int> GeneralCellNfaces(ngenc);
      readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces);
      size_t nside = sizeElementType(elg,"GeneralCellNsides");
      vector<int> GeneralCellNsides(nside);
      readElementType(elg,"GeneralCellNsides",GeneralCellNsides);
      size_t nnodes = sizeElementType(elg,"GeneralCellNodes");
      vector<int> GeneralCellNodes(nnodes);
      readElementType(elg,"GeneralCellNodes",GeneralCellNodes);
      int cnt1 = 0;
      int cnt2 = 0;
      int cnt = 0;
      for(size_t i=0;i<ngenc;++i) {
        bool blank = true;
        int nf = GeneralCellNfaces[i];
        for(int f=0;f<nf;++f) {
          int fs = GeneralCellNsides[cnt1++];
          for(int n=0;n<fs;++n) {
            int nd = GeneralCellNodes[cnt2++];
            if(iblank[nd] < 2)
              blank = false;
          }
        }
        if(!blank)
          cnt++;
        else
          gencIblanked += (int)i;
      }
      ngenc_b = cnt;
      if(ngenc-ngenc_b > 0)
        cout << ngenc-ngenc_b << " general cells iblanked" << endl;
    }
    H5Gclose(elg);
    Loci::hdf5CloseFile(file_id);

  }
  ntets_orig = ntets;
  nhexs_orig = nhexs;
  nprsm_orig = nprsm;
  npyrm_orig = npyrm;
  ngenc_orig = ngenc;

  ntets = ntets_b;
  nhexs = nhexs_b;
  nprsm = nprsm_b;
  npyrm = npyrm_b;
  ngenc = ngenc_b;
  ntetsIblank = ntets_orig-ntets;
  nhexsIblank = nhexs_orig-nhexs;
  nprsmIblank = nprsm_orig-nprsm;
  npyrmIblank = npyrm_orig-npyrm;
  ngencIblank = ngenc_orig-ngenc;

  // Now check for variables
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i];
    string filename = out_dir+'/' + varname + "_sca." + iteration + "_" + casename;
    struct stat tmpstat;
    if(stat(filename.c_str(),&tmpstat) == 0) {
      nodalScalarVars[varname] = filename;
    } else {
      filename = out_dir+'/' + varname + "_vec." + iteration + "_" + casename;
      if(stat(filename.c_str(),&tmpstat) == 0) {
        nodalVectorVars[varname] = filename;
      } else {
      }
    }
  }

  error = false;
} // End of volumePart()


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool volumePart::hasNodalScalarVar(string var) const
{
  map<string,string>::const_iterator mi=nodalScalarVars.find(var);
  return (mi != nodalScalarVars.end());
} // End of hasNodalScalarVar()


/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool volumePart::hasNodalVectorVar(string var) const {
  map<string,string>::const_iterator mi=nodalVectorVars.find(var);
  return (mi != nodalVectorVars.end());
}


/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> volumePart::getNodalScalarVars() const  {
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
vector<string> volumePart::getNodalVectorVars() const {
  vector<string> tmp;
  map<string,string>::const_iterator mi;
  for(mi=nodalVectorVars.begin();mi!=nodalVectorVars.end();++mi)
    tmp.push_back(mi->first);
  return tmp;
}


/** ****************************************************************************
 * @brief Read the position file and write all position vectors to pos.
 * @param pos[out] Array of position vectors from posFile
 ******************************************************************************/
void volumePart::getPos(vector<vector3d<float> > &pos) const
{
  pos.clear();
  string filename = posFile;
  hid_t  file_id  = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos");
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT);
#endif
  if(elg < 0) {
    H5Fclose(file_id);
    return;
  }
  size_t nsz = sizeElementType(elg,"data");
  if(nsz != nnodes) {
    H5Gclose(elg);
    H5Fclose(file_id);
    return;
  }

  vector<vector3d<float> > tmp(nsz);
  pos.swap(tmp);
  readElementType(elg,"data",pos);
  H5Gclose(elg);
  H5Fclose(file_id);
} // End of getPos()


/** ****************************************************************************
 * @brief Read the position file and write all position vectors to pos.
 * @param pos[out] Array of position vectors from posFile
 ******************************************************************************/
void volumePart::getPos(vector<vector3d<double> > &pos) const
{
  pos.clear();
  string filename = posFile;
  hid_t  file_id  = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;

#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos");
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT);
#endif

  if(elg < 0) {
    H5Fclose(file_id);
    return;
  }
  size_t nsz = sizeElementType(elg,"data");
  if(nsz != nnodes) {
    H5Gclose(elg);
    H5Fclose(file_id);
    return;
  }

  vector<vector3d<double> > tmp(nsz);
  pos.swap(tmp);
  readElementType(elg,"data",pos);
  H5Gclose(elg);
  H5Fclose(file_id);
}


/** ****************************************************************************
 * @brief
 * @param tets
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const
{
  tets.clear();
  if(ntets <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,ntets_orig-start);
  vector<Array<int,4> > tets_local(lsize);
  readElementTypeBlock(elg,"tetrahedra",tets_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = tetsIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    tets.swap(tets_local);
    return;
  }
  // Remove iblanked cells
  vector<Array<int,4> > tets_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    tets_new[cnt] = tets_local[cp];
    cnt++;
  }ENDFORALL;
  tets.swap(tets_new);
}


/** ****************************************************************************
 * @brief
 * @param tetids
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getTetIds(vector<int> &tetids, size_t start, size_t size) const
{
  tetids.clear();
  if(ntets <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,ntets_orig-start);
  vector<int > tets_local(lsize);
  readElementTypeBlock(elg,"tetrahedra_ids",tets_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = tetsIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    tetids.swap(tets_local);
    return;
  }
  // Remove iblanked cells
  vector<int > tets_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    tets_new[cnt] = tets_local[cp];
    cnt++;
  }ENDFORALL;
  tetids.swap(tets_new);
}


/** ****************************************************************************
 * @brief
 * @param pyrms
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const
{
  pyrms.clear();
  if(npyrm <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,npyrm_orig-start);
  vector<Array<int,5> > pyrms_local(lsize);
  readElementTypeBlock(elg,"pyramid",pyrms_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = pyrmIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    pyrms.swap(pyrms_local);
    return;
  }
  // Remove iblanked cells
  vector<Array<int,5> > pyrms_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    pyrms_new[cnt] = pyrms_local[cp];
    cnt++;
  }ENDFORALL;
  pyrms.swap(pyrms_new);
}

/** ****************************************************************************
 * @brief
 * @param pyrmids
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const
{
  pyrmids.clear();
  if(npyrm <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,npyrm_orig-start);
  vector<int > pyrms_local(lsize);
  readElementTypeBlock(elg,"pyramid_ids",pyrms_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = pyrmIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    pyrmids.swap(pyrms_local);
    return;
  }
  // Remove iblanked cells
  vector<int > pyrms_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    pyrms_new[cnt] = pyrms_local[cp];
    cnt++;
  }ENDFORALL;
  pyrmids.swap(pyrms_new);
}

/** ****************************************************************************
 * @brief
 * @param prsms
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const
{
  prsms.clear();
  if(nprsm <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,nprsm_orig-start);
  vector<Array<int,6> > prsms_local(lsize);
  readElementTypeBlock(elg,"prism",prsms_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = prsmIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    prsms.swap(prsms_local);
    return;
  }
  // Remove iblanked cells
  vector<Array<int,6> > prsms_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    prsms_new[cnt] = prsms_local[cp];
    cnt++;
  }ENDFORALL;
  prsms.swap(prsms_new);
}


/** ****************************************************************************
 * @brief
 * @param prsmids
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const
{
  prsmids.clear();
  if(nprsm <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,nprsm_orig-start);
  vector<int > prsms_local(lsize);
  readElementTypeBlock(elg,"prism_ids",prsms_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = prsmIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    prsmids.swap(prsms_local);
    return;
  }
  // Remove iblanked cells
  vector<int > prsms_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    prsms_new[cnt] = prsms_local[cp];
    cnt++;
  }ENDFORALL;
  prsmids.swap(prsms_new);
}


/** ****************************************************************************
 * @brief
 * @param hexs
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const
{
  hexs.clear();
  if(nhexs <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,nhexs_orig-start);
  vector<Array<int,8> > hexs_local(lsize);
  readElementTypeBlock(elg,"hexahedra",hexs_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = hexsIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    hexs.swap(hexs_local);
    return;
  }
  // Remove iblanked cells
  vector<Array<int,8> > hexs_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    hexs_new[cnt] = hexs_local[cp];
    cnt++;
  }ENDFORALL;
  hexs.swap(hexs_new);
}

/** ****************************************************************************
 * @brief
 * @param hexids
 * @param start
 * @param size
 ******************************************************************************/
void volumePart::getHexIds(vector<int> &hexids, size_t start, size_t size) const
{
  hexids.clear();
  if(nhexs <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = min(size,nhexs_orig-start);
  vector<int > hexs_local(lsize);
  readElementTypeBlock(elg,"hexahedra_ids",hexs_local,start,lsize);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = hexsIblanked & interval(start,start+lsize-1);
  if(iblank==EMPTY) {
    hexids.swap(hexs_local);
    return;
  }
  // Remove iblanked cells
  vector<int > hexs_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    hexs_new[cnt] = hexs_local[cp];
    cnt++;
  }ENDFORALL;
  hexids.swap(hexs_new);
}


/** ****************************************************************************
 * @brief
 * @param genCellNfaces
 * @param genCellNsides
 * @param genCellNodes
 ******************************************************************************/
void volumePart::getGenCell(vector<int> &genCellNfaces,
                            vector<int> &genCellNsides,
                            vector<int> &genCellNodes) const
{
  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;

  vector<int> GeneralCellNfaces(ngenc_orig);
  readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces);
  size_t nside = sizeElementType(elg,"GeneralCellNsides");
  vector<int> GeneralCellNsides(nside);
  readElementType(elg,"GeneralCellNsides",GeneralCellNsides);
  size_t nnodes = sizeElementType(elg,"GeneralCellNodes");
  vector<int> GeneralCellNodes(nnodes);
  readElementType(elg,"GeneralCellNodes",GeneralCellNodes);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);
  // If no general cells iblanked, then return
  if(ngenc_orig == ngenc) {
    genCellNfaces.swap(GeneralCellNfaces);
    genCellNsides.swap(GeneralCellNsides);
    genCellNodes.swap(GeneralCellNodes);
    return;
  }
  int currentFaceOffset = 0;
  int currentNodeOffset = 0;
  int skip_cells = 0;
  int skip_faces = 0;
  int skip_nodes = 0;
  for(size_t i=0;i<ngenc_orig;++i) {
    bool blank = gencIblanked.inSet(i);
    // nf is number of faces for this general cell
    int nf = GeneralCellNfaces[i];
    // nn is the number of nodes firthis general cell
    int nn = 0;
    for(int f=0;f<nf;++f) {
      nn += GeneralCellNsides[currentFaceOffset+f];
    }
    if(blank) {
      skip_cells += 1;
      skip_faces += nf;
      skip_nodes += nn;
    } else {
      if(skip_cells > 0) {
        GeneralCellNfaces[i-skip_cells] = GeneralCellNfaces[i];
        for(int f=0;f<nf;++f)
          GeneralCellNsides[currentFaceOffset-skip_faces+f] =
            GeneralCellNsides[currentFaceOffset+f];
        for(int n=0;n<nn;++n)
          GeneralCellNodes[currentNodeOffset-skip_nodes+n] =
            GeneralCellNodes[currentNodeOffset+n];
      }
    }
    currentFaceOffset += nf;
    currentNodeOffset += nn;
  }

  GeneralCellNfaces.resize(GeneralCellNfaces.size()-skip_cells);
  GeneralCellNsides.resize(GeneralCellNsides.size()-skip_faces);
  GeneralCellNodes.resize(GeneralCellNodes.size()-skip_nodes);
  genCellNfaces.swap(GeneralCellNfaces);
  genCellNsides.swap(GeneralCellNsides);
  genCellNodes.swap(GeneralCellNodes);
}

/** ****************************************************************************
 * @brief
 * @param genids
 ******************************************************************************/
void volumePart::getGenIds(vector<int> &genids) const
{
  genids.clear();
  if(ngenc <=0)
    return;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements");
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int lsize = ngenc_orig;
  vector<int > gens_local(lsize);
  readElementType(elg,"GeneralCell_ids",gens_local);
  H5Gclose(elg);
  Loci::hdf5CloseFile(file_id);

  entitySet iblank = gencIblanked & interval(0,lsize-1);
  if(iblank==EMPTY) {
    genids.swap(gens_local);
    return;
  }
  // Remove iblanked cells
  vector<int > gens_new(lsize-iblank.size());
  int cnt = 0;
  entitySet dom = (~(iblank)) & interval(0,lsize-1);
  FORALL(dom,cp) {
    gens_new[cnt] = gens_local[cp];
    cnt++;
  }ENDFORALL;
  genids.swap(gens_new);
}

/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void volumePart::getNodalScalar(string varname, vector<float> &vals) const
{
  vals.clear();
  map<string,string>::const_iterator mi = nodalScalarVars.find(varname);
  if(mi == nodalScalarVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
  string vname = getVarNameFromFile(file_id,varname);
#ifdef VERBOSE
  if(vname != varname) {
    cerr << "reading var '" << vname << "' from file '" << filename << "'" << endl;
  }
#endif
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,vname.c_str());
#else
  hid_t elg = H5Gopen(file_id,vname.c_str(),H5P_DEFAULT);

#endif
  if(elg < 0)
    return;
  int nsz = sizeElementType(elg,"data");
  vector<float> tmp(nsz);
  vals.swap(tmp);
  readElementType(elg,"data",vals);
  H5Gclose(elg);
  H5Fclose(file_id);
}

/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void volumePart::getNodalVector(string varname, vector<vector3d<float> > &vals) const
{
  vals.clear();
  map<string,string>::const_iterator mi = nodalVectorVars.find(varname);
  if(mi == nodalVectorVars.end())
    return;
  string filename = mi->second;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id < 0)
    return;
  string vname = getVarNameFromFile(file_id,varname);
#ifdef VERBOSE
  if(vname != varname) {
    cerr << "reading var '" << vname << "' from file '" << filename << "'" << endl;
  }
#endif
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,vname.c_str());
#else
  hid_t elg = H5Gopen(file_id,vname.c_str(),H5P_DEFAULT);
#endif
  if(elg < 0)
    return;
  int nsz = sizeElementType(elg,"data");
  vector<vector3d<float> > tmp(nsz);
  vals.swap(tmp);
  readElementType(elg,"data",vals);
  H5Gclose(elg);
  H5Fclose(file_id);
}

/** ****************************************************************************
 * @brief
 * @param blank
 ******************************************************************************/
void volumePart::getNodalIblank(vector<unsigned char> &blank) const
{
  if(has_iblank) {
    Loci::entitySet idom = iblank.domain();
    vector<unsigned char> tmp(idom.size());
    int cnt = 0;
    FORALL(idom,ii) {
      tmp[cnt] = iblank[ii];
      cnt++;
    } ENDFORALL;
    blank.swap(tmp);
  } else {
    blank.clear();
  }
}

/** ****************************************************************************
 * @brief
 * @param vars
 ******************************************************************************/
void volumePartDerivedVars::processDerivedVars(const vector<string> &vars)
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
volumePartDerivedVars::volumePartDerivedVars(volumePartP part,
                                             string output_dir,
                                             string casename, string iteration,
                                             vector<string> vars)
{
  error = part->fail();
  partName = part->getPartName();
  nnodes = part->getNumNodes();
  ntets = part->getNumTets();
  nhexs = part->getNumHexs();
  nprsm = part->getNumPrsm();
  npyrm = part->getNumPyrm();
  ngenc = part->getNumGenc();
  ntetsIblank = part->getNumTetsIblank();
  nhexsIblank = part->getNumHexsIblank();
  nprsmIblank = part->getNumPrsmIblank();
  npyrmIblank = part->getNumPyrmIblank();
  ngencIblank = part->getNumGencIblank();

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
      cerr << "Unable to open file " << filename << endl;
    }
  }

  processDerivedVars(vars);
}

/** ****************************************************************************
 * @brief
 * @param var
 * @return true
 * @return false
 ******************************************************************************/
bool volumePartDerivedVars::hasNodalScalarVar(string var) const
{
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
bool volumePartDerivedVars::hasNodalVectorVar(string var) const
{
  return shadowPart->hasNodalVectorVar(var);
}

/** ****************************************************************************
 * @brief
 * @return vector<string>
 ******************************************************************************/
vector<string> volumePartDerivedVars::getNodalScalarVars() const
{
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
vector<string> volumePartDerivedVars::getNodalVectorVars() const
{
  return shadowPart->getNodalVectorVars();
}

/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void volumePartDerivedVars::getPos(vector<vector3d<float> > &pos) const
{
  shadowPart->getPos(pos);
}

/** ****************************************************************************
 * @brief
 * @param pos
 ******************************************************************************/
void volumePartDerivedVars::getPos(vector<vector3d<double> > &pos) const
{
  shadowPart->getPos(pos);
}

/** ****************************************************************************
 * @brief
 * @param tets
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const
{
  shadowPart->getTetBlock(tets,start,size);
}

/** ****************************************************************************
 * @brief
 * @param tetids
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getTetIds(vector<int> &tetids, size_t start, size_t size) const
{
  shadowPart->getTetIds(tetids,start,size);
}

/** ****************************************************************************
 * @brief
 * @param pyrms
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const
{
  shadowPart->getPyrmBlock(pyrms,start,size);
}

/** ****************************************************************************
 * @brief
 * @param pyrmids
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const
{
  shadowPart->getPyrmIds(pyrmids,start,size);
}

/** ****************************************************************************
 * @brief
 * @param prsms
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const
{
  shadowPart->getPrsmBlock(prsms,start,size);
}

/** ****************************************************************************
 * @brief
 * @param prsmids
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const
{
  shadowPart->getPrsmIds(prsmids,start,size);
}

/** ****************************************************************************
 * @brief
 * @param hexs
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const
{
  shadowPart->getHexBlock(hexs,start,size);
}

/** ****************************************************************************
 * @brief
 * @param hexids
 * @param start
 * @param size
 ******************************************************************************/
void volumePartDerivedVars::getHexIds(vector<int> &hexids, size_t start, size_t size) const
{
  shadowPart->getHexIds(hexids,start,size);
}

/** ****************************************************************************
 * @brief
 * @param genCellNfaces
 * @param genCellNsides
 * @param genCellNodes
 ******************************************************************************/
void volumePartDerivedVars::getGenCell(vector<int> &genCellNfaces,
                                       vector<int> &genCellNsides,
                                       vector<int> &genCellNodes) const
{
  shadowPart->getGenCell(genCellNfaces,genCellNsides,genCellNodes);
}

/** ****************************************************************************
 * @brief
 * @param genids
 ******************************************************************************/
void volumePartDerivedVars::getGenIds(vector<int> &genids) const
{
  shadowPart->getGenIds(genids);
}


/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void volumePartDerivedVars::getNodalScalar(string varname, vector<float> &vals) const
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
} // End of getNodalScalar()

/** ****************************************************************************
 * @brief
 * @param varname
 * @param vals
 ******************************************************************************/
void volumePartDerivedVars::getNodalVector(string varname, vector<vector3d<float> > &vals) const
{
  shadowPart->getNodalVector(varname,vals);
}

/** ****************************************************************************
 * @brief
 * @param blank
 ******************************************************************************/
void volumePartDerivedVars::getNodalIblank(vector<unsigned char> &blank) const
{
  shadowPart->getNodalIblank(blank);
}

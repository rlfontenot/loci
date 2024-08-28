/** ****************************************************************************
 * @file      extract.cc
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

#include <dirent.h>
#include <Loci.h>
#include "extract.h"
#include "fileFunctions.h"
#include "hdf5_Functions.h"

using std::cerr;
using std::endl;
using std::cout;

#pragma  GCC diagnostic ignored "-Wunused-variable"

string output_dir;

/** ****************************************************************************
 * @brief
 * @param ac
 * @param av
 ******************************************************************************/
void Usage(int ac, char *av[])
{
  cerr << av[0] << ": Incorrect Usage" << endl;
  cout << "Usage:" << endl;
  cout << av[0] << " <package> [package options] <case_name> <time_step> <variable(s)>" << endl;
  cout << endl;
  cout << "where <package> may be:" << endl
       << "-2d         : extract for the 2dgv plotting package"             << endl
       << "-fv         : extract for the FieldView post-processing package" << endl
       << "-en         : extract for the Ensight post-processing package"   << endl
       << "-en_with_id : extract for the Ensight post-processing package"   << endl
       << "              with node id and element id"                       << endl
       << "-cgns       : extract for the CGNS post-processing package"      << endl
       << "-tec        : extract for the TecPlot post-procesing package"    << endl
       << "-vtk        : extract for the Paraview post-procesing package"   << endl
       << "-vtk64      : extract for the Paraview post-procesing package"   << endl
       << "              (for large cases, must use >= Paraview 3.98)"      << endl
       << "-vtk_surf   : extract boundary surface mesh for the Paraview"    << endl
       << "              post-procesing package"                            << endl
       << "-vtk_surf64 : extract boundary surface mesh for the Paraview"    << endl
       << "              post-procesing package (for large cases, must "    << endl
       << "              use >= Paraview 3.98)"                             << endl
       << "-ascii      : extract to an ascii file"                          << endl
       << "-surf       : extract boundary surface mesh"                     << endl
       << "-cut        : extract a cutting plane for the 2dgv plotting"     << endl
       << "              package"                                           << endl
       << "-mean       : generate mean and variance from a family of ouput" << endl
       << "              variables"                                         << endl
       << "-combine    : combine mean and variance from online averaging"   << endl
       << "              ouput"                                             << endl
       << "-fcombine   : combine favre mean and variance from online"       << endl
       << "              averaging ouput"                                   << endl
       << endl
       << "Variables are defined by the solver, but typically include: "    << endl
       << "r          - nodal density"                                      << endl
       << "p          - nodal log10 pressure"                               << endl
       << "P          - nodal absolute pressure"                            << endl
       << "pg         - nodal gage pressure"                                << endl
       << "u          - nodal velocity magnitude"                           << endl
       << "m          - nodal mach number"                                  << endl
       << "t          - nodal temperature"                                  << endl
       << "a          - nodal soundspeed"                                   << endl
       << "f<species> - nodal mass fractions for <species>"                 << endl
       << "v          - nodal velocity vector"                              << endl
       << "x          - x coordinate"                                       << endl
       << "y          - y coordinate"                                       << endl
       << "z          - z coordinate"                                       << endl
       << endl
       << "Boundary Variables:"                                             << endl
       << "qdot  - wall boundary heat flux"                                 << endl
       << "yplus - wall boundary y plus"                                    << endl
       << "tau   - viscous wall shear stress vector"                        << endl
       << "tw    - viscous wall boundary temperature"                       << endl
       << "pw    - viscous wall boundary pressure"                          << endl
       << "n     - boundary normal (-ascii only)"                           << endl
       << "area  - boundary area   (-ascii only)"                           << endl
       << endl;

  cout << "extra options for particle extraction"                           << endl;
  cout << "  -mp <n> : maximum particles to extract. Default (if not"       << endl;
  cout << "            specified) is all particles available"               << endl;

  cout << "extra options for the 2dgv postprocessing package"               << endl
       << "  -bc <boundary_tag>  : specify which boundary to extract for"   << endl
       << "                        (multiple -bc's are merged)"             << endl
       << "  -xy : project boundary onto z=0 plane (default)"               << endl
       << "  -yz : project boundary onto x=0 plane"                         << endl
       << "  -xz : project boundary onto y=0 plane"                         << endl
       << "  -xr : project boundary onto x,radius plane (for axisymmetric"  << endl
       << "        grids)"                                                  << endl
       << endl
       << "extra options for orienting cutting plane"                       << endl
       << "  -xy          : originate transformations from z=0 plane (default)" << endl
       << "  -yz          : originate transformations from x=0 plane"       << endl
       << "  -xz          : originate transformations from y=0 plane"       << endl
       << "  -Rx <amount> : rotate cutting plane about x-axis"              << endl
       << "  -Ry <amount> : rotate cutting plane about y-axis"              << endl
       << "  -Rz <amount> : rotate cutting plane about z-axis"              << endl
       << "  -Sx <amount> : translate cutting plane along x-axis"           << endl
       << "  -Sy <amount> : translate cutting plane along y-axis"           << endl
       << "  -Sz <amount> : translate cutting plane along z-axis"           << endl
       << endl
       << "extra options for averaging feature '-mean'"                     << endl
       << "  -end <value> : ending iteration number for averaging"          << endl
       << "  -inc <value> : value to increment between iterations for"      << endl
       << "                 averaging"                                      << endl
       << endl
       << "extra options for controlling extract"                           << endl
       << "   -dir <directory> : change extract directory from default"     << endl
       << "                      'output'"                                  << endl
       << "   -skipPartDirs    : do not process the data in the part "      << endl
       << "                      subdirectories, if present. Default is to" << endl
       << "                      process the data in the part sub-"         << endl
       << "                      directories, if present."                  << endl
       << endl
       << "example:  to extract OH species from time step 50 of ssme"       << endl
       << "          simulation for visualization with 2dgv use:"           << endl
       << "   " << av[0] << " -2d -bc 1 -xr ssme 50 fOH"                    << endl
       << endl
       << "example: to extract an ascii table of boundary heat flux and"    << endl
       << "         x locations:"                                           << endl
       << "   " << av[0] << " -ascii -bc 4 nozzle 0 x qdot"                 << endl
       << endl
       << "example: to extract a cutting plane with various transformations:" << endl
       << "   " << av[0] << " -cut -xz -Sy 1.5 -Rx 30 -Rz -15 nozzle 0 P t" << endl
       << endl
       << "example: to extract 5000 particles with associated temperature"  << endl
       << "         from time step 100 for visualization with Ensight:"     << endl
       << "   " << av[0] << " -en combustor 100 particle_temp -mp 5000"     << endl
       << endl
       << "example: to compute mean and variance values for velocity and"   << endl
       << "         temperature for iterations 1000-3000 output every 100"  << endl
       << "         iterations run:"                                        << endl
       << "   " << av[0] << " -mean -end 3000 -inc 100 nozzle 1000 t v"     << endl
       << endl
       << "NOTE: outputs variabes tMean, tVar, vMean, vVar, vCuv, vCuw,"    << endl
       << "      vCvw at iteration 1000."                                   << endl;
  exit(-1);
} // End of Usage()


/** ****************************************************************************
 * @brief
 * @param dval
 * @param var_name
 * @param casename
 * @param iteration
 ******************************************************************************/
void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) {
  if(var_name == "m") {
    string filename = output_dir+"/a_sca."+iteration + "_" + casename;

    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    fact_db facts;
    store<float> soundSpeed;
    readData(file_id,"a",soundSpeed.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    filename = output_dir+"/v_vec." + iteration +"_" + casename;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    store<vector3d<float> > u;
    readData(file_id,"v",u.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    entitySet dom = u.domain();
    int c = 0;
    FORALL(dom,nd) {
      float m = norm(u[nd])/soundSpeed[nd];
      dval[c++] = m;
    } ENDFORALL;
  } else if(var_name == "p" || var_name == "P") {
    string filename = output_dir+"/pg_sca."+iteration + "_" + casename;

    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    fact_db facts;
    store<float> pg;
    readData(file_id,"pg",pg.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    filename = output_dir+"/Pambient_par." + iteration +"_" + casename;

    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    param<float> Pambient;
    readData(file_id,"Pambient",Pambient.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    entitySet dom = pg.domain();
    bool log = (var_name == "p");
    int c = 0;
    FORALL(dom,nd) {
      float p = pg[nd]+Pambient[nd];
      if(log)
        dval[c++] = log10(p);
      else
        dval[c++] = p;
    } ENDFORALL;
  } else if(var_name == "u") {
    fact_db facts;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    store<vector3d<float> > u;
    readData(file_id,"v",u.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    entitySet dom = u.domain();
    int c = 0;
    FORALL(dom,nd) {
      float m = norm(u[nd]);
      dval[c++] = m;
    } ENDFORALL;
  } else if(var_name == "x" || var_name =="y" || var_name == "z") {
    store<vector3d<float> > pos;
    string posname = getPosFile(output_dir,iteration,casename);
    hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to get grid positions for iteration " << iteration
           << endl;
      cerr << "does file '" << posname << "' exist?" << endl;
      Loci::Abort();
      exit(-1);
    }

    fact_db facts;
    readData(file_id,"pos",pos.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);
    entitySet dom = pos.domain();
    int c = 0;
    if(var_name == "x") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x;
      } ENDFORALL;
    }
    if(var_name == "y") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].y;
      } ENDFORALL;
    }
    if(var_name == "z") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].z;
      } ENDFORALL;
    }
  } else if(var_name == "0" || var_name =="1" || var_name == "2") {
    fact_db facts;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      return;
    }

    store<vector3d<float> > u;
    readData(file_id,"v",u.Rep(),EMPTY,facts);
    Loci::hdf5CloseFile(file_id);

    entitySet dom = u.domain();
    int c = 0;
    if(var_name == "0") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].x;
      } ENDFORALL;
    }
    if(var_name == "1") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].y;
      } ENDFORALL;
    }
    if(var_name == "2") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].z;
      } ENDFORALL;
    }
  } else {
    cerr << "don't know how to get derived variable " << var_name << endl;
  }
}


/** ****************************************************************************
 * @brief
 * @param output_dir
 * @param iteration
 * @param casename
 * @return vector<string>
 ******************************************************************************/
vector<string> volumeSurfaceNames(string output_dir, string iteration,
                                  string casename) {
  string gridtopo = getTopoFileName(output_dir, casename, iteration);
  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries");
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT);
#endif
  hsize_t num_bcs = 0;
  H5Gget_num_objs(bndg,&num_bcs);
  vector<string>  bc_names;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024];
    memset(buf, '\0', 1024);
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf));
    buf[1023]='\0';
    bc_names.push_back(string(buf));
  }
  H5Gclose(bndg);
  H5Fclose(file_id);
  return bc_names;
}

/** ****************************************************************************
 * @brief
 * @param volSurface
 * @param vp
 * @param output_dir
 * @param iteration
 * @param casename
 * @param varlist
 ******************************************************************************/
void extractVolumeSurfaces(vector<surfacePartP> &volSurface,
                           volumePartP vp,
                           string output_dir,
                           string iteration,
                           string casename,
                           vector<string> varlist) {
  string gridtopo = getTopoFileName(output_dir, casename, iteration);

  cout << "extracting topology from '" << gridtopo << "'" << endl;

  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

  map<string,int> elementScalars;
  map<string,int> elementVectors;
  vector<string> variable_file(varlist.size());
  vector<map<int,int> > entityMap(varlist.size());
  vector<vector<int> > entityIds(varlist.size());
  vector<vector<float> > scalarElementVars(varlist.size());
  vector<vector<vector3d<float> > > vectorElementVars(varlist.size());

  for(size_t i = 0;i<varlist.size();++i) {
    string var = varlist[i];
    string filename = output_dir+'/' + var + "_bnd." + iteration + "_" + casename;
    struct stat tmpstat;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      elementScalars[var] = i;
      variable_file[i] = filename;
      continue;
    }
    filename = output_dir+'/' + var + "_bndvec." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      elementVectors[var] = i;
      variable_file[i] = filename;
      continue;
    }
  }

  map<string,int>::const_iterator mi;
  for(mi=elementScalars.begin();mi!=elementScalars.end();++mi) {
    int id = mi->second;
    string filename(variable_file[id]);
    string varname = varlist[id];
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);

    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      continue;
    }

#ifdef H5_USE_16_API
    hid_t di = H5Gopen(file_id,"dataInfo");
#else
    hid_t di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT);
#endif
    size_t nbel = sizeElementType(di,"entityIds");

    vector<int> elemIds(nbel);
    readElementType(di,"entityIds",elemIds);

    H5Gclose(di);
    vector<float> var(nbel);
    readElementType(file_id,varname.c_str(),var);

    entityIds[id] = elemIds;
    for(size_t i=0;i<nbel;++i) {
      entityMap[id][elemIds[i]] = int(i);
    }
    scalarElementVars[id] = var;
    H5Fclose(file_id);
  }

  for(mi=elementVectors.begin();mi!=elementVectors.end();++mi) {
    int id = mi->second;
    string filename(variable_file[id]);
    string varname = varlist[id];
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT);
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl;
      continue;
    }

#ifdef H5_USE_16_API
    hid_t di = H5Gopen(file_id,"dataInfo");
#else
    hid_t di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT);
#endif
    size_t nbel = sizeElementType(di,"entityIds");

    vector<int> elemIds(nbel);
    readElementType(di,"entityIds",elemIds);

    H5Gclose(di);
    vector<vector3d<float> > var(nbel);
    readElementType(file_id,varname.c_str(),var);

    entityIds[id] = elemIds;
    for(size_t i=0;i<nbel;++i) {
      entityMap[id][elemIds[i]] = int(i);
    }
    vectorElementVars[id] = var;
    H5Fclose(file_id);
  }

#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries");
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT);
#endif
  hsize_t num_bcs = 0;
  H5Gget_num_objs(bndg,&num_bcs);
  vector<string>  bc_names;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024];
    memset(buf, '\0', 1024);
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf));
    buf[1023]='\0';
    bc_names.push_back(string(buf));
  }
  vector<surfacePartCopy * > surfaceWork(num_bcs);
  for(hsize_t bc=0;bc<num_bcs;++bc) {
#ifdef H5_USE_16_API
    hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str());
#else
    hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str(),H5P_DEFAULT);
#endif

    size_t nquads = sizeElementType(bcg,"quads");
    size_t ntrias = sizeElementType(bcg,"triangles");
    size_t ngeneral = sizeElementType(bcg,"nside_sizes");

    vector<Array<int,3> > trias(ntrias);
    readElementType(bcg,"triangles",trias);
    vector<Array<int,4> > quads(nquads);
    readElementType(bcg,"quads",quads);

    vector<int> nside_sizes(ngeneral);
    readElementType(bcg,"nside_sizes",nside_sizes);
    size_t nside_nodes_size = sizeElementType(bcg,"nside_nodes");
    vector<int> nside_nodes(nside_nodes_size);
    readElementType(bcg,"nside_nodes",nside_nodes);

    vector<int > trias_id(ntrias);
    readElementType(bcg,"triangles_id",trias_id);
    vector<int > quads_id(nquads);
    readElementType(bcg,"quads_id",quads_id);
    vector<int > nside_id(ngeneral);
    readElementType(bcg,"nside_id",nside_id);

    vector<unsigned char> iblank;
    vp->getNodalIblank(iblank);

    if(iblank.size() > 0) {
      int cnt = 0;
      for(size_t i=0;i<ntrias;++i) {
        bool blank = true;
        for(int j=0;j<3;++j)
          if(iblank[trias[i][j]-1] < 2)
            blank = false;
        if(blank)
          cnt++;
        else
          if(cnt != 0) { // If there are some blanked copy into place
            trias[i-cnt]=trias[i];
            trias_id[i-cnt] = trias_id[i];
          }
      }
      if(cnt > 0) {
        size_t newsz = trias.size()-cnt;
        trias.resize(newsz);
        trias_id.resize(newsz);
      }
      cnt = 0;
      for(size_t i=0;i<nquads;++i) {
        bool blank = true;
        for(int j=0;j<4;++j)
          if(iblank[quads[i][j]-1] < 2)
            blank = false;
        if(blank)
          cnt++;
        else
          if(cnt != 0) { // If there are some blanked copy into place
            quads[i-cnt]=quads[i];
            quads_id[i-cnt] = quads_id[i];
          }
      }
      if(cnt > 0) {
        size_t newsz = quads.size()-cnt;
        quads.resize(newsz);
        quads_id.resize(newsz);
      }
      cnt  = 0;
      int cnt2 = 0;
      int nside_off = 0;
      for(size_t i=0;i<ngeneral;++i) {
        bool blank = true;
        for(int j=0;j<nside_sizes[i];++j) {
          if(iblank[nside_nodes[nside_off+j]-1] < 2)
            blank = false;
        }
        if(blank) {
          cnt++;
          cnt2 += nside_sizes[i];
        } else {
          if(cnt != 0) {
            nside_sizes[i-cnt] = nside_sizes[i];
            nside_id[i-cnt] = nside_id[i];
            for(int j=0;j<nside_sizes[i];++j)
              nside_nodes[nside_off-cnt2+j] = nside_nodes[nside_off+j];
          }
        }
        nside_off += nside_sizes[i];
      }
      if(cnt > 0) {
        size_t newsz = nside_sizes.size()-cnt;
        nside_sizes.resize(newsz);
        nside_id.resize(newsz);
        size_t newsz2 = nside_nodes.size()-cnt2;
        nside_nodes.resize(newsz2);
      }
    }
    surfaceWork[bc] =
      new surfacePartCopy(bc_names[bc],trias,trias_id, quads,quads_id, nside_sizes,nside_nodes, nside_id);

    for(mi=elementScalars.begin();mi!=elementScalars.end();++mi) {
      int id = mi->second;
      string varname = varlist[id];
      vector<float> qvals(quads.size());
      vector<float> tvals(trias.size());
      vector<float> gvals(nside_sizes.size());
      bool valid = true;
      for(size_t i=0;i<quads.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(quads_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
          break;
        } else {
          qvals[i] = scalarElementVars[id][ii->second];
        }
      }
      for(size_t i=0;i<trias.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(trias_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
          break;
        } else {
          tvals[i] = scalarElementVars[id][ii->second];
        }
      }
      for(size_t i=0;i<nside_sizes.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(nside_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
          break;
        } else {
          gvals[i] = scalarElementVars[id][ii->second];
        }
      }
      if(valid)
        surfaceWork[bc]->registerElementScalar(varname,qvals,tvals,gvals);
    }

    for(mi=elementVectors.begin();mi!=elementVectors.end();++mi) {
      int id = mi->second;
      string varname = varlist[id];
      vector<vector3d<float> > qvals(quads.size());
      vector<vector3d<float> > tvals(trias.size());
      vector<vector3d<float> > gvals(nside_sizes.size());
      bool valid = true;
      for(size_t i=0;i<quads.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(quads_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
        } else {
          qvals[i] = vectorElementVars[id][ii->second];
        }
      }
      for(size_t i=0;i<trias.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(trias_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
        } else {
          tvals[i] = vectorElementVars[id][ii->second];
        }
      }
      for(size_t i=0;i<nside_sizes.size();++i) {
        map<int,int>::const_iterator ii = entityMap[id].find(nside_id[i]);
        if(ii==entityMap[id].end()) {
          valid = false;
        } else {
          gvals[i] = vectorElementVars[id][ii->second];
        }
      }

      if(valid)
        surfaceWork[bc]->registerElementVector(varname,qvals,tvals,gvals);
    }

  }
  H5Gclose(bndg);
  H5Fclose(file_id);

  {
    vector<vector3d<double> > posvol;
    vp->getPos(posvol);
    for(hsize_t bc=0;bc<num_bcs;++bc) {
      surfaceWork[bc]->registerPos(posvol);
    }
  }
  {
    vector<string> vars = vp->getNodalScalarVars();
    for(size_t i=0;i<vars.size();++i) {
      vector<float> val;
      vp->getNodalScalar(vars[i],val);
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        surfaceWork[bc]->registerNodalScalar(vars[i],val);
      }
    }
  }
  {
    vector<string> vars = vp->getNodalVectorVars();
    for(size_t i=0;i<vars.size();++i) {
      vector<vector3d<float> > val;
      vp->getNodalVector(vars[i],val);
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        surfaceWork[bc]->registerNodalVector(vars[i],val);
      }
    }
  }
  volSurface = vector<surfacePartP>(num_bcs);
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    volSurface[bc] = surfaceWork[bc];
  }

}


namespace Loci {
  void disableDebugDir();
}

/** ****************************************************************************
 * @brief
 *
 * @param ac
 * @param av
 * @return int
 ******************************************************************************/
int main(int ac, char *av[])
{
  output_dir = "output";
  Loci::disableDebugDir();
  Loci::Init(&ac,&av);

  ePlot_type plot_type = NONE;

  bool           found_casename  = false;
  bool           found_iteration = false;
  bool           skip_part_dirs  = false;
  string         casename;
  string         iteration;
  vector<string> variables;
  vector<string> boundaries;
  float          xShift = 0.0;
  float          yShift = 0.0;
  float          zShift = 0.0;
  float          temp;
  int            view   = VIEWXY;
  affineMapping  transformMatrix;

  // record the maximum particle number to extract
  // a value < 0 means that there is no maximum particle
  // number limit, i.e., all particles are to be extracted
  // default is to extract all particles.  users can use
  // command line switch "-mp <n>" to set the maximum number
  // of particles to be extracted. if the requested particle
  // number is larger than the available particle number, then
  // all particles will be extracted.
  int  max_particles = -1;
  int  end_iter      = -1;
  int  inc_iter      = -1;
  bool id_required   = false; //ensight has the option to display node and element ids

  vector<string> partlist;
  bool novolume = false;
  for(int i=1;i<ac;++i)
  {
    if(av[i][0] == '-')
    {
      if(!strcmp(av[i],"-ascii"))         plot_type = ASCII;
      else if(!strcmp(av[i],"-surf"))     plot_type = SURFACE;
      else if(!strcmp(av[i],"-mean"))     plot_type = MEAN;
      else if(!strcmp(av[i],"-combine"))  plot_type = COMBINE;
      else if(!strcmp(av[i],"-fcombine")) plot_type = FCOMBINE;
      else if(!strcmp(av[i],"-2d"))       plot_type = TWODGV;
      else if(!strcmp(av[i],"-en"))       plot_type = ENSIGHT;
      else if(!strcmp(av[i],"-fv"))       plot_type = FIELDVIEW;
      else if(!strcmp(av[i],"-en_with_id"))
      {
        plot_type   = ENSIGHT;
        id_required = true;
      }else if(!strcmp(av[i],"-cgns"))
      {
        plot_type   = CGNS;
        id_required = true;
      }else if(!strcmp(av[i],"-tec"))
      {
        plot_type = TECPLOT;
#ifndef USE_NATIVE_TECPLOT
        cerr << "Note, This compiled version is using the older ASCII tecplot format." << endl
             << "If you are using a recent version you can configure Loci to use" << endl
             << "the native binary tecplot format that will be more effective for use with" << endl
             << "tecplot360. " << endl;
#endif
      }
      else if(!strcmp(av[i],"-vtk"))        plot_type = VTK;
      else if(!strcmp(av[i],"-vtk_surf"))   plot_type = VTK_SURFACE;
      else if(!strcmp(av[i],"-vtk64"))      plot_type = VTK64;
      else if(!strcmp(av[i],"-vtk_surf64")) plot_type = VTK_SURFACE64;
      else if(!strcmp(av[i],"-cut"))        plot_type = CUTTINGPLANE;
      else if(!strcmp(av[i],"-Sx"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> xShift).fail())
          Usage(ac, av);
      }else if(!strcmp(av[i],"-Sy"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> yShift).fail())
          Usage(ac, av);
      }else if(!strcmp(av[i],"-Sz"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> zShift).fail())
          Usage(ac, av);
      }else if(!strcmp(av[i],"-Rx"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> temp).fail())
          Usage(ac, av);
        transformMatrix.rotateX(-temp);
      }else if(!strcmp(av[i],"-Ry"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> temp).fail())
          Usage(ac, av);
        transformMatrix.rotateY(-temp);
      }else if(!strcmp(av[i],"-Rz"))
      {
        i++;
        std::istringstream iss(av[i]);
        if((iss >> std::dec >> temp).fail())
          Usage(ac, av);
        transformMatrix.rotateZ(-temp);
      }
      else if(!strcmp(av[i],"-xy"))  view=VIEWXY;
      else if(!strcmp(av[i],"-yz"))
      {
        view=VIEWYZ;
        transformMatrix.rotateY(90.0);
        transformMatrix.rotateZ(90.0);
      }else if(!strcmp(av[i],"-xz"))
      {
        view=VIEWXZ;
        transformMatrix.rotateX(90.0);
      }
      else if(!strcmp(av[i],"-xr"))       view=VIEWXR;
      else if(!strcmp(av[i],"-novolume")) novolume = true;
      else if(!strcmp(av[i],"-bc"))
      {
        i++;
        string v(av[i]);
        if(av[i][0] >='0' && av[i][0] <= '9')
          v = "BC_"+v;
        boundaries.push_back(v);
        partlist.push_back(v);
      }else if(!strcmp(av[i],"-part"))
      {
        i++;
        string v(av[i]);
        partlist.push_back(v);
      }else if(!strcmp(av[i],"-mp"))
      {
        // get the number of particles
        ++i;
        string n(av[i]);
        if(!valid_int(n))
        {
          cerr << "argument followed option '-mp' is not an integer,"
               << " used default value" << endl;
        }else
        {
          max_particles = str2int(n);
        }
      }else if(!strcmp(av[i],"-dir"))
      {
        ++i;
        output_dir = string(av[i]);
      }else if(!strcmp(av[i],"-inc"))
      {
        ++i;
        inc_iter = atoi(av[i]);
      }else if(!strcmp(av[i],"-end"))
      {
        ++i;
        end_iter = atoi(av[i]);
      } else if(!strcmp(av[i],"-skipPartDirs") ||
                !strcmp(av[i],"-onlyvolume") ) 
      {
        skip_part_dirs = true;
      }else 
      {
        cerr << "unknown option " << av[i] << endl;
        Usage(ac,av);
      }
    }else
    {
      if(found_iteration)
      {
        variables.push_back(string(av[i]));
      }else if(found_casename)
      {
        iteration       = string(av[i]);
        found_iteration = true;
      }else
      {
        casename       = string(av[i]);
        found_casename = true;
      }
    } // End If(av[i][0] == '-')
  } // End For(ac)

  if(boundaries.size() > 0)
    novolume = true;
  if(plot_type == NONE)
  {
    Usage(ac,av);
  }

  // if output directory doesn't exist, create one
  struct stat statbuf;
  if(stat(output_dir.c_str(), &statbuf))
  {
    mkdir(output_dir.c_str(), 0755);
  }else
  {
    if(!S_ISDIR(statbuf.st_mode))
    {
      cerr << "file '" << output_dir
           <<"' should be a directory!, rename 'output' and start again."
           << endl;
      Loci::Abort();
    }
  }

  // scan for parts
  if(partlist.size() == 0)
  {
    std::set<string> partfind;
    DIR *dp = opendir(output_dir.c_str());

    // Look in output directory and find all variables
    if(dp == 0)
    {
      cerr << "unable to open directory '" << output_dir << "'" << endl;
      exit(-1);
    }
    dirent *entry = readdir(dp);

    while(entry != 0)
    {
      string postfix;
      string vname;
      string vtype;
      string filename     = entry->d_name;
      string searchHeader = casename+"_SURF.";
      string filesub      = filename.substr(0,searchHeader.size());
      if(searchHeader == filesub && !skip_part_dirs)
      {
        size_t len       = filename.size()-searchHeader.size();
        string partname  = filename.substr(searchHeader.size(),len);
        string filecheck = output_dir + "/" + filename + "/topo_file." + iteration;
        struct stat tmpstat;

        if(stat(filecheck.c_str(),&tmpstat)== 0)
        {
          partfind.insert(partname);
          partlist.push_back(partname);
        }
      } // End If(searchHeader)
      entry = readdir(dp);
    } // End While(entry)

    closedir(dp);

    vector<string> vsurfs = volumeSurfaceNames(output_dir, iteration, casename);
    for(size_t i=0; i<vsurfs.size(); ++i)
    {
      if(partfind.find(vsurfs[i]) == partfind.end())
      {
        partlist.push_back(vsurfs[i]);
      }
    }
  } // End If(partlist.size)

  if(variables.size() == 0)
  {
    std::set<string>  varset;
    DIR              *dp = opendir(output_dir.c_str());
    // Look in output directory and find all variables
    if(dp == 0)
    {
      cerr << "unable to open directory '" << output_dir << "'" << endl;
      exit(-1);
    }

    dirent *entry = readdir(dp);
    string  tail  = iteration + "_" + casename;
    while(entry != 0)
    {
      string filename = entry->d_name;
      string postfix;
      string vname;
      string vtype;
      int    nsz = filename.size();
      int    dot = -1;
      for(int i=nsz-1;i>=0;--i)
      {
        if(filename[i] == '.')
        {
          dot = i;
          break;
        }
      }

      for(int i=dot+1;i<nsz;++i)
        postfix += filename[i];

      int und = -1;
      if(dot > 0)
      {
        for(int i=dot-1;i>=0;--i)
        {
          if(filename[i] == '_')
          {
            und = i;
            break;
          }
        }

        if(und > 0)
        {
          for(int i=und+1;i<dot;++i) vtype += filename[i];
          for(int i=0;i<und;++i)     vname += filename[i];
        }else
        {
          for(int i=0;i<dot;++i)
            vname += filename[i];
        }
      } // End If(dot)

      // Add derived variables
      if(dot>0 && und>0 && postfix == tail)
      {
        if(vtype == "sca"    || vtype == "vec"   || vtype == "bnd" ||
           vtype == "bndvec" || vtype == "ptsca" || vtype == "ptvec")
        {
          varset.insert(vname);
        }

        if(vtype == "sca" && vname == "pg")
        {
          varset.insert("P");
        }

        if(vtype == "sca" && vname == "a")
        {
          varset.insert("m");
        }
      } // End If(dot && und && postfix)
      entry = readdir(dp);
    } // End While(entry)

    closedir(dp);

    // Now check each part for variables
    if(partlist.size() > 0)
    {
      for(size_t i=0; i<partlist.size(); ++i)
      {
        string  dirname = output_dir+"/"+casename+"_SURF."+partlist[i];
        DIR    *dp      = opendir(dirname.c_str());

        // Look in output directory and find all variables
        if(dp == 0 || skip_part_dirs)
        {
          continue;
        }

        dirent *entry = readdir(dp);
        for( ; entry != 0; entry=readdir(dp))
        {
          string filename = entry->d_name;
          int    fsz      = filename.size();
          int    isz      = iteration.size();
          if(fsz <= isz)
            continue;

          string fiter = filename.substr(fsz-(isz+1),isz+1);

          if(fiter != string("."+iteration))
            continue;

          string remainder = filename.substr(0,fsz-(isz+1));
          int    remsz     = remainder.size();
          if(remsz <= 4)
            continue;

          string postfix = remainder.substr(remsz-4,4);
          if(postfix == "_sca" || postfix == "_vec")
          {
            string vname = remainder.substr(0,remsz-4);
            varset.insert(vname);
            continue;
          }

          if(remsz <= 5)
            continue;

          postfix = remainder.substr(remsz-5,5);
          if(postfix == "_bsca" || postfix == "_bvec")
          {
            string vname = remainder.substr(0,remsz-5);
            varset.insert(vname);
          }
        } // End For(entry)
        closedir(dp);
      } // End For(partlist.size)
    } // End If(partlist.size)

    std::set<string>::const_iterator vi;
    for(vi=varset.begin();vi!=varset.end();++vi)
      variables.push_back(*vi);

    if(variables.size() == 0)
    {
      variables.push_back("x");
      variables.push_back("y");
      variables.push_back("z");
    }

    cout << "extracting variables: ";
    for(size_t i=0; i<variables.size(); ++i)
    {
      cout << ' ' << variables[i];
    }
    cout << endl;
  } // End If(variables.size)

  // Find out variable types and variable files
  vector<int>    variable_type(variables.size());
  vector<string> variable_file(variables.size());
  bool           particle_info_requested = false;

  for(size_t i=0; i<variables.size(); ++i)
  {
    const string var(variables[i]);
    string       filename = output_dir+'/' + var + "_hdf5." + iteration;
    struct stat  tmpstat;
    if(stat(filename.c_str(), &tmpstat) == 0)
    {
      variable_type[i] = NODAL_SCALAR;
      variable_file[i] = filename;
      continue;
    }

    filename = output_dir+'/' + var + "_sca." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)== 0)
    {
      variable_type[i] = NODAL_SCALAR;
      variable_file[i] = filename;
      continue;
    }

    filename = output_dir+'/' + var + "_vec." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)== 0)
    {
      variable_type[i] = NODAL_VECTOR;
      variable_file[i] = filename;
      continue;
    }

    filename = output_dir+'/' + var + "_bnd." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)== 0)
    {
      variable_type[i] = BOUNDARY_SCALAR;
      variable_file[i] = filename;
      continue;
    }

    filename = output_dir+'/' + var + "_bndvec." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)== 0)
    {
      variable_type[i] = BOUNDARY_VECTOR;
      variable_file[i] = filename;
      continue;
    }

    filename = output_dir+'/' + var + "_ptsca." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)==0)
    {
      variable_type[i]        = PARTICLE_SCALAR;
      variable_file[i]        = filename;
      particle_info_requested = true;
      continue;
    }

    filename = output_dir+'/' + var + "_ptvec." + iteration + "_" + casename;
    if(stat(filename.c_str(),&tmpstat)==0)
    {
      variable_type[i]        = PARTICLE_VECTOR;
      variable_file[i]        = filename;
      particle_info_requested = true;
      continue;
    }

    if(plot_type == ASCII && boundaries.size() > 0)
    {
      if(var == "n")
      {
        variable_type[i] = BOUNDARY_DERIVED_VECTOR;
        continue;
      }

      if(var == "area" || var == "x" || var == "y" || var == "z")
      {
        variable_type[i] = BOUNDARY_DERIVED_SCALAR;
        continue;
      }
    } // End If(plot_type)

    if(var == "m")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    if(var == "p")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    if(var == "P")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    if(var == "u")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    if(var == "0" || var == "1" || var == "2")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    if(var == "x" || var == "y" || var == "z")
    {
      variable_type[i] = NODAL_DERIVED;
      continue;
    }

    // ... other derived variables here

    if(partlist.size() == 0)
      cerr << "Warning, variable '" << var
           << "' is unknown and will not be processed." << endl;
    variable_type[i] = UNDEFINED;
  } // End For(variables.size)


  // we will first check to see if particle position is present
  // in case of any particle information extraction
  if(particle_info_requested)
  {
    string      filename = output_dir +"/particle_pos." + iteration + "_" + casename;
    struct stat tmpstat;
    if(stat(filename.c_str(),&tmpstat)!=0)
    {
      cerr << "Warning: particle geometry '" << filename << "' must be "
           << "presented for any particle related variable extraction." << endl;
      Loci::Finalize();
      exit(-1);
    }
  } // End If(particle_info_requested)

#ifdef H5_USE_16_API
  H5Eset_auto(NULL,NULL);
#else
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);
#endif

  struct stat tmpstat;
  string      topoFile      = getTopoFileName(output_dir, casename, iteration);
  bool        timesyncerror = false;

  if(stat(topoFile.c_str(), &tmpstat) == 0)
  {
    struct stat gridstat;
    string      gridfile = casename + ".vog";

    if(stat(gridfile.c_str(), &gridstat) == 0)
    {
      if(gridstat.st_mtime > tmpstat.st_mtime)
        timesyncerror = true;
    }
  } // End If(topoFile)

  if(timesyncerror)
  {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl;
  } // End If(timesyncerror)

  if(plot_type == ASCII)
  {
    if(boundaries.size() == 0)
    {
      process_ascii_nodal(casename, iteration, variables, variable_type, variable_file);
      Loci::Finalize();
      exit(0);
    }else
    {
      process_ascii_bndry(casename, iteration, variables, variable_type, variable_file,
                          boundaries);
      Loci::Finalize();
      exit(0);
    }
  } // End If(plot_type)

  if(plot_type == MEAN) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
           << "       and option -inc to specify iteration increment value for iterations" << endl
           << "       to specify which files to average!" << endl;
      Loci::Finalize();
      exit(-1);
    }

    process_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter);
    Loci::Finalize();
    exit(0);
  }
  if(plot_type == COMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
           << "       and option -inc to specify iteration increment value for iterations" << endl
           << "       to specify which files to average!" << endl;
      Loci::Finalize();
      exit(-1);
    }

    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,false);
    Loci::Finalize();
    exit(0);
  }
  if(plot_type == FCOMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
           << "       and option -inc to specify iteration increment value for iterations" << endl
           << "       to specify which files to average!" << endl;
      Loci::Finalize();
      exit(-1);
    }

    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,true);
    Loci::Finalize();
    exit(0);
  }
  if(plot_type == TWODGV) {
    if(variables.size() != 1) {
      cerr << "2dgv extract can only extract one variable at a time."
           << endl;
      Usage(ac,av);
    }
    if(boundaries.size() == 0) {
      cerr << "2dgv extract must have the projected boundaries identified using the '-bc' flag" << endl;
      Usage(ac,av);
    }
    get_2dgv(casename,iteration,variables,variable_type,variable_file,
             boundaries,view);
    Loci::Finalize();
    exit(0);
  }
  if(plot_type == SURFACE) {
    if(boundaries.size() ==0) {
      cerr << "'extract -surf' must have one boundary surface identified using the '-bc' flag" << endl;
      Usage(ac,av);
    }
    if(iteration == "") {
      cerr << "'extract -surf' must specify iteration to extract from"
           << endl;
      Usage(ac,av);
    }
    get_surf(casename,iteration,variables,variable_type,variable_file,
             boundaries);
    Loci::Finalize();
    exit(0);
  }

  // New grid topology processor
  postProcessorP postprocessor = 0;

  switch(plot_type) {
  case ENSIGHT:
    postprocessor = new ensightPartConverter(id_required);
    break;
  case CGNS:
    postprocessor = new cgnsPartConverter(id_required);
    break;
  case FIELDVIEW:
    postprocessor = new fieldViewPartConverter;
    break;
  case TECPLOT:
    postprocessor = new tecplotPartConverter;
    break;
  case VTK:
    postprocessor = new vtkPartConverter(false);
    break;
  case VTK_SURFACE:
    postprocessor = new vtkSurfacePartConverter(false);
    break;
  case VTK64:
    postprocessor = new vtkPartConverter(true);
    break;
  case VTK_SURFACE64:
    postprocessor = new vtkSurfacePartConverter(true);
    break;
  case CUTTINGPLANE:
    postprocessor = new cuttingPlanePartConverter(transformMatrix, -xShift, -yShift, -zShift);
    break;
  default:
    cerr << "Unknown export method!" << endl;
    break;
  }

  if(postprocessor != 0) {
#ifdef H5_USE_16_API
    H5Eset_auto(NULL,NULL);
#else
    H5Eset_auto(H5E_DEFAULT,NULL,NULL);
#endif
    vector<surfacePartP> parts;


    // Check for derived variable input requirements
    std::set<string> varset;
    for(size_t i=0;i<variables.size();++i)
      varset.insert(variables[i]);
    if(varset.find("m") != varset.end()) {
      varset.insert("a");
      varset.insert("v");
    }
    if(varset.find("P") != varset.end()) {
      varset.insert("pg");
    }
    if(varset.find("p") != varset.end()) {
      varset.insert("pg");
    }
    if(varset.find("u") != varset.end()) {
      varset.insert("v");
    }
    if(varset.find("0") != varset.end()) {
      varset.insert("v");
    }
    if(varset.find("1") != varset.end()) {
      varset.insert("v");
    }
    if(varset.find("2") != varset.end()) {
      varset.insert("v");
    }
    vector<string> varlist;
    std::set<string>::const_iterator vi;
    for(vi=varset.begin();vi!=varset.end();++vi)
      varlist.push_back(*vi);
    variables.swap(varlist);

    if(postprocessor->processesSurfaceElements()) {
      std::set<string> volsearch;
      for(size_t i=0;i<partlist.size();++i) {
        string name = partlist[i];
        string dir  = output_dir + "/" + casename + "_SURF." + name;
        if(skip_part_dirs)
        {
          volsearch.insert(name);
        }else
        {
          surfacePartP sp = new surfacePart(name,dir,iteration,variables);
          if(sp->fail()) 
          {
            volsearch.insert(name);
          }else 
          {
            cout << "part: " << name << endl;
            parts.push_back(sp);
          }
        }
      }
      if(!volsearch.empty()) {
        vector<surfacePartP> volSurface;

        volumePartP vp =
          new volumePart(output_dir,iteration,casename,variables);
        if(!vp->fail())
          extractVolumeSurfaces(volSurface,vp,output_dir,iteration,
                                casename,variables);

        for(size_t i=0;i<volSurface.size();++i) {
          if(volsearch.find(volSurface[i]->getPartName()) != volsearch.end()) {
            cout << "part: " << volSurface[i]->getPartName() << endl;
            parts.push_back(volSurface[i]);
          }
        }
      }
    }
    volumePartP vp = 0;
    if(!novolume && postprocessor->processesVolumeElements()) {
      string testfile = getTopoFileName(output_dir, casename, iteration);
      struct stat tmpstat;
      cout << "checking " << testfile << endl;
      if(stat(testfile.c_str(),&tmpstat)==0) {
        // topo file exists, so there is a volume grid
        cout << "creating volume part" << endl;
        vp = new volumePart(output_dir,iteration,casename,variables);
        if(vp->fail()) {
          vp = 0;
        }
      }
    }
    particlePartP pp = 0;
    if(postprocessor->processesParticleElements()) {
      string testfile = output_dir + "/particle_pos."+iteration + "_" + casename;
      struct stat tmpstat;
      cout << "checking " << testfile << endl;
      if(stat(testfile.c_str(),&tmpstat)==0) {
        pp = new particlePart(output_dir,iteration,casename,variables,max_particles);
      }
    }

    if(parts.size() > 0) {
      vector<surfacePartP> modparts(parts.size());
      for(size_t i=0;i<parts.size();++i)
        modparts[i] = new surfacePartDerivedVars(parts[i],
                                                 output_dir,casename,
                                                 iteration, variables);
      postprocessor->addSurfaceParts(modparts);
    }

    if(vp!=0) {
      volumePartP vpn = new volumePartDerivedVars(vp,
                                                  output_dir,casename,
                                                  iteration,variables);
      postprocessor->addVolumePart(vpn);
    }
    if(pp!=0)
      postprocessor->addParticlePart(pp);
    postprocessor->exportPostProcessorFiles(casename,iteration);
    Loci::Finalize();
    exit(0);
  }


  if(timesyncerror) {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl;
  }
  Loci::Finalize();
} // End of main()


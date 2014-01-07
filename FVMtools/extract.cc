//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include <Loci.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;
using std::sort ;
using std::unique ;
using std::ifstream ;

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "extract.h"

string output_dir ;

    

void Usage(int ac, char *av[]) {
  cerr << av[0] << ": Incorrect Usage" << endl ;
  cout << "Usage:" << endl;
  cout << av[0] << " <package> [package options] <case_name> <time_step> <variable(s)>" << endl ;
  cout << endl ;
  cout << "where <package> may be:" << endl
       << "-2d :  extract for the 2dgv plotting package" << endl
       << "-fv :  extract for the FieldView post-processing package" << endl
       << "-en :  extract for the Ensight post-processing package" << endl
       << "-en_with_id :  extract for the Ensight post-processing package with node id and element id" << endl
       << "-tec:  extract for the TecPlot post-procesing package" << endl
       << "-vtk:   extract for the Paraview post-procesing package" << endl
       << "-vtk64: extract for the Paraview post-procesing package (for large cases, must use >= Paraview 3.98)" << endl
       << "-vtk_surf:  extract boundary surface mesh for the Paraview post-procesing package" << endl
       << "-vtk_surf64:  extract boundary surface mesh for the Paraview post-procesing package (for large cases, must use >= Paraview 3.98)" << endl
       << "-ascii: extract to an ascii file" << endl
       << "-surf: extract boundary surface mesh" << endl
       << "-cut:  extract a cutting plane for the 2dgv plotting package" << endl
       << "-mean: generate mean and variance from a family of ouput variables" << endl 
       << "-combine: combine mean and variance from online averaging ouput" << endl 
       << "-fcombine: combine favre mean and variance from online averaging ouput" << endl 
       << endl ;
  cout << "Variables are defined by the solver, but typically include: " << endl
       << "r     - nodal density" << endl 
       << "p     - nodal log10 pressure" << endl
       << "P     - nodal absolute pressure" << endl
       << "pg    - nodal gage pressure" << endl 
       << "u     - nodal velocity magnitude" << endl
       << "m     - nodal mach number" << endl
       << "t     - nodal temperature" << endl
       << "a     - nodal soundspeed" << endl
    //       << "G - extract Gamma" << endl
    //       << "R - extract Gas Constant R~" << endl
       << "f<species> - nodal mass fractions for <species>" << endl
       << "v     - nodal velocity vector"<<endl
       << "x     - x coordinate" << endl 
       << "y     - y coordinate" << endl 
       << "z     - z coordinate" << endl
       << "Boundary Variables:" << endl 
       << "qdot  - wall boundary heat flux" << endl
       << "yplus - wall boundary y plus" << endl
       << "tau   - viscous wall shear stress vector" << endl
       << "tw    - viscous wall boundary temperature" << endl
       << "pw    - viscous wall boundary pressure" << endl
       << "n     - boundary normal (-ascii only)" << endl
       << "area  - boundary area   (-ascii only)" << endl
       << endl ;

  cout << "extra options for particle extraction" << endl ;
  cout << "  -mp <n> : maximum particles to extract" << endl ;
  cout << "          : default (if not specified) is all particles available"
       << endl ;

  cout << "extra options for the 2dgv postprocessing package" << endl 
       << "  -bc <boundary_tag>  : specify which boundary to extract for (multiple -bc's"
       << endl 
       << "             are merged)" << endl 
       << "  -xy : project boundary onto z=0 plane (default)" << endl
       << "  -yz : project boundary onto x=0 plane" << endl
       << "  -xz : project boundary onto y=0 plane" << endl
       << "  -xr : project boundary onto x,radius plane (for axisymmetric grids)"
       << endl << endl ;
  cout << "extra options for orienting cutting plane" << endl
       << "  -xy : originate transformations from z=0 plane (default)" << endl
       << "  -yz : originate transformations from x=0 plane" << endl
       << "  -xz : originate transformations from y=0 plane" << endl
       << "  -Rx <amount> : rotate cutting plane about x-axis" << endl
       << "  -Ry <amount> : rotate cutting plane about y-axis" << endl
       << "  -Rz <amount> : rotate cutting plane about z-axis" << endl
       << "  -Sx <amount> : translate cutting plane along x-axis" << endl
       << "  -Sy <amount> : translate cutting plane along y-axis" << endl
       << "  -Sz <amount> : translate cutting plane along z-axis" << endl << endl;
  cout << "extra options for averaging feature '-mean'" << endl 
       << "  -end <value> : ending iteration number for averaging" << endl
       << "  -inc <value> : value to increment between iterations for averaging" << endl 
       << endl ;
  cout << "extra options for controlling extract" << endl
       << "   -dir <directory> : change extract directory from default 'output'"

       << endl << endl ;
  cout << "example:  to extract OH species from time step 50 of ssme simulation for visualization with 2dgv use:" << endl
       << av[0] << " -2d -bc 1 -xr ssme 50 fOH" << endl ;
  cout << "example: to extract an ascii table of boundary heat flux and x locations:"<<endl
       << av[0] << " -ascii -bc 4 nozzle 0 x qdot" << endl ;
  cout << "example: to extract a cutting plane with various transformations:" << endl
       << av[0] << " -cut -xz -Sy 1.5 -Rx 30 -Rz -15 nozzle 0 P t" << endl;
  cout << "example: to extract 5000 particles with associated temperature"
       << " from time step 100 for visualization with Ensight:" << endl
       << av[0] << " -en combustor 100 particle_temp -mp 5000" << endl ;

  cout << "example: to compute mean and variance values for velocity and temperature" 
       << " for iterations 1000-3000 output every 100 iterations run:"
       << endl
       << av[0] << " -mean -end 3000 -inc 100 nozzle 1000 t v" << endl
       << "NOTE: outputs variabes tMean, tVar, vMean, vVar, vCuv, vCuw, vCvw"
       << "      at iteration 1000." << endl ;
  exit(-1) ;
}

int  sizeElementType(hid_t group_id, const char *element_name) {
  hid_t dataset = H5Dopen(group_id,element_name) ;
  if(dataset < 0) {
    H5Eclear() ;
    return 0 ;
  }
  hid_t dspace = H5Dget_space(dataset) ;

  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  
  
  H5Dclose(dataset) ;
  return int(size) ;
  
}

string getPosFile(string output_dir,string iteration, string casename) {
  string posname = output_dir+"/grid_pos." + iteration + "_" + casename ;
  struct stat tmpstat ;
  if(stat(posname.c_str(),&tmpstat) != 0) {
    posname = output_dir+"/grid_pos." + casename ;
  } else if(tmpstat.st_size == 0) {
    posname = output_dir+"/grid_pos." + casename ;
  }
  return posname ;
}

void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) {
  if(var_name == "m") {
    string filename = output_dir+"/a_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> soundSpeed ;
    Loci::readContainer(file_id,"a",soundSpeed.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd])/soundSpeed[nd] ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "p" || var_name == "P") {
    string filename = output_dir+"/pg_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> pg ;
    Loci::readContainer(file_id,"pg",pg.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = output_dir+"/Pambient_par." + iteration +"_" + casename ;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    param<float> Pambient ;
    Loci::readContainer(file_id,"Pambient",Pambient.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = pg.domain() ;
    bool log = (var_name == "p") ;
    int c = 0 ;
    FORALL(dom,nd) {
      float p = pg[nd]+Pambient[nd] ;
      if(log)
        dval[c++] = log10(p) ;
      else
        dval[c++] = p ;
    } ENDFORALL ;
  } else if (var_name == "u") {
    fact_db facts ;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd]) ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "x" || var_name =="y" || var_name == "z") {
    store<vector3d<float> > pos ;
    string posname = getPosFile(output_dir,iteration,casename) ;
    hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to get grid positions for iteration " << iteration
           << endl ;
      cerr << "does file '" << posname << "' exist?" << endl ;
      Loci::Abort() ;
      exit(-1) ;
    }

    fact_db facts ;
    Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    entitySet dom = pos.domain() ;
    int c = 0 ;
    if(var_name == "x") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "y") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].y ;
      } ENDFORALL ;
    }
    if(var_name == "z") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].z ;
      } ENDFORALL ;
    }
  } else if(var_name == "0" || var_name =="1" || var_name == "2") {
    fact_db facts ;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    if(var_name == "0") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "1") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].y ;
      } ENDFORALL ;
    }
    if(var_name == "2") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].z ;
      } ENDFORALL ;
    }
  } else {
    cerr << "don't know how to get derived variable " << var_name << endl ;
  }
}

string getTopoFileName(string output_dir, string casename, string iteration) {
  string gridtopo = output_dir+"/" + casename +".topo" ;
  string toponamefile = output_dir + "/topo_file." + iteration + "_" + casename ;
  struct stat tmpstat ;
  if(stat(toponamefile.c_str(),&tmpstat)== 0) {
    ifstream tinput(toponamefile.c_str()) ;
    string name  ;
    tinput >> name ;
    name = output_dir + "/" + name ;
    if(stat(name.c_str(),&tmpstat)==0)
      gridtopo=name ;
  }
  return gridtopo ;
}

void setup_grid_topology(string casename, string iteration) {
  fact_db facts ;
  string file = casename + ".vog" ;
  struct stat tmpstat ;
  if(stat(file.c_str(),&tmpstat) != 0) {
    cerr << "unable to find vog file = " << file << endl ;
    Loci::Abort() ;
  }

  if(!Loci::setupFVMGrid(facts,file)) {
    cerr << "unable to read grid " << file << endl ;
  }
  // if output directory doesn't exist, create one
  struct stat statbuf ;

  if(stat(output_dir.c_str(),&statbuf))
    mkdir(output_dir.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file '"
           << output_dir << "' should be a directory!, rename and start again."
           << endl ;
      Loci::Abort() ;
    }
  
  string filename = output_dir+"/"+casename+".topo" ;
  if(stat(filename.c_str(),&tmpstat)!= 0) {
    Loci::createLowerUpper(facts) ;
    multiMap upper,lower,boundary_map,face2node ;
    Map ref ;
    store<string> boundary_names ;
    constraint geom_cells ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;
    face2node = facts.get_variable("face2node") ;
    ref = facts.get_variable("ref") ;
    boundary_names = facts.get_variable("boundary_names") ;
    geom_cells = facts.get_variable("geom_cells") ;
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;

    // If topology file does not exist, create it.
    Loci::parallelWriteGridTopology(filename.c_str(),
                                    upper.Rep(),lower.Rep(),boundary_map.Rep(),
                                    face2node.Rep(),
                                    ref.Rep(),
                                    boundary_names.Rep(),
                                    pos.Rep(),
                                    *geom_cells,
                                    facts) ;
  }

  store<vector3d<double> > pos ;
  pos = facts.get_variable("pos") ;
  filename = getPosFile(output_dir,iteration,casename) ;
  hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT) ;
  
  Loci::writeContainer(file_id,"pos",pos.Rep(),facts) ;
  
  Loci::hdf5CloseFile(file_id) ;

}


//read volume element ids only for ensight output
void extract_grid(string casename, string iteration,
                  grid_topo_handler *topo,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int max_particles, bool id_required) {
  
  
  vector<string> bnd_scalar_vars,bnd_scalar_filenames ;
  vector<string> bnd_vector_vars,bnd_vector_filenames ;
  vector<string> mfvars ;
  vector<string> pt_scalar_vars, pt_scalar_filenames ;
  vector<string> pt_vector_vars, pt_vector_filenames ;

  bool particle_extract = false ;
  {
    vector<string> vname ;
    vector<int> vtype;
    vector<string> vfile ;
  
    for(size_t i=0;i<variables.size();++i) {
      const string var_name(variables[i]) ;
      const int vt = variable_types[i] ;
      const string filename(variable_filenames[i]) ;
      switch(vt) {
      case NODAL_SCALAR:
      case NODAL_DERIVED:
      case NODAL_VECTOR:
        vname.push_back(var_name) ;
        vfile.push_back(filename) ;
        vtype.push_back(vt) ;
        break ;
      case NODAL_MASSFRACTION:
        mfvars.push_back(var_name) ;
        break ;
      case BOUNDARY_SCALAR:
        bnd_scalar_vars.push_back(var_name) ;
        bnd_scalar_filenames.push_back(filename) ;
        break ;
      case BOUNDARY_VECTOR:
        bnd_vector_vars.push_back(var_name) ;
        bnd_vector_filenames.push_back(filename) ;
        break ;
      case PARTICLE_SCALAR:
        pt_scalar_vars.push_back(var_name) ;
        pt_scalar_filenames.push_back(filename) ;
        particle_extract = true ;
        break ;
      case PARTICLE_VECTOR:
        pt_vector_vars.push_back(var_name) ;
        pt_vector_filenames.push_back(filename) ;
        particle_extract = true ;
        break ;
      default:
        cerr << "unable to process variable " << var_name << endl ;
        break ;
      }
    }
    for(size_t i=0;i<mfvars.size();++i) {
      vname.push_back(mfvars[i]) ;
      vfile.push_back(string("")) ;
      vtype.push_back(NODAL_MASSFRACTION) ;
    }
    for(size_t i=0;i<bnd_scalar_vars.size();++i) {
      vname.push_back(bnd_scalar_vars[i]) ;
      vfile.push_back(bnd_scalar_filenames[i]) ;
      vtype.push_back(BOUNDARY_SCALAR) ;
    }
    for(size_t i=0;i<bnd_vector_vars.size();++i) {
      vname.push_back(bnd_vector_vars[i]) ;
      vfile.push_back(bnd_vector_filenames[i]) ;
      vtype.push_back(BOUNDARY_VECTOR) ;
    }
    for(size_t i=0;i<pt_scalar_vars.size();++i) {
      vname.push_back(pt_scalar_vars[i]) ;
      vfile.push_back(pt_scalar_filenames[i]) ;
      vtype.push_back(PARTICLE_SCALAR) ;
    }
    for(size_t i=0;i<pt_vector_vars.size();++i) {
      vname.push_back(pt_vector_vars[i]) ;
      vfile.push_back(pt_vector_filenames[i]) ;
      vtype.push_back(PARTICLE_VECTOR) ;
    }

    variables = vname ;
    variable_filenames = vfile ;
    variable_types = vtype ;
  }

  Array<int,7> events ;
  topo->fileWritingSequence(events) ;
  FATAL(Loci::MPI_processes != 1) ;
  store<vector3d<float> > pos ;
  string posname = getPosFile(output_dir,iteration,casename) ;
  hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to get grid positions for iteration " << iteration
         << endl ;
    cerr << "does file '" << posname << "' exist?" << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }

  fact_db facts ;
  Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;
  int npnts = pos.domain().size() ;

  string iblankname = output_dir+"/grid_iblank." + iteration + "_" + casename ;
  store<unsigned char> iblank ;
  entitySet pdom = interval(1,npnts) ;
  iblank.allocate(pdom) ;
  struct stat tmpstat ;
  bool has_iblank = false ;
  if(stat(iblankname.c_str(),&tmpstat)== 0) {
    hid_t file_id = Loci::hdf5OpenFile(iblankname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to get iblank info for iteration " << iteration
           << endl ;
      cerr << "is file '" << iblankname << "' corrupted?" << endl ;
      Loci::Abort() ;
      exit(-1) ;
    }

    store<unsigned char> iblank_tmp ;
    Loci::readContainer(file_id,"iblank",iblank_tmp.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    entitySet dom = iblank_tmp.domain() ;
    int cnt = 1 ;
    FORALL(dom,nd) {
      iblank[cnt++] = iblank_tmp[nd] ;
    } ENDFORALL ;
    has_iblank = true ;
  } else {
    for(int i=1;i<=npnts;++i)
      iblank[i] = 0 ;
  }
  
  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;

  cout << "extracting topology from '" << gridtopo << "'" << endl;

  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;

  hid_t elg = H5Gopen(file_id,"elements") ;

  int ntets = sizeElementType(elg,"tetrahedra") ;
  int nhexs = sizeElementType(elg,"hexahedra") ;
  int nprsm = sizeElementType(elg,"prism") ;
  int npyrm = sizeElementType(elg,"pyramid") ;
  int ngenc = sizeElementType(elg,"GeneralCellNfaces") ;

  int ntets_b = ntets ;
  int nhexs_b = nhexs ;
  int nprsm_b = nprsm ;
  int npyrm_b = npyrm ;
  int ngenc_b = ngenc ;
  const int block_size=65536 ; // Size of blocking factor
  if(has_iblank) {
    // need to adjust the number of elements based on iblanking.
    if(ntets > 0) {
      int nblocks = (ntets-1)/block_size+1 ;
      int remain = ntets ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,4> > tets(size) ;
        readElementTypeBlock(elg,"tetrahedra",tets,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<4;++j)
            if(iblank[tets[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
        }
      }
      WARN(remain != 0) ;
      ntets_b = cnt ;
      if(ntets-ntets_b > 0)
        cout << ntets-ntets_b << " tetrahedra iblanked" << endl ;
    }
    if(nhexs > 0) {
      int nblocks = (nhexs-1)/block_size+1 ;
      int remain = nhexs ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        
        vector<Array<int,8> > hexs(size) ;
        readElementTypeBlock(elg,"hexahedra",hexs,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<8;++j)
            if(iblank[hexs[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
        }
      }
      WARN(remain != 0) ;
      nhexs_b = cnt ;
      if(nhexs-nhexs_b > 0)
        cout << nhexs-nhexs_b << " hexahedra iblanked" << endl ;
    }
    if(nprsm > 0) {
      int nblocks = (nprsm-1)/block_size+1 ;
      int remain = nprsm ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,6> > prsm(size) ;
        readElementTypeBlock(elg,"prism",prsm,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<6;++j)
            if(iblank[prsm[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
        }
      }
      WARN(remain != 0) ;
      nprsm_b = cnt ;
      if(nprsm-nprsm_b > 0)
        cout << nprsm-nprsm_b << " prisms iblanked" << endl ;
        
    }
    if(npyrm > 0) {
      int nblocks = (npyrm-1)/block_size+1 ;
      int remain = npyrm ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,5> > pyrm(size) ;
        readElementTypeBlock(elg,"pyramid",pyrm,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<5;++j)
            if(iblank[pyrm[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
        }
      }
      WARN(remain != 0) ;
      npyrm_b = cnt ;
      if(npyrm-npyrm_b > 0)
        cout << npyrm-npyrm_b << " pyramids iblanked" << endl ;
    }
    if(ngenc > 0) {
        vector<int> GeneralCellNfaces(ngenc) ;
        readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
        int nside = sizeElementType(elg,"GeneralCellNsides") ;
        vector<int> GeneralCellNsides(nside) ;
        readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
        int nnodes = sizeElementType(elg,"GeneralCellNodes") ;
        vector<int> GeneralCellNodes(nnodes) ;
        readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
        int cnt1 = 0 ;
        int cnt2 = 0 ;
        int cnt = 0 ;
        for(int i=0;i<ngenc;++i) {
          bool blank = true ;
          int nf = GeneralCellNfaces[i] ;
          for(int f=0;f<nf;++f) {
            int fs = GeneralCellNsides[cnt1++] ;
            for(int n=0;n<fs;++n) {
              int nd = GeneralCellNodes[cnt2++] ;
              if(iblank[nd] < 2)
                blank = false ;
            }
          }
          if(!blank)
            cnt++ ;
        }
        ngenc_b = cnt ;
        if(ngenc-ngenc_b > 0)
          cout << ngenc-ngenc_b << " general cells iblanked" << endl ;
    }
  }
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string>  bc_names ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
    bc_names.push_back(string(buf)) ;
  }
  
  cout << "npnts = " << npnts << ' '
       << "ntets = " << ntets_b << ' '
       << "npyrm = " << npyrm_b << ' '
       << "nprsm = " << nprsm_b << ' '
       << "nhexs = " << nhexs_b << ' '
       << "ngenc = " << ngenc_b << endl ;

  cout << "bcs = " ;
  for(size_t i=0;i<bc_names.size();++i)
    cout << bc_names[i] << ' ' ;
  cout << endl ;

  double time = 0 ;
  int ncycle = 0 ;
  string timefile = output_dir+"/timestep_txt." + iteration + "_" + casename ;
  std::ifstream timein(timefile.c_str(),ios::in) ;
  if(!timein.fail())
    timein >> ncycle >> time ;
  
  topo->open(casename,iteration,npnts,ntets_b,nprsm_b,npyrm_b,nhexs_b,ngenc_b,
             bc_names, variables,variable_types,time) ;


  for(int i=0;i<7;++i) {
    switch(events[i]) {
    case GRID_POSITIONS:
      {
        int minpos = pos.domain().Min() ;
        topo->create_mesh_positions(&pos[minpos],npnts) ;
        pos.allocate(EMPTY) ;
      }
      break ;
    case GRID_VOLUME_ELEMENTS:
      topo->create_mesh_elements() ;

      if(id_required && ntets > 0){
        int nblocks = (ntets-1)/block_size+1 ;
        int remain = ntets ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<int> tets_ids(size) ;
          vector<Array<int,4> > tets(size) ;
          readElementTypeBlock(elg,"tetrahedra_ids",tets_ids,start,size) ;
          readElementTypeBlock(elg,"tetrahedra",tets,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<4;++j)
                if(iblank[tets[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  {
                    tets[i-cnt]=tets[i] ;
                    tets_ids[i-cnt] = tets_ids[i];
                  }
              
            }
          }
          topo->write_tets_ids(&tets_ids[0],size-cnt,b,nblocks,ntets_b) ;
        }
      }


      if(ntets > 0) {
        int nblocks = (ntets-1)/block_size+1 ;
        int remain = ntets ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<Array<int,4> > tets(size) ;
          readElementTypeBlock(elg,"tetrahedra",tets,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<4;++j)
                if(iblank[tets[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  tets[i-cnt]=tets[i] ;
            }
          }
          topo->write_tets(&tets[0],size-cnt,b,nblocks,ntets_b) ;
        }
      }

      if(id_required && npyrm > 0){

        
        int nblocks = (npyrm-1)/block_size+1 ;
        int remain = npyrm ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<Array<int,5> > pyrm(size) ;
          vector<int > pyrm_ids(size) ;
          readElementTypeBlock(elg,"pyramid",pyrm,start,size) ;
          readElementTypeBlock(elg,"pyramid_ids",pyrm_ids,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<5;++j)
                if(iblank[pyrm[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  { pyrm[i-cnt]=pyrm[i] ;
                    pyrm_ids[i-cnt]=pyrm_ids[i] ;
                  }
            }
          }
          topo->write_pyrm_ids(&pyrm_ids[0],size-cnt,b,nblocks,npyrm_b) ;
        }

      }


        
      if(npyrm > 0) {
        int nblocks = (npyrm-1)/block_size+1 ;
        int remain = npyrm ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<Array<int,5> > pyrm(size) ;
          readElementTypeBlock(elg,"pyramid",pyrm,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<5;++j)
                if(iblank[pyrm[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  pyrm[i-cnt]=pyrm[i] ;
            }
          }
          topo->write_pyrm(&pyrm[0],size-cnt,b,nblocks,npyrm_b) ;
        }
      }


      if(id_required && nprsm > 0) {
        int nblocks = (nprsm-1)/block_size+1 ;
        int remain = nprsm ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<Array<int,6> > prsm(size) ;
          vector<int> prsm_ids(size);
          readElementTypeBlock(elg,"prism",prsm,start,size) ;
          readElementTypeBlock(elg,"prism_ids",prsm_ids,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<6;++j)
                if(iblank[prsm[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  {
                    prsm[i-cnt]=prsm[i] ;
                    prsm_ids[i-cnt]=prsm_ids[i] ;
                  }
            }
          }
          topo->write_prsm_ids(&prsm_ids[0],size-cnt,b,nblocks,nprsm_b) ;
          
        }
      }
      
      if(nprsm > 0) {
        int nblocks = (nprsm-1)/block_size+1 ;
        int remain = nprsm ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          vector<Array<int,6> > prsm(size) ;
          readElementTypeBlock(elg,"prism",prsm,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<6;++j)
                if(iblank[prsm[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  prsm[i-cnt]=prsm[i] ;
            }
          }
          topo->write_prsm(&prsm[0],size-cnt,b,nblocks,nprsm_b) ;
        }
      }
      if(id_required && nhexs > 0) {
        int nblocks = (nhexs-1)/block_size+1 ;
        int remain = nhexs ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          
          vector<Array<int,8> > hexs(size) ;
          vector<int> hexs_ids(size) ;
          readElementTypeBlock(elg,"hexahedra",hexs,start,size) ;
          readElementTypeBlock(elg,"hexahedra_ids",hexs_ids,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<8;++j)
                if(iblank[hexs[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  {
                    hexs[i-cnt]=hexs[i] ;
                     hexs_ids[i-cnt]=hexs_ids[i] ;
                  }
            }
          }
          topo->write_hexs_ids(&hexs_ids[0],size-cnt,b,nblocks,nhexs_b) ;
        }
      }
      
      if(nhexs > 0) {
        int nblocks = (nhexs-1)/block_size+1 ;
        int remain = nhexs ;
        int start = 0 ;
        for(int b=0;b<nblocks;++b) {
          int size = min(block_size,remain) ;
          
          vector<Array<int,8> > hexs(size) ;
          readElementTypeBlock(elg,"hexahedra",hexs,start,size) ;
          remain -= size ;
          start += size ;
          int cnt = 0 ;
          if(has_iblank) {
            for(int i=0;i<size;++i) {
              bool blank = true ;
              for(int j=0;j<8;++j)
                if(iblank[hexs[i][j]] < 2)
                  blank = false ;
              if(blank)
                cnt++ ;
              else 
                if(cnt != 0) // If there are some blanked copy into place
                  hexs[i-cnt]=hexs[i] ;
            }
          }
          topo->write_hexs(&hexs[0],size-cnt,b,nblocks,nhexs_b) ;
        }
      }
      if(id_required && ngenc > 0) {

        // still need to do general cell iblanking
        vector<int> GeneralCellNfaces(ngenc) ;
        readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
        int nside = sizeElementType(elg,"GeneralCellNsides") ;
        vector<int> GeneralCellNsides(nside) ;
        readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
        int nnodes = sizeElementType(elg,"GeneralCellNodes") ;
        vector<int> GeneralCellNodes(nnodes) ;
        readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
        vector<int> GeneralCell_ids(ngenc);
        readElementType(elg,"GeneralCell_ids",GeneralCell_ids) ;
        
        if(has_iblank) {
          int cnt1 = 0 ;
          int cnt2 = 0 ;
          int skip_cells = 0;
          int skip_faces = 0 ;
          int skip_nodes = 0 ;
          for(int i=0;i<ngenc;++i) {
            bool blank = true ;
            int nf = GeneralCellNfaces[i] ;
            int cnt1s = cnt1 ;
            int cnt2s = cnt2 ;
            for(int f=0;f<nf;++f) {
              int fs = GeneralCellNsides[cnt1++] ;
              for(int n=0;n<fs;++n) {
                int nd = GeneralCellNodes[cnt2++] ;
                if(iblank[nd] < 2)
                  blank = false ;
              }
            }
            if(blank) {
              skip_cells += 1 ;
              skip_faces += cnt1-cnt1s ;
              skip_nodes += cnt2-cnt2s ;
            } else {
              if(skip_cells > 0) {
                GeneralCellNfaces[i-skip_cells] = GeneralCellNfaces[i] ;
                GeneralCell_ids[i-skip_cells] = GeneralCell_ids[i] ;
                for(int j=0;j<cnt1-cnt1s;++j)
                  GeneralCellNsides[cnt1s+j-skip_faces] =
                    GeneralCellNsides[cnt1s+j] ;
                for(int j=0;j<cnt2-cnt2s;++j)
                  GeneralCellNodes[cnt2s+j-skip_nodes] =
                    GeneralCellNodes[cnt2s+j] ;
              }
            }
            
          }
          ngenc -= skip_cells ;
          nside -= skip_faces ;
          nnodes -= skip_nodes ;
        }
        topo->write_general_cell_ids(&GeneralCell_ids[0],ngenc) ;
      }

      if(ngenc > 0) {

        // still need to do general cell iblanking
        vector<int> GeneralCellNfaces(ngenc) ;
        readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
        int nside = sizeElementType(elg,"GeneralCellNsides") ;
        vector<int> GeneralCellNsides(nside) ;
        readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
        int nnodes = sizeElementType(elg,"GeneralCellNodes") ;
        vector<int> GeneralCellNodes(nnodes) ;
        readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
        if(has_iblank) {
          int cnt1 = 0 ;
          int cnt2 = 0 ;
          int skip_cells = 0;
          int skip_faces = 0 ;
          int skip_nodes = 0 ;
          for(int i=0;i<ngenc;++i) {
            bool blank = true ;
            int nf = GeneralCellNfaces[i] ;
            int cnt1s = cnt1 ;
            int cnt2s = cnt2 ;
            for(int f=0;f<nf;++f) {
              int fs = GeneralCellNsides[cnt1++] ;
              for(int n=0;n<fs;++n) {
                int nd = GeneralCellNodes[cnt2++] ;
                if(iblank[nd] < 2)
                  blank = false ;
              }
            }
            if(blank) {
              skip_cells += 1 ;
              skip_faces += cnt1-cnt1s ;
              skip_nodes += cnt2-cnt2s ;
            } else {
              if(skip_cells > 0) {
                GeneralCellNfaces[i-skip_cells] = GeneralCellNfaces[i] ;
                for(int j=0;j<cnt1-cnt1s;++j)
                  GeneralCellNsides[cnt1s+j-skip_faces] =
                    GeneralCellNsides[cnt1s+j] ;
                for(int j=0;j<cnt2-cnt2s;++j)
                  GeneralCellNodes[cnt2s+j-skip_nodes] =
                    GeneralCellNodes[cnt2s+j] ;
              }
            }
            
          }
          ngenc -= skip_cells ;
          nside -= skip_faces ;
          nnodes -= skip_nodes ;
        }
        topo->write_general_cell(&GeneralCellNfaces[0],ngenc,
                                 &GeneralCellNsides[0],nside,
                                 &GeneralCellNodes[0],nnodes) ;
      }
      topo->close_mesh_elements() ;
      H5Gclose(elg) ;
      break ;
    case GRID_BOUNDARY_ELEMENTS:
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str()) ;
        
        int nquads = sizeElementType(bcg,"quads") ;
        int ntrias = sizeElementType(bcg,"triangles") ;
        int ngeneral = sizeElementType(bcg,"nside_sizes") ;
        
        vector<Array<int,3> > trias(ntrias) ;
        readElementType(bcg,"triangles",trias) ;
        vector<Array<int,4> > quads(nquads) ;
        readElementType(bcg,"quads",quads) ;

        vector<int> nside_sizes(ngeneral) ;
        readElementType(bcg,"nside_sizes",nside_sizes) ;
        int nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
        vector<int> nside_nodes(nside_nodes_size) ;
        readElementType(bcg,"nside_nodes",nside_nodes) ;
        
        vector<int> node_set ;
        for(int i=0;i<nside_nodes_size;++i)
          node_set.push_back(nside_nodes[i]) ;
        
        for(int i=0;i<ntrias;++i) {
          node_set.push_back(trias[i][0]) ;
          node_set.push_back(trias[i][1]) ;
          node_set.push_back(trias[i][2]) ;
        }
        for(int i=0;i<nquads;++i) {
          node_set.push_back(quads[i][0]) ;
          node_set.push_back(quads[i][1]) ;
          node_set.push_back(quads[i][2]) ;
          node_set.push_back(quads[i][3]) ;
        }
        sort(node_set.begin(),node_set.end()) ;
        node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;
        
        map<int,int> nmap ;
        for(size_t i=0;i<node_set.size();++i) {
          nmap[node_set[i]] = i+1 ;
        }
        
        
        topo->create_boundary_part(bc_names[bc],&node_set[0],node_set.size()) ;
        
        for(int i=0;i<ntrias;++i) 
          for(int j=0;j<3;++j) 
            trias[i][j] = nmap[trias[i][j]] ;
        for(int i=0;i<nquads;++i) 
          for(int j=0;j<4;++j) 
            quads[i][j] = nmap[quads[i][j]] ;
        for(int i=0;i<nside_nodes_size;++i)
          nside_nodes[i] = nmap[nside_nodes[i]] ;
        
        vector<int > trias_id(ntrias) ;
        readElementType(bcg,"triangles_id",trias_id) ;
        vector<int > quads_id(nquads) ;
        readElementType(bcg,"quads_id",quads_id) ;
        vector<int > nside_id(ngeneral) ;
        readElementType(bcg,"nside_id",nside_id) ;
        
        int cnt = 0 ;
        for(int i=0;i<ntrias;++i) {
          bool blank = true ;
          for(int j=0;j<3;++j)
            if(iblank[node_set[trias[i][j]-1]] < 2)
              blank = false ;
          if(blank)
            cnt++ ;
          else 
            if(cnt != 0) { // If there are some blanked copy into place
              trias[i-cnt]=trias[i] ;
              trias_id[i-cnt] = trias_id[i] ;
            }
        }
        ntrias -=cnt ;
        

        topo->write_trias(&trias[0],&trias_id[0],ntrias) ;

        cnt = 0 ;
        for(int i=0;i<nquads;++i) {
          bool blank = true ;
          for(int j=0;j<4;++j)
            if(iblank[node_set[quads[i][j]-1]] < 2)
              blank = false ;
          if(blank)
            cnt++ ;
          else 
            if(cnt != 0) { // If there are some blanked copy into place
              quads[i-cnt]=quads[i] ;
              quads_id[i-cnt] = quads_id[i] ;
            }
        }
        nquads -=cnt ;
        


        topo->write_quads(&quads[0],&quads_id[0],nquads) ;
        topo->write_general_face(&nside_sizes[0], &nside_id[0], ngeneral,
                                 &nside_nodes[0], nside_nodes_size) ;
        H5Gclose(bcg) ;
        topo->close_boundary_part() ;
        
      }
      break ;
    case NODAL_VARIABLES:
      {
        
        for(size_t i=0;i<variables.size();++i) {
          const string var_name(variables[i]) ;
          const int vt = variable_types[i] ;
          const string filename(variable_filenames[i]) ;
          switch(vt) {
          case NODAL_SCALAR:
            {
              hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                                 H5F_ACC_RDONLY,
                                                 H5P_DEFAULT) ;
              if(file_id < 0) {
                cerr << "unable to open file '" << filename << "'!" << endl ;
                continue ;
              }

              fact_db facts ;
              store<float> scalar ;
              Loci::readContainer(file_id,var_name,scalar.Rep(),EMPTY,facts) ;
              entitySet dom = scalar.domain() ;
              int sz = dom.size() ;
              int min_val= dom.Min() ;
              topo->output_nodal_scalar(&scalar[min_val],sz,var_name) ;
              Loci::hdf5CloseFile(file_id) ;
            }
            break;
          case NODAL_DERIVED:
            {
              vector<float> dval(npnts) ;
              getDerivedVar(dval,var_name,casename,iteration) ;
              int sz = dval.size() ;
              topo->output_nodal_scalar(&dval[0],sz,var_name) ;
            }
            break;
          case NODAL_VECTOR:
            {
              hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                                 H5F_ACC_RDONLY,
                                                 H5P_DEFAULT) ;
              if(file_id < 0) {
                cerr << "unable to open file '" << filename << "'!" << endl ;
                continue ;
              }
              
              fact_db facts ;
              store<vector3d<float> > vec ;
              Loci::readContainer(file_id,var_name,vec.Rep(),EMPTY,facts) ;
              entitySet dom = vec.domain() ;
              int sz = dom.size() ;
              int min_val= dom.Min() ;
              topo->output_nodal_vector(&vec[min_val],sz,var_name) ;
              Loci::hdf5CloseFile(file_id) ;
            }
            break;
          case NODAL_MASSFRACTION:
            break ;
          case BOUNDARY_SCALAR:
          case BOUNDARY_VECTOR:
            break ;
          case PARTICLE_SCALAR:
          case PARTICLE_VECTOR:
            break ;
          default:
            cerr << "unable to process variable " << var_name << endl ;
            break ;
          }
        }

        if(mfvars.size() > 0) {
          string filename = output_dir+"/mix." + iteration + "_" + casename ;
          
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            Loci::Abort() ;
            exit(-1) ;
          }
          
          fact_db facts ;
          storeVec<float> mix ;
          Loci::readContainer(file_id,"mixture",mix.Rep(),EMPTY,facts) ;
          param<string> species_names ;
          Loci::readContainer(file_id,"species_names",species_names.Rep(),EMPTY,facts) ;
          Loci::hdf5CloseFile(file_id) ;
          
          map<string,int> smap ;
          std::istringstream iss(*species_names) ;
          for(int i=0;i<mix.vecSize();++i) {
            string s ;
            iss >> s ;
            smap[s] = i ;
          }
      
          entitySet dom = mix.domain() ;
          
          vector<float> vec(npnts) ;
          
          for(size_t i=0;i<mfvars.size();++i) {
            const string var_name(mfvars[i]) ;
            string sp = string(mfvars[i].c_str()+1) ;
            map<string,int>::const_iterator mi = smap.find(sp) ;
            if(mi == smap.end()) {
              cerr << "warning, species " << sp << " does not exist in dataset!"
                   << endl ;
            } else {
              const int ind = mi->second ;
              int c = 0 ;
              FORALL(dom,nd) {
                vec[c++] = mix[nd][ind] ;
              } ENDFORALL ;
              topo->output_nodal_scalar(&vec[0],npnts,var_name) ;
            }
          }

        }
      }
      break;
    case BOUNDARY_VARIABLES:
      if(bnd_scalar_vars.size() > 0) {
        for(size_t b=0;b<bnd_scalar_vars.size();++b) {
          const string filename(bnd_scalar_filenames[b]) ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }
          
          hid_t di = H5Gopen(file_id,"dataInfo") ;
          int nbel = sizeElementType(di,"entityIds") ;
          
          vector<int> elemIds(nbel) ;
          readElementType(di,"entityIds",elemIds) ;
          
          H5Gclose(di) ;
          vector<float> var(nbel) ;
          readElementType(file_id,bnd_scalar_vars[b].c_str(),var) ;
          
          topo->output_boundary_scalar(&var[0],&elemIds[0],nbel,
                                       bnd_scalar_vars[b]) ;
        }
      }
      
      if(bnd_vector_vars.size() > 0) {
        for(size_t b=0;b<bnd_vector_vars.size();++b) {
          const string filename(bnd_vector_filenames[b]) ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }
          
          hid_t di = H5Gopen(file_id,"dataInfo") ;
          int nbel = sizeElementType(di,"entityIds") ;
          
          vector<int> elemIds(nbel) ;
          readElementType(di,"entityIds",elemIds) ;
          
          H5Gclose(di) ;
          vector<vector3d<float> > var(nbel) ;
          readElementType(file_id,bnd_vector_vars[b].c_str(),var) ;
          
          topo->output_boundary_vector(&var[0],&elemIds[0],nbel,
                                       bnd_vector_vars[b]) ;
        }
      }
      break;
    case PARTICLE_POSITIONS:
      if(particle_extract) {
        string posname = output_dir+"/particle_pos." + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                           H5F_ACC_RDONLY,
                                           H5P_DEFAULT) ;
        if(file_id < 0) {
          cerr << "unable to get particle positions for iteration "
               << iteration << endl ;
          Loci::Abort() ;
          exit(-1) ;
        }
        
        int np = sizeElementType(file_id, "particle position") ;
        vector<vector3d<float> > ppos(np) ;
        readElementType(file_id, "particle position", ppos) ;
        
        topo->create_particle_positions(&ppos[0], np, max_particles) ;
      }
      break ;
    case PARTICLE_VARIABLES:
      for(size_t v=0;v<pt_scalar_vars.size();++v) {
        const string varname(pt_scalar_vars[v]) ;
        const string filename(pt_scalar_filenames[v]) ;
        hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                           H5F_ACC_RDONLY, H5P_DEFAULT) ;
        if(file_id < 0) {
          cerr << "unable to open file '" << filename << "'!" << endl ;
          continue ;
        }
        int np = sizeElementType(file_id, varname.c_str()) ;
        vector<float> scalar(np) ;
        readElementType(file_id, varname.c_str(), scalar) ;

        topo->output_particle_scalar(&scalar[0], np, max_particles, varname) ;
        
      }

      for(size_t v=0;v<pt_vector_vars.size();++v) {
        const string varname(pt_vector_vars[v]) ;
        const string filename(pt_vector_filenames[v]) ;
        hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                           H5F_ACC_RDONLY, H5P_DEFAULT) ;

        if(file_id < 0) {
          cerr << "unable to open file '" << filename << "'!" << endl ;
          continue ;
        }
        int np = sizeElementType(file_id, varname.c_str()) ;
        vector<vector3d<float> > vec(np) ;
        readElementType(file_id, varname.c_str(), vec) ;

        topo->output_particle_vector(&vec[0], np, max_particles, varname) ;
      }
      
      break ;
    default:
      cerr << "internal error, problem with topo sequence" << endl ;
    }
  }
  topo->close() ;
  H5Gclose(bndg) ;
  H5Fclose(file_id) ;
}
namespace Loci {
  void disableDebugDir() ;
}
int main(int ac, char *av[]) {
  output_dir = "output" ;
  Loci::disableDebugDir() ;
  Loci::Init(&ac,&av) ;

  enum {ASCII,TWODGV,ENSIGHT,FIELDVIEW,TECPLOT,VTK,VTK_SURFACE,VTK64,VTK_SURFACE64,CUTTINGPLANE, SURFACE, MEAN, COMBINE, FCOMBINE, NONE} plot_type = NONE ;

  string casename ;
  bool found_casename = false ;
  bool found_iteration = false ;
  string iteration ;
  vector<string> variables ;
  vector<string> boundaries ;
  float xShift, yShift, zShift, temp;

  xShift = yShift = zShift = 0.0;
  int view = VIEWXY ;
  affineMapping transformMatrix;

  // record the maximum particle number to extract
  // a value < 0 means that there is no maximum particle
  // number limit, i.e., all particles are to be extracted
  // default is to extract all particles.  users can use
  // command line switch "-mp <n>" to set the maximum number
  // of particles to be extracted. if the requested particle
  // number is larger than the available particle number, then
  // all particles will be extracted.
  int max_particles = -1 ;

  int end_iter = -1 ;
  int inc_iter = -1 ;
  bool id_required = false;//ensight has the option to display node and element ids

  
  for(int i=1;i<ac;++i) {
    if(av[i][0] == '-') {
      if(!strcmp(av[i],"-ascii"))
        plot_type = ASCII ;
      else if(!strcmp(av[i],"-surf"))
        plot_type = SURFACE ;
      else if(!strcmp(av[i],"-mean"))
        plot_type = MEAN ;
      else if(!strcmp(av[i],"-combine"))
        plot_type = COMBINE ;
      else if(!strcmp(av[i],"-fcombine"))
        plot_type = FCOMBINE ;
      else if(!strcmp(av[i],"-2d"))
        plot_type = TWODGV ;
      else if(!strcmp(av[i],"-en"))
        plot_type = ENSIGHT ;
      else if(!strcmp(av[i],"-en_with_id"))
        {
          plot_type = ENSIGHT ;
          id_required = true;
        }
      else if(!strcmp(av[i],"-fv"))
        plot_type = FIELDVIEW ;
      else if(!strcmp(av[i],"-tec")) {
        plot_type = TECPLOT ;
	cout << "*****************************************************************************"<< endl ;
	cout << "NOTICE!! Latest versions of tecplot360 will perform better using using the" << endl
	     << "         the Ensight importer.  It is recommended that you use extract -en " << endl
	     << "         instead of extract -tec." << endl ;
	cout << "*****************************************************************************"<< endl ;
      } 
      else if(!strcmp(av[i],"-vtk"))
        plot_type = VTK ;
      else if(!strcmp(av[i],"-vtk_surf"))
        plot_type = VTK_SURFACE ;
      else if(!strcmp(av[i],"-vtk64"))
        plot_type = VTK64 ;
      else if(!strcmp(av[i],"-vtk_surf64"))
        plot_type = VTK_SURFACE64 ;
      else if(!strcmp(av[i],"-cut"))
	plot_type = CUTTINGPLANE ;
      else if(!strcmp(av[i],"-Sx")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> xShift).fail())
	    Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Sy")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> yShift).fail())
	    Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Sz")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> zShift).fail())
	    Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Rx")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
	    Usage(ac, av) ;
        transformMatrix.rotateX(-temp) ;
      }
      else if(!strcmp(av[i],"-Ry")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
	    Usage(ac, av) ;
	transformMatrix.rotateY(-temp) ;
      }
      else if(!strcmp(av[i],"-Rz")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
	    Usage(ac, av) ;
	transformMatrix.rotateZ(-temp) ;
      }
      else if(!strcmp(av[i],"-xy")) 
        view=VIEWXY ;
      else if(!strcmp(av[i],"-yz")) {
        view=VIEWYZ ;
	transformMatrix.rotateY(90.0) ;
	transformMatrix.rotateZ(90.0) ;
      }
      else if(!strcmp(av[i],"-xz")) {
        view=VIEWXZ ;
	transformMatrix.rotateX(90.0) ;
      }
      else if(!strcmp(av[i],"-xr")) 
        view=VIEWXR ;
      else if(!strcmp(av[i],"-bc")) {
        i++ ;
        string v(av[i]) ;
        if(av[i][0] >='0' && av[i][0] <= '9')
          v = "BC_"+v ;
        boundaries.push_back(v) ;
      } else if(!strcmp(av[i],"-mp")) {
        // get the number of particles
        ++i ;
        string n(av[i]) ;
        if(!valid_int(n)) {
          cerr << "argument followed option '-mp' is not an integer,"
               << " used default value" << endl ;
        } else {
          max_particles = str2int(n) ;
        }
        
      } else if(!strcmp(av[i],"-dir")) {
        ++i ;
        output_dir = string(av[i]) ;
      } else if(!strcmp(av[i],"-inc")) {
        ++i ;
        inc_iter = atoi(av[i]) ;
      } else if(!strcmp(av[i],"-end")) {
        ++i ;
        end_iter = atoi(av[i]) ;
      } else {
        cerr << "unknown option " << av[i] << endl ;
        Usage(ac,av) ;
      }
      
    } else {
      if(found_iteration)
        variables.push_back(string(av[i])) ;
      else if(found_casename) {
        iteration = string(av[i]) ;
        found_iteration = true ;
      } else {
        casename = string(av[i]) ;
        found_casename = true ;
      }
    }
  }
  if(plot_type == NONE) {
    Usage(ac,av) ;
  }

  // if output directory doesn't exist, create one
  struct stat statbuf ;
  if(stat(output_dir.c_str(),&statbuf))
    mkdir(output_dir.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file '" << output_dir
           <<"' should be a directory!, rename 'output' and start again."
           << endl ;
      Loci::Abort() ;
    }
  
      

  if(variables.size() == 0) {
    DIR *dp = opendir(output_dir.c_str()) ;
    // Look in output directory and find all variables
    if(dp == 0) {
      cerr << "unable to open directory '" << output_dir << "'" << endl ;
      exit(-1) ;
    }
    dirent *entry = readdir(dp) ;
  
    string tail = iteration + "_" + casename ;
    while(entry != 0) {
      string filename = entry->d_name ;
      string postfix ;
      string vname ;
      string vtype ;
      int nsz = filename.size() ;
      int dot = -1 ;
      for(int i=nsz-1;i>=0;--i)
        if(filename[i] == '.') {
          dot = i ;
          break ;
        }
      for(int i=dot+1;i<nsz;++i)
        postfix += filename[i] ;
      int und = -1 ;
      if(dot > 0) {
        for(int i=dot-1;i>=0;--i)
          if(filename[i] == '_') {
            und = i ;
            break ;
          }
        if(und > 0) {
          for(int i=und+1;i<dot;++i)
            vtype += filename[i] ;
          for(int i=0;i<und;++i)
            vname += filename[i] ;
        } else
          for(int i=0;i<dot;++i)
            vname += filename[i] ;
        
        
      }
        
      if(dot>0 && und>0 && postfix == tail) {
        if(vtype == "sca" || vtype == "vec" || vtype == "bnd" ||
           vtype == "bndvec" || vtype == "ptsca" || vtype == "ptvec") {
          variables.push_back(vname) ;
        }
        if(vtype == "sca" && vname == "pg") {
          variables.push_back("P") ;
        }
        if(vtype == "sca" && vname == "a") {
          variables.push_back("m") ;
        }
      }
      if(dot > 0 && postfix == tail) {
        if(vname == "mix" && vtype == "") {
          string mfile = output_dir+"/mix." + iteration + "_" + casename ;
          
          hid_t file_id = Loci::hdf5OpenFile(mfile.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << mfile << "'!" << endl ;
            Loci::Abort() ;

            exit(-1) ;
          }
          fact_db facts ;
          
          param<string> species_names ;
          Loci::readContainer(file_id,"species_names",species_names.Rep(),EMPTY,facts) ;
          Loci::hdf5CloseFile(file_id) ;
          std::istringstream iss(*species_names) ;
          string s ;
          do {
            s = "" ;
            iss >> s ;
            if(s != "") {
              string v = "f"+s ;
              variables.push_back(v) ;
            }
          } while(!iss.eof() && s!= "") ;
            
        }
      }
      
       
      entry = readdir(dp) ;
    }
    if(variables.size() == 0) {
      variables.push_back("x") ;
      variables.push_back("y") ;
      variables.push_back("z") ;
    }
    cout << "extracting variables:" ;
    for(size_t i=0;i<variables.size();++i) {
      cout << ' ' << variables[i] ;
    }
    cout << endl ;
    closedir(dp) ;
      
  }
  
  vector<int> variable_type(variables.size()) ;
  vector<string> variable_file(variables.size()) ;

  bool particle_info_requested = false ;
  
  for(size_t i=0;i<variables.size();++i) {
    const string var(variables[i]) ;
    string filename = output_dir+'/' + var + "_hdf5." + iteration ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
      
    filename = output_dir+'/' + var + "_sca." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }

    filename = output_dir+'/' + var + "_vec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = output_dir+'/' + var + "_bnd." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = output_dir+'/' + var + "_bndvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }

    filename = output_dir+'/' + var + "_ptsca." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)==0) {
      variable_type[i] = PARTICLE_SCALAR ;
      variable_file[i] = filename ;
      particle_info_requested = true ;
      continue ;
    }

    filename = output_dir+'/' + var + "_ptvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)==0) {
      variable_type[i] = PARTICLE_VECTOR ;
      variable_file[i] = filename ;
      particle_info_requested = true ;
      continue ;
    }

    if(plot_type == ASCII && boundaries.size() > 0) {
      if(var == "n") {
        variable_type[i] = BOUNDARY_DERIVED_VECTOR ;
        continue ;
      }
      if(var == "area" || var == "x" || var == "y" || var == "z") {
        variable_type[i] = BOUNDARY_DERIVED_SCALAR ;
        continue ;
      }
    }
    if(var == "m") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "p") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "P") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "u") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "0" || var == "1" || var == "2") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "x" || var == "y" || var == "z") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
      
    if(var.size()>1 && var[0] == 'f') {
      variable_type[i] = NODAL_MASSFRACTION ;
      continue ;
    }

    // ... other derived variables here

    cerr << "Warning, variable '" << var << "' is unknown and will not be processed." << endl ;
    variable_type[i] = UNDEFINED ;
  }

  H5Eset_auto(NULL,NULL) ;

  // we will first check to see if particle position is present
  // in case of any particle information extraction
  if(particle_info_requested) {
    string filename = output_dir +"/particle_pos." + iteration + "_" + casename ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)!=0) {
      cerr << "Warning: particle geometry '" << filename << "' must"
           << " be presented for any"
           << " particle related variable extraction." << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
  }

  string filename = getTopoFileName(output_dir, casename, iteration) ;
  struct stat tmpstat ;
  string posfile = getPosFile(output_dir,iteration,casename) ;
  //  cout << "posfile = " << posfile << endl ;
  if(stat(filename.c_str(),&tmpstat)!= 0 ||
     stat(posfile.c_str(),&tmpstat) != 0) {
    cerr << "Warning, no grid topology information.  Will attempt to generate!"
         << endl ;
    setup_grid_topology(casename,iteration) ;
  }
  bool timesyncerror = false ;
  if(stat(filename.c_str(),&tmpstat)==0) {
    struct stat gridstat ;
    string gridfile = casename + ".vog" ;
    if(stat(gridfile.c_str(),&gridstat)==0) {
      if(gridstat.st_mtime > tmpstat.st_mtime)
        timesyncerror = true ;
    }
  }

  if(timesyncerror) {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl ;
  }

  
  if(plot_type == ASCII) {
    if(boundaries.size() == 0) {
      process_ascii_nodal(casename,iteration,
                          variables,variable_type,variable_file) ;
      Loci::Finalize() ;
      exit(0) ;
    } else {
      process_ascii_bndry(casename,iteration,
                          variables,variable_type,variable_file,
                          boundaries) ;
      Loci::Finalize() ;
      exit(0) ;
    }
  }

  if(plot_type == MEAN) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    process_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == COMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,false) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == FCOMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,true) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == TWODGV) {
    if(variables.size() != 1) {
      cerr << "2dgv extract can only extract one variable at a time."
           << endl ;
      Usage(ac,av) ;
    }
    if(boundaries.size() == 0) {
      cerr << "2dgv extract must have the projected boundaries identified using the '-bc' flag" << endl ;
      Usage(ac,av) ;
    }
    get_2dgv(casename,iteration,variables,variable_type,variable_file,
             boundaries,view) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == SURFACE) {
    if(boundaries.size() ==0) {
      cerr << "'extract -surf' must have one boundary surface identified using the '-bc' flag" << endl ;
      Usage(ac,av) ;
    }
    if(iteration == "") {
      cerr << "'extract -surf' must specify iteration to extract from"
           << endl ;
      Usage(ac,av) ;
    }
    get_surf(casename,iteration,variables,variable_type,variable_file,
             boundaries) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  
  // process grid topology
  grid_topo_handler *topo_out = 0 ;
 
  switch(plot_type) {
  case ENSIGHT:
    topo_out = new ensight_topo_handler( id_required);
    break ;
  case FIELDVIEW:
    topo_out = new fv_topo_handler ;
    break ;
  case TECPLOT:
    topo_out = new tecplot_topo_handler ;
    break ;
  case VTK:
    topo_out = new vtk_topo_handler(false) ;
    break ;
  case VTK_SURFACE:
    topo_out = new vtk_surf_topo_handler(boundaries,false) ;
    break ;
  case VTK64:
    topo_out = new vtk_topo_handler(true) ;
    break ;
  case VTK_SURFACE64:
    topo_out = new vtk_surf_topo_handler(boundaries,true) ;
    break ;
  case CUTTINGPLANE:
    topo_out = new cuttingplane_topo_handler(transformMatrix, -xShift, -yShift, -zShift) ;
    break ;
  default:
    cerr << "Unknown export method!" << endl ;
    break ;
  }

  if(topo_out != 0) {
    extract_grid(casename,iteration,
                 topo_out,variables,
                 variable_type,variable_file,max_particles, id_required) ;//read volume element  only for ensight output  
  }
  
  if(timesyncerror) {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl ;
  }
  Loci::Finalize() ;
}


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
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;

#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>

using namespace std ;
int TECINI(const char *problem_name,const char *variables_name,
           const char *tec_name,const char *dot,int *Debug,int *VIsDouble) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  cout << "Initializing file " << filename << endl ;
 
  ofstream ofile(filename.c_str(),ios::out) ;
  ofile << "variables = " << variables_name << endl ;
  ofile.close() ;

  return(0) ;
}

int TECZNE(const char *zone_name, int *npnts, int *ncells,const char *tec_name,
           const char *elem_type) {

  string filename = tec_name ;
  ofstream ofile(filename.c_str(),ios::app) ;
 
  ofile << "ZONE, T = \"" << zone_name << '"'
        << ", DATAPACKING=BLOCK" 
        << ", n = " << *npnts
        << ", e = " << *ncells
        << ", ZONETYPE = " << elem_type << endl ;
  ofile.close() ;

  return(0) ;
}

int TECDAT(int *npnts, float *vals, const char *tec_name) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  ofile.precision(16) ;

  for(int i = 0; i < *npnts; i++)
    ofile << vals[i] << endl ;

  ofile.close() ;

  return(0) ;
}

int TECNOD(int *nm, int *ncells, const char *tec_name) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  for(int i = 0; i < *ncells; i++) {
    int j = i*8 ;
    ofile << nm[j+0] << " " << nm[j+1] << " "
          << nm[j+2] << " " << nm[j+3] << " "
          << nm[j+4] << " " << nm[j+5] << " "
          << nm[j+6] << " " << nm[j+7] << endl ;
  }

  ofile.close() ;

  return(0) ;
}

int TECEND() {
  return(0)  ;
}

void tecplot_topo_handler::open(string casename, string iteration ,size_t inpnts,
                                size_t intets, size_t inprsm, size_t inpyrm,
                                size_t inhexs, size_t ingen,
                                const vector<string> &bc_names,
                                const vector<string> &variables,
                                const vector<int> &variable_types,
                                double time) {
  npnts = inpnts ;
  ntets = intets ;
  nprsm = inprsm ;
  npyrm = inpyrm ;
  nhexs = inhexs ;
  ngen = ingen ;
  nvars = 3 ;

  filename = "tec_"+casename+"_"+iteration+".dat" ;
  int ncells = ntets+nprsm+npyrm+nhexs ;
  string varstring = "\"x\", \"y\", \"z\"" ;
  for(size_t i=0;i<variables.size();++i) {
    if(variable_types[i] == NODAL_SCALAR ||
       variable_types[i] == NODAL_DERIVED ||
       variable_types[i] == NODAL_MASSFRACTION) {
      varstring += ", \"" + variables[i] + '"' ;
      nvars++ ;
    } else if(variable_types[i] == NODAL_VECTOR) {
      varstring += ", \"" + variables[i] +".x\"" ;
      varstring += ", \"" + variables[i] +".y\"" ;
      varstring += ", \"" + variables[i] +".z\"" ;
      nvars += 3;
    } else {
      cerr << "boundary variable " << variables[i] << "not currently supported by tecplot extractor" << endl ;
    }
  }
  int Debug = 1 ;
  int VIsDouble = 0 ;
  TECINI(casename.c_str(),varstring.c_str(),filename.c_str(),
         ".",&Debug,&VIsDouble) ;
  int npt = npnts ;
  TECZNE("VOLUME_MESH",&npt,&ncells,filename.c_str(),"FEBRICK") ;
  
}
void tecplot_topo_handler::close_mesh_elements() {
  int ncells = bricks.size() ;
  TECNOD(&bricks[0][0],&ncells, filename.c_str()) ;
}


void tecplot_topo_handler::close() {
}
void tecplot_topo_handler::create_mesh_positions(vector3d<float> pos[], size_t pts) {
  vector<float> pos_x(npnts),pos_y(npnts),pos_z(npnts) ;
  for(size_t i=0;i<npnts;++i) {
    pos_x[i] = pos[i].x ;
    pos_y[i] = pos[i].y ;
    pos_z[i] = pos[i].z ;
  }
  int npt = npnts ;
  TECDAT(&npt,&pos_x[0],filename.c_str()) ;
  TECDAT(&npt,&pos_y[0],filename.c_str()) ;
  TECDAT(&npt,&pos_z[0],filename.c_str()) ;
}

void tecplot_topo_handler::write_tets(Array<int,4> tets[], size_t ntets, int block, int nblocks, size_t tottets) {
  for(size_t i=0;i<ntets;++i) {
    Array<int,8> brick ;
    brick[0] = tets[i][0] ;
    brick[1] = tets[i][1] ;
    brick[2] = tets[i][2] ;
    brick[3] = brick[2] ;
    brick[4] = tets[i][3] ;
    brick[5] = brick[4] ;
    brick[6] = brick[4] ;
    brick[7] = brick[4] ;
    bricks.push_back(brick) ;
  }    
}
void tecplot_topo_handler::write_pyrm(Array<int,5> pyrm[], size_t npyrm, int block, int nblocks, size_t totpyrm) {
  for(size_t i=0;i<npyrm;++i) {
    Array<int,8> brick ;
    brick[0] = pyrm[i][0] ;
    brick[1] = pyrm[i][1] ;
    brick[2] = pyrm[i][2] ;
    brick[3] = pyrm[i][3] ;
    brick[4] = pyrm[i][4] ;
    brick[5] = brick[4] ;
    brick[6] = brick[4] ;
    brick[7] = brick[4] ;
    bricks.push_back(brick) ;
  }
  
}
void tecplot_topo_handler::write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int nblocks, size_t totprsm) {
  for(size_t i=0;i<nprsm;++i) {
    Array<int,8> brick ;
    brick[0] = prsm[i][0] ;
    brick[1] = prsm[i][1] ;
    brick[2] = prsm[i][2] ;
    brick[3] = prsm[i][2] ;
    brick[4] = prsm[i][3] ;
    brick[5] = prsm[i][4] ;
    brick[6] = prsm[i][5] ;
    brick[7] = prsm[i][5] ;
    bricks.push_back(brick) ;
  }
}
void tecplot_topo_handler::write_hexs(Array<int,8> hexs[], size_t nhexs, int block, int nblocks, size_t tothexs) {
  for(size_t i=0;i<nhexs;++i) {
    Array<int,8> brick ;
    for(int j=0;j<8;++j)
      brick[j] = hexs[i][j] ;
    bricks.push_back(brick) ;
  }
}

void tecplot_topo_handler::write_general_cell(int nfaces[], size_t nnfaces,
                                              int nsides[], size_t nnsides,
                                              int nodes[], size_t nnodes) {
  cerr << "tecplot extract module doesn't support general cells!" << endl ;
}

void tecplot_topo_handler::create_boundary_part(string name,int node_set[],
                                                size_t npnts) {
  boundary_name = name ;
  ordinary_faces.clear() ;
  node_ids.clear() ;
  for(size_t i=0;i<npnts;++i)
    node_ids.push_back(node_set[i]) ;
}
  

void tecplot_topo_handler::write_quads(Array<int,4> quads[],
                                       int quads_ids[], size_t nquads) {
  for(size_t i=0;i<nquads;++i) {
    ordinary_faces.push_back(quads[i]) ;
    elem_ids.push_back(quads_ids[i]) ;
  }
}
void tecplot_topo_handler::write_trias(Array<int,3> trias[],
                                       int trias_ids[], size_t ntrias) {
  for(size_t i=0;i<ntrias;++i) {
    Array<int,4> a ;
    a[0] = trias[i][0] ;
    a[1] = trias[i][1] ;
    a[2] = trias[i][2] ;
    a[3] = trias[i][2] ;
    
    ordinary_faces.push_back(a) ;
    elem_ids.push_back(trias_ids[i]) ;
  }
}

void tecplot_topo_handler::write_general_face(int nside_sizes[],
                                              int nside_ids[], size_t ngeneral,
                                              int nside_nodes[],
                                              size_t nside_nodes_size) {
  if(ngeneral != 0) {
    cerr << "general boundary facets ignored in tecplot extract"
         << endl ;
  }

  ofstream ofile(filename.c_str(),ios::app) ;

  int nelm = ordinary_faces.size() ;
  ofile << "ZONE T = \"" << boundary_name << '"'
        << ", N = " << npnts
        << ", E = " << nelm
        << ", ZONETYPE=FEQUADRILATERAL"
        << ",VARSHARELIST = ([1" ;
  for(int i=1;i<nvars;++i)
    ofile << ',' << i+1 ;
  ofile << "]=1)" << endl ;

  for(int i=0;i<nelm;++i)
    ofile << node_ids[ordinary_faces[i][0]-1] << ' '
          << node_ids[ordinary_faces[i][1]-1] << ' '
          << node_ids[ordinary_faces[i][2]-1] << ' '
          << node_ids[ordinary_faces[i][3]-1] << endl ;
}

void tecplot_topo_handler::close_boundary_part() {
  ordinary_faces.clear() ;
  node_ids.clear() ;
  boundary_name = "" ;
}

void tecplot_topo_handler::output_nodal_scalar(float val[], size_t npnts,
                                               string varname) {
  int npt = npnts ;
  TECDAT(&npt,&val[0],filename.c_str()) ;
}

void tecplot_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               size_t npnts, string varname) {
  vector<float> tmp(npnts) ;
  for(size_t i=0;i<npnts;++i)
    tmp[i] = val[i].x ;
  int npt = npnts ;
  TECDAT(&npt,&tmp[0],filename.c_str()) ;
  for(size_t i=0;i<npnts;++i)
    tmp[i] = val[i].y ;
  TECDAT(&npt,&tmp[0],filename.c_str()) ;
  for(size_t i=0;i<npnts;++i)
    tmp[i] = val[i].z ;
  TECDAT(&npt,&tmp[0],filename.c_str()) ;
}

void tecplot_topo_handler::output_boundary_scalar(float val[], int node_set[],
                                                  size_t nvals, string varname) {
}  

void tecplot_topo_handler::output_boundary_vector(vector3d<float> val[],
                                                  int node_set[],
                                                  size_t nvals, string varname) {
}  

bool tecplotPartConverter::processesVolumeElements() const {
  return true ;
}
bool tecplotPartConverter::processesSurfaceElements() const {
  return true ;
}
bool tecplotPartConverter::processesParticleElements() const {
  return false ;
}

void tecplotPartConverter::exportPostProcessorFiles(string casename, string iteration) const {
}

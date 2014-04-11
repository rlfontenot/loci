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
#ifndef EXTRACT_H
#define EXTRACT_H
#include <sstream>
#include <string>
#include <stdlib.h>

using std::list;

enum var_type {NODAL_SCALAR,NODAL_VECTOR,NODAL_DERIVED,NODAL_MASSFRACTION, BOUNDARY_SCALAR,BOUNDARY_VECTOR, BOUNDARY_DERIVED_SCALAR, BOUNDARY_DERIVED_VECTOR, PARTICLE_SCALAR, PARTICLE_VECTOR, UNDEFINED} ;
enum view_type {VIEWXY=0,VIEWYZ=1,VIEWXZ=2,VIEWXR=3}  ;


enum gridWritingPhases { GRID_POSITIONS, GRID_VOLUME_ELEMENTS, GRID_BOUNDARY_ELEMENTS, NODAL_VARIABLES,BOUNDARY_VARIABLES, PARTICLE_POSITIONS, PARTICLE_VARIABLES} ;

// convert a string to an integer
inline int str2int(string s) {
  std::stringstream ss ;
  ss << s ;
  int ret ;
  ss >> ret ;
  return ret ;
}

// determine whether range [b,e) contains only digits
// [b,e) must be a valid range and contains characters
template <class ForwardIter> inline bool
alldigits(ForwardIter b, ForwardIter e) {
  for(;b!=e;++b)
    if(!isdigit(*b))
      return false ;
                                                                                
  return true ;
}
// determine whether s contains a valid integer (including the sign)
inline bool
valid_int(const std::string& s) {
  if(s.empty())
    return false ;
  if( (s[0] == '-') || (s[0] == '+')) {
    if(s.size() > 1)
      return alldigits(s.begin()+1,s.end()) ;
    else
      return false ;
  }else
    return alldigits(s.begin(),s.end()) ;
}

class grid_topo_handler {
public:
  virtual ~grid_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_VOLUME_ELEMENTS ;
    sequence[2] = GRID_BOUNDARY_ELEMENTS ;
    sequence[3] = NODAL_VARIABLES ;
    sequence[4] = BOUNDARY_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) = 0 ;
  virtual void close() = 0 ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) = 0 ;
  virtual void create_mesh_elements() = 0 ;
  
  virtual void write_tets_ids(int tets_ids[], size_t ntets,int block, int numblocks,size_t tottets){} 
  virtual void write_pyrm_ids(int pyrm_ids[],size_t npyrm,int block, int numblocks,size_t totpyrm){}
  virtual void write_prsm_ids(int prsm_ids[],size_t nprsm,int block, int numblocks,size_t totprsm){}
  virtual void write_hexs_ids(int hex_ids[], size_t nhexs,int block, int numblocks,size_t tothexs){}
  virtual void write_general_cell_ids(int nfaces_ids[],size_t nnfaces){}
  
  virtual void write_tets(Array<int,4> tets[],  size_t ntets,int block, int numblocks,size_t tottets) = 0 ;
  virtual void write_pyrm(Array<int,5> pyrm[],size_t npyrm,int block, int numblocks,size_t totpyrm) = 0 ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int numblocks,size_t totprsm) = 0 ;
  virtual void write_hexs(Array<int,8> hexs[],  size_t nhexs,int block, int numblocks,size_t tothexs) = 0 ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) = 0 ;
  
  virtual void close_mesh_elements() = 0 ;

  virtual void create_boundary_part(string name,int node_set[],size_t npnts) = 0 ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) = 0 ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) = 0 ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) = 0 ;
  virtual void close_boundary_part() = 0 ;
  virtual void create_nodal_output() = 0 ;
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) = 0 ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) = 0 ;
  virtual void close_nodal_output() = 0 ;
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) = 0 ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) = 0 ;

  // maxp is the maximum particles to be extracted
  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) = 0 ;
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) = 0 ;
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) = 0 ;
} ;

class ensight_topo_handler : public grid_topo_handler {
  string dirname ;
  size_t npnts ;
  size_t ntets, nprsm, npyrm, nhexs, ngen ;
  FILE *OFP ;
  int part_id ;
  vector<vector<int> > part_nodes ;
  vector<vector<int> > part_tria_ids ;
  vector<vector<int> > part_quad_ids ;
  vector<vector<int> > part_nside_ids ;
  
  vector<vector3d<float> > positions ;

  bool particle_output ;
  string particle_geo_filename ;
  enum id_option {OFF, GIVEN, ASSIGN, IGNORE}; 
  id_option node_id_opt;
  id_option element_id_opt;
public:
  ensight_topo_handler(){  OFP=0;particle_output=false;node_id_opt = OFF; element_id_opt = OFF;}
  //constructor to set the node_id_opt and element_id_opt
  ensight_topo_handler(bool id_required){
    OFP=0;
    particle_output=false;
    if(id_required){
      node_id_opt = GIVEN;
      element_id_opt = GIVEN;
    }else{
      node_id_opt = OFF;
      element_id_opt = OFF;
    }
  }
  virtual ~ensight_topo_handler() {}
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) ;
  virtual void create_mesh_elements() {}

  virtual void write_tets(Array<int,4> tets[], size_t ntets,int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm(Array<int,5> prsm[], size_t npyrm,int block, int numblocks,size_t totprym) ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int numblocks,size_t totprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], size_t nhexs,int block, int numblocks,size_t tothexs) ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) ;

  virtual void write_tets_ids(int tets_ids[], size_t ntets,int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm_ids(int pyrm_ids[],size_t npyrm,int block, int numblocks,size_t totpyrm) ;
  virtual void write_prsm_ids(int prsm_ids[],size_t nprsm,int block, int numblocks,size_t totprsm)  ;
  virtual void write_hexs_ids(int hex_ids[], size_t nhexs,int block, int numblocks,size_t tothexs)  ;
  virtual void write_general_cell_ids(int nfaces_ids[],size_t nnfaces) ;

  
  virtual void close_mesh_elements() {}
  virtual void create_boundary_part(string name,int node_set[], size_t npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) ;
  
  virtual void create_particle_positions(vector3d<float> pos[],
					 size_t np, size_t maxp) ;
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) ;
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) ;
} ;

class tecplot_topo_handler : public grid_topo_handler {
  string filename ;
  size_t npnts ;
  size_t ntets, nprsm, npyrm, nhexs, ngen ;
  int nvars ;
  vector<Array<int, 8> > bricks ;
  string boundary_name ;
  vector<Array<int,4> > ordinary_faces ;
  vector<int> node_ids ;
  vector<int> elem_ids ;
public:
  tecplot_topo_handler(){}
  virtual ~tecplot_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = NODAL_VARIABLES ;
    sequence[2] = GRID_VOLUME_ELEMENTS ;
    sequence[3] = GRID_BOUNDARY_ELEMENTS ;
    sequence[4] = BOUNDARY_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) ;
  virtual void create_mesh_elements() {}
  virtual void write_tets(Array<int,4> tets[], size_t ntets, int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm(Array<int,5> prsm[], size_t npyrm, int block, int numblocks,size_t totpyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm, int block, int numblocks,size_t totprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], size_t nhexs, int block, int numblocks,size_t tothexs) ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) ;
  
  virtual void close_mesh_elements() ;
  virtual void create_boundary_part(string name,int node_set[], size_t npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) ;
    
  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) {}
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) {}
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) {}
} ;

class vtk_topo_handler : public grid_topo_handler {
  long long unsigned int npnts ;
  long long unsigned int ntets, nprsm, npyrm, nhexs, ngen ;
  long long unsigned int face_off ;
  long long unsigned int Offset ;
  long long unsigned int off ;
  long long unsigned int data_count ;
  string filename ;
  int nvars ;
  string boundary_name ;
  float * pos;
  float * data_store ;
  vector<int> data_size;
  vector<int> conn ;
  vector<int> cell_offsets ;
  vector<int> cell_faces ;
  vector<int> face_offsets ;
  vector<unsigned char> cell_types ;
  bool is_64bit;
  int int_size;
public:
  vtk_topo_handler(bool bit64)
  {
    Offset = 0 ; 
    data_count = 0 ;
    off = 0 ;
    face_off = 0 ;
    is_64bit = bit64;
    int_size = is_64bit ? sizeof(long long unsigned int) : sizeof(unsigned int);
  }
  virtual ~vtk_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_VOLUME_ELEMENTS ;
    sequence[2] = NODAL_VARIABLES ;
    sequence[3] = GRID_BOUNDARY_ELEMENTS ;
    sequence[4] = BOUNDARY_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) ;
  virtual void create_mesh_elements() {}
  virtual void write_tets(Array<int,4> tets[], size_t ntets, int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm(Array<int,5> prsm[], size_t npyrm, int block, int numblocks,size_t totpyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm, int block, int numblocks,size_t totprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], size_t nhexs, int block, int numblocks,size_t tothexs) ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) ;
  virtual void close_mesh_elements() ;
  virtual void create_boundary_part(string name,int node_set[], size_t npnts) {}
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) {}
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) {}
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t  nside_nodes_size) {}
  virtual void close_boundary_part() {}
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) {}
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) {}
    
  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) {}
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) {}
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) {}
} ;

class vtk_surf_topo_handler : public grid_topo_handler {
  string filename ;
  long long unsigned int npnts, ncells;
  long long unsigned int ntets, nprsm, npyrm, nhexs, ngen ;
  long long unsigned int elem_offset;
  int nvars ;
  map<int,int> bmap;
  int part_index;
  bool output_boundary;
  FILE *fid;
  vector<Array<int, 8> > bricks ;
  string boundary_name ;
  vector<int> node_ids ;
  vector<int> elem_ids ;
  vector<int> elem_conn ;
  vector<int> elem_offsets ;
  map<int,int> nmap ;
  vector<int> G2L;
  vector<unsigned char> elem_types ;
  vector<string> boundaries ;
  vector<int> data_size ;
  vector<float> elem_data ;
  vector<string> data_names ;
  float * position;
  bool is_64bit;
  int int_size;
public:
  vtk_surf_topo_handler(vector<string> &these_boundaries, bool bit64)
  {
    boundaries = these_boundaries;
    is_64bit = bit64;
    part_index = 0;
    npnts = 0; ncells = 0; elem_offset = 0; output_boundary = false;	  
    int_size = is_64bit ? sizeof(long long unsigned int) : sizeof(int);
  }
  virtual ~vtk_surf_topo_handler() {
    
  }
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_BOUNDARY_ELEMENTS ;
    sequence[1] = BOUNDARY_VARIABLES ;
    sequence[2] = GRID_POSITIONS ;
    sequence[3] = GRID_VOLUME_ELEMENTS ;
    sequence[4] = NODAL_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[],size_t npnts);
  virtual void create_mesh_elements() {}
  virtual void write_tets(Array<int,4> tets[], size_t ntets,int block, int numblocks,size_t tottets) {}
  virtual void write_pyrm(Array<int,5> prsm[], size_t npyrm, int block, int numblocks,size_t totpyrm) {}
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int numblocks,size_t totprsm) {}
  virtual void write_hexs(Array<int,8> hexs[], size_t nhexs,int block, int numblocks,size_t tothexs) {}
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) {}
  virtual void close_mesh_elements() {}
  virtual void create_boundary_part(string name,int node_set[], size_t npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], size_t npnts,string valname) {}
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) {}
  virtual void close_nodal_output() {}
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) ;
    
  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) {}
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string
				      valname) {}
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) {}
} ;

class fv_topo_handler : public grid_topo_handler {
  string dirname ;
  string filename ;
  string particleFilename ;
  size_t npnts ;
  size_t ntets, nprsm, npyrm, nhexs, ngen ;
  int part_id ;
  vector<Array<int,4> > ordinary_faces ;
  vector<vector<int> > part_nodes ;
  vector<int> elem_ids ;
  bool general_boundary ;
  bool first_var ;
  bool first_boundary ;
  FILE *OFP ;

  // particle variable buffer
  vector<string> particle_scalars_names ;
  vector<vector<float> > particle_scalars ;

  vector<string> particle_vectors_names ;
  vector<vector<vector3d<float> > > particle_vectors ;
public:
  fv_topo_handler(){OFP=0;first_var = true ; first_boundary=true ; general_boundary = false ;}
  virtual ~fv_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_BOUNDARY_ELEMENTS ;
    sequence[2] = GRID_VOLUME_ELEMENTS ;
    sequence[3] = NODAL_VARIABLES ;
    sequence[4] = BOUNDARY_VARIABLES ;
    // NOTE: in fieldview, we have to read/write the particle
    // variables first because we need to combine the scalar and
    // vector variables with the particle position into a single
    // file. we first read all the variables into a memory buffer
    // and then when we write the particle positions, we combine
    // them into a single file. therefore the order here is important.
    // we cannot start the particle position phase before the
    // particle variable phase.
    sequence[5] = PARTICLE_VARIABLES ;
    sequence[6] = PARTICLE_POSITIONS ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) ;
  virtual void create_mesh_elements() ;
  
  virtual void write_tets(Array<int,4> tets[], size_t ntets,int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm(Array<int,5> prsm[], size_t npyrm,int block, int numblocks,size_t totprym) ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int numblocks,size_t totprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], size_t nhexs,int block, int numblocks,size_t tothexs) ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) ;


  virtual void close_mesh_elements() ;
  virtual void create_boundary_part(string name,int node_set[], size_t npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) ;
    
  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) ;
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) ;
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) ;
} ;

struct affineMapping {
  float M[4][4] ;
  affineMapping() ;
  void Combine(affineMapping a) ;
  void translate(vector3d<float> tv) ;
  void rotateX(float theta) ;
  void rotateY(float theta) ;
  void rotateZ(float theta) ;
  vector3d<float> MapNode(vector3d<float> v) ;
};

class cuttingplane_topo_handler : public grid_topo_handler {
  template <int T>
    static bool arrLessThan(Array<int,T> arr1, Array<int,T> arr2) ;   // Used for sorting and searching through the Array datatype
  bool registerFace(int faceNode[], int nNodes, int cellNum) ;        // Finds cuts in the face and registers them into the vector of intersections
  void disambiguateFace(int faceNode[], int nNodes, int cellNum) ;     // Cuts a nonplanar face into planar triangles and sends each to registerFace()
  void checkLoop(int start, int end) ;     // Walks around the loop made by a cut cell and handles any double loops or errors
  int cellCount ;                          // Used to offset cell numbers to ensure unique numbering
  int numDisFaces ;                        // Total number of faces disambiguated
  int searchNodeMap(Array<int,2> arr1) ;   // Finds node arr1 in nodeMap and returns index
  vector<vector3d<float> > nodes ;         // Copy of all nodes and 3d coordinates
  vector<vector<float> > nodeVal ;         // Stores scalar property values for each node
  vector<string> variables;                // Stores property names in strings
  string strIter;                          // Stores iteration number
  vector<Array<int,5> > intersects ;       // List of all edges formed from plane intersection
  map<int, int> cellMap ;                  // Maps old 3d cell numbers to the new smaller set of 2d cells
  vector<Array<int,2> > nodeMap ;          // Maps old 3d node numbers to the new smaller set of 2d nodes
  list<vector<int> > disambiguatedFaces ;  // Stores all unmatched previously disambiguated faces and corresponding fabricated nodes
  list<vector<int> > resolvedFace;         // Stores all matched disambiguated faces
  affineMapping transMatrix ;              // Transformation matrix that transforms the grid topology to the desired position
  int tetsCut, prsmCut, pyrmCut, hexsCut ;
public:
  cuttingplane_topo_handler(affineMapping &transformMatrix,
			    float xShift, float yShift, float zShift);
  virtual ~cuttingplane_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = NODAL_VARIABLES ;
    sequence[1] = GRID_POSITIONS ;
    sequence[2] = GRID_VOLUME_ELEMENTS ;
    sequence[3] = GRID_BOUNDARY_ELEMENTS ;
    sequence[4] = BOUNDARY_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,size_t npnts,
                    size_t ntets, size_t nprsm, size_t npyrm, size_t nhexs, size_t ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], size_t npnts) ;
  virtual void create_mesh_elements() ;
  
  virtual void write_tets(Array<int,4> tets[], size_t ntets,int block, int numblocks,size_t tottets) ;
  virtual void write_pyrm(Array<int,5> pyrm[], size_t npyrm,int block, int numblocks,size_t totpyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int numblocks,size_t totprsm) ;
  virtual void write_hexs(Array<int,8> hexs[], size_t  nhexs,int block, int numblocks,size_t tothexs) ;
  virtual void write_general_cell(int nfaces[], size_t nnfaces,
                                  int nsides[], size_t nnsides,
                                  int nodes[], size_t nnodes) ;

  virtual void close_mesh_elements() ;

  virtual void create_boundary_part(string name,int node_set[],size_t npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           size_t nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           size_t ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], size_t ngeneral,
                             int nside_nodes[], size_t nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() ;
  virtual void output_nodal_scalar(float val[], size_t npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   size_t npnts, string valname) ;
  virtual void close_nodal_output() ;
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      size_t nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      size_t nvals, string valname) ;

  virtual void create_particle_positions(vector3d<float> pos[],
                                         size_t np, size_t maxp) {}
  virtual void output_particle_scalar(float val[], size_t np,
                                      size_t maxp, string valname) {}
  virtual void output_particle_vector(vector3d<float> val[], size_t np,
                                      size_t maxp, string valname) {}
} ;

void get_2dgv(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries,
              int view) ;

void get_surf(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries) ;

void process_ascii_nodal(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames) ;

void process_ascii_bndry(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames,
                         vector<string> boundaries) ;

void process_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter) ;

void combine_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter,
		  bool do_favre) ;

size_t  sizeElementType(hid_t group_id, const char *element_name) ;

template<class T> void readElementType(hid_t group_id, const char *element_name,
                                       vector<T> &v) {
  if(v.size() > 0) {
    hid_t dataset = H5Dopen(group_id,element_name) ;

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }
}

template<class T>
void readElementTypeBlock(hid_t group_id,
                          const char *element_name,
                          vector<T> &v,
                          size_t start_elem,
                          int block_size) {
  if(block_size > 0) {
    hid_t dataset = H5Dopen(group_id,element_name) ;

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;

    hsize_t count = block_size ;
    hsize_t start = start_elem ;
    hid_t dataspace = H5Dget_space(dataset) ;
    hsize_t stride = 1 ;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
    int rank = 1 ;
    hsize_t dimension = count ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    hid_t datatype = dp->get_hdf5_type() ;
    H5Dread(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&v[0]) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;
  }
}


template<class T> void writeElementType(hid_t group_id,
                                        const char *element_name,
                                        std::vector<T> &v) {
  hsize_t array_size = v.size() ;
  if(array_size == 0)
    return ;
  int rank = 1 ;
  hsize_t dimension = array_size ;

  hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

  typedef data_schema_traits<T> traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  
#ifdef H5_INTERFACE_1_6_4
  hsize_t start = 0 ;
#else
  hssize_t start = 0 ;
#endif
  hsize_t stride = 1 ;
  hsize_t count = v.size() ;
  hid_t datatype = dp->get_hdf5_type() ;
  hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                            dataspace, H5P_DEFAULT) ;
  if(count != 0) {
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                        &start, &stride, &count, NULL) ;
    hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
    H5Dwrite(dataset,datatype,memspace,dataspace,
             H5P_DEFAULT, &v[0]) ;
    H5Sclose(memspace) ;
  }
  H5Dclose(dataset) ;
  H5Sclose(dataspace) ;
  H5Tclose(datatype) ;
}
  
  
void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) ;

namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]) ;
  }
}

extern string output_dir ;

string getPosFile(string output_dir,string iteration, string casename) ;


class particlePartBase : public Loci::CPTR_type {
 protected:
  bool error ;
  string partName ;
  size_t numParticles ;
 public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  size_t getNumParticles() const { return numParticles ; }
  virtual bool hasScalarVar(string var) const = 0 ;
  virtual bool hasVectorVar(string var) const = 0 ;
  virtual std::vector<string> getScalarVars() const = 0 ;
  virtual std::vector<string> getVectorVars() const = 0 ;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const = 0;
  virtual void getParticleScalar(string varname, vector<float> &val) const = 0 ;
  virtual void getParticleVector(string varname, 
				 vector<vector3d<float> > &val) const  = 0 ;
} ;

typedef Loci::CPTR<particlePartBase> particlePartP ;

class particlePart: public particlePartBase {
  string directory ;
  string posfile ;
  map<string,string> scalarVars ;
  map<string,string> vectorVars ;
  int stride_size ;
 public:
  particlePart() { error = true ; }
  particlePart(string output_dir, string iteration, string casename,
	       vector<string> vars, int maxparticles) ;
  virtual bool hasScalarVar(string var) const ;
  virtual bool hasVectorVar(string var) const ;
  virtual std::vector<string> getScalarVars() const ;
  virtual std::vector<string> getVectorVars() const ;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const ;
  virtual void getParticleScalar(string varname, vector<float> &val) const ;
  virtual void getParticleVector(string varname, 
				 vector<vector3d<float> > &val) const ;
} ;

// Create abstraction for parts
class volumePartBase : public Loci::CPTR_type {
 protected:
  bool error ;
  string partName ;
  size_t nnodes ;
  size_t ntets, nhexs, nprsm, npyrm, ngenc ; 
  size_t ntetsIblank, nhexsIblank,nprsmIblank,npyrmIblank,ngencIblank ;
 public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  size_t getNumNodes() const { return nnodes ; }
  
  size_t getNumTets() const { return ntets ; }
  size_t getNumHexs() const { return nhexs ; }
  size_t getNumPrsm() const { return nprsm ; }
  size_t getNumPyrm() const { return npyrm ; }
  size_t getNumGenc() const { return ngenc ; }
  size_t getNumTetsIblank() const { return ntetsIblank ; }
  size_t getNumHexsIblank() const { return nhexsIblank ; }
  size_t getNumPrsmIblank() const { return nprsmIblank ; }
  size_t getNumPyrmIblank() const { return npyrmIblank ; }
  size_t getNumGencIblank() const { return ngencIblank ; }

  virtual bool hasNodalScalarVar(string var) const = 0 ;
  virtual bool hasNodalVectorVar(string var) const = 0 ;
  virtual std::vector<string> getNodalScalarVars() const = 0 ;
  virtual std::vector<string> getNodalVectorVars() const = 0 ;
  virtual void getPos(vector<vector3d<float> > &pos) const = 0 ;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const = 0 ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const = 0 ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const = 0 ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const = 0 ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const = 0 ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const = 0 ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const = 0 ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const = 0 ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const = 0;
  virtual void getGenIds(vector<int> &genids) const = 0 ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const = 0 ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0 ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const = 0 ;
} ;

typedef Loci::CPTR<volumePartBase> volumePartP ;

class volumePartDerivedVars : public volumePartBase {
  volumePartP shadowPart ;
  enum derivedVar_t { VAR_M, VAR_P, VAR_logp, VAR_U, VAR_0, VAR_1, VAR_2,
		      VAR_X, VAR_Y, VAR_Z } ;
  map<string,derivedVar_t> derivedVars ;
  float Pambient ;
  void processDerivedVars(const vector<string> &vars) ;
 public:
  volumePartDerivedVars() {error = true ;}
  volumePartDerivedVars(volumePartP part,
			string output_dir, string iteration, string casename,
			vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual std::vector<string> getNodalScalarVars() const ;
  virtual std::vector<string> getNodalVectorVars() const ;
  virtual void getPos(vector<vector3d<float> > &val) const ;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const ;
  virtual void getGenIds(vector<int> &genids) const ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const ;
} ;

class volumePart : public volumePartBase {
  string directory ;
  string topoFile ;
  string posFile ;
  // maps from variables to file name
  map<string,string> nodalScalarVars ;
  map<string,string> nodalVectorVars ;
  bool has_iblank ;
  store<unsigned char> iblank ;
  size_t ntets_orig, nhexs_orig, nprsm_orig, npyrm_orig, ngenc_orig ; 
  Loci::entitySet tetsIblanked, hexsIblanked, prsmIblanked, pyrmIblanked ;
  Loci::entitySet gencIblanked ;
 public:
  volumePart() {error = true ;}
  volumePart(string output_dir, string iteration, string casename,
	      vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual std::vector<string> getNodalScalarVars() const ;
  virtual std::vector<string> getNodalVectorVars() const ;
  virtual void getPos(vector<vector3d<float> > &val) const ;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const ;
  virtual void getGenIds(vector<int> &genids) const ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const ;
} ;

class surfacePartBase : public Loci::CPTR_type {
 protected:
  bool error ;
  string partName ;
  int nnodes, nquads, ntrias, ngenf ;
 public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  int getNumNodes() const { return nnodes ; }
  int getNumQuads() const { return nquads ; }
  int getNumTrias() const { return ntrias ; }
  int getNumGenfc() const { return ngenf ; }
  virtual bool hasNodalScalarVar(string var) const = 0 ;
  virtual bool hasNodalVectorVar(string var) const = 0 ;
  virtual bool hasElementScalarVar(string var) const = 0 ;
  virtual bool hasElementVectorVar(string var) const = 0 ;
  virtual std::vector<string> getNodalScalarVars() const = 0 ;
  virtual std::vector<string> getNodalVectorVars() const = 0 ;
  virtual std::vector<string> getElementScalarVars() const = 0 ;
  virtual std::vector<string> getElementVectorVars() const = 0 ;
  virtual void getQuads(vector<Array<int,4> > &quads) const = 0 ;
  virtual void getTrias(vector<Array<int,3> > &trias) const = 0 ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const =0;
  virtual void getPos(vector<vector3d<float> > &pos) const = 0;
  virtual void getNodalScalar(string varname, vector<float> &vals) const = 0 ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const =0 ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const = 0 ;
} ;
typedef Loci::CPTR<surfacePartBase> surfacePartP ;

class surfacePartDerivedVars : public surfacePartBase {
  surfacePartP shadowPart ;
  enum derivedVar_t { VAR_M, VAR_P, VAR_logp, VAR_U, VAR_0, VAR_1, VAR_2,
		      VAR_X, VAR_Y, VAR_Z } ;
  map<string,derivedVar_t> derivedVars ;
  float Pambient ;
  void processDerivedVars(const vector<string> &vars) ;
 public:
  surfacePartDerivedVars() { error = true ; shadowPart=0;} 
  surfacePartDerivedVars(surfacePartP part, string output_dir,
			 string casename ,
			 string iteration, 
			 vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;
    
      
class surfacePart : public surfacePartBase {
  string directory ;
  string topoFile ;
  string posFile ;
  // maps from variables to file name
  map<string,string> nodalScalarVars ;
  map<string,string> nodalVectorVars ;
  map<string,string> elementScalarVars ;
  map<string,string> elementVectorVars ;

  vector<int> quad_ord ;
  vector<int> tri_ord ;
  vector<int> gen_ord ;
  entitySet quadSet,triSet,genSet ;
 public:
  surfacePart() {error = true ;}
  surfacePart(string name, string directory, string iteration,
	      vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;

class surfacePartCopy : public surfacePartBase {
  vector<Array<int,3> > trifaces ;
  vector<Array<int,4> > quadfaces ;
  vector<int> nfacenodes, gennodes ;
  vector<int> nodemap ;
  vector<vector3d<float> > pos ;
  map<string,vector<float> > nodalScalars ;
  map<string,vector<vector3d<float> > > nodalVectors ;
  map<string,Array<vector<float>,3> > elementScalars ;
  map<string,Array<vector<vector3d<float> >,3> > elementVectors ;
 public:
  surfacePartCopy() {error = true ;}
  surfacePartCopy(string name, vector<Array<int,3> > &triangles,
                  vector<Array<int,4> > &quads,
                  vector<int> &genface2n, vector<int> &gnodes) ;
  void registerPos(const vector<vector3d<float> > &pos) ;
  void registerNodalScalar(string name,const vector<float> &val) ;
  void registerNodalVector(string name,const vector<vector3d<float> > &val) ;
  void registerElementScalar(string name, 
			     const vector<float> &qval,
			     const vector<float> &tval,
			     const vector<float> &gval) ;
  void registerElementVector(string name, 
			     const vector<vector3d<float> > &qval,
			     const vector<vector3d<float> > &tval,
			     const vector<vector3d<float> > &gval) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;

class postProcessorConvert : public Loci::CPTR_type {
 protected:
  vector<surfacePartP> surfacePartList ;
  vector<volumePartP> volumePartList ;
  vector<particlePartP> particlePartList ;
 public:
  void addSurfaceParts(const vector<surfacePartP> &list) {
    for(size_t i=0;i<list.size();++i)
      surfacePartList.push_back(list[i]) ;
  }
  void addVolumePart(volumePartP volpart) {
    volumePartList.push_back(volpart) ;
  }
  void addParticlePart(particlePartP particles) {
    particlePartList.push_back(particles) ;
  }
  virtual bool processesVolumeElements() const = 0 ;
  virtual bool processesSurfaceElements() const = 0 ;
  virtual bool processesParticleElements() const = 0 ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const = 0 ;
} ;

typedef Loci::CPTR<postProcessorConvert> postProcessorP ;

class ensightPartConverter : public postProcessorConvert {
 public:
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class tecplotPartConverter : public postProcessorConvert {
 public:
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class vtkPartConverter : public postProcessorConvert {
  bool bit64 ;
 public:
  vtkPartConverter() {bit64 = false; } 
  vtkPartConverter(bool input) {bit64 = input; } ;
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class vtkSurfacePartConverter : public postProcessorConvert {
  bool bit64 ;
 public:
  vtkSurfacePartConverter() {bit64=false; }
  vtkSurfacePartConverter(bool input) { bit64=input; }
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class fieldViewPartConverter : public postProcessorConvert {
 public:
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class cuttingPlanePartConverter : public postProcessorConvert 
{
  affineMapping transformMatrix ;
  float xShift, yShift, zShift ;
 public:
  cuttingPlanePartConverter(const affineMapping &m,
			    float xs,float ys,float zs) {
    transformMatrix = m ;
    xShift=xs ;
    yShift=ys ;
    zShift=zs ;
  }
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

#endif

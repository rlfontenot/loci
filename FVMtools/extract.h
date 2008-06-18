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

using std::list;

enum var_type {NODAL_SCALAR,NODAL_VECTOR,NODAL_DERIVED,NODAL_MASSFRACTION, BOUNDARY_SCALAR,BOUNDARY_VECTOR, BOUNDARY_DERIVED_SCALAR, BOUNDARY_DERIVED_VECTOR, PARTICLE_SCALAR, PARTICLE_VECTOR, UNDEFINED} ;
enum view_type {VIEWXY=0,VIEWYZ=1,VIEWXZ=2,VIEWXR=3}  ;


enum gridWritingPhases { GRID_POSITIONS, GRID_VOLUME_ELEMENTS, GRID_BOUNDARY_ELEMENTS, NODAL_VARIABLES,BOUNDARY_VARIABLES, PARTICLE_POSITIONS, PARTICLE_VARIABLES} ;

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
  virtual void open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types) = 0 ;
  virtual void close() = 0 ;
  virtual void create_mesh_positions(vector3d<float> pos[], int npnts) = 0 ;
  virtual void create_mesh_elements() = 0 ;
  virtual void write_tets(Array<int,4> tets[], int ntets) = 0 ;
  virtual void write_pyrm(Array<int,5> pyrm[], int npyrm) = 0 ;
  virtual void write_prsm(Array<int,6> prsm[], int nprsm) = 0 ;
  virtual void write_hexs(Array<int,8> hexs[], int nhexs) = 0 ;
  virtual void write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[], int nnsides,
                                  int nodes[], int nnodes) = 0 ;
  virtual void close_mesh_elements() = 0 ;

  virtual void create_boundary_part(string name,int node_set[],int npnts) = 0 ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           int nquads) = 0 ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           int ntrias) = 0 ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
                             int nside_nodes[], int nside_nodes_size) = 0 ;
  virtual void close_boundary_part() = 0 ;
  virtual void create_nodal_output() = 0 ;
  virtual void output_nodal_scalar(float val[], int npnts, string valname) = 0 ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   int npnts, string valname) = 0 ;
  virtual void close_nodal_output() = 0 ;
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      int nvals, string valname) = 0 ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      int nvals, string valname) = 0 ;

  virtual void create_particle_positions(vector3d<float> pos[], int np) = 0 ;
  virtual void output_particle_scalar(float val[],
                                      int np, string valname) = 0 ;
  virtual void output_particle_vector(vector3d<float> val[],
                                      int np, string valname) = 0 ;
} ;

class ensight_topo_handler : public grid_topo_handler {
  string dirname ;
  int npnts ;
  int ntets, nprsm, npyrm, nhexs, ngen ;
  FILE *OFP ;
  int part_id ;
  vector<vector<int> > part_nodes ;
  vector<vector<int> > part_tria_ids ;
  vector<vector<int> > part_quad_ids ;
  vector<vector<int> > part_nside_ids ;
  
  vector<vector3d<float> > positions ;

  bool particle_output ;
  string particle_geo_filename ;
public:
  ensight_topo_handler(){OFP=0;particle_output=false;}
  virtual ~ensight_topo_handler() {}
  virtual void open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], int npnts) ;
  virtual void create_mesh_elements() {}
  virtual void write_tets(Array<int,4> tets[], int ntets) ;
  virtual void write_pyrm(Array<int,5> prsm[], int npyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], int nprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], int nhexs) ;
  virtual void write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[], int nnsides,
                                  int nodes[], int nnodes) ;
  virtual void close_mesh_elements() {}
  virtual void create_boundary_part(string name,int node_set[], int npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           int nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           int ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
                             int nside_nodes[], int nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], int npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   int npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      int nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      int nvals, string valname) ;
  
  virtual void create_particle_positions(vector3d<float> pos[], int np) ;
  virtual void output_particle_scalar(float val[], int np, string valname) ;
  virtual void output_particle_vector(vector3d<float> val[],
                                      int np, string valname) ;
} ;

class tecplot_topo_handler : public grid_topo_handler {
  string filename ;
  int npnts ;
  int ntets, nprsm, npyrm, nhexs, ngen ;
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
  virtual void open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], int npnts) ;
  virtual void create_mesh_elements() {}
  virtual void write_tets(Array<int,4> tets[], int ntets) ;
  virtual void write_pyrm(Array<int,5> prsm[], int npyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], int nprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], int nhexs) ;
  virtual void write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[], int nnsides,
                                  int nodes[], int nnodes) ;
  virtual void close_mesh_elements() ;
  virtual void create_boundary_part(string name,int node_set[], int npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           int nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           int ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
                             int nside_nodes[], int nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], int npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   int npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      int nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      int nvals, string valname) ;
    
  virtual void create_particle_positions(vector3d<float> pos[], int np) {}
  virtual void output_particle_scalar(float val[],
                                      int np, string valname) {}
  virtual void output_particle_vector(vector3d<float> val[],
                                      int np, string valname) {}
} ;

class fv_topo_handler : public grid_topo_handler {
  string filename ;
  int npnts ;
  int ntets, nprsm, npyrm, nhexs, ngen ;
  int part_id ;
  vector<Array<int,4> > ordinary_faces ;
  vector<vector<int> > part_nodes ;
  vector<int> elem_ids ;
  bool general_boundary ;
  bool first_var ;
  bool first_boundary ;
  FILE *OFP ;
public:
  fv_topo_handler(){OFP=0;first_var = true ; first_boundary=true ; general_boundary = false ;}
  virtual ~fv_topo_handler() {}
  virtual void fileWritingSequence(Array<int,7> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_BOUNDARY_ELEMENTS ;
    sequence[2] = GRID_VOLUME_ELEMENTS ;
    sequence[3] = NODAL_VARIABLES ;
    sequence[4] = BOUNDARY_VARIABLES ;
    sequence[5] = PARTICLE_POSITIONS ;
    sequence[6] = PARTICLE_VARIABLES ;
  }
  virtual void open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], int npnts) ;
  virtual void create_mesh_elements() ;
  virtual void write_tets(Array<int,4> tets[], int ntets) ;
  virtual void write_pyrm(Array<int,5> prsm[], int npyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], int nprsm)  ;
  virtual void write_hexs(Array<int,8> hexs[], int nhexs) ;
  virtual void write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[], int nnsides,
                                  int nodes[], int nnodes) ;
  virtual void close_mesh_elements() ;
  virtual void create_boundary_part(string name,int node_set[], int npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           int nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           int ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
                             int nside_nodes[], int nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() {}
  virtual void output_nodal_scalar(float val[], int npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   int npnts, string valname) ;
  virtual void close_nodal_output() {} ; 
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      int nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      int nvals, string valname) ;
    
  virtual void create_particle_positions(vector3d<float> pos[], int np) {}
  virtual void output_particle_scalar(float val[],
                                      int np, string valname) {}
  virtual void output_particle_vector(vector3d<float> val[],
                                      int np, string valname) {}
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
  virtual void open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types) ;
  virtual void close() ;
  virtual void create_mesh_positions(vector3d<float> pos[], int npnts) ;
  virtual void create_mesh_elements() ;
  virtual void write_tets(Array<int,4> tets[], int ntets) ;
  virtual void write_pyrm(Array<int,5> pyrm[], int npyrm) ;
  virtual void write_prsm(Array<int,6> prsm[], int nprsm) ;
  virtual void write_hexs(Array<int,8> hexs[], int nhexs) ;
  virtual void write_general_cell(int nfaces[], int nnfaces,
                                  int nsides[], int nnsides,
                                  int nodes[], int nnodes) ;
  virtual void close_mesh_elements() ;

  virtual void create_boundary_part(string name,int node_set[],int npnts) ;
  virtual void write_quads(Array<int,4> quads[], int quad_ids[],
                           int nquads) ;
  virtual void write_trias(Array<int,3> trias[], int tria_ids[],
                           int ntrias) ;
  virtual void write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,
                             int nside_nodes[], int nside_nodes_size) ;
  virtual void close_boundary_part() ;
  virtual void create_nodal_output() ;
  virtual void output_nodal_scalar(float val[], int npnts, string valname) ;
  virtual void output_nodal_vector(vector3d<float> val[],
                                   int npnts, string valname) ;
  virtual void close_nodal_output() ;
  virtual void output_boundary_scalar(float val[], int node_set[],
                                      int nvals, string valname) ;
  virtual void output_boundary_vector(vector3d<float> val[], int node_set[],
                                      int nvals, string valname) ;

  virtual void create_particle_positions(vector3d<float> pos[], int np) ;
  virtual void output_particle_scalar(float val[],
                                      int np, string valname) ;
  virtual void output_particle_vector(vector3d<float> val[],
                                      int np, string valname) ;
} ;

void get_2dgv(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries,
              int view) ;
void process_ascii_nodal(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames) ;

void process_ascii_bndry(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames,
                         vector<string> boundaries) ;
int  sizeElementType(hid_t group_id, const char *element_name) ;

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

void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) ;

namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]) ;
  }
}


#endif

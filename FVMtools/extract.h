#ifndef EXTRACT_H
#define EXTRACT_H

enum var_type {NODAL_SCALAR,NODAL_VECTOR,NODAL_DERIVED,NODAL_MASSFRACTION, BOUNDARY_SCALAR,BOUNDARY_VECTOR, BOUNDARY_DERIVED_SCALAR, BOUNDARY_DERIVED_VECTOR, UNDEFINED} ;
enum view_type {VIEWXY=0,VIEWYZ=1,VIEWXZ=2,VIEWXR=3}  ;


enum gridWritingPhases { GRID_POSITIONS, GRID_VOLUME_ELEMENTS, GRID_BOUNDARY_ELEMENTS, NODAL_VARIABLES,BOUNDARY_VARIABLES } ;

class grid_topo_handler {
public:
  virtual ~grid_topo_handler() {}
  virtual void fileWritingSequence(Array<int,5> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_VOLUME_ELEMENTS ;
    sequence[2] = GRID_BOUNDARY_ELEMENTS ;
    sequence[3] = NODAL_VARIABLES ;
    sequence[4] = BOUNDARY_VARIABLES ;
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
public:
  ensight_topo_handler(){OFP=0;}
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
    
} ;

class tecplot_topo_handler : public grid_topo_handler {
  string filename ;
  int npnts ;
  int ntets, nprsm, npyrm, nhexs, ngen ;
  int part_id ;
  vector<Array<int, 8> > bricks ;
public:
  tecplot_topo_handler(){}
  virtual ~tecplot_topo_handler() {}
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
    
} ;

class fv_topo_handler : public grid_topo_handler {
  string filename ;
  int npnts ;
  int ntets, nprsm, npyrm, nhexs, ngen ;
  int part_id ;
  vector<Array<int,4> > ordinary_faces ;
  vector<vector<int> > part_nodes ;
  bool first_var ;
  bool first_boundary ;
  FILE *OFP ;
public:
  fv_topo_handler(){OFP=0;first_var = true ; first_boundary=true ;}
  virtual ~fv_topo_handler() {}
  virtual void fileWritingSequence(Array<int,5> &sequence) {
    sequence[0] = GRID_POSITIONS ;
    sequence[1] = GRID_BOUNDARY_ELEMENTS ;
    sequence[2] = GRID_VOLUME_ELEMENTS ;
    sequence[3] = NODAL_VARIABLES ;
    sequence[4] = BOUNDARY_VARIABLES ;
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

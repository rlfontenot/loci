#include <Loci.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>

#include <Tools/xdr.h>

#include <sys/types.h>
#include <sys/stat.h>
using std::istream ;
using std::ostream ;
using std::endl ;
using std::cin ;
using std::cout ;
using std::cerr ;
using std::ios ;
using std::ofstream ;
using std::ifstream ;

using std::istringstream ;
using std::ostringstream ;
using std::string ;
using std::char_traits ;
#include <vector>
using std::vector ;
#include <map>
using std::map ;

void ENUsage(char *s) {
  cerr << "Usage: " << endl 
       << string(s) << " <problem_name> <time_step> <key> " << endl << endl ;
  cerr << "where <key> may be: " << endl
       << "r - extract density " << endl
       << "p - extract log10 pressure" << endl
       << "P - extract actual pressure" << endl 
       << "v - extract velocity vector uvel, vvel, wvel" << endl
       << "U - extract the magnitude of velocity " << endl  
       << "m - extract mach number" << endl
       << "t - extract temperature" << endl
       << "a - extract mach number" << endl
       << "G - extract Gamma" << endl
       << "R - extract mixture gas constant R~" << endl
       << "q - extract heat flux at the boundary" << endl
       << "y - extract y_plus at the boundary" << endl
       << "w - extract wall stress ws_x, ws_y, ws_z at the boundary" << endl
       << "wp - extract wall pressure at the boundary" << endl
       << "wt - extract wall temperature at the boundary" << endl
       << "f<species> - extract mass fractions for <species>" << endl << endl ;
  
  cerr << "example:  to extract OH species from time step 800 of ssme simulation use:" << endl
       << string(s) << "-fv ssme 200 fOH" << endl ;
}

namespace {
bool compare_face(const std::vector<int> &face1, const std::vector<int> &face2) {
  if(face1.size() != face2.size())
 return false ;
  vector<int>::size_type fsz = face1.size() ;
  vector<int>::size_type ind ;
  for(ind=0;ind<fsz;++ind)
    if(face2[ind] == face1[0])
      break ;
  if(ind == face2.size())
    return false ;
  vector<int>::size_type l ;
  for( l=1;l<fsz;++l) {
    int li = (fsz+ind-l)%fsz ;
    if(face1[l] != face2[li])
      break ;
  }
  if(l == fsz)
    return true ;
  for( l=1;l<fsz;++l) {
    int li = (ind+l)%fsz ;
    if(face1[l] != face2[li])
      break ;
  }
  return (l == fsz) ;
}
}

// data structure to held the boundary info
struct qf_boundary {
  vector<int> facenode ;
  float q, y, wp, wt ;
  float w[3] ;
  qf_boundary(float q=0.0,float y=0.0,
              float wp=0.0,float wt=0.0,
              float w_x=0.0,float w_y=0.0,
              float w_z=0.0):q(q),y(y),wp(wp),
                             wt(wt) {
    w[0]=w_x,w[1]=w_y,w[2]=w_z ;
  }
} ;

int EnsightExtract(int ac, char *av[]) {
  if(ac <= 1) {
    ENUsage(av[0]) ;
    exit(0) ;
  }
  char bstr[] = "grid" ;
  char *ncyc = bstr ;
  char tmp_buf[80] ;
  char *problem_name = av[1] ;
  char filename[512] ;
  ac-- ;
  av++ ;
  if(ac >= 2) {
    ncyc = av[1] ; 
    av++ ;
    ac-- ; 
  }
  char case_name[256], geo_name[256], grid_name[256], name[256], dir_name[256] ;
  int extra_vars = 0 ;
  int extra_vec_vars = 0 ;
  int num_extra_vars = 0 ;
  sprintf(dir_name, "%s_ensight", problem_name) ;
  sprintf(case_name, "%s/%s.case", dir_name, problem_name) ;
  sprintf(geo_name, "%s/%s.geo", dir_name, problem_name) ;
  sprintf(name, "%s.geo", problem_name) ;
  sprintf(grid_name, "%s.xdr", problem_name) ;
  int num_vars = 0 ;
  int num_bvars = 0 ;
  char** var_buf = new char*[20] ;
  char** vector_var_files = new char*[20] ;
  char** scalar_var_files = new char*[20] ;
  // buffers for boundary variables
  char** bvar_buf = new char*[20] ;
  char** bvector_var_files = new char*[20] ;
  char** bscalar_var_files = new char*[20] ;
  
  struct stat statbuf ;
  if(stat(dir_name,&statbuf))
    mkdir(dir_name,0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file " << dir_name << " should be a directory!, rename 'ensight' and start again."
           << endl ;
      exit(-1) ;
    }
  for(int i = 0; i < 20; ++i) {
    var_buf[i] = new char[80] ;
    vector_var_files[i] = new char[80] ;
    scalar_var_files[i] = new char[80] ;
  }
  for(int i=0;i<20;++i) {
    bvar_buf[i] = new char[80] ;
    bvector_var_files[i] = new char[80] ;
    bscalar_var_files[i] = new char[80] ;
  }

  // unused variables
  //char** vector_bvar_files = new char*[8] ;
  //char** scalar_bvar_files = new char*[8] ;
  //char** bvar_buf = new char*[8] ;
  int tmp = 0 ;
  int scalar_file_tmp = 0 ;
  int vector_file_tmp = 0 ;
  int tmp_b = 0 ;
  int bscalar_file_tmp = 0 ;
  int bvector_file_tmp = 0 ;
  entitySet vec_indices, scalar_indices ;
  entitySet bvec_indices, bscalar_indices ;
  dMap vec_map ;
  // flag to indicate whether to read qf file
  bool read_boundary_value = false ; 
  if(ac == 1) {
    sprintf(var_buf[0], "m") ;
    sprintf(scalar_var_files[0], "%s/%s.m", dir_name, problem_name) ;
    scalar_indices += 0 ;
  }
  while(ac > 1) {
    if(!strcmp(av[1],"v")) {
      sprintf(var_buf[tmp], "v") ;
      sprintf(vector_var_files[vector_file_tmp], "%s/%s.v", dir_name, problem_name) ;
      vec_indices += tmp ;
      vec_map[tmp++] = vector_file_tmp++ ;
    } else if(!strcmp(av[1],"w")) {
      sprintf(bvar_buf[tmp_b],"w") ;
      sprintf(bvector_var_files[bvector_file_tmp],"%s/%s.w",dir_name,problem_name) ;
      bvec_indices += tmp_b ;
      tmp_b++ ;
      bvector_file_tmp++ ;
      
      if(!read_boundary_value)
        read_boundary_value = true ;
    } else if (!strcmp(av[1],"q") ||
             !strcmp(av[1],"y") ||
             !strcmp(av[1],"wp")||
             !strcmp(av[1],"wt")) {
      sprintf(bvar_buf[tmp_b],av[1]) ;
      sprintf(bscalar_var_files[bscalar_file_tmp], "%s/%s.%s", dir_name, problem_name, av[1]) ;
      bscalar_indices += tmp_b ;
      tmp_b++ ;
      bscalar_file_tmp++ ;
      
      if(!read_boundary_value)
        read_boundary_value = true ;
    } else {
      var_buf[tmp] = av[1] ;
      if(var_buf[tmp][0] == 'n') {
        std::string s = std::string(av[1] + 1) ;
	sprintf(scalar_var_files[scalar_file_tmp++], "%s/%s.%s", dir_name, problem_name, s.c_str()) ;
	scalar_indices += tmp++ ;
	num_extra_vars++ ;
      } else if(var_buf[tmp][0] == 'f') {
	std::string s = string(var_buf[tmp] + 1) ;
	sprintf(scalar_var_files[scalar_file_tmp++], "%s/%s.%s", dir_name, problem_name, s.c_str()) ;
	scalar_indices += tmp++ ;
      } else if(var_buf[tmp][0] == 'x') {
        std::string s = string(var_buf[tmp] + 1) ;
	sprintf(vector_var_files[vector_file_tmp], "%s/%s.%s", dir_name, problem_name, s.c_str()) ;
        vec_indices += tmp ;
        vec_map[tmp++] = vector_file_tmp++ ;
        //	vec_indices += tmp++ ;
	num_extra_vars++ ;
      } 
      else {
	sprintf(scalar_var_files[scalar_file_tmp++], "%s/%s.%s", dir_name, problem_name, av[1]) ;
	scalar_indices += tmp++ ;
      }
    }
    av++ ; 
    ac-- ;
  }
  
  if(tmp > 0)
    num_vars = tmp ;
  if(tmp_b > 0)
    num_bvars = tmp_b ;
  int *var_map = new int[num_extra_vars] ;
  cout << "Total number of variables = " << num_vars+num_bvars << endl ;
  cout << " Number of variables to be read from other files = " << num_extra_vars << endl ;
  for(int i = 0; i < scalar_file_tmp ; ++i)
    cout << "scalar variable files to be written " << scalar_var_files[i] << endl ;
  for(int i=0;i<bscalar_file_tmp;++i)
    cout << "boundary scalar variable files to be written "
         << bscalar_var_files[i] << endl ;
  for(int i = 0; i < vector_file_tmp ; ++i)
    cout << "vector variable files to be written " << vector_var_files[i] << endl ;
  for(int i=0;i<bvector_file_tmp;++i)
    cout << "boundary vector variable files to be written "
         << bvector_var_files[i] << endl ;
  ofstream of(case_name, ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << name << endl ;
  if(scalar_file_tmp || vector_file_tmp ||
     bscalar_file_tmp || bvector_file_tmp)
    of << "VARIABLE" << endl ;
  entitySet::const_iterator ei = scalar_indices.begin() ;
  for(int i = 0; i < scalar_file_tmp; ++i) {
    memset(tmp_buf, '\0', 80) ;
    if(var_buf[*ei][0] == 'n' || var_buf[*ei][0] == 'f') {
      std::string s = std::string(var_buf[*ei] + 1) ;
      sprintf(tmp_buf, "%s.%s", problem_name, s.c_str()) ;
      of << "scalar per node:\t" << s.c_str() << '\t' << tmp_buf << endl ;
    }
    else {
      sprintf(tmp_buf, "%s.%s", problem_name, var_buf[*ei]) ;
      of << "scalar per node:\t" << var_buf[*ei]  << '\t' << tmp_buf << endl ;
    }
    ++ei ;
  }
  ei = bscalar_indices.begin() ;
  for(int i=0;i<bscalar_file_tmp;++i) {
    memset(tmp_buf,'\0',80) ;
    sprintf(tmp_buf,"%s.%s",problem_name,bvar_buf[*ei]) ;
    of << "scalar per element:\t" << bvar_buf[*ei]
       << '\t' << tmp_buf << endl ;
    ++ei ;
  }
  ei = vec_indices.begin() ;
  for(int i = 0; i < vector_file_tmp; ++i) {
    memset(tmp_buf, '\0', 80) ;
    if(var_buf[*ei][0] == 'x') {
      string s = std::string(var_buf[*ei] + 1) ;
      sprintf(tmp_buf, "%s.%s", problem_name, s.c_str()) ;
      of << "vector per node:\t" << s.c_str() << '\t' << tmp_buf << endl ;
    }
    else {
      sprintf(tmp_buf, "%s.%s", problem_name, var_buf[*ei]) ;
      of << "vector per node:\t" <<  var_buf[*ei] << '\t' << tmp_buf << endl ;
    }
    ++ei ;
  }
  ei = bvec_indices.begin() ;
  for(int i=0;i<bvector_file_tmp;++i) {
    memset(tmp_buf,'\0',80) ;
    sprintf(tmp_buf,"%s.%s",problem_name,bvar_buf[*ei]) ;
    of << "vector per element:\t" << bvar_buf[*ei]
       << '\t' << tmp_buf << endl ;
    ++ei ;
  }

  cout << "opening " << geo_name << endl ;
  FILE *OFP = fopen(geo_name, "wb") ;
  if(OFP == NULL) {
    cerr << "can't open file " << geo_name << endl ;
    exit(-1) ;
  }
  
  FILE *IFP_GRD = fopen(grid_name, "r") ;
  if(IFP_GRD == NULL) {
    cerr << "open failed on " << grid_name << endl ;
    perror(grid_name) ;
    exit(-1) ;
  }
  XDR xdr_grd ;
  xdrstdio_create(&xdr_grd, IFP_GRD, XDR_DECODE) ;
  
  xdr_int(&xdr_grd, &tmp) ;
  xdr_int(&xdr_grd, &tmp) ;
  xdr_int(&xdr_grd, &tmp) ;
  int npnts, nfaces, ncells ;
  xdr_int(&xdr_grd, &npnts) ;
  xdr_int(&xdr_grd, &nfaces) ;
  xdr_int(&xdr_grd, &ncells) ;
  xdr_int(&xdr_grd, &tmp) ;
  xdr_int(&xdr_grd, &tmp) ;
  float* tmp_pos ;
  double tmp_double ;
  tmp_pos = new float[3*npnts] ; 
  for(int j = 0; j < 3*npnts; ++j) {
    xdr_double(&xdr_grd, &tmp_double) ;
    tmp_pos[j] = float(tmp_double) ;
  }
  cout << "Finished reading the positions " << endl ;
  entitySet faces = interval(1,nfaces) ;
  entitySet cells = interval(1, ncells) ; 
  Map cl, cr ;
  multiMap face2node ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;
  store<int> count ;
  count.allocate(faces) ;
  std::vector<int> offset ;
  int *off_cl_cr ;
  off_cl_cr = new int[3*faces.size() + 1] ;
  for(int i = 0; i < (faces.size() * 3) + 1; ++i)
    xdr_int(&xdr_grd, &off_cl_cr[i]) ;
  tmp = 0 ;
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) {
    offset.push_back(off_cl_cr[tmp++]) ;
    cl[*ei] = off_cl_cr[tmp++] ;
    cr[*ei] = off_cl_cr[tmp++] ;
    if(cl[*ei] < 0) 
      if(cr[*ei] < 0) {
	cerr << " boundary condition on both sides of a face?" << endl ;
	exit(1) ;
      } else {
	int tmp_swap = cr[*ei] ;
	cr[*ei] = cl[*ei] ;
	cl[*ei] = tmp_swap ;
      }
  }
  offset.push_back(off_cl_cr[tmp]) ;
  delete [] off_cl_cr ;
  entitySet::const_iterator ii = faces.begin() ;
  for(vector<int>::size_type i = 1; i < offset.size(); ++i) {
    count[*ii] = offset[i] - offset[i-1] ;
    ++ii ;
  } 
  face2node.allocate(count) ;
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) 
    for(int i = 0; i < count[*ei]; ++i) {
      xdr_int(&xdr_grd, &face2node[*ei][i]) ;
      face2node[*ei][i] += 1 ;
      
    }
  cout << " Done with creating cl, cr, face2node " << endl ; 
  entitySet cri = Loci::MapRepP(cr.Rep())->image(faces) ;
  entitySet orig_boundaries =  cri & interval(Loci::UNIVERSE_MIN,-1) ;
  tmp = orig_boundaries.size() ;
  
  memset(tmp_buf, '\0', 80) ; 
  sprintf(tmp_buf, "C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Ensight model geometry description") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Grid file used is %s", grid_name) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "node id given") ;
  //sprintf(tmp_buf, "node id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "element id given") ;
  //sprintf(tmp_buf, "element id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  tmp = 1 ;
  fwrite(&tmp, sizeof(int), 1, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Entire") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&npnts, sizeof(int), 1, OFP) ;
  cout << "npnts = " << npnts << endl ;
  for(int i = 1 ; i <= npnts; ++i)
    fwrite(&i, sizeof(int), 1, OFP) ;
  //Coordinate information 
  for(int i = 0; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  
  for(int i = 1; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  
  for(int i = 2; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  //delete [] tmp_pos ;

  store<int> counter ;
  counter.allocate(cells) ;
  for(entitySet::const_iterator ei=cells.begin();ei!=cells.end();++ei)
    counter[*ei] = 0 ;
  
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) {
    counter[cl[*ei]]++ ;
    if(cr[*ei] > 0)
      counter[cr[*ei]]++ ;
  }
  
  multiMap cell2faces ;
  cell2faces.allocate(counter) ;
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) {
    cell2faces[cl[*ei]][--counter[cl[*ei]]] = *ei ;
    if(cr[*ei] > 0)
      cell2faces[cr[*ei]][--counter[cr[*ei]]] = *ei ;
  }
  // modified to size = 5
  int sizes[5] ;
  for(int i = 0; i < 5; ++i)
    sizes[i] = 0 ; 
  entitySet unclassified_cells ;
  int num_triangles, num_quads ;
  store<int> elem_type ;
  elem_type.allocate(cells) ;
  for(entitySet::const_iterator ei = cells.begin(); ei != cells.end(); ++ei) {
    num_triangles = 0 ;
    num_quads = 0 ;
    int num_others = 0 ;
    int triangle_nodes[3][2] ;
    for(const Entity *ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
        if(count[*ej] == 3) {
            if(num_triangles < 2) {
              triangle_nodes[0][num_triangles] = face2node[*ej][0] ;
              triangle_nodes[1][num_triangles] = face2node[*ej][1] ;
              triangle_nodes[2][num_triangles] = face2node[*ej][2] ;
            }
            num_triangles++ ;
        } else if(count[*ej] == 4)
	num_quads++ ;
      else
        num_others++ ;
    }
    bool prism_test = false ;

    if((num_triangles == 2) && (num_quads == 3) && (num_others == 0)) {
        prism_test = true ;
        for(int i=0;i<3;++i)
          for(int j=0;j<3;++j)
            if(triangle_nodes[i][0] == triangle_nodes[j][1])
              prism_test = false ;
    }

    bool hex_test = false ;
    if( (num_triangles == 0) && (num_quads == 6) && (num_others == 0)) {
        const Entity *ef = cell2faces.begin(*ei) ;
        int count = 0 ;
        for(const Entity *ej = ef+1; ej != cell2faces.end(*ei); ++ej) {
            bool find = false ;
            for(int i=0;i<4;++i)
              for(int j=0;j<4;++j)
                if(face2node[*ef][i] == face2node[*ej][j])
                  find = true ;
            if(find)
              count++ ;
        }
        if(count == 4)
          hex_test = true ;
    }
          
    
    // new classification code
    if( (num_triangles == 4) && (num_quads == 0) && (num_others == 0)) {
      sizes[0]++ ; //Increment the number of tets
      elem_type[*ei] = 0 ;
    }
    else if( hex_test ) {
      sizes[1]++ ; // Increment the number of hexs
      elem_type[*ei] = 1 ;
    }
    else if( prism_test ) {
      sizes[2]++ ; // Increment the number of prisms
      elem_type[*ei] = 2 ;
    }
    else if( (num_triangles == 4) && (num_quads == 1) && (num_others == 0)) {
      sizes[3]++ ; // Increment the number of pyramids
      elem_type[*ei] = 3 ;
    }
    else {
      sizes[4]++ ;
      elem_type[*ei] = 4 ;
    }
    
  }
  cout << " Number of tets = " << sizes[0] << endl ;
  cout << " Number of hexs = " << sizes[1] << endl ;
  cout << " Number of prisms = " << sizes[2] << endl ;
  cout << " Number of pyramids = " << sizes[3] << endl ;

  cout << " Number of nfaced cells = " << sizes[4] << endl ;
  
  MapVec<4> tet_elems ;
  tet_elems.allocate(Loci::interval(0, sizes[0]-1)) ;
  MapVec<5> pyramid_elems ;
  pyramid_elems.allocate(Loci::interval(0, sizes[3]-1)) ;
  MapVec<6> prism_elems ;
  prism_elems.allocate(Loci::interval(0, sizes[2]-1)) ;
  MapVec<8> hex_elems ;
  hex_elems.allocate(Loci::interval(0, sizes[1]-1)) ;
  vector<vector<vector<int> > > nfaced_elems(sizes[4]) ;
  std::vector<int> tet_id(sizes[0]), pyramid_id(sizes[3]), prism_id(sizes[2]), hex_id(sizes[1]), nfaced_id(sizes[4]) ;
  int tet_no = 0 ;
  int hex_no = 0 ;
  int pyramid_no = 0 ;
  int prism_no = 0 ;
  int nfaced_no = 0 ;
  for(entitySet::const_iterator ei = cells.begin(); ei != cells.end(); ++ei) {
    int quad_id[6], triangle_id[2], common[2] ;
    int quad_count = 0 , triangle_count = 0, tmp_count = 0  ;
    entitySet match_pair_1, match_pair_2, nodes_1, nodes_2, tmp_nodes ;
    const Entity *ej ;
    if(elem_type[*ei] == 0) {
      ej = cell2faces.begin(*ei) ;
      for(int i=0;i<3;++i) 
        tet_elems[tet_no][i] = face2node[*ej][i] ;
      ++ej ;
      tet_elems[tet_no][3] = -1 ;
      for(int i=0;i<3;++i) {
        bool found = true ;
        for(int j=0;j<3;++j) {
          if(tet_elems[tet_no][j] == face2node[*ej][i])
            found = false ;
        }
        if(found) {
          tet_elems[tet_no][3] = face2node[*ej][i] ;
        }
      }
      if(tet_elems[tet_no][3] == -1) {
        cerr << "unable to reconstruct tet" << endl ;
      }
      tet_id[tet_no] = *ei ;
      tet_no++ ;
    }
    else if(elem_type[*ei] == 1) {
      for(ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) 
	quad_id[quad_count++] = *ej ;
      bool degenerate = face2node[quad_id[0]][0] == face2node[quad_id[0]][3];
      for(int j=0;j<3;++j) 
        if(face2node[quad_id[0]][j] == face2node[quad_id[0]][j+1])
          degenerate = true ;
      if(degenerate) {
        for(int i=1;i<6;++i) {
          degenerate = face2node[quad_id[i]][0] == face2node[quad_id[i]][3];
          for(int j=0;j<3;++j) 
            if(face2node[quad_id[i]][j] == face2node[quad_id[i]][j+1])
              degenerate = true ;
          if(!degenerate) {
            std::swap(quad_id[i],quad_id[0]) ;
            break ;
          }
        }
      }
      hex_elems[hex_no][0] = face2node[quad_id[0]][0] ;
      hex_elems[hex_no][1] = face2node[quad_id[0]][1] ;
      hex_elems[hex_no][2] = face2node[quad_id[0]][2] ;
      hex_elems[hex_no][3] = face2node[quad_id[0]][3] ;
      for(int i = 0; i < 4; i+=2) {
        int n1 = hex_elems[hex_no][i] ;
        int n2 = hex_elems[hex_no][i+1] ;
        
        int cnt = 0 ;
        for(int j=1; j<quad_count; ++j) {
          for(int k=0;k<4;++k) {
            int f1 = face2node[quad_id[j]][k] ;
            int f2 = face2node[quad_id[j]][(k+1)%4] ;
            if((f1 == n1 && f2 == n2)) {
              hex_elems[hex_no][i+4] = face2node[quad_id[j]][(k-1+4)%4] ;
              hex_elems[hex_no][i+1+4] = face2node[quad_id[j]][(k+2)%4] ;
              cnt++ ;
            }
            if((f1 == n2 && f2 == n1)) {
              hex_elems[hex_no][i+4] = face2node[quad_id[j]][(k+2)%4] ;
              hex_elems[hex_no][i+1+4] = face2node[quad_id[j]][(k-1+4)%4] ;
              cnt++ ;
            }
          }
        }
        if(cnt != 1) {
          cerr << " Element = " << *ei << endl ;
          cerr << "Error: Hex elem ordering screwed up " <<  endl ;
          for(const Entity *ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
            for(const Entity *fc=face2node.begin(*ej);fc!=face2node.end(*ej);++fc)
              cerr << *fc << ' ' ;
            cerr << endl ;
          }
        }
      }
      hex_id[hex_no] = *ei ;
      hex_no++ ;
    } else if(elem_type[*ei] == 2) {
      for(ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
	
	if(count[*ej] == 4) {
	  quad_id[quad_count] = *ej ;  
	  quad_count++ ;
	} else if(count[*ej] == 3) {
	  triangle_id[triangle_count] = *ej ;
	  triangle_count++ ;
	} 
      } 
      prism_elems[prism_no][0] = face2node[triangle_id[0]][0] ;
      prism_elems[prism_no][1] = face2node[triangle_id[0]][1] ;
      prism_elems[prism_no][2] = face2node[triangle_id[0]][2] ;
      
      for(int i = 0; i < 3; ++i) {
	tmp_count = 0 ;
	for(int j = 0; j < quad_count; ++j) {
	  for(int k = 0; k < 4; ++k)
	    if(prism_elems[prism_no][i] == face2node[quad_id[j]][k]) 
	      common[tmp_count++] = quad_id[j] ;
	}
	nodes_1 = EMPTY ;
	nodes_2 = EMPTY ;
	for(int k = 0; k < 4; ++k)
	  nodes_1 += face2node[common[0]][k] ;
	for(int k = 0; k < 4; ++k)
	  nodes_2 += face2node[common[1]][k] ;
	match_pair_2 = nodes_1 & nodes_2 ;
	match_pair_2 -= interval(prism_elems[prism_no][i], prism_elems[prism_no][i]) ;
	if(match_pair_2.size() != 1) {
            cerr << " Element = " << *ei << endl ;
            cerr << "match_pair_2.size() = " << match_pair_2.size() << endl ;
            cerr << "Error: Prism elem ordering screwed up " <<  endl ;
            for(const Entity *ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
                for(const Entity *fc=face2node.begin(*ej);fc!=face2node.end(*ej);++fc)
                  cerr << *fc << ' ' ;
                cerr << endl ;
            }
            cerr << endl ;
            
	}
	prism_elems[prism_no][i+3] = *match_pair_2.begin() ;
      }
      
      prism_id[prism_no] = *ei ;
      prism_no++ ;
    } else if(elem_type[*ei] == 3) {
      entitySet total_nodes ;
      entitySet final_nodes ;

      for(ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
	if(count[*ej] == 4)
	  quad_id[0] = *ej ;
	for(int i = 0; i < count[*ej]; ++i)
	  total_nodes += face2node[*ej][i] ;
      }
      
      for(int i = 0; i < count[quad_id[0]]; ++i) {
	pyramid_elems[pyramid_no][i] = face2node[quad_id[0]][i] ;
	final_nodes += pyramid_elems[pyramid_no][i] ;
      }  

      total_nodes -= final_nodes ;
      if(total_nodes.size() != 1)
	cerr << "Pyramid elements ordering screwed up " << endl ;
      pyramid_elems[pyramid_no][4] = *total_nodes.begin() ;
      pyramid_id[pyramid_no] = *ei ;
      pyramid_no++ ;
    }else if(elem_type[*ei] == 4) {
      vector<vector<int> > faces ;
      for(ej = cell2faces.begin(*ei); ej != cell2faces.end(*ei); ++ej) {
        vector<int> facenode ;
	for(int i = 0; i < count[*ej]; ++i)
	  facenode.push_back(face2node[*ej][i]) ;
        faces.push_back(facenode) ;
      }
      nfaced_elems[nfaced_no] = faces ;
      ++nfaced_no ;
    }else {
      cerr <<  "ERROR:  Unknown element type " << endl ;
      return 0;	
    }
  }
  if(tet_no) {
    entitySet nodeSet_temp ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "tetra4") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&tet_no, sizeof(int), 1, OFP) ;
    for(int i = 0; i < tet_no; ++i)
      fwrite(&tet_id[i], sizeof(int), 1, OFP) ;
    for(int i = 0; i < tet_no; ++i)
      fwrite(tet_elems[i], sizeof(int), 4, OFP) ;
  }
  if(hex_no) {
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "hexa8") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&hex_no, sizeof(int), 1, OFP) ;
    for(int i = 0; i < hex_no; ++i)
      fwrite(&hex_id[i], sizeof(int), 1, OFP) ;
    for(int i = 0; i < hex_no; ++i)
      fwrite(hex_elems[i], sizeof(int), 8, OFP) ;
  }
  if(prism_no) {
    sprintf(tmp_buf, "penta6") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&prism_no, sizeof(int), 1, OFP) ;
    for(int i = 0; i < prism_no; ++i)
      fwrite(&prism_id[i], sizeof(int), 1, OFP) ;
    for(int i = 0; i < prism_no; ++i)
      fwrite(prism_elems[i], sizeof(int), 6, OFP) ;
  }
  if(pyramid_no) {
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "pyramid5") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&pyramid_no, sizeof(int), 1, OFP) ;
    for(int i = 0; i < pyramid_no; ++i)
      fwrite(&pyramid_id[i], sizeof(int), 1, OFP) ;
    for(int i = 0; i < pyramid_no; ++i)
      fwrite(pyramid_elems[i], sizeof(int), 5, OFP) ;
  }
  if(nfaced_no > 0) {
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "nfaced") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&nfaced_no, sizeof(int), 1, OFP) ;
    // writing out element id
    for(int i = 0; i < nfaced_no; ++i)
      fwrite(&nfaced_id[i], sizeof(int), 1, OFP) ;
    // writing out face number for each element
    for(int i = 0; i < nfaced_no; ++i) {
      int facenum = nfaced_elems[i].size() ;
      fwrite(&facenum, sizeof(int), 1, OFP) ;
    }
    // writing out node number for each face in each element
    for(int i = 0; i < nfaced_no; ++i) {
      vector<vector<int> >& faces = nfaced_elems[i] ;
      int size = faces.size() ;
      for(int j=0;j<size;++j) {
        int nodenum = faces[j].size() ;
        fwrite(&nodenum, sizeof(int), 1, OFP) ;
      }
    }
    // finally write out the node indices for each face in each element
    for(int i = 0; i < nfaced_no; ++i) {
      vector<vector<int> >& faces = nfaced_elems[i] ;
      int size = faces.size() ;
      for(int j=0;j<size;++j) {
        vector<int>& facenode = faces[j] ;
        int nodesize = facenode.size() ;
        for(int k=0;k<nodesize;++k) {
          int nodeindex = facenode[k] ;
          fwrite(&nodeindex, sizeof(int), 1, OFP) ;
        }
      }
    }
    
  }

  // writing out the boundary information
  entitySet crcell = Loci::MapRepP(cr.Rep())->image(faces) ;
  entitySet rev_boundaries =  crcell & interval(Loci::UNIVERSE_MIN,-1) ;
  std::set<int,std::greater<int> > boundaries ;
  for(entitySet::const_iterator ei=rev_boundaries.begin();
      ei!=rev_boundaries.end();++ei)
    boundaries.insert(*ei) ;

  int part_id = 2 ;
  map<int,entitySet> boundary_face_map ;
  std::map<int,vector<int> > part_node_table ; // used in variable output
  for(std::set<int,std::greater<int> >::const_iterator bi=boundaries.begin();
      bi!=boundaries.end();++bi) {
    // get the bc id
    std::stringstream ss ;
    ss << -(*bi) ;
    std::string s = ss.str() ;
    // first compute those faces that are on this boundary
    entitySet bnd ; bnd += *bi ;
    std::pair<entitySet,entitySet> pimage =
      Loci::MapRepP(cr.Rep())->preimage(bnd) ;
    entitySet boundary_faces = pimage.first ;
    boundary_face_map[part_id] = boundary_faces ;
    int face_num = boundary_faces.size() ;
    if(face_num == 0) continue ;
    
    entitySet nodes ;
    Map face_node_num ;
    face_node_num.allocate(boundary_faces) ;
    for(entitySet::const_iterator fi=boundary_faces.begin();
        fi!=boundary_faces.end();++fi) {
      face_node_num[*fi] = count[*fi] ;
      for(int i=0;i<count[*fi];++i)
        nodes += face2node[*fi][i] ;
    }
    int node_num = nodes.size() ;

    Map node_offset ;
    node_offset.allocate(nodes) ;
    int nof = 1 ;
    vector<int> node_local ;
    for(entitySet::const_iterator ei=nodes.begin();
        ei!=nodes.end();++ei) {
      node_offset[*ei] = nof++ ;
      node_local.push_back(*ei) ;
    }
    part_node_table[part_id] = node_local ;
    
    std::map<int,vector<int> > face_node_id ;
    for(entitySet::const_iterator fi=boundary_faces.begin();
        fi!=boundary_faces.end();++fi) {
      vector<int> tmp_nodes ;
      for(int i=0;i<count[*fi];++i) {
        int n = face2node[*fi][i] ;
        tmp_nodes.push_back(node_offset[n]) ;
      }
      face_node_id[*fi] = tmp_nodes ;
    }

    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&part_id, sizeof(int), 1, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "BC_%s", s.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&node_num, sizeof(int), 1, OFP) ;
    cout << "number of nodes for BC_" << s << " = " << node_num << endl ;
    cout << "number of faces for BC_" << s << " = " << face_num << endl ;
    for(entitySet::const_iterator ni=nodes.begin();
        ni!=nodes.end(); ++ni) {
      int nid = *ni ;
      fwrite(&nid, sizeof(int), 1, OFP) ;
    }
    //Coordinate information
    // the index is offset by 1 since node is numbered from 1
    for(entitySet::const_iterator ni=nodes.begin();
        ni!=nodes.end();++ni) {
      int index = 3*(*ni-1) ;
      fwrite(&tmp_pos[index],sizeof(float),1,OFP) ;
    }
    for(entitySet::const_iterator ni=nodes.begin();
        ni!=nodes.end();++ni) {
      int index = 3*(*ni-1)+1 ;
      fwrite(&tmp_pos[index],sizeof(float),1,OFP) ;
    }
    for(entitySet::const_iterator ni=nodes.begin();
        ni!=nodes.end();++ni) {
      int index = 3*(*ni-1)+2 ;
      fwrite(&tmp_pos[index],sizeof(float),1,OFP) ;
    }
    // fill in boundary face info
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "nsided") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&face_num, sizeof(int), 1, OFP) ;
    for(entitySet::const_iterator fi=boundary_faces.begin();
        fi!=boundary_faces.end();++fi) {
      int face_no = *fi ;
      fwrite(&face_no,sizeof(int),1,OFP) ;
    }
    for(entitySet::const_iterator fi=boundary_faces.begin();
        fi!=boundary_faces.end();++fi) {
      //int node_no = face_node_num[*fi] ;
      fwrite(&face_node_num[*fi],sizeof(int),1,OFP) ;
    }
    for(entitySet::const_iterator fi=boundary_faces.begin();
        fi!=boundary_faces.end();++fi) {
      std::map<int,vector<int> >::const_iterator mi ;
      mi = face_node_id.find(*fi) ;
      FATAL(mi == face_node_id.end()) ;
      for(vector<int>::const_iterator vi=mi->second.begin();
          vi!=mi->second.end();++vi) {
        int node_of = *vi ;
        fwrite(&node_of,sizeof(int),1,OFP) ;
      }
    }

    ++part_id ;
  }
  delete [] tmp_pos ;
  /////////////
  fclose(IFP_GRD) ;
  fclose(OFP) ;
  OFP = NULL ;
  // read qf file here
  vector<qf_boundary> qf_boundary_vec ;
  map<int,qf_boundary> matched_qf_boundary ;
  entitySet matched_faces ;
  vector<int> qf_exist_part ;
  if(read_boundary_value) {
    cout << " Reading in boundary values " << endl ;

    char qf_name[512] ;
    sprintf(qf_name, "output/qf.%s", ncyc) ; 
    ifstream ifile(qf_name, ios::in) ;
    if(ifile.fail())  {
      cerr << " Fail to open file: " << qf_name << endl ;
      exit(-1) ;
    }
    int read_count = 0 ;
    while(!ifile.eof()) {
      read_count++ ;
      qf_boundary a_bound ;
      int dummy ;
      double d ;
      int invalid_id = 0 ;
      for(int i = 0; i < 4; ++i) {
        ifile >> dummy ;
        if(ifile.eof()) {
          cerr << qf_name << " format error." << endl ;
          exit(-1) ;
        }
        dummy++ ;
        if(dummy != 0)
          a_bound.facenode.push_back(dummy) ;
        else {
          if(++invalid_id > 1) {
            cerr << qf_name << " format error "
                 << "(a face must have at least 3 nodes defined)."
                 << endl ;
            exit(-1) ;
          }
        }
      }
      // std::sort(a_bound.facenode.begin(),a_bound.facenode.end()) ;
      // variable orders in the qf file
      //ofile << xc << " " << yc << " " << zc << " " 
      //    << tx << " " << ty << " " << tz << " " 
      //    << qf << " " << p << ' ' << T  << ' '<< yp_w << endl ;
      // skip the next three value
      for(int i = 0; i < 3; ++i) {
        ifile >> d ;
        if(ifile.eof()) {
          cerr << qf_name << " format error." << endl ;
          exit(-1) ;
        }
      }
      // read in the rest values
      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.w[0] = float(d) ;
      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.w[1] = float(d) ;
      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.w[2] = float(d) ;
      
      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.q = float(d) ;

      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.wp = float(d) ;

      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.wt = float(d) ;

      ifile >> d ;
      if(ifile.eof()) {
        cerr << qf_name << " format error." << endl ;
        exit(-1) ;
      }
      a_bound.y = float(d) ;

      // add to the vector collection
      qf_boundary_vec.push_back(a_bound) ;
      // consume the rest of the line
      while(!ifile.eof()) {
        int c = ifile.peek();
        if(!isspace(c) && (c!=EOF))
          break ;
        c = ifile.get() ;
      }
    } // end while(!ifile.eof())
    ifile.close() ;
    cout << "Finished reading in boundary values " << endl ;
    cout << "Number of values read from the qf file = "
         << read_count << endl ;

    cout << " Matching boundary value to faces..." << endl ;

    // first get all the boundary faces
    entitySet crcell = Loci::MapRepP(cr.Rep())->image(faces) ;
    entitySet boundaries = crcell & interval(Loci::UNIVERSE_MIN,-1) ;

    std::pair<entitySet,entitySet> pimage =
      Loci::MapRepP(cr.Rep())->preimage(boundaries) ;
    entitySet boundary_faces = pimage.first ;

    // then build a node2face map
    map<int,entitySet> node2face ;
    for(entitySet::const_iterator ei=face2node.domain().begin();
        ei!=face2node.domain().end();++ei) {
      int sz = face2node.end(*ei) - face2node.begin(*ei) ;
      for(int i=0;i<sz;++i)
	node2face[face2node[*ei][i]] += *ei ;
    }

    // then we begin matching
    int match_count = 1 ;
    int unmatched_number = 0 ;
    //timeval t1,t2 ;
    //gettimeofday(&t1,NULL) ;
    for(vector<qf_boundary>::const_iterator vi=qf_boundary_vec.begin();
        vi!=qf_boundary_vec.end();++vi,match_count++) {
      const vector<int>& facenode = (vi->facenode) ;
      entitySet face = interval(Loci::UNIVERSE_MIN,Loci::UNIVERSE_MAX) ;
      // we only need to intersect 3 times
      for(int ni=0;ni<3;++ni) {
        map<int,entitySet>::const_iterator mi ;
        mi = node2face.find(facenode[ni]) ;
        if(mi==node2face.end()) {
          cerr << "Data inconsistent, matching failed "
               << "(in lookup node infomation.)." << endl ;
          exit(-1) ;
        }
        face &= mi->second ;
      }
      // by now, face should only contain one entity
      // and that is the match result
      if(face.size() != 1) {
        cerr << "Data inconsistent, matching failed "
             << "(in the " << match_count << "th matching, "
             << "face.size() = " << face.size() << ")." << endl ;
        unmatched_number++ ;
        //exit(-1) ;
      } else {
        // then we store the result
        int face_id = *(face.begin()) ;
        matched_faces += face_id ;
        matched_qf_boundary[face_id] = *vi ;
      }
    }
    //gettimeofday(&t2,NULL) ;
    // matching finished, we no longer need the
    // memory allocated for qf_boundary_vec, so
    // we clean it to release some resourse.
    qf_boundary_vec.clear() ;

    // now we need to calculate which boundary faces
    // are matched in which part
    if(unmatched_number > 0)
      cerr << "There are " << unmatched_number
           << " faces unmatched." << endl ;
    if( (matched_faces - boundary_faces) != EMPTY)
      cerr << "Some faces in " << qf_name
           << " are not on boundaries." << endl ;

    entitySet matched_faces_copy = matched_faces ;
    for(map<int,entitySet>::const_iterator mi=boundary_face_map.begin();
        mi!=boundary_face_map.end();++mi) {
      entitySet common = (mi->second) & matched_faces_copy ;
      entitySet diff = (mi->second) - matched_faces_copy ;
      if(common == EMPTY)
        continue ;
      qf_exist_part.push_back(mi->first) ;
      matched_faces_copy -= common ;
      // fill in what's missing with 0.0
      if(diff != EMPTY) {
        cerr << " Not all values for BC are defined." << endl ;
      }
      qf_boundary zero_bound ;
      for(entitySet::const_iterator ei=diff.begin();
          ei!=diff.end();++ei)
        matched_qf_boundary[*ei] = zero_bound ;
    }
    cout << " Finished matching boundary values." << endl ;
  } // end if(read_boundary_value)
  
  sprintf(filename, "output/qn_xdr.%s", ncyc) ;
  int t = 0 ;
  for(int n = 0; n < num_vars; ++n) {
    float* val_vec = 0 ; // initialize it
    int type = *var_buf[n] ;
    int count = 0 ;
    memset(tmp_buf, '\0', 80) ;
    if(scalar_indices.inSet(n)) {
      if(var_buf[n][0] == 'f') {
	std::string s = string(var_buf[n] + 1) ;
	cout << "Writing out scalar variable " << s << "from qn file " <<  endl ;
      } else if(var_buf[n][0] == 'n') {
	std::string s = string(var_buf[n] + 1) ;
	extra_vars = 1 ;
	var_map[count++] = n ;
      } else
        if(var_buf[n][0] == 'x') {
	std::string s = string(var_buf[n] + 1) ;
	cout << "Writing out additional vector variable " << s << endl ;
	extra_vec_vars = 1 ;
	var_map[count++] = n ;
      } 

      cout << "opening scalar var file " << scalar_var_files[t] << endl ;
      OFP = fopen(scalar_var_files[t++], "wb") ;
      if(OFP == NULL) {
	cerr << "can't open file " << scalar_var_files[t-1] << endl ;
	exit(-1) ;
      }
      sprintf(tmp_buf, "Variable file %s", scalar_var_files[t-1]) ;
    }
    else {
      val_vec = new float[3*npnts] ;
      if(var_buf[n][0] == 'x') {
          std::string s = string(var_buf[n] + 1) ;
          extra_vec_vars = 1 ;
          var_map[count++] = n ;

          cout << "opening additional vector var file "
               << vector_var_files[vec_map[n]] << endl ;
          
          OFP = fopen(vector_var_files[vec_map[n]], "wb") ;
          if(OFP == NULL) {
              cerr << "can't open file " << vector_var_files[vec_map[n]] << endl ;
              exit(-1) ;
          }
      } else {
          cout << "opening vector var file " <<
              vector_var_files[vec_map[n]] << endl ;

          OFP = fopen(vector_var_files[vec_map[n]], "wb") ;
          if(OFP == NULL) {
              cerr << "can't open file " << vector_var_files[vec_map[n]] << endl ;
              exit(-1) ;
          }
      }
      sprintf(tmp_buf, "Variable file %s", vector_var_files[vec_map[n]]) ;

    }
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    tmp = 1 ;
    fwrite(&tmp, sizeof(int), 1, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;

    if(!extra_vars && !extra_vec_vars) {
      FILE *IFP_QN = fopen(filename, "r") ;
      if(IFP_QN == NULL) {
	cerr << "open failed on " << filename << endl ;
	exit(-1) ;
      }
      XDR xdr_qn ;
      xdrstdio_create(&xdr_qn, IFP_QN, XDR_DECODE) ;
      int ns ;
      char tmp_char[256] ;
      int tmp_size = 0 ;
      xdr_int(&xdr_qn, &ns) ;
      int numvars = 7 + ns ;
      int nv = numvars ;
      string *species = new string[ns];
      for(int i = 0; i < ns; ++i) {
	xdr_int(&xdr_qn, &tmp_size) ;
	for(int j = 0; j < tmp_size; ++j)
	  xdr_char(&xdr_qn, &tmp_char[j]) ;
	for(int k = 0; k < tmp_size; ++k)
	  species[i] += tmp_char[k] ;
      }
      int nnodes = 0;
      xdr_int(&xdr_qn, &nnodes) ;
      if(nnodes != npnts) {
	cerr << " Error : number of nodes from the output file and the grid file differ " << endl ;
	exit(-1) ;
      }
      int squery =- 1 ;
      if(var_buf[n][0] == 'f') {
	string s = string(var_buf[n]+1) ;
	for(int i = 0; i < ns; ++i)
	  if(species[i] == s) {
	    squery = i ;
	    cout << " species[ " << i << " ] = " << species[i] << " s = " << s <<  endl ;
	  }
	if(squery == -1) {
	  cerr << "species " << s << " not in data " << endl ;
	  exit(-1) ;
	}
      }
      double *qb = new double[nv] ;
      tmp = 0 ;
      float* scalar_node_val = new float[nnodes] ;
      for(int i = 0; i < nnodes; ++i) {
	for(int ii = 0; ii < nv; ++ii) 
	  xdr_double(&xdr_qn, &qb[ii]) ;
	
	double r =  qb[0] ;
	double u =  qb[1] ;
	double v =  qb[2] ;
	double w =  qb[3] ;
	double a =  qb[4] ;
	double T =  qb[5] ;
	double P =  qb[6] ;
	double U = sqrt(u*u+v*v+w*w) ;
	double M = U/a ;
        double Gamma = a*a*r/P ;
        double Rt = P/(r*T) ;
        
	int t ;
	float val = 0.0 ;
	if(scalar_indices.inSet(n)) {
	  switch(type) {
	  case 'r':
	    val = r ;
	    break ;
	  case 'P':
	    val = P ;
	    break ;
	  case 'p':
	    val = log10(P)  ;
	    break ;
	  case 'U':
	    val = U ;
	    break ;
	  case 'm':
	    val = M ;
	    break ;
          case 'a':
            val = a ;
            break ;
          case 'G':
            val = Gamma ;
            break ;
          case 'R':
            val = Rt ;
            break ;
	  case 't':
	    val = T ;
	    break ;
	  case 'f':
	    t = squery+7 ;
	    val =  qb[t] ;
	    break ;	
	  default:
	    cerr << "wrong type"<< endl ;
	    exit(0) ;
	  }
          scalar_node_val[i] = val ;
          fwrite(&val, sizeof(float), 1, OFP) ;
	}
	else {
	  switch(type) {
	  case 'v':
	    val_vec[tmp++] = u ;
	    val_vec[tmp++] = v ;
	    val_vec[tmp++] = w ;
	    break ;
	  }
	}
      }
      if(!scalar_indices.inSet(n)) {
	for(int i = 0; i < 3*npnts; i+=3) 
	  fwrite(&val_vec[i], sizeof(float), 1, OFP) ;
	
	for(int i = 1; i < 3*npnts; i+=3) 
	  fwrite(&val_vec[i], sizeof(float), 1, OFP) ;
	
	for(int i = 2; i < 3*npnts; i+=3) 
	  fwrite(&val_vec[i], sizeof(float), 1, OFP) ;
	//delete [] val_vec ;
      }
      // write out the boundary part
      for(int part_num=2;part_num<part_id;++part_num) {
        memset(tmp_buf, '\0', 80) ;
        sprintf(tmp_buf, "part") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
        fwrite(&part_num, sizeof(int), 1, OFP) ;
        memset(tmp_buf, '\0', 80) ;
        sprintf(tmp_buf, "coordinates") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
        // get all the node in this part
        std::map<int,vector<int> >::const_iterator mi ;
        mi = part_node_table.find(part_num) ;
        FATAL(mi == part_node_table.end()) ;
        const vector<int>& part_node = mi->second ;
        if(scalar_indices.inSet(n)) {
          for(vector<int>::const_iterator vi=part_node.begin();
              vi!=part_node.end();++vi) {
            // we need to offset by 1 since node is numbered from 1
            fwrite(&scalar_node_val[*vi-1],sizeof(float),1,OFP) ;
          }
        }else {
          vector<int>::const_iterator vi ;
          // index offset by 1 since node is numbered from 1
          for(vi=part_node.begin();vi!=part_node.end();++vi) {
            int index = 3*(*vi-1) ;
            fwrite(&val_vec[index],sizeof(float),1,OFP) ;
          }
          for(vi=part_node.begin();vi!=part_node.end();++vi) {
            int index = 3*(*vi-1)+1 ;
            fwrite(&val_vec[index],sizeof(float),1,OFP) ;
          }
          for(vi=part_node.begin();vi!=part_node.end();++vi) {
            int index = 3*(*vi-1)+2 ;
            fwrite(&val_vec[index],sizeof(float),1,OFP) ;
          }
        }
      }
      delete [] scalar_node_val ;
      delete [] val_vec ;
      delete [] qb ;
      xdr_destroy(&xdr_qn) ;
      fclose(IFP_QN) ;
      fclose(OFP) ;
      OFP = NULL ;
    }
    else if(extra_vars) {
      char tmp_name[512] ;
      char *tmp_char = var_buf[var_map[count-1]] + 1 ;
      sprintf(tmp_name, "output/%s_hdf5.%s",tmp_char, ncyc) ;
      if(stat(tmp_name,&statbuf)==0) {
        cout << "opening file " << tmp_name << endl ;
        store<float> read_scalar ;
        hid_t  file_id, group_id;
        file_id = H5Fopen(tmp_name, H5F_ACC_RDONLY, H5P_DEFAULT);
        group_id = H5Gopen(file_id, tmp_char) ;
        //        entitySet dom = Loci::interval(0, npnts-1) ;
        entitySet dom = EMPTY ;
        Loci::read_container(group_id, read_scalar.Rep(), dom) ;
        float val ;
        float* scalar_node_val = new float[npnts] ;
        if(dom.size() != npnts) {
          cerr << " domains don't match in size for variable " << tmp_name
               << endl ;
        }
        entitySet::const_iterator ei = dom.begin() ;
        for(int i = 0; i < npnts; ++i,++ei) {
          val = read_scalar[*ei] ;
          scalar_node_val[i] = val ;
          fwrite(&val, sizeof(float), 1, OFP) ;
        }

        // write out the boundary part
        for(int part_num=2;part_num<part_id;++part_num) {
          memset(tmp_buf, '\0', 80) ;
          sprintf(tmp_buf, "part") ;
          fwrite(tmp_buf, sizeof(char), 80, OFP) ;
          fwrite(&part_num, sizeof(int), 1, OFP) ;
          memset(tmp_buf, '\0', 80) ;
          sprintf(tmp_buf, "coordinates") ;
          fwrite(tmp_buf, sizeof(char), 80, OFP) ;
          // get all the node in this part
          std::map<int,vector<int> >::const_iterator mi ;
          mi = part_node_table.find(part_num) ;
          FATAL(mi == part_node_table.end()) ;
          const vector<int>& part_node = mi->second ;
          for(vector<int>::const_iterator vi=part_node.begin();
              vi!=part_node.end();++vi) {
            // we need to offset by 1 since node is numbered from 1
            fwrite(&scalar_node_val[*vi-1],sizeof(float),1,OFP) ;
          }
        }
        delete [] scalar_node_val ;
        
        extra_vars = 0 ;
        H5Gclose(group_id) ;
        H5Fclose(file_id) ;
        fclose(OFP) ;
        OFP = NULL ;
      } else {
        sprintf(tmp_name, "output/%s_%s.xdr",tmp_char, ncyc) ;
        cout << "opening file " << tmp_name << endl ;
        FILE *IFP_QN = fopen(tmp_name, "r") ;
        if(IFP_QN == NULL) {
          cerr << "open failed on " << tmp_name << endl ;
          exit(-1) ;
        }
        XDR xdr_qn ;
        xdrstdio_create(&xdr_qn, IFP_QN, XDR_DECODE) ;
        int nnodes  =0;
        xdr_int(&xdr_qn, &nnodes) ;
        if(nnodes != npnts) {
          cerr << " Error : number of nodes from the output file and the grid file differ " << endl ;
          exit(-1) ;
        }
        double tmp_double ;
        float val ;
        float* scalar_node_val = new float[nnodes] ;
        for(int i=0;i<nnodes;++i) {
          xdr_double(&xdr_qn, &tmp_double) ;
          val = float(tmp_double) ;
          scalar_node_val[i] = val ;
          fwrite(&val, sizeof(float), 1, OFP) ;
        }
        // write out the boundary part
        for(int part_num=2;part_num<part_id;++part_num) {
          memset(tmp_buf, '\0', 80) ;
          sprintf(tmp_buf, "part") ;
          fwrite(tmp_buf, sizeof(char), 80, OFP) ;
          fwrite(&part_num, sizeof(int), 1, OFP) ;
          memset(tmp_buf, '\0', 80) ;
          sprintf(tmp_buf, "coordinates") ;
          fwrite(tmp_buf, sizeof(char), 80, OFP) ;
          // get all the node in this part
          std::map<int,vector<int> >::const_iterator mi ;
          mi = part_node_table.find(part_num) ;
          FATAL(mi == part_node_table.end()) ;
          const vector<int>& part_node = mi->second ;
          for(vector<int>::const_iterator vi=part_node.begin();
              vi!=part_node.end();++vi) {
            // we need to offset by 1 since node is numbered from 1
            fwrite(&scalar_node_val[*vi-1],sizeof(float),1,OFP) ;
          }
        }
        delete [] scalar_node_val ;

        xdr_destroy(&xdr_qn) ;
        fclose(IFP_QN) ;
      }
    }
    else if(extra_vec_vars) {
      char tmp_name[512] ;
      char *tmp_char = var_buf[var_map[count-1]] + 1 ;
      sprintf(tmp_name, "output/%s_hdf5.%s",tmp_char, ncyc) ;
      cout << "opening file " << tmp_name << endl ;
      store<vector3d<float> > read_vector ;
      hid_t  file_id, group_id;
      file_id = H5Fopen(tmp_name, H5F_ACC_RDONLY, H5P_DEFAULT);
      group_id = H5Gopen(file_id, tmp_char) ;
      entitySet dom = EMPTY ;
      Loci::read_container(group_id, read_vector.Rep(), dom) ;
      if(dom.size() != npnts) {
        cerr << "domains don't match in size for variable " << tmp_name
             << endl ;
      }
      float val[3] ;
      float* vec_node_val = new float[3*npnts] ;
      entitySet::const_iterator ei = dom.begin() ;
      for(int i = 0; i < npnts; ++i,++ei) {
        val[0] = read_vector[*ei].x ;
        val[1] = read_vector[*ei].y ;
        val[2] = read_vector[*ei].z ;
             
	for(int j = 0; j < 3; ++j) {
          vec_node_val[3*i+j] = val[j] ;
        }
      }
      for(int i = 0; i < 3*npnts; i+=3) 
        fwrite(&vec_node_val[i], sizeof(float), 1, OFP) ;
	
      for(int i = 1; i < 3*npnts; i+=3) 
        fwrite(&vec_node_val[i], sizeof(float), 1, OFP) ;
      
      for(int i = 2; i < 3*npnts; i+=3) 
        fwrite(&vec_node_val[i], sizeof(float), 1, OFP) ;
      // write out the boundary part
      for(int part_num=2;part_num<part_id;++part_num) {
        memset(tmp_buf, '\0', 80) ;
        sprintf(tmp_buf, "part") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
        fwrite(&part_num, sizeof(int), 1, OFP) ;
        memset(tmp_buf, '\0', 80) ;
        sprintf(tmp_buf, "coordinates") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
        // get all the node in this part
        std::map<int,vector<int> >::const_iterator mi ;
        mi = part_node_table.find(part_num) ;
        FATAL(mi == part_node_table.end()) ;
        const vector<int>& part_node = mi->second ;
        vector<int>::const_iterator vi ;
        // index offset by 1 since node is numbered from 1
        for(vi=part_node.begin();vi!=part_node.end();++vi) {
          int index = 3*(*vi-1) ;
          fwrite(&vec_node_val[index],sizeof(float),1,OFP) ;
        }
        for(vi=part_node.begin();vi!=part_node.end();++vi) {
          int index = 3*(*vi-1)+1 ;
          fwrite(&vec_node_val[index],sizeof(float),1,OFP) ;
        }
        for(vi=part_node.begin();vi!=part_node.end();++vi) {
          int index = 3*(*vi-1)+2 ;
          fwrite(&vec_node_val[index],sizeof(float),1,OFP) ;
        }
      }
      delete [] vec_node_val ;
      
      extra_vec_vars = 0 ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
      fclose(OFP) ;
      OFP = NULL ;
    } 
  }
  // finally write out boundary info in qf. file
  int bscalar_file_index = 0 ;
  int bvector_file_index = 0 ;
  for(int n = 0; n < num_bvars; ++n) {
    memset(tmp_buf, '\0', 80) ;
    if(bscalar_indices.inSet(n)) {
      cout << "Writing out variable " << bvar_buf[n]
           << " from scalar file " <<  endl ;
      cout << "opening scalar boundary var file " << bscalar_var_files[bscalar_file_index] << endl ;
      OFP = fopen(bscalar_var_files[bscalar_file_index], "wb") ;
      if(OFP == NULL) {
	cerr << "can't open file "
             << bscalar_var_files[bscalar_file_index]
             << endl ;
	exit(-1) ;
      }
      sprintf(tmp_buf, "Variable file %s",
              bscalar_var_files[bscalar_file_index]) ;
      bscalar_file_index++ ;
    }
    else {
        cout << "Writing out vector variable " << bvar_buf[n] << endl ;
        cout << "opening boundary vector variable file "
             << bvector_var_files[bvector_file_index] << endl ;
      OFP = fopen(bvector_var_files[bvector_file_index], "wb") ;
      if(OFP == NULL) {
	cerr << "can't open file "
             << bvector_var_files[bvector_file_index]
             << endl ;
	exit(-1) ;
      }
      sprintf(tmp_buf, "Variable file %s",
              bvector_var_files[bvector_file_index]) ;
      bvector_file_index++ ;
    }
    // write out 
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    // we only need to write out the needed part
    for(vector<int>::const_iterator vi=qf_exist_part.begin();
        vi!=qf_exist_part.end();++vi) {
      sprintf(tmp_buf, "part") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int part_id = *vi ;
      fwrite(&part_id, sizeof(int), 1, OFP) ;
      memset(tmp_buf, '\0', 80) ;
      sprintf(tmp_buf, "nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      // retrieve the boundary faces
      map<int,entitySet>::const_iterator mi ;
      mi=boundary_face_map.find(part_id) ;
      if(mi==boundary_face_map.end()) {
        cerr << "Data inconsistent." << endl ;
        exit(-1) ;
      }
      const entitySet& bfaces = mi->second ;
      entitySet::const_iterator ei ;
      if(!strcmp(bvar_buf[n],"q")) {
        // write out q value
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.q),sizeof(float),1,OFP) ;
        }
      } else if(!strcmp(bvar_buf[n],"y")) {
        // write out y value
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.y),sizeof(float),1,OFP) ;
        }
      } else if(!strcmp(bvar_buf[n],"wp")) {
        // write out wp value
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.wp),sizeof(float),1,OFP) ;
        }
      } else if(!strcmp(bvar_buf[n],"wt")) {
        // write out wt value
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.wt),sizeof(float),1,OFP) ;
        }
      } else if(!strcmp(bvar_buf[n],"w")) {
        // write out the vector value "w"
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.w[0]),sizeof(float),1,OFP) ;
        }
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.w[1]),sizeof(float),1,OFP) ;
        }
        for(ei=bfaces.begin();ei!=bfaces.end();++ei) {
          const qf_boundary& qfb = matched_qf_boundary[*ei] ;
          fwrite(&(qfb.w[2]),sizeof(float),1,OFP) ;
        }
      } else {
        cerr << "Could not extract variable: "
             << bvar_buf[n] << endl ;
        exit(-1) ;
      }
    }
    fclose(OFP) ;
    OFP = NULL ;
  }
    
  cout << "Done " << endl ;
  return 0 ;
}

 

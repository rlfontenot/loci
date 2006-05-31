#include <Loci.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
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
#include "fv_extract.h"
#include <vector>
using std::vector ;

/*
** fv_encode_elem_header:  return an encoded binary element header
**
** Input:
**    elem_type:  integer element type as shown in fv_reader_tags.h
**    wall_info:  array of integer "wall" flags, one for each face of
**                the element.  The wall flags are used during streamline
**                calculation.  Currently, the only meaningful values are
**                A_WALL and NOT_A_WALL as shown in fv_reader_tags.h.
**                Streamlines are forced away from faces marked as
**                "A_WALL", by limiting velocity and position very near
**                the wall.
** Output:
**    Function return value is the encoded binary element header.
*/

#ifdef __STDC__
unsigned int fv_encode_elem_header (int elem_type, int wall_info[])
#else
unsigned int fv_encode_elem_header (elem_type, wall_info)
int elem_type;
int wall_info[];
#endif
{
  unsigned int header;
  int i, nfaces;
  switch (elem_type)
    {
    case FV_TET_ELEM_ID:
      header = (1 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 4;
            break;
    case FV_HEX_ELEM_ID:
      header = (4 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 6;
      break;
    case FV_PRISM_ELEM_ID:
      header = (3 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 5;
      break;
    case FV_PYRA_ELEM_ID:
      header = (2 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 5;
      break;
    default:
      fprintf(stderr, "ERROR:  Unknown element type\n");
      return 0;
    }
  
  for (i = 0; i < nfaces; i++)
    {
      unsigned int u = wall_info[i];
      if (u > A_WALL)
        {
	  fprintf(stderr, "ERROR:  Bad wall value\n");
	  return 0;
        }
      header |= (u << (i*BITS_PER_WALL));
    }
  return header;
}

void FVUsage(char *s) {
  cerr << "Usage: " << endl 
       << string(s) << " -fv <problem_name> <time_step> <key> " << endl << endl ;
  cerr << "where <key> may be: " << endl
       << "r - extract density " << endl
       << "p - extract log10 pressure" << endl
       << "P - extract actual pressure" << endl 
       << "u - extract velocity vector uvel, vvel, wvel" << endl
       << "U - extract the magnitude of velocity " << endl  
       << "m - extract mach number" << endl
       << "t - extract temperature" << endl
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
inline bool compare_face(const std::vector<int> &face1, const std::vector<int> &face2) {
  if(face1.size() != face2.size())
    return false ;
  size_t fsz = face1.size() ;
  size_t ind ;
  for(ind=0;ind<fsz;++ind)
    if(face2[ind] == face1[0])
      break ;
  if(ind == face2.size())
    return false ;
  size_t l ;
  for( l=1;l<fsz;++l) {
    int li = (fsz+ind-l)%fsz ;
    if(face1[l] != face2[li])
      break ;
  }
  if(l == fsz)
    return true ;
  for( l=1;l<fsz;++l) {
    size_t li = (ind+l)%fsz ;
    if(face1[l] != face2[li])
      break ;
  }
  return (l == fsz) ;
}
}
int FieldViewExtract(int ac, char *av[]) {
//int main(int ac, char *av[]) {
  if(ac <= 1) {
    FVUsage(av[0]) ;
    exit(0) ;
  }
   
  static int hex_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
				NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  static int tet_walls[4] = {  NOT_A_WALL, NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  static int pyramid_walls[5] = {  NOT_A_WALL, NOT_A_WALL, NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  /*
    static int hex_walls[6] = { A_WALL, A_WALL, A_WALL,
    A_WALL, A_WALL, A_WALL };
    static int tet_walls[4] = {  A_WALL, A_WALL, A_WALL, A_WALL };
    static int pyramid_walls[5] = { A_WALL, A_WALL, A_WALL, A_WALL, A_WALL };
  */
  
  char *ncyc = "grid" ;
  int ibuf[10] ;
  char fv[80] ;
  char name_buf[500] ;
  char *problem_name = av[1] ;
  char filename[512] ;
  ac-- ;
  av++ ;
  if(ac >= 2) {
    ncyc = av[1] ; 
    av++ ;
    ac-- ; 
  }
  char fv_name[256] ;
  int extra_vars = 0 ;
  int extra_vec_vars = 0 ;
  int num_extra_vars = 0 ;
  sprintf(fv_name, "fv_%s.bin", ncyc) ;
  FILE *OFP = fopen(fv_name, "wb") ;
  if(OFP == NULL) {
    cerr << "can't open file " << fv_name << endl ;
    exit(-1) ;
  }
  // Write out the magic number and the fieldview version info 
  ibuf[0] = FV_MAGIC ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  sprintf(fv, "FIELDVIEW") ;
  ibuf[0] = 2 ;
  ibuf[1] = 5 ;
  fwrite(&fv, sizeof(char), 80, OFP) ;
  fwrite(ibuf, sizeof(int), 2,OFP) ;
  
  float TIME, FSMACH, ALPHA, RE ;
  TIME = 0.0 ;
  FSMACH = 0.0 ;
  ALPHA = 0.0 ;
  RE = 0.0 ;
  fwrite(&TIME, sizeof(float), 1, OFP) ;
  fwrite(&FSMACH, sizeof(float), 1, OFP) ;
  fwrite(&ALPHA, sizeof(float), 1, OFP) ;
  fwrite(&RE, sizeof(float), 1, OFP) ;
  
  ibuf[0] = 1 ; // Number of grids 
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  int num_vars = 0 ;
  int num_bvars = 0 ;
  char** var_buf ;
  char** bvar_buf ;
  var_buf = new char*[20] ;
  bvar_buf = new char*[8] ;
  for(int i = 0; i < 20; ++i)
    var_buf[i] = new char[80] ;
  for(int i = 0; i < 8; ++i)
    bvar_buf[i] = new char[80] ;
  int tmp = 0 ;
  int b_tmp = 0 ;
  int index[7] ;
  if(ac == 1)
    sprintf(var_buf[0], "m") ;
  while(ac > 1) {
    if(!strcmp(av[1],"v")) {
      sprintf(var_buf[tmp++], "0") ;
      sprintf(var_buf[tmp++], "1") ;
      sprintf(var_buf[tmp++], "2") ;
    } else if(!strcmp(av[1],"y")) {
      bvar_buf[b_tmp] = "yplus_w" ;
      index[b_tmp++] = 6 ;
    }
    else if(!strcmp(av[1],"W")) {
      bvar_buf[b_tmp++] = "wstress_mag" ;
    }
    else if(!strcmp(av[1],"q"))  {
      bvar_buf[b_tmp] = "heatf" ;
      index[b_tmp++] = 3 ;
    }
    else if(!strcmp(av[1],"w")) {
      index[b_tmp] = 0 ;
      sprintf(bvar_buf[b_tmp++], "0") ;
      index[b_tmp] = 1 ;
      sprintf(bvar_buf[b_tmp++], "1") ;
      index[b_tmp] = 2 ;
      sprintf(bvar_buf[b_tmp++], "2") ;
    }
    else if(!strcmp(av[1],"wp")) {
      index[b_tmp] = 4 ;
      sprintf(bvar_buf[b_tmp++], "wall_pressure") ;
    } else if(!strcmp(av[1],"wt")) {
      index[b_tmp] = 5 ;
      sprintf(bvar_buf[b_tmp++], "wall_temperature") ;
    } 
    else {
      var_buf[tmp++] = av[1] ;
    }
    av++ ; 
    ac-- ;
  }
  if(tmp > 0)
    num_vars = tmp ;
  if(b_tmp > 0)
    num_bvars = b_tmp ;
  //else
  //num_vars = 0 ;
  /*
  if(b_tmp > 0) {
    if((b_tmp == 2) || (b_tmp > 3) || (num_vars > 1)) {
      cout << "Only one boundary variable can be extracted at a time " << endl ;
      cout << "Make sure there are no other variables " << endl ;
      cout << " Use the following variables one at a time: " << endl ;
    cout << " y ----- yplus " << endl ;
    cout << " q ----- heat flux " << endl ;
    cout << " w ----- shear stress vector" << endl ;
    exit(-1) ;
  }
  else if(b_tmp == 1) {
    sprintf(var_buf[0], "m") ;
    num_vars = 1 ;
    num_bvars = 1 ;
  }
  else if(b_tmp == 3) {
    tmp = 0 ;
    sprintf(var_buf[tmp++], "0") ;
    sprintf(var_buf[tmp++], "1") ;
    sprintf(var_buf[tmp++], "2") ;
    num_bvars = 3 ;
    num_vars = 3 ;
  }
  }
  */
  for(int i = 0; i < num_vars; ++i)
    if((var_buf[i][0] == 'n') || (var_buf[i][0] == 'x'))
      num_extra_vars++ ;
  int *var_map = new int[num_extra_vars] ;
  cout << "Total number of variables = " << num_vars << endl ;
  cout << "Number of boundary variables = " << num_bvars << endl ;
  cout << " Number of variables to be read from other files = " << num_extra_vars << endl ;
  sprintf(filename, "output/qn_xdr.%s", ncyc) ;
  sprintf(name_buf,"%s.xdr",problem_name) ;
  FILE *IFP_GRD = fopen(name_buf, "r") ;
  if(IFP_GRD == NULL) {
    cerr << "open failed on " << name_buf << endl ;
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
  for(size_t i = 1; i < offset.size(); ++i) {
    count[*ii] = offset[i] - offset[i-1] ;
    ++ii ;
  } 
  face2node.allocate(count) ;
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) 
    for(int i = 0; i < count[*ei]; ++i) {
      xdr_int(&xdr_grd, &face2node[*ei][i]) ;
      face2node[*ei][i] += 1 ;
      
    }
  //cout << " face2node = " << face2node << endl ;
  cout << " Done with creating cl, cr, face2node " << endl ; 
  entitySet cri = Loci::MapRepP(cr.Rep())->image(faces) ;
  entitySet orig_boundaries =  cri & interval(Loci::UNIVERSE_MIN,-1) ;
  ibuf[0] = orig_boundaries.size() ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  cout <<  "Num face types = " << ibuf[0] << endl ;
  char tmp_buf[80] ;
  char qf_name[512] ;
  sprintf(qf_name, "output/qf.%s", ncyc) ; 
  ifstream ifile(qf_name, ios::in) ;
  dmultiMap qf_f2n ;
  std::vector<float> vec_float ;
  std::vector<entitySet> vec_set ;
  dMap face_map ;
  if((!ifile.fail()) && num_bvars)  {
    cout << " Reading in boundary values " << endl ;
    tmp = 0 ;
    do{
      std::vector<int> tmp_vec ;
      int dummy ;
      char c ; 
      double d ;
      for(int i = 0; i < 4; ++i) {
	ifile >> dummy ;
	dummy += 1 ;
	if(dummy != 0)
	  tmp_vec.push_back(dummy) ;
      }
      qf_f2n[tmp] = tmp_vec ;
      tmp++ ;
      for(int i = 0; i < 3; ++i) 
	ifile >> d ;
      for(int i = 3; i < 10; ++i) {
	ifile >> d ;
	vec_float.push_back(float(d)) ;
      }
      c = ifile.get();
      while(!ifile.eof() && isspace(c)) 
	c = ifile.get() ;
      
      if(!isspace(c)) {
	ifile.putback(c);
      }
    } while(!ifile.eof()) ;
    cout << "Finished reading in boundary values " << endl ;
    cout << "Number of values read from the qf file = " << tmp << endl ;
    HASH_MAP(int, entitySet) node2face ;
    entitySet qf_dom = qf_f2n.domain() ;
    ifile.close() ;
    
    qf_dom = face2node.domain() ;
    for(entitySet::const_iterator ei = qf_dom.begin(); ei != qf_dom.end(); ++ei) {
      int sz = face2node.end(*ei) - face2node.begin(*ei) ;
      for(int i = 0; i < sz; ++i)
	node2face[face2node[*ei][i]] += *ei ;
    }
    qf_dom =  qf_f2n.domain() ;
    entitySet tmp_face ;
    entitySet left_out ;
    for(entitySet::const_iterator ei = qf_dom.begin(); ei != qf_dom.end(); ++ei) {
      std::vector<int> tmp_vec = qf_f2n[*ei] ;
      std::sort(tmp_vec.begin(), tmp_vec.end()) ;
      int nd1 = tmp_vec[0] ;
      bool match = false ;
      //for(int j = 0; j < node2face[nd1].size(); ++j) {
      for(entitySet::const_iterator pi = node2face[nd1].begin(); pi != node2face[nd1].end(); ++pi) {	
	std::vector<int> fn ;
	for(const int* k = face2node.begin(*pi); k != face2node.end(*pi); ++k)
	  fn.push_back(*k) ;
	std::sort(fn.begin(), fn.end()) ;
	if(compare_face(tmp_vec,fn)) {
	  if(match) {
	    cerr << "found two matches!" << endl ;
	  } 
	  match = true ;
	  tmp_face += *pi ;
	  face_map[*pi] = *ei ;
	  left_out += *ei ;
	}
      } 
    }
    cout << "Finished matching the faces " << endl ;   
    cout << "Number of faces found by matching = " << tmp_face.size() << endl ;     
    left_out = qf_f2n.domain() - left_out ;
    cout << " left_out = " << left_out << endl ;
    cout << "Number of faces left out = " << left_out.size() << endl ;
    
    for(entitySet::const_iterator ei = left_out.begin(); ei != left_out.end(); ++ei) {
      cout << "qf_f2n[ " << *ei << " ] = " ;
      for(size_t i = 0; i < qf_f2n[*ei].size(); ++i)
	cout << qf_f2n[*ei][i] << "  " ;
      cout << endl ;
    } 
    
    entitySet rem ;
    entitySet finally = tmp_face  ;
    for(entitySet::const_iterator ei = orig_boundaries.begin(); ei != orig_boundaries.end(); ++ei) {
      entitySet tmp_image = cr.preimage(entitySet(interval(*ei, *ei))).first ;
      rem = tmp_image & tmp_face ;
      if(rem != EMPTY) {
	ibuf[0] = 1 ;
	vec_set.push_back(rem) ;
	finally -= rem ;
      }
      else
	ibuf[0] = 0 ;
      ibuf[1] = 0 ;
      sprintf(tmp_buf,"BC_%d",-*ei) ; 
      cout << ibuf[0] <<  "  " << ibuf[1] << "  " << tmp_buf << endl ;
      fwrite(ibuf, sizeof(int), 2, OFP) ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    }
    if(finally != EMPTY) {
      cerr << "Some boundary faces left unmatched " << endl ;
      cerr << "Left out faces = " << finally << endl ;
      cerr << " Number of faces left out = " << finally.size() << endl ;
    }
    
    fwrite(&num_vars, sizeof(int), 1, OFP) ;
    for(int i = 0; i < num_vars; ++i) {
      int type = *var_buf[i] ;
      if(type >= '0' && type <= '9') {
	memset(tmp_buf, '\0', 80) ;
	if(type == '0') {
	  sprintf(tmp_buf, "uvel;v") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	  //cout << "variable =   " << tmp_buf << endl ;
	}
	if(type == '1') {
	  memset(tmp_buf, '\0', 80) ;
	  sprintf(tmp_buf, "vvel") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	  //cout << "variable =  " << tmp_buf << endl ;
	}
	if(type == '2') {
	  memset(tmp_buf, '\0', 80) ;
	  sprintf(tmp_buf, "wvel") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ; 
	  //cout << "variable = " << tmp_buf << endl ;
	}
      } 
      else {
	fwrite(var_buf[i], sizeof(char), 80, OFP) ;
	//cout << "variable = " << var_buf[i] << endl ;
      }
    }
    fwrite(&num_bvars, sizeof(int), 1, OFP) ;
    for(int i = 0; i < num_bvars; ++i) {
      int type = *bvar_buf[i] ;
      if(type >= '0' && type <= '9') {
	memset(tmp_buf, '\0', 80) ;
	if(type == '0') {
	  sprintf(tmp_buf, "ws_x;wall_stress") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	  //cout << "variable =   " << tmp_buf << endl ;
	}
	if(type == '1') {
	  memset(tmp_buf, '\0', 80) ;
	  sprintf(tmp_buf, "ws_y") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	  //cout << "variable =  " << tmp_buf << endl ;
	}
	if(type == '2') {
	  memset(tmp_buf, '\0', 80) ;
	  sprintf(tmp_buf, "ws_z") ;
	  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	  //cout << "variable = " << tmp_buf << endl ;
	}
      }
      else 
	fwrite(bvar_buf[i], sizeof(char), 80, OFP) ;
      //cout << "variable = " << bvar_buf[i] << endl ;
    }
  }
  else { 
    cerr << "No Boundary values present " << endl ;
    FORALL(orig_boundaries, bc) {
      ibuf[0] = 0 ;
      ibuf[1] = 0 ;
      sprintf(tmp_buf,"BC_%d",-bc) ; 
      fwrite(ibuf, sizeof(int), 2, OFP) ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    } ENDFORALL ;
    fwrite(&num_vars, sizeof(int), 1, OFP) ;
    for(int i = 0; i < num_vars; ++i) {
      int type = *var_buf[i] ;
       if(type >= '0' && type <= '9') {
	 memset(tmp_buf, '\0', 80) ;
	 if(type == '0') {
	   sprintf(tmp_buf, "uvel;v") ;
	   fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	   //cout << "variable =   " << tmp_buf << endl ;
	 }
	 if(type == '1') {
	   memset(tmp_buf, '\0', 80) ;
	   sprintf(tmp_buf, "vvel") ;
	   fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	   //cout << "variable =  " << tmp_buf << endl ;
	 }
	 if(type == '2') {
	   memset(tmp_buf, '\0', 80) ;
	   sprintf(tmp_buf, "wvel") ;
	   fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	   //cout << "variable = " << tmp_buf << endl ;
	 }
       }
       else {
	 fwrite(var_buf[i], sizeof(char), 80, OFP) ;
	 //cout << "variable = " << var_buf[i] << endl ;
       }
    }
    ibuf[0] = 0 ; // number of boundary variables 
    fwrite(ibuf, sizeof(int), 1, OFP) ;
  }
  ibuf[0] = FV_NODES ; 
  ibuf[1] = npnts ;
  fwrite(ibuf, sizeof(int), 2, OFP) ;
  
  /* Coordinate information */
  for(int i = 0; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  
  for(int i = 1; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  
  for(int i = 2; i < 3*npnts; i+=3) 
    fwrite(&tmp_pos[i], sizeof(float), 1, OFP) ;
  delete [] tmp_pos ;
  
  tmp = 1 ;
  entitySet unclassified_faces, unclassified_cells ; 
  for(entitySet::const_iterator ei = orig_boundaries.begin(); ei != orig_boundaries.end(); ++ei) {
    entitySet tmp_image = cr.preimage(entitySet(interval(*ei, *ei))).first ;
    ibuf[0] = FV_FACES;
    ibuf[1] = tmp;	/* first boundary type */
    ibuf[2] = tmp_image.size() ;
    cout << " FV_FACES = " << ibuf[0] << "  boundary_type = " << ibuf[1] << "  number of faces = " << ibuf[2] << endl ;
    fwrite(ibuf, sizeof(int), 3, OFP) ;
    for(entitySet::const_iterator ej = tmp_image.begin(); ej != tmp_image.end(); ++ej) {
      if((count[*ej] != 4) && (count[*ej] != 3)) {
	cerr << "ERROR: Number of nodes per face is equal to  " << count[*ej] << "  : Currently not supported " << endl ;  
	unclassified_faces += *ej ; 
      }
      if(count[*ej] == 4) {
	fwrite(&face2node[*ej][0], sizeof(int), 4, OFP) ;
      } else {
	fwrite(&face2node[*ej][0], sizeof(int), count[*ej], OFP) ;
	ibuf[0] = 0 ;
	fwrite(ibuf, sizeof(int), 1, OFP) ;
      }
    }
    tmp++ ;
  }
  store<entitySet> cell2faces ;
  cell2faces.allocate(cells) ;
  for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ++ei) {
    cell2faces[cl[*ei]] += *ei ;
    if(cr[*ei] > 0)
      cell2faces[cr[*ei]] += *ei ;
  }
  int sizes[4] ;
  for(int i = 0; i < 4; ++i)
    sizes[i] = 0 ;
  char efile_name[512] ;
  sprintf(efile_name, "%s_elem_file.bin", problem_name) ; 
  struct stat stat_elem, stat_grid ;
  int success1 = stat(efile_name, &stat_elem) ;
  int success2 =stat(name_buf, &stat_grid) ;
  FILE *EFILE = NULL;
  if(success2 != 0) {
    cerr << "unable to stat " << name_buf<< endl ;
  }
  if(success1 == 0 && stat_elem.st_mtime > stat_grid.st_mtime) 
    EFILE = fopen(efile_name, "rb") ;
  bool regenerate = 0 ;
  fpos_t before_pos ;
  fgetpos(OFP, &before_pos) ;
  if(EFILE != NULL) {
    cout << " READING IN ELEMENT INFORMATION " << endl ;
    unsigned int header ;
    int tmp_count ;
    fread(ibuf, sizeof(int), 5, EFILE) ;
    int total_elems = ibuf[1] + ibuf[2] + ibuf[3] + ibuf[4] ;
    fwrite(ibuf, sizeof(int), 5, OFP) ; 
    if(ncells == total_elems) {
      for(int i = 0; i < total_elems; ++i) {
	if(feof(EFILE)) {
	  regenerate = 1 ;
	  cerr << "Premature end of element file " << endl ;
	  cerr << "Generating a new element file " << endl ;
	  i = total_elems-1 ;
	}
	fread(&header, sizeof(unsigned int), 1, EFILE) ;
	fread(&tmp_count, sizeof(int), 1, EFILE) ;
	fread(ibuf, sizeof(int), tmp_count, EFILE) ;
	fwrite(&header, sizeof(unsigned int), 1, OFP) ;
	fwrite(ibuf, sizeof(int), tmp_count, OFP) ;
      }
    }
    else {
      cerr << "There seems to be a mismatch in the number of elements" <<endl ;
      cerr << "Number of elements in the grid = " << ncells << endl ;
      cerr << " Number of elements in the element file = " << total_elems << endl ;
      cerr << "Generating a new element file " << endl ;
      regenerate = 1 ;
    }
  }
  if((regenerate) || (EFILE == NULL)) {
    const fpos_t after_pos = before_pos ;
    fsetpos(OFP, &after_pos) ;
    FILE *OFILE = fopen(efile_name, "wb") ;
    int num_triangles, num_quads ;
    store<int> elem_type ;
    elem_type.allocate(cells) ;
    for(entitySet::const_iterator ei = cells.begin(); ei != cells.end(); ++ei) {
      num_triangles = 0 ;
      num_quads = 0 ;
      for(entitySet::const_iterator ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
	if(count[*ej] == 3)
	  num_triangles++ ;
	else if(count[*ej] == 4)
	  num_quads++ ;
      }
      if(!num_quads) {
	sizes[0]++ ; //Increment the number of tets
	elem_type[*ei] = FV_TET_ELEM_ID ;
      }
      else if(!num_triangles) {
	sizes[1]++ ; // Increment the number of hexs
	elem_type[*ei] = FV_HEX_ELEM_ID ;
      }
      else if(num_quads > num_triangles) {
	sizes[2]++ ; // Increment the number of prisms
	elem_type[*ei] = FV_PRISM_ELEM_ID ;
      }
      else if(num_triangles > num_quads) {
	sizes[3]++ ; // Increment the number of pyramids
	elem_type[*ei] = FV_PYRA_ELEM_ID ;
      }
      else {
	unclassified_cells += *ei ;
	cerr << " Cell " << *ei << " does not belong to the basic type " << endl ;
      }
    }
    ibuf[0] = FV_ELEMENTS ;
    unsigned int elem_header ;
    ibuf[1] =  sizes[0] ;// tet_count ;
    ibuf[2] =  sizes[1] ;// hex_count ;
    ibuf[3] =  sizes[2] ;// prism_count ;
    ibuf[4] =  sizes[3] ;// pyramid_count ;
    fwrite(ibuf, sizeof(int), 5, OFP) ;
    fwrite(ibuf, sizeof(int), 5, OFILE) ; 
    cout << " Number of tets = " << sizes[0] << endl ;
    cout << " Number of hexs = " << sizes[1] << endl ;
    cout << " Number of prisms = " << sizes[2] << endl ;
    cout << " Number of pyramids = " << sizes[3] << endl ;
    cout << " Number of unclassified faces = " << unclassified_faces.size() << endl ;
    cout << " Number of unclassified cells = " << unclassified_cells.size() << endl ;
    for(entitySet::const_iterator ei = cells.begin(); ei != cells.end(); ++ei) {
      entitySet total_nodes, final_nodes ;
      int quad_id[6], triangle_id[2], common[2] ;
      int quad_count = 0 , triangle_count = 0, tmp_count = 0  ;
      entitySet match_pair_1, match_pair_2, nodes_1, nodes_2, tmp_nodes ;
      entitySet::const_iterator ej ;
      switch(elem_type[*ei]) {
      case FV_TET_ELEM_ID:
	for( ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
	  for(int i = 0; i < count[*ej]; ++i)
	    total_nodes += face2node[*ej][i] ;
	  if(total_nodes.size() == 4)
	    break ;
	}
	ej = cell2faces[*ei].begin() ;
	for(int i = 0; i < count[*ej]; ++i) {
	  ibuf[i] = face2node[*ej][i] ;
	  final_nodes += ibuf[i] ;
	}
	total_nodes -= final_nodes ;
	if(total_nodes.size() > 1)
	  cout << "Error : Tet element ordering screwed up " << endl ;
	ej = total_nodes.begin() ;
	ibuf[3] = *ej ;
	elem_header = fv_encode_elem_header(FV_TET_ELEM_ID, tet_walls);
	if (elem_header == 0)
	  {
	    cerr <<  "fv_encode_elem_header failed for tets. " << endl ;
	    exit(-1);
	  }
	fwrite(&elem_header, sizeof(unsigned int), 1, OFP) ;
	fwrite(&elem_header, sizeof(unsigned int), 1, OFILE) ;
	tmp_count = 4 ;
	fwrite(&tmp_count, sizeof(int), 1, OFILE) ;
	fwrite(ibuf, sizeof(int), 4, OFP) ;
	fwrite(ibuf, sizeof(int), 4, OFILE) ;
	break ;
      case FV_HEX_ELEM_ID:
	for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) 
	  quad_id[quad_count++] = *ej ;
	
	ibuf[0] = face2node[quad_id[0]][0] ;
	ibuf[1] = face2node[quad_id[0]][1] ;
	ibuf[2] = face2node[quad_id[0]][3] ;
	ibuf[3] = face2node[quad_id[0]][2] ;
	for(int i = 0; i < 4; ++i) {
	  tmp_count = 0 ;
	  for(int j = 1; j < quad_count; ++j) {
	    for(int k = 0; k < 4; ++k)
	      if(ibuf[i] == face2node[quad_id[j]][k]) 
		common[tmp_count++] = quad_id[j] ;
	  }
	  nodes_1 = EMPTY ;
	  nodes_2 = EMPTY ;
	  for(int k = 0; k < 4; ++k)
	    nodes_1 += face2node[common[0]][k] ;
	  for(int k = 0; k < 4; ++k)
	    nodes_2 += face2node[common[1]][k] ;
	  match_pair_2 = nodes_1 & nodes_2 ;
	  match_pair_2 -= interval(ibuf[i], ibuf[i]) ;
	  if(match_pair_2.size() != 1) {
	    cout << " Element = " << *ei << endl ;
	    cout << "Error: Hex elem ordering screwed up " <<  endl ;
	  }
	  ej = match_pair_2.begin() ;
	  ibuf[i+4] = *ej ;
	}
	elem_header = fv_encode_elem_header(FV_HEX_ELEM_ID, hex_walls) ;
	if (elem_header == 0)
	  {
	    cerr << "fv_encode_elem_header failed for  hex. " << endl ;
	    exit(-1);
	  }
	fwrite(&elem_header, sizeof(unsigned int), 1, OFP);
	fwrite(&elem_header, sizeof(unsigned int), 1, OFILE) ;
	fwrite(ibuf, sizeof(int), 8, OFP) ;
	tmp_count = 8 ;
	fwrite(&tmp_count, sizeof(int), 1, OFILE) ;
	fwrite(ibuf, sizeof(int), 8, OFILE) ;
	break ;
      case FV_PRISM_ELEM_ID:
	for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
	  if(count[*ej] == 4) {
	    quad_id[quad_count] = *ej ;  
	    quad_count++ ;
	  } else if(count[*ej] == 3) {
	    triangle_id[triangle_count] = *ej ;
	    triangle_count++ ;
	  } 
	} 
	
	//for(int i = 0; i < 4; ++i)
	//cout <<  "face2node[ " << i << " ] = " << face2node[quad_id[0]][i] << endl ;
	for(int j=0; j < 3; ++j)
	  nodes_1 += face2node[triangle_id[0]][j] ;
	for(int j=0; j < 3; ++j)
	  nodes_2 += face2node[triangle_id[1]][j] ;
	ibuf[0] = face2node[quad_id[0]][0] ;
	ibuf[1] = face2node[quad_id[0]][1] ;
	ibuf[2] = face2node[quad_id[0]][2] ;
	ibuf[3] = face2node[quad_id[0]][3] ;
        tmp = 1 ;
	while(!(nodes_1.inSet(ibuf[0]) && nodes_2.inSet(ibuf[1]))) {
	  for(int i = 0; i < 4; ++i)
	    ibuf[i] = face2node[quad_id[0]][((i+tmp)%4)] ;
	  tmp++ ;
	}
      	
	match_pair_1 += ibuf[0] ;
	match_pair_1 += ibuf[3] ;
	
	match_pair_2 += ibuf[1] ;
	match_pair_2 += ibuf[2] ;
	
	if((match_pair_2 & nodes_1) != EMPTY) {
	  tmp_nodes = nodes_1 - match_pair_2 ;
	  if(tmp_nodes.size() != 1)
	    cerr << "Error : Prism ordering screwed up: element =  " << *ei  << endl ;
	  ej = tmp_nodes.begin() ;
	  ibuf[4] = *ej ;
	  
	} else if((match_pair_2 & nodes_2) != EMPTY) {
	  tmp_nodes = nodes_2 - match_pair_2 ;
	  if(tmp_nodes.size() != 1)
	    cerr << "Error : Prism ordering screwed up: element =  " << *ei  << endl ;
	  ej = tmp_nodes.begin() ;
	  ibuf[4] = *ej ;
	}
	if((match_pair_1 & nodes_1) != EMPTY) {
	  tmp_nodes = nodes_1 - match_pair_1 ;
	  if(tmp_nodes.size() != 1)
	    cerr << "Error : Prism ordering screwed up: element =  " << *ei  << endl ;
	  ej = tmp_nodes.begin() ;
	  ibuf[5] = *ej ;
	  
	} else if((match_pair_1 & nodes_2) != EMPTY) {
	  tmp_nodes = nodes_2 - match_pair_1 ;
	  if(tmp_nodes.size() != 1)
	    cerr << "Error : Prism ordering screwed up: element =  " << *ei  << endl ;
	  ej = tmp_nodes.begin() ;
	  ibuf[5] = *ej ;
	}
	elem_header = fv_encode_elem_header(FV_PRISM_ELEM_ID, pyramid_walls);
	if(elem_header == 0)
	  {
	    cerr << "fv_encode_elem_header failed for  prism. " << endl ;
	    exit(-1);
	  }
	fwrite (&elem_header, sizeof(unsigned int), 1, OFP);
	fwrite(ibuf, sizeof(int), 6, OFP) ;
	fwrite(&elem_header, sizeof(unsigned int), 1, OFILE) ;
	tmp_count = 6 ;
	fwrite(&tmp_count, sizeof(int), 1, OFILE) ;
	fwrite(ibuf, sizeof(int), 6, OFILE) ;
	//cout << "  Nodes are  " ;
	//for(int i = 0; i < 6; ++i)
	//cout <<  ibuf[i] << "  " ;
	//cout << endl ;
	break ; 
      case FV_PYRA_ELEM_ID:
	//cout << "pyramid elem " << *ei << endl ;
	for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
	  if(count[*ej] == 4)
	    quad_id[0] = *ej ;
	  for(int i = 0; i < count[*ej]; ++i)
	    total_nodes += face2node[*ej][i] ;
	}
	//cout << " total_nodes=  " << total_nodes << endl ;
	
	for(int i = 0; i < count[quad_id[0]]; ++i) {
	  ibuf[i] = face2node[quad_id[0]][i] ;
	  final_nodes += ibuf[i] ;
	}  
	//cout << " final_nodes = " << final_nodes << endl ;
	total_nodes -= final_nodes ;
	if(total_nodes.size() != 1)
	  cerr << "Pyramid elements ordering screwed up " << endl ;
	ej = total_nodes.begin() ;
	ibuf[4] = *ej ;
	elem_header = fv_encode_elem_header(FV_PYRA_ELEM_ID, pyramid_walls);
	if (elem_header == 0)
	  {
	    cerr << "fv_encode_elem_header failed for  pyramid. " << endl ;
	    exit(-1);
	  }
	fwrite (&elem_header, sizeof(unsigned int), 1, OFP) ;
	fwrite(ibuf, sizeof(int), 5, OFP) ;
	fwrite(&elem_header, sizeof(unsigned int), 1, OFILE) ;
	tmp_count = 5 ;
	fwrite(&tmp_count, sizeof(int), 1, OFILE) ;
	fwrite(ibuf, sizeof(int), 5, OFILE) ;
	break ; 
      default:
	cerr <<  "ERROR:  Unknown element type " << endl ;
	return 0;	
      }
    }   
    fclose(OFILE) ;
  }
  cout << " Done with grouping the element types " << endl ;
  ibuf[0] = FV_VARIABLES ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  for(int n = 0; n < num_vars; ++n) {
    int type = *var_buf[n] ;
    int count = 0 ;
    if(type >= '0' && type <= '9') {
      if(type == '0') {
	memset(tmp_buf, '\0', 80) ;
	sprintf(tmp_buf, "uvel") ;
	cout << "Writing out variable =   " << tmp_buf << endl ;
      }
      if(type == '1') {
	memset(tmp_buf, '\0', 80) ;
	sprintf(tmp_buf, "vvel") ;
	cout << "Writing out variable =  " << tmp_buf << endl ;
      }
      if(type == '2') {
	memset(tmp_buf, '\0', 80) ;
	sprintf(tmp_buf, "wvel") ;
	cout << "Writing out variable = " << tmp_buf << endl ;
      }
    }
    else {
      if(var_buf[n][0] == 'f') {
	std::string s = string(var_buf[n] + 1) ;
	cout << "Writing out variable " << s << endl ;
      } else if(var_buf[n][0] == 'n') {
	std::string s = string(var_buf[n] + 1) ;
	cout << "Writing out variable " << s << endl ;
	extra_vars = 1 ;
	var_map[count++] = n ;
      } else if(var_buf[n][0] == 'x') {
	std::string s = string(var_buf[n] + 1) ;
	cout << "Writing out variable " << s << endl ;
	extra_vec_vars = 1 ;
	var_map[count++] = n ;
      }
      else
	cout << "Writing out variable " << var_buf[n] << endl ;
    }
    if((!extra_vars) && (!extra_vec_vars)) {
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
      for(int i=0;i<ns;++i) {
	xdr_int(&xdr_qn, &tmp_size) ;
	for(int j = 0; j < tmp_size; ++j)
	  xdr_char(&xdr_qn, &tmp_char[j]) ;
	for(int k = 0; k < tmp_size; ++k)
	  species[i] += tmp_char[k] ;
      }
      int nnodes  =0;
      xdr_int(&xdr_qn, &nnodes) ;
      if(nnodes != npnts) {
	cerr << " Error : number of nodes from the output file and the grid file differ " << endl ;
	cerr <<  "This should not happen: See whether the run is from a plot-3d grid or not " << endl ;
	cerr << "Currently Plot3d grid runs cannot be visualized with the xdr format generated from cobalt grid  as the number of nodes don't match " << endl ;
	exit(-1) ;
      }
      int squery=-1 ;
      if(var_buf[n][0] == 'f') {
	string s = string(var_buf[n]+1) ;
	for(int i=0;i<ns;++i)
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
      for(int i=0;i<nnodes;++i) {
	for(int ii=0;ii<nv;++ii) 
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
	case 't':
	  val = T ;
	  break ;
	case '0':
	  val = u ;
	  break ;
	case '1': 
	  val = v ;
	  break ;
	case '2':
	  val = w ;
	  break ;
        case 'a':
          val = a ;
          break ;
        case 'g':
          val = Gamma ;
          break ;
        case 'R':
          val = Rt ;
          break ;
	case 'f':
	  t = squery+7 ;
	 val =  qb[t] ;
	 break ;	
	default:
	  cerr << "wrong type"<< endl ;
	  exit(0) ;
	}
	fwrite(&val, sizeof(float), 1, OFP) ;
	//cout << "  variable = " << var_buf[n] << " node number = " << i <<"  value " << val << endl ;
      }
      delete [] qb ;
      xdr_destroy(&xdr_qn) ;
      fclose(IFP_QN) ;
    }
    else if(extra_vars) {
      char tmp_name[512] ;
      char *tmp_char = var_buf[var_map[count-1]] + 1 ;
      sprintf(tmp_name, "output/%s_hdf5.%s",tmp_char, ncyc) ;
      cout << "opening file " << tmp_name << endl ;
      store<double> read_scalar ;
      hid_t  file_id, group_id;
      file_id = H5Fopen(tmp_name, H5F_ACC_RDONLY, H5P_DEFAULT);
      group_id = H5Gopen(file_id, tmp_char) ;
      //      entitySet dom = Loci::interval(0, npnts-1) ;
      entitySet tmpdom = EMPTY ;
      Loci::read_container(group_id, read_scalar.Rep(), tmpdom) ;
      if(tmpdom.size() != npnts) {
        cerr << " domains don't match in size for variable " << tmp_char
             << endl ;
      }
      float val ;
      for(entitySet::const_iterator ei=tmpdom.begin();
          ei!=tmpdom.end();++ei) {
	val = float(read_scalar[*ei]) ;
	fwrite(&val, sizeof(float), 1, OFP) ;
      }
      extra_vars = 0 ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
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
      //      entitySet dom = Loci::interval(0, npnts-1) ;
      entitySet tmpdom = EMPTY ;
      Loci::read_container(group_id, read_vector.Rep(), tmpdom) ;
      float val[3] ;
      if(tmpdom.size() != npnts) {
        cerr << " domains don't match in size for variable " << tmp_name
            << endl ;
      }
      for(entitySet::const_iterator ei=tmpdom.begin();
          ei!=tmpdom.end();++ei) {
        val[0] = read_vector[*ei].x ;
        val[1] = read_vector[*ei].y ;
        val[2] = read_vector[*ei].z ;
	fwrite(val, sizeof(float), 3, OFP) ;
      }
      extra_vec_vars = 0 ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
    } 
  }
  ibuf[0] = FV_BNDRY_VARS ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  for(int i = 0; i < num_bvars; ++i) {
    if(!strcmp(bvar_buf[i], "wstress_mag")) {
      cout << "Writing out variable " << bvar_buf[i] << endl ;
      for(size_t j = 0; j < vec_set.size(); ++j)
	for(entitySet::const_iterator ei = vec_set[j].begin(); ei != vec_set[j].end(); ++ei) {
	  float val = 0.0 ;
	  double mag = 0.0 ;
	  for(int p = 0; p < 3; ++p)
	    mag += vec_float[7*face_map[*ei] + p] * vec_float[7*face_map[*ei] + p] ;
	  val = sqrt(mag) ;
	  fwrite(&val, sizeof(float), 1, OFP) ;
	}
    }
    else {
      //ibuf[0] = flag.domain().size() ;
      int type = *bvar_buf[i] ;
      if(type == '0')
	cout << "Writing out x component of wall stress ws_x " << endl ;
      else if(type == '1')
	cout << "Writing out y component of wall stress ws_y " << endl ;
      else if(type == '2')
	cout << "Writing out z component of wall stress ws_z " << endl ;
      else 
	cout << "Writing out variable " << bvar_buf[i] << endl ;
      for(size_t j = 0; j < vec_set.size(); ++j)
	for(entitySet::const_iterator ei = vec_set[j].begin(); ei != vec_set[j].end(); ++ei)
	  fwrite(&vec_float[7*face_map[*ei] + index[i]], sizeof(float), 1, OFP) ;
    }
  }
  cout << "Done "<< endl ;
  xdr_destroy(&xdr_grd) ;
  fclose(IFP_GRD) ;
  fclose(OFP) ; 
  return(0) ;
}
  
 

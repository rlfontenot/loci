
#include <Loci.h>
#include <stdlib.h>
#include <math.h>
#include <string>
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

#include "tec_extract.h"
#include <vector>
using std::vector ;

void TECUsage(char *s) {
  cerr << "Usage: " << endl 
       << string(s) << " -tec <problem_name> <time_step> <key> " << endl << endl ;
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
       << string(s) << "-tec ssme 200 fOH" << endl ;
}

  int TECINI(char *problem_name,char *variables_name,
             char *tec_name,char *dot,int *Debug,int *VIsDouble) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  cout << "Initializing file " << filename << endl ;
 
  ofstream ofile(filename.c_str(),ios::out) ;
  ofile << "variables = " << variables_name << endl ;
  ofile.close() ;

  return(0) ;
}

  int TECZNE(char *zone_name, int *npnts, int *ncells,char *tec_name,
             char *block_type, char *elem_type) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
 
  ofile << "Zone, T = " << zone_name << ", n = " << *npnts << ", e = " <<
            *ncells << ", F = " << block_type << ", ET = " << elem_type << endl ;
  ofile.close() ;

  return(0) ;
}

  int TECDAT(int *npnts, float *vals, char *tec_name) {

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

  int TECNOD(int *nm, int *ncells, char *tec_name) {

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

int TecplotExtract(int ac, char *av[]) {

  if(ac <= 1) {
    TECUsage(av[0]) ;
    exit(0) ;
  }
   
  char *ncyc = "grid" ;
  int ibuf[10] ;
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
  
  int num_vars = 0 ;
  int extra_vars = 0 ;
  int extra_vec_vars = 0 ;
  int num_extra_vars = 0 ;
  char** var_buf ;
  var_buf = new char*[20] ;
  for(int i = 0; i < 20; ++i)
    var_buf[i] = new char[80] ;
  int tmp = 0 ;
  if(ac == 1)
    sprintf(var_buf[0], "m") ;
  while(ac > 1) {
    if(!strcmp(av[1],"v")) {
      sprintf(var_buf[tmp++], "0") ;
      sprintf(var_buf[tmp++], "1") ;
      sprintf(var_buf[tmp++], "2") ;
    } else {
      var_buf[tmp++] = av[1] ;
    }
    av++ ; 
    ac-- ;
  }
  if(tmp > 0)
    num_vars = tmp ;

for(int i = 0; i < num_vars; ++i)
    if((var_buf[i][0] == 'n') || (var_buf[i][0] == 'x'))
      num_extra_vars++ ;
  int *var_map = new int[num_extra_vars] ;
  cout << "Total number of variables = " << num_vars << endl ;
  cout << " Number of variables to be read from other files = " << num_extra_vars << endl ;
  sprintf(filename, "output/qn_xdr.%s", ncyc) ;
  sprintf(name_buf,"%s.xdr",problem_name) ;
  FILE *IFP_GRD = fopen(name_buf, "r") ;
cout << "IFP_GRD: " << IFP_GRD << endl ;
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
  float* tmp_pos_x,*tmp_pos_y,*tmp_pos_z ;
  double tmp_double ;
  tmp_pos_x = new float[npnts] ; 
  tmp_pos_y = new float[npnts] ;
  tmp_pos_z = new float[npnts] ;
  for(int j = 0; j < npnts; ++j) {
    xdr_double(&xdr_grd, &tmp_double) ;
    tmp_pos_x[j] = float(tmp_double) ;
    xdr_double(&xdr_grd, &tmp_double) ;
    tmp_pos_y[j] = float(tmp_double) ;
    xdr_double(&xdr_grd, &tmp_double) ;
    tmp_pos_z[j] = float(tmp_double) ;
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
  cout << "before read offset " << endl ;
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
  cout <<  "Num face types = " << ibuf[0] << endl ;
  char tmp_buf[80] ;
  char qf_name[512] ;
  sprintf(qf_name, "output/qf.%s", ncyc) ; 
  ifstream ifile(qf_name, ios::in) ;
  dmultiMap qf_f2n ;
  std::vector<float> vec_float ;
  std::vector<entitySet> vec_set ;
  dMap face_map ;

  cerr << "No Boundary values present " << endl ;
  FORALL(orig_boundaries, bc) {
    sprintf(tmp_buf,"BC_%d",-bc) ; 
  } ENDFORALL ;

  char variables_name[1000] ; variables_name[0]='\0' ;
  strcat(variables_name,"\"x\", \"y\", \"z\"") ;
  for(int i = 0; i < num_vars; ++i) {
    int type = *var_buf[i] ;
    if(type >= '0' && type <= '9') {
      if(type == '0') {
        strcat(variables_name,", \"uvel\"") ;
      }
      if(type == '1') {
        strcat(variables_name,", \"vvel\"") ;
      }
      if(type == '2') {
        strcat(variables_name,", \"wvel\"") ;
      }
    }
    else if(var_buf[i][0]=='x') {
      strcat(variables_name,", \"") ; strcat(variables_name,&var_buf[i][1]) ;
      strcat(variables_name,"X") ; strcat(variables_name,"\"") ;
      strcat(variables_name,", \"") ; strcat(variables_name,&var_buf[i][1]) ;
      strcat(variables_name,"Y") ; strcat(variables_name,"\"") ;
      strcat(variables_name,", \"") ; strcat(variables_name,&var_buf[i][1]) ;
      strcat(variables_name,"Z") ; strcat(variables_name,"\"") ;
    }
    else if(var_buf[i][0]=='n' || var_buf[i][0]=='f') {
      strcat(variables_name,", \"") ;
      strcat(variables_name,var_buf[i]+1) ;
      strcat(variables_name,"\"") ;
    }
    else {
      strcat(variables_name,", \"") ;
      strcat(variables_name,var_buf[i]) ;
      strcat(variables_name,"\"") ;
    }
  }

  int Debug = 1 ;
  int VIsDouble = 0 ;
  char tec_name[256] ;
  sprintf(tec_name, "tec_%s.dat", ncyc) ;
//Initialize Tecplot ascii output
  int err = TECINI(problem_name,variables_name,
                   tec_name,".",&Debug,&VIsDouble) ;

  err = TECZNE("SINGLE_ZONE", &npnts, &ncells,tec_name,
                   "FEBLOCK","BRICK") ;
  
  /* Coordinate information */
  err = TECDAT(&npnts,&tmp_pos_x[0],tec_name) ;
  err = TECDAT(&npnts,&tmp_pos_y[0],tec_name) ;
  err = TECDAT(&npnts,&tmp_pos_z[0],tec_name) ;

  delete [] tmp_pos_x ;
  delete [] tmp_pos_y ;
  delete [] tmp_pos_z ;
  
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

//get connectivity list NM
  int *NM = new int[ncells*8] ; //for Brick element only

  {
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
	elem_type[*ei] = TEC_TET_ELEM_ID ;
      }
      else if(!num_triangles) {
	sizes[1]++ ; // Increment the number of hexs
	elem_type[*ei] = TEC_HEX_ELEM_ID ;
      }
      else if(num_triangles < num_quads) {
        sizes[2]++ ; // Increment the number of pyramids
        elem_type[*ei] = TEC_PRISM_ELEM_ID ;
      }
      else if(num_triangles > num_quads) {
        sizes[3]++ ; // Increment the number of pyramids
        elem_type[*ei] = TEC_PYRA_ELEM_ID ;
      }
      else {
        cerr << "triangles= " << num_triangles << " quads= " << num_quads << endl ;
	cerr << " Cell " << *ei << " does not belong to the basic type " << endl ;
      }
    }
    ibuf[1] =  sizes[0] ;// tet_count ;
    ibuf[2] =  sizes[1] ;// hex_count ;
    ibuf[3] =  sizes[2] ;// prism_count ;
    ibuf[4] =  sizes[3] ;// pyramid_count ;
    cout << " Number of tets = " << sizes[0] << endl ;
    cout << " Number of hexs = " << sizes[1] << endl ;
    cout << " Number of prisms = " << sizes[2] << endl ;
    cout << " Number of pyramids = " << sizes[3] << endl ;

    int cn = 0 ;
    for(entitySet::const_iterator ei = cells.begin(); ei != cells.end(); ++ei) {
      entitySet total_nodes, final_nodes ;
      int quad_id[6], triangle_id[3],common[2] ;
      int quad_count = 0 , triangle_count = 0, tmp_count = 0  ;
      entitySet match_pair_1, match_pair_2, nodes_1, nodes_2, tmp_nodes ;
      entitySet::const_iterator ej ;
      switch(elem_type[*ei]) {
      case TEC_HEX_ELEM_ID:
        for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej)
          quad_id[quad_count++] = *ej ;
	
//changed from FV order
	ibuf[0] = face2node[quad_id[0]][0] ;
	ibuf[1] = face2node[quad_id[0]][3] ;
	ibuf[2] = face2node[quad_id[0]][2] ;
	ibuf[3] = face2node[quad_id[0]][1] ;
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
        
        for(int i=0;i<8;++i)
          NM[8*cn+i] = ibuf[i ] ;
        cn++ ;
	break ;
      case TEC_TET_ELEM_ID:
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
        ibuf[4] = *ej ;

        ibuf[3] = ibuf[2] ;
        ibuf[5] = ibuf[4] ;
        ibuf[6] = ibuf[4] ;
        ibuf[7] = ibuf[4] ;

        for(int i=0;i<8;++i)
          NM[8*cn+i] = ibuf[i ] ;
        cn++ ;
        break ;
      case TEC_PRISM_ELEM_ID:
#ifdef OLDSTUFF
        for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
          if(count[*ej] == 4) {
            quad_id[quad_count] = *ej ;
            quad_count++ ;
          } else if(count[*ej] == 3) {
            triangle_id[triangle_count] = *ej ;
            triangle_count++ ;
          }
        }

        for(int j=0; j < 3; ++j)
          nodes_1 += face2node[triangle_id[0]][j] ;
        for(int j=0; j < 3; ++j)
          nodes_2 += face2node[triangle_id[1]][j] ;
        ibuf[0] = face2node[quad_id[0]][0] ;
        ibuf[1] = face2node[quad_id[0]][3] ;
        ibuf[2] = face2node[quad_id[0]][2] ;
        ibuf[3] = face2node[quad_id[0]][1] ;
        if(nodes_1.inSet(ibuf[0])) {
          nodes_1 -= ibuf[0] ;
          nodes_1 -= ibuf[1] ;
          if(nodes_1.size()>1)
            cerr << "something is wrong with the Prism cells\n" ;
          ibuf[4] = *(nodes_1.begin()) ; 
          ibuf[5] = ibuf[4] ;
          nodes_2 -= ibuf[2] ;
          nodes_2 -= ibuf[3] ;
          if(nodes_2.size()>1)
            cerr << "something is wrong with the Prism cells\n" ;
          ibuf[6] = *(nodes_2.begin()) ;
          ibuf[7] = ibuf[6] ;
        } else {
          nodes_2 -= ibuf[0] ;
          nodes_2 -= ibuf[1] ;
          if(nodes_2.size()>1)
            cerr << "something is wrong with the Prism cells\n" ;
          ibuf[4] = *(nodes_2.begin()) ;
          ibuf[5] = ibuf[4] ; 
          nodes_1 -= ibuf[2] ;
          nodes_1 -= ibuf[3] ;
          if(nodes_1.size()>1)
            cerr << "something is wrong with the Prism cells\n" ;
          ibuf[6] = *(nodes_1.begin()) ;
          ibuf[7] = ibuf[6] ;
        }
#endif
	for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
	  if(count[*ej] == 4) {
	    quad_id[quad_count] = *ej ;  
	    quad_count++ ;
	  } else if(count[*ej] == 3) {
	    triangle_id[triangle_count] = *ej ;
	    triangle_count++ ;
	  } 
	} 
	
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
        ibuf[7] = ibuf[5] ;
        ibuf[6] = ibuf[4] ;
        ibuf[5] = ibuf[6] ;
        ibuf[4] = ibuf[7] ;
        
        for(int i=0;i<8;++i)
          NM[8*cn+i] = ibuf[i ] ;
        cn++ ;
        break ;
      case TEC_PYRA_ELEM_ID:
        for(ej = cell2faces[*ei].begin(); ej != cell2faces[*ei].end(); ++ej) {
          if(count[*ej] == 4)
            quad_id[0] = *ej ;
          for(int i = 0; i < count[*ej]; ++i)
            total_nodes += face2node[*ej][i] ;
        }
        ibuf[0] = face2node[quad_id[0]][0] ;
        ibuf[1] = face2node[quad_id[0]][3] ;
        ibuf[2] = face2node[quad_id[0]][2] ;
        ibuf[3] = face2node[quad_id[0]][1] ;
        for(int i = 0; i < count[quad_id[0]]; ++i) 
          final_nodes += ibuf[i] ;
        total_nodes -= final_nodes ;
        if(total_nodes.size() != 1)
          cerr << "Pyramid elements ordering screwed up " << endl ;
        ej = total_nodes.begin() ;

        ibuf[4] = *ej ;
        ibuf[5] = ibuf[4] ;
        ibuf[6] = ibuf[4] ;
        ibuf[7] = ibuf[4] ;

        for(int i=0;i<8;++i)
          NM[8*cn+i] = ibuf[i ] ;
        cn++ ;
        break ;
      default:
	cerr <<  "ERROR:  only brick element type is supported" << endl ;
	return 0;	
      }
    }   
  }
  cout << " Done with grouping the element types " << endl ;

  float *values = new float[npnts] ;

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
      } else
	cout << "Writing out variable " << var_buf[n] << endl ;
    }

    if ((!extra_vars) && (!extra_vec_vars)) {
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
        cerr << "Currently Plot3d grid runs cannot be visualized with the xdr fo rmat generated from cobalt grid  as the number of nodes don't match " << endl ;
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
        case 'f':
          t = squery+7 ;
         val =  qb[t] ;
         break ;
        default:
          cerr << "wrong type"<< endl ;
          exit(0) ;
        }
        values[i] = val ;
      }
      delete [] qb ;
      xdr_destroy(&xdr_qn) ;
      fclose(IFP_QN) ;
      err = TECDAT(&nnodes,&values[0],tec_name) ;
    }
    else if(extra_vars) {
      char tmp_name[512] ;
      char *tmp_char = var_buf[var_map[count-1]] + 1 ;
      sprintf(tmp_name, "output/%s_hdf5.%s",tmp_char, ncyc) ;
      cout << "opening file " << tmp_name << endl ;
      store<float> read_scalar ;
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
      ostringstream oss ;
      oss << tec_name ;
      string filename = oss.str() ;
      ofstream ofile(filename.c_str(),ios::app) ;
      //      ofile.precision(16) ;

      float val ;
      for(entitySet::const_iterator ei=tmpdom.begin();
          ei!=tmpdom.end();++ei) {
          val = read_scalar[*ei] ;
          ofile << val << endl ;
//        fwrite(&val, sizeof(float), 1, OFP) ;
      }
      extra_vars = 0 ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
      ofile.close() ;
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
      if(tmpdom.size() != npnts) {
        cerr << " domains don't match in size for variable " << tmp_name
            << endl ;
      }
      ostringstream oss ;
      oss << tec_name ;
      string filename = oss.str() ;
      ofstream ofile(filename.c_str(),ios::app) ;
      float *xVal=new float[npnts],*yVal=new float[npnts],*zVal=new
        float[npnts] ;
      int count=0 ;
      for(entitySet::const_iterator ei=tmpdom.begin();
          ei!=tmpdom.end();++ei,++count) {
        xVal[count] = read_vector[*ei].x ;
        yVal[count] = read_vector[*ei].y ;
        zVal[count] = read_vector[*ei].z ;
//      fwrite(val, sizeof(float), 3, OFP) ;
      }
      for(int i=0;i<npnts;++i) ofile << xVal[i] << endl ;
      for(int i=0;i<npnts;++i) ofile << yVal[i] << endl ;
      for(int i=0;i<npnts;++i) ofile << zVal[i] << endl ;
      delete [] xVal ; delete [] yVal ; delete [] zVal ;
      extra_vec_vars = 0 ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
      ofile.close() ;
    }
  }

//write out nodal connectivity
  err = TECNOD(NM, &ncells, tec_name) ;
  cout << "Done "<< endl ;

  xdr_destroy(&xdr_grd) ;
  fclose(IFP_GRD) ;

  err = TECEND() ;

  return(0) ;
} 




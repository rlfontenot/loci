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

#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>


/* Numeric tags (codes) for FIELDVIEW binary file format. */
 
#define FV_MAGIC        0x00010203      /* decimal 66051 */
 
/* Content of the file (grid only, results only or combined). */
#define FV_GRIDS_FILE           1
#define FV_RESULTS_FILE         2
#define FV_COMBINED_FILE        3
 
#define FV_NODES                1001
#define FV_FACES                1002
#define FV_ELEMENTS             1003
#define FV_VARIABLES            1004
#define FV_BNDRY_VARS           1006
#define FV_ARB_POLY_FACES       1007
#define FV_ARB_POLY_ELEMENTS    1008
#define FV_ARB_POLY_BNDRY_VARS  1009
 
#define FV_TET_ELEM_ID          1
#define FV_HEX_ELEM_ID          2
#define FV_PRISM_ELEM_ID        3
#define FV_PYRA_ELEM_ID         4
#define FV_ARB_POLY_ELEM_ID     5
 
/* Values for "wall_info" array (see comments in fv_encode_elem_header). */
#ifdef __STDC__
#define A_WALL         (07u)
#define NOT_A_WALL     (0u)
#else
#define A_WALL         (07)
#define NOT_A_WALL     (0)
#endif

/* Don't change these - used by fv_encode_elem_header ! */
#define MAX_NUM_ELEM_FACES     6
#define BITS_PER_WALL  3
#define ELEM_TYPE_BIT_SHIFT    (MAX_NUM_ELEM_FACES*BITS_PER_WALL)
 
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
 
unsigned int fv_encode_elem_header (int elem_type, int wall_info[])
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

string convert_fv_compatible(string var) {
  if(var == "i")
    return string("i_") ;
  if(var == "j")
    return string("j_") ;
  if(var == "k")
    return string("k_") ;
  if(var == "x")
    return string("x_") ;
  if(var == "y")
    return string("y_") ;
  if(var == "z")
    return string("z_") ;
  if(var == "R")
    return string("R_") ;
  return var ;
}

void fv_topo_handler::open(string casename, string iteration ,int inpnts,
                           int intets, int inprsm, int inpyrm,
                           int inhexs, int ingen,
                           const vector<string> &bc_names,
                           const vector<string> &variables,
                           const vector<int> &variable_types) {
  npnts = inpnts ;
  ntets = intets ;
  nprsm = inprsm ;
  npyrm = inpyrm ;
  nhexs = inhexs ;
  ngen = ingen ;
  part_id = 1 ;
  filename = casename + "_fv_"+iteration + ".bin";
  OFP = fopen(filename.c_str(), "wb") ;
  if(OFP == NULL) {
    cerr << "can't open file " << filename << endl ;
    exit(-1) ;
  }
  int ibuf[10] ;
  char fv[80] ;
  // Write out the magic number and the fieldview version info 
  ibuf[0] = FV_MAGIC ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;
  memset(fv,'\0',80) ;
  sprintf(fv, "FIELDVIEW") ;
  fwrite(&fv, sizeof(char), 80, OFP) ;
  ibuf[0] = 3 ;
  ibuf[1] = 0 ;
  fwrite(ibuf, sizeof(int), 2,OFP) ;

  ibuf[0] = FV_COMBINED_FILE ;
  fwrite(ibuf,sizeof(int),1,OFP) ;

  ibuf[0] = 0 ;
  fwrite(ibuf,sizeof(int),1,OFP) ;

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

  vector<string> nlist ;
  vector<string> blist ;
  for(size_t i=0;i<variables.size();++i) {
    int vt = variable_types[i] ;
    const string var(convert_fv_compatible(variables[i])) ;
    if(vt == NODAL_SCALAR || vt == NODAL_DERIVED || vt == NODAL_MASSFRACTION)
      nlist.push_back(var) ;
    if(vt == NODAL_VECTOR) {
      string v = var ;
      nlist.push_back(v + "x ; " + v) ;
      nlist.push_back(v + "y") ;
      nlist.push_back(v + "z") ;
    }
    if(vt == BOUNDARY_SCALAR || vt == BOUNDARY_DERIVED_SCALAR) {
      blist.push_back(var) ;
    }
    if(vt == BOUNDARY_VECTOR || vt == BOUNDARY_DERIVED_VECTOR) {
      string v = var ;
      blist.push_back(v + "x ; " + v) ;
      blist.push_back(v + "y") ;
      blist.push_back(v + "z") ;
    }
  }
  // Number of face types (boundary flags)
  ibuf[0] = bc_names.size() ;
  fwrite(ibuf, sizeof(int), 1, OFP) ;

  for(size_t i=0;i<bc_names.size();++i) {
    ibuf[0] = 0 ; // zero, no boundary variables
    if(blist.size() > 0)
      ibuf[0] = 1 ;
    ibuf[1] = 1 ; // Normals should be consistent
    fwrite(ibuf,sizeof(int),2,OFP) ;
    memset(fv,'\0',80) ;
    sprintf(fv, "%s",bc_names[i].c_str()) ;
    fwrite(&fv, sizeof(char), 80, OFP) ;
  }

  // Nodal variables
  ibuf[0] = nlist.size() ;
  fwrite(ibuf,sizeof(int),1,OFP) ;
  for(size_t i=0;i<nlist.size();++i) {
    memset(fv,'\0',80) ;
    sprintf(fv, "%s",nlist[i].c_str()) ;
    fwrite(&fv, sizeof(char), 80, OFP) ;
  }
  // boundary variables 
  ibuf[0] = blist.size() ;
  fwrite(ibuf,sizeof(int),1,OFP) ;
  for(size_t i=0;i<blist.size();++i) {
    memset(fv,'\0',80) ;
    sprintf(fv, "%s",blist[i].c_str()) ;
    fwrite(&fv, sizeof(char), 80, OFP) ;
  }
}
void fv_topo_handler::close() {
  if(first_boundary) {
    int ibuf[2] ;
    ibuf[0] = FV_BNDRY_VARS ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_boundary=false ;
  }
  fclose(OFP) ;
}
void fv_topo_handler::create_mesh_positions(vector3d<float> pos[], int pts) {
  int ibuf[2] ;
  ibuf[0] = FV_NODES ;
  ibuf[1] = pts ;
  fwrite(ibuf,sizeof(int),2,OFP) ;
  for(int i=0;i<pts;++i) {
    float x = pos[i].x;
    fwrite(&x,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<pts;++i) {
    float y = pos[i].y ;
    fwrite(&y,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<pts;++i) {
    float z = pos[i].z;
    fwrite(&z,sizeof(float),1,OFP) ;
  }
}


void fv_topo_handler::create_mesh_elements() {
  int ibuf[5] ;
  ibuf[0] = FV_ELEMENTS ;
  ibuf[1] = ntets ;
  ibuf[2] = nhexs ;
  ibuf[3] = nprsm ;
  ibuf[4] = npyrm ;
  fwrite(ibuf,sizeof(int),5,OFP) ;
}
void fv_topo_handler::write_tets(Array<int,4> tets[], int ntets) {
  static int tet_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                              NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  
  unsigned int elem_header = fv_encode_elem_header(FV_TET_ELEM_ID,
                                                   tet_walls) ;
  for(int i=0;i<ntets;++i) {
    fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
    fwrite(&tets[i][0],sizeof(int),4,OFP) ;
  }
}
void fv_topo_handler::write_pyrm(Array<int,5> pyrm[], int npyrm) {
  static int pyrm_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                              NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  
  unsigned int elem_header = fv_encode_elem_header(FV_PYRA_ELEM_ID,
                                                   pyrm_walls) ;
  for(int i=0;i<npyrm;++i) {
    fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
    fwrite(&pyrm[i][0],sizeof(int),5,OFP) ;
  }
}
void fv_topo_handler::write_prsm(Array<int,6> prsm[], int nprsm) {
  static int prsm_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                              NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  
  unsigned int elem_header = fv_encode_elem_header(FV_PRISM_ELEM_ID,
                                                   prsm_walls) ;
  for(int i=0;i<nprsm;++i) {
    fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
    Array<int,6> prsml ;
    prsml[0] = prsm[i][0] ;
    prsml[1] = prsm[i][3] ;
    prsml[2] = prsm[i][4] ;
    prsml[3] = prsm[i][1] ;
    prsml[4] = prsm[i][5] ;
    prsml[5] = prsm[i][2] ;
    fwrite(&prsml[0],sizeof(int),6,OFP) ;
  }

}
void fv_topo_handler::write_hexs(Array<int,8> hexs[], int nhexs) {
  static int hex_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                              NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };

  unsigned int elem_header = fv_encode_elem_header(FV_HEX_ELEM_ID,
                                                   hex_walls) ;
  for(int i=0;i<nhexs;++i) {
    fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
    Array<int,8> hlocal ;
    hlocal[0] = hexs[i][0] ;
    hlocal[1] = hexs[i][1] ;
    hlocal[2] = hexs[i][3] ;
    hlocal[3] = hexs[i][2] ;
    hlocal[4] = hexs[i][4] ;
    hlocal[5] = hexs[i][5] ;
    hlocal[6] = hexs[i][7] ;
    hlocal[7] = hexs[i][6] ;
    fwrite(&hlocal[0],sizeof(int),8,OFP) ;
  }
  
}

void fv_topo_handler::write_general_cell(int nfaces[], int nnfaces,
                                         int nsides[], int nnsides,
                                         int nodes[], int nnodes) {
  int ibuf[4] ;
  ibuf[0] = FV_ARB_POLY_ELEMENTS ;
  ibuf[1] = nnfaces ;
  fwrite(ibuf,sizeof(int),2,OFP) ;
  int sc = 0 ;
  int nc = 0 ;
  for(int i=0;i<nnfaces;++i) {
    int nf = nfaces[i] ;
    int nnd = 0 ;
    for(int j=0;j<nf;++j)
      nnd+=nsides[sc+j] ;
    vector<int> tmp(nnd) ;
    for(int j=0;j<nnd;++j)
      tmp[j] = nodes[nc+j] ;
    sort(tmp.begin(),tmp.end()) ;
    tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end()) ;
    
    ibuf[0] = nf ;
    ibuf[1] = tmp.size() ;
    ibuf[2] = -1 ;
    fwrite(ibuf,sizeof(int),3,OFP) ;
    for(int j=0;j<nf;++j) {
      ibuf[0] = A_WALL ;
      fwrite(ibuf,sizeof(int),1,OFP) ;
      ibuf[0] = nsides[sc] ;
      fwrite(ibuf,sizeof(int),1,OFP) ;
      fwrite(&nodes[nc],sizeof(int),nsides[sc],OFP) ;
      nc += nsides[sc] ; // next node starting place
      sc++ ; // next face
      ibuf[0] = 0 ;
      fwrite(ibuf,sizeof(int),1,OFP) ;
    }
  }

  if(nc != nnodes) {
    cerr << " something wrong!" << endl ;
  }
  if(sc != nnsides) {
    cerr << "something worng with sides!" << endl ;
  }
}

void fv_topo_handler::close_mesh_elements() {
}
void fv_topo_handler::create_boundary_part(string name,int node_set[],
                                           int npnts) {
  ordinary_faces.clear() ;
  part_nodes.push_back(vector<int>(npnts)) ;
  for(int i=0;i<npnts;++i)
    part_nodes[part_id-1][i] = node_set[i] ;

}
  

void fv_topo_handler::write_quads(Array<int,4> quads[],
                                  int quads_ids[], int nquads) {
  for(int i=0;i<nquads;++i) {
    ordinary_faces.push_back(quads[i]) ;
    elem_ids.push_back(quads_ids[i]) ;
  }
}
void fv_topo_handler::write_trias(Array<int,3> trias[],
                                  int trias_ids[], int ntrias) {
  for(int i=0;i<ntrias;++i) {
    Array<int,4> a ;
    a[0] = trias[i][0] ;
    a[1] = trias[i][1] ;
    a[2] = trias[i][2] ;
    a[3] = 0 ;
    
    ordinary_faces.push_back(a) ;
    elem_ids.push_back(trias_ids[i]) ;
  }
}

void fv_topo_handler::write_general_face(int nside_sizes[],
                                         int nside_ids[], int ngeneral,
                                         int nside_nodes[],
                                         int nside_nodes_size) {
  int ibuf[4] ;
  ibuf[1] = part_id ;
  ibuf[2] = ngeneral+ordinary_faces.size() ;
  for(int i=0;i<ngeneral;++i) {
    elem_ids.push_back(nside_ids[i]) ;
  }
  if(ngeneral == 0) {
    ibuf[0] = FV_FACES ;
    fwrite(ibuf,sizeof(int),3,OFP) ;
    for(size_t i=0;i<ordinary_faces.size();++i) {
      ibuf[0]= part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
      ibuf[1]= part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
      ibuf[2]= part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
      if(ordinary_faces[i][3] == 0) {
        ibuf[3]= 0 ;
      } else {
        ibuf[3]= part_nodes[part_id-1][ordinary_faces[i][3]-1] ;
      }
      fwrite(ibuf,sizeof(int),4,OFP) ;
    }
  } else {
    ibuf[0] = FV_ARB_POLY_FACES ;
    fwrite(ibuf,sizeof(int),3,OFP) ;
    for(size_t i=0;i<ordinary_faces.size();++i) {
      if(ordinary_faces[i][3] == 0) {
        ibuf[0] = 3 ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        ibuf[0]= part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
        ibuf[1]= part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
        ibuf[2]= part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
        fwrite(ibuf,sizeof(int),3,OFP) ;
      } else {
        ibuf[0] = 4 ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        ibuf[0]= part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
        ibuf[1]= part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
        ibuf[2]= part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
        ibuf[3]= part_nodes[part_id-1][ordinary_faces[i][3]-1] ;
        fwrite(ibuf,sizeof(int),4,OFP) ;
      }
    }

    int cnt = 0 ;
    for(int i=0;i<ngeneral;++i) {
      ibuf[0] = nside_sizes[i] ;
      fwrite(ibuf,sizeof(int),1,OFP) ;
      for(int j=0;j<nside_sizes[i];++j) {
        ibuf[0] = part_nodes[part_id-1][nside_nodes[cnt++]-1] ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
      }
    }
  }
}

void fv_topo_handler::close_boundary_part() {
  part_id++ ;
  ordinary_faces.clear() ;
}

void fv_topo_handler::output_nodal_scalar(float val[], int npnts,
                                          string varname) {
  if(first_var) {
    int ibuf[1] ;
    ibuf[0] = FV_VARIABLES ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_var = false ;
  }
  fwrite(val,sizeof(float),npnts,OFP) ;
}

void fv_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               int npnts, string varname) {
  if(first_var) {
    int ibuf[1] ;
    ibuf[0] = FV_VARIABLES ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_var = false ;
  }
  for(int i=0;i<npnts;++i) {
    float d = val[i].x ;
    fwrite(&d,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<npnts;++i) {
    float d = val[i].y ;
    fwrite(&d,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<npnts;++i) {
    float d = val[i].z ;
    fwrite(&d,sizeof(float),1,OFP) ;
  }
}

void fv_topo_handler::output_boundary_scalar(float val[], int node_set[],
                                             int nvals, string varname) {
  if(first_var) {
    int ibuf[1] ;
    ibuf[0] = FV_VARIABLES ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_var = false ;

  }
  if(first_boundary) {
    int ibuf[2] ;
    ibuf[0] = FV_BNDRY_VARS ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_boundary=false ;
  }
  vector<float> valout(elem_ids.size()) ;
  dMap einv ;
  for(size_t i=0;i<elem_ids.size();++i) {
    einv[elem_ids[i]] = i ;
    valout[i] = 0.0 ;
  }
  for(int i=0;i<nvals;++i) {
    valout[einv[node_set[i]]] = val[i] ;
  }
  fwrite(&valout[0],sizeof(float),valout.size(),OFP) ;
}  

void fv_topo_handler::output_boundary_vector(vector3d<float> val[],
                                                  int node_set[],
                                                  int nvals, string varname) {
  if(first_var) {
    int ibuf[1] ;
    ibuf[0] = FV_VARIABLES ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_var = false ;
  }
  if(first_boundary) {
    int ibuf[2] ;
    ibuf[0] = FV_BNDRY_VARS ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    first_boundary=false ;
  }
  vector<float> valout(elem_ids.size()) ;
  dMap einv ;
  for(size_t i=0;i<elem_ids.size();++i) {
    einv[elem_ids[i]] = i ;
    valout[i] = 0.0 ;
  }
  for(int i=0;i<nvals;++i) 
    valout[einv[node_set[i]]] = val[i].x ;
  fwrite(&valout[0],sizeof(float),valout.size(),OFP) ;

  for(size_t i=0;i<elem_ids.size();++i) 
    valout[i] = 0.0 ;
  for(int i=0;i<nvals;++i) 
    valout[einv[node_set[i]]] = val[i].y ;
  fwrite(&valout[0],sizeof(float),valout.size(),OFP) ;

  for(size_t i=0;i<elem_ids.size();++i) 
    valout[i] = 0.0 ;
  for(int i=0;i<nvals;++i) 
    valout[einv[node_set[i]]] = val[i].z ;
  fwrite(&valout[0],sizeof(float),valout.size(),OFP) ;
}  

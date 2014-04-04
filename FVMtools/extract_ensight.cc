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
#include <stdio.h>
#include <strings.h> 
#include <Loci.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <set> 
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::set ;
using std::ofstream ;
using std::ios ;

#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>

void ensight_topo_handler::open(string casename, string iteration ,size_t inpnts,
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
  part_id = 1 ;
  dirname = casename + "_ensight."+iteration ;
  struct stat statbuf ;
  if(stat(dirname.c_str(),&statbuf))
    mkdir(dirname.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file " << dirname << " should be a directory!, rename 'ensight' and start again."
           << endl ;
      exit(-1) ;
    }

  string geo_filename =  casename + ".geo" ;
  string case_filename = dirname + "/" + casename + ".case" ;
  string particle_filename ;
  
  // scan to see if particle geometry is included
  for(size_t i=0;i<variable_types.size();++i) {
    if( (variable_types[i] == PARTICLE_SCALAR) ||
        (variable_types[i] == PARTICLE_VECTOR)) {
      particle_output = true ;
      particle_geo_filename = casename + "_particles.geo" ;
    }
  }
  
  ofstream of(case_filename.c_str(),ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << geo_filename << endl ;
  if(particle_output) {
    of << "measured:  " << particle_geo_filename << endl ;
    particle_geo_filename = dirname + "/" + particle_geo_filename ;
  }

  geo_filename = dirname + "/"+geo_filename ;
  if(variables.size() > 0) {
    of << "VARIABLE" << endl ;
    for(size_t i=0;i<variables.size();++i) {
      const int vt = variable_types[i] ;

      if(vt == NODAL_SCALAR || vt == NODAL_DERIVED || vt == NODAL_MASSFRACTION) {
        of << "scalar per node:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
      if(vt == NODAL_VECTOR) {
        of << "vector per node:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
      if(vt == BOUNDARY_SCALAR) {
        of << "scalar per element:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
      if(vt == BOUNDARY_VECTOR) {
        of << "vector per element:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
      if(vt == PARTICLE_SCALAR) {
        of << "scalar per measured node:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
       if(vt == PARTICLE_VECTOR) {
        of << "vector per measured node:\t " << variables[i] << '\t'
           << variables[i] << endl ;
      }
    }
  }
  of.close() ;
  
  
  OFP = fopen(geo_filename.c_str(),"wb") ;
  if(OFP == NULL) {
    cerr << "unable to open file '" << geo_filename << "' for writing!" ;
    exit(-1) ;
  }
  char tmp_buf[80] ;
  
  memset(tmp_buf, '\0', 80) ; 
  snprintf(tmp_buf,80, "C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "Ensight model geometry description") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "Grid file used is %s", casename.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  
  if(node_id_opt == GIVEN) snprintf(tmp_buf, 80, "node id given") ;
  else snprintf(tmp_buf, 80, "node id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  
  memset(tmp_buf, '\0', 80) ;
  
  if(element_id_opt == GIVEN)  snprintf(tmp_buf,80, "element id given") ;
  else snprintf(tmp_buf, 80,"element id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  
}
void ensight_topo_handler::close() {
  fclose(OFP) ;
}
void ensight_topo_handler::create_mesh_positions(vector3d<float> pos[], size_t pts) {
  positions.reserve(pts) ;
  for(size_t i=0;i<pts;++i)
    positions.push_back(pos[i]) ;
  
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80,"part") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&part_id, sizeof(int), 1, OFP) ;
  part_id++ ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "Entire") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  int pnts =  pts ;
  fwrite(&pnts, sizeof(int), 1, OFP) ;

  if(node_id_opt == GIVEN || node_id_opt == IGNORE){
    for(size_t i = 0; i < pts; i++){
      int nid = i+1;
      fwrite(&nid, sizeof(int), 1, OFP) ;
    }
  }

  for(size_t i=0;i<pts;++i) {
    float x = pos[i].x;
    fwrite(&x,sizeof(float),1,OFP) ;
  }
  for(size_t i=0;i<pts;++i) {
    float y = pos[i].y ;
    fwrite(&y,sizeof(float),1,OFP) ;
  }
  for(size_t i=0;i<pts;++i) {
    float z = pos[i].z;
    fwrite(&z,sizeof(float),1,OFP) ;
  }
}

void ensight_topo_handler::write_tets(Array<int,4> tets[], size_t ntets, int block,int nblocks, size_t tottets) {
  if(tottets > 0) {
    if(block == 0 && (element_id_opt != GIVEN && element_id_opt != IGNORE)) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "tetra4") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = tottets ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
    }
    
    if(ntets > 0)
      fwrite(tets,sizeof(Array<int,4>),ntets,OFP) ;
  }
}
void ensight_topo_handler::write_tets_ids( int tets_ids[], size_t ntets, int block,int nblocks, size_t tottets) {
  if(element_id_opt == GIVEN || element_id_opt == IGNORE) {
    if(tottets > 0) {
      if(block == 0) {
        char tmp_buf[80] ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf, 80, "tetra4") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	int tot = tottets ;
        fwrite(&tot, sizeof(int), 1, OFP) ;
      }
    
    if(ntets > 0)
      fwrite(tets_ids,sizeof(int),ntets,OFP) ;
    }
  }
}

void ensight_topo_handler::write_pyrm(Array<int,5> pyrm[], size_t npyrm,int block, int nblocks, size_t totpyrm) {
  if(totpyrm > 0) {
    if(block==0&& (element_id_opt != GIVEN && element_id_opt != IGNORE) ) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "pyramid5") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = totpyrm ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
    }
    if(npyrm > 0)
      fwrite(pyrm, sizeof(Array<int,5>), npyrm, OFP) ;
  }
}

void ensight_topo_handler::write_pyrm_ids(int pyrm_ids[], size_t npyrm,int block, int nblocks, size_t totpyrm) {
 if(element_id_opt == GIVEN || element_id_opt == IGNORE) {
  if(totpyrm > 0) {
    if(block==0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "pyramid5") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = totpyrm ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
    }
    if(npyrm > 0)
      fwrite(pyrm_ids, sizeof(int), npyrm, OFP) ;
  }
 }
}


void ensight_topo_handler::write_prsm(Array<int,6> prsm[], size_t nprsm,int block, int nblocks, size_t totprsm) {
  if(totprsm > 0) {
    if(block==0 && (element_id_opt != GIVEN && element_id_opt != IGNORE) ) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "penta6") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = totprsm ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
    }
    if(nprsm > 0)
      fwrite(prsm, sizeof(Array<int,6> ), nprsm, OFP) ;
  }
}

void ensight_topo_handler::write_prsm_ids(int prsm_ids[], size_t nprsm,int block, int nblocks, size_t totprsm) {
   if(element_id_opt == GIVEN || element_id_opt == IGNORE) {
     if(totprsm > 0) {
      if(block==0) {
        char tmp_buf[80] ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf, 80, "penta6") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	int tot = totprsm ;
        fwrite(&tot, sizeof(int), 1, OFP) ;
      }
      if(nprsm > 0)
        fwrite(prsm_ids, sizeof(int), nprsm, OFP) ;
     
     }
   }
}

void ensight_topo_handler::write_hexs(Array<int,8> hexs[], size_t nhexs, int block, int nblocks, size_t tothexs) {
  if(tothexs > 0) {
    if(block == 0 &&(element_id_opt != GIVEN && element_id_opt != IGNORE)) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "hexa8") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = tothexs ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
    }
    if(nhexs > 0)
      fwrite(hexs, sizeof(Array<int,8>), nhexs, OFP) ;
  }
}

void ensight_topo_handler::write_hexs_ids( int hexs_ids[], size_t nhexs, int block, int nblocks, size_t tothexs) {
  if(element_id_opt == GIVEN || element_id_opt == IGNORE) {
    if(tothexs > 0) {
      if(block == 0) {
        char tmp_buf[80] ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "hexa8") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
	int tot = tothexs ;
        fwrite(&tot, sizeof(int), 1, OFP) ;
      }
    if(nhexs > 0)
      fwrite(hexs_ids, sizeof(int), nhexs, OFP) ;
    }
  }
}



void ensight_topo_handler::write_general_cell(int nfaces[],size_t nnfaces,
                                               int nsides[], size_t nnsides,
                                               int nodes[], size_t nnodes) {
  if(nnfaces > 0) {
    if(element_id_opt != GIVEN && element_id_opt != IGNORE){
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "nfaced") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nnf = nnfaces ;
      fwrite(&nnf, sizeof(int), 1, OFP) ;
    }
    fwrite(nfaces, sizeof(int), nnfaces, OFP) ;
    fwrite(nsides, sizeof(int),nnsides,OFP) ;
    fwrite(nodes,sizeof(int),nnodes,OFP) ;
  }
}

void ensight_topo_handler::write_general_cell_ids(int nfaces_ids[], size_t nnfaces) {
  if(element_id_opt == GIVEN || element_id_opt == IGNORE) {
    if(nnfaces > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "nfaced") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nnf = nnfaces ;
      fwrite(&nnf, sizeof(int), 1, OFP) ;
      fwrite(nfaces_ids, sizeof(int), nnfaces, OFP) ;
    }
  }
}

void ensight_topo_handler::create_boundary_part(string name,int node_set[],
                                                size_t npnts) {
  part_nodes.push_back(vector<int>(npnts)) ;
  part_tria_ids.push_back(vector<int>(0)) ;
  part_quad_ids.push_back(vector<int>(0)) ;
  part_nside_ids.push_back(vector<int>(0)) ;
  
  for(size_t i=0;i<npnts;++i)
    part_nodes[part_id-2][i] = node_set[i] ;

  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&part_id, sizeof(int), 1, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "%s", name.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  int npt = npnts ;
  fwrite(&npt, sizeof(int), 1, OFP) ;

  if(node_id_opt == GIVEN || node_id_opt == IGNORE){
    fwrite(&node_set[0], sizeof(int), npnts,OFP) ;
  }
  for(size_t i=0;i<npnts;++i) {
    float x = positions[node_set[i]-1].x;
    fwrite(&x,sizeof(float),1,OFP) ;
  }
  for(size_t i=0;i<npnts;++i) {
    float y = positions[node_set[i]-1].y ;
    fwrite(&y,sizeof(float),1,OFP) ;
  }
  for(size_t i=0;i<npnts;++i) {
    float z = positions[node_set[i]-1].z;
    fwrite(&z,sizeof(float),1,OFP) ;
  }
}
  

void ensight_topo_handler::write_quads(Array<int,4> quads[],
                                       int quads_ids[], size_t nquads) {
  if(nquads > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "quad4") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    int nq = nquads ;
    fwrite(&nq, sizeof(int),1,OFP) ;
    if(element_id_opt == GIVEN || element_id_opt == IGNORE ){
      fwrite(&quads_ids[0], sizeof(int),nquads,OFP) ;
    }
    fwrite(&quads[0],sizeof(Array<int,4>),nquads,OFP) ;
    vector<int> qid(nquads) ;
    for(size_t i=0;i<nquads;++i) 
      qid[i] = quads_ids[i] ;
    part_quad_ids[part_id-2].swap(qid) ;
  }
}
void ensight_topo_handler::write_trias(Array<int,3> trias[],
                                       int trias_ids[], size_t ntrias) {
  if(ntrias > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "tria3") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    int nt = ntrias ;
    fwrite(&nt, sizeof(int),1,OFP) ;
    if(element_id_opt == GIVEN || element_id_opt == IGNORE){
      fwrite(&trias_ids[0], sizeof(int),ntrias,OFP) ;
    }
    fwrite(&trias[0],sizeof(Array<int,3>),ntrias,OFP) ;
    vector<int> qid(ntrias) ;
    for(size_t i=0;i<ntrias;++i) 
      qid[i] = trias_ids[i] ;
    part_tria_ids[part_id-2].swap(qid) ;
  }
}

void ensight_topo_handler::write_general_face(int nside_sizes[],
                                         int nside_ids[], size_t ngeneral,
                                         int nside_nodes[],
                                         size_t nside_nodes_size) {
  if(ngeneral > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "nsided") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    int ngen = ngeneral ;
    fwrite(&ngen, sizeof(int),1,OFP) ;
    if(element_id_opt == GIVEN ||element_id_opt == IGNORE ){
      fwrite(&nside_ids[0], sizeof(int),ngeneral,OFP) ;
    }
    size_t tot = 0 ;
    for(size_t i=0;i<ngeneral;++i)
      tot += nside_sizes[i] ;
    if(nside_nodes_size != tot) {
      cerr << "mismatch in node size and faces size " << nside_nodes_size
           << " was " << tot << endl ;
    }
    fwrite(&nside_sizes[0],sizeof(int),ngeneral,OFP) ;
    fwrite(&nside_nodes[0],sizeof(int),nside_nodes_size,OFP) ;
    vector<int> qid(ngeneral) ;
    for(size_t i=0;i<ngeneral;++i) 
      qid[i] = nside_ids[i] ;
    part_nside_ids[part_id-2].swap(qid) ;
  }
}

void ensight_topo_handler::close_boundary_part() {
  part_id++ ;
}

void ensight_topo_handler::output_nodal_scalar(float val[], size_t npnts,
                                               string varname) {
  string filename = dirname + '/' + varname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename << "' for writing variable info!"
         << endl ;

    return ;
  }
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  int tmp = 1 ;
  fwrite(&tmp, sizeof(int), 1, FP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  fwrite(val,sizeof(float),npnts,FP) ;
  for(size_t i=0;i<part_nodes.size();++i) {
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int tmp = i+2 ;
    fwrite(&tmp, sizeof(int), 1, FP) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    for(size_t j=0;j<part_nodes[i].size();++j) {
      int nd = part_nodes[i][j]-1 ;
      fwrite(&val[nd],sizeof(float),1,FP) ;
    }
  }
  fclose(FP) ;
    
}

void ensight_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               size_t npnts, string varname) {
  string filename = dirname + '/' + varname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename << "' for writing variable info!"
         << endl ;

    return ;
  }
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  int tmp = 1 ;
  fwrite(&tmp, sizeof(int), 1, FP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  for(size_t i=0;i<npnts;++i) {
    float d = val[i].x ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
  for(size_t i=0;i<npnts;++i) {
    float d = val[i].y ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
  for(size_t i=0;i<npnts;++i) {
    float d = val[i].z ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
    
  for(size_t i=0;i<part_nodes.size();++i) {
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int tmp = i+2 ;
    fwrite(&tmp, sizeof(int), 1, FP) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    for(size_t j=0;j<part_nodes[i].size();++j) {
      int nd = part_nodes[i][j]-1 ;
      fwrite(&val[nd].x,sizeof(float),1,FP) ;
    }
    for(size_t j=0;j<part_nodes[i].size();++j) {
      int nd = part_nodes[i][j]-1 ;
      fwrite(&val[nd].y,sizeof(float),1,FP) ;
    }
    for(size_t j=0;j<part_nodes[i].size();++j) {
      int nd = part_nodes[i][j]-1 ;
      fwrite(&val[nd].z,sizeof(float),1,FP) ;
    }
  }

  fclose(FP) ;
  
}

void ensight_topo_handler::output_boundary_scalar(float val[], int node_set[],
                                                  size_t nvals, string varname) {
  map<int,int> nmap ;
  map<int,int>::const_iterator mi ;

  for(size_t i=0;i<nvals;++i)
    nmap[node_set[i]] = i ;
  
  string filename = dirname + '/' + varname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename << "' for writing variable info!"
         << endl ;

    return ;
  }
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80 ,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;


  int nparts = part_tria_ids.size() ;
  for(int i=0;i<nparts;++i) {
    int id = -1 ;
    if(part_tria_ids[i].size()>0)
      id = part_tria_ids[i][0] ;
    if(part_quad_ids[i].size()>0)
      id = part_quad_ids[i][0] ;
    if(part_nside_ids[i].size()>0)
      id = part_nside_ids[i][0] ;
    bool inpart = (nmap.find(id) != nmap.end()) ;
    if(!inpart) {
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int partid = i+2 ;
    fwrite(&partid,sizeof(int),1,FP) ;
    if(part_tria_ids[i].size() > 0 ) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"tria3") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_tria_ids[i].size();++j) {
        mi = nmap.find(part_tria_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second] ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
    }
    if(part_quad_ids[i].size() > 0) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"quad4") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_quad_ids[i].size();++j) {
        mi = nmap.find(part_quad_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second] ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
    }
    if(part_nside_ids[i].size() > 0) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_nside_ids[i].size();++j) {
        mi = nmap.find(part_nside_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second] ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
    }
  }
  fclose(FP) ;
}  

void ensight_topo_handler::output_boundary_vector(vector3d<float> val[],
                                                  int node_set[],
                                                  size_t nvals, string varname) {
  map<int,int> nmap ;
  map<int,int>::const_iterator mi ;

  for(size_t i=0;i<nvals;++i)
    nmap[node_set[i]] = i ;
  
  string filename = dirname + '/' + varname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename << "' for writing variable info!"
         << endl ;

    return ;
  }
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;


  int nparts = part_tria_ids.size() ;
  for(int i=0;i<nparts;++i) {
    int id = -1 ;
    if(part_tria_ids[i].size()>0)
      id = part_tria_ids[i][0] ;
    if(part_quad_ids[i].size()>0)
      id = part_quad_ids[i][0] ;
    if(part_nside_ids[i].size()>0)
      id = part_nside_ids[i][0] ;
    bool inpart = (nmap.find(id) != nmap.end()) ;
    if(!inpart) {
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int partid = i+2 ;
    fwrite(&partid,sizeof(int),1,FP) ;
    if(part_tria_ids[i].size() > 0 ) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"tria3") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_tria_ids[i].size();++j) {
        mi = nmap.find(part_tria_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].x ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_tria_ids[i].size();++j) {
        mi = nmap.find(part_tria_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].y ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_tria_ids[i].size();++j) {
        mi = nmap.find(part_tria_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].z ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
    }
    if(part_quad_ids[i].size() > 0) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"quad4") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_quad_ids[i].size();++j) {
        mi = nmap.find(part_quad_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].x ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_quad_ids[i].size();++j) {
        mi = nmap.find(part_quad_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].y ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_quad_ids[i].size();++j) {
        mi = nmap.find(part_quad_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].z ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }

    }
    if(part_nside_ids[i].size() > 0) {
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80,"nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, FP) ;
      for(size_t j=0;j!=part_nside_ids[i].size();++j) {
        mi = nmap.find(part_nside_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].x ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_nside_ids[i].size();++j) {
        mi = nmap.find(part_nside_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].y ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
      for(size_t j=0;j!=part_nside_ids[i].size();++j) {
        mi = nmap.find(part_nside_ids[i][j]) ;
        float tmp = 0 ;
        if(mi!=nmap.end()) 
          tmp = val[mi->second].z ;
        fwrite(&tmp,sizeof(float),1,FP) ;
      }
    }
  }
  fclose(FP) ;
}

void
ensight_topo_handler::create_particle_positions
(vector3d<float> pos[], size_t pts, size_t maxp) {

  int m = pts ;
  if(maxp >= 0 && maxp < pts)
    m = maxp ;

  int jump ;
  if(m>0)
    jump = pts / m ;
  else
    jump = 1 ;
  
  FILE *FP = 0 ;
  FP = fopen(particle_geo_filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << particle_geo_filename
         << "' for writing particle geometry info!" << endl ;
    return ;
  }

  char tmp_buf[80] ;

  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80,"C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "particle positions") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "particle coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  // number of points
  fwrite(&m, sizeof(int), 1, FP) ;
  // point ids
  for(int i=1;i<m+1;++i)
    fwrite(&i, sizeof(int), 1, FP) ;
  // point coordinates
  int k = 0 ;
  for(int i=0;i<m;++i,k+=jump) {
    float x = pos[k].x ;
    float y = pos[k].y ;
    float z = pos[k].z ;
    fwrite(&x, sizeof(float), 1, FP) ;
    fwrite(&y, sizeof(float), 1, FP) ;
    fwrite(&z, sizeof(float), 1, FP) ;
  }
  
  fclose(FP) ;
}

void
ensight_topo_handler::output_particle_scalar(float val[], size_t np,
                                             size_t maxp, string valname) {
  int m = np ;
  if(maxp >= 0 && maxp < np)
    m = maxp ;

  int jump ;
  if(m>0)
    jump = np / m ;
  else
    jump = 1 ;

  string filename = dirname + '/' + valname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename
         << "' for writing variable info!" << endl ;
    return ;
  }

  char tmp_buf[80] ;

  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "Per particle scalar: %s", valname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  int k = 0 ;
  for(int i=0;i<m;++i,k+=jump)
    fwrite(&val[k], sizeof(float), 1, FP) ;

  fclose(FP) ;
}

void
ensight_topo_handler::output_particle_vector(vector3d<float> val[], size_t np,
                                             size_t maxp, string valname) {
  int m = np ;
  if(maxp >= 0 && maxp < np)
    m = maxp ;

  int jump ;
  if(m>0)
    jump = np / m ;
  else
    jump = 1 ;

  string filename = dirname + '/' + valname ;
  FILE *FP = 0 ;
  FP = fopen(filename.c_str(), "wb") ;
  if(FP==0) {
    cerr << "can't open file '" << filename
         << "' for writing variable info!" << endl ;
    return ;
  }

  char tmp_buf[80] ;

  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "Per particle vector: %s", valname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  int k = 0 ;
  for(int i=0;i<m;++i,k+=jump) {
    float x = val[k].x ;
    float y = val[k].y ;
    float z = val[k].z ;
    
    fwrite(&x, sizeof(float), 1, FP) ;
    fwrite(&y, sizeof(float), 1, FP) ;
    fwrite(&z, sizeof(float), 1, FP) ;
  }

  fclose(FP) ;
}

void ensightPartConverter::exportPostProcessorFiles(string casename,
						    string iteration) {
  string dirname = casename + "_ensight."+iteration ;
  struct stat statbuf ;
  if(stat(dirname.c_str(),&statbuf))
    mkdir(dirname.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file " << dirname << " should be a directory!, rename 'ensight' and start again."
           << endl ;
      exit(-1) ;
    }
  
  string geo_filename =  casename + ".geo" ;
  string case_filename = dirname + "/" + casename + ".case" ;
  set<string> nodal_scalars ;
  set<string> nodal_vectors ;
  set<string> element_scalars ;
  set<string> element_vectors ;
  for(size_t i =0;i<surfacePartList.size();++i) {
    vector<string> nscalars = surfacePartList[i].getNodalScalarVars() ;
    for(size_t j=0;j<nscalars.size();++j) 
      nodal_scalars.insert(nscalars[j]) ;
    vector<string> nvectors = surfacePartList[i].getNodalVectorVars() ;
    for(size_t j=0;j<nvectors.size();++j) 
      nodal_vectors.insert(nvectors[j]) ;
    vector<string> escalars = surfacePartList[i].getElementScalarVars() ;
    for(size_t j=0;j<escalars.size();++j) 
      element_scalars.insert(escalars[j]) ;
    vector<string> evectors = surfacePartList[i].getElementVectorVars() ;
    for(size_t j=0;j<evectors.size();++j) 
      element_vectors.insert(evectors[j]) ;
  }


  //write out case file
  ofstream of(case_filename.c_str(),ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << geo_filename << endl ;
  
  geo_filename = dirname + "/"+geo_filename ;
  
  if(nodal_scalars.size()+nodal_vectors.size()+
     element_scalars.size()+element_vectors.size() > 0) {
    of << "VARIABLE" << endl ;
    set<string>::const_iterator si ;
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
      of << "scalar per node:\t " << *si << '\t' << *si << endl ;
    }
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      of << "vector per node:\t " << *si << '\t' << *si << endl ;
    }
    for(si=element_scalars.begin();si!=element_scalars.end();++si) {
      of << "scalar per element:\t " << *si << '\t' << *si << endl ;
    }
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      of << "vector per element:\t " << *si << '\t' << *si << endl ;
    }
  }
  of.close() ;
  
  FILE *OFP ;
  OFP = fopen(geo_filename.c_str(),"wb") ;
  if(OFP == NULL) {
    cerr << "unable to open file '" << geo_filename << "' for writing!" ;
    exit(-1) ;
  }
  char tmp_buf[80] ;
  
  memset(tmp_buf, '\0', 80) ; 
  snprintf(tmp_buf,80, "C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "Ensight model geometry description") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "Grid file used is %s", casename.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "node id off") ;  
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80,"element id off") ;  
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;

  int part_id = 0 ;
  vector<int> partnums(surfacePartList.size()) ;
  // Loop over parts, write out their geometry
  for(size_t i =0;i<surfacePartList.size();++i) {
    part_id++ ;
    partnums[i] = part_id ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&part_id, sizeof(int), 1, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    string name = surfacePartList[i].getPartName();
    snprintf(tmp_buf,80, "%s", name.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    int npt = surfacePartList[i].getNumNodes() ;
    fwrite(&npt, sizeof(int), 1, OFP) ;
    vector<vector3d<float> > pos ;
    surfacePartList[i].getPos(pos) ;
    for(int j=0;j<npt;++j) {
      float x = pos[j].x ;
      fwrite(&x,sizeof(float),1,OFP) ;
    }
    for(int j=0;j<npt;++j) {
      float y = pos[j].y ;
      fwrite(&y,sizeof(float),1,OFP) ;
    }
    for(int j=0;j<npt;++j) {
      float z = pos[j].z ;
      fwrite(&z,sizeof(float),1,OFP) ;
    }
    int nquads = surfacePartList[i].getNumQuads() ;
    if(nquads > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "quad4") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nq = nquads ;
      fwrite(&nq, sizeof(int),1,OFP) ;
      vector<Array<int,4> > quads ;
      surfacePartList[i].getQuads(quads) ;
      fwrite(&quads[0],sizeof(Array<int,4>),nquads,OFP) ;
    }
    int ntrias = surfacePartList[i].getNumTrias() ;
    if(ntrias > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "tria3") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nt = ntrias ;
      fwrite(&nt, sizeof(int),1,OFP) ;
      vector<Array<int,3> > trias ;
      surfacePartList[i].getTrias(trias) ;
      fwrite(&trias[0],sizeof(Array<int,3>),ntrias,OFP) ;
    }    
    int ngeneral = surfacePartList[i].getNumGenfc() ;
    if(ngeneral > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int ngen = ngeneral ;
      fwrite(&ngen, sizeof(int),1,OFP) ;
      
      vector<int> nside_sizes,nside_nodes ;
      surfacePartList[i].getGenf(nside_sizes,nside_nodes) ;
      
      int tot = 0 ;
      for(int j=0;j<ngeneral;++j)
	tot += nside_sizes[j] ;
      int nside_nodes_size = nside_nodes.size() ;
      if(nside_nodes_size != tot) {
	cerr << "mismatch in node size and faces size " << nside_nodes_size
	     << " was " << tot << endl ;
      }
      fwrite(&nside_sizes[0],sizeof(int),ngeneral,OFP) ;
      fwrite(&nside_nodes[0],sizeof(int),nside_nodes_size,OFP) ;
    }
  }

  fclose(OFP) ;
  
  // Finished writing out the the geo file, now write out the variables
  set<string>::const_iterator si ;
  // write out nodal scalars
  for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;

      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i].hasNodalScalarVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<float> vals ;
	surfacePartList[i].getNodalScalar(varname,vals) ;
	fwrite(&vals[0],sizeof(float),vals.size(),FP) ;
      }
    }
    fclose(FP) ;
  }
  // write out nodal vector variables
  for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;

      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      memset(tmp_buf, '\0', 80) ;
      if(surfacePartList[i].hasElementVectorVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<vector3d<float> > vals ;
	surfacePartList[i].getNodalVector(varname,vals) ;
	int nvals = vals.size() ;
	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].x,sizeof(float),1,FP) ;

	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].y,sizeof(float),1,FP) ;

	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].z,sizeof(float),1,FP) ;
      }
    }
    fclose(FP) ;
  }
  // write out element scalars
  for(si=element_scalars.begin();si!=element_scalars.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i].hasElementScalarVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
        vector<float> qvals, tvals, gvals ;
        surfacePartList[i].getElementScalar(varname,qvals,tvals,gvals) ;

        if(qvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "quad4") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&qvals[0],sizeof(float),qvals.size(),FP) ;
        }
        if(tvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "tria3") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&tvals[0],sizeof(float),tvals.size(),FP) ;
        }
        if(gvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "nsided") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&gvals[0],sizeof(float),gvals.size(),FP) ;
        }
      }
    }
    fclose(FP) ;
  }
  // write out element vectors
  for(si=element_vectors.begin();si!=element_vectors.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i].hasElementVectorVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
        vector<vector3d<float> > qvals, tvals, gvals ;
        surfacePartList[i].getElementVector(varname,qvals,tvals,gvals) ;
        
        if(qvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "quad4") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].z,sizeof(float),1,FP) ;
        }
        if(tvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "tria3") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].z,sizeof(float),1,FP) ;
        }
        if(gvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "nsided") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].z,sizeof(float),1,FP) ;
        }
      }
    }
    fclose(FP) ;
  }

}
// ... the end ...



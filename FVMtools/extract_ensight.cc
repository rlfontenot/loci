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

void ensight_topo_handler::open(string casename, string iteration ,int inpnts,
                                int intets, int inprsm, int inpyrm,
                                int inhexs, int ingen,
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
  sprintf(tmp_buf, "C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Ensight model geometry description") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Grid file used is %s", casename.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "node id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "element id off") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  
}
void ensight_topo_handler::close() {
  fclose(OFP) ;
}
void ensight_topo_handler::create_mesh_positions(vector3d<float> pos[], int pts) {
  positions.reserve(pts) ;
  for(int i=0;i<pts;++i)
    positions.push_back(pos[i]) ;
  
  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&part_id, sizeof(int), 1, OFP) ;
  part_id++ ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "Entire") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&pts, sizeof(int), 1, OFP) ;

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

void ensight_topo_handler::write_tets(Array<int,4> tets[], int ntets) {
  if(ntets > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "tetra4") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&ntets, sizeof(int), 1, OFP) ;

    fwrite(tets,sizeof(Array<int,4>),ntets,OFP) ;
  }
}
void ensight_topo_handler::write_pyrm(Array<int,5> pyrm[], int npyrm) {
  if(npyrm > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "pyramid5") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&npyrm, sizeof(int), 1, OFP) ;
    fwrite(pyrm, sizeof(Array<int,5>), npyrm, OFP) ;
  }
}
void ensight_topo_handler::write_prsm(Array<int,6> prsm[], int nprsm) {
  if(nprsm > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "penta6") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&nprsm, sizeof(int), 1, OFP) ;
    fwrite(prsm, sizeof(Array<int,6> ), nprsm, OFP) ;
  }
}
void ensight_topo_handler::write_hexs(Array<int,8> hexs[], int nhexs) {
  if(nhexs > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "hexa8") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&nhexs, sizeof(int), 1, OFP) ;
    fwrite(hexs, sizeof(Array<int,8>), nhexs, OFP) ;
  }
}

void ensight_topo_handler::write_general_cell(int nfaces[], int nnfaces,
                                               int nsides[], int nnsides,
                                               int nodes[], int nnodes) {
  if(nfaces > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "nfaced") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&nnfaces, sizeof(int), 1, OFP) ;
    fwrite(nfaces, sizeof(int), nnfaces, OFP) ;
    fwrite(nsides, sizeof(int),nnsides,OFP) ;
    fwrite(nodes,sizeof(int),nnodes,OFP) ;
  }
}

void ensight_topo_handler::create_boundary_part(string name,int node_set[],
                                                int npnts) {
  part_nodes.push_back(vector<int>(npnts)) ;
  part_tria_ids.push_back(vector<int>(0)) ;
  part_quad_ids.push_back(vector<int>(0)) ;
  part_nside_ids.push_back(vector<int>(0)) ;
  
  for(int i=0;i<npnts;++i)
    part_nodes[part_id-2][i] = node_set[i] ;

  char tmp_buf[80] ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&part_id, sizeof(int), 1, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "%s", name.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  fwrite(&npnts, sizeof(int), 1, OFP) ;
    
  for(int i=0;i<npnts;++i) {
    float x = positions[node_set[i]-1].x;
    fwrite(&x,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<npnts;++i) {
    float y = positions[node_set[i]-1].y ;
    fwrite(&y,sizeof(float),1,OFP) ;
  }
  for(int i=0;i<npnts;++i) {
    float z = positions[node_set[i]-1].z;
    fwrite(&z,sizeof(float),1,OFP) ;
  }
}
  

void ensight_topo_handler::write_quads(Array<int,4> quads[],
                                       int quads_ids[], int nquads) {
  if(nquads > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "quad4") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&nquads, sizeof(int),1,OFP) ;
    fwrite(&quads[0],sizeof(Array<int,4>),nquads,OFP) ;
    vector<int> qid(nquads) ;
    for(int i=0;i<nquads;++i) 
      qid[i] = quads_ids[i] ;
    part_quad_ids[part_id-2].swap(qid) ;
  }
}
void ensight_topo_handler::write_trias(Array<int,3> trias[],
                                       int trias_ids[], int ntrias) {
  if(ntrias > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "tria3") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&ntrias, sizeof(int),1,OFP) ;
    fwrite(&trias[0],sizeof(Array<int,3>),ntrias,OFP) ;
    vector<int> qid(ntrias) ;
    for(int i=0;i<ntrias;++i) 
      qid[i] = trias_ids[i] ;
    part_tria_ids[part_id-2].swap(qid) ;
  }
}

void ensight_topo_handler::write_general_face(int nside_sizes[],
                                         int nside_ids[], int ngeneral,
                                         int nside_nodes[],
                                         int nside_nodes_size) {
  if(ngeneral > 0) {
    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "nsided") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&ngeneral, sizeof(int),1,OFP) ;
    int tot = 0 ;
    for(int i=0;i<ngeneral;++i)
      tot += nside_sizes[i] ;
    if(nside_nodes_size != tot) {
      cerr << "mismatch in node size and faces size " << nside_nodes_size
           << " was " << tot << endl ;
    }
    fwrite(&nside_sizes[0],sizeof(int),ngeneral,OFP) ;
    fwrite(&nside_nodes[0],sizeof(int),nside_nodes_size,OFP) ;
    vector<int> qid(ngeneral) ;
    for(int i=0;i<ngeneral;++i) 
      qid[i] = nside_ids[i] ;
    part_nside_ids[part_id-2].swap(qid) ;
  }
}

void ensight_topo_handler::close_boundary_part() {
  part_id++ ;
}

void ensight_topo_handler::output_nodal_scalar(float val[], int npnts,
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
  sprintf(tmp_buf,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  int tmp = 1 ;
  fwrite(&tmp, sizeof(int), 1, FP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  fwrite(val,sizeof(float),npnts,FP) ;
  for(size_t i=0;i<part_nodes.size();++i) {
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int tmp = i+2 ;
    fwrite(&tmp, sizeof(int), 1, FP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    for(size_t j=0;j<part_nodes[i].size();++j) {
      int nd = part_nodes[i][j]-1 ;
      fwrite(&val[nd],sizeof(float),1,FP) ;
    }
  }
  fclose(FP) ;
    
}

void ensight_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               int npnts, string varname) {
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
  sprintf(tmp_buf,"variable : %s",varname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "part") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  int tmp = 1 ;
  fwrite(&tmp, sizeof(int), 1, FP) ;
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "coordinates") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  for(int i=0;i<npnts;++i) {
    float d = val[i].x ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
  for(int i=0;i<npnts;++i) {
    float d = val[i].y ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
  for(int i=0;i<npnts;++i) {
    float d = val[i].z ;
    fwrite(&d,sizeof(float),1,FP) ;
  }
    
  for(size_t i=0;i<part_nodes.size();++i) {
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int tmp = i+2 ;
    fwrite(&tmp, sizeof(int), 1, FP) ;
    memset(tmp_buf, '\0', 80) ;
    sprintf(tmp_buf, "coordinates") ;
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
                                                  int nvals, string varname) {
  map<int,int> nmap ;
  map<int,int>::const_iterator mi ;

  for(int i=0;i<nvals;++i)
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
  sprintf(tmp_buf,"variable : %s",varname.c_str()) ;
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
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int partid = i+2 ;
    fwrite(&partid,sizeof(int),1,FP) ;
    if(part_tria_ids[i].size() > 0 ) {
      memset(tmp_buf, '\0', 80) ;
      sprintf(tmp_buf,"tria3") ;
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
      sprintf(tmp_buf,"quad4") ;
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
      sprintf(tmp_buf,"nsided") ;
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
                                                  int nvals, string varname) {
  map<int,int> nmap ;
  map<int,int>::const_iterator mi ;

  for(int i=0;i<nvals;++i)
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
  sprintf(tmp_buf,"variable : %s",varname.c_str()) ;
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
    sprintf(tmp_buf, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    int partid = i+2 ;
    fwrite(&partid,sizeof(int),1,FP) ;
    if(part_tria_ids[i].size() > 0 ) {
      memset(tmp_buf, '\0', 80) ;
      sprintf(tmp_buf,"tria3") ;
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
      sprintf(tmp_buf,"quad4") ;
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
      sprintf(tmp_buf,"nsided") ;
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
(vector3d<float> pos[], int pts, int maxp) {

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
  sprintf(tmp_buf,"C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;
  
  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "particle positions") ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  memset(tmp_buf, '\0', 80) ;
  sprintf(tmp_buf, "particle coordinates") ;
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
ensight_topo_handler::output_particle_scalar(float val[], int np,
                                             int maxp, string valname) {
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
  sprintf(tmp_buf, "Per particle scalar: %s", valname.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, FP) ;

  int k = 0 ;
  for(int i=0;i<m;++i,k+=jump)
    fwrite(&val[k], sizeof(float), 1, FP) ;

  fclose(FP) ;
}

void
ensight_topo_handler::output_particle_vector(vector3d<float> val[], int np,
                                             int maxp, string valname) {
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
  sprintf(tmp_buf, "Per particle vector: %s", valname.c_str()) ;
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

// ... the end ...



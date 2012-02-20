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

void vtk_topo_handler::open(string casename, string iteration ,int inpnts,
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
  nvars = 0 ;
  filename = "vtk_"+casename+"_"+iteration+".vtu" ;
  int ncells = ntets+nprsm+npyrm+nhexs ;
  string varstring = "\"x\", \"y\", \"z\"" ;
  for(size_t i=0;i<variables.size();++i) {
    if(variable_types[i] == NODAL_SCALAR ||
       variable_types[i] == NODAL_DERIVED ||
       variable_types[i] == NODAL_MASSFRACTION) {
      nvars++ ;
    } else if(variable_types[i] == NODAL_VECTOR) {
      nvars += 3;
    } else {
      cerr << "Variable type for " << variables[i] << " is not currently supported by vtk extractor" << endl ;
    }
  }
  data_store = new float[nvars * npnts];
  FILE * fid = fopen(filename.c_str(), "w"); 
  fprintf(fid,"<?xml version='1.0'?>\n");
  fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
  fprintf(fid,"  <UnstructuredGrid>\n");
  fprintf(fid,"    <Piece NumberOfPoints='%d' NumberOfCells='%d'>\n",npnts,ncells);
  fprintf(fid,"      <Points>\n");
  fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%d'/>\n",Offset);
  fprintf(fid,"      </Points>\n");
  Offset += 3 * npnts * sizeof(float) + sizeof(int) ;
  fclose(fid);
}
void vtk_topo_handler::close_mesh_elements() {
  int ncells = bricks.size() ;
  int *nm = &bricks[0][0];
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"      <Cells>\n");
  fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  int counter=0 ;
  conn = new int[ntets*4 + nhexs*8 + nprsm*6 + npyrm*5] ;
  cell_offsets = new int[ntets+nhexs+nprsm+npyrm] ;
  cell_types = new unsigned char[ntets+nhexs+nprsm+npyrm] ;
  unsigned int tet = 10, pyrm = 14, prsm = 13, hex = 12;
  for(int i = 0; i < ncells; i++) {
    int j = i*8 ;
    for (int k=0;k<4;k++) conn[counter+k] = nm[j+k];
    if      (nm[j+7] > -1) {
      for (int k=4;k<8;k++) conn[counter+k] = nm[j+k];
      counter += 8 ;
      cell_types[i] = (unsigned char)hex ;
    }
    else if (nm[j+5] > -1) {
      conn[counter+4] = nm[j+4];
      conn[counter+5] = nm[j+5];
      counter += 6 ;
      cell_types[i] = (unsigned char)prsm ;
    }
    else if (nm[j+4] > -1) {
      conn[counter+4] = nm[j+4];
      counter += 5 ;
      cell_types[i] = (unsigned char)pyrm ;
    }
    else if (nm[j+3] > -1) {
      counter += 4 ;
      cell_types[i] = (unsigned char)tet ;
    }
    else printf("vtk extract doesn't understand cell connectivity\n") ;
    cell_offsets[i] = counter ;
  }
  Offset += counter * sizeof (int) + sizeof(int) ;
  fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  Offset += ncells * sizeof (int) + sizeof(int) ;
  fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  Offset += ncells * sizeof (unsigned char) + sizeof(int) ;
  fprintf(fid,"      </Cells>\n");
  fprintf(fid,"      <PointData>\n");
  fclose(fid);
}

void vtk_topo_handler::close() {
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"      </PointData>\n");
  fprintf(fid,"    </Piece>\n");
  fprintf(fid,"  </UnstructuredGrid>\n");
  fprintf(fid,"  <AppendedData encoding='raw'>\n");
  fprintf(fid,"_");
  int Scalar = npnts * sizeof (float), Vector = 3 * Scalar;
  int Cells = (ntets+nprsm+nhexs+npyrm) * sizeof(int);
  int CellChars = (ntets+nprsm+nhexs+npyrm) * sizeof(unsigned char);
  int Conn = (ntets*4+nprsm*6+nhexs*8+npyrm*5) * sizeof(int);
  fwrite((const char *) (&Vector), 4, 1, fid) ;
  fwrite((const char *) (&pos[0]), sizeof (float), 3 * npnts, fid) ;
  if (pos) delete [] pos ;
  fwrite((const char *) (&Conn), 4, 1, fid) ;
  fwrite((const char *) &conn[0], sizeof (int), 8*nhexs+6*nprsm+5*npyrm+4*ntets, fid) ;
  if (conn) delete [] conn ;
  fwrite((const char *) (&Cells), 4, 1, fid) ;
  fwrite((const char *) &cell_offsets[0], sizeof (int), nhexs+nprsm+npyrm+ntets, fid) ;
  if (cell_offsets) delete [] cell_offsets ;
  fwrite((const char *) (&CellChars), 4, 1, fid) ;
  fwrite((const char *) &cell_types[0], sizeof (unsigned char), nhexs+nprsm+npyrm+ntets, fid) ;
  if (cell_types) delete [] cell_types ;
  int curr_loc = 0;
  for (int i=0;i<(int)data_size.size();i++) {
    int ndata = data_size[i], ndata_size = ndata * sizeof (float); 
    fwrite((const char *) (&ndata_size), 4, 1, fid) ;
    fwrite((const char *) (&data_store[curr_loc]), sizeof (float), ndata, fid) ;
    curr_loc += ndata;
  }
  if (data_store) delete [] data_store;
  vector<int>().swap(data_size);
  fprintf(fid,"  </AppendedData>\n");
  fprintf(fid,"</VTKFile>\n");
  fclose(fid);
}
void vtk_topo_handler::create_mesh_positions(vector3d<float> position[], int pts) {
  pos = new float[3*npnts];
  for(int i=0;i<npnts;++i) {
    pos[3*i]   = position[i].x ;
    pos[3*i+1] = position[i].y ;
    pos[3*i+2] = position[i].z ;
  }
}

void vtk_topo_handler::write_tets(Array<int,4> tets[], int ntets, int block, int nblocks, int tottets) {
  for(int i=0;i<ntets;++i) {
    Array<int,8> brick ;
    brick[0] = tets[i][0]-1 ;
    brick[1] = tets[i][1]-1 ;
    brick[2] = tets[i][2]-1 ;
    brick[3] = tets[i][3]-1 ;
    brick[4] = -1 ;
    brick[5] = -1 ;
    brick[6] = -1 ;
    brick[7] = -1 ;
    bricks.push_back(brick) ;
  }    
}
void vtk_topo_handler::write_pyrm(Array<int,5> pyrm[], int npyrm, int block, int nblocks, int totpyrm) {
  for(int i=0;i<npyrm;++i) {
    Array<int,8> brick ;
    brick[0] = pyrm[i][0]-1 ;
    brick[1] = pyrm[i][1]-1 ;
    brick[2] = pyrm[i][2]-1 ;
    brick[3] = pyrm[i][3]-1 ;
    brick[4] = pyrm[i][4]-1 ;
    brick[5] = -1 ;
    brick[6] = -1 ;
    brick[7] = -1 ;
    bricks.push_back(brick) ;
  }
}
void vtk_topo_handler::write_prsm(Array<int,6> prsm[], int nprsm,int block, int nblocks, int totprsm) {
  for(int i=0;i<nprsm;++i) {
    Array<int,8> brick ;
    brick[0] = prsm[i][0]-1 ;
    brick[1] = prsm[i][1]-1 ;
    brick[2] = prsm[i][2]-1 ;
    brick[3] = prsm[i][3]-1 ;
    brick[4] = prsm[i][4]-1 ;
    brick[5] = prsm[i][5]-1 ;
    brick[6] = -1 ;
    brick[7] = -1 ;
    bricks.push_back(brick) ;
  }
}
void vtk_topo_handler::write_hexs(Array<int,8> hexs[], int nhexs, int block, int nblocks, int tothexs) {
  for(int i=0;i<nhexs;++i) {
    Array<int,8> brick ;
    for(int j=0;j<8;++j) brick[j] = hexs[i][j]-1 ;
    bricks.push_back(brick) ;
  }
}

void vtk_topo_handler::write_general_cell(int nfaces[], int nnfaces,
                                              int nsides[], int nnsides,
                                              int nodes[], int nnodes) {
  cerr << "vtk extract module doesn't support general cells!" << endl ;
}

void vtk_topo_handler::output_nodal_scalar(float val[], int npnts,
                                               string varname) {
  FILE * fid = fopen(filename.c_str(),"a") ;
  fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='1' format='appended' offset='%d'/>\n",varname.c_str(),Offset) ;
  Offset += npnts * sizeof(float) + sizeof(int) ;
  int ndata = npnts;
  data_size.push_back(ndata);
  for(int i = 0; i < npnts; i++) data_store[data_count++] = val[i] ;
  fclose(fid);
}

void vtk_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               int npnts, string varname) {
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='3' format='appended' offset='%d'/>\n",varname.c_str(),Offset);
  Offset += 3 * npnts * sizeof(float) + sizeof(int) ;
  int vec_size = 3 * npnts;
  data_size.push_back(vec_size);
  for(int i = 0; i < npnts; i++) { 
    data_store[data_count++] = val[i].x ;
    data_store[data_count++] = val[i].y ;
    data_store[data_count++] = val[i].z ;
  }
  fclose(fid);
}

void vtk_surf_topo_handler::open(string casename, string iteration ,int npnts,
                    int ntets, int nprsm, int npyrm, int nhexs, int ngen,
                    const vector<string> &bc_names,
                    const vector<string> &variables,
                    const vector<int> &variable_types,
                    double time)
{
  string filename = "vtk_surf_"+casename+"_"+iteration+".vtu" ;
  fid = fopen(filename.c_str(),"w");
  fprintf(fid,"<?xml version='1.0'?>\n");
  fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
  fprintf(fid,"  <UnstructuredGrid>\n");
}

void vtk_surf_topo_handler::close()
{
  npnts = node_ids.size();
  ncells = elem_types.size();
  int Offset = 0;
  fprintf(fid,"    <Piece NumberOfPoints='%d' NumberOfCells='%d'>\n",npnts,ncells);
  fprintf(fid,"      <Points>\n");
  fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%d'/>\n",Offset);
  fprintf(fid,"      </Points>\n");
  fprintf(fid,"      <Cells>\n");
  Offset += 3 * npnts * sizeof(float) + sizeof(int) ;
  fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  Offset += (int)elem_conn.size() * sizeof (int) + sizeof(int) ;
  fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  Offset += ncells * sizeof (int) + sizeof(int) ;
  fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%d'/>\n",Offset);
  Offset += ncells * sizeof (unsigned char) + sizeof(int) ;
  fprintf(fid,"      </Cells>\n");
  fprintf(fid,"      <CellData>\n");
  for (int i=0;i<(int)data_names.size();i++) { 
    int comp=-1;
    if      (data_size[i] ==     ncells) comp = 1;
    else if (data_size[i] == 3 * ncells) comp = 3;
    else { cout << "Wrong size" << endl; exit(1); }
    fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%d'/>\n",data_names[i].c_str(),comp,Offset) ;
    Offset += comp * ncells * sizeof (float) + sizeof(int); 
  }
  fprintf(fid,"      </CellData>\n");
  fprintf(fid,"    </Piece>\n");
  fprintf(fid,"  </UnstructuredGrid>\n");
  fprintf(fid,"  <AppendedData encoding='raw'>\n");
  fprintf(fid,"_");
  int Scalar = (int) node_ids.size() * sizeof (float), Vector = 3 * Scalar;
  int Cells = (int) elem_types.size() * sizeof(int);
  int CellChars = (int) elem_types.size() * sizeof(unsigned char);
  int Conn = (int) elem_conn.size() * sizeof(int);
  fwrite((const char *) (&Vector), 4, 1, fid) ;
  fwrite((const char *) (&position[0]), sizeof (float), 3 * (int)node_ids.size(), fid) ;
  fwrite((const char *) (&Conn), 4, 1, fid) ;
  fwrite((const char *) (&elem_conn[0]), sizeof (int), (int) elem_conn.size(), fid) ;
  fwrite((const char *) (&Cells), 4, 1, fid) ;
  fwrite((const char *) (&elem_offsets[0]), sizeof (int), (int) elem_offsets.size(), fid) ;
  fwrite((const char *) (&CellChars), 4, 1, fid) ;
  fwrite((const char *) (&elem_types[0]), sizeof (unsigned char), (int) elem_types.size(), fid) ;
  int curr_loc = 0;
  for (int i=0;i<(int)this->data_size.size();i++) {
    int ndata = this->data_size[i], ndata_size = ndata * sizeof (float);
    fwrite((const char *) (&ndata_size), 4, 1, fid) ;
    fwrite((const char *) (&elem_data[curr_loc]), sizeof (float), ndata, fid) ;
    curr_loc += ndata;
  }
  fprintf(fid,"  </AppendedData>\n");
  fprintf(fid,"</VTKFile>\n");
  fclose(fid);
}

void vtk_surf_topo_handler::create_boundary_part(string name,int node_set[], int np)
{
  output_boundary = false;
  for (int i=0;i<(int)boundaries.size();i++) {
    if (name == boundaries[i]) {
      output_boundary = true;
      for (int j=0;j<np;j++) {
	int node = node_set[j];
        node_ids.push_back(node);
	nmap[j+1+part_index] = node;
      }
    }
  }
}

void vtk_surf_topo_handler::create_mesh_positions(vector3d<float> Pos[],int np)
{
  position = new float[3*node_ids.size()];

  vector<int> G2L(np,-1);

  for (int i=0;i<(int)node_ids.size();i++) {
    int j=node_ids[i]-1;
    position[3*i]   = Pos[j].x; 
    position[3*i+1] = Pos[j].y; 
    position[3*i+2] = Pos[j].z; 
    G2L[j] = i;
  }
  for (int i=0;i<(int)elem_conn.size();i++) elem_conn[i] = G2L[elem_conn[i]];
  vector<int>().swap(G2L);
}

void vtk_surf_topo_handler::write_quads(Array<int,4> quads[], int quad_ids[],int nquads) 
{
  if (output_boundary) {
    int type = 9;
    unsigned char char_type = (unsigned char)type;
    for (int i=0;i<nquads;i++) {
      elem_ids.push_back(quad_ids[i]);
      for (int j=0;j<4;j++) {
	int id = nmap[quads[i][j]+part_index]-1;
        elem_conn.push_back(id);
      }
      elem_offset += 4;
      elem_offsets.push_back(elem_offset);
      elem_types.push_back(char_type);
    }
  }
}

void vtk_surf_topo_handler::write_trias(Array<int,3> trias[], int tria_ids[],int ntrias) 
{
  if (output_boundary) {
    int type = 5;
    unsigned char char_type = (unsigned char)type;
    for (int i=0;i<ntrias;i++) {
      elem_ids.push_back(tria_ids[i]);
      for (int j=0;j<3;j++) {
	int id = nmap[trias[i][j]+part_index]-1;
        elem_conn.push_back(id);
      }
      elem_offset += 3;
      elem_offsets.push_back(elem_offset);
      elem_types.push_back(char_type);
    }
  }
}

void vtk_surf_topo_handler::write_general_face(int nside_sizes[], int nside_ids[], int ngeneral,int nside_nodes[], int nside_nodes_size) 
{
  if (output_boundary) {
    int type = 7,cnt=0;
    unsigned char char_type = (unsigned char)type;
    for (int i=0;i<ngeneral;i++) {
      elem_ids.push_back(nside_ids[i]);
      for (int j=0;j<nside_sizes[i];j++) {
	int id = nmap[nside_nodes[cnt++]+part_index]-1;
        elem_conn.push_back(id);
      }
      elem_offset += nside_sizes[i];
      elem_offsets.push_back(elem_offset);
      elem_types.push_back(char_type);
    }
  }
}

void vtk_surf_topo_handler::close_boundary_part()
{
  if (output_boundary) part_index += node_ids.size();
  output_boundary = false;
}

void vtk_surf_topo_handler::output_boundary_scalar(float val[], int node_set[],int nvals, string valname) 
{
  data_names.push_back(valname);
  int size = elem_ids.size();
  data_size.push_back(size);
  
  vector<float> valout(elem_ids.size(),0);  
  
  if (!bmap.size()) for (size_t i=0;i<elem_ids.size();i++) bmap[elem_ids[i]] = i;

  for (int i=0;i<nvals;i++) valout[bmap[node_set[i]]] = val[i];

  for (int i=0;i<(int)elem_ids.size();i++) elem_data.push_back(valout[i]);
}

void vtk_surf_topo_handler::output_boundary_vector(vector3d<float> val[], int node_set[],int nvals, string valname) 
{
  data_names.push_back(valname);
  int size = 3 * elem_ids.size();
  data_size.push_back(size);
  
  vector<float> xvalout(elem_ids.size(),0);  
  vector<float> yvalout(elem_ids.size(),0);  
  vector<float> zvalout(elem_ids.size(),0);  
  
  if (!bmap.size()) for (size_t i=0;i<elem_ids.size();i++) bmap[elem_ids[i]] = i;

  for (int i=0;i<nvals;i++) {
    int j = bmap[node_set[i]];
    xvalout[j] = val[i].x;
    yvalout[j] = val[i].y;
    zvalout[j] = val[i].z;
  }

  for (int i=0;i<(int)elem_ids.size();i++) {
    elem_data.push_back(xvalout[i]);
    elem_data.push_back(yvalout[i]);
    elem_data.push_back(zvalout[i]);
  }
}
    

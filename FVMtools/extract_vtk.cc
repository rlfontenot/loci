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
  if (is_64bit) filename = "vtk64_"+casename+"_"+iteration+".vtu" ;
  else filename = "vtk_"+casename+"_"+iteration+".vtu" ;
  long long unsigned int ncells = ntets+nprsm+npyrm+nhexs+ngen;
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
  if (is_64bit) fprintf(fid,"<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>\n");
  else fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
  fprintf(fid,"  <UnstructuredGrid>\n");
  fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",npnts,ncells);
  fprintf(fid,"      <Points>\n");
  fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
  fprintf(fid,"      </Points>\n");
  Offset += 3 * npnts * sizeof(float) + int_size ;
  fclose(fid);
}
void vtk_topo_handler::close_mesh_elements() {
  unsigned int ncells = ntets+nprsm+npyrm+nhexs+ngen;
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"      <Cells>\n");
  fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += off * sizeof (int) + int_size ;
  fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += ncells * sizeof (int) + int_size ;
  fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += ncells * sizeof (unsigned char) + int_size ;
  if (ngen) { // include faces and faceoffsets
    fprintf(fid,"        <DataArray type='Int32' Name='faces' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += (long long unsigned int) cell_faces.size() * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='faceoffsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += (long long unsigned int) face_offsets.size() * sizeof (int) + int_size ;
  }
  fprintf(fid,"      </Cells>\n");
  fprintf(fid,"      <PointData>\n");
  fclose(fid);
}

void vtk_topo_handler::close() {
  /*
  cout << "Int:                    " << std::numeric_limits<int>::max() << endl;
  cout << "Unsigned Int:           " << std::numeric_limits<unsigned int>::max() << endl;
  cout << "Long Long Unsigned Int: " << std::numeric_limits<long long unsigned int>::max() << endl;
  */
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"      </PointData>\n");
  fprintf(fid,"    </Piece>\n");
  fprintf(fid,"  </UnstructuredGrid>\n");
  fprintf(fid,"  <AppendedData encoding='raw'>\n");
  fprintf(fid,"_");
  long long unsigned int Scalar = npnts * sizeof (float), Vector = 3 * Scalar;
  long long unsigned int Cells = (long long unsigned int) cell_types.size() * sizeof(int);
  long long unsigned int CellChars = (long long unsigned int) cell_types.size() * sizeof(unsigned char);
  long long unsigned int Conn = (long long unsigned int) conn.size() * sizeof(int);
  if (!is_64bit) {
    if ( Scalar > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
    if ( Vector > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
    if ( Cells > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
    if ( Conn > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
  }
  fwrite((const char *) (&Vector), int_size, 1, fid) ;
  fwrite((const char *) (&pos[0]), sizeof (float), 3 * npnts, fid) ;
  if (pos) delete [] pos ;
  fwrite((const char *) (&Conn), int_size, 1, fid) ;
  fwrite((const char *) &conn[0], sizeof (int), (long long unsigned int) conn.size(), fid) ;
  vector<int>().swap(conn) ;
  fwrite((const char *) (&Cells), int_size, 1, fid) ;
  fwrite((const char *) &cell_offsets[0], sizeof (int), (long long unsigned int) cell_offsets.size(), fid) ;
  vector<int>().swap(cell_offsets) ;
  fwrite((const char *) (&CellChars), int_size, 1, fid) ;
  fwrite((const char *) &cell_types[0], sizeof (unsigned char), (long long unsigned int) cell_types.size(), fid) ;
  vector<unsigned char>().swap(cell_types) ;
  if (ngen) {
    long long unsigned int Faces = (long long unsigned int) cell_faces.size() * sizeof(int) ;
    long long unsigned int FaceConn = (long long unsigned int) face_offsets.size() * sizeof(int) ;
    if (!is_64bit) {
      if ( Faces > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
      if ( FaceConn > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
    }
    fwrite((const char *) (&Faces), int_size, 1, fid) ;
    fwrite((const char *) &cell_faces[0], sizeof (int), (long long unsigned int) cell_faces.size(), fid) ;
    fwrite((const char *) (&FaceConn), int_size, 1, fid) ;
    fwrite((const char *) &face_offsets[0], sizeof (int), (long long unsigned int) face_offsets.size(), fid) ;
  }
  long long unsigned int curr_loc = 0;
  for (long long unsigned int i=0;i<(long long unsigned int)data_size.size();i++) {
    long long unsigned int ndata = data_size[i], ndata_size = ndata * sizeof (float); 
    if (!is_64bit) {
      if ( ndata_size > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
    }
    fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
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
  for(long long unsigned int i=0;i<npnts;++i) {
    long long unsigned int j=3*i;
    pos[j]   = position[i].x ;
    pos[j+1] = position[i].y ;
    pos[j+2] = position[i].z ;
  }
}

void vtk_topo_handler::write_tets(Array<int,4> cells[], int ncells, int block, int nblocks, int tottets) {
  unsigned int type = 10 ;
  unsigned char ctype = (unsigned char)type ;
  int nnpc = 4;

  const int nf=4,nnpf[nf]={3,3,3,3};
  vector< vector<int> > ff(nf);
  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
  ff[0][0] = 0; ff[0][1] = 1; ff[0][2] = 3;
  ff[1][0] = 1; ff[1][1] = 2; ff[1][2] = 3;
  ff[2][0] = 2; ff[2][1] = 0; ff[2][2] = 3;
  ff[3][0] = 0; ff[3][1] = 2; ff[3][2] = 1;

  if (ngen) {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
      
      cell_faces.push_back(nf); face_off++;
      for (int j=0;j<nf;j++) {
        cell_faces.push_back(nnpf[j]); face_off++;
        for (int k=0;k<nnpf[j];k++) {
	  int nn = cells[i][ff[j][k]] - 1;
	  cell_faces.push_back(nn); face_off++;
	}
      }
      face_offsets.push_back(face_off) ;
    }
  }
  else {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
    }
  }
}

void vtk_topo_handler::write_pyrm(Array<int,5> cells[], int ncells, int block, int nblocks, int totpyrm) {
  unsigned int type = 14 ;
  unsigned char ctype = (unsigned char)type ;
  int nnpc = 5;

  const int nf=5,nnpf[nf]={4,3,3,3,3};
  vector< vector<int> > ff(nf);
  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
  ff[0][0] = 0; ff[0][1] = 3; ff[0][2] = 2; ff[0][3] = 1;
  ff[1][0] = 0; ff[1][1] = 1; ff[1][2] = 4;
  ff[2][0] = 1; ff[2][1] = 2; ff[2][2] = 4;
  ff[3][0] = 2; ff[3][1] = 3; ff[3][2] = 4;
  ff[4][0] = 3; ff[4][1] = 0; ff[4][2] = 4; 

  if (ngen) {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
      
      cell_faces.push_back(nf); face_off++;
      for (int j=0;j<nf;j++) {
        cell_faces.push_back(nnpf[j]); face_off++;
        for (int k=0;k<nnpf[j];k++) {
	  int nn = cells[i][ff[j][k]] - 1;
	  cell_faces.push_back(nn); face_off++;
	}
      }
      face_offsets.push_back(face_off) ;
    }
  }
  else {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
    }
  }
}

void vtk_topo_handler::write_prsm(Array<int,6> cells[], int ncells, int block, int nblocks, int totprsm) {
  unsigned int type = 13 ;
  unsigned char ctype = (unsigned char)type ;
  int nnpc = 6;

  const int nf=5,nnpf[nf]={4,4,4,3,3};
  vector< vector<int> > ff(nf);
  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
  ff[0][0] = 0; ff[0][1] = 2; ff[0][2] = 5; ff[0][3] = 3;
  ff[1][0] = 2; ff[1][1] = 1; ff[1][2] = 4; ff[1][3] = 5;
  ff[2][0] = 0; ff[2][1] = 3; ff[2][2] = 4; ff[2][3] = 1;
  ff[3][0] = 0; ff[3][1] = 1; ff[3][2] = 2;
  ff[4][0] = 3; ff[4][1] = 5; ff[4][2] = 4; 

  if (ngen) {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
      
      cell_faces.push_back(nf); face_off++;
      for (int j=0;j<nf;j++) {
        cell_faces.push_back(nnpf[j]); face_off++;
        for (int k=0;k<nnpf[j];k++) {
	  int nn = cells[i][ff[j][k]] - 1;
	  cell_faces.push_back(nn); face_off++;
	}
      }
      face_offsets.push_back(face_off) ;
    }
  }
  else {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
    }
  }
}

void vtk_topo_handler::write_hexs(Array<int,8> cells[], int ncells, int block, int nblocks, int tothexs) {
  unsigned int type = 12 ;
  unsigned char ctype = (unsigned char)type ;
  int nnpc = 8;

  const int nf=6,nnpf[nf]={4,4,4,4,4,4};
  vector< vector<int> > ff(nf);
  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
  ff[0][0] = 0; ff[0][1] = 1; ff[0][2] = 5; ff[0][3] = 4;
  ff[1][0] = 1; ff[1][1] = 2; ff[1][2] = 6; ff[1][3] = 5;
  ff[2][0] = 2; ff[2][1] = 3; ff[2][2] = 7; ff[2][3] = 6;
  ff[3][0] = 3; ff[3][1] = 0; ff[3][2] = 4; ff[3][3] = 7;
  ff[4][0] = 4; ff[4][1] = 5; ff[4][2] = 6; ff[4][3] = 7;
  ff[5][0] = 0; ff[5][1] = 3; ff[5][2] = 2; ff[5][3] = 1;

  if (ngen) {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
      
      cell_faces.push_back(nf); face_off++;
      for (int j=0;j<nf;j++) {
        cell_faces.push_back(nnpf[j]); face_off++;
        for (int k=0;k<nnpf[j];k++) {
	  int nn = cells[i][ff[j][k]] - 1;
	  cell_faces.push_back(nn); face_off++;
	}
      }
      face_offsets.push_back(face_off) ;
    }
  }
  else {
    for(int i=0;i<ncells;++i) {
      for (int j=0;j<nnpc;j++) { 
        int nd = cells[i][j]-1 ;
        conn.push_back(nd) ;
        off++ ;
      }
      cell_types.push_back(ctype);
      cell_offsets.push_back(off);
    }
  }
}

void vtk_topo_handler::write_general_cell(int nfaces[], int nnfaces,
                                              int nsides[], int nnsides,
                                              int nodes[], int nnodes) {
  unsigned int type = 42 ;
  unsigned char ctype = (unsigned char)type ;
  int face = 0, node = 0 ;
  vector<int> tmp;
  
  for(int i=0;i<nnfaces;++i) {
    int tsz=0,nf = nfaces[i] ;
    for (int j=0;j<nf;j++) tsz += nsides[face++] ; 
    tmp.resize(tsz) ; face-=nf ; tsz=0 ;
    cell_faces.push_back(nf) ; face_off++ ;
    for (int j=0;j<nf;j++) { 
      int ne = nsides[face++] ; 
      cell_faces.push_back(ne) ; face_off++ ;
      for (int k=0;k<ne;k++) { 
        int nd = nodes[node++]-1 ;
	cell_faces.push_back(nd) ; face_off++ ;
	tmp[tsz++] = nd ;
      }
    }

    sort(tmp.begin(),tmp.end()) ;
    tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end()) ;
    
    for (size_t j=0;j<tmp.size();j++) conn.push_back(tmp[j]) ;  
    off += (long long unsigned int) tmp.size() ;
    
    cell_types.push_back(ctype) ;
    cell_offsets.push_back(off) ;
    face_offsets.push_back(face_off) ;
  }

  if (nnodes  != node) cerr << "Problem in write_general_cell" << endl ;
  if (nnsides != face) cerr << "Problem in write_general_cell" << endl ;
}

void vtk_topo_handler::output_nodal_scalar(float val[], int npnts,
                                               string varname) {
  FILE * fid = fopen(filename.c_str(),"a") ;
  fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='1' format='appended' offset='%llu'/>\n",varname.c_str(),Offset) ;
  Offset += (long long unsigned int) npnts * sizeof(float) + int_size ;
  long long unsigned int ndata = (long long unsigned int) npnts;
  data_size.push_back(ndata);
  for(long long unsigned int i = 0; i < (long long unsigned int)npnts; i++) data_store[data_count++] = val[i] ;
  fclose(fid);
}

void vtk_topo_handler::output_nodal_vector(vector3d<float> val[],
                                               int npnts, string varname) {
  FILE * fid = fopen(filename.c_str(),"a");
  fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='3' format='appended' offset='%llu'/>\n",varname.c_str(),Offset);
  Offset += 3 * (long long unsigned) npnts * sizeof(float) + int_size ;
  long long unsigned int vec_size = 3 * (long long unsigned int)npnts;
  data_size.push_back(vec_size);
  for(long long unsigned int i = 0; i < (long long unsigned int) npnts; i++) { 
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
  string filename;
  if (is_64bit) filename = "vtk_surf64_"+casename+"_"+iteration+".vtu" ;
  else filename = "vtk_surf_"+casename+"_"+iteration+".vtu" ;
  fid = fopen(filename.c_str(),"w");
  fprintf(fid,"<?xml version='1.0'?>\n");
  if (is_64bit) fprintf(fid,"<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>\n");
  else fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
  fprintf(fid,"  <UnstructuredGrid>\n");
}

void vtk_surf_topo_handler::close()
{
  npnts = node_ids.size();
  ncells = elem_types.size();
  long long unsigned int Offset = 0;
  fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",npnts,ncells);
  fprintf(fid,"      <Points>\n");
  fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
  fprintf(fid,"      </Points>\n");
  fprintf(fid,"      <Cells>\n");
  Offset += 3 * npnts * sizeof(float) + int_size ;
  fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += (long long unsigned int)elem_conn.size() * sizeof (int) + int_size ;
  fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += ncells * sizeof (int) + int_size ;
  fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
  Offset += ncells * sizeof (unsigned char) + int_size ;
  fprintf(fid,"      </Cells>\n");
  fprintf(fid,"      <CellData>\n");
  for (long long unsigned int i=0;i<(long long unsigned int)data_names.size();i++) { 
    int comp=-1;
    if      ( (int) data_size[i] ==     (int) ncells) comp = 1;
    else if ( (int) data_size[i] == 3 * (int) ncells) comp = 3;
    else { cout << "Wrong size" << endl; exit(1); }
    fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",data_names[i].c_str(),comp,Offset) ;
    Offset += comp * ncells * sizeof (float) + int_size; 
  }
  fprintf(fid,"      </CellData>\n");
  fprintf(fid,"    </Piece>\n");
  fprintf(fid,"  </UnstructuredGrid>\n");
  fprintf(fid,"  <AppendedData encoding='raw'>\n");
  fprintf(fid,"_");
  long long unsigned int Scalar = (long long unsigned int) node_ids.size() * sizeof (float), Vector = 3 * Scalar;
  long long unsigned int Cells = (long long unsigned int) elem_types.size() * sizeof(int);
  long long unsigned int CellChars = (long long unsigned int) elem_types.size() * sizeof(unsigned char);
  long long unsigned int Conn = (long long unsigned int) elem_conn.size() * sizeof(int);
  fwrite((const char *) (&Vector), int_size, 1, fid) ;
  fwrite((const char *) (&position[0]), sizeof (float), 3 * (long long unsigned int)node_ids.size(), fid) ;
  fwrite((const char *) (&Conn), int_size, 1, fid) ;
  fwrite((const char *) (&elem_conn[0]), sizeof (int), (long long unsigned int) elem_conn.size(), fid) ;
  fwrite((const char *) (&Cells), int_size, 1, fid) ;
  fwrite((const char *) (&elem_offsets[0]), sizeof (int), (long long unsigned int) elem_offsets.size(), fid) ;
  fwrite((const char *) (&CellChars), int_size, 1, fid) ;
  fwrite((const char *) (&elem_types[0]), sizeof (unsigned char), (long long unsigned int) elem_types.size(), fid) ;
  long long unsigned int curr_loc = 0;
  for (long long unsigned int i=0;i<(long long unsigned int)this->data_size.size();i++) {
    long long unsigned int ndata = this->data_size[i], ndata_size = ndata * sizeof (float);
    fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
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
    

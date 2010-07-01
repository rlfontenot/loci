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
#include <iostream>
#include <fstream>
#include "Tools/xdr.h"
#include <algorithm>
#include <map> 
#include <vector> 
using std::cout ;
using std::endl ;
using std::cerr ;
using std::map ;
int main(int ac, char* av[]) {
  char* filename = av[1] ; 
  cout << " Dumping into new binary grid format " << endl ;
  char buf[512] ;
  sprintf(buf,"%s.cog",filename) ;
  FILE* IFP = fopen(buf, "r") ;
  if(IFP == NULL) {
    cerr << "can't open '" << buf << "'" << endl ;
    exit(-1) ;
  }
  const int MAX = 50000 ;
  int npatch, maxppf, maxfpc ;
  int ndim, nzones ;
  int npnts, nfaces, ncells ;
  if(fscanf(IFP, "%d%d%d", &ndim, &nzones, &npatch)!=3) {
    cerr << "fscanf failed" << endl ;
    exit(-1) ;
  }
  if(ndim != 3) {
    cerr << "currently only supports 3-D grid files" << endl ;
    exit(-1) ;
  }
  if(nzones != 1) {
    cerr << "currently only supports single zone files" << endl ;
    exit(-1) ;
  }
  if(fscanf(IFP, "%d%d%d%d%d", &npnts, &nfaces, &ncells, &maxppf, &maxfpc)!=5){
    cerr << "scanf failed" << endl ;
    exit(-1) ;
  }
  cout << " ndim = " << ndim << endl ;
  cout << " npnts = " << npnts << endl ;
  cout << " nfaces = " << nfaces << endl ;
  cout << "ncells = " << ncells << endl ;
  
  char out_buf[512] ;
  sprintf(out_buf,"%s.xdr",filename) ;
  FILE *FP = fopen(out_buf, "w") ;
  if(FP == NULL) {
    cerr << "can't open " << out_buf <<  endl ;
    return(-1);
  }
  XDR xdr_handle ;
  xdrstdio_create(&xdr_handle, FP, XDR_ENCODE) ;
  xdr_int(&xdr_handle, &ndim) ;
  xdr_int(&xdr_handle, &nzones) ;
  xdr_int(&xdr_handle, &npatch) ;
  xdr_int(&xdr_handle, &npnts) ;
  xdr_int(&xdr_handle, &nfaces) ;
  xdr_int(&xdr_handle, &ncells) ;
  xdr_int(&xdr_handle, &maxppf) ;
  xdr_int(&xdr_handle, &maxfpc) ;
  int node_reads, face_reads ;
  int rem_nodes, rem_faces ;
  node_reads = npnts / MAX ;
  rem_nodes = npnts % MAX ;
  face_reads = nfaces / MAX ;
  rem_faces = nfaces % MAX ;
  double* tmp_pos ;
  if(node_reads)
    tmp_pos = new double[3*MAX] ;
  else
    tmp_pos = new double[3*rem_nodes] ;
  for(int i = 0; i < node_reads; ++i) {
    int tmp = 0 ;
    for(int j = 0; j < 3*MAX; ++j) {
      if(fscanf(IFP, "%lf", &tmp_pos[tmp])!=1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      xdr_double(&xdr_handle, &tmp_pos[tmp++]) ;
    }
  }
  int tmp = 0 ;
  for(int j = 0; j < 3*rem_nodes; ++j) {
    if(fscanf(IFP, "%lf", &tmp_pos[tmp]) != 1) {
      cerr << "scanf failed" ;
      exit(-1) ;
    }
    xdr_double(&xdr_handle, &tmp_pos[tmp++]) ;
  }
  
  fpos_t after_pos ;
  fgetpos(IFP, &after_pos) ;
  int local_size ;
  int dummy_node ; 
  std::vector<int> offset(3) ;
  int off = 0 ;
  map<int, std::vector<std::vector<int> > > hm_int ;
  
  for(int i = 0; i < face_reads; ++i) {
    for(int j = 0; j < MAX; ++j) {
      if(fscanf(IFP, "%d", &local_size)!= 1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      for(int k = 0; k < local_size; ++k)
	if(fscanf(IFP, "%d", &dummy_node) !=1) {
	  cerr << "scanf failed" << endl ; 
	  exit(-1) ;
	}
      if(fscanf(IFP, "%d%d", &offset[1], &offset[2]) !=2 ) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      if((offset[1] > 0) && (offset[2] > 0)) {
	offset[0] = off ;
	xdr_int(&xdr_handle, &offset[0]) ;
	xdr_int(&xdr_handle, &offset[1]) ;
	xdr_int(&xdr_handle, &offset[2]) ;
	off += local_size ;
      }
      else {
	int loc_min = std::min(offset[1], offset[2]) ;
	offset[0] = local_size ;
	(hm_int[loc_min]).push_back(offset) ;
      }
    }
  }
  for(int j = 0; j < rem_faces; ++j) {
    if(fscanf(IFP, "%d", &local_size)!= 1) {
      cerr << "scanf failed" << endl ;
      exit(-1);
    }
    for(int k = 0; k < local_size; ++k)
      if(fscanf(IFP, "%d", &dummy_node)!=1) {
	cerr << "scanf failed" << endl ; 
	exit(-1) ;
      }
    if(fscanf(IFP, "%d%d", &offset[1], &offset[2]) != 2) {
      cerr << "scanf failed" << endl ;
      exit(-1); 
    }
    if((offset[1] > 0) && (offset[2] > 0)) {
      offset[0] = off ;
      xdr_int(&xdr_handle, &offset[0]) ;
      xdr_int(&xdr_handle, &offset[1]) ;
      xdr_int(&xdr_handle, &offset[2]) ;
      off += local_size ;
    }
    else {
      int loc_min = std::min(offset[1], offset[2]) ;
      offset[0] = local_size ;
      hm_int[loc_min].push_back(offset) ; 
    }
  }
  std::map<int, std::vector<std::vector<int> > >::const_iterator hmi ;
  for(hmi = hm_int.begin(); hmi != hm_int.end(); ++hmi) {
    for(std::vector<std::vector<int> >::const_iterator vvi = hmi->second.begin(); vvi != hmi->second.end(); ++vvi) {
      std::vector<int> tmp_vec = *vvi ; 
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &tmp_vec[1]) ;
      xdr_int(&xdr_handle, &tmp_vec[2]) ;
      off += tmp_vec[0] ;
    }
  }
  xdr_int(&xdr_handle, &off) ;
  int* node_vec = new int[maxppf] ;
  int dnode[2] ;
  std::map<int, std::vector<std::vector<int> > >mv_int ;
  const fpos_t tmp_position = after_pos ;
  fsetpos(IFP, &tmp_position) ;
  for(int i = 0; i < face_reads; ++i) {
    for(int j = 0; j < MAX; ++j) {
      if(fscanf(IFP, "%d", &local_size)!=1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      int tmp = 0 ;
      for(int k = 0; k < local_size; ++k) {
	if(fscanf(IFP, "%d", &node_vec[tmp])!=1) {
	  cerr << "scanf failed" << endl ;
	  exit(-1) ;
	}
	node_vec[tmp] -= 1 ;
	tmp++ ;
	//xdr_int(&xdr_handle, &node_vec[tmp++]) ;
      }
      if(fscanf(IFP, "%d", &dnode[0])!=1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      if(fscanf(IFP, "%d", &dnode[1])!=1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      if((dnode[0] > 0) && (dnode[1] > 0)) {
	for(int k = 0; k < local_size; ++k)
	  xdr_int(&xdr_handle, &node_vec[k]) ;
      }
      else {
	int loc_min = std::min(dnode[0], dnode[1]) ;
	std::vector<int> tmp_vec ;
	for(int k = 0; k < local_size; ++k)
	  tmp_vec.push_back(node_vec[k]) ;
	mv_int[loc_min].push_back(tmp_vec) ;
      }
    }
  }
  for(int j = 0; j < rem_faces; ++j) {
    if(fscanf(IFP, "%d", &local_size)!=1) {
      cerr << "scanf failed" << endl ;
      exit(-1) ;
    }
    int tmp = 0 ;
    for(int k = 0; k < local_size; ++k) {
      if(fscanf(IFP, "%d", &node_vec[tmp])!=1) {
	cerr << "scanf failed" << endl ;
	exit(-1) ;
      }
      node_vec[tmp] -= 1 ;
      tmp++ ;
      //xdr_int(&xdr_handle, &node_vec[tmp++]) ;
    } 
    if(fscanf(IFP, "%d", &dnode[0]) !=1) {
      cerr << "scanf failed" << endl ;
      exit(-1) ;
    }
    if(fscanf(IFP, "%d", &dnode[1]) != 1) {
      cerr << "scanf failed" << endl ;
      exit(-1) ;
    }

    if((dnode[0] > 0) && (dnode[1] > 0)) {
      for(int k = 0; k < local_size; ++k)
	xdr_int(&xdr_handle, &node_vec[k]) ;
    }
    else {
      int loc_min = std::min(dnode[0], dnode[1]) ;
      std::vector<int> tmp_vec ;
      for(int k = 0; k < local_size; ++k)
	tmp_vec.push_back(node_vec[k]) ;
      mv_int[loc_min].push_back(tmp_vec) ;
    }
  }
  std::map<int, std::vector<std::vector<int> > >::const_iterator mvi ;
  for(mvi = mv_int.begin(); mvi != mv_int.end(); ++mvi) {
    for(std::vector<std::vector<int> >::const_iterator vvi = mvi->second.begin(); vvi != mvi->second.end(); ++vvi) {
      std::vector<int> tmp_vec = *vvi ; 
      for(size_t k = 0; k < tmp_vec.size(); ++k)
	xdr_int(&xdr_handle, &tmp_vec[k]) ;
    }
  }
  fclose(FP);
  fclose(IFP) ;
  xdr_destroy(&xdr_handle) ;
  delete [] tmp_pos ;
  return(0) ;
  
}
  

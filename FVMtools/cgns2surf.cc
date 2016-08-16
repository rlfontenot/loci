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

/*cgns library can be compiled with 64-bit or 32-bit,
  which will decide cgsize_t is long int or int
  Loci library can also be compiled with 64-bit or 32-bit, which will affect gEntity
  static_cast is used in case the size of gEntity and cgsize_t are different
*/

/*
  Only read in Grid specification part, i.e., grid coordinates and element connectivity.   
  1. All data located in one file. The file name is provided by user 
  2. Only read in the first zone of first base
  4. Only unstructured grid. Vertex dimension is 3D, and physical dimension is 3D.
  5. The data read only contains core vertex coordinates and core elements, and boundary conditions
  6. Allowed element type: TRI_3, QUAD_4,
  TETRA_4,  PYRA_5, PYRA_14,
  PENTA_6,  HEXA_8,
*/

/*write vog
  1. read in boundary surface and boundary names
  2. remove transparent faces
  3. post scale unit
  
*/
/*questions:
  1. Does vog file deal with structured grid?
  2. Is boundary condition info need to written into another function?
  3. Does cgns know anything about transparent faces?
*/


#include <strings.h>
#include <string.h>
#include <stdio.h>
//#include "vogtools.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include <ctype.h>
#include <iostream>
#ifdef _WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif
#include "cgnslib.h"
#include "hash.h"

#ifndef CGNS_ENUMV
# define CGNS_ENUMV(V) V
# define CGNS_ENUMT(T) T
#endif

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;
using std::string ;
using std::vector ;
using std::istringstream ;
using std::pair;
using std::ofstream;
using std::min;
using std::max;


typedef struct {
  cgsize_t id;
  int bcnum;
  int nnodes;
  cgsize_t nodes[4];
} FACE;



/*-------------------------------------------------------------------*/

static int compare_faces (void *v1, void *v2)
{
  FACE *f1 = (FACE *)v1;
  FACE *f2 = (FACE *)v2;
  int i, k;
  cgsize_t id, nn, n1[4], n2[4];

  if (f1->nnodes != f2->nnodes)
    return (f1->nnodes - f2->nnodes);

  for (i = 0; i < f1->nnodes; i++) {
    id = f1->nodes[i];
    for (k = 0; k < i; k++) {
      if (n1[k] > id) {
        nn = n1[k];
        n1[k] = id;
        id = nn;
      }
    }
    n1[i] = id;
  }
  for (i = 0; i < f2->nnodes; i++) {
    id = f2->nodes[i];
    for (k = 0; k < i; k++) {
      if (n2[k] > id) {
        nn = n2[k];
        n2[k] = id;
        id = nn;
      }
    }
    n2[i] = id;
  }

  for (i = 0; i < f1->nnodes; i++) {
    if (n1[i] != n2[i])
      return (int)(n1[i] - n2[i]);
  }
  return 0;
}
/*-------------------------------------------------------------------*/
void ListHash(HASH hash, FACE ** Faces){
  HASH_TAB *tabp = (HASH_TAB *)hash;
  int ne = 0;
  for (size_t i = 0; i < tabp->size; i++) {
    for (BUCKET *p = tabp->table[i]; p != NULL; p = p->next){
      Faces[ne++] = (FACE *)(p->entry);
    }
  }
}

/*-------------------------------------------------------------------*/

static size_t hash_face (void *v)
{
  FACE *f = (FACE *)v;
  int n;
  size_t hash = 0;

  for (n = 0; n < f->nnodes; n++)
    hash += (size_t)f->nodes[n];
  return hash;
}

/*-------------------------------------------------------------------*/
static void error_exit(const char *func)
{
  printf("CGNSlib ERROR:");
  if (func != NULL && *func)
    printf("%s:", func);
  printf("%s\n", cg_get_error());
  exit(1);
}
/*-------------------------------------------------------------------*/

static FACE *new_face (cgsize_t faceid)
{
  FACE *face = (FACE *)malloc(sizeof(FACE));
  if (face == NULL)
    error_exit( "malloc failed for a new face");
  face->id = faceid;
  face->bcnum = 0;
  face->nnodes = 0;
  return face;
}

/*--------------------------------------------------------------------*/

static int sort_faces (const void *v1, const void *v2)
{
  FACE **f1 = (FACE **)v1;
  FACE **f2 = (FACE **)v2;

  return (int)((*f1)->id - (*f2)->id);
}
  
/*--------------------------------------------------------------------*/

static int sort_points (const void *v1, const void *v2)
{
  return (int)(*((cgsize_t *)v1) - *((cgsize_t *)v2));
}

/*--------------------------------------------------------------------*/

static int find_point (cgsize_t id, cgsize_t np, cgsize_t *pts)
{
  cgsize_t lo = 0, hi = np - 1, mid;

  if (!np || id < pts[0]) return 0;
  if (id == pts[0]) return 1;
  if (!hi || id > pts[hi]) return 0;
  if (id == pts[hi]) return 1;

  while (lo <= hi) {
    mid = (lo + hi) >> 1;
    if (id == pts[mid]) return 1;
    if (id < pts[mid])
      hi = mid - 1;
    else
      lo = mid + 1;
  }
  return 0;
}





struct surf_face
{
  int bc_id;
  vector<cgsize_t> trias;
  vector<cgsize_t> quads;
};

/*--------------------------------------------------------------------*/


static surf_face* get_bc(vector<surf_face>& bc_data, int id){
  for(size_t i = 0; i < bc_data.size(); i++){
    if(bc_data[i].bc_id == id) return &bc_data[i];
  }
  return NULL;
}
  


/*--------------------------------------------------------------------*/
/* set the bcnum of faces */ 
static void unstructured_boundary (int nb, cgsize_t nFaces, FACE** Faces,
                                   CGNS_ENUMT(PointSetType_t) ptype,
                                   CGNS_ENUMT(GridLocation_t) location,
                                   cgsize_t np, cgsize_t *ptset)
{
  cgsize_t nf, is, ie;
  int n;

  if (ptype == CGNS_ENUMV(PointRange) ||
      ptype == CGNS_ENUMV(PointList)) {
    if (location == CGNS_ENUMV(FaceCenter) ||
        location == CGNS_ENUMV(CellCenter)) {
      ptype = (ptype == CGNS_ENUMV(PointRange) ?
               CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
    }
    else if (location != CGNS_ENUMV(Vertex)) {
      return;
    }
  }

  if (ptype == CGNS_ENUMV(PointRange)) {
    if (ptset[0] < ptset[1]) {
      is = ptset[0];
      ie = ptset[1];
    }
    else {
      is = ptset[1];
      ie = ptset[0];
    }
    for (nf = 0; nf < nFaces; nf++) {
      for (n = 0; n < Faces[nf]->nnodes; n++) {
        if (Faces[nf]->nodes[n] < is ||
            Faces[nf]->nodes[n] > ie) break;
      }
      if (n == Faces[nf]->nnodes)
        Faces[nf]->bcnum = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(PointList)) {
    qsort(ptset, (size_t)np, sizeof(cgsize_t), sort_points);
    for (nf = 0; nf < nFaces; nf++) {
      for (n = 0; n < Faces[nf]->nnodes; n++) {
        if (!find_point(Faces[nf]->nodes[n], np, ptset))
          break;
      }
      if (n == Faces[nf]->nnodes)
        Faces[nf]->bcnum = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(ElementRange)) {
    if (ptset[0] < ptset[1]) {
      is = ptset[0];
      ie = ptset[1];
    }
    else {
      is = ptset[1];
      ie = ptset[0];
    }
    for (nf = 0; nf < nFaces; nf++) {
      if (Faces[nf]->id >= is && Faces[nf]->id <= ie)
        Faces[nf]->bcnum = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(ElementList)) {
    qsort(ptset, (size_t)np, sizeof(cgsize_t), sort_points);
    for (nf = 0; nf < nFaces; nf++) {
      if (find_point(Faces[nf]->id, np, ptset))
        Faces[nf]->bcnum = nb;
    }
    return;
  }
}
/* read in boundary conditions, set the bcbum for each face, and store the boundary names in bc_ids */
static void boundary_conditions (cgsize_t nFaces,  FACE** Faces, cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                 vector<pair<int,string> >& bc_ids)
{
  int nb, ib, nrmlindex[3];
  cgsize_t is;
  cgsize_t np; //num_points
  cgsize_t *ptset;//point  set
  char name[33];
  CGNS_ENUMT(BCType_t) bctype;
  CGNS_ENUMT(PointSetType_t) ptype;
  CGNS_ENUMT(GridLocation_t) location;
  CGNS_ENUMT(DataType_t) datatype;
  int dim =  1; //for unstructured mesh
  int nBocos = 0;
  if (cg_nbocos (cgFile, cgBase, cgZone, &nBocos))
    error_exit ("cg_nbocos");

  if (nBocos) {
    // Bocos = (BOCO *)malloc(nBocos * sizeof(BOCO));
    // if (Bocos == NULL) error_exit( "malloc failed for bocos");
 
    for (nb = 1; nb <= nBocos; nb++) {
      if (cg_boco_info(cgFile, cgBase, cgZone, nb, name,
                       &bctype, &ptype, &np, nrmlindex, &is,
                       &datatype, &ib))
        error_exit("cg_boco_info");
      //Bocos[nb-1].type = bctype;
      //strcpy(Bocos[nb-1].name, name);
      string bc_name = string(name);
      bc_ids.push_back(pair<int, string>(nb, bc_name));
      if (cg_boco_gridlocation_read(cgFile, cgBase, cgZone,
                                    nb, &location))
        error_exit("cg_boco_gridlocation_read");
      ptset = (cgsize_t *)malloc((size_t)(np * dim) * sizeof(cgsize_t));
      if (ptset == NULL)
        error_exit( "malloc failed for boco ptset");
      if (cg_boco_read(cgFile, cgBase, cgZone, nb, ptset, 0))
        error_exit("cg_boco_read");

      //set up the bcnum of faces
      unstructured_boundary(nb, nFaces, Faces, ptype, location, np, ptset);
      
      free(ptset);
    }
  }
  //why?
  for (np = 0; np < nFaces; np++) {
    if (Faces[np]->bcnum < 0){
      cout<< " negative bcnum" << Faces[np]->bcnum<<endl;
      Faces[np]->bcnum = nBocos - Faces[np]->bcnum;
    }
  }
}



/*read in and hash the faces, set up boundary number and boundary names*/ 
static void unstructured_faces (HASH& facehash, cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                const vector<string>& boundaries,vector<surf_face>& bc_data
                                ){

  
  
  vector<pair<int,string> > surf_ids;
  int nsect = 0, nn = 0, ip = 0;
  
  CGNS_ENUMT(ElementType_t) elemtype, et;
  char name[33];
  cgsize_t i, n, is, ie, ne;
  cgsize_t size, *conn;
  FACE face, *pf;
  
  if (cg_nsections (cgFile, cgBase, cgZone, &nsect))
    error_exit ("cg_nsections");
  for (int ns = 1; ns <= nsect; ns++) {
    int nf = 0;
    if (cg_section_read (cgFile, cgBase, cgZone, ns,
                         name, &elemtype, &is, &ie, &nn, &ip))
      error_exit ("cg_section_read");
    if (elemtype < CGNS_ENUMV(TRI_3) ||
        (elemtype > CGNS_ENUMV(QUAD_9) &&
         elemtype != CGNS_ENUMV(MIXED))) continue;
    ne = ie - is + 1;
    if (cg_ElementDataSize (cgFile, cgBase, cgZone, ns, &size))
      error_exit ("cg_ElementDataSize");
    conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
    if (conn == NULL)
      error_exit ("malloc failed for element connectivity");
    if (cg_elements_read (cgFile, cgBase, cgZone, ns, conn, NULL))
      error_exit ("cg_elements_read");

    et = elemtype;
    for (i = 0, n = 0; n < ne; n++) {
      if (elemtype == CGNS_ENUMV(MIXED))
        et = (CGNS_ENUMT(ElementType_t))conn[i++];
      switch (et) {
      case CGNS_ENUMV(TRI_3):
      case CGNS_ENUMV(TRI_6):
        nf = 3;
      break;
      case CGNS_ENUMV(QUAD_4):
      case CGNS_ENUMV(QUAD_8):
      case CGNS_ENUMV(QUAD_9):
        nf = 4;
      break;
      default:
        nf = 0;
        break;
      }
      if (nf) {
        face.nnodes = nf;
        for (nn = 0; nn < face.nnodes; nn++)
          face.nodes[nn] = conn[i+nn];
        pf = (FACE *)HashFind(facehash, &face);
        if (NULL == pf) {
          pf = new_face(0);
          pf->nnodes = face.nnodes;
          for (nn = 0; nn < face.nnodes; nn++)
            pf->nodes[nn] = face.nodes[nn];
          HashAdd(facehash, pf);
        }
        pf->id = is + n;
        pf->bcnum = -ns;
      }
      cg_npe (et, &nn);
      i += nn;
    }
    free (conn);
  }

  cgsize_t nFaces = HashSize(facehash);
  FACE ** Faces = (FACE **)malloc(nFaces * sizeof(FACE *));
  if (Faces == NULL)
    error_exit("malloc failed for face list");
  
  ListHash(facehash, Faces);
  HashDestroy(facehash, NULL);
    
  qsort(Faces, nFaces, sizeof(FACE *), sort_faces);

 
  boundary_conditions (nFaces, Faces, cgFile, cgBase, cgZone, surf_ids );
 
  {
    if(boundaries.empty()){
       bc_data.resize(surf_ids.size());
      for(size_t i = 0; i < surf_ids.size(); i++){
        bc_data[i].bc_id = surf_ids[i].first;
      }
    }else{
      vector<int> selected_ids;
      for(size_t i = 0 ; i < boundaries.size(); i++){
       
        for(size_t j = 0; j < surf_ids.size(); j++){
          if(boundaries[i] == surf_ids[j].second)selected_ids.push_back(surf_ids[j].first);
        }
      }
      sort(selected_ids.begin(),selected_ids.end()) ;
      selected_ids.erase(unique(selected_ids.begin(),selected_ids.end()),selected_ids.end()) ;
      bc_data.resize(selected_ids.size());
      for(size_t i = 0; i < selected_ids.size(); i++){
        bc_data[i].bc_id = selected_ids[i];
      }
    }
  }
  
  
  
  int previous_id = Faces[0]->bcnum;
  cgsize_t nf = 0;
  surf_face* the_surf = get_bc(bc_data, previous_id);
  
   while( nf < nFaces) {
     while(nf < nFaces && Faces[nf]->bcnum ==previous_id){
      if(the_surf != NULL){
        if (Faces[nf]->nnodes == 4){
          for(int j = 0; j < 4; j++)
            the_surf->quads.push_back(Faces[nf]->nodes[j]);
        }else{
          for(int j = 0; j < 3; j++)
            the_surf->trias.push_back(Faces[nf]->nodes[j] );
        }
      }
      nf++;
     }
    if(nf < nFaces){
      previous_id = Faces[nf]->bcnum; 
      the_surf = get_bc(bc_data, previous_id);
    }
   }
  
  //clean up
  for(cgsize_t i = 0; i < nFaces; i++){
    free(Faces[i]);
  }
  free(Faces);
}




/*read from cgns file,
  The coordinates are read in as double precision, the format in file can be in single or double precision 
*/
void readCGNS(string filename, int cgbase, int cgzone,  const vector<string>& boundaries,
              vector<double > &pos, vector<surf_face>& bc_data
              ) {

 
  int celldim = 3, phydim = 3;
 
  //sizes: Number of vertices, cells, and boundary vertices in each (index)-dimension.
  //
  // Note that for unstructured grids, the number of cells is the number of highest order elements.
  //Thus, in three dimensions it's the number of 3-D cells, and in two dimensions it's the number of 2-D cells.

  //Also for unstructured grids, if the nodes are sorted between internal nodes and boundary nodes,
  //the optional parameter NBoundVertex must be set equal to the number of boundary nodes.
  //By default, NBoundVertex equals zero, meaning that the nodes are unsorted.

  //Note that a non-zero value for NBoundVertex only applies to unstructured grids.
  //For structured grids, the NBoundVertex parameter always equals 0 in all directions.
  
  cgsize_t sizes[3]; 
  cgsize_t start_index, end_index;
  
  int num_bases=0, num_zones = 0;
  int index_file = 0, index_base=0, index_zone=0;
  cgsize_t num_nodes=0;
  
  cgsize_t errs = 0;
  char bname[33], zname[33];
  CGNS_ENUMT(ZoneType_t) ztype;
    
  cout <<"start reading cgns file" << endl;
  
  if(cg_open (filename.c_str(), CG_MODE_READ, &index_file)) error_exit(" unable to open CGNS grid file ");
  if(cg_nbases (index_file, &num_bases))error_exit("error reading number of bases");
  if(num_bases != 1 && cgbase <= num_bases){
    index_base = cgbase;
  }else{
    index_base =1; //assume only one base and its index is 1
  }
  if(cg_base_read (index_file, index_base, bname, &celldim, &phydim))error_exit("error reading base information");
  if(celldim != 3 || phydim != 3){
    cg_close(index_file);
    cerr << "only 3D cell and physical dimensions are allowed in CGNS file" << endl;
    exit(-1);
  }
  if(cg_nzones (index_file, index_base, &num_zones)) error_exit("error reading number of zones");
  if(num_zones != 1 && cgzone <= num_zones){
    index_zone = cgzone;
  }else{
    index_zone = 1;//assume only one zone and its index is 1
  }
  if(cg_zone_type (index_file, index_base, index_zone, &ztype))error_exit("error reading zone type");
  if (ztype != Unstructured) error_exit("can only handle unstructured grid in CGNS file");
  if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))error_exit("error reading zone information");
    
  //3D unstructured sizes:		NVertex, NCell3D, NBoundVertex
  num_nodes = sizes[0];
  if(num_nodes < 1){
    cg_close(index_file);
    cerr << "number of nodes < 1 " << endl;
    exit(-1);
  }
  
  cgsize_t mxsize = num_nodes;
  pos.resize(mxsize*3) ;
  if(num_nodes > 0){ // Read in positions
   
    start_index = 1;
    end_index = num_nodes;
        
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateX", RealDouble,
                          &start_index, &end_index, &pos[0]);
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateY", RealDouble,
                          &start_index, &end_index, &pos[mxsize]);
    
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateZ", RealDouble,
                          &start_index, &end_index, &pos[2*mxsize]); 
   

    
  }
  
 
//ElementRange contains the index of the first and last elements defined
    //in ElementConnectivity. The elements are indexed with a global numbering
    //system, starting at 1, for all element sections under a given Zone_t data
    //structure. The global numbering insures that each element, whether it's a
    //cell, a face, or an edge, is uniquely identified by its number. They are
    //also listed as a continuous list of element numbers within any single
    //element section. Therefore the number of elements in a section is:

    //ElementSize = ElementRange.end - ElementRange.start + 1
     
 
  
  

  //create a a hash  
  HASH facehash;
  facehash = HashCreate( 127,
                        compare_faces, hash_face);
  if (NULL == facehash)
    error_exit("hash table creation failed");
      
      
 
  
  //read in trias and quads, then hash and sort faces, get boundary information
  
  cout << "start  unstructured_faces " << endl;
  unstructured_faces (facehash, index_file, index_base, index_zone,
                      boundaries,
                      bc_data
                      );

  cout << "end  unstructured_faces " << endl; 
  cout <<"finish reading cgns file" << endl;
        
}



  
void write_surf(const string& filename,
                const vector<double>& pos,
                vector<surf_face>& bc_data
                ){

  vector<int> node_set;
  for(size_t bc=0;bc<bc_data.size();++bc) {
    for(size_t i=0;i<bc_data[bc].trias.size();++i) {
      node_set.push_back(bc_data[bc].trias[i]) ;
    }
    for(size_t i=0;i<bc_data[bc].quads.size();++i) {
      node_set.push_back(bc_data[bc].quads[i]) ;
    }
  }

 

  if(node_set.size() == 0) {
    cerr << "The input boundary name should match the BOCO names in cgns file !" << endl;
    cerr << "nothing to process in surf extract! Exiting."
         << endl ;
    exit(-1) ;
  }
  
  sort(node_set.begin(),node_set.end()) ;
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;

  int mx = node_set[0] ;
  int mn = node_set[0] ;
  for(size_t i=0;i<node_set.size();++i) {
    mx = max(mx,node_set[i]) ;
    mn = min(mn,node_set[i]) ;
  }

  
  vector<int> nmap(mx+1, -1);
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1 ; //assume original numbering start with 1?
  }

  for(size_t bc=0;bc<bc_data.size();++bc) {
    for(size_t i=0;i<bc_data[bc].trias.size();++i) {
      if(nmap[bc_data[bc].trias[i]] < 0)cerr << "nmap invalid for  " << bc_data[bc].trias[i] << endl ;
      bc_data[bc].trias[i] = nmap[bc_data[bc].trias[i]]; 
    }
  
    for(size_t i=0;i<bc_data[bc].quads.size();++i) {
      if(nmap[bc_data[bc].quads[i]] < 0)cerr << "nmap invalid for  " << bc_data[bc].quads[i] << endl ;
      bc_data[bc].quads[i] = nmap[bc_data[bc].quads[i]]; 
    }   
  }

  



  ofstream ssfile(filename.c_str(),ios::out) ;
  ssfile.precision(15) ;

  int ntri = 0 ;
  int nqua = 0 ;
  for(size_t bc=0;bc<bc_data.size();++bc) {
    ntri += bc_data[bc].trias.size()/3;
    nqua += bc_data[bc].quads.size()/4;
  }
 
  ssfile << ntri << ' ' << nqua << ' ' 
         << node_set.size() << endl ;
  
  
  double normal_spacing = 0 ;
  int offset = pos.size()/3;
  for(size_t i=0;i<node_set.size();++i) {
    
    int nd = node_set[i]-1 ;
    ssfile << pos[nd] << ' ' << pos[nd+offset] << ' ' << pos[nd+2*offset]
           << ' '<< normal_spacing << endl ;
  }
  // output triangles
  for(size_t bc=0;bc<bc_data.size();++bc) {
    for(size_t i=0;i<bc_data[bc].trias.size()/3;++i) {
       for(int j = 0 ; j < 3; j++){
         ssfile << bc_data[bc].trias[i*3+j] << ' ';
       }
       ssfile  <<bc_data[bc].bc_id << " 0 1" << endl;
    }
    
  }
// output quads
  for(size_t bc=0;bc<bc_data.size();++bc) {
    for(size_t i=0;i<bc_data[bc].quads.size()/4;++i) {
      for(int j = 0 ; j < 4; j++){
        ssfile << bc_data[bc].quads[i*4+j] << ' ';
      }
      ssfile << bc_data[bc].bc_id << " 0 1" << endl;
    }
  }
 
  ssfile.close();
}


 int main(int ac, char* av[]) {
    const char *filename ; 
    std::string tmp_str ;
    
   int cgbase = 1;
   int cgzone = 1;
   vector<string> bcnames;
   while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.

     if(ac >= 3 && !strcmp(av[1],"-B")){
      cgbase = atoi(av[2]);
      ac -= 2 ;
      av += 2 ;
    }else if(ac >= 3 && !strcmp(av[1],"-Z")){
      cgzone = atoi(av[2]);
      ac -= 2 ;
      av += 2 ;
     }else if(ac >= 3 && !strcmp(av[1],"-bc")){
       bcnames.push_back(string(av[2]));
       ac -= 2 ;
       av += 2 ;
     }
   }
   if(ac == 2) {
     tmp_str.append(av[1]) ;
   }else{
     cerr << "Usage: cgns2surf <options> <file>" << endl
          << "Where options are listed below and <file> is the filename sans postfix" << endl
          << "flags:" << endl
          <<" -bc <string> : boundary name to extract, this option can repeat itself " << endl
          << "  -B <int> : base number in cgns file" << endl
          << "  -Z <int> : zone number in cgns file" << endl
          << endl ;
     exit(-1) ;
   }
   cout << " cgbase " << cgbase << " cgzone " << cgzone << endl;  

  int loc = 0;
  loc = tmp_str.find('.') ;
  std::string new_str = tmp_str.substr(0, loc) ;
  filename = new_str.c_str() ;
  char buf[512] ;
  bzero(buf,512) ;
  snprintf(buf,511,"%s.cgns",filename) ;
  string infile = buf ;
  string outfile = string(filename) + string(".surf") ;
  if(bcnames.size() ==1)outfile = bcnames[0]+string(".surf") ;
  vector<double> pos ;
  vector<surf_face> bc_data;
  readCGNS(infile, cgbase,  cgzone,   bcnames,
           pos,  bc_data); 
  
  //write surf
  write_surf(outfile, pos, bc_data);
 }

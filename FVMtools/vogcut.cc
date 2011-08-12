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
#include "vogtools.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <limits>

using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ifstream;
using Loci::vector3d ;
using Loci::entitySet;
typedef vector3d<double> vect3d;
using std::ofstream;

const int NEW_BC_ID = std::numeric_limits<int>::min();


unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

bool readVolTags(hid_t input_fid,
                 vector<pair<string,Loci::entitySet> > &volDat) {
  using namespace Loci ;
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  vector<pair<string,entitySet> > volTags ;
  hid_t cell_info = H5Gopen(input_fid,"cell_info") ;
  if(cell_info > 0) {
    vector<string> vol_tag ;
    vector<entitySet> vol_set ;
    vector<int> vol_id ;
    
    hsize_t num_tags = 0 ;
    H5Gget_num_objs(cell_info,&num_tags) ;
    for(hsize_t tg=0;tg<num_tags;++tg) {
      char buf[1024] ;
      memset(buf, '\0', 1024) ;
      H5Gget_objname_by_idx(cell_info,tg,buf,sizeof(buf)) ;
      buf[1023]='\0' ;
      
      string name = string(buf) ;
      hid_t vt_g = H5Gopen(cell_info,buf) ;
      hid_t id_a = H5Aopen_name(vt_g,"Ident") ;
      int ident ;
      H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
      H5Aclose(id_a) ;
      entitySet dom ;
      HDF5_ReadDomain(vt_g,dom) ;
      vol_tag.push_back(name) ;
      vol_set.push_back(dom) ;
      vol_id.push_back(ident) ;
      H5Gclose(vt_g) ;
    }
    int maxi =0 ;
    for(size_t i=0;i<vol_id.size();++i)
      maxi = max(vol_id[i],maxi) ;
    vector<pair<string,entitySet> > tmp(maxi+1) ;
    volTags.swap(tmp) ;
    for(size_t i=0;i<vol_id.size();++i) {
      volTags[vol_id[i]].first = vol_tag[i] ;
      volTags[vol_id[i]].second = vol_set[i] ;
    }
  } else {
    hid_t file_info = H5Gopen(input_fid,"file_info") ;
    long numCells = readAttributeLong(file_info,"numCells") ;
    volTags.push_back(pair<string,entitySet>
                      (string("Main"),
                       entitySet(interval(0,numCells-1)))) ;
    H5Gclose(file_info) ;
  }
  
  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);
  volDat.swap(volTags) ;
  return true ;
}

namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern entitySet faceCluster(const multiMap &face2node,
                               const Map &cl, const Map &cr, entitySet faces,
                               vector<unsigned char> &cluster_info,
                               vector<unsigned short> &cluster_sizes);
  extern int fillClusterFaceSizes(unsigned char *cluster, int *sizes) ;
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Map &cl, Map &cr, int face_base) ;

}


//if a node is outside the sphere
inline bool mark_node(const vect3d& p, const vect3d& center, double r){
  return (dot(p-center, p-center) >= r*r);
}

//if a face is outside the sphere:
//if all the nodes of a face are outside the sphere, then the face is outside the sphere
inline bool mark_face(const vector<bool>& nodestat, const Loci::multiMap& face2node, Loci::Entity f){
  
  int fsz = face2node[f].size() ;
  bool isOut = true;
  for(int i = 0; i < fsz; i++){
    if(!nodestat[face2node[f][i]]){
      isOut = false;
      break;
    }
  }
  return isOut;
}



//read in the positions of the nodes 
bool readNodes(hid_t input_fid, 
               vector<vect3d>& pos_dat //pos of all nodes
               ) {
      
  // read in positions
  hid_t fid = H5Gopen(input_fid,"file_info") ;
  unsigned long numNodes = readAttributeLong(fid,"numNodes") ;
  H5Gclose(fid) ;
  
  hsize_t  count = numNodes ;

#ifdef H5_INTERFACE_1_6_4
  hsize_t lstart = 0 ;
#else
  hssize_t lstart = 0 ;
#endif 

  // Read in pos data from file i
  pos_dat.resize(numNodes) ;

  hid_t node_g = H5Gopen(input_fid,"node_info") ;
  hid_t dataset = H5Dopen(node_g,"positions") ;
  hid_t dspace = H5Dget_space(dataset) ;
  
  hsize_t stride = 1 ;
   
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
  int rank = 1 ;
  hsize_t dimension = count ;
  hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
  typedef Loci::data_schema_traits<vect3d > traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  hid_t datatype = dp->get_hdf5_type() ;
  hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
		      &pos_dat[0]) ;
  if(err < 0) {
    cerr << "unable to read positions "<< endl ;
    exit(-1) ;
  }
  H5Sclose(dspace) ;
  H5Dclose(dataset) ;
  H5Gclose(node_g) ;
  return true;
}

bool readTags(string tagfile, vector<bool>& tags){
  
  ifstream inFile;
  inFile.open(tagfile.c_str());
  if(!inFile){
    cerr <<"can not open " << tagfile<< " for input" << endl;
    return false;
  }
  short int  tag; 
  for(unsigned int i = 0; i < tags.size(); i++){
    if(inFile >> tag) tags[i] = (tag != 0);
    else tags[i] = false;
  }
  return true;
}
// bool writeTags(string tagfile, const vector<bool>& tags){
  
//   ofstream outFile;
//   outFile.open(tagfile.c_str());
//   if(!outFile){
//     cerr <<"can not open " << tagfile<< " for output" << endl;
//     return false;
//   }
   
//   for(unsigned int i = 0; i < tags.size(); i++){
//     if(tags[i]) outFile << "1 ";
//     else outFile<<"0 ";
//   }
//   return true;
// }

bool createMaps(hid_t input_fid,
                const vector<bool> nodestat, //the IN/OUT state of all nodes in the original grid
                vector<int>& nodeMap, //the index of nodes
                vector<int>& cellMap, //the index of cells
                unsigned long &newNumNodes, 
                unsigned long &newNumFaces,
                unsigned long &newNumCells,
                vector<vector<int> >& face_data, //the faces of new grids
                Loci::entitySet& bc_ids //the boundary ids of the new grid
                ){
  
  hid_t fid = H5Gopen(input_fid,"file_info") ;
  unsigned long numNodes = readAttributeLong(fid,"numNodes");
  unsigned long numCells = readAttributeLong(fid,"numCells");
  H5Gclose(fid) ;
  
  //initialize maps to -2
  {
    vector<int>(numNodes, -2).swap(nodeMap);
    vector<int>(numCells, -2).swap(cellMap);
  }

    //first time read through face_info to select cells to keep
   
    hid_t face_g  = H5Gopen(input_fid,"face_info") ;
    hid_t  dataset = H5Dopen(face_g,"cluster_sizes") ;
    hid_t dspace = H5Dget_space(dataset) ;
    hsize_t stride = 1 ;
    hsize_t size = 0 ;
    H5Sget_simple_extent_dims(dspace,&size,NULL) ;
    vector<unsigned short> csizes(size) ;
    hsize_t  dimension = size ;
 
#ifdef H5_INTERFACE_1_6_4
    hsize_t start = 0 ;
#else
    hssize_t start = 0 ;
#endif
    hsize_t  count = size ;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
    int rank = 1 ;
    hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
    hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
    if(err < 0) {
      cerr << "unable to read cluster sizes from file" << endl ;
      exit(-1) ;
    }
    H5Dclose(dataset) ;
    H5Sclose(memspace) ;
  
  {
    // Read in clusters 
    dataset = H5Dopen(face_g,"cluster_info") ;
    dspace = H5Dget_space(dataset) ;
    start = 0 ;

    for(size_t c=0;c<size;++c) { // Loop over clusters
      count = csizes[c] ;
      vector<unsigned char> cluster(count) ;
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      dimension = count ;
      memspace = H5Screate_simple(rank,&dimension,NULL) ;
      err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster[0]) ;
      if(err < 0) {
        cerr << "unable to read cluster from file " << endl ;
      }
      start += count ;
      // Now scan cluster into local buffers

      int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
      Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
      Loci::store<int> fcnts ;
      fcnts.allocate(fclust) ;
      Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
      multiMap face2node ;
      Map cl, cr ;
      face2node.allocate(fcnts) ;
      cl.allocate(fclust) ;
      cr.allocate(fclust) ;
      Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0);
      H5Sclose(memspace) ;
    
      // Now loop over faces in cluster to select cells to keep
      for(int f=0;f<nfaces;++f) {
        bool isOut = mark_face(nodestat, face2node, f);
        if(!isOut){//if any node of the face is in, both cl and cr is kept
          cellMap[cl[f]] = -1;
          if(cr[f]>=0) cellMap[cr[f]] = -1;
        }
      }//end loop over faces
    }//end loop over cluster
  
  
  }

  
  //second time read through face_info to select nodes and faces to keep
 

  {
    // Read in clusters 
    dataset = H5Dopen(face_g,"cluster_info") ;
    dspace = H5Dget_space(dataset) ;
    start = 0 ;
    for(size_t c=0;c<size;++c) { // Loop over clusters
      count = csizes[c] ;
      vector<unsigned char> cluster(count) ;
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      dimension = count ;
      memspace = H5Screate_simple(rank,&dimension,NULL) ;
      err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster[0]) ;
      if(err < 0) {
        cerr << "unable to read cluster from file " << endl ;
      }
      start += count ;
      // Now scan cluster into local buffers

      int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
      Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
      Loci::store<int> fcnts ;
      fcnts.allocate(fclust) ;
      Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
      multiMap face2node ;
      Map cl,cr ;
      face2node.allocate(fcnts) ;
      cl.allocate(fclust) ;
      cr.allocate(fclust) ;
      Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0);
      H5Sclose(memspace) ;
    
      // Now loop over faces in cluster to select faces and nodes to keep
      for(int f=0;f<nfaces;++f) {
        if(cr[f] < 0 && cellMap[cl[f]] == -1){ //boundary face in old grid
          bc_ids += -cr[f];
        
          int sz = face2node.num_elems(f); 
          vector<int> tmpVec(sz+2);
          tmpVec[0] = cl[f];
          tmpVec[1] = cr[f];
          for(int i = 0; i < sz; i++){
            tmpVec[i+2] = face2node[f][i];
            nodeMap[face2node[f][i]] = -1;
          }
          face_data.push_back(tmpVec);
        }else if(cr[f] >= 0 &&(cellMap[cl[f]] == -1 || cellMap[cr[f]]==-1)){//inner face in old grid
          if(cellMap[cl[f]] == -1 && cellMap[cr[f]] == -1){//still inner face in new grid
            int sz = face2node.num_elems(f); 
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cl[f];
            tmpVec[1] = cr[f];
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }else if(cellMap[cl[f]] == -1 && cellMap[cr[f]] == -2){//boundary face in new grid, no flipping
          
            int sz = face2node.num_elems(f); 
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cl[f];
            tmpVec[1] = NEW_BC_ID;
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }else if(cellMap[cl[f]] == -2 && cellMap[cr[f]] == -1){//boundary face in new grid,  flipping
            int sz = face2node.num_elems(f);
            for(int n=0;n<sz/2;++n) {
              std::swap(face2node[f][n],face2node[f][sz-1-n]) ;
            }
          
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cr[f];
            tmpVec[1] = NEW_BC_ID;
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }
        }
    
      }//end loop over faces
    }//end loop over cluster
    H5Gclose(face_g) ;
  }
   
  // Count nodes 
  int nnodes = 0 ;
  for(unsigned long i = 0; i < numNodes; i++){
    if(nodeMap[i]==-1) nodeMap[i] = nnodes++;
  }
  newNumNodes = nnodes;
  
  int ncells = 0;
  for(unsigned long i = 0; i < numCells; i++){
    if(cellMap[i]==-1) cellMap[i] = ncells++;
  }
  newNumCells = ncells;
  newNumFaces = face_data.size();

  
  if(nnodes==0 ||ncells==0||newNumFaces==0)return false;  
  else return true; 
}



void Usage() {
  cout << "******************************************************************************" << endl ;

  cout << "vogcut is a utility for extracting  a sphere region from a grid." << endl
       << endl ;
  cout << "Usage: " << endl
       << "  vogcut <options>" << endl
       << "where <options> may be given by:" << endl 
       <<  endl ;
  cout << "  -g <inputfile>           | Specify input grid filename " << endl
       << "  -o <outputfile>          | Specify output grid filename" << endl
       << "  -center <x0 y0 z0>       | Specify center of sphere" << endl
       << "  -r <value>               | Specify radius of sphere" << endl
       << "  -tag <tagfile>           | Specify tag filename " <<endl<<endl;
       
  cout << "example1: vogcut -g grid.vog -o grid_cut.vog -center 0 0 0 -r 1"<<endl;
  cout << "example2: vogcut -g grid.vog -o grid_cut.vog -tag grid.tag"<<endl<<endl;
  cout << "******************************************************************************" << endl ;
}
  
int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << "vogcut is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  string output_file="";
  string input_file="" ;
  string tag_file="";
  vect3d center = vect3d(0.0, 0.0, 0.0);//center of sphere
  double radius  = 1;//radius of sphere
  bool tag_opt = false;
  
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
    if(opt == "-g") {
      i++ ;
      if(i<ac) {
        input_file = string(av[i]) ;
        parsed = true ;
      }
    }
    
    if(opt == "-center") {
      i++ ;
      if((i+2)<ac) {
        center = vect3d(atof(av[i]), atof(av[i+1]), atof(av[i+2])) ;
        parsed = true ;
        i = i+2;
      }
    }
    if(opt == "-r") {
      i++ ;
      if(i<ac) {
        radius = atof(av[i]) ;
        parsed = true ;
      }
    }
    
    if(opt == "-tag") {
      i++ ;
      if(i<ac) {
        tag_file = string(av[i]);
        parsed = true ;
        tag_opt = true;
      }
    }  
           
    if(opt == "-o") {
      i++ ;
      if(i<ac) {
        output_file = string(av[i]) ;
        parsed = true ;
      }
    }
    
    if(!parsed) {
      cerr << "unable to parse command line argument '" << av[i] << "'" << endl ;
      Usage() ;
      exit(-1) ;
    }
  }

   
  if(input_file.empty()) {
    cerr << "no input grid files to process!" << endl ;
    Usage() ;
    exit(-1) ;
  }


  

  bool fail = false ;

#define DEBUG
#ifndef DEBUG
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);
#endif

  {
    size_t p = input_file.find('.');
    if(p == string::npos)input_file = input_file+".vog";
  }
    
  if(output_file==""){
    size_t p = input_file.find('.');
    output_file = input_file.substr(0, p)+"_cut.vog";
  }
  
  
  hid_t input_fid= H5Fopen(input_file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(input_fid <= 0) {
    cerr << "unable to open file '" << input_file << "'"<< endl ;
    exit(-1) ;
    fail = true ;
  }
  if(fail)
    exit(-1) ; 
  
  vector<vect3d> pos; //position of nodes in input grid
  //read in positions
  readNodes(input_fid, 
            pos);
 
  //mark the nodes
  unsigned long inNumNodes = pos.size();
  vector<bool> nodestat(inNumNodes, false);
  if(tag_opt)readTags(tag_file, nodestat);
  else{
    for(unsigned long i = 0; i < inNumNodes; i++)
      nodestat[i] = mark_node(pos[i], center, radius);
    //writeTags(string("grid.tag"), nodestat); 
  }
  
  vector<int> nodeMap;
  vector<int> cellMap;
  unsigned long newNumCells, newNumFaces, newNumNodes;
  vector<vector<int> > face_data;
  Loci::entitySet bc_ids;
  if(!createMaps(input_fid,
                 nodestat,
                 nodeMap,
                 cellMap,
                 newNumNodes,
                 newNumFaces,
                 newNumCells,
                 face_data,
                 bc_ids
                 )){
    cerr<<"ERROR: nothing in the sphere"<< endl;
    return 0;
  }

   
  
  //get new positions
  store<vector3d<double> > act_pos_dat;
  entitySet act_nodes = interval(0, newNumNodes-1);
  act_pos_dat.allocate(act_nodes) ;
  size_t local_count = 0; 
  for(size_t ii=0;ii<pos.size();++ii){
    if(nodeMap[ii]!= -2){
      act_pos_dat[local_count++] = pos[ii] ;
    }
  } 
    
  //define the new boudnary id
    
  int new_bc_id = 1;
  if(bc_ids.size()>0)new_bc_id = bc_ids.Max()+1;
    
  //create maps
  multiMap face2node;
  Map cl, cr;
  store<int> fcount;
  entitySet faces = interval(0, newNumFaces-1);
  cl.allocate(faces);
  cr.allocate(faces);
  fcount.allocate(faces);
  bool new_bc_id_used = false;
    
  for(unsigned long f = 0; f < newNumFaces; f++){
    fcount[f] = face_data[f].size()-2;
    cl[f] = cellMap[face_data[f][0]];
    if(face_data[f][1]>=0) cr[f] =cellMap[face_data[f][1]];
    else cr[f] = face_data[f][1];
    if(cr[f]==NEW_BC_ID){
      cr[f] = -new_bc_id;
      new_bc_id_used = true;
    }
  }
  face2node.allocate(fcount);
    
  for(unsigned long f = 0; f < newNumFaces; f++){
    for(int i = 0; i < fcount[f]; i++){
      face2node[f][i] = nodeMap[face_data[f][i+2]];
    }
  }
  VOG::colorMatrix(act_pos_dat,cl,cr,face2node) ;
    
  // Setup volume tags
    
  vector<pair<string,entitySet> > volTags ;
  vector<pair<string,entitySet> > volDat ;
  readVolTags(input_fid,volDat) ;
  if(!volDat.empty()){
    for(size_t j=0;j<volDat.size();++j){
      string name = volDat[j].first;
      entitySet domain = volDat[j].second;
      entitySet image;
      for(entitySet::const_iterator ei = domain.begin(); ei != domain.end(); ei++){
        if(cellMap[*ei] >=0)image += cellMap[*ei];
      }
      if(image != EMPTY) volTags.push_back(make_pair<string, entitySet>(volDat[j].first, image));
    }
        

    cout << "volume tags =" << endl  ;
    for(size_t i=0;i<volTags.size();++i)
      cout << i << " - " << volTags[i].first << volTags[i].second << endl ;
  }
  //read in boundary names 
  vector<pair<int,string> > boundary_ids;
  Loci::readBCfromVOG(input_file,boundary_ids) ;
  std::map<int,string> bc_map;
  for(unsigned int i = 0; i < boundary_ids.size(); i++){
    bc_map[boundary_ids[i].first] = boundary_ids[i].second;
  }
        
  vector<pair<int,string> > surf_ids ;
  if(new_bc_id !=1){
    bool bc_map_empty = (boundary_ids.size()==0);
    for(entitySet::const_iterator ei=bc_ids.begin(); ei != bc_ids.end(); ei++){
        
      if(!bc_map_empty){
        int id = *ei ;
        std::map<int, string>::const_iterator mi= bc_map.find(id) ;
        if(mi != bc_map.end()){
          surf_ids.push_back(pair<int,string>(*ei,mi->second)) ;
          
        }else{
          int id = *ei ;
          char bcname[512] ;
          sprintf(bcname,"BC_%d",id) ;
          string bcstr(bcname) ;
          surf_ids.push_back(pair<int,string>(*ei,bcstr)) ;
        }
      }else{
        int id = *ei ;
        char bcname[512] ;
        sprintf(bcname,"BC_%d",id) ;
        string bcstr(bcname) ;
        surf_ids.push_back(pair<int,string>(*ei,bcstr)) ;
      }
    }
  }
 
  if(new_bc_id_used){     
    char bcname[512] ;
    sprintf(bcname,"BC_%d",new_bc_id) ;
    string bcstr(bcname) ;
    surf_ids.push_back(pair<int,string>(new_bc_id, bcstr)) ; 
  }
  Loci::writeVOG(output_file, act_pos_dat, cl, cr, face2node,surf_ids, volTags) ;
  // Close everything up
  H5Fclose(input_fid) ;
  Loci::Finalize() ;
  return 0 ;
}
    

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
#include <Loci>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using Loci::vector3d ;
using std::ofstream ;
using std::ios ;


unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}
int  sizeElementType(hid_t group_id, const char *element_name) {
  hid_t dataset = H5Dopen(group_id,element_name) ;
  if(dataset < 0) {
    H5Eclear() ;
    return 0 ;
  }
  hid_t dspace = H5Dget_space(dataset) ;

  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  
  
  H5Dclose(dataset) ;
  return int(size) ;
  
}

template<class T> void readElementType(hid_t group_id, const char *element_name,
                                       vector<T> &v) {
  if(v.size() > 0) {
    hid_t dataset = H5Dopen(group_id,element_name) ;

    typedef Loci::data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }
}





struct surface_info {
  string name ;
  int id ;
  vector<Loci::Array<int,3> > trias ;
  vector<Loci::Array<int,4> > quads ;
  vector<int> nside_sizes;
  vector<int> nside_nodes;
  vector<vector3d<double> > pos;
} ;

//read in a boundary surface
void readSurface(string bc_name,
                 string casename,
                 string iteration,
                 surface_info& surf
                 ) {
  //open the topology file
  string gridtopo;
  if(iteration=="") gridtopo = "output/" + bc_name +"/"+casename+"_bndry.topo" ;
  else  gridtopo = "output/" + bc_name +"/"+casename+"_"+iteration+"_bndry.topo" ;
  cout << "reading topo file : " << gridtopo << endl;

  hid_t fi ;
  hid_t file_id ; 

  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id <= 0) {
    cerr << "unable to open file '" << gridtopo << "'"<< endl ;
  }
  fi = H5Gopen(file_id,bc_name.c_str()) ;


  //get sizes
  int numQuads = sizeElementType(fi,"quads") ;
  int numTrias = sizeElementType(fi,"triangles") ;
  int numGnrls = sizeElementType(fi,"nside_sizes") ;
  int nside_nodes_size = sizeElementType(fi,"nside_nodes") ;
  
  //get surface is
  unsigned long surf_id = readAttributeLong(fi,"id");

  //read in vectors
  surf.trias.resize(numTrias);
  surf.quads.resize(numQuads);
  surf.nside_sizes.resize(numGnrls);
  surf.nside_nodes.resize(nside_nodes_size);
  
  readElementType(fi,"triangles",surf.trias) ;
  readElementType(fi,"quads",surf.quads) ;
  readElementType(fi,"nside_sizes",surf.nside_sizes) ;
  readElementType(fi,"nside_nodes",surf.nside_nodes) ; 
  surf.id = surf_id;  
  surf.name = bc_name;

  //close the file
  H5Gclose(fi) ;
  H5Fclose(file_id) ;

  //open position file
  string gridpos;
  if(iteration=="")gridpos = "output/" + bc_name +"/"+casename+"_bndry.pos" ;
  else gridpos = "output/" + bc_name +"/"+casename+"_"+iteration+ "_bndry.pos"; 
  cout << "reading position file : " << gridpos << endl;
  
  file_id = H5Fopen(gridpos.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id <= 0) {
    cerr << "unable to open file '" << gridpos << "'"<< endl ;
  }
  fi = H5Gopen(file_id,bc_name.c_str()) ;
  
  //read in positions
  int numNodes = sizeElementType(fi,"positions") ;
  surf.pos.resize(numNodes);
  readElementType(fi,"positions",surf.pos) ;
  //close file
  H5Gclose(fi) ;
  H5Fclose(file_id) ;
}
  
//read in a cutting plane
void readSurface_cut(string casename,
                     string iteration,
                     int id,
                     surface_info& surf
                     ) {
  //open the cutting plane file
  string gridtopo;
  if(iteration=="")gridtopo = "output/cut/" +casename+".cut" ;
  else gridtopo = "output/cut/" +casename+"_"+iteration+".cut" ;
  
  hid_t fi ;
  hid_t file_id ; 

  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id <= 0) {
    cerr << "unable to open file '" << gridtopo << "'"<< endl ;
  }
  fi = H5Gopen(file_id,"cuttingplane") ;

  //read in the sizes
  int numGnrls = sizeElementType(fi,"nside_sizes") ;
  int nside_nodes_size = sizeElementType(fi,"nside_nodes") ;

  //read in general faces 
  surf.trias.clear();
  surf.quads.clear();
  surf.nside_sizes.resize(numGnrls);
  surf.nside_nodes.resize(nside_nodes_size);
  
  readElementType(fi,"nside_sizes",surf.nside_sizes) ;
  readElementType(fi,"nside_nodes",surf.nside_nodes) ; 

  surf.id = id;  //assign an id to cutting plane
  surf.name = "cuttingplane"; //assign a name to the cutting plane
  
  //read in num of edge nodes
  int numNodes = sizeElementType(fi,"positions") ;

  //read in edge nodes positions
  surf.pos.resize(numNodes);
  readElementType(fi,"positions",surf.pos) ;

  //read in num of inner_edge nodes
  int numInnerNodes = sizeElementType(fi,"positions2") ;
  
  //read in inner_edge nodes positions
  if(numInnerNodes > 0){
    vector<vector3d<double> > pos2(numInnerNodes);
    readElementType(fi,"positions2",pos2);
    surf.pos.resize(numNodes+numInnerNodes);
    for(int j = 0; j <numInnerNodes; j++){
      surf.pos[j+numNodes] = pos2[j]; 
    }
  }
  
  H5Gclose(fi) ;
  H5Fclose(file_id) ;
}

void writeEnsight(string casename,
                  string iteration,
                  string ofile,
                  vector<surface_info>& surf_list
                  ) {
  enum id_option {OFF, GIVEN, ASSIGN, IGNORE}; 
  id_option node_id_opt = OFF;
  id_option element_id_opt = OFF;
  int  part_id = 1 ;
  FILE *OFP ;
  
  //mkdir {$casename}_ensight.{$iteration}
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
  

  
  string geo_filename =  ofile + ".geo" ;
  string case_filename = dirname + "/" + ofile + ".case" ;
  
  //write case file
  ofstream of(case_filename.c_str(),ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << geo_filename << endl ;
  of.close() ;

  //open geo file
  geo_filename =  dirname + "/" + geo_filename;
  OFP = fopen(geo_filename.c_str(),"wb") ;
  if(OFP == NULL) {
    cerr << "unable to open file '" << geo_filename << "' for writing!" ;
    exit(-1) ;
  }
  
  //write format
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

  //for each surface
  for(unsigned int s = 0; s <surf_list.size(); s++){
    const surface_info& surf = surf_list[s];
    string name = surf.name;
    int npnts = surf.pos.size();

    //create a part
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&part_id, sizeof(int), 1, OFP) ;
  
    //write name
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "%s", name.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  

    //write coordinates
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&npnts, sizeof(int), 1, OFP) ;
    for(int i=0;i<npnts;++i) {
      float x = surf.pos[i].x;
      fwrite(&x,sizeof(float),1,OFP) ;
    }
    for(int i=0;i<npnts;++i) {
      float y = surf.pos[i].y ;
      fwrite(&y,sizeof(float),1,OFP) ;
    }
    for(int i=0;i<npnts;++i) {
      float z = surf.pos[i].z;
      fwrite(&z,sizeof(float),1,OFP) ;
    }
    //write trias
    int ntrias = surf.trias.size();
    if(ntrias > 0){
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "tria3");
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      fwrite(&ntrias, sizeof(int),1,OFP) ;
      fwrite(&surf.trias[0], sizeof(int),ntrias*3,OFP) ;
    }

    //write quads
    int nquads = surf.quads.size();
    if(nquads > 0){
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "quad4");
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      fwrite(&nquads, sizeof(int),1,OFP) ;
      fwrite(&surf.quads[0], sizeof(int),nquads*4,OFP) ;
    }

  
    //write general
    int ngeneral = surf.nside_sizes.size();
    if(ngeneral > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      fwrite(&ngeneral, sizeof(int),1,OFP) ;
    
      int tot = 0 ;
      for(int i=0;i<ngeneral;++i)
        tot += surf.nside_sizes[i] ;
      if(int(surf.nside_nodes.size()) != tot) {
        cerr << "mismatch in node size and faces size " << surf.nside_nodes.size()
             << " was " << tot << endl ;
      }
      fwrite(&surf.nside_sizes[0],sizeof(int),ngeneral,OFP) ;
      fwrite(&surf.nside_nodes[0],sizeof(int),tot,OFP) ;
    }
  
    part_id++ ;
  }
  fclose(OFP) ;
}




// void writeSurf(string filename_surf,
//                const surface_info &tmp_surf,
//                const vector<vector3d<double> > &tmp_p) {
  
//   std::ofstream ofile(filename_surf.c_str(),std::ios::out) ;
  
//   if(ofile.fail()) {
//     cerr << "unable to open output file '" << filename_surf << "'" << endl ;
//     Usage() ;
//   }
 
//   int ntri=tmp_surf.trias.size(), nqua=tmp_surf.quads.size(),ngen=tmp_surf.nside_sizes.size();
 
//   ofile << ntri << ' ' << nqua << ' ' << tmp_p.size() << endl ;
//   ofile.precision(14) ;
  
//   double normal_spacing = 0 ;
//   for(size_t i=0;i<tmp_p.size();++i) {
//     ofile << tmp_p[i].x << ' ' << tmp_p[i].y << ' ' << tmp_p[i].z
// 	  << ' ' << normal_spacing << endl ;
//   }
  
 
//   for(size_t j=0;j<tmp_surf.trias.size();++j){
//     ofile << tmp_surf.trias[j][0] << ' '
//           << tmp_surf.trias[j][1] << ' '
//           << tmp_surf.trias[j][2] << ' '
//           << tmp_surf.id << ' '
//           << 0 << ' '  // reconnection flag
//           << 0 << endl ; // bc flag
//   }
//   // output quad faces
//   for(size_t j=0;j<tmp_surf.quads.size();++j){
//     ofile << tmp_surf.quads[j][0] << ' '
//           << tmp_surf.quads[j][1] << ' '
//           << tmp_surf.quads[j][2] << ' '
//           << tmp_surf.quads[j][3] << ' '
//           << tmp_surf.id << ' '
//           << 0 << ' '  // reconnection flag
//           << 0 << endl ; // bc flag
    
//   }
//   // Now write out general faces
//   if(ngen > 0) {
//     ofile << ngen << endl ;
//     int offset = 0;
//     for(size_t j=0;j<tmp_surf.nside_sizes.size();++j) {
//       size_t nf = tmp_surf.nside_sizes[j] ;
//       ofile << nf ;
//       for(size_t k=0;k<nf;++k)
//         ofile << ' ' << tmp_surf.nside_nodes[offset+k] ;
//       offset += nf;
//       ofile << endl ;
//     }
//   }
//   ofile.close();
// }

void Usage() {
  cout << "Utility for extracting surfaces from chem output" << endl
       << endl ;
  cout << "Usage: " << endl
       << "  extract_chem   -case <casename> [-bc <bc_name>] [-iter <iteration>] [-cut ] [-o <out>]" << endl
       <<  endl ;
  cout << " This program can extract one or more boundary surfaces, and/or a cutting plane" << endl 
       << " The output file is in casename_ensight./casename.case "<< endl;
  exit(-1) ;
}

int main(int ac, char *av[]) {
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  
  
  if(ac < 2){
    Usage();
    exit(0) ;
  }
  
  vector<string> boundaries;
  string bc_name="";
  string casename;
  string iteration="";
  bool cut = false;
  string ofile = "";
  
  
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
    if(opt == "-bc") {
      i++ ;
      if(i<ac) {
        boundaries.push_back(string(av[i])) ;
        parsed = true ;
      }
    }else  if(opt == "-case") {
      i++ ;
      if(i<ac) {
        casename = string(av[i]) ;
        parsed = true ;
      }
    }else  if(opt == "-iter") {
      i++ ;
      if(i<ac) {
        iteration = string(av[i]) ;
        parsed = true ;
      }
    }else  if(opt == "-cut") {
      cut = true;
      parsed = true ;
    
    }else  if(opt == "-o") {
      i++ ;
      if(i<ac) {
        ofile = string(av[i]) ;
        parsed = true ;
      }

    }else{
      parsed = false ;
    }
    
    if(!parsed) {
      cerr << "unable to parse command line argument '" << av[i] << "'" << endl ;
      Usage() ;
    }
  }
  if(casename==""){
    cout<<"ERROR: case name is not specified" << endl;
    Usage() ;
    return -1;
  }
  
 
  int num_boundaries = boundaries.size();
  if(cut) num_boundaries++;
  if(num_boundaries==0){
    Usage() ;
    return -1;
  }
  
  vector<surface_info> surf_list(num_boundaries);
  
 

  for(unsigned int i = 0; i < boundaries.size(); i++){
    string bc_name = boundaries[i];
    readSurface(bc_name,
                casename,
                iteration,
                surf_list[i]);
  }

  if(cut){
    int id = num_boundaries;
    readSurface_cut(casename,
                    iteration,
                    id,
                    surf_list[num_boundaries-1]
                    );
  }
  
  
  if(ofile=="") ofile = casename;

  cout << "writing ensight file: " << ofile+".case and "<< ofile+".geo"<< endl;  
  
  
  writeEnsight(casename, iteration, ofile, surf_list);
  
  //  writeSurf(outfile, surf, pos);
  
 
}

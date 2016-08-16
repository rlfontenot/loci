//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#include <stdlib.h>
#include <math.h>
#include "cgnslib.h"
#include <stdio.h>
#include <string.h>
#include <iostream> 
#include <fstream>
#include <string>
#include <set> 
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
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
using std::ifstream ;
using std::string;
using std::set;
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


string output_dir ;

void Usage(int ac, char *av[]) {
  cerr << av[0] << ": Incorrect Usage" << endl ;
  cout << "Usage:" << endl;
  cout << av[0] << "  <case_name> <iteration> <variable(s)>" << endl ;
  cout << endl ;
  cout << "Variables are defined by the solver, but typically include: " << endl
       << "r     - nodal density" << endl 
       << "p     - nodal log10 pressure" << endl
       << "P     - nodal absolute pressure" << endl
       << "pg    - nodal gage pressure" << endl 
       << "u     - nodal velocity magnitude" << endl
       << "m     - nodal mach number" << endl
       << "t     - nodal temperature" << endl
       << "a     - nodal soundspeed" << endl
       << "f<species> - nodal mass fractions for <species>" << endl
       << "v     - nodal velocity vector"<<endl
       << "x     - x coordinate" << endl 
       << "y     - y coordinate" << endl 
       << "z     - z coordinate" << endl
       << "Boundary Variables:" << endl 
       << "qdot  - wall boundary heat flux" << endl
       << "yplus - wall boundary y plus" << endl
       << "tau   - viscous wall shear stress vector" << endl
       << "tw    - viscous wall boundary temperature" << endl
       << "pw    - viscous wall boundary pressure" << endl
       << "n     - boundary normal (-ascii only)" << endl
       << "area  - boundary area   (-ascii only)" << endl
       << endl ;
  cout << "extra options for controlling extract" << endl
       << "   -dir <directory> : change extract directory from default 'output'"

       << endl << endl ;
  exit(-1) ;
}

struct Zone{
  map<string, int> var2sol;
  bool isVolume;
  int num_quad;
  int num_tria;
  int num_genfc;
};
  


 
//if user-specified variables is empty, collect all variables in file
//otherwise collect the intersection of user-specified variables and file variables
void collect_variables(const vector<string>& variables, //user-specified variables
                       int index_file, int index_base,
                       set<string>& nodal_scalars,
                       set<string>& nodal_vectors,
                       set<string>& element_scalars,
                       set<string>& element_vectors,
                       vector<Zone*>& zone_data 
                       ) {

  char  sname[80], fname[80] , zname[80];
  cgsize_t sizes[3];
  CGNS_ENUMT(GridLocation_t) location;
  CGNS_ENUMT(DataType_t) dtype;
  
  bool collect_all = false;
  if(variables.empty()) collect_all = true;
  int num_zones = 0;
  if(cg_nzones (index_file, index_base, &num_zones)) cg_error_exit();
  cout<<" num_zones: " << num_zones << endl;
  
  zone_data.resize(num_zones);
  for(int index_zone = 1; index_zone <= num_zones; index_zone++){
    if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))cg_error_exit();
    cout << " index_zone : " << index_zone << " zname " << zname << " sizes " << sizes[0] << " " << sizes[1] << " " << sizes[2] << endl;
    Zone* azone = new Zone();
    zone_data[index_zone-1] = azone;
    azone->isVolume = false;
    string partName(zname);
    if(partName.find("Volume") != std::string::npos)azone->isVolume = true; 
    
    map<string, int>& var2sol = azone->var2sol;
    int nsols = 0;
    if(cg_nsols(index_file, index_base, index_zone, &nsols))cg_error_exit();
    for(int index_sol = 1; index_sol <= nsols; index_sol++){
      bool isNodal = true;
      if(cg_sol_info(index_file, index_base, index_zone, index_sol, sname,
                     &location))cg_error_exit();
      if(location != CGNS_ENUMV(Vertex))isNodal = false;
      int nfields = 0;
      if(cg_nfields(index_file, index_base, index_zone, index_sol, &nfields))cg_error_exit();
      for(int index_field = 1; index_field <= nfields; index_field++){
        if(cg_field_info(index_file, index_base, index_zone, index_sol, index_field,
                         &dtype, fname))cg_error_exit();
        string varname = string(fname);
        std::string::reverse_iterator rit=varname.rbegin(); 
        
        bool isVector = false;
        if(*rit=='X' || *rit=='Y' || *rit=='Z'){
          isVector = true;
          varname = varname.substr(0, varname.length()-1);
        }
        if(collect_all){
          if(isNodal){
            if(isVector)nodal_vectors.insert(varname);
            else nodal_scalars.insert(varname);
          }else{
            if(isVector)element_vectors.insert(varname);
            else element_scalars.insert(varname);
          }
          var2sol[varname] = index_sol;
       
        }else{
          for(size_t i = 0; i < variables.size(); i++){
            if(varname == variables[i]){
              if(isNodal){
                if(isVector)nodal_vectors.insert(varname);
                else nodal_scalars.insert(varname);
              }else{
                if(isVector)element_vectors.insert(varname);
                else element_scalars.insert(varname);
              }
              var2sol[varname] = index_sol;
              break;
            }
          }
        }
       
      }
    }
  }
}


      
   
  
/*read grid information and variables from cgns file,
  and then write into enisght file
*/
void get_ensight(const vector<string>& variables,
                 string casename, string iteration,
                 bool id_required) {
  //prepare to write
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
  string cgns_filename = casename +"_" +iteration +".cgns";

  int celldim = 3, phydim = 3;
  cgsize_t sizes[3]; 
  cgsize_t start_index, end_index;
  int nbndry, parent_flag;
  int num_bases=0, num_zones = 0, num_sections = 0;
  int index_file = 0, index_base=0,  index_sect = 0;
   
  char bname[80], zname[80], cname[80], sname[80];
 
  
  CGNS_ENUMT(ZoneType_t) ztype;
  CGNS_ENUMT(ElementType_t) etype;
  CGNS_ENUMT(GridLocation_t) location;
  CGNS_ENUMT(DataType_t) dtype = CGNS_ENUMV(RealSingle);
    
  set<string> nodal_scalars ;
  set<string> nodal_vectors ;
  set<string> element_scalars ;
  set<string> element_vectors ;
  vector<Zone*> zone_data; 
  cout <<"start reading cgns file" << endl;
  if(cg_open (cgns_filename.c_str(), CG_MODE_READ, &index_file)) cg_error_exit();
  if(cg_nbases (index_file, &num_bases))cg_error_exit();
  if(num_bases != 1){
    cout<<" there are " << num_bases << " bases"<< endl;
    cout<< "only read the first one" << endl;
  }
  index_base =1; //assume only one base and its index is 1
  if(cg_base_read (index_file, index_base, bname, &celldim, &phydim))cg_error_exit();
  if(celldim != 3 || phydim != 3)cg_error_exit();
  if(cg_nzones (index_file, index_base, &num_zones)) cg_error_exit();
  //collect variables
  collect_variables(  variables,
                      index_file, index_base,
                      nodal_scalars,
                      nodal_vectors,
                      element_scalars,
                      element_vectors,
                      zone_data);
 
  //write out case file
  ofstream of(case_filename.c_str(),ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << geo_filename << endl ;
  
  geo_filename = dirname + "/"+geo_filename ;
  
  if(nodal_scalars.size()+nodal_vectors.size()+
     element_scalars.size()+element_vectors.size()
     ) {
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
  if(id_required)snprintf(tmp_buf,80, "element id given") ;
  else snprintf(tmp_buf, 80,"element id off") ;  
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;

  
  int part_id = 0 ;
  // Volume Parts Output
  map<int, int> zone2part;
  for(int index_zone = 1; index_zone <= num_zones; index_zone++)
    {
      if(zone_data[index_zone-1]->isVolume){
        if(cg_zone_type (index_file, index_base, index_zone, &ztype))cg_error_exit();
        if (ztype != Unstructured) cg_error_exit();
        if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))cg_error_exit();
        
      
        //3D unstructured sizes:		NVertex, NCell3D, NBoundVertex
        int  num_nodes = sizes[0];
        if(num_nodes < 1){
          cg_close(index_file);
          cerr << "number of nodes < 1 " << endl;
          exit(-1);
        }
        
        //output coordinates
        part_id++ ;
        zone2part[index_zone] = part_id;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf, 80, "part") ;
        fwrite(tmp_buf,sizeof(char), 80, OFP) ;
        fwrite(&part_id,sizeof(int),1,OFP) ;
        fwrite(zname, sizeof(char), 80, OFP) ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf, 80, "coordinates") ;
        fwrite(tmp_buf, sizeof(char), 80, OFP) ;
        fwrite(&num_nodes,sizeof(int),1,OFP) ;
        //read in x
        start_index = 1;
        end_index = num_nodes;
        vector<float> buf(num_nodes);   
        if(cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateX", dtype,
                          &start_index, &end_index, &buf[0]))cg_error_exit();
        fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;
      
        if(cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateY", dtype,
                          &start_index, &end_index, &buf[0]))cg_error_exit();
        fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;
      
        if(cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateZ", dtype,
                          &start_index, &end_index, &buf[0]))cg_error_exit(); 
        fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;
      
      
      
        if(cg_nsections (index_file, index_base, index_zone,
                         &num_sections))cg_error_exit();
      
        if(num_sections <1){
          cg_close(index_file);
          cerr << "number of section < 1 " << endl;
          exit(-1);
        }


        //read in connectivity
        vector<int> tets;
        vector<int> pyramids;
        vector<int > prisms;
        vector<int> hexs;
        vector<int> genface; 
        vector<int> gencell; 
      
       
        for (index_sect = 1; index_sect <= num_sections; ++index_sect)
          {
            
            cgsize_t size = 0;
            cgsize_t *conn = NULL;
            if(cg_section_read (index_file, index_base, index_zone,
                                index_sect, cname, &etype,
                                &start_index, &end_index,
                                &nbndry, &parent_flag))cg_error_exit();
        
            if (cg_ElementDataSize (index_file, index_base, index_zone, index_sect, &size))
              cg_error_exit ();
            
            
            conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
            if (conn == NULL)
              cg_error_exit ();
            if (cg_elements_read (index_file, index_base, index_zone, index_sect, conn, NULL))
              cg_error_exit ();
 
        
            switch (etype) {
            case CGNS_ENUMV(TETRA_4):
              for(cgsize_t i = 0 ; i < size; i++) tets.push_back(conn[i]);
              break;
            case CGNS_ENUMV(PYRA_5):
              for(cgsize_t i = 0 ; i < size; i++) pyramids.push_back(conn[i]);
              break;
            case CGNS_ENUMV(PENTA_6):
              for(cgsize_t i = 0 ; i < size; i++) prisms.push_back(conn[i]);
              break;
            case CGNS_ENUMV(HEXA_8):
              for(cgsize_t i = 0 ; i < size; i++) hexs.push_back(conn[i]);
              break;
            case  CGNS_ENUMV(NGON_n):
              for(cgsize_t i = 0 ; i < size; i++) genface.push_back(conn[i]);
              break;
            case  CGNS_ENUMV(NFACE_n):
              for(cgsize_t i = 0 ; i < size; i++) gencell.push_back(conn[i]);
              break; 
            default:
              // cerr<<"Warning: unknown element type in voulem part: " <<etype <<  endl;
              break;
            }
          }

        cout << " num tets " << tets.size() << endl;
        cout << " num prramids " << pyramids.size() << endl;
        cout << " num prism : "<<  prisms.size() << endl;
        cout << " num hex " <<  hexs.size()<< endl;
        cout << " num genc " << gencell.size() << endl;; 
   
        //write out connectivity
        vector<int> genCellNfaces, genCellNsides,genCellNodes ; 
        char name[80];
        if(tets.size() > 0){
          int num_elem = tets.size()/4;
          memset(name, '\0', 80) ;
          snprintf(name, 80, "tetra4") ;
          fwrite(name, sizeof(char), 80, OFP) ;
          fwrite(&num_elem, sizeof(int), 1, OFP) ;
          if(id_required){//write fake ids
            vector<int> ids(num_elem, 1);
            fwrite(&ids[0],sizeof(int),num_elem,OFP); 
          }
          fwrite(&tets[0],sizeof(int),tets.size(),OFP);
        }
        if(pyramids.size() > 0){
          int num_elem = pyramids.size()/5;
          memset(name, '\0', 80) ;
          snprintf(name,80, "pyramid5") ;
          fwrite(name, sizeof(char), 80, OFP) ;
          fwrite(&num_elem, sizeof(int), 1, OFP) ;
          if(id_required){//write fake ids
            vector<int> ids(num_elem, 1);
            fwrite(&ids[0],sizeof(int),num_elem,OFP); 
          }       
          fwrite(&pyramids[0],sizeof(int),pyramids.size(),OFP);
        }
        if(prisms.size() > 0){
          int num_elem = prisms.size()/6;
          memset(name, '\0', 80) ;
          snprintf(name,80, "penta6") ; 
          fwrite(name, sizeof(char), 80, OFP) ;
          fwrite(&num_elem, sizeof(int), 1, OFP) ;
          if(id_required){//write fake ids
            vector<int> ids(num_elem, 1);
            fwrite(&ids[0],sizeof(int),num_elem,OFP); 
          }       
          fwrite(&prisms[0],sizeof(int),prisms.size(),OFP);
        }
        if(hexs.size() > 0){
          int num_elem = hexs.size()/8;
          memset(name, '\0', 80) ;
          snprintf(name,80, "hexa8") ;
          fwrite(name, sizeof(char), 80, OFP) ;
          fwrite(&num_elem, sizeof(int), 1, OFP) ;
          if(id_required){//write fake ids
            vector<int> ids(num_elem, 1);
            fwrite(&ids[0],sizeof(int),num_elem,OFP); 
          }       
          fwrite(&hexs[0],sizeof(int),hexs.size(),OFP);
        }
        if(gencell.size() > 0){
          size_t face_index = 0;
          vector<int> nnodes;
          vector<int> face_offset;
          while(face_index < genface.size()){
            int num_nodes = genface[face_index];
            nnodes.push_back(num_nodes);
            face_offset.push_back(face_index);
            face_index += num_nodes +1;
          }
           
           
          size_t cell_index = 0;
          face_index = 0;
          while(cell_index < gencell.size()){
            int num_face = gencell[cell_index];
            cell_index++;
            genCellNfaces.push_back(num_face);
            for(int i = 0; i < num_face; i++){
              int face_index = face_offset[gencell[cell_index+i]];
              int num_nodes = genface[face_index];
              face_index++;
              genCellNsides.push_back(num_nodes);
              for(int j = 0; j < num_nodes; j++){
                genCellNodes.push_back(genface[face_index+j]);
              }
            }
           
            cell_index += num_face;
          }
             
          
          memset(name, '\0', 80) ;
          snprintf(name,80, "nfaced") ;
          fwrite(name, sizeof(char), 80, OFP) ;
          int nnf = genCellNfaces.size() ;
          fwrite(&nnf, sizeof(int), 1, OFP) ;
          if(id_required){//write fake ids
            vector<int> ids(nnf, 1);
            fwrite(&ids[0],sizeof(int),nnf,OFP); 
          }
          fwrite(&genCellNfaces[0], sizeof(int), nnf, OFP) ;
          int nnsides = genCellNsides.size() ;
          fwrite(&genCellNsides[0], sizeof(int),nnsides,OFP) ;
          int num_nodes = genCellNodes.size() ;
          fwrite(&genCellNodes[0],sizeof(int),num_nodes,OFP) ;
        }
      }
    }
  
  for(int index_zone = 1; index_zone <= num_zones; index_zone++){ //other zones are surface part
    Zone* z = zone_data[index_zone -1];
    if(!(z->isVolume)){
      if(cg_zone_type (index_file, index_base, index_zone, &ztype))cg_error_exit();
      if (ztype != Unstructured) cg_error_exit();
      if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))cg_error_exit();

    
      //3D unstructured sizes:		NVertex, NCell3D, NBoundVertex
      int num_nodes = sizes[0];
      if(num_nodes < 1){
        cg_close(index_file);
        cerr << "number of nodes < 1 " << endl;
        exit(-1);
      }

      part_id++ ;
      
      zone2part[index_zone] = part_id;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "part") ;
      fwrite(tmp_buf,sizeof(char), 80, OFP) ;
      fwrite(&part_id,sizeof(int),1,OFP) ;
      fwrite(zname, sizeof(char), 80, OFP) ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "coordinates") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      fwrite(&num_nodes,sizeof(int),1,OFP) ;
      //read in x
      start_index = 1;
      end_index = num_nodes;
      
      vector<float> buf(num_nodes);   
      if(cg_coord_read (index_file, index_base, index_zone,
                        "CoordinateX", dtype,
                        &start_index, &end_index, &buf[0]))cg_error_exit();
       

      fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;

      if(cg_coord_read (index_file, index_base, index_zone,
                        "CoordinateY", dtype,
                        &start_index, &end_index, &buf[0]))cg_error_exit();
    
      fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;
      if(cg_coord_read (index_file, index_base, index_zone,
                        "CoordinateZ", dtype,
                        &start_index, &end_index, &buf[0]))cg_error_exit(); 
      fwrite(&buf[0],sizeof(float),num_nodes,OFP) ;

    
    
      if(cg_nsections (index_file, index_base, index_zone,
                       &num_sections))cg_error_exit();
    
      if(num_sections <1){
        cg_close(index_file);
        cerr << "number of section < 1 " << endl;
        exit(-1);
      }

      vector<int > quads;
      vector<int > trias;
      vector<int > genface;
      vector<int> genCellNfaces, genCellNsides,genCellNodes ;
    
      for (index_sect = 1; index_sect <= num_sections; ++index_sect)
        {
          cgsize_t size = 0;
          cgsize_t* conn = NULL;
          if(cg_section_read (index_file, index_base, index_zone,
                              index_sect, cname, &etype,
                              &start_index, &end_index,
                              &nbndry, &parent_flag))cg_error_exit();
         
          if (cg_ElementDataSize (index_file, index_base, index_zone, index_sect, &size))
            cg_error_exit ();
          conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
          if (conn == NULL)
            cg_error_exit ();
          if (cg_elements_read (index_file, index_base, index_zone, index_sect, conn, NULL))
            cg_error_exit ();
 
      
     
          //int num_elem = end_index -  start_index + 1;
          switch (etype) {
          case CGNS_ENUMV(QUAD_4):
            for(int i = 0 ; i < size; i++) quads.push_back(conn[i]);
            break;
          case CGNS_ENUMV(TRI_3):
            for(int i = 0 ; i < size; i++)trias.push_back(conn[i]);
            break;
          case CGNS_ENUMV(NGON_n):
            for(int i = 0 ; i < size; i++) genface.push_back(conn[i]);
            break;
          default:
            cerr<<"Warning: unknown element type in zone " <<index_zone << endl;
            break;
          }
        }
      //read element id

      
      vector<int> quads_ids;
      vector<int> trias_ids;
      vector<int> genface_ids;
      char name[80];
      CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(Integer);
      if(id_required) {
        //if(cg_goto(index_file, index_base, index_zone, "ElementVariables","ElementIds", "end"))cg_error_exit();
      
        int nsol = 0;
        if(cg_nsols(index_file, index_base, index_zone, &nsol))cg_error_exit();
        for(int index_sol = 1; index_sol<= nsol; index_sol++){
          if(cg_sol_info(index_file, index_base, index_zone, index_sol, sname,
                         &location))cg_error_exit();
          if(location==CGNS_ENUMV(CellCenter) && string(sname) == string("ElementVariables")){
            int data_dim;
            cgsize_t dim_vals;
            if(cg_sol_size(index_file, index_base, index_zone, index_sol, &data_dim,
                           &dim_vals))cg_error_exit();
            vector<int> ids(dim_vals);
            cgsize_t range_min = 1 , range_max = dim_vals;
            memset(name, '\0', 80) ;
            sprintf(name, "ElementIds");
            if(cg_field_read(index_file, index_base, index_zone, index_sol, name,
                             datatype, &range_min, &range_max,
                             &ids[0]))cg_error_exit();
            
            size_t ei = 0;
            for(size_t i = 0; i < quads.size(); i++){
              quads_ids.push_back(ids[ei++]);
            }
            for(size_t i = 0; i < trias.size(); i++){
              trias_ids.push_back(ids[ei++]);
            }
            for(size_t i = 0; i < genface.size(); i++){
              genface_ids.push_back(ids[ei++]);
            }
          }
        }
      }

    
      if(quads.size() > 0){
        int num_elem = quads.size()/4;
        memset(name, '\0', 80) ;
        snprintf(name, 80, "quad4") ;
        fwrite(name, sizeof(char), 80, OFP) ;
        fwrite(&num_elem, sizeof(int), 1, OFP) ;
        if(id_required){//write fake ids
          fwrite(&quads_ids[0],sizeof(int),num_elem,OFP); 
        }
        fwrite(&quads[0],sizeof(int),quads.size(),OFP);
      }
      if(trias.size() > 0) {
        int num_elem = trias.size()/3;
        memset(name, '\0', 80) ;
        snprintf(name, 80, "tria3") ;
        fwrite(name, sizeof(char), 80, OFP) ;
        fwrite(&num_elem, sizeof(int),1,OFP) ;
        if(id_required){
          fwrite(&trias_ids[0], sizeof(int),num_elem,OFP) ;
        }
        fwrite(&trias[0],sizeof(int), trias.size(),OFP) ;
      }
      z->num_genfc = 0;
      if(genface.size() > 0) {
        vector<int> nside_sizes,nside_nodes ;
        size_t ei = 0;
        while(ei < genface.size()){
          int num_nodes = genface[ei];
          nside_sizes.push_back(num_nodes);
          ei++;
          for(int i = 0; i < num_nodes; i++){
            nside_nodes.push_back(genface[ei+i]);
          }
          ei += num_nodes;
        }
        z->num_genfc = nside_sizes.size();
        int nside_nodes_size  = nside_nodes.size();
        memset(name, '\0', 80) ;
        snprintf(name, 80, "nsided") ;
        fwrite(name, sizeof(char), 80, OFP) ;
        int ngen = nside_sizes.size() ;
        fwrite(&ngen, sizeof(int),1,OFP) ;
        if(id_required){
          fwrite(&genface_ids[0], sizeof(int),ngen,OFP) ;
        }
        fwrite(&nside_sizes[0],sizeof(int),ngen,OFP) ;
        fwrite(&nside_nodes[0],sizeof(int),nside_nodes_size,OFP) ;
      }
      z->num_quad = quads.size()/4;
      z->num_tria = trias.size()/3;
   
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
    // Loop over volume parts and write out variables for each part if they
    // exist
    

    for(int index_zone = 1; index_zone <= num_zones; index_zone++){
      map<string, int>& var2sol = zone_data[index_zone-1]->var2sol;
      if(var2sol.find(varname) != var2sol.end()){
        int index_sol = var2sol[varname];
        
        int data_dim = 0;
        cgsize_t size = 0;
        if(cg_sol_size( index_file, index_base, index_zone, index_sol, &data_dim,
                        &size))cg_error_exit();
        
        if(data_dim != 1 || size <= 0)cg_error_exit();
        cgsize_t range_min =1, range_max=size;
        vector<float> vals(size);

        
        char name[80];
        memset(name, '\0', 80) ;
        snprintf(name,80,"%s",varname.c_str()) ;
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &vals[0]))cg_error_exit();

       
          
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "part") ;
        fwrite(tmp_buf, sizeof(char), 80, FP) ;
        int tmp = zone2part[index_zone] ;
        fwrite(&tmp, sizeof(int), 1, FP) ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "coordinates") ;
        fwrite(tmp_buf, sizeof(char), 80, FP) ;
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
    // Loop over volume parts and write out variables for each part if they
    // exist

    char name[80];
    for(int index_zone = 1; index_zone <= num_zones; index_zone++){
      map<string, int>& var2sol = zone_data[index_zone-1]->var2sol;
      if(var2sol.find(varname) != var2sol.end()){
        int index_sol = var2sol[varname];
        int data_dim = 0;
        cgsize_t size = 0;
        if(cg_sol_size( index_file, index_base, index_zone, index_sol, &data_dim,
                        &size))cg_error_exit();
        if(data_dim != 1 || size <= 0)cg_error_exit();
        cgsize_t range_min =1, range_max=size;
        vector<float> vals(size); 
          
    
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "part") ;
        fwrite(tmp_buf, sizeof(char), 80, FP) ;
        int tmp = zone2part[index_zone];
        fwrite(&tmp, sizeof(int), 1, FP) ;
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "coordinates") ;
        fwrite(tmp_buf, sizeof(char), 80, FP) ;
        string nameX = varname + "X";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameX.c_str()); 
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &vals[0]))cg_error_exit();
          
        fwrite(&vals[0],sizeof(float),size,FP) ;
          
        string nameY = varname + "Y";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameY.c_str()); 
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &vals[0]))cg_error_exit();
        fwrite(&vals[0],sizeof(float),size,FP) ;
          
        string nameZ = varname + "Z";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameZ.c_str()); 
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &vals[0]))cg_error_exit();
        fwrite(&vals[0],sizeof(float),size,FP) ;
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
    for(int index_zone = 1; index_zone <= num_zones; index_zone++){
      map<string, int>& var2sol = zone_data[index_zone-1]->var2sol;
      if((!zone_data[index_zone-1]->isVolume) && var2sol.find(varname) != var2sol.end()){
        int index_sol = var2sol[varname];
        int data_dim = 0;
        cgsize_t size = 0;
        if(cg_sol_size( index_file, index_base, index_zone, index_sol, &data_dim,
                        &size))cg_error_exit();
        if(data_dim != 1 || size <= 0)cg_error_exit();
        cgsize_t range_min =1, range_max =size;
        vector<float> vals(size);
        char name[80];
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", varname.c_str());
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &vals[0]))cg_error_exit();
        vector<float> qvals, tvals, gvals;

        int count = 0;
        for(int i = 0; i < zone_data[index_zone-1]->num_quad; i++){
          qvals[i] = vals[count++];
        }
        for(int i = 0; i < zone_data[index_zone-1]->num_tria; i++){
          tvals[i] = vals[count++];
        }
        for(int i = 0; i < zone_data[index_zone-1]->num_genfc; i++){
          gvals[i] = vals[count++];
        }
        
        memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = zone2part[index_zone] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;

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

    for(int index_zone = 1; index_zone <= num_zones; index_zone++){
      map<string, int>& var2sol = zone_data[index_zone-1]->var2sol;
      if((!zone_data[index_zone -1]->isVolume) && var2sol.find(varname) != var2sol.end()){
        int index_sol = var2sol[varname];
        int data_dim = 0;
        cgsize_t size = 0;
        if(cg_sol_size( index_file, index_base, index_zone, index_sol, &data_dim,
                        &size))cg_error_exit();
        if(data_dim != 1 || size <= 0)cg_error_exit();
        cgsize_t range_min =1, range_max = size;
        vector<float> xvals(size);
        vector<float> yvals(size);
        vector<float> zvals(size);
        char name[80];
        string nameX = varname+"X";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameX.c_str());
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &xvals[0]))cg_error_exit();

        string nameY = varname+"Y";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameY.c_str());
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &yvals[0]))cg_error_exit();

        string nameZ = varname+"Z";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameZ.c_str());
        if(cg_field_read(index_file, index_base, index_zone, index_sol ,name,
                         dtype, &range_min, &range_max,
                         &zvals[0]))cg_error_exit();


        
        memset(tmp_buf, '\0', 80) ;
        snprintf(tmp_buf,80, "part") ;
        fwrite(tmp_buf, sizeof(char), 80, FP) ;
        int tmp = zone2part[index_zone];
        fwrite(&tmp, sizeof(int), 1, FP) ;
        vector<float> qvals_x, tvals_x, gvals_x ;
        vector<float> qvals_y, tvals_y, gvals_y ;
        vector<float> qvals_z, tvals_z, gvals_z ;
        int count = 0;
        for(int i = 0; i < zone_data[index_zone-1]->num_quad; i++){
          qvals_x[i] = xvals[count];
          qvals_y[i] = yvals[count];
          qvals_z[i] = zvals[count];
          count++;
        }
        for(int i = 0; i < zone_data[index_zone-1]->num_tria; i++){
          tvals_x[i] = xvals[count];
          tvals_y[i] = yvals[count];
          tvals_y[i] = zvals[count];
          count++;
        }

        
        for(int i = 0; i < zone_data[index_zone-1]->num_genfc; i++){
          gvals_x[i] = xvals[count];
          gvals_y[i] = yvals[count];
          gvals_z[i] = zvals[count];
          count++;
        }
        
        
        if(qvals_x.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "quad4") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&qvals_x[0],sizeof(float),qvals_x.size(),FP); 
          fwrite(&qvals_y[0],sizeof(float),qvals_y.size(),FP);  
          fwrite(&qvals_z[0],sizeof(float),qvals_z.size(),FP);  
        }
        if(tvals_x.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "tria3") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&tvals_x[0],sizeof(float),tvals_x.size(),FP); 
          fwrite(&tvals_y[0],sizeof(float),tvals_y.size(),FP);  
          fwrite(&tvals_z[0],sizeof(float),tvals_z.size(),FP); 
        }
        if(gvals_x.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "nsided") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&gvals_x[0],sizeof(float),gvals_x.size(),FP); 
          fwrite(&gvals_y[0],sizeof(float),gvals_y.size(),FP);  
          fwrite(&gvals_z[0],sizeof(float),gvals_z.size(),FP); 
        }
      }
    }
    fclose(FP) ;
  }
  
  

  for(size_t i = 0; i < zone_data.size(); i++)delete zone_data[i];
        
}


int main(int ac, char *av[]) {
  output_dir = "output" ;

 
  string casename ;
  bool found_casename = false ;
  bool found_iteration = false ;
  string iteration ;
  vector<string> variables ;
  

  
 
  bool id_required = true;//ensight has the option to display node and element ids

  vector<string> partlist ;
 
  for(int i=1;i<ac;++i) {
    if(av[i][0] == '-') {
      if(!strcmp(av[i],"-dir")) {
        ++i ;
        output_dir = string(av[i]) ;
      } else {
        cerr << "unknown option " << av[i] << endl ;
        Usage(ac,av) ;
      }
      
    } else {
      if(found_iteration)
        variables.push_back(string(av[i])) ;
      else if(found_casename) {
        iteration = string(av[i]) ;
        found_iteration = true ;
      } else {
        casename = string(av[i]) ;
        found_casename = true ;
      }
    }
  }
  

  // Check for derived variable input requirements
  std::set<string> varset ;
  for(size_t i=0;i<variables.size();++i)
    varset.insert(variables[i]) ;
  if(varset.find("m") != varset.end()) {
    varset.insert("a") ;
    varset.insert("v") ;
  }
  if(varset.find("P") != varset.end()) {
    varset.insert("pg") ;
  }
  if(varset.find("p") != varset.end()) {
    varset.insert("pg") ;
  }
  if(varset.find("u") != varset.end()) {
    varset.insert("v") ;
  }
  if(varset.find("0") != varset.end()) {
    varset.insert("v") ;
  }
  if(varset.find("1") != varset.end()) {
    varset.insert("v") ;
  }
  if(varset.find("2") != varset.end()) {
    varset.insert("v") ;
  }
  vector<string> varlist ;
  std::set<string>::const_iterator vi ;
  for(vi=varset.begin();vi!=varset.end();++vi)
    varlist.push_back(*vi) ;
  variables.swap(varlist) ;
      
  get_ensight(variables,
              casename, iteration,
              id_required);

    
 
} 
    


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

using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using Loci::vector3d ;

struct file_info {
  unsigned long numNodes ;
  unsigned long numFaces ;
  unsigned long numCells ;
  vector<short> cluster_sizes ;
} ;

unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

struct affineMapping {
  double M[4][4] ;
  double determinant() { return
      (M[0][0]*M[1][1]*M[2][2]+
       M[1][0]*M[2][1]*M[0][2]+
       M[2][0]*M[0][1]*M[1][2]) -
      (M[0][0]*M[2][1]*M[1][2]+
       M[1][0]*M[0][1]*M[2][2]+
       M[2][0]*M[1][1]*M[0][2]) ;
  }

  bool leftHanded() {return (determinant() < 0) ; }
  affineMapping() {
    for(int i=0;i<4;++i) {
      for(int j=0;j<4;++j)
        M[i][j] = 0 ;
      M[i][i] = 1 ;
    }
  }
  void Combine(affineMapping a) {
    double Mtmp[4][4] ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) {
        double mtmp = 0 ;
        for(int k=0;k<4;++k)
          mtmp += a.M[i][k]*M[k][j] ;
        Mtmp[i][j] = mtmp ;
      }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        M[i][j] = Mtmp[i][j] ;
  }
  void translate(vector3d<double> tv) {
    affineMapping tmp ;
    tmp.M[0][3] = tv.x ;
    tmp.M[1][3] = tv.y ;
    tmp.M[2][3] = tv.z ;
    Combine(tmp) ;
  }
  void scale(vector3d<double> tv) {
    affineMapping tmp ;
    tmp.M[0][0] = tv.x ;
    tmp.M[1][1] = tv.y ;
    tmp.M[2][2] = tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  sth ;
    tmp.M[2][1] = -sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][2] = -sth ;
    tmp.M[2][0] =  sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][1] =  sth ;
    tmp.M[1][0] = -sth ;
    tmp.M[1][1] =  cth ;
    Combine(tmp) ;
  }
  vector3d<double> Map(vector3d<double> v) {
    double tmp[4] ;
    tmp[0] = v.x ;
    tmp[1] = v.y ;
    tmp[2] = v.z ;
    tmp[3] = 1. ;
    double res[4] ;
    for(int i=0;i<4;++i)
      res[i] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        res[i] += M[i][j]*tmp[j] ;
    vector3d<double> r(res[0],res[1],res[2]) ;
    return r ;
  }
   
} ;

namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern int fillClusterFaceSizes(unsigned char *cluster, int *sizes) ;
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Map &cl, Map &cr, int face_base) ;
  vector<unsigned char>
  encode_face_cluster(const multiMap &face2node,
                      const Map &cl, const Map &cr,
                      entitySet fcluster,
                      entitySet nodeSet,
                      entitySet cellSet) ;
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

void Usage() {
  cout << "Utility for merging vog grids into a single vog file." << endl
       << endl ;
  cout << "Usage: " << endl
       << "  vogmerge <options>" << endl
       << "where <options> may be given by:" << endl 
       <<  endl ;
  cout << "  -g <gridfile>           | Specify input grid filename" << endl
       << "  -o <gridfile>           | Specify output grid filename" << endl
       << "  -xshift <value>         | translate grid in x-coordinate" << endl
       << "  -yshift <value>         | translate grid in y-coordinate" << endl
       << "  -zshift <value>         | translate grid in z-coordinate" << endl
       << "  -xscale <value>         | scale grid in x-coordinate" << endl
       << "  -yscale <value>         | scale grid in y-coordinate" << endl
       << "  -zscale <value>         | scale grid in z-coordinate" << endl
       << "  -scale  <value>         | scale grid in all coordinates" << endl
       << "  -xrotate <value>        | rotate scale about x-axis (degrees)" << endl
       << "  -yrotate <value>        | rotate scale about y-axis (degrees)" << endl
       << "  -zrotate <value>        | rotate scale about z-axis (degrees)" << endl
       << "  -bc <oldname>,<newname> | rename boundary surface" << endl
       << "  -tag <name>             | specify volume tag for input grid" << endl
       << endl ;
}
  
int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << "vogmerge is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  string output_file = "merged.vog" ;
  vector<string> input_files ;
  vector<map<string,string> > bc_rename ;
  vector<vector<pair<string,entitySet> > > volTag ;
  vector<bool> applyMap ;
  vector<affineMapping> gridXform ;
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
    if(opt == "-g") {
      i++ ;
      if(i<ac) {
        input_files.push_back(string(av[i])) ;
        bc_rename.push_back(map<string,string>()) ;
        volTag.push_back(vector<pair<string,entitySet> >()) ;
        applyMap.push_back(false) ;
        gridXform.push_back(affineMapping()) ;
        parsed = true ;
      }
    }
    if(input_files.size() > 0) {
      if(opt == "-xshift") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(atof(av[i]),0,0) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yshift") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(0,atof(av[i]),0) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zshift") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(0,0,atof(av[i])) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-xrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateX(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateY(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateZ(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-xscale") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(atof(av[i]),1.,1.) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yscale") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(1.,atof(av[i]),1.) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zscale") {
        i++ ;
        if(i<ac) {
          vector3d<double> tv(1.,1.,atof(av[i])) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-scale") {
        i++ ;
        if(i<ac) {
          double v = atof(av[i]) ;
          vector3d<double> tv(v,v,v) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }      
      
      if(opt== "-bc") {
        i++ ;
        if(i<ac) {
          char *p = index(av[i],',') ;
          if(p==0) {
            cerr << "-bc option, comma missing in argument" << endl ;
            Usage() ;
            exit(-1) ;
          }
          *p = '\0' ;
          
          string first(av[i]), second(p+1) ;
          bc_rename.back()[first]=second ;
          parsed = true ;
        }
      }
      if(opt=="-tag") {
        i++ ;
        if(i<ac) {
          vector<pair<string,entitySet> > vp ;
          vp.push_back(pair<string,entitySet>(string(av[i]),Loci::EMPTY)) ;
          volTag.back() = vp ;
          parsed = true ;
        }
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

  int num_inputs = input_files.size() ;
  if(num_inputs == 0) {
    cerr << "no input grid files to process!" << endl ;
    Usage() ;
    exit(-1) ;
  }
  vector<hid_t> input_fid(num_inputs) ;

  bool fail = false ;

  hid_t output_fid = 0 ;
  output_fid = H5Fcreate(output_file.c_str(),
                         H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT) ;
  if(output_fid <= 0) {
    cerr << "unable to open output file '" << output_file << "'" << endl ;
    fail = true ;
  }
#define DEBUG
#ifndef DEBUG
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);
#endif
  
  for(int i=0;i<num_inputs;++i) 
    input_fid[i] = H5Fopen(input_files[i].c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  for(int i=0;i<num_inputs;++i)  
    if(input_fid[i] <= 0) {
      cerr << "unable to open file '" << input_files[i] << "'"<< endl ;
      fail = true ;
    }
  if(fail)
    exit(-1) ;
  vector<file_info> fileData(num_inputs) ;
  for(int i=0;i<num_inputs;++i) {
    hid_t fi = H5Gopen(input_fid[i],"file_info") ;

    fileData[i].numNodes = readAttributeLong(fi,"numNodes") ;
    fileData[i].numFaces = readAttributeLong(fi,"numFaces") ;
    fileData[i].numCells = readAttributeLong(fi,"numCells") ;
    
    H5Gclose(fi) ;
  }

  // Setup volume tags
  for(int i=0;i<num_inputs;++i) {
    entitySet compSet = Loci::interval(0,fileData[i].numCells-1) ;
    if(volTag[i].size() == 0) {
      vector<pair<string,entitySet> > volDat ;
      readVolTags(input_fid[i],volDat) ;
      for(size_t j=0;j<volDat.size();++j)
        volTag[i].push_back(volDat[j]) ;
    } else {
      if(volTag[i].size() != 1)
        cerr << "Grid '" << input_files[i] << "' has more than one volume tag!" << endl ;
      volTag[i][0].second = compSet ;
    }
  }
  
  
  vector<vector<pair<int,string> > > boundary_ids(num_inputs) ;
  for(int i=0;i<num_inputs;++i) {
    Loci::readBCfromVOG(input_files[i],boundary_ids[i]) ;
    int sz = boundary_ids[i].size() ;
    map<string,string>::const_iterator mi ;
    for(int bc=0;bc<sz;++bc) {
      if((mi = bc_rename[i].find(boundary_ids[i][bc].second)) != bc_rename[i].end())
        boundary_ids[i][bc].second = mi->second ;
    }
  }    

  map<string,int> newbcs ;
  int bcid = 1 ;
  for(int i=0;i<num_inputs;++i) {
    int bndsi = boundary_ids[i].size() ;
    for(int j=0;j<bndsi;++j) {
      string bnd_name = boundary_ids[i][j].second ;
      if(newbcs.find(bnd_name) == newbcs.end()) {
        newbcs[bnd_name] = bcid++ ;
      }
    }
  }
  
  unsigned long numNodes = 0 ;
  unsigned long numFaces = 0 ;
  unsigned long numCells = 0 ;
  for(int i=0;i<num_inputs;++i) {
    numNodes += fileData[i].numNodes ;
    numFaces += fileData[i].numFaces ;
    numCells += fileData[i].numCells ;
  }


  // Begin writing output file

  {
    // Write out file info
    hid_t file_id = H5Gcreate(output_fid,"file_info",0) ;
    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    hid_t att_id = H5Acreate(file_id,"numNodes",H5T_STD_I64BE,
                             dataspace_id,H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numNodes) ;
    H5Aclose(att_id) ;
    att_id = H5Acreate(file_id,"numFaces", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numFaces) ;
    H5Aclose(att_id) ;
    att_id = H5Acreate(file_id,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numCells) ;
    H5Aclose(att_id) ;
    H5Gclose(file_id) ;
    cout << "numNodes = " << numNodes
         << ", numFaces = " << numFaces
         << ", numCells = " << numCells << endl ;
  }

  {
    // Write out node information
    hid_t node_id = H5Gcreate(output_fid,"node_info",0) ;

    
    
    int rank = 1 ;
    hsize_t dimension = numNodes ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

#ifdef H5_INTERFACE_1_6_4
    hsize_t start = 0 ;
    hsize_t lstart = 0 ;
#else
    hssize_t start = 0 ;
    hssize_t lstart = 0 ;
#endif
    hsize_t stride = 1 ;
    typedef Loci::data_schema_traits<Loci::vector3d<double> > traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    hid_t dataseto = H5Dcreate(node_id,"positions",dp->get_hdf5_type(),
                               dataspace,H5P_DEFAULT) ;
    // Write out node info
    for(int i=0;i<num_inputs;++i) {
      hsize_t count = fileData[i].numNodes ;
      
      // Read in pos data from file i
      vector<Loci::vector3d<double> > pos_dat(fileData[i].numNodes) ;
      hid_t node_g = H5Gopen(input_fid[i],"node_info") ;
      hid_t dataset = H5Dopen(node_g,"positions") ;
      hid_t dspace = H5Dget_space(dataset) ;
      
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
      int rank = 1 ;
      hsize_t dimension = count ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                          &pos_dat[0]) ;
      if(err < 0) {
        cerr << "unable to read positions from '" << input_files[i] << "'" << endl ;
        exit(-1) ;
      }
      H5Sclose(dspace) ;
      H5Dclose(dataset) ;
      H5Gclose(node_g) ;

      if(applyMap[i]) {
        for(size_t ii=0;ii<pos_dat.size();++ii)
          pos_dat[ii] = gridXform[i].Map(pos_dat[ii]) ;
      }
      // Now write out this file's part of the positions
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, &start,&stride,&count,NULL)  ;
      dimension = count ;
      rank = 1 ;
      err = H5Dwrite(dataseto,datatype,memspace,dataspace,
                     H5P_DEFAULT,&pos_dat[0]) ;
      if(err < 0) {
        cerr << "unable to write positions to '" << output_file << "'" << endl ;
        exit(-1) ;
      }
      
      start += count ;
      
      H5Sclose(memspace) ;
    }
    
    H5Sclose(dataspace) ;
    H5Dclose(dataseto) ;
    H5Gclose(node_id) ;
  }

  {
    // Write Face Info

    // First read in clusters from each file and write modified
    // clusters to a scratch file
    FILE *scratch = tmpfile() ;

    // Total output file cluster sizes
    vector<unsigned short> cluster_sizes ;

    int nodeOffset = 0 ;
    int cellOffset = 0 ;
    // Loop over files
    for(int i=0;i<num_inputs;++i) {

      // Generate bc remap ;
      map<int,int> bcr ;
      int bndsi = boundary_ids[i].size() ;
      for(int j=0;j<bndsi;++j) {
        string bnd_name = boundary_ids[i][j].second ;
        int bnd_tag = boundary_ids[i][j].first ;
        bcr[-bnd_tag] = -newbcs[bnd_name] ;
      }        

      
      hid_t face_g = H5Gopen(input_fid[i],"face_info") ;

      // Read cluster sizes
      hid_t dataset = H5Dopen(face_g,"cluster_sizes") ;
      hid_t dspace = H5Dget_space(dataset) ;
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      vector<unsigned short> csizes(size) ;
      hsize_t dimension = size ;
      hsize_t stride = 1 ;
 #ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t count = size ;
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      int rank = 1 ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
      hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
      if(err < 0) {
        cerr << "unable to read cluster sizes from file '" <<
          input_files[i] << "'" << endl ;
        exit(-1) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(memspace) ;

      // Read in clusters and transform
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
          cerr << "unable to read cluster from file '" << input_files[i]
               << "'" << endl ;
        }
        start += count ;

        // Now scan cluster into local buffers

        int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
        Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
        Loci::store<int> fcnts ;
        fcnts.allocate(fclust) ;
        Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
        Loci::multiMap face2node ;
        Loci::Map cl,cr ;
        face2node.allocate(fcnts) ;
        cl.allocate(fclust) ;
        cr.allocate(fclust) ;
        Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0) ;
        if(gridXform[i].leftHanded()) {
          // If coordinate system changes handedness, we need to flip face
          // ordering to get correct normal orientation
          for(int f=0;f<nfaces;++f) {
            int fsz = face2node[f].size() ;
            for(int n=0;n<fsz/2;++n) {
              std::swap(face2node[f][n],face2node[f][fsz-1-n]) ;
            }
          }
        }
        for(int f=0;f<nfaces;++f) {
          int fsz = face2node[f].size() ;
          for(int n=0;n<fsz;++n)
            face2node[f][n] += nodeOffset ;
          cl[f] += cellOffset ;
          if(cr[f] >=0)
            cr[f] += cellOffset ;
          else {
            map<int,int>::const_iterator mi ;
            if((mi=bcr.find(cr[f])) != bcr.end()) {
              cr[f] = mi->second ;
            } else {
              int id = -cr[f] ;
              char bcname[512] ;
              sprintf(bcname,"BC_%d",id) ;
              string bcstr(bcname) ;

              map<string,string>::const_iterator mis ;
              
              if((mis = bc_rename[i].find(bcstr)) != bc_rename[i].end())
                bcstr = mis->second ;
              if(newbcs.find(bcstr) == newbcs.end())
                newbcs[bcstr] = bcid++ ;
              bcr[-id] = -newbcs[bcstr] ;
              cr[f] = bcr[cr[f]] ;
            }
          }
        }
        entitySet cellSet = cl.image(fclust) +cr.image(fclust) ;
        entitySet nodeSet = Loci::MapRepP(face2node.Rep())->image(fclust) ;
        if(cellSet.size() > 256 || nodeSet.size() > 256) {
          cerr << "problem encoding cluster" << endl ;
        }

        vector<unsigned char> clusterout =
          Loci::encode_face_cluster(face2node,cl,cr,fclust,nodeSet,cellSet) ;
        // Write cluster to tmp file
        fwrite(&clusterout[0],clusterout.size(),1,scratch) ;
        unsigned short clsz = clusterout.size() ;
        cluster_sizes.push_back(clsz) ;
      }

      H5Gclose(face_g) ;
      nodeOffset += fileData[i].numNodes ;
      cellOffset += fileData[i].numCells ;
    }

    rewind(scratch) ;

    hid_t face_id = H5Gcreate(output_fid,"face_info",0) ;

    // Write cluster sizes
    int rank = 1 ;
    hsize_t dimension = cluster_sizes.size() ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
#ifdef H5_INTERFACE_1_6_4
    hsize_t start = 0 ;
#else
    hssize_t start = 0 ;
#endif

    hsize_t stride = 1;
    hid_t dataset = H5Dcreate(face_id,"cluster_sizes",H5T_NATIVE_USHORT,
                              dataspace,H5P_DEFAULT) ;
    hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
    H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&dimension,NULL) ;
    hid_t err = H5Dwrite(dataset,H5T_NATIVE_USHORT,memspace,dataspace,H5P_DEFAULT,
                   &cluster_sizes[0]) ;
    if(err<0) {
      cerr << "unable to write cluster sizes to '" << output_file << "'" << endl ;
      exit(-1) ;
    }
    H5Sclose(memspace) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;

    // Write cluster_info

    long cluster_info_size = 0 ;
    for(size_t s=0;s<cluster_sizes.size();++s) 
      cluster_info_size += cluster_sizes[s] ;

    rank = 1 ;
    dimension = cluster_info_size ;
    dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    start = 0 ;
    stride = 1 ;
    dataset = H5Dcreate(face_id,"cluster_info",H5T_NATIVE_UCHAR,dataspace,
                        H5P_DEFAULT) ;
    
    vector<unsigned char > data(8096) ;
    while(cluster_info_size > 0) {
      int sz = 8096 ;
      if(cluster_info_size < 8096)
        sz = cluster_info_size ;
      
      fread(&data[0],sz,1,scratch) ;

      hsize_t dim = sz ;
      hsize_t count = sz ;
      memspace = H5Screate_simple(rank,&dim,NULL) ;
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      err = H5Dwrite(dataset,H5T_NATIVE_UCHAR,memspace,dataspace,H5P_DEFAULT,
                     &data[0]) ;
      if(err < 0) {
        cerr << "unable to write cluster_info" << endl ;
        exit(-1) ;
      }
      H5Sclose(memspace) ;
      cluster_info_size -= sz ;
      start += sz ;
    }

    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;
    H5Gclose(face_id) ;
    fclose(scratch) ;
  }

  cout << "boundary conditions in output are:" << endl ;
  map<string,int>::const_iterator mi ;
  for(mi=newbcs.begin();mi!=newbcs.end();++mi) 
    cout << mi->first <<  ' ' << mi->second << endl ;

  // Write out new boundary conditions
  hid_t surf_id = H5Gcreate(output_fid,"surface_info",0) ;
  for(mi=newbcs.begin();mi!=newbcs.end();++mi) {
    hid_t bc_id = H5Gcreate(surf_id,mi->first.c_str(),0) ;
    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
    hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                             dataspace_id, H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_INT,&(mi->second)) ;
    H5Aclose(att_id) ;
    H5Gclose(bc_id) ;
  }
  H5Gclose(surf_id) ;


  // Compute Volume Tags
  map<string,entitySet> volTagMap ;
  int cellcnt = 0 ;
  for(int i=0;i<num_inputs;++i) {
    for(size_t j=0;j<volTag[i].size();++j) {
      string tag = volTag[i][j].first ;
      entitySet tagSet = (volTag[i][j].second >> cellcnt) ;
      volTagMap[tag] += tagSet ;
    }
    cellcnt+= fileData[i].numCells ;
  }
  vector<pair<string,entitySet> > volTags ;

  if(volTagMap["Main"] != Loci::EMPTY)
    volTags.push_back(pair<string,entitySet>(string("Main"),volTagMap["Main"])) ;
  
                    
  map<string,entitySet>::const_iterator vtmi ;
  for(vtmi=volTagMap.begin();vtmi!=volTagMap.end();++vtmi) {
    if(vtmi->first != "Main")
      volTags.push_back(pair<string,entitySet>(vtmi->first,vtmi->second)) ;
  }

  cout << "volume tags =" << endl  ;
  for(size_t i=0;i<volTags.size();++i)
    cout << i << " - " << volTags[i].first << volTags[i].second << endl ;

  
  hid_t cell_info = H5Gcreate(output_fid,"cell_info", 0) ;

  for(size_t i=0;i<volTags.size();++i) {
    hid_t vol_id = H5Gcreate(cell_info,volTags[i].first.c_str(),0) ;
    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
    hid_t att_id = H5Acreate(vol_id,"Ident", H5T_NATIVE_INT,
                             dataspace_id, H5P_DEFAULT) ;
    int num = int(i) ;
    H5Awrite(att_id,H5T_NATIVE_INT,&num) ;
    H5Aclose(att_id) ;
    Loci::HDF5_WriteDomain(vol_id,volTags[i].second) ;
    H5Gclose(vol_id) ;
  }
  H5Gclose(cell_info) ;
  
  // Close everything up
  for(int i=0;i<num_inputs;++i) 
    H5Fclose(input_fid[i]) ;
  H5Fclose(output_fid) ;
}

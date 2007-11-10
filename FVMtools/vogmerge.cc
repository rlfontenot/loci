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

struct file_info {
  unsigned long numNodes ;
  unsigned long numFaces ;
  unsigned long numCells ;
  vector<short> cluster_sizes ;
} ;

unsigned long readAttributeLong(hid_t group, char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

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

int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  Loci::Init(&ac, &av) ;
  string output_file = "out" ;
  vector<string> input_files ;
  vector<map<string,string> > bc_rename ;
  vector<vector3d<double> > translate ;
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    if(opt == "-g") {
      i++ ;
      if(i<ac) {
        input_files.push_back(string(av[i])) ;
        translate.push_back(vector3d<double>(0,0,0)) ;
        bc_rename.push_back(map<string,string>()) ;
      }
    }
    if(input_files.size() > 0) {
      if(opt == "-xshift") {
        i++ ;
        if(i<ac)
          translate.back().x += atof(av[i]) ;
      }
      if(opt == "-yshift") {
        i++ ;
        if(i<ac)
          translate.back().y += atof(av[i]) ;
      }
      if(opt == "-zshift") {
        i++ ;
        if(i<ac)
          translate.back().z += atof(av[i]) ;
      }
      if(opt== "-bc") {
        i++ ;
        if(i<ac) {
          char *p = index(av[i],',') ;
          if(p==0) {
            cerr << "-bc option, comma missing in argument" << endl ;
            exit(-1) ;
          }
          *p = '\0' ;
          
          string first(av[i]), second(p+1) ;
          bc_rename.back()[first]=second ;
        }
      }
    }
        
    if(opt == "-o") {
      i++ ;
      if(i<ac)
        output_file = string(av[i]) ;
    }
  }

  int num_inputs = input_files.size() ;
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

      for(size_t ii=0;ii<pos_dat.size();++ii)
        pos_dat[ii] += translate[i] ;
      
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
      for(unsigned int c=0;c<size;++c) { // Loop over clusters
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
          encode_face_cluster(face2node,cl,cr,fclust,nodeSet,cellSet) ;
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
    hsize_t lstart = 0 ;
#else
    hssize_t start = 0 ;
    hssize_t lstart = 0 ;
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
    lstart = 0 ;
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
  
  // Close everything up
  for(int i=0;i<num_inputs;++i) 
    H5Fclose(input_fid[i]) ;
  H5Fclose(output_fid) ;
}

#include <hdf5.h>

#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include <rpc/xdr.h>
#include <rpc/rpc.h>
#include "sciTypes.h"
#include "defines.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;
using std::vector;
using Loci::MPI_rank;
using Loci::MPI_processes;
namespace Loci{
 entitySet faceCluster(const multiMap &face2node,
                        const Map &cl, const Map &cr, entitySet faces,
                        vector<unsigned char> &cluster_info,
                       vector<unsigned short> &cluster_sizes) ;
   bool readBCfromVOG(string filename,
                    vector<pair<int,string> > &boundary_ids);


  hid_t writeVOGOpen(string filename);
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids);
  void writeVOGClose(hid_t file_id) ;
}

void writeVOGNode(hid_t file_id,
                Loci::storeRepP &pos,
                  const_store<Loci::FineNodes> &inner_nodes);


void colorMatrix(Map &cl, Map &cr, multiMap &face2node);
namespace Loci{
hid_t writeVOGOpen(string filename) {
    hid_t file_id = 0 ;
    if(MPI_rank==0) 
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
    return file_id ;
  }
   void writeVOGClose(hid_t file_id) {
    if(MPI_rank == 0) H5Fclose(file_id) ;
  }
  
       
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids) {
    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      if(surface_ids.size() != 0) {
        group_id = H5Gcreate(file_id,"surface_info",0) ;
        for(size_t i=0;i<surface_ids.size();++i) {
          hid_t bc_id = 0 ;
          bc_id = H5Gcreate(group_id,surface_ids[i].second.c_str(),0) ;
          hsize_t dims = 1 ;
          hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
          
          hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                                   dataspace_id, H5P_DEFAULT) ;
          H5Awrite(att_id,H5T_NATIVE_INT,&surface_ids[i].first) ;
          H5Aclose(att_id) ;
          H5Gclose(bc_id) ;
        }
        H5Gclose(group_id) ;
      }
    }
  }


}
  void writeVOGFace(hid_t file_id, Map &cl, Map &cr, multiMap &face2node) {
    // Compute cell set
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
  
    Map tmp_cl, tmp_cr;
    multiMap tmp_face2node;

  
    long long local_num_faces = face2node.domain().size()  ;

  
    long long num_cells = geom_cells.size() ;
    long long num_faces = 0 ;

    // Reduce these variables
    MPI_Allreduce(&local_num_faces,&num_faces,1,MPI_LONG_LONG_INT,
                  MPI_SUM,MPI_COMM_WORLD) ;

    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      group_id = H5Gopen(file_id,"file_info") ;

      cerr << "num_cells = " << num_cells << endl
           << "num_faces = " << num_faces << endl ;

      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
      hid_t att_id = H5Acreate(group_id,"numFaces", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_faces) ;
      H5Aclose(att_id) ;
      att_id = H5Acreate(group_id,"numCells", H5T_STD_I64BE,
                         dataspace_id, H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_cells) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;
      group_id = H5Gcreate(file_id,"face_info",0) ;
    }
  
    entitySet faces = face2node.domain() ;
    vector<pair<pair<int,int>, int> > f_ord(faces.size()) ;
    int i = 0 ;
    // For small number of cells, sort to keep bc groupings
    if(num_cells<100000) {
      FORALL(faces,fc) {
        f_ord[i].first.first = cr[fc] ;
        f_ord[i].first.second = cl[fc] ;
        f_ord[i].second = fc ;
        i++ ;
      } ENDFORALL ;
      sort(f_ord.begin(),f_ord.end()) ;
    } else {
      FORALL(faces,fc) {
        f_ord[i].first.first = cl[fc] ;
        f_ord[i].first.second = cr[fc] ;
        f_ord[i].second = fc ;
        i++ ;
      } ENDFORALL ;
    }

    i=0 ;
    store<int> count ;
    count.allocate(faces) ;
    FORALL(faces,fc) {
      int nfc = f_ord[i].second ;
      count[fc] = face2node[nfc].size() ;
      i++ ;
    } ENDFORALL ;
    tmp_face2node.allocate(count) ;
    tmp_cl.allocate(faces) ;
    tmp_cr.allocate(faces) ;
    i=0 ;
  
    int mc = (geom_cells).Min() ;
    // Nodes should be adjusted to start from zero also... for the general case
    FORALL(faces,fc) {
      int nfc = f_ord[i].second ;
      tmp_cl[fc] = cl[nfc]-mc ;
      tmp_cr[fc] = cr[nfc] ;
      if(tmp_cr[fc] >= 0)
        tmp_cr[fc] -= mc ;
      for(int j=0;j<count[fc];++j)
        tmp_face2node[fc][j] = face2node[nfc][j] ;
      i++ ;
    } ENDFORALL ;

    vector<unsigned char> cluster_info ;
    vector<unsigned short> cluster_sizes ;
    while(faces != EMPTY) {
      entitySet fcluster = faceCluster(tmp_face2node,tmp_cl,tmp_cr,faces,
                                       cluster_info,cluster_sizes) ;
      faces -= fcluster ;
    }

    Loci::writeUnorderedVector(group_id,"cluster_sizes",cluster_sizes) ;
    Loci::writeUnorderedVector(group_id,"cluster_info",cluster_info) ;
  
  
    if(MPI_rank == 0) {
      H5Gclose(group_id) ;
    }
  }



// void writeVOG(string filename,store<vector3d<double> > &pos,
//               store<Loci::FineNodes> & inner_nodes,
//               Map &cl, Map &cr, multiMap &face2node,
//               vector<pair<int,string> > surface_ids) {
//     // write grid file
//   hid_t file_id = writeVOGOpen(filename) ;
//   writeVOGSurf(file_id,surface_ids) ;
//   writeVOGNode(file_id,pos) ;
//   writeVOGFace(file_id,cl,cr,face2node) ;
//   writeVOGClose(file_id) ;
// }




class node_output_file : public pointwise_rule {

  store<bool> node_output ;
  const_param<string> outfile_par ;
  const_param<string> meshfile_par ;
  const_store<Loci::FineNodes> inner_nodes;
  
  
  
public:
  node_output_file(){
    name_store("node_output", node_output);
    name_store("outfile_par", outfile_par);
    name_store("meshfile_par", meshfile_par);
     
    name_store("inner_nodes", inner_nodes);
    
 
    
    input("outfile_par");
    input("meshfile_par");
    input("inner_nodes");
 
    
    output("node_output");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    hid_t file_id = Loci::writeVOGOpen(*outfile_par) ;
    vector<pair<int,string> > boundary_ids;
    
    Loci::readBCfromVOG(*meshfile_par, boundary_ids);
    
    Loci::writeVOGSurf(file_id,boundary_ids);
    Loci::storeRepP pos = Loci::exec_current_fact_db->get_variable("pos");
    writeVOGNode(file_id, pos, inner_nodes);
    Loci::writeVOGClose(file_id) ;
  }
};
register_rule<node_output_file> register_node_output_file;

class face_output_file : public pointwise_rule {
 
  store<bool> face_output ;
  const_param<string> outfile_par ;
  const_store<Loci::FineFaces> fine_faces;

  
public:
  face_output_file(){
    name_store("face_output", face_output);
    name_store("outfile_par", outfile_par);
    name_store("fine_faces", fine_faces);
    
    
    input("outfile_par");
    input("fine_faces");
    
    output("face_output");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    //first put fine_faces into a maps
    hid_t file_id = 0 ;

       long numNodes = 0;

     if(MPI_rank==0){
       
       //read numNodes
       long long tmp_numNodes;  
       file_id = H5Fopen((*outfile_par).c_str(),H5F_ACC_RDWR,H5P_DEFAULT) ;
       
       hid_t    group_id = H5Gopen(file_id,"file_info");
       
       hid_t attr = H5Aopen_name(group_id,"numNodes");
       //data type and variable has to match
       hid_t ret  = H5Aread(attr, H5T_NATIVE_LLONG, &tmp_numNodes);
       
       ret =  H5Aclose(attr);
       ret = H5Gclose(group_id);
       numNodes = tmp_numNodes;
     }
     //problem here; MPI has no long long 
     if(MPI_processes > 1) MPI_Bcast(&numNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD);
     
     //compute numFaces
      int local_num_face = 0;
     for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
       local_num_face += fine_faces[*ei].size();
     }

     std::vector<int> face_sizes= Loci::all_collect_sizes(local_num_face);

     long long numFaces = 0;
     for(int i =0; i < MPI_processes; i++) numFaces += face_sizes[i];

     long long cell_min = numNodes + numFaces +1;
     long long face_min = numNodes+1;
     for(int i =0; i < MPI_rank; i++) face_min += face_sizes[i];
     long long face_max = face_min + face_sizes[MPI_rank] -1;

  
     Map cl, cr;
     multiMap face2node;
     store<int> count;
     entitySet faces = interval(face_min, face_max);
     cl.allocate(faces);
     cr.allocate(faces);
     count.allocate(faces);
     
     entitySet::const_iterator fid = faces.begin();
     for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
       for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
         cl[*fid] = fine_faces[*ei][i][0] + cell_min;
         if(fine_faces[*ei][i][1]>=0) cr[*fid] = fine_faces[*ei][i][1] + cell_min;
         else  cr[*fid] = fine_faces[*ei][i][1];
         count[*fid] = fine_faces[*ei][i].size()-2;
         fid++;
       }
     }
     
     face2node.allocate(count);

     fid = faces.begin();
     for(sequence::const_iterator ei = seq.begin(); ei != seq.end(); ei++){
       for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
         for(int j = 0; j < count[*fid]; j++){
           //vog file node index start with 0
           face2node[*fid][j] = fine_faces[*ei][i][j+2]-1;
         }
         fid++;
       }
     }


     colorMatrix(cl, cr, face2node);
     writeVOGFace(file_id, cl, cr, face2node);
     Loci::writeVOGClose(file_id);
  }
};

register_rule<face_output_file> register_face_output_file;

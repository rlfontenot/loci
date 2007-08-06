#include <Loci.h>

#include <vector>
#include <string>
namespace VOG {

  struct BC_descriptor {
    std::string name ;
    int id ;
    bool BC,Visc,Recon,Source,Trans,Rebuild ;
  } ;

  extern std::vector<BC_descriptor> readTags(std::string filename) ;

  using std::vector ;
  // Utility routine for sample sort
  template <class T> void parGetSplitters(vector<T> &splitters,
                                          const vector<T> &input,
                                          MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    splitters = vector<T>(p-1) ;
    vector<T> allsplits(p*(p-1)) ;

    int nlocal = input.size() ;
    if(nlocal < p) {
      std::cerr << "sample sort needs at least p elements per processor"
                << std::endl ;
    }
    for(int i=1;i<p;++i) 
      splitters[i-1] = input[(i*nlocal)/p] ;

    int tsz = sizeof(T) ;
    MPI_Allgather(&splitters[0],(p-1)*tsz,MPI_BYTE,
                  &allsplits[0],(p-1)*tsz,MPI_BYTE,comm) ;
    
    sort(allsplits.begin(),allsplits.end()) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }

  template <class T> void parSampleSort(vector<T> &list, MPI_Comm comm) {
    // First sort list locally
    sort(list.begin(),list.end()) ;

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;
    int me = 0 ;
    MPI_Comm_rank(comm,&me) ;

    vector<T> splitters ;
    parGetSplitters(splitters,list,comm) ;

    int s=0 ;
    vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || list[i] < splitters[s]) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !(list[i] < splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(T) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(T) ;

    vector<T> sorted_pnts(result_size) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    list.swap(sorted_pnts) ;
    sort(list.begin(),list.end()) ;
    return ;
  }
  // Optimize indicies of mesh to increase locality
  extern void optimizeMesh(store<vector3d<double> > &pos,
                           Map &cl, Map &cr, multiMap &face2node) ;
  // Establish geometrically consistent face orientation
  extern void orientFaces(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern void colorMatrix(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern void writeVOG(std::string filename,store<vector3d<double> > &pos,
                       Map &cl, Map &cr, multiMap &face2node) ;
  extern vector<int> simplePartitionVec(int mn, int mx, int p) ;

}


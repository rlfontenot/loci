#include <Loci.h>
#include <parSampleSort.h>

#include <vector>
#include <string>
namespace VOG {

  struct BC_descriptor {
    std::string name ;
    int id ;
    bool Trans ;
  } ;

  extern std::vector<BC_descriptor> readTags(std::string filename) ;

  using std::vector ;

  // Optimize indicies of mesh to increase locality
  extern void optimizeMesh(store<vector3d<double> > &pos,
                           Map &cl, Map &cr, multiMap &face2node) ;
  // Establish geometrically consistent face orientation
  extern void orientFaces(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern void colorMatrix(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern std::vector<int> simplePartitionVec(int mn, int mx, int p) ;

}


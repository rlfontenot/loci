#include <Loci.h>
#include <iostream>
#include <fstream>
#include "vogtools.h"
#include "Tools/debug.h"

#include <map>

#include <sstream>
using std::istringstream ;
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;
using std::ifstream ;
using std::ios ;
#include <ctype.h>

bool checkToken(string s1, string s2) {
  int s1l = s1.length() ;
  int s2l = s2.length() ;
  if(s1l < s2l)
    return false ;
  for(int i=0;i<s2l;++i)
    if(toupper(s1[i]) != toupper(s2[i]))
      return false ;
  return true ;
}
struct quadFace {
  Array<int,4> nodes ;
  int cell ;
  bool left ;
} ;

struct triaFace {
  Array<int,3> nodes ;
  int cell ;
  bool left ;
} ;

inline bool quadCompare(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.nodes[3] < f2.nodes[3])) ;
}

inline bool triaCompare(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2])) ;
}


inline bool quadEqual(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] == f2.nodes[0]) &&
          (f1.nodes[1] == f2.nodes[1]) &&
          (f1.nodes[2] == f2.nodes[2]) &&
          (f1.nodes[3] == f2.nodes[3]));
}

inline bool triaEqual(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] == f2.nodes[0]) &&
          (f1.nodes[1] == f2.nodes[1]) &&
          (f1.nodes[2] == f2.nodes[2]));

}

int main(int ac, char *av[]) {
  using namespace Loci ;
  string Lref = "NOSCALE" ;
  bool optimize = true ;
  Loci::Init(&ac,&av) ;

  if(Loci::MPI_processes != 1) {
    cerr << "msh2vog converter can only run serial at present" << endl ;
  }

  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
	Loci::Finalize() ;
	exit(0) ;
      }
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-in")) {
      Lref = "1 inch" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-ft")) {
      Lref = "1 foot" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-cm")) {
      Lref = "1 centimeter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-m")) {
      Lref = "1 meter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mm")) {
      Lref = "1 millimeter" ;
      ac-- ;
      av++ ;
   } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }

  if(Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units." << endl ;
    exit(-1) ;
  }


  if(ac != 2) {
    cerr << "Usage: cobalt2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -m : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;
      exit(-1) ;
  }

  if(Lref == "")
    Lref = "1 meter" ;
  
  if(!isdigit(Lref[0])) {
    Lref = string("1") + Lref ;
  }

  string filename = av[1] ;
  filename += ".msh" ;

  string outfile = av[1] ;
  outfile += ".vog" ;

  ifstream gridfile(filename.c_str(),ios::in) ;
  if(gridfile.fail()) {
    cerr << "open failed on '" << filename << "'" << endl ;
    exit(-1) ;
  }
  bool binary_format = false ;
  std::map<int,string> surfIds ;
  std::map<int,string> voltags ;
  std::map<int,int> physicalSurfaceTag ;
  std::map<int,int> physicalVolumeTag ;
  string text ;
  vector<vector3d<double> > vpos ;
  int ncells = 0 ;
  int cellid = 1 ;
  vector<size_t> posid ;
  vector<triaFace> tria ;
  vector<quadFace> quad ;

  size_t numNodes=0, minNodeTag=0,maxNodeTag=0 ;

  while(true) {
    getline(gridfile,text) ;
    if(gridfile.eof()) {
      break ;
    }
    if(checkToken(text,"$MeshFormat")) {
      double gmshver = 0 ;
      int bin=0, dsz=0 ;
      gridfile>>gmshver >> bin >> dsz ;
      if(gmshver < 4.0 || gmshver >=5.0) {
	cerr << "gmesh format converter only works for version 4 file format, file version=" << gmshver  << endl ;
	Loci::Abort() ;
      }
      cout << "GMesh Format: " << gmshver << endl ;
      if(bin != 0) {
	binary_format = true ;
	cerr << "binary format not currently implemented" <<endl ;
	Loci::Abort() ;
      }
      
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndMeshFormat")) ;
    }
    if(checkToken(text,"$PhysicalNames")) {
      cout << "---- Begin Physical Names Section ----" << endl ;
      int numNames = 0 ;
      gridfile >> numNames ;
      for(int i=0;i<numNames;++i) {
	int dim=0,tag=0 ;
	string tagname ;
	gridfile >> dim >> tag >> tagname ;
	string dropquotes ;
	for(int i=0;i<tagname.size();++i) {
	  if(tagname[i] != '"')
	    dropquotes += tagname[i] ;
	}
	if(dim == 2) {
	  surfIds[tag] = dropquotes ;
	}
	if(dim == 3) {
	  voltags[tag] = dropquotes ;
	}
      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndPhysicalNames")) ;
      cout << "Named Surface Boundaries:" << endl ;
      for(std::map<int,string>::const_iterator mi=surfIds.begin();mi!=surfIds.end();++mi) {
	cout << mi->second << ":" <<mi->first << endl ;
      }
      if(voltags.begin() != voltags.end()) {
	cout << "Named Volumes:" << endl ;
	for(std::map<int,string>::const_iterator mi=voltags.begin();mi!=voltags.end();++mi) {
	  cout << mi->second << ":" <<mi->first << endl ;
	}
      }
      cout << "---- End Physical Names ----" << endl ;
    }
    if(checkToken(text,"$Entities")) {
      size_t npnts=0,ncurves=0,nsurf=0,nvol=0 ;
      gridfile >> npnts >> ncurves >> nsurf >> nvol ;
      for(size_t i=0;i<npnts;++i) {
	int ptag ;
	double x,y,z ;
	int nphystags=0 ;
	gridfile >> ptag >> x >> y >> z >> nphystags ;
	for(int j=0;j<nphystags;++j) {
	  int val ;
	  gridfile >> val ;
	}
      }
      for(size_t i=0;i<ncurves;++i) {
	int ctag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;
	gridfile >> ctag >> minX >> minY >> minZ ;
	gridfile >> maxX >> maxY >> maxZ >> numPhysicalTags ;
	for(int j=0;j<numPhysicalTags;++j) {
	  int tags ;
	  gridfile >> tags ;
	}
	int nbpts ;
	gridfile >> nbpts ;
	for(int i=0;i<nbpts;++i){
	  int val ;
	  gridfile >> val ;
	}
      }
      for(size_t i=0;i<nsurf;++i) {
	int surfTag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;
	gridfile >> surfTag >> minX >> minY>>minZ>> maxX >> maxY >> maxZ
		 >> numPhysicalTags ;
	physicalSurfaceTag[surfTag] = 0 ;
	for(size_t j=0;j<numPhysicalTags;++j) {
	  int val ;
	  gridfile >> val ;
	  physicalSurfaceTag[surfTag] = val ;
	}
	size_t numBoundingSurfaces ;
	gridfile >> numBoundingSurfaces ;
	for(size_t j=0;j<numBoundingSurfaces;++j) {
	  int val ;
	  gridfile >> val ;
	}
      }
      for(size_t i=0;i<nvol;++i) {
	int volTag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;
	gridfile >> volTag >> minX >> minY>>minZ>> maxX >> maxY >> maxZ
		 >> numPhysicalTags ;
	physicalVolumeTag[volTag] = 0 ;
	for(size_t j=0;j<numPhysicalTags;++j) {
	  int val ;
	  gridfile >> val ;
	  physicalVolumeTag[volTag] = val ;
	}
	size_t numBoundingSurfaces ;
	gridfile >> numBoundingSurfaces ;
	for(size_t j=0;j<numBoundingSurfaces;++j) {
	  int val ;
	  gridfile >> val ;
	}
	
      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndEntities")) ;
      //      cout << "Found Surface Entities:" << endl ;
      //      for(std::map<int,int>::const_iterator mi=physicalSurfaceTag.begin();
      //	  mi!=physicalSurfaceTag.end();++mi) {
      //	if(mi->second >0)
      //	  cout << mi->first << " " << surfIds[mi->second] << endl ;
      //      }
    }    
    if(checkToken(text,"$Nodes")) {
      size_t numEntityBlocks=0 ;
      gridfile >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag ;
      cout << "Nodes: " << numNodes << " [" << minNodeTag << "," << maxNodeTag << "]" << endl ;
      if(maxNodeTag-minNodeTag+1 != numNodes ||
	 minNodeTag != 1) {
	cerr << "Converter expecting continuous node labels starting from 1"
	     << endl ;
	Loci::Abort() ;
      }
      vector<vector3d<double> > tmp(numNodes) ;
      vpos.swap(tmp) ;
      vector<size_t> itmp(numNodes) ;
      posid.swap(itmp) ;
      size_t cnt = 0 ;
      for(size_t i=0;i<numEntityBlocks;++i) {
	int dim=0,etag=0, parametric=0, blksze=0 ;
	gridfile >> dim >> etag >> parametric>> blksze ;
	for(int j=0;j<blksze;++j){
	  size_t nid = 0 ;
	  gridfile >> nid ;
	  posid[cnt+j] = nid ;
	}
	for(int j=0;j<blksze;++j){
	  double x=0,y=0,z=0 ;
	  gridfile >> x >> y >> z ;
	  vpos[cnt+j] = vector3d<double>(x,y,z) ;
	}
	cnt += blksze ;
      }
      if(cnt != numNodes) {
	cerr << "did not read in numNodes cnt=" << cnt << " numNodes=" <<numNodes << endl ;
      }
      //      bool idsok=true ;
      //      for(size_t i=0;i<numNodes;++i)
      //	if(posid[i] != i+1)
      //	  idsok = false ;
      //      if(idsok) {
      //	cout << "node tags consecutive" << endl ;
      //      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndNodes")) ;
    }
    if(checkToken(text,"$Elements")) {
      //      cout << "Elements" << endl ;
      size_t numEntityBlocks=0,numElements=0,minElemTag=0,maxElemTag=0 ;
      gridfile >> numEntityBlocks >> numElements >> minElemTag >> maxElemTag ;
      cout << "Elements: " << numElements << " [" << minElemTag << "," << maxElemTag << "]" << endl ;
      if(maxElemTag-minElemTag+1 != numElements ||
	 minElemTag != 1) {
	cerr << "Converter expecting continuous cell labels starting from 1"
	     << endl ;
	Loci::Abort() ;
      }
      for(size_t i=0;i<numEntityBlocks;++i) {
	int dim=0,tag=0,type=0 ;
	size_t nelem = 0 ;
	gridfile >> dim >> tag >> type >> nelem ;

	for(size_t j=0;j<nelem;++j) {
	  switch(type) {
	  case 1: // 2-node line
	    {
	      int id,pt1,pt2 ;
	      gridfile >> id >> pt1 >> pt2 ;
	    }
	    break ;
	  case 2: // 3-node triangle
	    {
	      int id, pt1,pt2,pt3 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 ;
	      int psid = physicalSurfaceTag[tag] ;
	      if(psid > 0) { // boundary surface
		triaFace tf ;
		tf.nodes[0] = pt1 ;
		tf.nodes[1] = pt2 ;
		tf.nodes[2] = pt3 ;
		tf.left = true ;
		tf.cell = -psid ;
		tria.push_back(tf) ;
	      } else {
		cerr << "id=" << id << ",psid=" << psid << endl ;
	      }
	      
	    }
	    break ;
	  case 3: // 4-node quadrangle
	    {
	      int id, pt1,pt2,pt3,pt4 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 >>pt4 ;
	      int psid = physicalSurfaceTag[tag] ;
	      if(psid > 0) { // boundary surface
		quadFace qf ;
		qf.nodes[0] = pt1 ;
		qf.nodes[1] = pt2 ;
		qf.nodes[2] = pt3 ;
		qf.nodes[3] = pt4 ;
		qf.left = true ;
		qf.cell = -psid ;
		quad.push_back(qf) ;
	      }
	    }
	    break ;
	  case 4: // 4 node tetrahedron
	    {
	      int id, pt1,pt2,pt3,pt4 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 >>pt4 ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt2 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt4 ;
	      tf.nodes[2] = pt3 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt2 ;
	      tf.nodes[2] = pt4 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt2 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt4 ;
	      tria.push_back(tf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 5: // 8 node hexahedron
	    {
	      int id, pt1,pt2,pt3,pt4 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 >>pt4 ;
	      int pt5,pt6,pt7,pt8 ;
	      gridfile >> pt5 >> pt6 >> pt7 >> pt8 ;

	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt3 ;
	      qf.nodes[3] = pt2 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt2 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt5 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt5 ;
	      qf.nodes[2] = pt8 ;
	      qf.nodes[3] = pt4 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt3 ;
	      qf.nodes[2] = pt4 ;
	      qf.nodes[3] = pt8 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt8 ;
	      qf.nodes[2] = pt5 ;
	      qf.nodes[3] = pt6 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt6 ;
	      qf.nodes[2] = pt2 ;
	      qf.nodes[3] = pt3 ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 6: // 6 node prism
	    {
	      int id, pt1,pt2,pt3,pt4,pt5,pt6 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 >>pt4>>pt5 >> pt6 ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt2 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt4 ;
	      tf.nodes[1] = pt5 ;
	      tf.nodes[2] = pt6 ;
	      tria.push_back(tf) ;
	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt2 ;
	      qf.nodes[2] = pt5 ;
	      qf.nodes[3] = pt4 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt3 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt2 ;
	      qf.nodes[1] = pt3 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt5 ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 7: // 5 node pyramid
	    {
	      int id, pt1,pt2,pt3,pt4,pt5 ;
	      gridfile >> id >> pt1 >> pt2 >> pt3 >>pt4>>pt5  ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt2 ;
	      tf.nodes[2] = pt5 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt2 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt3 ;
	      tf.nodes[1] = pt4 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt4 ;
	      tf.nodes[1] = pt1 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt3 ;
	      qf.nodes[3] = pt2 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 15: // 1 point node
	    {
	      int id, pt1 ;
	      gridfile >> id >> pt1 ;
	      cerr << "1 point node not supported" << endl ;
	    }
	  default:
	    cout << "element type= " << type << " unsupported by converter."
		 << endl ;
	    Loci::Abort() ;
	  }
	}
      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndElements")) ;
      cout << "Read " << ncells << " Cells" << endl ;
    }
  }

  // prepare triangle faces (sort them)
  for(size_t i=0;i<tria.size();++i) {
    // pos numbers nodes from zero
    tria[i].nodes[0] -= 1 ;
    tria[i].nodes[1] -= 1 ;
    tria[i].nodes[2] -= 1 ;

    if(tria[i].nodes[0] > tria[i].nodes[1]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[1]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[0] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[1] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[1],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
  }
  //sort tria
  sort(tria.begin(),tria.end(),triaCompare) ;

  //  if(MPI_rank==0)cerr<<" preparing quad faces" << endl;
  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<quad.size();++i) {
    // pos numbers nodes from zero
    quad[i].nodes[0] -=1 ;
    quad[i].nodes[1] -=1 ;
    quad[i].nodes[2] -=1 ;
    quad[i].nodes[3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad[i].nodes[0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad[i].nodes[j]) {
        vs = quad[i].nodes[j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad[i].nodes[(j+nv)&0x3] ;
    // next make orientation so that it will match other face
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[j] ;
    else {
      for(size_t j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[(4 - j) &0x3 ] ;
      quad[i].left = !quad[i].left ;
    }
  }
  sort(quad.begin(),quad.end(),quadCompare) ;

  int matchcnt = 0 ;
  int unmatchcnt = 0 ;
  for(size_t i=0;i+1<tria.size();i++) {
    if(triaEqual(tria[i],tria[i+1])) {
      i++ ;
      matchcnt++ ;
    } else {
      unmatchcnt++ ;
    }
    
  }
  //  cout << "tria matched" << matchcnt << ", unmatched =" << unmatchcnt << endl ;
  
  matchcnt = 0 ;
  unmatchcnt = 0 ;
  for(size_t i=0;i+1<quad.size();i++) {
    if(quadEqual(quad[i],quad[i+1])) {
      i++ ;
      matchcnt++ ;
    } else {
      unmatchcnt++ ;
    }
  }
  //  cout << "quad matched" << matchcnt << ", unmatched =" << unmatchcnt << endl ;
  
  //  cout << "pos.size()=" << vpos.size() << endl ;
  //  cout << "tria.size()=" << tria.size() << endl ;
  //  cout << "quad.size()=" << quad.size() << endl ;

  //due to edge split, trias and quads can turn into general faces
  int ntria = tria.size()/2 ;
  int nquad = quad.size()/2 ;
  int nfaces = ntria+nquad ;
  //  if(MPI_rank==0)cerr<<" creating face2node, cl, cr" << endl;


  int facebase = std::numeric_limits<int>::min()+2048 ;

  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;

  if(Loci::MPI_rank == 0) {
    for(int i=0;i<MPI_processes;++i) 
      if(facesizes[i] == 0) {
	cerr << "Run msh2vog with fewer than " << i << " processors for this mesh!" << endl ;
	Loci::Abort() ;
	break ;
      }
  }
  entitySet faces = interval(facebase,facebase+nfaces-1) ;

  if(nfaces == 0) {
    faces = EMPTY ;
  }
  store<int> count ;
  Map cl,cr ;
  
  cl.allocate(faces) ;
  cr.allocate(faces) ;

  count.allocate(faces) ;
  int fc = facebase ;

  for(int i=0;i<ntria;i++)
    count[fc++] = 3 ;

  for(int i=0;i<nquad;i++)
    count[fc++] = 4 ;


  multiMap face2node ;
  face2node.allocate(count) ;
  fc = facebase ;
  int ccerror = 0 ;
  
  for(int i=0;i<ntria;++i) {
    int corner[3];
    for(int j=0;j<3;++j) {
      corner[j] = tria[i*2].nodes[j] ;
      FATAL(tria[i*2].nodes[j] != tria[i*2+1].nodes[j]) ;
    }
    int c1 = tria[i*2].cell ;
    int c2 = tria[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(tria[i*2+1].left)
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2+1].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2+1].nodes[2-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(tria[i*2].left)
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2].nodes[2-j] ;
    } else {
      if(!tria[i*2].left)
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((tria[i*2].left && tria[i*2+1].left) ||
         (!tria[i*2].left && !tria[i*2+1].left)) {
        cerr << "consistency error" << endl ;
	ccerror++ ;
      }
    }

    for(int j=0;j<3;++j) {
      face2node[fc][j] = corner[j] ;
    }

    fc++ ;
  }

  for(int i=0;i<nquad;++i) {
    int corner[4];
    for(int j=0;j<4;++j) {
      corner[j] = quad[i*2].nodes[j] ;
      FATAL(quad[i*2].nodes[j] != quad[i*2+1].nodes[j]) ;
    }
    int c1 = quad[i*2].cell ;
    int c2 = quad[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }

    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(quad[i*2+1].left)
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2+1].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2+1].nodes[3-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(quad[i*2].left)
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2].nodes[3-j] ;
    } else {
      if(!quad[i*2].left)
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((quad[i*2].left && quad[i*2+1].left) ||
         (!quad[i*2].left && !quad[i*2+1].left)) {
        cerr << "consistency error" << endl ;
	ccerror++ ;
      }
    }

    for(int j=0;j<4;++j) {
      face2node[fc][j] = corner[j] ;
    }
    fc++ ;
  }
  if(ccerror > 0) 
    cerr << "consistency error #=" << ccerror << endl ;
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  store<vector3d<double> > pos ;
  entitySet nodes = interval(0,numNodes-1) ;
  pos.allocate(nodes) ;
  for(size_t i=0;i<numNodes;++i)
    pos[i] = vpos[i] ;

  vector<pair<int,string> > surf_ids ;
  for(std::map<int,string>::const_iterator mi=surfIds.begin();mi!=surfIds.end();++mi) {
    surf_ids.push_back(pair<int,string>(mi->first,mi->second)) ;
  }
  
  VOG::colorMatrix(pos,cl,cr,face2node) ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  //  if(MPI_rank == 0)
    //    cerr << "done  coloring" << endl ;
  if(optimize) {
    if(MPI_rank == 0)
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;
  Loci::Finalize() ;
  return 0 ;
}  

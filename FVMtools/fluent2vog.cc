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

/**

fluent2vog converter program:

This program converts fluent mesh files (with a '.msh') extension into
a VOG format file that is suitable for use with Loci programs.  To use this
program enter:

fluent2vog <units option> <-o> <casename>

where units option can be:

-m : meters
-cm : centemeters
-mm : millimeters
-in : inches
-ft : feet
-Lref "reference length" : user specified scaling factor

The "-o" option, if provided, will disable node reordering optimizations.

The <casename> is the base name without extension of the fluent mesh file.

This converter is able to convert both 2-D and 3-D meshes.  For 2-D
meshes, the converter automatically extrudes the mesh once cell in the
Z coordinate direction.  Symbolic boundary names are extracted from the
Fluent mesh file and used as symbolic boundary face names in the VOG file.
Volume tagged cells are ignored.  Periodic boundary maps, hanging nodes, and
other features are ignored by the converter.

 **/

#include <Loci.h>
#include "vogtools.h"
#include <iostream>
#include <cstdlib>
#include <string>

double extrude_dist = 0.01 ; // In 2-D grids, extrude 1 cm

using namespace std ;


istream &getOpenParen(istream &s) {
  while(s.peek() != '(' && !s.eof())
    s.get() ;
  s.get() ;
  return s ;
}

istream &getCloseParen(istream &s) {
  int count = 0 ;
  while(s.peek() != ')' && !s.eof()) {
    int c = s.get() ;
    if(c<0){
      cerr << "failed to find closing parenthesis" << endl ;
      exit(-1)  ;
    }
    if(c == '(') {
      count++ ;
      while(count > 0) {
	c = s.get() ;
        if(c<0) {
          cerr << "failed to find closing parenthesis" << endl;
          exit(-1) ;
        }

	if(c == '(')
	  count++ ;
	if(c == ')') 
	  count-- ;
      }
    }
  }
  s.get() ;
  return s ;
}

int getDecimal(istream &s) {
  int val = 0 ;
  while(s.peek() == ' ' || s.peek() == '\t' || 
	s.peek() == '\r' || s.peek() == '\n') 
    s.get() ;
  while(s.peek() >= '0' && s.peek() <= '9') {
    char c = s.get() ;
    int d = c - '0' ;
    val = val*10 + d ;
  }
  return val ;
}

int getHex(istream &s) {
  int val = 0 ;
  while(s.peek() == ' ' || s.peek() == '\t' || 
	s.peek() == '\r' || s.peek() == '\n') 
    s.get() ;
  while((s.peek() >= '0' && s.peek() <= '9') ||
	(s.peek() >= 'a' && s.peek() <= 'f') ||
	(s.peek() >= 'A' && s.peek() <= 'F')) {
    char c = s.get() ;
    int h = c - '0' ;
    if(c >= 'a' && c <= 'f') 
      h = (c - 'a') + 10 ;
    if(c >= 'A' && c <= 'F')
      h = (c - 'A') + 10 ;
    val = val*16+h ;
  }
  return val ;
}


string getString(istream &s) {
  while(s.peek() == ' ' || s.peek() == '\t' || 
	s.peek() == '\r' || s.peek() == '\n') 
    s.get() ;
  string str ;
  if(s.peek() == '"') {
    s.get() ;
    while(s.peek() != '"') {
      char c = s.get() ;
      str += c ;
    }
    s.get() ;
  }
  return str ;
}

string getName(istream &s) {
  while(s.peek() == ' ' || s.peek() == '\t' || 
	s.peek() == '\r' || s.peek() == '\n') 
    s.get() ;
  string str ;
  while((s.peek() >= 'a' && s.peek() <='z')||
	(s.peek() >= 'A' && s.peek() <='Z')||
	(s.peek() >= '0' && s.peek() <='9')||
 	s.peek() == '_' || s.peek() == '-'){
    char c = s.get() ;
    if(c == '-')
      c = '_' ;
    str += c ;
  }

  return str ;
}    

int getZoneInfo(istream &s,int &start, int &end, int &type, int &dim) {
  getOpenParen(s) ;
  int zone = getHex(s) ;
  start = getHex(s) ;
  end = getHex(s) ;
  type = getHex(s) ;
  dim = getHex(s) ;
  getCloseParen(s) ;
  return zone ;
}

int scanFluentFile(string filename,
                   store<vector3d<double> > &pos,
                   Map &cl, Map &cr, multiMap &face2node,
                   vector<VOG::BC_descriptor> &bcs) {
  int maxzoneid = 0 ;
  
  map<int,int> zone_map ;
  store<int> count ;
  store<Loci::Array<int,4> > face_info ;
  ifstream s(filename.c_str(),ios::in) ;
  if(s.fail()) {
    return -1 ;
  }
  vector<pair<pair<int,int>,pair<int,int> > > edge_list ;

  int n2dpos = 0;
  int dimension = 0 ;
  while(!s.eof()) {
    getOpenParen(s) ;
    int code = getDecimal(s) ;
    switch(code) {
    case 0: // comment
      getString(s) ;
      getCloseParen(s) ;
      break ;
    case 1: // Header
      getString(s) ;
      getCloseParen(s) ;
      break ;
    case 2: // Dimension
      dimension = getDecimal(s) ;
      if(dimension == 2) {
        cout << "Converting 2-D mesh, extruding creating two symmetry planes."
             << endl ;
      } else if(dimension != 3) {
        cerr << "Don't know how to process a "
             << dimension << " dimensional mesh." << endl ;
        exit(-1) ;
      } 
      break ;
    case 10:
      {
	int start,end,type,dim ;
	int zone = getZoneInfo(s,start,end,type,dim) ;
	if(type != 1) {
	  cerr << "can only handle type 1 nodes" << endl;
          //	  exit(-2) ;
	}
        if(zone == 0) { // allocate pos ;
          if(dimension == 2) {
            n2dpos = end-start+1 ;
            entitySet dom = interval(start-1,end-1+n2dpos) ;
            pos.allocate(dom) ;
            cout << "reading " << end << " nodes." << endl ;
          } else {
            entitySet dom = interval(start-1,end-1) ;
            pos.allocate(dom) ;
            cout << "reading " << end << " nodes." << endl ;
          }
        } else {
          if(pos.domain() == EMPTY) {
            cerr << "No zone 0 to set node size!" << endl ;
            return -1 ;
          }
          getOpenParen(s) ;
          if(dimension == 2) {
            for(int i=start-1;i<end;++i) {
              s >> pos[i].x >> pos[i].y ;
              pos[i].z = -0.5*extrude_dist ;
            }
            for(int i=start-1;i<end;++i) {
              pos[i+n2dpos]= pos[i] ;
              pos[i+n2dpos].z = 0.5*extrude_dist;
            }
          } else {
            for(int i=start-1;i<end;++i) {
              s >> pos[i].x >> pos[i].y >> pos[i].z ;
            }
          }
          getCloseParen(s) ;
        }
        
	getCloseParen(s) ;
      }
      break ;
    case 18: // Periodic BC
      cerr << "warning, periodic zone ignored" << endl ;
      getCloseParen(s) ;
      break ;
    case 12: // Cells
      {
	int start, end, type, dim ;
	int zone = getZoneInfo(s,start,end,type,dim) ;
	getCloseParen(s) ;
        if(zone == 0)
          cout << "reading " << end-start+1 << " cells." << endl ;
      }
      break ;
    case 13: // Faces
      {
	int start, end, type, dim ;
	int zone = getZoneInfo(s,start,end,type,dim) ;
        maxzoneid = max(maxzoneid,zone) ;
        if(zone == 0) {
          if(dimension == 2) {
            cout << "reading " << end-start+1 << " edges." << endl ;
          } else {
            entitySet dom = interval(start-1,end-1) ;
            cl.allocate(dom) ;
            cr.allocate(dom) ;
            count.allocate(dom) ;
            face_info.allocate(dom) ;
          }
        } else {
          if(dimension == 2) {
            getOpenParen(s) ;
            if(dim == 0) { // mixed element type (better be edges)
              for(int i=start-1;i<end;++i) {
                int nfaces = getHex(s) ;
                if(nfaces != 2) {
                  cerr << "faces should be edges in 2-D mesh." << endl ;
                  exit(-1) ;
                }
                int n1 = getHex(s)-1 ;
                int n2 = getHex(s)-1 ;

                int c1 = getHex(s) ; // cell left side
                int c2 = getHex(s) ; // cell right side
                if(c2 == 0) {
                  c2 = -zone ;
                  zone_map[zone] = 1 ;
                } else if(c1 == 0) {
                  c1 = -zone ;
                  zone_map[zone] = 1 ;
                  std::swap(c1,c2) ;
                  std::swap(n1,n2) ;
                }
                pair<int,int> nc = pair<int,int>(n1,n2) ;
                pair<int,int> cc = pair<int,int>(c1,c2) ;
                edge_list.push_back(pair<pair<int,int>,pair<int,int> >(nc,cc)) ;
              }
            } else if(dim == 2) { // edges
              for(int i=start-1;i<end;++i) {
                int n1 = getHex(s)-1 ;
                int n2 = getHex(s)-1 ;

                int c1 = getHex(s) ; // cell left side
                int c2 = getHex(s) ; // cell right side
                if(c2 == 0) {
                  c2 = -zone ;
                  zone_map[zone] = 1 ;
                } else if(c1 == 0) {
                  c1 = -zone ;
                  zone_map[zone] = 1 ;
                  std::swap(c1,c2) ;
                  std::swap(n1,n2) ;
                }
                pair<int,int> nc = pair<int,int>(n1,n2) ;
                pair<int,int> cc = pair<int,int>(c1,c2) ;
                edge_list.push_back(pair<pair<int,int>,pair<int,int> >(nc,cc)) ;
              }
            } else {
              cerr << "unsupported face type " << endl ;
              return -1 ;
            }
            getCloseParen(s) ;
          } else {
            getOpenParen(s) ;
            if(dim == 0) { // mixed element type
              for(int i=start-1;i<end;++i) {
                int nfaces = getHex(s) ;
                for(int j=0;j<nfaces;++j)
                  face_info[i][j] = getHex(s) ;
                cl[i] = getHex(s) ;
                cr[i] = getHex(s) ;
                if(cr[i] == 0) {
                  cr[i] = -zone ;
                  zone_map[zone] = 1 ;
                } else if(cl[i] == 0) {
                  cl[i] = -zone ;
                  zone_map[zone] = 1 ;
                  std::swap(cl[i],cr[i]) ;
                  for(int j=0;j<nfaces/2;++j)
                    std::swap(face_info[i][j],face_info[i][nfaces-j-1]) ;
                }
                count[i] = nfaces ;
              }
            } else if(dim == 3) { // triangles
              for(int i=start-1;i<end;++i) {
                int nfaces = 3 ;
                for(int j=0;j<nfaces;++j)
                  face_info[i][j] = getHex(s) ;
                cl[i] = getHex(s) ;
                cr[i] = getHex(s) ;
                if(cr[i] == 0) {
                  cr[i] = -zone ;
                  zone_map[zone] = 1 ;
                } else if(cl[i] == 0) {
                  cl[i] = -zone ;
                  zone_map[zone] = 1 ;
                  std::swap(cl[i],cr[i]) ;
                  for(int j=0;j<nfaces/2;++j)
                    std::swap(face_info[i][j],face_info[i][nfaces-j-1]) ;
                }
                count[i] = nfaces ;
              }
            } else if(dim == 4) { // quads
              for(int i=start-1;i<end;++i) {
                int nfaces = 4 ;
                for(int j=0;j<nfaces;++j)
                  face_info[i][j] = getHex(s) ;
                cl[i] = getHex(s) ;
                cr[i] = getHex(s) ;
                if(cr[i] == 0) {
                  cr[i] = -zone ;
                  zone_map[zone] = 1 ;
                } else if(cl[i] == 0) {
                  cl[i] = -zone ;
                  zone_map[zone] = 1 ;
                  std::swap(cl[i],cr[i]) ;
                  for(int j=0;j<nfaces/2;++j)
                    std::swap(face_info[i][j],face_info[i][nfaces-j-1]) ;
                }
                count[i] = nfaces ;
              }
            } else {
              cerr << "unsupported face type " << endl ;
              return -1 ;
            }
            getCloseParen(s) ;
          }
        }
	getCloseParen(s) ;
      }
      break ;
    case 45: // Boundary names
      {
	getOpenParen(s) ;
	int zoneId = getDecimal(s) ;
        maxzoneid = max(maxzoneid,zoneId) ;
	string type = getName(s) ;
	string name = getName(s) ;
        if(zone_map.find(zoneId) != zone_map.end()) {
          VOG::BC_descriptor tmp ;
          tmp.name = name ;
          tmp.id = zoneId ;
          tmp.Trans = false ;
          bcs.push_back(tmp) ;
          cout << "Boundary Face: " << zoneId << " " << name << endl ;
        }
	getCloseParen(s) ;
	getCloseParen(s) ;
      }
      break ;
    default:
      cerr << "unknown section " << code << endl ;
      cerr << "ignored" << endl ;
      getCloseParen(s) ;
      break ;
    }
  }

  if(dimension == 2) { // Need to extrude mesh
    // Create names for two symmetry planes ;
    int s1 = maxzoneid+1 ;
    int s2 = maxzoneid+2 ;
    VOG::BC_descriptor tmp ;
    tmp.name = "SymmetryZ1" ;
    tmp.id = s1 ;
    tmp.Trans = false ;
    bcs.push_back(tmp) ;
    cout << "Boundary Face: " <<tmp.id << " " << tmp.name << endl ;
    tmp.name = "SymmetryZ2" ;
    tmp.id = s2 ;
    bcs.push_back(tmp) ;
    cout << "Boundary Face: " <<tmp.id << " " << tmp.name << endl ;
    
    vector<pair<int,pair<int,int> > > etmp ;
    for(size_t i=0;i<edge_list.size();++i) {
      if(edge_list[i].second.first >=0) {
        int cell = edge_list[i].second.first ;
        pair<int,int> nl = pair<int,int>(edge_list[i].first.second,
                                         edge_list[i].first.first) ;
        etmp.push_back(pair<int,pair<int,int> > ( cell,nl) );
      }
      if(edge_list[i].second.second >=0) {
        int cell = edge_list[i].second.second ;
        pair<int,int> nl = pair<int,int>(edge_list[i].first.first,
                                         edge_list[i].first.second) ;
        etmp.push_back(pair<int,pair<int,int> > ( cell,nl) );
      }
    }
    sort(etmp.begin(),etmp.end()) ;

    // number of cells
    int cnt = 1 ;
    int c = etmp[0].first ;
    for(size_t i=1;i<etmp.size();++i) {
      if(c != etmp[i].first) {
        cnt++ ;
        c = etmp[i].first ;
      }
    }
    int cellcnt = cnt ;
    int nfaces = edge_list.size() + cellcnt*2 ;

    entitySet fdom = interval(0,nfaces-1) ;
    cl.allocate(fdom) ;
    cr.allocate(fdom) ;
    store<int> count ;
    count.allocate(fdom) ;

    cnt = 0 ;

    c = -1 ;
    for(size_t i=0;i<etmp.size();++i) {
      if(c != etmp[i].first) { // found new face
        c = etmp[i].first ;
        int xc = 1 ;
        for(size_t j=i+1;j<etmp.size();++j)
          if(etmp[j].first == c)
            xc++ ;
          else
            break ;
        i += xc-1 ;
        count[cnt] = xc ;
        cl[cnt] = c ;
        cr[cnt] = -s1 ;
        count[cnt+cellcnt] = xc ;
        cl[cnt+cellcnt] = c ;
        cr[cnt+cellcnt] = -s2 ;
        cnt++ ;
      }
    }
    if(cnt != cellcnt) {
      cerr << "counts don't match!" << endl;
    }
    cnt += cellcnt ;
    
    for(size_t i=0;i<edge_list.size();++i) {
      cl[cnt] = edge_list[i].second.first ;
      cr[cnt] = edge_list[i].second.second ;
      count[cnt] = 4 ;
      cnt++ ;
    }

    face2node.allocate(count) ;
    cnt = 0 ;
    
    for(size_t i=0;i<etmp.size();i+=count[cnt++]) {
      face2node[cnt][0] = etmp[i].second.first ;
      face2node[cnt][1] = etmp[i].second.second ;
      int sv = face2node[cnt][1] ;
      for(int k=2;k<count[cnt];++k) {
        int nsv = -1 ;
        for(int l=1;l<count[cnt];++l)
          if(etmp[i+l].second.first == sv)
            nsv = etmp[i+l].second.second ;
        if(nsv == -1) {
          cerr << "unable to follow edges around face!" << endl ;
          exit(-1) ;
        }
        face2node[cnt][k] = nsv ;
        sv = nsv ;
      }
      // Build face on opposite side
      for(int k=0;k<count[cnt];++k) {
        face2node[cnt+cellcnt][count[cnt]-k-1] = face2node[cnt][k]+n2dpos ;
      }
    }
    cnt = cellcnt*2 ;
    
    for(size_t i=0;i<edge_list.size();++i) {
      face2node[cnt][0] = edge_list[i].first.first ;
      face2node[cnt][1] = edge_list[i].first.second ;
      face2node[cnt][2] = edge_list[i].first.second + n2dpos ;
      face2node[cnt][3] = edge_list[i].first.first + n2dpos ;
      cnt++ ;
    }
                                          
    return 0 ;
  } else {
    face2node.allocate(count) ;
    entitySet dom = face2node.domain() ;
    FORALL(dom,ii) {
      int fsz = face2node[ii].size() ;
      for(int j=0;j<fsz;++j)
        face2node[ii][j] = face_info[ii][fsz-j-1]-1 ;
    } ENDFORALL ;
    return 0 ;
  }
}


int main(int ac, char *av[]) {
  using namespace Loci ;
  using namespace VOG ;

  bool optimize = true ;
  
  Loci::Init(&ac,&av) ;

  if(Loci::MPI_processes > 1) {
    cerr << "fluent2vog does not run multiprocessor" << endl ;
    Loci::Abort() ;
  }
  string Lref = "NOSCALE" ;
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
      Lref = "1 mm" ;
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
    cerr << "Usage: fluent2vog <options> <file>" << endl
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

  Loci::UNIT_type tp ;
  istringstream iss(Lref) ;
  iss >> tp ;
  double posScale = tp.get_value_in("meter") ;

  if(Loci::MPI_rank == 0) {
    cout << "input grid file units = " << tp ;
    if(posScale != 1.0) 
      cout << " = " << posScale << " meters " ;
    cout << endl ;
  }

  string filename = av[1] ;
  filename += ".msh" ;

  string outfile = av[1] ;
  outfile += ".vog" ;

  store<vector3d<double> > pos;
  Map cl, cr;
  multiMap face2node;
  vector<BC_descriptor> bcs ;

  // Read in Mesh
  int success = scanFluentFile(filename,pos,cl,cr,face2node,bcs) ;

  if(success != 0) {
    cerr << "Error reading '" << filename << "'" << endl ;
    Loci::Abort() ;
  }

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }
  
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0) 
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }
  
  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;


  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;
  
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  return 0 ;
}


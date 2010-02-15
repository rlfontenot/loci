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
#include <Loci.h>
#include "vogtools.h"
#include <iostream>
#include <cstdlib>
#include <string>

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
    char c = s.get() ;
    if(c == '(') {
      count++ ;
      while(count > 0) {
	c = s.get() ;
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

  map<int,int> zone_map ;
  store<int> count ;
  store<Loci::Array<int,4> > face_info ;
  ifstream s(filename.c_str(),ios::in) ;
  if(s.fail()) {
    return -1 ;
  }
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
      if(dimension != 3) {
	cerr << "only three dimensional mesh supported." << endl ;
	return -1;
      }
      break ;
    case 10:
      {
	int start,end,type,dim ;
	int zone = getZoneInfo(s,start,end,type,dim) ;
	if(type != 1) {
	  cerr << "can only handle type 1 nodes" << endl;
	  exit(-2) ;
	}
        if(zone == 0) { // allocate pos ;
          entitySet dom = interval(start-1,end-1) ;
          pos.allocate(dom) ;
          cout << "reading " << end << " nodes." << endl ;
        } else {
          if(pos.domain() == EMPTY) {
            cerr << "No zone 0 to set node size!" << endl ;
            return -1 ;
          }
          getOpenParen(s) ;
          for(int i=start-1;i<end;++i) {
            s >> pos[i].x >> pos[i].y >> pos[i].z ;
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
          cout << "reading " << end-start << " cells." << endl ;
      }
      break ;
    case 13: // Faces
      {
	int start, end, type, dim ;
	int zone = getZoneInfo(s,start,end,type,dim) ;
        if(zone == 0) {
          entitySet dom = interval(start-1,end-1) ;
          cl.allocate(dom) ;
          cr.allocate(dom) ;
          count.allocate(dom) ;
          face_info.allocate(dom) ;
          cout << "reading " << end-start << " faces." << endl ;
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
                zone_map[i] = 1 ;
              } else if(cl[i] == 0) {
                cl[i] = -zone ;
                zone_map[zone] = 1 ;
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
                zone_map[i] = 1 ;
              } else if(cl[i] == 0) {
                cl[i] = -zone ;
                zone_map[zone] = 1 ;
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
	getCloseParen(s) ;
      }
      break ;
    case 45: // Boundary names
      {
	getOpenParen(s) ;
	int zoneId = getDecimal(s) ;
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

  face2node.allocate(count) ;
  entitySet dom = face2node.domain() ;
  FORALL(dom,ii) {
    int fsz = face2node[ii].size() ;
    for(int j=0;j<fsz;++j)
      face2node[ii][j] = face_info[ii][fsz-j-1]-1 ;
  } ENDFORALL ;
  return 0 ;
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


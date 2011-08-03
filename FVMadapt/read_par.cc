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
#include <queue>
#include "read_par.h"
#include "defines.h"
#include <limits>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::ifstream;

double get_distance(const vect3d& p, const vect3d& p1, const vect3d& p2){
  if( norm(p1-p2) < NORMALIZE_ZERO_THRESHOLD) return norm(p - p1);
  double dotp1 = dot(p2-p1, p-p1);
  double dotp2 = dot(p1-p2, p-p2);
  if(dotp1 > 0 && dotp2 > 0) return norm(cross(p2-p1, p-p1))/norm(p2-p1);
  if(dotp1 <= 0) return norm(p-p1);
  if(dotp2 <= 0) return norm(p-p2);
  cerr<<"WARNING: get_distance reach dummy code" << endl;
  return norm(p-0.5*(p1+p2));
}
  
void readPar(string filename, vector<source_par>& sources){
  sources.clear();
 
  ifstream inFile;
  
  inFile.open(filename.c_str());
  if(!inFile){
    cerr <<"can not open " << filename << " for input" << endl;
    Loci::Abort();
  }
  
  int numSources = 0;
  inFile >> numSources;
  if(numSources <= 0){
    cerr<<"can not get numSources " << endl;
    Loci::Abort();
  }
  sources.resize(numSources);
  
 
  
  for(int i = 0; i < numSources; i++){

    source_par s;
    vect3d p1, p2;
    inFile >> p1.x >> p1.y >>p1.z;
    inFile >> p2.x >> p2.y >> p2.z;
    inFile >> s.r0 >> s.s0 >> s.r1 >> s.s1 >> s.a;
    
    s.p1 = p1;
    s.p2 = p2;
    sources[i] = s;
  }
  inFile.close();
    // s1.p1= vect3d(1.5, 0, 3);
//   s1.p2= vect3d(1.5, 1, 3);
//   s1.r0 = 0.2;
//   s1.s0 = 0.2;
//   s1.r1 = 0.5;
//   s1.s1 = 0.1; 
//   s1.a = 0.2;
 
   cerr<<"reading parameters" << endl;  
}

double get_spacing(const vect3d& p, const source_par& s){
  double sp;
  double r = get_distance(p, s.p1, s.p2);
  if(r <= s.r0) sp = s.s0;
  else if(r <= s.r1) sp = s.s1;
  else sp = s.s1*pow(r/s.r1, s.a);
  return sp;
}
double get_min_spacing(const vector<Node*>& nodes, const vector<source_par>& ss){
  //if(ss.size() == 1) return 0.25;
  double spacing = std::numeric_limits<double>::max(); 
  for(unsigned int np = 0; np < nodes.size(); np++){
    for(unsigned int i = 0; i < ss.size(); i++){
      spacing = min(spacing, get_spacing(nodes[np]->p, ss[i]));
    }
  }
  return spacing;
}


bool tag_cell(const vector<Node*>& nodes, const vector<source_par>& sources, double min_edge_len){


  //get bounding box of the cell
  double minx, miny, minz, maxx, maxy, maxz;
  minx=miny=minz=std::numeric_limits<float>::max();
  maxx=maxy=maxz=std::numeric_limits<float>::min();
  for(unsigned int i = 0; i < nodes.size(); i++){
    minx = min(minx, (nodes[i]->p).x);
    miny = min(miny, (nodes[i]->p).y);
    minz = min(minz, (nodes[i]->p).z);
    maxx = max(maxx, (nodes[i]->p).x);
    maxy = max(maxy, (nodes[i]->p).y);
    maxz = max(maxz, (nodes[i]->p).z);
  }
 
  for(unsigned int i = 0; i<sources.size(); i++){
    //if any feature is in the bounding box
    vect3d p1 = sources[i].p1;
    vect3d p2 = sources[i].p2;
    if((p1.x <= maxx && p1.x >= minx &&
        p1.y <= maxy && p1.y >= miny &&
        p1.z <= maxz && p1.z >= minz
        )||(p2.x <= maxx && p2.x >= minx &&
            p2.y <= maxy && p2.y >= miny &&
            p2.z <= maxz && p2.z >= minz)){
      if(min_edge_len>=2*sources[i].r1)return true;
    }
  }
      
  double min_spacing = get_min_spacing(nodes, sources);
  if(min_edge_len>=2*min_spacing)return true;
  else return false;
  
}

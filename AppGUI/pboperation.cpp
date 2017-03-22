#include "pboperation.h"
#define PI 3.14159265358979323846264338327950
#include <stdlib.h>
vect3d get_wireframe_center(const vector<vect3d> &pos,
                         const vector<int>  &trias)
{
  if(pos.size()==0) return vect3d(0.0, 0.0, 0.0);
  //find all edges of surface id, and it boundary tag is initialied to -2
  map<pair<int, int>, int> all_edges;
  int p[3];
  for(size_t i = 0; i < trias.size()/3; i++){

     p[0] = trias[3*i];
     p[1] = trias[3*i+1];
     p[2] = trias[3*i+2];
    for(int j =0; j < 3; j++){
      pair<int, int> edge;
      edge.first = min(p[j],p[(j+1)%3]);
      edge.second = max(p[j],p[(j+1)%3]);
      all_edges[edge]=-2;
    }
  }

  


  //for each edge in all_edges, reset its boundary tag
  int fid =0; 
  for(size_t i = 0; i < trias.size()/3; i++,fid++){
    p[0] = trias[3*i];
    p[1] = trias[3*i+1];
    p[2] = trias[3*i+2];
    for(int j =0; j < 3; j++){
      pair<int, int> edge;
      edge.first = min(p[j],p[(j+1)%3]);
      edge.second = max(p[j],p[(j+1)%3]);
      
      if(all_edges[edge] == -2){
         all_edges[edge] = fid;
      }else if(all_edges[edge] == fid){
        cerr << " same edge2face appear twice in find_bdry_nodes()" << endl;
        exit(0);
      }else{
        all_edges[edge] = -1;
      }
        
    }
  }
     
  //collect all boundary edges 
  vector<pair<pair<int,int>,int> > bdry_edges;
  for(map<pair<int, int>, int>::const_iterator p = all_edges.begin(); p != all_edges.end(); p++){
    if(p->second == -2){
      cerr<< "edge never be visited in find_bdry_nodes(): "<< (p->first).first << " " << (p->first).second << endl;
      exit(0);
    }else if(p->second != -1){
      bdry_edges.push_back(pair<pair<int, int>, int>(p->first, p->second));
      
    }
  }

  double total = 0;
  vect3d totalVec = vect3d(0.0, 0.0, 0.0);
  
  //compute the wireframe center
  for(unsigned int i = 0; i < bdry_edges.size(); i++){
    double len = norm(pos[bdry_edges[i].first.first]-pos[bdry_edges[i].first.second]);
    vect3d center = 0.5*(pos[bdry_edges[i].first.first]+pos[bdry_edges[i].first.second]);
    total+=len;
    totalVec += len*center;
  }

  return totalVec/total;
}


vect3d get_average_normal(const vector<vect3d> &pos,
                          const vector<int>  &trias){

  double total  = 0;
  vect3d totalVec = vect3d(0.0, 0.0, 0.0);
  vect3d p[3];
  for(size_t i = 0; i < trias.size()/3; i++){
    p[0] =pos[trias[3*i]];
    p[1] = pos[trias[3*i+1]];
    p[2] = pos[trias[3*i+2]];
    vect3d p01 = p[1]-p[0];
    vect3d p02 = p[2]-p[0];
    vect3d normal = cross(p01, p02);
    if(norm(normal) < 1e-34)continue;
    normal = normal/norm(normal);
    double area = dot(p01, p02);
    total += area;
    totalVec += area*normal;
  }

  return totalVec/total;
}
    
bool angleBetween(positions3d v1,positions3d v2, double& angle, vect3d& axis) {
  double norm1 = norm(v1);
  double norm2 = norm(v2);
  angle = 0;
  axis = positions3d(0, 0, 1);
  
  if(norm1 > 1e-34 && norm2 > 1e-34){
        
    // turn vectors into unit vectors 
    positions3d n1 = v1/norm1;
    positions3d n2 = v2/norm2; 
    angle = acos(dot(n1,n2)); 
	
    if( fabs(angle) > 1e-5 && fabs(angle-PI)>1e-5){
      axis = cross(n1, n2);
    }else if(fabs(angle) <= 1e-5){
      // if no noticable rotation is available return zero rotation
      // this way we avoid Cross product artifacts 
      axis = positions3d(0, 0, 1);
      angle = 0;
    }else if(fabs(angle-PI) < 1e-5){
      // in this case there are 2 lines on the same axis 
      // there are an infinite number of normals 
      // in this case. Anyone of these normals will be 
      // a valid rotation (180 degrees). so I pick one axis that is perpendicular to n1
      angle = PI;
      if(fabs(fabs(n1.z)-1.0)<1e-5)axis=positions3d(1, 0, 0);
      else axis = positions3d(-1*n1.y, n1.x, 0);
    }
  }else{
    return false;
  }
  return true;
}

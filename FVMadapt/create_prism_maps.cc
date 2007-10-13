#include <algorithm>
#include <cmath>
#include <iostream>
#include <Loci.h>
#include "hex_defines.h"

using std::cerr;
using std::cout;
using std::endl;
using Loci::storeRepP;


Array<Entity, 5> collect_prism_faces( const Entity* lower,
                                      const Entity* upper,
                                      const Entity* boundary_map,
                                      const Array<char, 5>& prism2face);
                                     

/* create_prismmaps
   * Create an ordering for the faces and nodes of a prism cell.
   *
   * prism2face	This is an array of 5 integers which give an ordering
   *		to the entities which appear in the cell2face map.
   * prism2node	This is an array of 6 integers which give an ordering
   *		to the entities which appear in the face2node maps
   *		of prism2triface[0] and prism2triface[1].  The first 3 entries
   *		refer to nodes of prism2face[0]  while the last
   *		3 entries refer to nodes of prism2triface[1] (zeta=1).
   *		coordinates of the prism cell as follows.
   */
  class create_prism_maps : public pointwise_rule {
    /** lower:	Mapping from cells to faces */
    const_multiMap lower;
    /** upper:	Mapping from cells to faces */
    const_multiMap upper;
    /** boundary_map: Mapping from cells to faces */
    const_multiMap boundary_map;
    /** face2node	Mapping from faces to nodes */
    const_multiMap face2node;
    store<Array<char,5> > orientCode;
    store<Array<char,5> > prism2face;
    store<Array<char,6> > prism2node;
  public:
    create_prism_maps() { 
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("face2node", face2node);
      name_store("prismOrientCode", orientCode);
      name_store("prism2face", prism2face);
      name_store("prism2node", prism2node);
      
      input("(upper,lower,boundary_map)->face2node");

      output("prismOrientCode");
      output("prism2face");
      output("prism2node");
      
      constraint("prisms");
    }
    virtual void compute(const sequence & seq) {
     
      do_loop(seq, this);
    }
    void calculate(Entity cc) {
      // Collect the entity designations of all nodes in the nodes array
      Entity nodes[6];
      // Collect the entity designations of all faces in the faces array
      Entity trifaces[2];
      Entity quadfaces[3];
      // Determine the relative orientations of the faces
      char orient_tri[2]; // value 'l' or 'r'
      char orient_quad[3]; //value 'l', or 'r'
      Array<char, 2> trirot;
      Array<char, 3> quadrot;
      Array<char, 2> prism2triface;
      Array<char, 3> prism2quadface;
      
      // Initialize topology deduction data structures
      int ntf=0;
      int nqf = 0;
      for (int f=0; f<lower.num_elems(cc); f++) {

        if(face2node.num_elems(lower[cc][f]) == 3){
          prism2triface[ntf] = f;
          orient_tri[ntf] = 'r';
          trirot[ntf] = 0;
          trifaces[ntf++] = lower[cc][f];
        }
        else if (face2node.num_elems(lower[cc][f]) == 4){
          prism2quadface[nqf] = f;
          orient_quad[nqf] = 'r';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = lower[cc][f];
        }
        
          
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }
     
      
      for (int f=0; f<upper.num_elems(cc); f++) {
        if(face2node.num_elems(upper[cc][f]) == 3){
          prism2triface[ntf] = f+6;
          orient_tri[ntf] = 'l';
          trirot[ntf] = 0;
          trifaces[ntf++] = upper[cc][f];
        }
        else if (face2node.num_elems(upper[cc][f]) == 4){
          prism2quadface[nqf] = f+6;
          orient_quad[nqf] = 'l';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = upper[cc][f];
        }
        
        
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }
      for (int f=0; f<boundary_map.num_elems(cc); f++) {

        if(face2node.num_elems(boundary_map[cc][f]) == 3){
          prism2triface[ntf] = f+12;
          orient_tri[ntf] = 'l';
          trirot[ntf] = 0;
          trifaces[ntf++] = boundary_map[cc][f];
        }
        else if (face2node.num_elems(boundary_map[cc][f]) == 4){
          prism2quadface[nqf] = f+12;
          orient_quad[nqf] = 'l';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = boundary_map[cc][f];
        }
        
        
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }

      if(ntf != 2 || nqf != 3){
        cerr<<"WARNING: error in clollecting prism faces" << endl;
        Loci::Abort();
      }
      //order the faces and everything according node_remap here

      
      // The first node of trifaces[0] will be nodes[0]
      nodes[0] = face2node[trifaces[0]][0];
      prism2node[cc][0] = 0;
      
      
      
      // Determine the orientation of the first face.
      // If orient[0]=='r', then faces[0] has an inward facing normal.
      if (orient_tri[0] == 'r') {
	// The second node of trifaces[0] will be nodes[1]
	nodes[1] = face2node[trifaces[0]][1];
	prism2node[cc][1] = 1;
	// The third node of faces[0] will be nodes[2]
	nodes[2] = face2node[trifaces[0]][2];
	prism2node[cc][2] = 2;
	// Set the rotation parameter to zero (no rotations, no flip)
	trirot[0] = 0;
      }
      // Otherwise, it's an outward facing normal, and we need to
      // reverse the direction of the first 3 nodes.
      else {
	// The 3rd node of faces[0] will be nodes[1]
	nodes[1] = face2node[trifaces[0]][2];
	prism2node[cc][1] = 2;
	// The 2nd node of faces[0] will be nodes[2]
	nodes[2] = face2node[trifaces[0]][1];
	prism2node[cc][2] = 1;
	// Set the rotation parameter to four (no rotations, one flip)
	trirot[0] = 4;
      }



      // Next, let's flag the other two faces which share nodes[0].
      // First, quadfaces[0]  will share nodes[0] and nodes[1]
      for (int f=0; f<3; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[quadfaces[f]][n] == nodes[0]) {
	    // If the next node in the face2node map is nodes[1]
	    // then this face has an outward facing normal.
	    if (face2node[quadfaces[f]][(n+1)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'l') {
		cerr << "WARNING[1a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // (n rotations)
             quadrot[0] = n ; // no flip 
	      // Make this face, quadfaces[0]
	      std::swap(prism2quadface[0], prism2quadface[f]);
	      std::swap(quadfaces[0], quadfaces[f]);
	      std::swap(orient_quad[0], orient_quad[f]);
	      // Store the other corners of the face in the node list
	      nodes[4] = face2node[quadfaces[0]][(n+2)%4];
	      nodes[3] = face2node[quadfaces[0]][(n+3)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	    // If the previous node in the face2node map is nodes[1]
	    // then this face has an inward facing normal.
	    else if (face2node[quadfaces[f]][(n+3)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'r') {
		cerr << "WARNING[1A]: Face orientation test failed!" << endl;
		throw -1;
	      }
              quadrot[0] = n + 4; //qh: if 'r', flip
	      // Make this face, faces[1]
	      std::swap(prism2quadface[0], prism2quadface[f]);
	      std::swap(quadfaces[0], quadfaces[f]);
	      std::swap(orient_quad[0], orient_quad[f]);
	      // Store the other corners of the face in the node list
	      nodes[4] = face2node[quadfaces[0]][(n+2)%4];
	      nodes[3] = face2node[quadfaces[0]][(n+1)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	  }
	}
      }
      // Next, quadfaces[2](2, 0, 3, 5)  will share nodes[2] and nodes[0]
      for (int f=1; f<3; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[quadfaces[f]][n] == nodes[2]) {
	    // If the next node in the face2node map is nodes[0]
	    // then this face has an outward facing normal.
	    if (face2node[quadfaces[f]][(n+1)%4] == nodes[0]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'l') {
		cerr << "WARNING[2a]: Face orientation test failed!" << endl;
		throw -1;
	      }
              // (n rotations,   no flip)
	      quadrot[2] = n;
	      // Make this face, faces[2]
	      std::swap(prism2quadface[2], prism2quadface[f]);
	      std::swap(quadfaces[2], quadfaces[f]);
	      std::swap(orient_quad[2], orient_quad[f]);
	      // Store the remaining corner of this face in the node list
              nodes[3] = face2node[quadfaces[2]][(n+2)%4];
              nodes[5] = face2node[quadfaces[2]][(n+3)%4];
              // Expedite the exit from these loops
              n = 4;
	      f = 3;
	    }
            // If the previous node in the face2node map is nodes[0]
	    // then this face has an inward facing normal.
            else if (face2node[quadfaces[f]][(n+3)%4] == nodes[0]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'r') {
		cerr << "WARNING[2A]: Face orientation test failed!" << endl;
		throw -1;
	      }
              
	      // (n rotations, one flip)
              quadrot[2] = n + 4;
	      // Make this face, faces[2]
	      std::swap(prism2quadface[2], prism2quadface[f]);
	      std::swap(quadfaces[2], quadfaces[f]);
	      std::swap(orient_quad[2], orient_quad[f]);
	      // Store the remaining corner of this face in the node list
              nodes[3] = face2node[quadfaces[2]][(n+2)%4];
              nodes[5] = face2node[quadfaces[2]][(n+1)%4];
              
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	  }
        }
      }
    
      // And trifaces[1]  will share nodes[3] and nodes[4] and nodes[5]
      for (int n=0; n<4; n++) {
        if (face2node[trifaces[1]][n] == nodes[3]) {
          // If the next node in the face2node map is nodes[4]
          // then this face has an outward facing normal
          if (face2node[trifaces[1]][(n+1)%3] == nodes[4]) {
            // Double check our deduction
            if (orient_tri[1] != 'l') {
              cerr << "WARNING[4a]: Face orientation test failed!" << endl;
              throw -1;
            }
            // nodes[5] is also shared with trifaces[1]
            if (face2node[trifaces[1]][(n+2)%3] != nodes[5]) {
              cerr << "WARNING[4b]: Coincident node test failed!" << endl;
              throw -1;
            }
              
            // (n rotations, no flip)
            trirot[1] = n;
            prism2node[cc][3] =  n;//qh: changed next 4 lines
            prism2node[cc][4] =  (n+1)%3;
            prism2node[cc][5] = (n+2)%3;
            n = 4;
          }
          // Otherwise the normal is inward facing
          else if (face2node[trifaces[1]][(n+1)%3] == nodes[5]) {
            // Double check our deduction
            if (orient_tri[1] != 'r') {
              cerr << "WARNING[4A]: Face orientation test failed!" << endl;
              throw -1;
            }
            // nodes[4] is also shared with trifaces[1]
            if (face2node[trifaces[1]][(n+2)%3] != nodes[4]) {
              cerr << "WARNING[4B]: Coincident node test failed!" << endl;
              throw -1;
            }
            
            // (n rotations, one flip)
	     
            trirot[1] = n + 4; 
            prism2node[cc][3] =  n;//qh: changed next 4 lines
            prism2node[cc][4] =  (n+2)%3;
            prism2node[cc][5] = (n+1)%3;
            n = 4;
            
          }
        }
      }
      // And finally quadfaces[1] should be the only one remaining
      for (int n=0; n<4; n++) {
	// quadfaces[1] should share nodes[1], nodes[2], nodes[5], and nodes[4]
	if (face2node[quadfaces[1]][n] == nodes[1]) {
	  // If the next node in the face2node map is nodes[2]
	  // then this face has an outward facing normal
	  if (face2node[quadfaces[1]][(n+1)%4] == nodes[2]) {
	    // Double check our deduction
	    if (orient_quad[1] != 'l') {
	      cerr << "WARNING[5a]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[5] is also shared with quadfaces[1]
	    if (face2node[quadfaces[1]][(n+2)%4] != nodes[5]) {
	      cerr << "WARNING[5b]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[4] is also shared with faces[3]
	    if (face2node[quadfaces[1]][(n+3)%4] != nodes[4]) {
	      cerr << "WARNING[5c]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, no flip)
	    quadrot[1] = n;
	    
          
	  }
	  // Otherwise the normal is inward facing
	  else if (face2node[quadfaces[1]][(n+1)%4] == nodes[4]) {
	    // Double check our deduction
	    if (orient_quad[1] != 'r') {
	      cerr << "WARNING[5A]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[6] is also shared with faces[2]
	    if (face2node[quadfaces[1]][(n+2)%4] != nodes[5]) {
	      cerr << "WARNING[5B]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[7] is also shared with faces[3]
	    if (face2node[quadfaces[1]][(n+3)%4] != nodes[2]) {
	      cerr << "WARNING[5C]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, one flip)
	   quadrot[1] = n + 4;
          }
	}
        
      }
      
      for(int i = 0; i < 2; i++){
        prism2face[cc][i] = prism2triface[i];
        orientCode[cc][i] = trirot[i];
      }
      for(int i = 0; i < 3; i++){
        prism2face[cc][i+2] = prism2quadface[i];
        orientCode[cc][i+2] = quadrot[i];
      }  


      //for print out
      
      if(false){
        cout << "cc " << cc << endl;
        cout << " prism2node: " ;
        for (int i = 0; i < 6; i++){
          cout << char(prism2node[cc][i] + '0') << " ";
        }
        cout << endl;
        
        for(int i = 0; i < 2; i++){
          cout <<  int(prism2triface[i]+ int(0)) << "; " << int(trirot[i]+int(0)) <<  "           ";
          
        }
        cout << endl;
        for(int i = 0; i < 3; i++){
          cout <<  int(prism2quadface[i]) << "; " << int(quadrot[i]) <<  "           ";
          
        }
        cout << endl;
        
        
      }
      
    }//calculate
    
  };

register_rule<create_prism_maps> register_create_prism_maps;


class determine_prism_fr : public pointwise_rule {
  const_Map cr;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > prism2face;
  const_blackbox<storeRepP>  node_remap;
  Map node_l2f;
  store<char> fr;
public:
  determine_prism_fr() {
      name_store("cr", cr);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("prism2face", prism2face);
      name_store("iface_remap", node_remap);
      name_store("fr", fr);
      
    
      input("cr->prism2face");
      input("cr->(lower,upper,boundary_map)");
      input("iface_remap");
      output("fr");
       constraint("cr->prisms");
    }
    virtual void compute(const sequence & seq) {
      if(seq.size()!=0){
        node_l2f = *node_remap;
        do_loop(seq, this);
      }
    }
  void calculate(Entity ff) {
    // Collect the entity designations of all faces in the faces array
    Array<Entity, 5> faces = collect_prism_faces(lower[cr[ff]].begin(), upper[cr[ff]].begin(),
                                                 boundary_map[cr[ff]].begin(), prism2face[cr[ff]]);
    
    // Store the original location of each face entity in faces
    int f;
    for ( f=0; f<5; f++) {
      if (node_l2f[ff] == node_l2f[faces[f]]) {
        fr[ff] = f;
        break;
      }
    }
      if(f ==5){
        cerr << "WARNING: can not find fr" << endl;
        Loci::Abort();
      }      
      //for print out
     //  if(false){
//         cout << "ff: " << ff << "fr[ff]: " <<char(fr[ff] + '0') <<endl ;
//       }
  }
};

register_rule<determine_prism_fr> register_determine_prism_fr;

class determine_prism_fl : public pointwise_rule {
  const_Map cl;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > prism2face;
   const_blackbox<storeRepP>  node_remap;
   Map node_l2f;
  store<char> fl;
public:
  determine_prism_fl() {
      name_store("cl", cl);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("prism2face", prism2face);
      name_store("face_remap", node_remap);
      name_store("fl", fl);
      
      input("face_remap");
      input("cl->prism2face");
      input("cl->(lower,upper,boundary_map)");
      
      output("fl");
      constraint("cl->prisms");
  }
  virtual void compute(const sequence & seq) {
    if(seq.size()!=0){
      node_l2f = *node_remap;
      do_loop(seq, this);
    }
  }
  void calculate(Entity ff) {
      // Collect the entity designations of all faces in the faces array
    Array<Entity, 5> faces = collect_prism_faces(lower[cl[ff]].begin(), upper[cl[ff]].begin(),
                                                 boundary_map[cl[ff]].begin(), prism2face[cl[ff]]);
    
    // Store the original location of each face entity in faces
    int f;
    for ( f=0; f<5; f++) {
      if (node_l2f[ff] == node_l2f[faces[f]]) {
        fl[ff] = f;
          break;
      }
    }
    if(f ==5){
      cerr << "WARNING: can not find fl" << endl;
        Loci::Abort();
    }
    // if(false){
//       cout << "ff: " << ff << "fr[ff]: " <<char(fl[ff] + '0') <<endl;
//     } 
  }
};

register_rule<determine_prism_fl> register_determine_prism_fl;


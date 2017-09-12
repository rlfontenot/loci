//#############################################################################
//#
//# Copyright 2015, Mississippi State University
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
#include <GLoci.h>
#include "defines.h"
using std::cout ;
using std::vector;
namespace Loci{
  int classify_cell(const gEntitySet& faces,const_gMultiMap &face2node);

  void parallelClassifyCell(fact_db &facts) {
    //get variables
    const_gMultiMap face2node ;
    face2node = facts.get_gvariable("face2node") ;
    const_gMultiMap upper ;
    upper = facts.get_gvariable("upper") ;
    const_gMultiMap lower ;
    lower = facts.get_gvariable("lower") ;
    const_gMultiMap boundary_map ;
    boundary_map = facts.get_gvariable("boundary_map") ;

    gConstraint geom_cells;
    geom_cells = facts.get_gvariable("geom_cells");
    gConstraint faces;
    faces = facts.get_gvariable("faces");
    
    //find actually face_dom
    gEntitySet face_dom = upper.image() + lower.image() + boundary_map.image();
    //expand map face2node
    if(Loci::MPI_processes > 1){
      vector<gEntitySet> init_ptn = g_all_collect_vectors<gEntity>(*faces);
      face2node.setRep(face2node.expand(face_dom, init_ptn));
    }

    gEntitySet hexcell;
    gEntitySet prism;
    gEntitySet gnrlcell;
    int elem_type;


  
    // Classify Cells
    GFORALL(*geom_cells, cc) {
      gEntitySet img = upper.image(cc) + lower.image(cc) + boundary_map.image(cc);
      elem_type = classify_cell(img,face2node) ;
      switch(elem_type) {
      case 1:
        hexcell += cc ; break ;
      case 2:
        prism += cc ; break ;
      default:
        gnrlcell += cc ;
      }
    } ENDGFORALL ;
     
    gConstraint hexcells;
    *hexcells = hexcell;
     
    gConstraint prisms;
    *prisms = prism;
     
    gConstraint gnrlcells;
    *gnrlcells = gnrlcell;

    
    facts.create_gfact("hexcells", hexcells, upper.get_domain_space());
    facts.create_gfact("prisms", prisms, upper.get_domain_space());
    facts.create_gfact("gnrlcells", gnrlcells, upper.get_domain_space());
        
    gEntitySet quadface;
    GFORALL(*faces, ff){
      if(face2node.num_elems(ff)==4)quadface += ff;
    }ENDGFORALL;
    gConstraint quadrangles;
    *quadrangles = quadface;
    facts.create_gfact("quadrangles", quadrangles, face2node.get_domain_space());
  }
  // int classify_cell(Entity *faces,int nfaces,const_multiMap &face2node);

  // void parallelClassifyCell(fact_db &facts) {
  //   //get variables
  //   const_multiMap face2node ;
  //   face2node = facts.get_gvariable("face2node") ;
  //   const_multiMap upper ;
  //   upper = facts.get_gvariable("upper") ;
  //   const_multiMap lower ;
  //   lower = facts.get_gvariable("lower") ;
  //   const_multiMap boundary_map ;
  //   boundary_map = facts.get_gvariable("boundary_map") ;

  //   constraint geom_cells;
  //   geom_cells = facts.get_gvariable("geom_cells");
  //   constraint faces;
  //   faces = facts.get_gvariable("faces");

  //   //find actually face_dom
  //   entitySet face_dom;
    
  //   FORALL(*geom_cells, cc) {
  //     for(int i=0;i<upper[cc].size();++i) face_dom += upper[cc][i];
  //     for(int i=0;i<lower[cc].size();++i) face_dom += lower[cc][i];
  //     for(int i=0;i<boundary_map[cc].size();++i) face_dom += boundary_map[cc][i];
  //   }ENDFORALL;
    
  //   //expand map face2node
  //   if(Loci::MPI_processes > 1){
  //     vector<entitySet> init_ptn = facts.get_init_ptn() ;
  //     entitySet out_of_dom = face_dom - face2node.domain();
  //     face2node.setRep(MapRepP(face2node.Rep())->expand(out_of_dom, init_ptn));
  //   }

  //   entitySet hexcell;
  //   entitySet prism;
  //   entitySet gnrlcell;
  //   int elem_type;


  
  //    // Classify Cells
  //    FORALL(*geom_cells, cc) {
  //      int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
  //      tmp_array<Entity> faces(nfaces) ;
  //      int cnt = 0 ;
  //      for(int i=0;i<upper[cc].size();++i)
  //        faces[cnt++] = upper[cc][i] ;
  //      for(int i=0;i<lower[cc].size();++i)
  //     faces[cnt++] = lower[cc][i] ;
  //      for(int i=0;i<boundary_map[cc].size();++i)
  //        faces[cnt++] = boundary_map[cc][i] ;
  //      elem_type = classify_cell(faces,nfaces,face2node) ;
  //      switch(elem_type) {
  //      case 1:
  //        hexcell += cc ; break ;
  //      case 2:
  //        prism += cc ; break ;
  //      default:
  //        gnrlcell += cc ;
  //      }
  //    } ENDFORALL ;
     
  //    constraint hexcells;
  //    *hexcells = hexcell;
     
  //    constraint prisms;
  //    *prisms = prism;
     
  //    constraint gnrlcells;
  //    *gnrlcells = gnrlcell;
     
  //    facts.create_fact("hexcells", hexcells);
  //    facts.create_fact("prisms", prisms);
  //    facts.create_fact("gnrlcells", gnrlcells);

  //    entitySet quadface;
  //    FORALL(*faces, ff){
  //      if(face2node.num_elems(ff)==4)quadface += ff;
  //    }ENDFORALL;

     
  //    constraint quadrangles;
  //    *quadrangles = quadface;
  //    facts.create_fact("quadrangles", quadrangles);
   
     
  // }


}




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
//#include <Loci.h>
#include <GLoci.h>
#include "defines.h"
using std::cout ;
using std::vector;
namespace Loci{
  int classify_cell(const gEntitySet& faces,const_gMultiMap &face2node);

  void parallelClassifyCell(gfact_db &facts) {
    //get variables
    const_gMultiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    const_gMultiMap upper ;
    upper = facts.get_variable("upper") ;
    const_gMultiMap lower ;
    lower = facts.get_variable("lower") ;
    const_gMultiMap boundary_map ;
    boundary_map = facts.get_variable("boundary_map") ;

    gConstraint geom_cells;
    geom_cells = facts.get_variable("geom_cells");
    gConstraint faces;
    faces = facts.get_variable("faces");
    
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

    
    facts.create_fact("hexcells", hexcells, upper.get_domain_space());
    facts.create_fact("prisms", prisms, upper.get_domain_space());
    facts.create_fact("gnrlcells", gnrlcells, upper.get_domain_space());
        
    gEntitySet quadface;
    GFORALL(*faces, ff){
      if(face2node.num_elems(ff)==4)quadface += ff;
    }ENDGFORALL;
    gConstraint quadrangles;
    *quadrangles = quadface;
    facts.create_fact("quadrangles", quadrangles, face2node.get_domain_space());
   }
}



